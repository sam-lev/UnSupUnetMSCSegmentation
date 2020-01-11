
import os
import sys
from time import time
import shutil

#for computing msc of images
import subprocess
from topoml.topology.read_msc import MSC
from topoml.topology.geometric_msc import GeoMSC

import datetime

import numpy as np

from skimage import io
from skimage.util import invert

from PIL import Image

import json
from networkx.readwrite import json_graph

from topoml.ml.utils import gaussian_fit
from topoml.NeuronTracer import NeuronTracer
from topoml.topology.utils import (
    get_pixel_values_from_arcs,
    is_ridge_arc,
    is_valley_arc,
)
""" give option to compute features from given selection or
    pre defined graph"""
from topoml.image.feature import (
    mean_filter,
    variance_filter,
    median_filter,
    minimum_filter,
    maximum_filter,
    gaussian_blur_filter,
    difference_of_gaussians_filter,
    sobel_filter,
    laplacian_filter,
    neighbor_filter,
)
from topoml.ui.colors import blue, light_green
from topoml.ui.vis import show_image

from topoml.topology.utils import build_msc
from topoml.topology.utils import build_geometric_msc
from topoml.ui.ArcSelector import ArcSelector
#To build MSC and select arcs
from topoml.ArcNeuronTracer import ArcNeuronTracer

from topoml.image.utils import blur_and_save, scale_intensity
from topoml.image.utils import augment_channels as aug_chan

class mscnn_segmentation:
    def __init__(self):
        self.image = None
        self.raw_image = None
        self.msc_selector = None
        self.fname_raw = None
        self.msc = None
        self.arcs = None
        self.msc_segmentaion = ArcNeuronTracer(construct_msc = False)
        if self.msc:
            self.msc.arcs = None
            self.msc.in_arcs = None
        self.msc = None
        self.msc_selector = None
        self.fname_raw = None
        self.in_arcs = set()
        self.out_arcs = set()
        self.arcs = None

    def clear_msc(self):
        self.msc_segmentaion = ArcNeuronTracer(construct_msc = False)
        if self.msc:
            self.msc.arcs = None
            self.msc.in_arcs = None
        self.msc = None
        self.msc_selector = None
        self.fname_raw = None
        self.in_arcs = set()
        self.out_arcs = set()
        self.arcs = None

    def compute_msc(self, image='./data/optic1.tif', write_path = '',fname_raw = None, persistence = 1, blur = True, blur_sigma = 2, delete_msc_files = True, scale_intensities = False):

        if not os.path.exists(os.path.join(write_path, 'raw_images')):
            os.mkdir(os.path.join(write_path, 'raw_images'))
        im = image
        img_name =  image.rsplit('/', 1)[1].rsplit('.', 1)[0]
        if fname_raw:
            fname = image
            fname_raw = os.path.join(write_path,'raw_images', fname.rsplit('.', 1)[0] + '.raw')
            image = io.imread(fname, as_gray=True, flatten=True)
            image.astype('float32').tofile(fname_raw)
        
        if blur:
            raw_image = io.imread(image, as_gray=True)
            if write_path:
                raw_path = os.path.join(write_path.rsplit(".", 1)[0],'raw_images',img_name)
            else:
                raw_path = image.rsplit(".", 1)[0]
            image, fname_raw = blur_and_save(raw_image,
                                                  raw_path+'PERS'+str(persistence), blur_sigma)
        if scale_intensities:
            img = io.imread(im, as_gray=True)
            image, fname_raw = scale_intensity(img,
                                               raw_path+'PERS'+str(persistence), scale_range=(0,255))
                
        proc = subprocess.Popen(['sh',
                                 os.path.join('..','MSCvisus','MSCTest2DViewer','test_2d_viewer'),
                                    fname_raw,
                                    str(image.shape[0]),
                                    str(image.shape[1]),
                                    "1",
                                    "1",
                                    str(persistence)],
                                   stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        #print(stdout.decode("utf-8"))
        #print(stderr.decode("utf-8"))
        msc = MSC()
        msc.read_from_file(fname_raw)
        self.msc = msc

        if delete_msc_files:
            raw_folder = os.path.join(write_path.rsplit(".", 1)[0],'raw_images')
            for the_file in os.listdir(raw_folder):
                file_path = os.path.join(raw_folder, the_file)
                try:
                    if os.path.isfile(file_path):
                        os.unlink(file_path)
                        #elif os.path.isdir(file_path): shutil.rmtree(file_path)
                except Exception as e:
                    print(e)
        return msc

    def compute_geomsc(self, image_filename='./data/optic1.tif'
                       , x=10, y=10
                       , persistence = 1
                       , blur = True, blur_sigma = 2
                       , write_path = '', delete_msc_files = False
                       , fname_raw = None
                       , invert_image = True
                       , grey_scale = True
                       , scale_intensities = False
                       , augment_channels = []):
        
        if not os.path.exists(os.path.join(write_path, 'raw_images')):
            os.mkdir(os.path.join(write_path, 'raw_images'))
        img_name =  image_filename.rsplit('/', 1)[1].rsplit('.', 1)[0]

        image = image_filename
        if fname_raw:
            fname = image_filename
            fname_raw = os.path.join(write_path,'raw_images', fname.rsplit('.', 1)[0] + '.raw')
            image = io.imread(fname, as_gray=grey_scale, flatten=True)
            if invert_image:
                image = invert(image)
            image.astype('float32').tofile(fname_raw)

        if scale_intensities:
            raw_image = io.imread(image, as_gray=grey_scale)
            if invert_image:
                raw_image = invert(raw_image)
            if write_path:
                raw_path = os.path.join(write_path.rsplit(".", 1)[0],'raw_images',img_name)
            else:
                raw_path = image.rsplit(".", 1)[0]
            image, fname_raw = scale_intensity(raw_image,
                                               raw_path+'PERS'+str(persistence), scale_range=(0,255))

        if len(augment_channels) > 0:
            im_name = image
            im = Image.open(image)#cv2.imread(image)
            im = np.array(im)
            #im = cv2.cvtColor(im, cv2.COLOR_BGR2HSV)
            write_dir =  os.path.join(write_path.rsplit(".", 1)[0],'augmented_images')
            if write_path:
                raw_path = os.path.join(write_dir,img_name)
            else:
                raw_path = image.rsplit(".", 1)[0]
            if not os.path.exists(write_dir):
                os.makedirs(write_dir)
            im, image = aug_chan(im, raw_path+'.'+im_name.split('.')[1], augment_channels)
            
        if blur:
            raw_image = io.imread(image, as_gray=grey_scale)
            if invert_image:
                raw_image = invert(raw_image)
            if write_path:
                raw_path = os.path.join(write_path.rsplit(".", 1)[0],'raw_images',img_name)
            else:
                raw_path = image.rsplit(".", 1)[0]
            image, fname_raw = blur_and_save(raw_image,raw_path+'PERS'+str(persistence),  blur_sigma=blur_sigma, grey_scale=grey_scale)
        else:
            raw_image = io.imread(image, as_gray=grey_scale)
            if invert_image:
                raw_image = invert(raw_image)
            if write_path:
                raw_path = os.path.join(write_path.rsplit(".", 1)[0],'raw_images',img_name)
            else:
                raw_path = image.rsplit(".", 1)[0]
            fname_raw = raw_path+'PERS'+str(persistence)+".raw"
            image = raw_image.astype("float32")
            image.tofile(fname_raw)
        # Sanity Check
        ##print("SHAPE ##### : ", image.shape)
        ##print("fname raw ###### : ", fname_raw)
        ##io.imshow(image)
        ##from matplotlib import pyplot as plt
        ##plt.show()
        
        if fname_raw is None:
            fname_raw = image_filename
        # hard coded path, different for other systems
        # will remain like this while still editing
        # geometric morse smale complex. After which
        # we will use pybind to eliminate the use of the shell
        # utilizing Attila's script directly.
        #/home/sam/Documents/PhD/Research/GradIntegrator/build/extract2dridgegraph
        #print("%%%%%%%%%%%% image file path %%%%%%: ", fname_raw)
        
        starting_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "..", "..","..")
        starting_dir = os.path.join( starting_dir
                                     ,"GradIntegrator"
                                     ,"build"
                                     ,"extract2dridgegraph")

        #print("list ex dir  %%%%%%%%%%", os.listdir(starting_dir))

        proc = subprocess.Popen(
            [   os.path.join('.',starting_dir,"extract2dridgegraph"), #note: '.' needed to run executable.
                fname_raw,
                str(image.shape[1]),
                str(image.shape[0]),
                str(persistence),
            ],
            shell=False,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stdout, stderr = proc.communicate()
        msc = GeoMSC()
        msc.read_from_file(fname_raw)
        self.msc = msc
        
        if delete_msc_files:
            raw_folder = os.path.join(write_path.rsplit(".", 1)[0],'raw_images')
            for the_file in os.listdir(raw_folder):
                file_path = os.path.join(raw_folder, the_file)
                try:
                    if os.path.isfile(file_path):
                        os.unlink(file_path)
                        #elif os.path.isdir(file_path): shutil.rmtree(file_path)
                except Exception as e:
                    print(e)
        return msc

    # msc for seg and image pre-computed and
    # this method overlays the valley/ridge arcs t0
    # create a segmentation
    def construct_msc_from_image(self, image ="./data/optic1.tif",  write_path = '.', persistence = 10, blur_sigma = 3, msc = None):
        #selector = ArcSelector(tracer.raw_image, tracer.msc, valley=False)
        image = image
        raw_image = io.imread(image, as_gray=True)
        p = persistence
        img_name =  image.rsplit('/', 1)[1].rsplit('.', 1)[0]
        #msc_w_path =  os.path.join(cwd,'data','msc_graphs')
        #group_msc_file =  os.path.join(msc_w_path,load_msc+'-MSC.csv')

        #ground truth segmenttion
        raw_path = os.path.join(write_path, 'raw_images')
        
        if not os.path.exists(os.path.join(write_path, 'ground_truth_seg')):
            os.mkdir(os.path.join(write_path, 'ground_truth_seg'))
        gt_msc_seg_path = os.path.join(write_path, 'ground_truth_seg')
        
        if not os.path.exists(os.path.join(write_path, 'msc_seg')):
            os.mkdir(os.path.join(write_path, 'msc_seg'))

        msc_seg_path = os.path.join(write_path, 'msc_seg')
        if not os.path.exists(os.path.join(msc_seg_path, 'persistence_'+str(persistence))):
            os.mkdir(os.path.join(msc_seg_path, 'persistence_'+str(persistence)))
        msc_seg_path = os.path.join(msc_seg_path, 'persistence_'+str(persistence))
            
        seg_img =  os.path.join(write_path, 'ground_truth_seg', img_name + '_seg.gif')
        img_msc =  os.path.join(msc_seg_path, img_name + '_pers-'+str(persistence)+'-MSC.tif')

        # if we want to compute arcneurontracer however compute msc in this class
        # obtain the saved or compute the msc for image and ground truth
        self.msc_segmentaion = ArcNeuronTracer(construct_msc = False)
        # (notice Switch in valley and ridge)
        #segmentation =  ArcNeuronTracer(seg_img, valley=True, ridge=False)
        
        #arce select but provides msc arcs found
        msc_selector = ArcSelector(image = raw_image, msc = self.msc,
                                          valley=True, ridge = False, invert=True)
        #seg_msc_selector =  ArcSelector(segmentation.raw_image, segmentation.msc, valley=False, ridge = True)
        unsup_seg_in_arcs, unsup_seg_out_arcs, unsup_seg_pixels = msc_selector.draw_binary_segmentation(img_msc, msc = self.msc, invert=False)#write_image(img_msc)
        
        """seg_img_msc =  seg_img.rsplit('.', 1)[0] + '_msc.tif'
        seg_in_arc, seg_out_arcs, seg_pixels = seg_msc_selector.draw_binary_segmentation(seg_img_msc) #write_image(seg_img_msc)
        """

    def construct_geomsc_from_image(self, image_filename ="./data/optic1.tif",  write_path = '.', persistence = 0, blur_sigma = 2, msc = None):
        #selector = ArcSelector(tracer.raw_image, tracer.msc, valley=False)
        
        raw_image = io.imread(image_filename, as_gray=True)
        p = persistence
        img_name =  image_filename.rsplit('/', 1)[1].rsplit('.', 1)[0]
        #msc_w_path =  os.path.join(cwd,'data','msc_graphs')
        #group_msc_file =  os.path.join(msc_w_path,load_msc+'-MSC.csv')

        #ground truth segmenttion
        raw_path = os.path.join(write_path, 'raw_images')
        
        if not os.path.exists(os.path.join(write_path, 'ground_truth_seg')):
            os.mkdir(os.path.join(write_path, 'ground_truth_seg'))
        gt_msc_seg_path = os.path.join(write_path, 'ground_truth_seg')
        
        if not os.path.exists(os.path.join(write_path, 'msc_seg')):
            os.mkdir(os.path.join(write_path, 'msc_seg'))

        msc_seg_path = os.path.join(write_path, 'msc_seg')
        if not os.path.exists(os.path.join(msc_seg_path, 'persistence_'+str(persistence))):
            os.mkdir(os.path.join(msc_seg_path, 'persistence_'+str(persistence)))
        msc_seg_path = os.path.join(msc_seg_path, 'persistence_'+str(persistence))
            
        seg_img =  os.path.join(write_path, 'ground_truth_seg', img_name + '_seg.gif')
        img_msc =  os.path.join(msc_seg_path, img_name + '_pers-'+str(persistence)+'-MSC.tif')

        # if we want to compute arcneurontracer however compute msc in this class
        # obtain the saved or compute the msc for image and ground truth
        #self.msc_segmentaion = ArcNeuronTracer(construct_msc = True, geometric_msc = True)

        # (notice Switch in valley and ridge)
        #segmentation =  ArcNeuronTracer(seg_img, valley=True, ridge=False)
        
        #arce select but provides msc arcs found
        msc_selector = ArcSelector(image = raw_image
                                   , msc = self.msc
                                   , valley=True
                                   , ridge = False
                                   , invert=True)
        #seg_msc_selector =  ArcSelector(segmentation.raw_image, segmentation.msc, valley=False, ridge = True)
        unsup_seg_in_arcs, unsup_seg_out_arcs, unsup_seg_pixels = msc_selector.draw_binary_segmentation(img_msc                                                                        , msc = self.msc                                                                , invert=True
                                                                                                        , reshape_out=False
                                                                                                        ,dpi = 158)#write_image(img_msc)
        
        """seg_img_msc =  seg_img.rsplit('.', 1)[0] + '_msc.tif'
        seg_in_arc, seg_out_arcs, seg_pixels = seg_msc_selector.draw_binary_segmentation(seg_img_msc) #write_image(seg_img_msc)"""
