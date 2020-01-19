import sys
import os
from skimage import io
import imageio

from .mscnn_segmentation import mscnn_segmentation
from .ml.MSCSample import SupervisedMSC
from .ui.ArcSelector import ArcSelector

#class with local file paths
from .ui.LocalSetup import LocalSetup

sys.path.append(os.getcwd())


LocalSetup = LocalSetup()
paths_for_multivax = LocalSetup.paths_for_multivax
# Paths for Multivax
project_base_path = LocalSetup.project_base_path
training_data_path = LocalSetup.training_data_path
training_seg_data_path  = LocalSetup.training_seg_data_path
testing_data_path = LocalSetup.testing_data_path
train_write_path = LocalSetup.train_write_path
test_write_path = LocalSetup.test_write_path
stare_training_data_path = LocalSetup.stare_training_data_path
stare_train_write_path = LocalSetup.stare_train_write_path

class MSCSegmentation:

    def __init__(self, msc=None, geomsc=None, labeled_segmentation=None, ridge=True, valley=False):
        self.msc = [] if msc is None else [msc]
        self.geomsc = [] if geomsc is None else [geomsc]
        self.labeled_segmentation = labeled_segmentation
        # msc after labeling arcs (updating arc.label_accuracy
        self.labeled_msc = []


    # use MSCTrainingSet to label msc segmentation for supervised learning
    # applying labels from hand segmented ground truth images.
    def label_msc(self, msc=None, geomsc=None, labeled_segmentation=None, invert=True):
        supMSC = SupervisedMSC(msc=msc, geomsc=geomsc, labeled_segmentation=labeled_segmentation)
        labeled_msc = supMSC.map_labeling(msc=msc, geomsc=geomsc, labeled_segmentation=labeled_segmentation,
                                          invert=invert)
        return labeled_msc


    # take raw 'image' used to compute msc and draw
    # msc valley and/or ridge edges displayed as
    # binary map if no labeling or color map if labeled based
    # off of ground truth hand segmentation determined by
    # percent overlap with labeled segmentation
    def draw_msc(self, image=None, image_filename=None, save_name=None, write_path='.'
                 ,msc=None, valley=True, ridge=False, invert=True, persistence=0):
        # check needed folders present else make
        p = persistence
        if image_filename is not None:
            image = io.imread(image_filename)#, as_gray=False)
            if save_name is None:
                save_name = image_filename.rsplit('/', 1)[1].rsplit('.', 1)[0]

        # ground truth segmentation
        raw_path = os.path.join(write_path, 'raw_images')
        # write folder path
        if not os.path.exists(os.path.join(project_base_path, train_write_path,'labeled_msc')):
            os.mkdir(os.path.join(project_base_path, train_write_path,'labeled_msc'))
        msc_seg_path = os.path.join(project_base_path, train_write_path,'labeled_msc')
        # persistence specific folder path
        if not os.path.exists(os.path.join(msc_seg_path, 'persistence_' + str(persistence))):
            os.mkdir(os.path.join(msc_seg_path, 'persistence_' + str(persistence)))
        msc_seg_path = os.path.join(msc_seg_path, 'persistence_' + str(persistence))

        # path name of msc to save to
        save_path_name = os.path.join(msc_seg_path, save_name + '_pers-' + str(persistence) + '-MSC.tif')

        msc_selector = ArcSelector(image=image, msc=msc, valley=True, ridge=False, invert=True)
        unsup_seg_in_arcs, unsup_seg_out_arcs, unsup_seg_pixels = msc_selector.draw_segmentation(save_path_name,
                                                                                                        msc=msc,
                                                                                                        invert=True)



    # iterate over various persistence values and compute the
    # morse smale segmentations for a number of images located at data_path
    # will create a folder at write path with raw images and another
    # fodler with msc segmentations for wach image at each persistence
    def msc_segment_images(self, persistence_values = [1], blur_sigmas = [3], data_path = '.'
                           , write_path = '.'):
        # check needed folders present else make
        if not os.path.exists(os.path.join(write_path, 'raw_images')):
            os.mkdir(os.path.join(write_path, 'raw_images'))
        # iterate through images and compute msc for each image
        # at various persistence values
        images = [f for f in os.listdir(data_path) if os.path.isfile(os.path.join(data_path, f)) and any(image_type in f.rsplit('.', 1)[1] for image_type in ['tif','gif','jpg','png','ppm'])]

        for img in images:
            for blur_sigma in blur_sigmas:
                for pers in persistence_values:
                    # construct msc object
                    mscnn = mscnn_segmentation()
                    mscnn.clear_msc()
                    msc = mscnn.compute_msc(image =  os.path.join(data_path,img), persistence = pers,
                                            blur_sigma=blur_sigma, write_path = write_path)
                    mscnn.construct_msc_from_image(image = os.path.join(data_path,img),
                                                   write_path = write_path, persistence = pers,
                                                   blur_sigma = blur_sigma)
                    self.msc.append(msc)


    # uses mscnn_segmentation class to construct geometric MSC
    # over given images and given persistences. Geometric in that
    # vertices added at edge intersections including ridge lines
    # and not only , min/max
    def geomsc_segment_images(self, persistence_values = [1], blur_sigmas = [3], data_path = '.',
                              write_path = '.', labeled_segmentation=None, label=True, save_binary_seg=False):
        # check needed folders present else make
        if not os.path.exists(os.path.join(write_path, 'raw_images')):
            os.mkdir(os.path.join(write_path, 'raw_images'))
        # iterate through images and compute msc for each image
        # at various persistence values
        images = sorted([f for f in os.listdir(data_path) if os.path.isfile(os.path.join(data_path, f)) and any(image_type in f.rsplit('.', 1)[1] for image_type in ['tif','gif','jpg','png','ppm'])])
        seg_images = sorted([f for f in os.listdir(training_seg_data_path)])

        msc_segmentations = []
        for im_num, img in enumerate(images):
            for blur_sigma in blur_sigmas:
                for pers in persistence_values:
                    # construct msc object
                    mscnn = mscnn_segmentation()
                    mscnn.clear_msc()
                    #compute geometric msc
                    msc = mscnn.compute_geomsc(image_filename =  os.path.join(data_path,img)
                                               , persistence = pers
                                               , blur=True
                                               , blur_sigma=blur_sigma
                                               , grey_scale=True
                                               , write_path = write_path
                                               , scale_intensities=False
                                               , augment_channels = [])

                    if labeled_segmentation is None and label:
                        labeled_segmentation = imageio.imread(os.path.join(training_seg_data_path,seg_images[im_num]))
                    if label:
                        labeled_msc = self.label_msc(geomsc=msc,labeled_segmentation=labeled_segmentation)
                        mscnn.msc = labeled_msc
                        msc = labeled_msc
                    # compute geomsc over image

                    mscnn.construct_geomsc_from_image(image_filename = os.path.join(data_path,img)
                                                      , write_path = write_path
                                                      , persistence = pers
                                                      , blur_sigma = blur_sigma
                                                      ,binary=save_binary_seg)

                    self.msc.append(msc)
                    msc_segmentations.append(msc)
        return msc_segmentations




                #draw_msc(image = None, image_filename = os.path.join(data_path,img), save_name=None,
                #         write_path = '.',msc = labeled_msc, valley = True, ridge = False,
                #         invert=True, persistence=pers)




    # declarations for running examples
    """persistence_values = [0.01]#[10, 12, 15, 20 , 23, 25, 30] # below 1 for GeoMSC
    blur_sigmas = [2.0] """
    #traininig_data_path = "/Users/multivax/Documents/PhD/4spring19/DeepLearning/DeepLearning/final_project/results/neuron_msc"
    #train_write_path = "/Users/multivax/Documents/PhD/4spring19/DeepLearning/DeepLearning/final_project/results/"

    # run the MSC segmentation over images
    """msc_segment_images(persistence_values = persistence_values, blur_sigma = blur_sigma,
                       data_path = training_data_path, write_path = train_write_path)"""

    #msc_segment_images(persistence_values = persistence_values, blur_sigma = blur_sigma,
    #                   data_path = stare_training_data_path, write_path = train_write_path)

    #msc_segment_images(persistence_values = persistence_values, blur_sigma = blur_sigma,
    #                   data_path = testing_data_path, write_path = test_write_path)

    # run the Geometric MSC segmentation over images
    """geomsc_segment_images(persistence_values = persistence_values
                          , blur_sigma = blur_sigmas
                          , data_path = training_data_path
                          , write_path = train_write_path)"""
