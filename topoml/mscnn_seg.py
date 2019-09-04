import sys
import os

from topoml.mscnn_segmentation import mscnn_segmentation

sys.path.append(os.getcwd())



training_data_path = "/Users/multivax/Documents/PhD/4spring19/DeepLearning/DeepLearning/HW3/datasets/drive/DRIVE/training/images"
#training_data_path = "/Users/multivax/Documents/PhD/4spring19/DeepLearning/DeepLearning/final_project/results/neuron_msc"

testing_data_path = "/Users/multivax/Documents/PhD/4spring19/DeepLearning/DeepLearning/HW3/datasets/drive/DRIVE/test/images"
train_write_path = "/Users/multivax/Documents/PhD/4spring19/DeepLearning/DeepLearning/HW3/datasets/drive/DRIVE/training/" # "/Users/multivax/Documents/PhD/4spring19/DeepLearning/DeepLearning/final_project/results/" #

test_write_path = "/Users/multivax/Documents/PhD/4spring19/DeepLearning/DeepLearning/HW3/datasets/drive/DRIVE/test/"


stare_training_data_path = "/Users/multivax/Documents/PhD/4spring19/DeepLearning/DeepLearning/HW3/datasets/stare/images"

stare_train_write_path = "/Users/multivax/Documents/PhD/4spring19/DeepLearning/DeepLearning/HW3/datasets/stare/"



# iterate over various persistence values and compute the
# morse smale segmentations for a number of images located at data_path
# will create a folder at write path with raw images and another
# fodler with msc segmentations for wach image at each persistence
def msc_segment_images(persistence_values = [1], blur_sigma = 3, data_path = '.', write_path = '.'):
    # check needed folders present else make
    if not os.path.exists(os.path.join(write_path, 'raw_images')):
        os.mkdir(os.path.join(write_path, 'raw_images'))
    # iterate through images and compute msc for each image
    # at various persistence values
    images = [f for f in os.listdir(data_path) if os.path.isfile(os.path.join(data_path, f)) and any(image_type in f.rsplit('.', 1)[1] for image_type in ['tif','gif','jpg','png','ppm'])]
        
    for img in images:
        for pers in persistence_values:

            # construct msc object
            mscnn = mscnn_segmentation()
            mscnn.clear_msc()
            
            msc = mscnn.compute_msc(image =  os.path.join(data_path,img), persistence = pers, blur_sigma=blur_sigma, write_path = write_path)

            mscnn.construct_msc_from_image(image = os.path.join(data_path,img),  write_path = write_path, persistence = pers, blur_sigma = blur_sigma)


persistence_values = [7]#[10, 12, 15, 20 , 23, 25, 30]
blur_sigma = 2.2

#traininig_data_path = "/Users/multivax/Documents/PhD/4spring19/DeepLearning/DeepLearning/final_project/results/neuron_msc"
#train_write_path = "/Users/multivax/Documents/PhD/4spring19/DeepLearning/DeepLearning/final_project/results/"
msc_segment_images(persistence_values = persistence_values, blur_sigma = blur_sigma,
                   data_path = training_data_path, write_path = train_write_path)

#msc_segment_images(persistence_values = persistence_values, blur_sigma = blur_sigma,
#                   data_path = stare_training_data_path, write_path = train_write_path)

#msc_segment_images(persistence_values = persistence_values, blur_sigma = blur_sigma,
#                   data_path = testing_data_path, write_path = test_write_path)
