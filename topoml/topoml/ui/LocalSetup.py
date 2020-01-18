import os

class LocalSetup:
        def __init__(self):
            self.paths_for_multivax=""" training_data_path = "/Users/multivax/Documents/PhD/4spring19/DeepLearning/DeepLearning/HW3/datasets/drive/DRIVE/training/images"
            #training_data_path = "/Users/multivax/Documents/PhD/4spring19/DeepLearning/DeepLearning/final_project/results/neuron_msc"
            testing_data_path = "/Users/multivax/Documents/PhD/4spring19/DeepLearning/DeepLearning/HW3/datasets/drive/DRIVE/test/images"
            train_write_path = "/Users/multivax/Documents/PhD/4spring19/DeepLearning/DeepLearning/HW3/datasets/drive/DRIVE/training/" # "/Users/multivax/Documents/PhD/4spring19/DeepLearning/DeepLearning/final_project/results/" #
            test_write_path = "/Users/multivax/Documents/PhD/4spring19/DeepLearning/DeepLearning/HW3/datasets/drive/DRIVE/test/"
            stare_training_data_path = "/Users/multivax/Documents/PhD/4spring19/DeepLearning/DeepLearning/HW3/datasets/stare/images"
            stare_train_write_path = "/Users/multivax/Documents/PhD/4spring19/DeepLearning/DeepLearning/HW3/datasets/stare/"
            """
            # Paths for Multivax
            self.project_base_path = "/home/sam/Documents/PhD/Research/GeoMSCxML/"
            self.training_data_path = os.path.join(self.project_base_path,"datasets","optics","drive","DRIVE","training","images")
            self.training_seg_data_path = os.path.join(self.project_base_path,"datasets","optics","drive","DRIVE","training","1st_manual")
            #training_data_path = "/Users/multivax/Documents/PhD/4spring19/DeepLearning/DeepLearning/final_project/results/neuron_msc"
            self.testing_data_path = os.path.join(self.project_base_path,"datasets","optics","drive","DRIVE","test","images")
            self.train_write_path = os.path.join(self.project_base_path,"datasets","optics","drive","DRIVE","training")
            # "/Users/multivax/Documents/PhD/4spring19/DeepLearning/DeepLearning/final_project/results/" #
            self.test_write_path = os.path.join(self.project_base_path,"datasets","optics","drive","DRIVE","test")
            self.stare_training_data_path = os.path.join(self.project_base_path,"datasets","optics","stare","images")
            self.stare_train_write_path = os.path.join(self.project_base_path,"datasets","optics","stare")
