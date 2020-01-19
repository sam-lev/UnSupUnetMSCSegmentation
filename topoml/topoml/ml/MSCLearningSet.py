
#external imports
#from torch.utils.data import Dataset
import numpy as np
#internal imports
from topoml.topology.MSCSegmentation import MSCSegmentation
from topoml.ml.utils import get_split
from topoml.ui.LocalSetup import LocalSetup


class MSCLearningSets:
        def __init__(self, persistence_values=[], blur_sigmas=[], training_data_path=None, validation_data_path=None, test_data_path=None
                    ,training_write_path = None, validation_write_path = None, test_write_path = None
                     , training_set=None, validation_set=None, test_set=None):

            self.validation_set = validation_set
            self.training_set = training_set
            self.test_set = test_set


            self.validation_set = MSCSegmentation.geomsc_segment_images(persistence_values=persistence_values, blur_sigmas=blur_sigmas
                                                  , data_path=training_data_path, write_path=training_write_path)
            self.training_set = MSCSegmentation.geomsc_segment_images(persistence_values=persistence_values,blur_sigmas=blur_sigmas
                                                                      , data_path=validation_data_path, write_path=validation_write_path)
            self.test_set = MSCSegmentation.geomsc_segment_images(persistence_values=persistence_values,blur_sigmas=blur_sigmas
                                                                      , data_path=test_data_path, write_path=test_write_path)


# Dataset class for the retina dataset
# each item of the dataset is a tuple with three items:
# - the first element is the input image to be segmented
# - the second element is the segmentation ground truth image
# - the third element is a mask to know what parts of the input image should be used (for training and for scoring)
class RetinaDataset():#Dataset):
    def transpose_first_index(self, x, with_hand_seg=False):
        if not with_hand_seg:
            x2 = (np.transpose(x[0], [2, 0, 1]), np.transpose(x[1], [2, 0, 1]), np.transpose(x[2], [2, 0, 1]))
        else:
            x2 = (np.transpose(x[0], [2, 0, 1]), np.transpose(x[1], [2, 0, 1]), np.transpose(x[2], [2, 0, 1]),
                  np.transpose(x[3], [2, 0, 1]))
        return x2

    def __init__(self, retina_array, split='train', do_transform=False, with_hand_seg=False):
        self.with_hand_seg = with_hand_seg
        indexes_this_split = get_split(np.arange(len(retina_array), dtype=np.int), split)
        self.retina_array = [self.transpose_first_index(retina_array[i], self.with_hand_seg) for i in
                             indexes_this_split]
        self.split = split
        self.do_transform = do_transform

    def __getitem__(self, index):
        sample = [x for x in self.retina_array[index]] #torch.FloatTensor(x)
        """if self.do_transform:
            v_gen = RandomVerticalFlipGenerator()
            h_gen = RandomHorizontalFlipGenerator()
            t = Compose([
                v_gen,
                RandomVerticalFlip(gen=v_gen),
                h_gen,
                RandomHorizontalFlip(gen=h_gen),
            ])
            sample = t(sample)"""
        return sample

    def __len__(self):
        return len(self.retina_array)