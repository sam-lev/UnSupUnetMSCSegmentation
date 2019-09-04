"""#!/home/sci/samlev/bin/bin/python3
"""
#SBATCH --time=21-00:00:00 # walltime, abbreviated by -t
#SBATCH --mem=200G
#SBATCH -o slurm-%j.out-%N # name of the stdout, using the job number (%j) and the first node (%N)

#SBATCH -e slurm-%j.err-%N # name of the stderr, using the job number (%j) and the first node (%N)

#SBATCH --gres=gpu:1
# needed for adjusting gpu allocation when running on medusa @ sci
import sys
import os
from time import time
sys.path.append(os.getcwd())
slurm = False

from topoml.msc_gnn import msc_gnn

# construct msc gnn handler
gnn = msc_gnn()
# select msc arcs or load msc from file




cvt_1 = 'cvt_ridge_arcs_n=2_k=20_train'#load training graphs
cvt_2 = 'cvt_ridge_arcs_n=4_k=10'
cvt_3 = 'cvt_ridge_arcs_n=20_k=5'
n2k30 = 'n-2+k-30'
n2k30v2 = 'n-2+k-30+v-2'
n1k50 = 'n-1_k-50'
n50k1 = 'n-50+k-1'
full_msc = "test_ridge_arcs"

cvt= (2, 30) #for cvt sampling, (number samples, rings in group(hops))

density_arc = 3.92 #nodes/hop
density_bg = 6.2 #nodes/hop

# train/test model by loading labeled msc or existing gnn object
# must have labeled train graph subset to test graph
train_group_name = n2k30v2
test_group_name = full_msc
saved_msc=full_msc

# hyper-param for gnn
learning_rate = 0.3
polarity = 9
weight_decay = 0.001
epochs = 10
depth = 3

# Random walks used to determine node pairs for unsupervised loss.
# to make a random walk collection use ./topoml/graphsage/utils.py
#example run: python3 topoml/graphsage/utils.py ./data/json_graphs/test_ridge_arcs-G.json ./data/random_walks/full_msc_n-1_k-40
##load_walks='full_msc_n-1_k-40'
##load_walks='full_msc_n-1_k-50_Dp-10_nW-100'
load_walks='full_msc_n-1_k-50_Dp-10_nW-50'

# file name for embedding of trained model 

embedding_name = 'n-'+str(cvt[0])+'_k-'+str(cvt[1])+'_lr-'+str(learning_rate)+'Plrty-'+str(polarity)+'_epochs-'+str(epochs)+'_depth-'+str(depth)+'trainGrp-'+train_group_name

print('embedding will be named: ',embedding_name)

# Do hand selection of arcs with window
"""gnn.select_msc( train_or_test='train', group_name= train_group_name, bounds=[[0,450],[0,300]], write_json=True, write_msc=True, load_msc=False)
"""

# Load msc from pre-labeled (train/test) file
"""gnn.select_msc( train_or_test='train', group_name= test_group_name, load_msc = saved_msc, write_json =True, write_msc = False, bounds=False)
"""

# Perform cvt sampling for training set  and load full msc groundtruth for testing. If train graph not saved can use as class variables in gnn object
gnn.select_msc( group_name= train_group_name, load_msc = saved_msc,
                write_json = True, write_msc = False, bounds=False,
                cvt=cvt)

print('training/validation set selected')

# construct unsupervised gnn model
aggregators = ['graphsage_maxpool', 'gcn', 'graphsage_seq',
               'graphsage_maxpool', 'graphsage_meanpool',
               'graphsage_seq', 'n2v']
aggregator = 'graphsage_meanpool'

print('... Beginning training with aggregator: ', aggregator)

gnn.unsupervised(aggregator=aggregator, slurm = slurm)


#------------------------TRAINING-----------------------#
start = time()

# no param needed, uses selected however can pass all or some of graph data structures needed for training.
"""gnn.train()
"""

# if loading labled training graph from file
gnn.train(train_prefix = train_group_name, embedding_name = embedding_name, load_walks=load_walks, learning_rate=learning_rate, epochs=epochs , weight_decay=weight_decay, polarity=polarity, depth=depth)


# if from cvt selection saved in class variables of msc_gnn
"""gnn.train( embedding_name = embedding_name, load_walks=load_walks, learning_rate=learning_rate, epochs=epochs, weight_decay=weight_decay, depth=depth)
"""

end = time()
print("-------TRAINING TIME: ", str(end-start)," . Begining Classification...")


#----------------------CLASSIFICATION------------------#
start = time()

print("training done")
print("beggining classification of unknown portion of MSC")
# if embedding graph made with test/train set the same (and named the same)
"""gnn.classify( embedding_prefix=embedding_name, aggregator = aggregator, learning_rate=learning_rate)
"""
#Can also pass graph if contained in gnn
"""gnn.classify()
"""

# Hand select test graph
"""gnn.select_msc( train_or_test='test', group_name= train_group_name,bounds=[[451,700],[301,600]], write_json=True, write_msc=True, load_msc=False)
"""

#if loading graph (test/train) for classification
gnn.classify( test_prefix=train_group_name, trained_prefix=train_group_name, embedding_prefix=embedding_name, aggregator = aggregator, learning_rate=learning_rate)

end =time()
print("------CLASSIFICATION TIME: ",str(end-start))


# Show labaling assigned by trained model
gnn.show_gnn_classification(pred_graph_prefix=embedding_name, train_view=False)#embedding_name)

# see train and val sets, must put in directory log-dir and make new folder
# with appropriate name of train graph, e.g. looks for graph in log-dir
gnn.show_gnn_classification(pred_graph_prefix=train_group_name, train_view=True)
