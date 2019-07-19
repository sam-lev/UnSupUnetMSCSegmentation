#!/home/sci/samlev/bin/bin/python3                                          
#SBATCH --time=21-00:00:00 # walltime, abbreviated by -t                       
#SBATCH --mem=60G                                                             
#SBATCH -o slurm-%j.out-%N # name of the stdout, using the job number (%j) and the first node (%N)                                                            
#SBATCH -e slurm-%j.err-%N # name of the stderr, using the job number (%j) and the first node (%N)                                                            
#SBATCH --gres=gpu:1

from __future__ import division
from __future__ import print_function

import os
import time
import tensorflow as tf
import numpy as np

from .models import SampleAndAggregate, SAGEInfo, Node2VecModel
from .minibatch import EdgeMinibatchIterator
from .neigh_samplers import UniformNeighborSampler
from .utils import load_data


#
# ii. Could update loss given msc behaviours with relation to 'negative' and
#       'positive' arcs-loss dependendt on aggregation of previous to current
#       both negative and positive samples.
# ii. Could add layers to allow convergence faster. Currently model is a
#       two layer multiperceptron
# i.  Could add adaptive hyperparamters e.g. loss drop given some epoch
#       number or exponential/other decay
# iv. Persistence affording subsets of graphs could be to a benefit in
#       explanding training set size or allowing training over set with
#       persistence dependence e.g. birth/death of msc added as hyperparmeter
#       during training. Start simple with high and have training set change by
#       lowering persitence iteratively.
# v.  Construct an aggregator that propogates along arcs depending on
#       persistence. Persistence Weighted Aggregator which is more likely to
#       move along high more persistent arcs. Weighted random walk based on
#       persistence.
# i.  Increase/diversify validation set 
#
# Geometric attributes:
#   geomretric no overlap each pixel part of arc is unique arc.
#   instead of each saddle with 4 incident, instead heres arc heres two nodes
#   at end.
#   extend across lines, cosine similarity of normals
#   laplacian adds dark spot on each line so tangential maxima connect in order to connect to minumum and not loose info for bounded region w/ in max 
#   send classification image
#   train on full image classify full image
#   Look into neighborhood properties based on geomretry
#
# -change weight initialization (currently xavier which is good w/ l2 regularization so would need to cater weight to feature attributes)
# - normalize features around zero instead of mean
# -add layers(?)
# -add identity element for features
# -try to overfit small set
# -play with negative bc loss meant to diverge pos from neg samples in embedding
#  - number neg samples plays large role cross entropy log of class seems to be log of number negative used. 

class unsupervised:
    def __init__(self, aggregator = 'graphsage_mean', slurm = False):
        self.aggregator =  aggregator
        self.graph = None
        self.features = None
        self.id_map = None
        self.node_classes = None
        self.slurm = slurm

    def train(self,  G=None, feats=None, id_map=None, walks=None, class_map=None, train_prefix='', load_walks=False, number_negative_samples = None, number_positive_samples=None, embedding_file_out = '', learning_rate = None, depth = 3, epochs = 200, positive_arcs = [], negative_arcs = [], weight_decay = 0.001, polarity = 6):
        slurm = self.slurm
        if slurm == False:
            os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"

        # Set random seed
        seed = 123
        np.random.seed(seed)
        tf.set_random_seed(seed)

        # Settings
        self.flags = tf.app.flags
        self.FLAGS = self.flags.FLAGS

        self.number_negative_samples = 20
        if number_negative_samples:
            self.number_negative_samples = number_negative_samples

        self.model_name = embedding_file_out

        if not learning_rate:
            learning_rate = 0.001
        self.learning_rate = learning_rate

        self.epochs = epochs

        self.depth = depth

        # define graph embedding dimensionality
        # dimension 2x used value with concat
        dim = int(474./10.)

        concat = False #mean aggregator only one to perform concat
        self.dim_feature_space = int((dim+1)/2) if concat else dim

        # Collect training data
        # train from saved file, assumes pre-labeled train/test nodes
        if train_prefix and G is None:
            print("loading graph data for gnn training")
            self.train_prefix = train_prefix
           
        
            train_data = load_data(train_prefix, load_walks=False, scheme_required = True, train_or_test='train')

            self.number_negative_samples = train_data[len(train_data)-2]
            number_positive_samples = train_data[len(train_data)-1]
            number_samples = self.number_negative_samples + number_positive_samples
            proportion_negative = int(number_samples/float(self.number_negative_samples))

        # train from passed graph, assumed pre-labeled(/processed)
        # graph with test/train nodes
        elif G is not None and feats is not None and  id_map is not None and class_map is not None and not train_prefix:
            train_prefix = 'nNeg-'+str(number_negative_samples)+'nPos-'+str(number_positive_samples)
            print("using pre-processed graph data for gnn training")
            self.number_negative_samples = number_negative_samples
            number_samples= number_negative_samples+number_positive_samples
            proportion_negative = int(number_samples/float(number_negative_samples))
            train_data = (G, feats, id_map, walks, class_map, [], [])
            
        # train from cvt sampled graph and respective in/out arcs as train
        elif positive_arcs and negative_arcs:
            train_data = load_data(positive_arcs, negative_arcs, load_walks=False, scheme_required = True, train_or_test='train')
            self.number_negative_samples =  len(negative_arcs)
            number_samples = len(positive_arcs)+len(negative_arcs)
            proportion_negative = int(number_samples/float(self.number_negative_samples))
        #keep labeled (test/train) graph for later use in testing
        self.graph = train_data[0]
        self.features = train_data[1]
        self.id_map = train_data[2]
        self.node_classes = train_data[4]
        
        tf.app.flags.DEFINE_boolean('log_device_placement', False,
                                    """Whether to log device placement.""")
        #core params..
        self.flags.DEFINE_string('model', self.aggregator, 'model names. See README for possible values.') #mean aggregator does not perform concat 
        self.flags.DEFINE_float('learning_rate', self.learning_rate, 'initial learning rate.')
        self.flags.DEFINE_integer('drop_1', 2, 'epoch to reduce learning rate first time  by a tenth')
        self.flags.DEFINE_integer('drop_2', 175, 'epoch to reduce learning rate for the second time by a tenth')
        self.flags.DEFINE_string("model_size", "big", "Can be big or small; model specific def'ns")
        self.flags.DEFINE_string('train_prefix', train_prefix, 'name of the object file that stores the training data. must be specified.')
        
        self.flags.DEFINE_string('model_name', self.model_name, 'name of the embedded graph model file is created.')

        self.flags.DEFINE_integer('depth', self.depth, 'epoch to reduce learning rate for the second time by a tenth') #I added this, journal advocates depth of 2 but loss seems to improve with more
        

        # left to default values in main experiments 
        self.flags.DEFINE_integer('epochs', self.epochs, 'number of epochs to train.')
        self.flags.DEFINE_float('dropout', 0.0, 'dropout rate (1 - keep probability).')
        self.flags.DEFINE_float('weight_decay', weight_decay, 'weight for l2 loss on embedding matrix.')
        self.flags.DEFINE_integer('max_degree', 64*3, 'maximum node degree.')
        self.flags.DEFINE_integer('samples_1', 16, 'number of samples in layer 1') #neighborhood sample size-currently set to whats used in paper, list of samples of variable hops away for convolving at each layer Length #layers +1
        self.flags.DEFINE_integer('samples_hidden', 16, 'number of samples in hidden layers')
        self.flags.DEFINE_integer('samples_2', 16,'number of users samples in layer 2') # neighborhoos sample size
        self.flags.DEFINE_integer('dim_1',32 , 'Size of output dim (final is 2x this, if using concat)')#mean aggregator does not perform concat else do, list of dimensions of the hidden representations from the input layer to the final latyer, length #Layers+1
        self.flags.DEFINE_integer('dim_hidden', 32, 'Size of output dim (final is 2x this, if using concat)')
        self.flags.DEFINE_integer('dim_2', 32, 'Size of output dim (final is 2x this, if using concat)')
        self.flags.DEFINE_boolean('random_context', False, 'Whether to use random context or direct edges')
        self.flags.DEFINE_integer('neg_sample_size', polarity, 'number of negative samples') # paper hard set to twenty rather than actual negative. defines the 'weight' on which neighboring negative nodes have on the loss function allowing a spread in the embedding space of positive and negative samples.
        self.flags.DEFINE_integer('batch_size', 8, 'minibatch size.')
        self.flags.DEFINE_integer('n2v_test_epochs', 1, 'Number of new SGD epochs for n2v.') # node to vector
        self.flags.DEFINE_integer('identity_dim', 0, 'Set to positive value to use identity embedding features of that dimension. Default 0.')

        #logging, saving, validation settings etc.
        self.flags.DEFINE_boolean('save_embeddings', True, 'whether to save embeddings for all nodes after training')
        self.flags.DEFINE_string('base_log_dir', './log-dir', 'base directory for logging and saving embeddings')
        self.flags.DEFINE_integer('validate_iter', 5000, "how often to run a validation minibatch.")
        self.flags.DEFINE_integer('validate_batch_size', 4, "how many nodes per validation sample.")
        self.flags.DEFINE_integer('gpu', 0, "which gpu to use.")
        self.flags.DEFINE_integer('print_every', 25, "How often to print training info.")
        self.flags.DEFINE_integer('max_total_steps', 1000000000, "Maximum total number of iterations")

        if slurm == False:
            os.environ["CUDA_VISIBLE_DEVICES"]=str(self.FLAGS.gpu)

            #self.GPU_MEM_FRACTION = 0.8

        #begin training
        print('Begin GNN training')
        print('')
        self._train(train_data[:-2])

    def log_dir(self):
        log_dir = self.FLAGS.base_log_dir + "/" + self.FLAGS.model_name+"-unsup-json_graphs" #+ self.FLAGS.train_prefix.split("/")[-2]
        log_dir += "/{model:s}_{model_size:s}_{lr:0.6f}/".format(
                model=self.FLAGS.model,
                model_size=self.FLAGS.model_size,
                lr=self.FLAGS.learning_rate)
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)
        return log_dir

    # Define model evaluation function
    def evaluate(self, sess, model, minibatch_iter, size=None):
        t_test = time.time()
        feed_dict_val = minibatch_iter.val_feed_dict(size)
        outs_val = sess.run([model.loss, model.ranks, model.mrr], 
                            feed_dict=feed_dict_val)
        return outs_val[0], outs_val[1], outs_val[2], (time.time() - t_test)

    def incremental_evaluate(self, sess, model, minibatch_iter, size):
        t_test = time.time()
        finished = False
        val_losses = []
        val_mrrs = []
        iter_num = 0
        while not finished:
            feed_dict_val, finished, _ = minibatch_iter.incremental_val_feed_dict(size, iter_num)
            iter_num += 1
            outs_val = sess.run([model.loss, model.ranks, model.mrr], 
                                feed_dict=feed_dict_val)
            val_losses.append(outs_val[0])
            val_mrrs.append(outs_val[2])
        return np.mean(val_losses), np.mean(val_mrrs), (time.time() - t_test)

    def save_val_embeddings(self,sess, model, minibatch_iter, size, out_dir, mod=""):
        val_embeddings = []
        finished = False
        seen = set([])
        nodes = []
        iter_num = 0
        name = "val"
        while not finished:
            feed_dict_val, finished, edges = minibatch_iter.incremental_embed_feed_dict(size, iter_num)
            iter_num += 1
            outs_val = sess.run([model.loss, model.mrr, model.outputs1], 
                                feed_dict=feed_dict_val)
            #ONLY SAVE FOR embeds1 because of planetoid
            for i, edge in enumerate(edges):
                if not edge[0] in seen:
                    val_embeddings.append(outs_val[-1][i,:])
                    nodes.append(edge[0])
                    seen.add(edge[0])
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        val_embeddings = np.vstack(val_embeddings)
        np.save(out_dir + name + mod + ".npy",  val_embeddings)
        with open(out_dir + name + mod + ".txt", "w") as fp:
            fp.write("\n".join(map(str,nodes)))

    def construct_placeholders(self):
        # Define placeholders
        placeholders = {
            'batch1' : tf.placeholder(tf.int32, shape=(None), name='batch1'),
            'batch2' : tf.placeholder(tf.int32, shape=(None), name='batch2'),
            # negative samples for all nodes in the batch
            'neg_samples': tf.placeholder(tf.int32, shape=(None,),
                name='neg_sample_size'),
            'dropout': tf.placeholder_with_default(0., shape=(), name='dropout'),
            'batch_size' : tf.placeholder(tf.int32, name='batch_size'),
        }
        return placeholders

    def _train(self, train_data, test_data=None):
        G = train_data[0]
        features = train_data[1]
        id_map = train_data[2]
        
        if not features is None:
            # pad with dummy zero vector
            features = np.vstack([features, np.zeros((features.shape[1],))])

        context_pairs = train_data[3] if self.FLAGS.random_context else None
        placeholders = self.construct_placeholders()
        minibatch = EdgeMinibatchIterator(G, 
                id_map,
                placeholders, batch_size=self.FLAGS.batch_size,
                max_degree=self.FLAGS.max_degree, 
                num_neg_samples=self.FLAGS.neg_sample_size,
                context_pairs = context_pairs)
        adj_info_ph = tf.placeholder(tf.int32, shape=minibatch.adj.shape)
        adj_info = tf.Variable(adj_info_ph, trainable=False, name="adj_info")

        if self.FLAGS.model == 'graphsage_mean':
            # Create model
            # for more layers add layers to MLP in models.py as well as
            # add SAGEInfo nodes for more layers [layer name, neighbor sampler,
            #               number neighbors sampled, out dim]
            sampler = UniformNeighborSampler(adj_info)
            layer_infos = [SAGEInfo("node", sampler, self.FLAGS.samples_1, self.FLAGS.dim_1),
                                SAGEInfo("node", sampler, self.FLAGS.samples_2, self.FLAGS.dim_2)]

            model = SampleAndAggregate(placeholders, 
                                         features,
                                         adj_info,
                                         minibatch.deg,
                                         layer_infos=layer_infos,
                                         depth = self.depth,
                                         model_size=self.FLAGS.model_size,
                                         identity_dim = self.FLAGS.identity_dim,
    
                                         logging=True)
            """callbacks=[
                ModelSaver(), # Record state graph at intervals during epochs
                InferenceRunner(dataset_train,
                                [ScalarStats('cost'), ClassificationError()]), #Compare to validation set
                ScheduledHyperParamSetter('learning_rate',
                                          [(1, 0.1), (args.drop_1, 0.01), (args.drop_2, 0.001)]) # denote current hyperparameters
            ],"""
        elif self.FLAGS.model == 'gcn':
            # Create model
            sampler = UniformNeighborSampler(adj_info)
            layer_infos = [SAGEInfo("node", sampler, self.FLAGS.samples_1, 2*self.FLAGS.dim_1)]
            for l in range(self.depth-2):
                layer = SAGEInfo("node", sampler, self.FLAGS.samples_hidden, 2*self.FLAGS.dim_hidden)
                layer_infos.append(layer)
            layer_infos.append( SAGEInfo("node", sampler, self.FLAGS.samples_2, 2*self.FLAGS.dim_2))
            layer_infos = [SAGEInfo("node", sampler, self.FLAGS.samples_1, 2*self.FLAGS.dim_1),
                                SAGEInfo("node", sampler, self.FLAGS.samples_2, 2*self.FLAGS.dim_2)]

            model = SampleAndAggregate(placeholders, 
                                         features,
                                         adj_info,
                                         minibatch.deg,
                                         layer_infos=layer_infos, 
                                         aggregator_type="gcn",
                                         model_size=self.FLAGS.model_size,
                                         identity_dim = self.FLAGS.identity_dim,
                                         concat=False,
                                         logging=True)

        elif self.FLAGS.model == 'graphsage_seq':
            sampler = UniformNeighborSampler(adj_info)
            layer_infos = [SAGEInfo("node", sampler, self.FLAGS.samples_1, self.FLAGS.dim_1),
                                SAGEInfo("node", sampler, self.FLAGS.samples_2, self.FLAGS.dim_2)]

            model = SampleAndAggregate(placeholders, 
                                         features,
                                         adj_info,
                                         minibatch.deg,
                                         layer_infos=layer_infos, 
                                         identity_dim = self.FLAGS.identity_dim,
                                         aggregator_type="seq",
                                         model_size=self.FLAGS.model_size,
                                         logging=True)

        elif self.FLAGS.model == 'graphsage_maxpool':
            sampler = UniformNeighborSampler(adj_info)
            layer_infos = [SAGEInfo("node", sampler, self.FLAGS.samples_1, self.FLAGS.dim_1)]
            for l in range(self.depth-2):
                layer = SAGEInfo("node", sampler, self.FLAGS.samples_hidden, 2*self.FLAGS.dim_hidden)
                layer_infos.append(layer)
            layer_infos.append(SAGEInfo("node", sampler, self.FLAGS.samples_2, self.FLAGS.dim_2))
            #layer_infos = [SAGEInfo("node", sampler, self.FLAGS.samples_1, self.FLAGS.dim_1),
            #                    SAGEInfo("node", sampler, self.FLAGS.samples_2, self.FLAGS.dim_2)]
            model = SampleAndAggregate(placeholders, 
                                        features,
                                        adj_info,
                                        minibatch.deg,
                                         layer_infos=layer_infos, 
                                         aggregator_type="maxpool",
                                         model_size=self.FLAGS.model_size,
                                         identity_dim = self.FLAGS.identity_dim,
                                         logging=True)
            """callbacks = [
                ScheduledHyperParamSetter('learning_rate',
                                          [(1, 0.001), (self.FLAGS.drop_1, 0.0001), (self.FLAGS.drop_2, 0.00001)])
                ]"""
        elif self.FLAGS.model == 'graphsage_meanpool':
            sampler = UniformNeighborSampler(adj_info)
            layer_infos = [SAGEInfo("node", sampler, self.FLAGS.samples_1, self.FLAGS.dim_1)]
            for l in range(self.depth-2):
                layer = SAGEInfo("node", sampler, self.FLAGS.samples_hidden, self.FLAGS.dim_hidden)
                layer_infos.append(layer)
            layer_infos.append(SAGEInfo("node", sampler, self.FLAGS.samples_2, self.FLAGS.dim_2))
            #layer_infos = [SAGEInfo("node", sampler, self.FLAGS.samples_1, self.FLAGS.dim_1),
            #                    SAGEInfo("node", sampler, self.FLAGS.samples_2, self.FLAGS.dim_2)]

            model = SampleAndAggregate(placeholders, 
                                        features,
                                        adj_info,
                                        minibatch.deg,
                                         layer_infos=layer_infos, 
                                         aggregator_type="meanpool",
                                         model_size=self.FLAGS.model_size,
                                         identity_dim = self.FLAGS.identity_dim,
                                         logging=True)

        elif self.FLAGS.model == 'n2v':
            model = Node2VecModel(placeholders, features.shape[0],
                                           minibatch.deg,
                                           #2x because graphsage uses concat
                                           nodevec_dim=2*self.FLAGS.dim_1,
                                           lr=self.FLAGS.learning_rate)
        else:
            raise Exception('Error: model name unrecognized.')

        config = tf.ConfigProto(log_device_placement=self.FLAGS.log_device_placement)
        config.gpu_options.allow_growth = True
        #config.gpu_options.per_process_gpu_memory_fraction = GPU_MEM_FRACTION
        config.allow_soft_placement = True

        # Initialize session
        sess = tf.Session(config=config)
        merged = tf.summary.merge_all()
        summary_writer = tf.summary.FileWriter(self.log_dir(), sess.graph)

        # Init variables
        sess.run(tf.global_variables_initializer(), feed_dict={adj_info_ph: minibatch.adj})

        # Train model

        train_shadow_mrr = None
        shadow_mrr = None

        total_steps = 0
        avg_time = 0.0
        epoch_val_costs = []

        train_adj_info = tf.assign(adj_info, minibatch.adj)
        val_adj_info = tf.assign(adj_info, minibatch.test_adj)
        for epoch in range(self.FLAGS.epochs): 
            minibatch.shuffle() 

            iter = 0
            print('...')
            print('Epoch: %04d' % (epoch + 1))
            print('...')

            epoch_val_costs.append(0)
            while not minibatch.end():
                # Construct feed dictionary
                feed_dict = minibatch.next_minibatch_feed_dict()
                feed_dict.update({placeholders['dropout']: self.FLAGS.dropout})

                t = time.time()
                # Training step
                outs = sess.run([merged, model.opt_op, model.loss, model.ranks, model.aff_all, 
                        model.mrr, model.outputs1], feed_dict=feed_dict)
                train_cost = outs[2]
                train_mrr = outs[5]
                if train_shadow_mrr is None:
                    train_shadow_mrr = train_mrr#
                else:
                    train_shadow_mrr -= (1-0.99) * (train_shadow_mrr - train_mrr)

                if iter % self.FLAGS.validate_iter == 0:
                    # Validation
                    sess.run(val_adj_info.op)
                    val_cost, ranks, val_mrr, duration  = self.evaluate(sess, model, minibatch, size=self.FLAGS.validate_batch_size)
                    sess.run(train_adj_info.op)
                    epoch_val_costs[-1] += val_cost
                if shadow_mrr is None:
                    shadow_mrr = val_mrr
                else:
                    shadow_mrr -= (1-0.99) * (shadow_mrr - val_mrr)

                if total_steps % self.FLAGS.print_every == 0:
                    summary_writer.add_summary(outs[0], total_steps)

                # Print results
                avg_time = (avg_time * total_steps + time.time() - t) / (total_steps + 1)

                if total_steps % self.FLAGS.print_every == 0:
                    print("Iter:", '%04d' % iter, 
                          "train_loss=", "{:.5f}".format(train_cost),
                          "train_mrr=", "{:.5f}".format(train_mrr), 
                          "train_mrr_ema=", "{:.5f}".format(train_shadow_mrr), # exponential moving average
                          "val_loss=", "{:.5f}".format(val_cost),
                          "val_mrr=", "{:.5f}".format(val_mrr), 
                          "val_mrr_ema=", "{:.5f}".format(shadow_mrr), # exponential moving average
                          "time=", "{:.5f}".format(avg_time))

                iter += 1
                total_steps += 1

                if total_steps > self.FLAGS.max_total_steps:
                    break

            if total_steps > self.FLAGS.max_total_steps:
                    break

        print("Optimization Finished.")
        if self.FLAGS.save_embeddings:
            sess.run(val_adj_info.op)

            self.save_val_embeddings(sess, model, minibatch, self.FLAGS.validate_batch_size, self.log_dir())

            if self.FLAGS.model == "n2v":
                # stopping the gradient for the already trained nodes
                train_ids = tf.constant([[id_map[n]] for n in G.nodes_iter() if not G.node[n]['val'] and not G.node[n]['test']],
                        dtype=tf.int32)
                test_ids = tf.constant([[id_map[n]] for n in G.nodes_iter() if G.node[n]['val'] or G.node[n]['test']], 
                        dtype=tf.int32)
                update_nodes = tf.nn.embedding_lookup(model.context_embeds, tf.squeeze(test_ids))
                no_update_nodes = tf.nn.embedding_lookup(model.context_embeds,tf.squeeze(train_ids))
                update_nodes = tf.scatter_nd(test_ids, update_nodes, tf.shape(model.context_embeds))
                no_update_nodes = tf.stop_gradient(tf.scatter_nd(train_ids, no_update_nodes, tf.shape(model.context_embeds)))
                model.context_embeds = update_nodes + no_update_nodes
                sess.run(model.context_embeds)

                # run random walks
                from .utils import run_random_walks
                nodes = [n for n in G.nodes_iter() if G.node[n]["val"] or G.node[n]["test"]]
                start_time = time.time()
                pairs = run_random_walks(G, nodes, num_walks=50)
                walk_time = time.time() - start_time

                test_minibatch = EdgeMinibatchIterator(G, 
                    id_map,
                    placeholders, batch_size=self.FLAGS.batch_size,
                    max_degree=self.FLAGS.max_degree, 
                    num_neg_samples=self.FLAGS.neg_sample_size,
                    context_pairs = pairs,
                    n2v_retrain=True,
                    fixed_n2v=True)

                start_time = time.time()
                print("Doing test training for n2v.")
                test_steps = 0
                for epoch in range(self.FLAGS.n2v_test_epochs):
                    test_minibatch.shuffle()
                    while not test_minibatch.end():
                        feed_dict = test_minibatch.next_minibatch_feed_dict()
                        feed_dict.update({placeholders['dropout']: self.FLAGS.dropout})
                        outs = sess.run([model.opt_op, model.loss, model.ranks, model.aff_all, 
                            model.mrr, model.outputs1], feed_dict=feed_dict)
                        if test_steps % self.FLAGS.print_every == 0:
                            print("Iter:", '%04d' % test_steps, 
                                  "train_loss=", "{:.5f}".format(outs[1]),
                                  "train_mrr=", "{:.5f}".format(outs[-2]))
                        test_steps += 1
                train_time = time.time() - start_time
                save_val_embeddings(sess, model, minibatch, self.FLAGS.validate_batch_size, log_dir(), mod="-test")
                print("Total time: ", train_time+walk_time)
                print("Walk time: ", walk_time)
                print("Train time: ", train_time)


if __name__ == '__main__':
    tf.app.run()
