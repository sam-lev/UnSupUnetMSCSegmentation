from __future__ import print_function

import numpy as np
import random
import json
import sys
import os

import networkx as nx
from networkx.readwrite import json_graph
version_info = list(map(int, nx.__version__.split('.')))
major = version_info[0]
minor = version_info[1]
assert (major <= 1) and (minor <= 11), "networkx major version > 1.11"

WALK_LEN=10
N_WALKS=50

def load_data(prefix='', positive_arcs = [], negative_arcs = [], normalize=True, load_walks=False, train_or_test = '', scheme_required = True):
    
    print('loading json graph family')

    require_scheme = scheme_required
    label_scheme = train_or_test
    
    G_data = json.load(open(prefix + "-G.json"))
    G = json_graph.node_link_graph(G_data)
    
    if isinstance(G.nodes()[0], int):
        conversion = lambda n : int(n)
    else:
        conversion = lambda n : n
        
    negative_sample_count = 0
    positive_sample_count = 0
    broken_count = 0
    total_nodes = 0
    val_count = 0
    unlabeled_nodes = 0

    #if not train_or_test:
    #    sys.exit("must specify graph as train or test set")

    if positive_arcs and negative_arcs:
        cvt_indices = positive_arcs + negative_arcs

    for node in G.nodes():
        total_nodes+=1
        before = broken_count
        if not 'val' in G.node[node] or not 'test' in G.node[node] or not 'label' in G.node[node]:# or not bool(np.sum(G.node[node]["label"])):
            G.remove_node(node)
            broken_count += 1
            
        if before == broken_count and "label" in G.node[node]:
            if G.node[node]["label"][0] == 1 and G.node[node]['train']:
                negative_sample_count+=1
            if G.node[node]["label"][1] == 1 and G.node[node]['train']:
                positive_sample_count+=1
            if G.node[node]['val']:
                val_count+=1
 
        if before == broken_count  and train_or_test == 'test' and require_scheme and not bool(np.sum(G.node[node]["label"])):
            G.remove_node(node)
            unlabeled_nodes+=1

    if os.path.exists(prefix + "-feats.npy"):
        feats = np.load(prefix + "-feats.npy")
    
    else:
        print("No features present.. Only identity features will be used.")
        feats = None

    
    print('feature shape: ', feats.shape)
    
    id_map = json.load(open(prefix + "-id_map.json"))
    id_map = {conversion(k):int(v) for k,v in id_map.items()}
    walks = []
    class_map = json.load(open(prefix + "-class_map.json"))
    if isinstance(list(class_map.values())[0], list):
        lab_conversion = lambda n : n
    else:
        lab_conversion = lambda n : int(n)

    class_map = {conversion(k):lab_conversion(v) for k,v in class_map.items()}# if bool(np.sum(v))==require_scheme}


    print("Removed {:d} nodes that lacked proper annotations".format(broken_count))
    print("Total Nodes: ",total_nodes)
    print("Negative Samples: ", negative_sample_count)
    print("Positive Samples: ", positive_sample_count)
    print("Validation Samples: ", val_count)
    print("Unlabeled samples: ", unlabeled_nodes )
    
    ## Make sure the graph has edge train_removed annotations
    ## (some datasets might already have this..)
    print("Loaded data.. now preprocessing..")
    for edge in G.edges():
        if (G.node[edge[0]]['val'] or G.node[edge[1]]['val'] or
            G.node[edge[0]]['test'] or G.node[edge[1]]['test']):
            G[edge[0]][edge[1]]['train_removed'] = True
        else:
            G[edge[0]][edge[1]]['train_removed'] = False

    if normalize and not feats is None:
        from sklearn.preprocessing import StandardScaler
        train_ids = np.array([id_map[n] for n in G.nodes() if G.node[n]['train']])
        train_feats = feats[train_ids]
        scaler = StandardScaler()
        scaler.fit(train_feats)
        feats = scaler.transform(feats)
    
    if load_walks:
        with open('./data/random_walks/'+load_walks + "-walks.txt") as fp:
            for line in fp:
                walks.append(map(conversion, line.split()))

    return G, feats, id_map, walks, class_map, negative_sample_count, positive_sample_count

def format_data(dual=None, features=None, node_id=None, id_map=None, node_classes=None, train_or_test = '', scheme_required = True, load_walks=False, normalize=True):

    #if not train_or_test:
    #sys.exit("must specify graph as train or test set")
            
    if dual:
        G = json_graph.node_link_graph(dual)
        if isinstance(G.nodes()[0], int):
            conversion = lambda n : int(n)
        else:
            conversion = lambda n : n

    feats = features

    if node_id:
        id_map = node_id
        id_map = {conversion(k):int(v) for k,v in id_map.items()}

    class_map = node_classes
    if isinstance(list(class_map.values())[0], list):
        lab_conversion = lambda n : n
    else:
        lab_conversion = lambda n : int(n)
    class_map = {conversion(k):lab_conversion(v) for k,v in class_map.items()}


    require_scheme = scheme_required
    if isinstance(G.nodes()[0], int):
        conversion = lambda n : int(n)
    else:
        conversion = lambda n : n
        
    negative_sample_count = 0
    positive_sample_count = 0
    broken_count = 0
    total_nodes = 0
    unlabeled_nodes = 0
    val_count = 0
    for node in G.nodes():
        total_nodes+=1
        before = broken_count
        
        if not 'val' in G.node[node] or not 'test' in G.node[node] or not 'train' in G.node[node]:
            G.remove_node(node)
            broken_count += 1
        if before == broken_count and "label" in G.node[node]:
            if G.node[node]["label"][0] > 0 and G.node[node]['train']:
                negative_sample_count+=1
            if G.node[node]["label"][1] > 0 and G.node[node]['train']:
                positive_sample_count+=1
            if G.node[node]['val']:
                val_count+=1
        #if before == broken_count and train_or_test and not G.node[node]['val']:
        #G.node[node]["train"] = label_scheme == 'train'
        #G.node[node]["test"] = label_scheme == 'test'
        if before==broken_count and require_scheme and not bool(np.sum(G.node[node]["label"])):
            G.remove_node(node)
            unlabeled_nodes+=1


    walks = []
    if isinstance(list(class_map.values())[0], list):
        lab_conversion = lambda n : n
    else:
        lab_conversion = lambda n : int(n)


    print("Removed {:d} nodes that lacked proper annotations".format(broken_count))
    print("..total nodes: ",total_nodes)
    print("..positive samples: ",positive_sample_count)
    print("..negative Samples: ", negative_sample_count)
    print("..validation samples: ", val_count)
    print("..unlabeled samples: ", unlabeled_nodes )
    ## Make sure the graph has edge train_removed annotations
    ## (some datasets might already have this..)
    print("Loaded data.. now preprocessing..")
    for edge in G.edges():
        if (G.node[edge[0]]['val'] or G.node[edge[1]]['val'] or
            G.node[edge[0]]['test'] or G.node[edge[1]]['test']):
            G[edge[0]][edge[1]]['train_removed'] = True
        else:
            G[edge[0]][edge[1]]['train_removed'] = False

    if normalize and not feats is None:
        from sklearn.preprocessing import StandardScaler
        train_ids = np.array([id_map[n] for n in G.nodes() if G.node[n]['train']])
        train_feats = feats[train_ids]
        scaler = StandardScaler()
        scaler.fit(train_feats)
        feats = scaler.transform(feats)
    
    if load_walks:
        print('loading walks...')
        with open('./data/random_walks/'+load_walks + "-walks.txt") as fp:
            for line in fp:
                walks.append(map(conversion, line.split()))

    return G, feats, id_map, walks, class_map, negative_sample_count, positive_sample_count

def run_random_walks(G, nodes, num_walks=N_WALKS):
    pairs = []
    for count, node in enumerate(nodes):
        if G.degree(node) == 0:
            continue
        for i in range(num_walks):
            curr_node = node
            for j in range(WALK_LEN):
                next_node = random.choice(G.neighbors(curr_node))
                # self co-occurrences are useless
                if curr_node != node:
                    pairs.append((node,curr_node))
                curr_node = next_node
        if count % 1000 == 0:
            print("Done walks for", count, "nodes")
    return pairs

if __name__ == "__main__":
    """ Run random walks """
    #example run: python3 topoml/graphsage/utils.py ./data/json_graphs/test_ridge_arcs-G.json ./data/random_walks/full_msc_n-1_k-40
    graph_file = sys.argv[1]
    out_file = sys.argv[2]
    G_data = json.load(open(graph_file))
    G = json_graph.node_link_graph(G_data)
    nodes = [n for n in G.nodes() if 'train' in G.node[n] or 'val' in G.node[n] and (G.node[n]['train'] or G.node[n]['val'])]
    print(len(nodes))
    G = G.subgraph(nodes)
    pairs = run_random_walks(G, nodes)
    with open(out_file+'-walks.txt', "w") as fp:
        fp.write("\n".join([str(p[0]) + "\t" + str(p[1]) for p in pairs]))
