topoml: Segmentation of MSC from Image, conversion to graph data structure
	graph sage learning
	
	mscnn_seg.py: 
	    - callable script e.g. python mscnn_seg.py
	    - use of msc_segmentation class to obtain msc or geometric msc from images
		      converted to graph data structure
		- compute geometric msc from images for graph data structure and learning
		Uses:  topoml.mscnn_segmentation import mscnn_segmentation     
		        which uses: topoml.topology.read_msc import MSC and topoml.topology.geometric_msc import GeoMSC       

	test_gnn: compute msc of passed image into graph structure or take msc_segementation class output.
		  run graph neural network

	mscnn_segmentation: compute msc segmentation to use by graph neural nets
			    encodes in graph data structure needed and saves example to file
			            compute geometric morse smale segmentation
			            
	topoml.arcneurontracer: encodes msc in graph structure of arcs nodes pixels
	                        inherets from neuron tracer which constructs the msc or geomsc using topoml.topology.utils.build_msc (or build_geometric_msc)

	topoml.neurontracer: parent class to constructing msc. uses topoml.topology.utils to construct msc
	
	graph_learning.classify_feats_graphs_v2 class graph_learners: employs various graph learning strategies
	    (CAMLP, Deep Walk, Confidence Label Prediction, Random Forests, locally and globally consistency) using msc saved as csv

    topoml.ui.ArcSelector: compute or use computed MSC to create graph
        -  msc_selector.draw_binary_segmentation  draw MSC graph as binary image
        
     topoml.topology.geometric_msc import GeoMSC 
        - read from file
        
GradIntegrator:
    extract2dridgegraph.cxx: compute geometric msc of 2d input. write to 
                             files filename"_vertices.txt" contents: vertex_id , x , y
                                   filename"_edges.txt" contents: edge_id, adjacent vertex 1 id, vertex 2 id, x_1, y_1, ... , x_i , y_i
                                   
Unsupervised_UNet: 
    compute msc of input image and use segmentation for training UNet.
