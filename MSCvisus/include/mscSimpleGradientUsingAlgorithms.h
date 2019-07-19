#ifndef MSC_SIMPLE_GRADIENT_USING_ALGORITHMS
#define MSC_SIMPLE_GRADIENT_USING_ALGORITHMS

#include "mscIndexTypes.h"
#include "mscBasicDataHandler.h"
#include "mscBasicMeshHandler.h"
#include "mscBasicMeshFunction.h"
#include "mscBasicGradientField.h"
#include "mscArrayFactory.h"

#include <vector>
#include <queue>
#include <map>
#include <set>




using namespace std;

template<typename dtype>
class mscSimpleGradientUsingAlgorithms {
protected:

	mscBasicMeshFunction<dtype>* my_mesh_function;
	mscBasicMeshHandler* my_mesh_handler;
	mscBasicGradientField* my_grad_field;
	mscArrayFactory* my_array_factory;

public:

	mscSimpleGradientUsingAlgorithms(
		mscBasicMeshFunction<dtype>* mesh_function,
		mscBasicMeshHandler* mesh_handler,
		mscBasicGradientField* grad_field,
		mscArrayFactory* array_factory) : 
	my_mesh_function(mesh_function),
		my_mesh_handler(mesh_handler), 
		my_grad_field(grad_field),
		my_array_factory(array_factory) {
	}

	// trace "down" in gradient and fill in the result vector
	// with all cells that are found

	virtual void count_critical_points(int dim) {
		int* counts = new int[dim];
		for (int i =0; i < dim; i++) counts[i] = 0;

		for (CELL_INDEX_TYPE i = 0; i < my_mesh_handler->num_cells(); i++) {
			if (my_grad_field->get_critical(i))
				counts[my_mesh_handler->dimension(i)]++;
		}

		for (int i =0; i < dim; i++) 
			printf("index-%d=%d\n", i, counts[i]);
	}

	virtual void trace_down_cells(const CELL_INDEX_TYPE& cellid, 
		vector<CELL_INDEX_TYPE>& result) {

			queue<CELL_INDEX_TYPE> cell_queue;
			cell_queue.push(cellid);

			result.clear();
			set<CELL_INDEX_TYPE> cell_visited;

			while (! cell_queue.empty()) {
				CELL_INDEX_TYPE current = cell_queue.front();
				cell_queue.pop();

				cell_visited.insert(current);
				result.push_back(current);

				cellIterator it;
				iteratorOperator& facets = my_mesh_handler->facets(current, it);
				for (facets.begin(it); facets.valid(it); facets.advance(it)) {
					CELL_INDEX_TYPE temp_id = facets.value(it);

					if (my_grad_field->get_critical(temp_id) &&
						cell_visited.count(temp_id) == 0) {
							result.push_back(temp_id);
							cell_visited.insert(temp_id);
					} else if (cell_visited.count(temp_id) == 0) {
						CELL_INDEX_TYPE pair = my_grad_field->get_pair(temp_id);
						result.push_back(temp_id);
						result.push_back(pair);
						cell_visited.insert(temp_id);
						cell_visited.insert(pair);
						cell_queue.push(pair);
					}		
				}
			}

	}

		virtual void trace_up_cells(const CELL_INDEX_TYPE& cellid, 
		vector<CELL_INDEX_TYPE>& result) {

			queue<CELL_INDEX_TYPE> cell_queue;
			cell_queue.push(cellid);

			result.clear();
			set<CELL_INDEX_TYPE> cell_visited;

			while (! cell_queue.empty()) {
				CELL_INDEX_TYPE current = cell_queue.front();
				cell_queue.pop();

				cell_visited.insert(current);
				result.push_back(current);

				cellIterator it;
				iteratorOperator& cofacets = my_mesh_handler->cofacets(current, it);
				for (cofacets.begin(it); cofacets.valid(it); cofacets.advance(it)) {
					CELL_INDEX_TYPE temp_id = cofacets.value(it);

					if (my_grad_field->get_critical(temp_id) &&
						cell_visited.count(temp_id) == 0) {
							result.push_back(temp_id);
							cell_visited.insert(temp_id);
					} else if (cell_visited.count(temp_id) == 0) {
						CELL_INDEX_TYPE pair = my_grad_field->get_pair(temp_id);
						result.push_back(temp_id);
						result.push_back(pair);
						cell_visited.insert(temp_id);
						cell_visited.insert(pair);
						cell_queue.push(pair);
					}		
				}
			}

	}


	virtual void trace_down_cells_restricted(const CELL_INDEX_TYPE& cellid, 
		vector<CELL_INDEX_TYPE>& result) {

			queue<CELL_INDEX_TYPE> cell_queue;
			cell_queue.push(cellid);

			DIM_TYPE temp_dim = my_grad_field->get_dim_asc_man(cellid) + 1;
			result.clear();
			set<CELL_INDEX_TYPE> cell_visited;

			while (! cell_queue.empty()) {
				CELL_INDEX_TYPE current = cell_queue.front();
				cell_queue.pop();

				cell_visited.insert(current);
				result.push_back(current);

				cellIterator it;
				iteratorOperator& facets = my_mesh_handler->facets(current, it);
				for (facets.begin(it); facets.valid(it); facets.advance(it)) {
					CELL_INDEX_TYPE temp_id = facets.value(it);

					if (my_grad_field->get_critical(temp_id) &&
						cell_visited.count(temp_id) == 0) {
							result.push_back(temp_id);
							cell_visited.insert(temp_id);
					} else if (cell_visited.count(temp_id) == 0 &&
						my_grad_field->get_dim_asc_man(temp_id) == temp_dim) {
							CELL_INDEX_TYPE pair = my_grad_field->get_pair(temp_id);
							result.push_back(temp_id);
							result.push_back(pair);
							cell_visited.insert(temp_id);
							cell_visited.insert(pair);
							cell_queue.push(pair);
					}		
				}
			}

	}
	virtual void trace_down_cells_restricted_counting(const CELL_INDEX_TYPE& cellid, 
		vector<CELL_INDEX_TYPE>& result, vector<int>& counts) {

			queue<CELL_INDEX_TYPE> cell_queue;
			cell_queue.push(cellid);

			DIM_TYPE temp_dim = my_grad_field->get_dim_asc_man(cellid) + 1;
			result.clear();
			counts.clear();
			set<CELL_INDEX_TYPE> cell_visited;

			// build the graph
			map<CELL_INDEX_TYPE, set<CELL_INDEX_TYPE> > node_graph;
			map<CELL_INDEX_TYPE, int > visit_counts;

			while (! cell_queue.empty()) {
				CELL_INDEX_TYPE current = cell_queue.front();
				cell_queue.pop();

				set<CELL_INDEX_TYPE> neighbors;

				cell_visited.insert(current);

				cellIterator it;
				iteratorOperator& facets = my_mesh_handler->facets(current, it);
				for (facets.begin(it); facets.valid(it); facets.advance(it)) {
					CELL_INDEX_TYPE temp_id = facets.value(it);

					if (my_grad_field->get_critical(temp_id)) {
						neighbors.insert(temp_id);
						if (visit_counts.count(temp_id) == 0) {
							visit_counts[temp_id] = 1;
						} else {
							visit_counts[temp_id]++;
						}
						
						cell_visited.insert(temp_id);
					} else if (my_grad_field->get_dim_asc_man(temp_id) == temp_dim) {
						CELL_INDEX_TYPE pair = my_grad_field->get_pair(temp_id);
						if (current == pair) continue;

						neighbors.insert(pair);
						if (visit_counts.count(pair) == 0) {
							visit_counts[pair] = 1;
						} else {
							visit_counts[pair]++;
						}
						if (cell_visited.count(pair) == 0) {
							cell_queue.push(pair);
						}
						cell_visited.insert(temp_id);
						cell_visited.insert(pair);

					}		
				}
				node_graph[current].insert(neighbors.begin(), neighbors.end());
			}
			//print graph
			printf("\ngraph of %d:\n", cellid);
			for(map<CELL_INDEX_TYPE, set<CELL_INDEX_TYPE> >::iterator mit = node_graph.begin();
				mit != node_graph.end(); mit++) {
					CELL_INDEX_TYPE tempid = (*mit).first;
					printf(" n=%d\n", tempid);
					for(set<CELL_INDEX_TYPE>::iterator sit = (*mit).second.begin();
						sit != (*mit).second.end(); sit++)
						printf("  -->%d\n", *sit);
			}
			// traverse graph from root
			cell_queue.push(cellid);
			while (! cell_queue.empty()) {
				CELL_INDEX_TYPE current = cell_queue.front();
				cell_queue.pop();
				result.push_back(current);
				counts.push_back(0);

				for (set<CELL_INDEX_TYPE>::iterator it = node_graph[current].begin();
					it != node_graph[current].end(); it++) {
						CELL_INDEX_TYPE tempid = *it;
						visit_counts[tempid]--;
						if (visit_counts[tempid] == 0) {
							cell_queue.push(tempid);
						}
				}
			}

			// the base case, 1 path from cell to itself
			visit_counts[cellid] = 1;
			for (int i = 0; i < result.size(); i++) {
				CELL_INDEX_TYPE current = result[i];
				int temp_count = visit_counts[current];
				counts[i] = temp_count;
				for (set<CELL_INDEX_TYPE>::iterator it = node_graph[current].begin();
					it != node_graph[current].end(); it++) {
						CELL_INDEX_TYPE tempid = *it;
						visit_counts[tempid] += temp_count;
				}
			}
	}


	void rec_man_trace_up(CELL_INDEX_TYPE& cellid, set<CELL_INDEX_TYPE>& res) {
		res.insert(cellid);
		CELL_INDEX_TYPE current = cellid;
		DIM_TYPE cdim = this->my_mesh_handler->dimension(cellid);
		cellIterator it;
		iteratorOperator& cofacets = my_mesh_handler->cofacets(current, it);
		for (cofacets.begin(it); cofacets.valid(it); cofacets.advance(it)) {
			CELL_INDEX_TYPE temp_id = cofacets.value(it);
			if (this->my_grad_field->get_critical(temp_id) || ! my_grad_field->get_assigned(temp_id)) continue; 
			
			CELL_INDEX_TYPE temp_pair = my_grad_field->get_pair(temp_id);

			if (temp_pair == cellid) continue;

			if (my_mesh_handler->dimension(temp_pair) != cdim) continue;
			
			rec_man_trace_up(temp_pair, res);
		}
	}



};

#endif