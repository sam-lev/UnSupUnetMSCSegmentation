#ifndef MSC_REGULAR_2D_GRID_IMPLICIT_MESH_HANDLER
#define MSC_REGULAR_2D_GRID_IMPLICIT_MESH_HANDLER

#include "mscIndexTypes.h"
#include "mscBasicIterator.h"
#include "mscRegularGridImplicitMeshHandler.h"
#include "mscBasicIterator.h"

#include <map>

using namespace std;


class mscRegularGrid2DImplicitMeshHandler : public mscRegularGridImplicitMeshHandler {

protected:
	CELL_INDEX_TYPE my_num_cells;
	CELL_INDEX_TYPE X, Y; // the number of vertices in X, Y, Z
	CELL_INDEX_TYPE dX, dY; // the number of cells in X, Y, Z
	CELL_INDEX_TYPE extent_list[4];

	CELL_INDEX_TYPE num_d_cells[3];
	CELL_INDEX_TYPE offset_list[4];
	//CELL_INDEX_TYPE vertex_offsets[4][5];
	CELL_INDEX_TYPE facet_position_list[8][7] ;
	CELL_INDEX_TYPE cofacet_position_list[8][7] ;
	DIM_TYPE facet_dim_list[8][7];
	DIM_TYPE cofacet_dim_list[8][7];


	CELL_INDEX_TYPE& sizes(DIM_TYPE val) {
		if (val == 0) return dX;
		return dY;
	}

	map<CELL_INDEX_TYPE, unsigned char> offset_2_loc;
	

	//void set_vertex_offsets() {
	//	
	//	vertex_offsets[0][0] = 1; vertex_offsets[0][1] = 0; 
	//	vertex_offsets[1][0] = 2; vertex_offsets[1][1] = 1; vertex_offsets[1][2] = -1;
	//	vertex_offsets[2][0] = 2; vertex_offsets[2][1] = X; vertex_offsets[2][2] = -X;
	//	vertex_offsets[3][0] = 4; vertex_offsets[3][1] = X+1; vertex_offsets[3][2] = X-1; 
	//	vertex_offsets[3][3] = -X+1; vertex_offsets[3][4] = -X-1;

	//}

	void initialize_values() {
		// num cells
		dX = 2*X-1;
		dY = 2*Y-1;
		my_num_cells = dX * dY;

		//numbers of cells
		num_d_cells[0] = X*Y;
		CELL_INDEX_TYPE xdir = (X-1)*Y;
		CELL_INDEX_TYPE ydir = X*(Y-1);
		
		num_d_cells[1] =  xdir + ydir;
		xdir = X*(Y-1);
		ydir = (X-1)*Y;
		num_d_cells[2] = xdir + ydir;


		//offset_list
		offset_list[0] = (1) + (0) * dX ;
		offset_list[1] = (-1) + (0) * dX;
		offset_list[2] = (0) + (1) * dX;
		offset_list[3] = (0) + (-1) * dX;

		extent_list[0]= dX-1;
		extent_list[1]= 0;
		extent_list[2]= dY-1;
		extent_list[3]= 0;

		//set_vertex_offsets();

		//face every other since start at point mod 2  is 1 for face
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
					int id = i + 2*j;
					DIM_TYPE dim = i + j;
					
					int place = 1;
					if (i%2 == 1) {
						facet_dim_list[id][place] = 0;
						facet_position_list[id][place++] = offset_list[0];
						facet_dim_list[id][place] = 0;
						facet_position_list[id][place++] = offset_list[1];

					}
					if (j%2 == 1) {
						facet_dim_list[id][place] = 1;
						facet_position_list[id][place++] = offset_list[2];
						facet_dim_list[id][place] = 1;
						facet_position_list[id][place++] = offset_list[3];
					}

					facet_position_list[id][0] = place - 1;
			
			}
		}
		//coface mod 2 is zero, start at point 0 then face then coface ect..
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
					int id = i + 2*j ;
					DIM_TYPE dim = i + j ;
					
					int place = 1;
					if (i%2 == 0) {
						cofacet_dim_list[id][place] = 0;
						cofacet_position_list[id][place++] = offset_list[0];
						cofacet_dim_list[id][place] = 0;
						cofacet_position_list[id][place++] = offset_list[1];
					}
					if (j%2 == 0) {
						cofacet_dim_list[id][place] = 1;
						cofacet_position_list[id][place++] = offset_list[2];
						cofacet_dim_list[id][place] = 1;
						cofacet_position_list[id][place++] = offset_list[3];
					}

					cofacet_position_list[id][0] = place - 1;
				}
			
		}

		for (int i = 0; i < 4; i++) {
			offset_2_loc[offset_list[i]] = i;
			//printf("TEST: %d -> %d -> %d\n", i, offset_list[i], offset_2_loc[offset_list[i]]);
		}
	}

	
	class allCellIteratorOperator : public iteratorOperator {

	public: 
		void begin(cellIterator& it) {
			it.my_location = 0;
			it.my_total_size = it.my_mesh_handler->num_cells();
		}			
	};

	allCellIteratorOperator my_all_cells_iterator_operator;

	class dCellIteratorOperator : public allCellIteratorOperator {

	public: 
		void begin(cellIterator& it) {
			it.my_location = 0;
			it.my_total_size = it.my_mesh_handler->num_cells();
			if (it.my_mesh_handler->dimension(it.my_location) != it.my_dim)
				advance(it);
		}	
		void advance(cellIterator& it) {
			it.my_location++;
			while (it.my_location < it.my_total_size &&
				it.my_mesh_handler->dimension(it.my_location) != it.my_dim) it.my_location++;
		}			
	};
	dCellIteratorOperator my_d_cells_iterator_operator;

	class dFacetIterator : public iteratorOperator {
	public: 
		virtual bool valid(cellIterator& it) {
			return it.my_location <= it.my_total_size;
		}

		virtual void begin(cellIterator& it) {
			
			it.my_location = 1;
			((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->cellid_2_coords(it.my_base, it.my_coords);
			it.my_dim = it.my_coords[0]%2 + 2 * (it.my_coords[1]%2) ;
			it.my_total_size = 
				((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->facet_position_list[it.my_dim][0];
			//neighbor 
			it.my_neighbor = it.my_base + 
				((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->facet_position_list[it.my_dim][it.my_location];	
		}

		void advance(cellIterator& it) {
			it.my_location++;
			it.my_neighbor = it.my_base + 
				((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->facet_position_list[it.my_dim][it.my_location];	
		}	

		virtual CELL_INDEX_TYPE value(cellIterator& it) {
			return it.my_neighbor;
		}
		
		virtual CELL_INDEX_TYPE index(cellIterator& it) {
			CELL_INDEX_TYPE temp = 
				((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->facet_position_list[it.my_dim][it.my_location];
			return ((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->offset_2_loc[temp];
		}

	};

	class dBoundaryFacetIterator : public dFacetIterator {
	public: 

		virtual void begin(cellIterator& it) {
			
			it.my_location = 0;
			((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->cellid_2_coords(it.my_base, it.my_coords);
			it.my_dim = it.my_coords[0]%2 + 2 * (it.my_coords[1]%2) ;
			it.my_total_size = 
				((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->facet_position_list[it.my_dim][0];
			advance(it);
		}

		void advance(cellIterator& it) {
			it.my_location++;
			

			// if invalid, return
			if (it.my_location > it.my_total_size) return;

			// get the coords of this point
			mscRegularGrid2DImplicitMeshHandler* temp_mesh_handler = 
				((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler);

			while (it.my_location <= it.my_total_size) {
				DIM_TYPE temp_dim = temp_mesh_handler->facet_dim_list[it.my_dim][it.my_location]; 
				if (it.my_location % 2 == 0 &&
					it.my_coords[temp_dim] != 0) {
						it.my_neighbor = it.my_base + 
							temp_mesh_handler->facet_position_list[it.my_dim][it.my_location];
						return;
				}
				if (it.my_location % 2 == 1 &&
					it.my_coords[temp_dim] != temp_mesh_handler->sizes(temp_dim) - 1) {
						it.my_neighbor = it.my_base + 
							temp_mesh_handler->facet_position_list[it.my_dim][it.my_location];
						return;

				}
				it.my_location++;
			}
		}

	};


		class dCoFacetIterator : public iteratorOperator {
	public: 
		virtual bool valid(cellIterator& it) {
			return it.my_location <= it.my_total_size;
		}

		virtual void begin(cellIterator& it) {
			
			it.my_location = 1;
			((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->cellid_2_coords(it.my_base, it.my_coords);
			it.my_dim = it.my_coords[0]%2 + 2 * (it.my_coords[1]%2);
			it.my_total_size = 
				((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->cofacet_position_list[it.my_dim][0];
			//neighbor 
			it.my_neighbor = it.my_base + 
				((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->cofacet_position_list[it.my_dim][it.my_location];	
		}

		void advance(cellIterator& it) {
			it.my_location++;
			it.my_neighbor = it.my_base + 
				((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->cofacet_position_list[it.my_dim][it.my_location];	
		}	

		virtual CELL_INDEX_TYPE value(cellIterator& it) {
			return it.my_neighbor;
		}
		
		virtual CELL_INDEX_TYPE index(cellIterator& it) {
			CELL_INDEX_TYPE temp = 
				((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->cofacet_position_list[it.my_dim][it.my_location];
			return ((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->offset_2_loc[temp];
		}

	};

	class dBoundaryCoFacetIterator : public dCoFacetIterator {
	public: 
		virtual void begin(cellIterator& it) {
			
			it.my_location = 0;
			((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->cellid_2_coords(it.my_base, it.my_coords);
			it.my_dim = it.my_coords[0]%2 + 2 * (it.my_coords[1]%2) ;
			it.my_total_size = 
				((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->cofacet_position_list[it.my_dim][0];
			advance(it);
		}

		void advance(cellIterator& it) {
			it.my_location++;
			
			// if invalid, return
			if (it.my_location > it.my_total_size) return;

			// get the coords of this point
			mscRegularGrid2DImplicitMeshHandler* temp_mesh_handler = 
				((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler);

			while (it.my_location <= it.my_total_size) {
				DIM_TYPE temp_dim = temp_mesh_handler->cofacet_dim_list[it.my_dim][it.my_location];
				//if not at at start or on coface return
				if (it.my_location % 2 == 0 &&
					it.my_coords[temp_dim] != 0) {
						it.my_neighbor = it.my_base + 
							temp_mesh_handler->cofacet_position_list[it.my_dim][it.my_location];
						return;
				}
				//if on face and not at end return
				if (it.my_location % 2 == 1 &&
					it.my_coords[temp_dim] != temp_mesh_handler->sizes(temp_dim) - 1) {
						it.my_neighbor = it.my_base + 
							temp_mesh_handler->cofacet_position_list[it.my_dim][it.my_location];
						return;

				}
				//since not at start or end and on face iterate to next
				it.my_location++;
			}
		}

	};
		class dNeighborIterator : public iteratorOperator {
	public: 
		virtual bool valid(cellIterator& it) {
			return it.my_location < 4;
		}

		virtual void begin(cellIterator& it) {
			
			it.my_location = 0;
			it.my_neighbor = it.my_base + 
				((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->offset_list[it.my_location] * 2;	
		}

		void advance(cellIterator& it) {
			it.my_location++;
			it.my_neighbor = it.my_base + 
				((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->offset_list[it.my_location] * 2;		
		}	

		virtual CELL_INDEX_TYPE value(cellIterator& it) {
			return it.my_neighbor;
		}
		
		virtual CELL_INDEX_TYPE index(cellIterator& it) {
			return it.my_location;
		}

	};


	class dBoundaryNeighborIterator : public iteratorOperator {

	public: 
		virtual bool valid(cellIterator& it) {
			//printf("%d asdfasdf %d\n", it.my_location, it.my_base);
			return it.my_location < 4;
		}
		virtual void begin(cellIterator& it) {
			
			it.my_location = 0;
			((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->cellid_2_coords(it.my_base, it.my_coords);		
			if (it.my_coords[0] + 2 > ((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->dX - 1) it.my_location = 1;
			it.my_neighbor = it.my_base + 
				((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->offset_list[it.my_location] * 2;	

		}

		virtual void advance(cellIterator& it) {
			it.my_location++;
			

			// if invalid, return
			if (it.my_location > 6) return;

			// get the coords of this point
			mscRegularGrid2DImplicitMeshHandler* temp_mesh_handler = 
				((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler);
			while (it.my_location < 6) {
				CELL_INDEX_TYPE tVal = it.my_coords[it.my_location / 2];
				if (it.my_location % 2 == 0) {
					// even location so upper extent
					if ( tVal + 2 > temp_mesh_handler->extent_list[it.my_location]) {
						it.my_location++;
					} else {
						it.my_neighbor = it.my_base + 
							temp_mesh_handler->offset_list[it.my_location] * 2;	
						return;
					}
				} else {
					if ( tVal - 2 < 0) {
						it.my_location++;
					} else {
						it.my_neighbor = it.my_base + 
							temp_mesh_handler->offset_list[it.my_location] * 2;	
						return;
					}
				}


			}
		}
		virtual CELL_INDEX_TYPE value(cellIterator& it) {
			return it.my_neighbor;
		}
		
		virtual CELL_INDEX_TYPE index(cellIterator& it) {
			return it.my_location;
		}
	};

	//class dVertexIterator : public iteratorOperator {
	//public: 
	//	virtual bool valid(cellIterator& it) {
	//		return it.my_location <= it.my_total_size;
	//	}

	//	virtual void begin(cellIterator& it) {
	//		
	//		it.my_location = 1;
	//		it.my_neighbor = it.my_base + 
	//			((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->vertex_offsets[it.my_dim][it.my_location];	
	//	}
		class dNeighborCellIterator : public iteratorOperator {
	public: 
		virtual bool valid(cellIterator& it) {
			return it.my_location < 4;
		}

		virtual void begin(cellIterator& it) {
			
			it.my_location = 0;
			it.my_neighbor = it.my_base + 
				((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->offset_list[it.my_location] ;	
		}

		void advance(cellIterator& it) {
			it.my_location++;
			it.my_neighbor = it.my_base + 
				((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->offset_list[it.my_location] ;		
		}	

		virtual CELL_INDEX_TYPE value(cellIterator& it) {
			return it.my_neighbor;
		}
		
		virtual CELL_INDEX_TYPE index(cellIterator& it) {
			return it.my_location;
		}

	};


	class dBoundaryNeighborCellIterator : public iteratorOperator {

	public: 
		virtual bool valid(cellIterator& it) {
			//printf("%d asdfasdf %d\n", it.my_location, it.my_base);
			return it.my_location < 4;
		}
		virtual void begin(cellIterator& it) {
			
			it.my_location = 0;
			((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->cellid_2_coords(it.my_base, it.my_coords);		
			if (it.my_coords[0] + 1 > ((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->dX - 1) it.my_location = 1;
			it.my_neighbor = it.my_base + 
				((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->offset_list[it.my_location] ;	

		}

		virtual void advance(cellIterator& it) {
			it.my_location++;
			

			// if invalid, return
			if (it.my_location > 4) return;

			// get the coords of this point
			mscRegularGrid2DImplicitMeshHandler* temp_mesh_handler = 
				((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler);
			while (it.my_location < 4) {
				CELL_INDEX_TYPE tVal = it.my_coords[it.my_location / 2];
				if (it.my_location % 2 == 0) {
					// even location so upper extent
					if ( tVal + 1 > temp_mesh_handler->extent_list[it.my_location]) {
						it.my_location++;
					} else {
						it.my_neighbor = it.my_base + 
							temp_mesh_handler->offset_list[it.my_location] ;	
						return;
					}
				} else {
					if ( tVal == 0) {
						it.my_location++;
					} else {
						it.my_neighbor = it.my_base + 
							temp_mesh_handler->offset_list[it.my_location] ;	
						return;
					}
				}


			}
		}
		virtual CELL_INDEX_TYPE value(cellIterator& it) {
			return it.my_neighbor;
		}
		
		virtual CELL_INDEX_TYPE index(cellIterator& it) {
			return it.my_location;
		}
	};
	//	void advance(cellIterator& it) {
	//		it.my_location++;
	//		it.my_neighbor = it.my_base + 
	//			((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->vertex_offsets[it.my_dim][it.my_location];	
	//	}	

	//	virtual CELL_INDEX_TYPE value(cellIterator& it) {
	//		return it.my_neighbor;
	//	}
	//	
	//	virtual CELL_INDEX_TYPE index(cellIterator& it) {
	//		return it.my_location;
	//	}

	//};
	//class dCoDCellIterator : public iteratorOperator {
	//public: 
	//	virtual bool valid(cellIterator& it) {
	//		return it.my_location < 4;
	//	}



	//	virtual void begin(cellIterator& it) {
	//		
	//		it.my_location = 0;
	//		
	//		i


	//		it.my_neighbor = it.my_base + 
	//			((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->vertex_offsets[it.my_dim][it.my_location];	
	//	}

	//	void advance(cellIterator& it) {
	//		it.my_location++;
	//		it.my_neighbor = it.my_base + 
	//			((mscRegularGrid2DImplicitMeshHandler*) it.my_mesh_handler)->vertex_offsets[it.my_dim][it.my_location];	
	//	}	

	//	virtual CELL_INDEX_TYPE value(cellIterator& it) {
	//		return it.my_neighbor;
	//	}
	//	
	//	virtual CELL_INDEX_TYPE index(cellIterator& it) {
	//		return it.my_location;
	//	}

	//};
	//dVertexIterator my_vertex_iterator_op;
	dBoundaryNeighborIterator my_b_neighbor_it_op;	
	dNeighborIterator my_neighbor_iter_op;
	dBoundaryNeighborCellIterator my_b_neighbor_cell_it_op;	
	dNeighborCellIterator my_neighbor_cell_iter_op;
	

	dCoFacetIterator my_cofacet_iter_op;
	dFacetIterator   my_facet_iter_op;
	dBoundaryCoFacetIterator my_b_cofacet_iter_op;
	dBoundaryFacetIterator my_b_facet_iter_op;



public:

	mscRegularGrid2DImplicitMeshHandler(CELL_INDEX_TYPE x, CELL_INDEX_TYPE y) :
	  X(x), Y(y) {
		  //printf("I GET CALLSED: %llu %llu %llu %llu\n", X, Y, x, y);
		  this->initialize_values();
	  }


	 virtual  ~mscRegularGrid2DImplicitMeshHandler() {
		  		printf("delete: mscRegularGrid2DImplicitMeshHandler \n");
	  }

	CELL_INDEX_TYPE offset_2_pair(unsigned char offset) {
		return offset_list[offset];
	}

	unsigned char pair_2_offset(CELL_INDEX_TYPE diff) {
		return offset_2_loc[diff];
	}

	// General mesh information
	  CELL_INDEX_TYPE num_cells() {
		  return my_num_cells;
	  }

	  CELL_INDEX_TYPE num_cells(DIM_TYPE dim) {
		  return num_d_cells[dim];
	  }

	  DIM_TYPE max_dim() { 
		  return 2;
	  }





	  virtual iteratorOperator& all_cells_iterator(cellIterator& it) {
		  it.my_mesh_handler = this;
		  return  my_all_cells_iterator_operator;
	  }


	  virtual iteratorOperator& d_cells_iterator(DIM_TYPE dim, cellIterator& it) {
		  it.my_mesh_handler = this;
		  it.my_dim = dim;
		  return my_d_cells_iterator_operator;
	  }

	  // Queries regarding individual cells
	  inline void cellid_2_coords(CELL_INDEX_TYPE cellid, CELL_INDEX_TYPE* coords) {
		  coords[0] = cellid % dX;
		  coords[1] = (cellid / dX);
	  }	   
	  inline CELL_INDEX_TYPE coords_2_cellid(const CELL_INDEX_TYPE* coords) {
		  return coords[0] + coords[1]* dX;
	  }

	  inline CELL_INDEX_TYPE num_cells_axis(int i) {
		  if (i == 0) return dX;
		  if (i == 1) return dY;
		  return 0;
	  }
	  BOUNDARY_TYPE boundary_value(CELL_INDEX_TYPE cellid) {
		  CELL_INDEX_TYPE coords[3];
		  cellid_2_coords(cellid, coords);
		  return (BOUNDARY_TYPE) (coords[0] == 0) || (coords[0] == dX-1) ||
			  (coords[1] == 0) || (coords[1] == dY-1) ;
	  }

	  DIM_TYPE dimension(const CELL_INDEX_TYPE& cellid) {
			CELL_INDEX_TYPE coords[3];
			  cellid_2_coords(cellid, coords);	
			  return (coords[0] % 2 + coords[1] % 2);
	  }

	  iteratorOperator& facets(CELL_INDEX_TYPE cellid, cellIterator& it) {
		  it.my_base = cellid;		 
		  it.my_mesh_handler = this;

		  if (this->boundary_value(cellid) == 0) {
			  return this->my_facet_iter_op;
		  } else {
			  //printf("HERHE\n");
			  return this->my_b_facet_iter_op;
		  }
	  }
	  iteratorOperator& cofacets(CELL_INDEX_TYPE cellid, cellIterator& it) {
		  it.my_base = cellid;
		  it.my_mesh_handler = this;

		  if (this->boundary_value(cellid) == 0) {
			  return this->my_cofacet_iter_op;
		  } else {
			  return this->my_b_cofacet_iter_op;
		  }
	  }

	  iteratorOperator& neighbor_vertices(CELL_INDEX_TYPE cellid, cellIterator& it) {
		  it.my_base = cellid;
		  it.my_mesh_handler = this;

		  CELL_INDEX_TYPE coords[3];
		  cellid_2_coords(cellid, coords);	

		  bool nearb = 
			  coords[0] <= 1 || coords[0] >= this->extent_list[0] - 2 ||
			  coords[1] <= 1 || coords[1] >= this->extent_list[2] - 2; // this is right

		  if (! nearb) {
			  return this->my_neighbor_iter_op;
		  } else {
			  return this->my_b_neighbor_it_op;
		  }
	  }
	  iteratorOperator& neighbors(CELL_INDEX_TYPE cellid, cellIterator& it) {
		  it.my_base = cellid;
		  it.my_mesh_handler = this;

		  CELL_INDEX_TYPE coords[3];
		  cellid_2_coords(cellid, coords);	

		  bool nearb = 
			  coords[0] == 0 || coords[0] == this->extent_list[0]  ||
			  coords[1] == 0 || coords[1] == this->extent_list[2] ; // this is right

		  if (! nearb) {
			  return this->my_neighbor_cell_iter_op;
		  } else {
			  return this->my_b_neighbor_cell_it_op;
		  }
	  }

	  //virtual iteratorOperator& vertices(CELL_INDEX_TYPE cellid, cellIterator& it) {
			//it.my_base = cellid;
			//it.my_mesh_handler = this;

			// CELL_INDEX_TYPE coords[3];
			// cellid_2_coords(cellid, coords);	
			//it.my_dim = 0;
			//if (coords[0] % 2 == 1) it.my_dim = it.my_dim | 1;
			//if (coords[1] % 2 == 1) it.my_dim = it.my_dim | 2;

			//it.my_total_size = vertex_offsets[it.my_dim][0];
			//return my_vertex_iterator_op;
	  //}
	
	  virtual void centroid(CELL_INDEX_TYPE cellid, float* coords) {
		  CELL_INDEX_TYPE icoords[3];
		  cellid_2_coords(cellid, icoords);	
		coords[0] = (float) icoords[0];
		coords[1] = (float) icoords[1];
	}
	  friend class mscRegularGrid3DGradientField;
	  
};

#endif
