#ifndef MSC_REGULAR_3D_GRID_IMPLICIT_MESH_HANDLER
#define MSC_REGULAR_3D_GRID_IMPLICIT_MESH_HANDLER

#include "mscIndexTypes.h"
#include "mscBasicIterator.h"
#include "mscRegularGridImplicitMeshHandler.h"
#include "mscBasicIterator.h"

#include <map>

using namespace std;


class mscRegularGrid3DImplicitMeshHandler : public mscRegularGridImplicitMeshHandler {

protected:
	CELL_INDEX_TYPE my_num_cells;
	CELL_INDEX_TYPE X, Y, Z; // the number of vertices in X, Y, Z
	CELL_INDEX_TYPE dX, dY, dZ; // the number of cells in X, Y, Z
	CELL_INDEX_TYPE extent_list[6];

	CELL_INDEX_TYPE num_d_cells[4];
	CELL_INDEX_TYPE offset_list[6];
	CELL_INDEX_TYPE facet_position_list[8][7] ;
	CELL_INDEX_TYPE cofacet_position_list[8][7] ;
	DIM_TYPE facet_dim_list[8][7];
	DIM_TYPE cofacet_dim_list[8][7];


	CELL_INDEX_TYPE& sizes(DIM_TYPE val) {
		if (val == 0) return dX;
		if (val == 1) return dY;
		return dZ;
	}

	map<CELL_INDEX_TYPE, unsigned char> offset_2_loc;

	void initialize_values() {
		// num cells
		dX = 2*X-1;
		dY = 2*Y-1;
		dZ = 2*Z-1;
		my_num_cells = dX * dY * dZ;


		//numbers of cells
		num_d_cells[0] = X*Y*Z;
		CELL_INDEX_TYPE xdir = (X-1)*Y*Z;
		CELL_INDEX_TYPE ydir = X*(Y-1)*Z;
		CELL_INDEX_TYPE zdir = X*Y*(Z-1);
		num_d_cells[1] =  xdir + ydir + zdir;
		xdir = X*(Y-1)*(Z-1);
		ydir = (X-1)*Y*(Z-1);
		zdir = (X-1)*(Y-1)*Z;
		num_d_cells[2] = xdir + ydir + zdir;
		num_d_cells[3] = (X-1)*(Y-1)*(Z-1);


		//offset_list
		offset_list[0] = (1) + (0) * dX + (0) * dX * dY;
		offset_list[1] = (-1) + (0) * dX + (0) * dX * dY;
		offset_list[2] = (0) + (1) * dX + (0) * dX * dY;
		offset_list[3] = (0) + (-1) * dX + (0) * dX * dY;
		offset_list[4] = (0) + (0) * dX + (1) * dX * dY;
		offset_list[5] = (0) + (0) * dX + (-1) * dX * dY;
		
		extent_list[0]= dX-1;
		extent_list[1]= 0;
		extent_list[2]= dY-1;
		extent_list[3]= 0;
		extent_list[4]= dZ-1;
		extent_list[5]= 0;


		
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				for (int k = 0; k < 2; k++) {
					int id = i + 2*j + 4*k;
					DIM_TYPE dim = i + j + k;
					
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
					if (k%2 == 1) {
						facet_dim_list[id][place] = 2;
						facet_position_list[id][place++] = offset_list[4];
						facet_dim_list[id][place] = 2;
						facet_position_list[id][place++] = offset_list[5];
					}
					facet_position_list[id][0] = place - 1;
				}
			}
		}
		
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				for (int k = 0; k < 2; k++) {
					int id = i + 2*j + 4*k;
					DIM_TYPE dim = i + j + k;
					
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
					if (k%2 == 0) {
						cofacet_dim_list[id][place] = 2;
						cofacet_position_list[id][place++] = offset_list[4];
						cofacet_dim_list[id][place] = 2;
						cofacet_position_list[id][place++] = offset_list[5];
					}
					cofacet_position_list[id][0] = place - 1;
				}
			}
		}

		for (int i = 0; i < 6; i++) {
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
			((mscRegularGrid3DImplicitMeshHandler*) it.my_mesh_handler)->cellid_2_coords(it.my_base, it.my_coords);
			it.my_dim = it.my_coords[0]%2 + 2 * (it.my_coords[1]%2) + 4* (it.my_coords[2]%2);
			it.my_total_size = 
				((mscRegularGrid3DImplicitMeshHandler*) it.my_mesh_handler)->facet_position_list[it.my_dim][0];
			//neighbor 
			it.my_neighbor = it.my_base + 
				((mscRegularGrid3DImplicitMeshHandler*) it.my_mesh_handler)->facet_position_list[it.my_dim][it.my_location];	
		}

		void advance(cellIterator& it) {
			it.my_location++;
			it.my_neighbor = it.my_base + 
				((mscRegularGrid3DImplicitMeshHandler*) it.my_mesh_handler)->facet_position_list[it.my_dim][it.my_location];	
		}	

		virtual CELL_INDEX_TYPE value(cellIterator& it) {
			return it.my_neighbor;
		}
		
		virtual CELL_INDEX_TYPE index(cellIterator& it) {
			CELL_INDEX_TYPE temp = 
				((mscRegularGrid3DImplicitMeshHandler*) it.my_mesh_handler)->facet_position_list[it.my_dim][it.my_location];
			return ((mscRegularGrid3DImplicitMeshHandler*) it.my_mesh_handler)->offset_2_loc[temp];
		}

	};

	class dBoundaryFacetIterator : public dFacetIterator {
	public: 

		virtual void begin(cellIterator& it) {
			
			it.my_location = 0;
			((mscRegularGrid3DImplicitMeshHandler*) it.my_mesh_handler)->cellid_2_coords(it.my_base, it.my_coords);
			it.my_dim = it.my_coords[0]%2 + 2 * (it.my_coords[1]%2) + 4* (it.my_coords[2]%2);
			it.my_total_size = 
				((mscRegularGrid3DImplicitMeshHandler*) it.my_mesh_handler)->facet_position_list[it.my_dim][0];
			advance(it);
		}

		void advance(cellIterator& it) {
			it.my_location++;
			

			// if invalid, return
			if (it.my_location > it.my_total_size) return;

			// get the coords of this point
			mscRegularGrid3DImplicitMeshHandler* temp_mesh_handler = 
				((mscRegularGrid3DImplicitMeshHandler*) it.my_mesh_handler);

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
			((mscRegularGrid3DImplicitMeshHandler*) it.my_mesh_handler)->cellid_2_coords(it.my_base, it.my_coords);
			it.my_dim = it.my_coords[0]%2 + 2 * (it.my_coords[1]%2) + 4* (it.my_coords[2]%2);
			it.my_total_size = 
				((mscRegularGrid3DImplicitMeshHandler*) it.my_mesh_handler)->cofacet_position_list[it.my_dim][0];
			//neighbor 
			it.my_neighbor = it.my_base + 
				((mscRegularGrid3DImplicitMeshHandler*) it.my_mesh_handler)->cofacet_position_list[it.my_dim][it.my_location];	
		}

		void advance(cellIterator& it) {
			it.my_location++;
			it.my_neighbor = it.my_base + 
				((mscRegularGrid3DImplicitMeshHandler*) it.my_mesh_handler)->cofacet_position_list[it.my_dim][it.my_location];	
		}	

		virtual CELL_INDEX_TYPE value(cellIterator& it) {
			return it.my_neighbor;
		}
		
		virtual CELL_INDEX_TYPE index(cellIterator& it) {
			CELL_INDEX_TYPE temp = 
				((mscRegularGrid3DImplicitMeshHandler*) it.my_mesh_handler)->cofacet_position_list[it.my_dim][it.my_location];
			return ((mscRegularGrid3DImplicitMeshHandler*) it.my_mesh_handler)->offset_2_loc[temp];
		}

	};

	class dBoundaryCoFacetIterator : public dCoFacetIterator {
	public: 
		virtual void begin(cellIterator& it) {
			
			it.my_location = 0;
			((mscRegularGrid3DImplicitMeshHandler*) it.my_mesh_handler)->cellid_2_coords(it.my_base, it.my_coords);
			it.my_dim = it.my_coords[0]%2 + 2 * (it.my_coords[1]%2) + 4* (it.my_coords[2]%2);
			it.my_total_size = 
				((mscRegularGrid3DImplicitMeshHandler*) it.my_mesh_handler)->cofacet_position_list[it.my_dim][0];
			advance(it);
		}

		void advance(cellIterator& it) {
			it.my_location++;
			
			// if invalid, return
			if (it.my_location > it.my_total_size) return;

			// get the coords of this point
			mscRegularGrid3DImplicitMeshHandler* temp_mesh_handler = 
				((mscRegularGrid3DImplicitMeshHandler*) it.my_mesh_handler);

			while (it.my_location <= it.my_total_size) {
				DIM_TYPE temp_dim = temp_mesh_handler->cofacet_dim_list[it.my_dim][it.my_location]; 
				if (it.my_location % 2 == 0 &&
					it.my_coords[temp_dim] != 0) {
						it.my_neighbor = it.my_base + 
							temp_mesh_handler->cofacet_position_list[it.my_dim][it.my_location];
						return;
				}
				if (it.my_location % 2 == 1 &&
					it.my_coords[temp_dim] != temp_mesh_handler->sizes(temp_dim) - 1) {
						it.my_neighbor = it.my_base + 
							temp_mesh_handler->cofacet_position_list[it.my_dim][it.my_location];
						return;

				}
				it.my_location++;
			}
		}

	};
	
	class dNeighborIterator : public iteratorOperator {
	public: 
		virtual bool valid(cellIterator& it) {
			return it.my_location < 6;
		}

		virtual void begin(cellIterator& it) {
			
			it.my_location = 0;
			it.my_neighbor = it.my_base + 
				((mscRegularGrid3DImplicitMeshHandler*) it.my_mesh_handler)->offset_list[it.my_location] * 2;	
		}

		void advance(cellIterator& it) {
			it.my_location++;
			it.my_neighbor = it.my_base + 
				((mscRegularGrid3DImplicitMeshHandler*) it.my_mesh_handler)->offset_list[it.my_location] * 2;		
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
			return it.my_location < 6;
		}
		virtual void begin(cellIterator& it) {
			
			it.my_location = 0;
			((mscRegularGrid3DImplicitMeshHandler*) it.my_mesh_handler)->cellid_2_coords(it.my_base, it.my_coords);		
			if (it.my_coords[0] + 2 > ((mscRegularGrid3DImplicitMeshHandler*) it.my_mesh_handler)->dX - 1) it.my_location = 1;
			it.my_neighbor = it.my_base + 
				((mscRegularGrid3DImplicitMeshHandler*) it.my_mesh_handler)->offset_list[it.my_location] * 2;	

		}

		virtual void advance(cellIterator& it) {
			it.my_location++;
			

			// if invalid, return
			if (it.my_location > 6) return;

			// get the coords of this point
			mscRegularGrid3DImplicitMeshHandler* temp_mesh_handler = 
				((mscRegularGrid3DImplicitMeshHandler*) it.my_mesh_handler);
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


	dBoundaryNeighborIterator my_b_neighbor_it_op;	
	dNeighborIterator my_neighbor_iter_op;
	dCoFacetIterator my_cofacet_iter_op;
	dFacetIterator   my_facet_iter_op;
	dBoundaryCoFacetIterator my_b_cofacet_iter_op;
	dBoundaryFacetIterator my_b_facet_iter_op;

public:

	mscRegularGrid3DImplicitMeshHandler(CELL_INDEX_TYPE x, CELL_INDEX_TYPE y, CELL_INDEX_TYPE z) :
	  X(x), Y(y), Z(z) {
		  //printf("I GET CALLSED: %llu %llu %llu %llu %llu %llu\n", X, Y, Z, x, y, z);
		  printf("I GET CALLSED: %d %d %d %d %d %d\n", X, Y, Z, x, y, z);
		  this->initialize_values();
	  }


	 virtual  ~mscRegularGrid3DImplicitMeshHandler() {
		  		printf("delete: mscRegularGrid3DImplicitMeshHandler \n");
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
		  return 3;
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
		  coords[1] = (cellid / dX) % dY;
		  coords[2] = (cellid / (dX*dY));
	  }	   

	  CELL_INDEX_TYPE coords_2_cellid(CELL_INDEX_TYPE* coords) {
		  return coords[0] + coords[1]*dX + coords[2]*dX*dY;
	  }

	  BOUNDARY_TYPE boundary_value(CELL_INDEX_TYPE cellid) {
		  CELL_INDEX_TYPE coords[3];
		  cellid_2_coords(cellid, coords);
		  //return (BOUNDARY_TYPE) (coords[0] == 0) + (coords[0] == dX-1) +
			 // (coords[1] == 0) + (coords[1] == dY-1) +
			 // (coords[2] == 0) + (coords[2] == dZ-1);
		  return (BOUNDARY_TYPE) (coords[0] == 0) || (coords[0] == dX-1) ||
			  (coords[1] == 0) || (coords[1] == dY-1) ||
			  (coords[2] == 0) || (coords[2] == dZ-1);
	  }

	  DIM_TYPE dimension(const CELL_INDEX_TYPE& cellid) {
			CELL_INDEX_TYPE coords[3];
			  cellid_2_coords(cellid, coords);	
			  return (coords[0] % 2 + coords[1] % 2 + coords[2] % 2);
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
			  coords[0] <= 1 || coords[0] >= this->extent_list[0] - 1 ||
			  coords[1] <= 1 || coords[1] >= this->extent_list[2] - 1 ||
			  coords[2] <= 1 || coords[2] >= this->extent_list[4] - 1;

		  if (! nearb) {
			  return this->my_neighbor_iter_op;
		  } else {
			  return this->my_b_neighbor_it_op;
		  }
	  }
	  
	 virtual void centroid(CELL_INDEX_TYPE cellid, float* coords) {
		  CELL_INDEX_TYPE icoords[3];
		  cellid_2_coords(cellid, icoords);	
		coords[0] = (float) icoords[0];
		coords[1] = (float) icoords[1];
		coords[2] = (float) icoords[2];

	}

	 CELL_INDEX_TYPE extent(DIM_TYPE i) {
		return extent_list[i*2];
	 }

	  friend class mscRegularGrid3DGradientField;
	  
};

#endif