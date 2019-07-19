#ifndef MSC_REGULAR_GRID_3D_DECOMPOSITION
#define MSC_REGULAR_GRID_3D_DECOMPOSITION

#include "mscIndexTypes.h"
#include "mscBasicMeshDecomposition.h"
#include "mscSimplePointerDataHandler.h"
#include "mscRegularGrid3DMeshFunction.h"
#include "mscRegularGrid3DGradientField.h"


class mscRegularGrid3DGradientField;

template<typename dtype>
class mscRegularGrid3DDecomposition : public mscBasicMeshDecomposition<dtype> {
protected:

	struct localExtents {
		INT_TYPE starts[3];
		INT_TYPE sizes[3];
		bool loaded;
		dtype* local_data;
		mscBasicDataHandler<dtype>* data_handler;
		mscBasicMeshFunction<dtype>* mesh_function;
		mscBasicMeshHandler* mesh_handler;
		mscBasicGradientField* grad;
	};

	localExtents globalSizes;
	localExtents desiredBlock;
	INT_TYPE mNumInDir[3];
	INT_TYPE mNumBlocks;
	char* mFilename;
	char* mOutGradFileName;
	char* mOutDatFileName;
	localExtents* blocks;

	bool isGlobalBoundary(int nx, int ny, int nz, localExtents& l) {
		int x = nx + l.starts[0]*2;
		int y = ny + l.starts[1]*2;
		int z = nz + l.starts[2]*2;
		return x == 0 || x == ((globalSizes.sizes[0]*2-1)-1) ||
			y == 0 || y == ((globalSizes.sizes[1]*2-1)-1) ||
			z == 0 || z == ((globalSizes.sizes[2]*2-1)-1);
	}

	virtual dtype* loadSubBlockFromFile(INT_TYPE block_id) {
		localExtents& l = blocks[block_id];
		// load values from file
		CELL_INDEX_TYPE num_data = l.sizes[0] * l.sizes[1] * l.sizes[2];
		//printf("reading %d data\n", num_data);
		dtype* local_data = new dtype[num_data];

		FILE* fin = fopen( mFilename, "rb");
		// seek to start
		CELL_INDEX_TYPE t_file_start = l.starts[2] * globalSizes.sizes[0] * 
			globalSizes.sizes[1];

		fseek(fin, sizeof(dtype) * t_file_start, SEEK_SET);
		dtype* curr_loc = local_data;
		printf("get here 1 : %d\n", t_file_start);

		int readcount = 0;
		for (int z = 0; z < l.sizes[2]; z++) {
				fseek(fin, sizeof(dtype) * l.starts[1] * globalSizes.sizes[0], SEEK_CUR);
			for(int y = 0; y < l.sizes[1]; y++) {
				//skip empty y space

				//for(int x = 0; x < l.sizes[0]; x++) {
				// skip empty x space
					fseek(fin, sizeof(dtype) * l.starts[0], SEEK_CUR);
					fread(curr_loc, sizeof(dtype), l.sizes[0], fin);
					readcount += l.sizes[0];

					curr_loc += /*sizeof(dtype) **/ l.sizes[0];
					fseek(fin, sizeof(dtype) * 
						(globalSizes.sizes[0] - (l.starts[0] + l.sizes[0])), SEEK_CUR);
				//}

			}
				fseek(fin, sizeof(dtype) * globalSizes.sizes[0] *
					(globalSizes.sizes[1]- (l.starts[1] + l.sizes[1])), SEEK_CUR);
		}
		fclose(fin);
		printf("load returned: read %d items\n", readcount);
		return local_data;
	}


public:

	mscRegularGrid3DDecomposition(
		INT_TYPE global_X,
		INT_TYPE global_Y,
		INT_TYPE global_Z,
		INT_TYPE size_X,
		INT_TYPE size_Y,
		INT_TYPE size_Z) {
			globalSizes.starts[0] = 0;
			globalSizes.starts[1] = 0;
			globalSizes.starts[2] = 0;
			globalSizes.sizes[0] = global_X;
			globalSizes.sizes[1] = global_Y;
			globalSizes.sizes[2] = global_Z;
			desiredBlock.sizes[0] = size_X;
			desiredBlock.sizes[1] = size_Y;
			desiredBlock.sizes[2] = size_Z;
	}

	void printBlockInfo(INT_TYPE block_id) {
		printf("%d:(%d, %d, %d)->(%d, %d, %d)\n",
			block_id,
			blocks[block_id].starts[0],
			blocks[block_id].starts[1],
			blocks[block_id].starts[2],
			blocks[block_id].starts[0] + blocks[block_id].sizes[0],
			blocks[block_id].starts[1] + blocks[block_id].sizes[1],
			blocks[block_id].starts[2] + blocks[block_id].sizes[2]
		);
	}

	void decompose() {
		// figure out number in each dimension
		for (int i = 0; i < 3; i++) {
			mNumInDir[i] = globalSizes.sizes[i] / desiredBlock.sizes[i] +
				(globalSizes.sizes[i] % desiredBlock.sizes[i] ? 1 : 0);
		}

		mNumBlocks = mNumInDir[0] * mNumInDir[1] * mNumInDir[2];
		blocks = new localExtents[mNumBlocks];

		// now simply create the list
		for (int z = 0; z < mNumInDir[2]; z++) {
			for (int y = 0; y < mNumInDir[1]; y++) {
				for (int x = 0; x < mNumInDir[0]; x++) {
					localExtents& l = blocks[x + y * mNumInDir[0] + z * mNumInDir[0] * mNumInDir[1]];
					l.loaded = false;

					l.starts[0] = x * desiredBlock.sizes[0];
					l.starts[1] = y * desiredBlock.sizes[1];
					l.starts[2] = z * desiredBlock.sizes[2];
					l.sizes[0] = desiredBlock.sizes[0] + 1;
					if (x == mNumInDir[0] - 1) 
						l.sizes[0] = globalSizes.sizes[0] - l.starts[0];

					l.sizes[1] = desiredBlock.sizes[1] + 1;
					if (y == mNumInDir[1] - 1) 
						l.sizes[1] = globalSizes.sizes[1] - l.starts[1];

					l.sizes[2] = desiredBlock.sizes[2] + 1;
					if (z == mNumInDir[2] - 1) 
						l.sizes[2] = globalSizes.sizes[2] - l.starts[2];
				}
			}
		}

	}

	void setInputFileName(char* filename) {
		mFilename = filename;
	}	   
	void setOutputFileNames(char* grad_name, char* dat_name) {
		mOutGradFileName = grad_name;
		mOutDatFileName = dat_name;
	}

	int numBlocks() { return mNumBlocks; }




	void loadBlock(INT_TYPE block_id) {

		localExtents& l = blocks[block_id];
		if (l.loaded) return;

		// load values from file
		CELL_INDEX_TYPE num_data = l.sizes[0] * l.sizes[1] * l.sizes[2];
		dtype* local_data = this->loadSubBlockFromFile(block_id);
		mscSimplePointerDataHandler<dtype>* data_handler = 
			new mscSimplePointerDataHandler<dtype>();
		data_handler->set_data(local_data);
		//data_handler->logify();

		printf("SIZES: %d %d %d\n", l.sizes[0], l.sizes[1], l.sizes[2]);
		mscRegularGrid3DImplicitMeshHandler* mesh_handler = 
			new mscRegularGrid3DImplicitMeshHandler(l.sizes[0], l.sizes[1], l.sizes[2]);

		mscArrayFactory* array_factory = new mscArrayFactory(REGULAR_ARRAY);
		mscRegularGrid3DMeshFunction<dtype>* mesh_function = 
			new mscRegularGrid3DMeshFunction<dtype>(data_handler, mesh_handler, array_factory);
		mesh_function->initialize();

		mscRegularGrid3DGradientField* gradient = 
			new mscRegularGrid3DGradientField(mesh_handler, array_factory);


		l.data_handler = data_handler;
		l.mesh_function = mesh_function;
		l.mesh_handler = mesh_handler;
		l.grad = gradient;
		
		l.loaded = true;
		// now we read all sub-blocks
		// create the data function

	}
	void unloadBlock(INT_TYPE block_id) {
		printf("deleteing %d\n", block_id);
		localExtents& l = blocks[block_id];
		if (! l.loaded) return;


		delete l.data_handler;
		delete l.mesh_function;
		delete l.mesh_handler;
		delete l.grad;
		l.loaded = false;
		printf("done\n");

	}


	void writeOutputs(INT_TYPE block_id) {
		struct bitfield {
			unsigned char assigned : 1;
			unsigned char flag : 1;
			//unsigned char critical : 1;
			//unsigned char insorter : 1;
			//unsigned char dimA : 3;
			unsigned char pair : 3;
			unsigned char ldir : 3;
			//unsigned char empty_flag : 1;
		};
		localExtents& l = blocks[block_id];
		// load values from file
		CELL_INDEX_TYPE num_data = 
			(l.sizes[0]*2 - 1) * 
			(l.sizes[1]*2 - 1) * 
			(l.sizes[2]*2 -1 );

		FILE* fgrad = fopen( mOutGradFileName, "rb+");
		FILE* fdat = fopen( mOutDatFileName, "rb+");

		// seek to start
		CELL_INDEX_TYPE t_file_start = 
			(l.starts[2]*2) * 
			(globalSizes.sizes[0]*2 - 1) * 
			(globalSizes.sizes[1]*2 - 1);
		printf("z seek+ :%d\n", t_file_start);
		fseek(fgrad, sizeof(bitfield) * t_file_start, SEEK_SET);
		fseek(fdat, sizeof(dtype) * t_file_start, SEEK_SET);
		CELL_INDEX_TYPE current_cell = 0;

		CELL_INDEX_TYPE counter = 0;
		CELL_INDEX_TYPE offsetcounter = t_file_start;
		for (int z = 0; z < (l.sizes[2]*2-1); z++) {
				fseek(fgrad, sizeof(bitfield) * (l.starts[1]*2) * (globalSizes.sizes[0]*2-1), SEEK_CUR);
				fseek(fdat, sizeof(dtype) * (l.starts[1]*2) * (globalSizes.sizes[0]*2-1), SEEK_CUR);
				offsetcounter +=sizeof(dtype) * (l.starts[1]*2) * (globalSizes.sizes[0]*2-1);		
				
				for(int y = 0; y < (l.sizes[1]*2-1); y++) {
				//printf("y seek+ :%d\n",(l.starts[1]*2) * (globalSizes.sizes[0]*2-1));


				//for(int x = 0; x < l.sizes[0]; x++) {

				//printf("x seek+ :%d",(l.starts[0]*2));
					fseek(fgrad, sizeof(bitfield) * (l.starts[0]*2), SEEK_CUR);
					fseek(fdat, sizeof(dtype) * (l.starts[0]*2), SEEK_CUR);
					offsetcounter+=sizeof(dtype) * (l.starts[0]*2);


					//printf("writing %d\n", sizeof(bitfield) * (l.sizes[0]*2-1));
					for (int k = 0; k < (l.sizes[0]*2-1); k++) {

						if (l.grad->get_assigned(current_cell) == 0) {
							printf("ERROR: not assigned\n");
						}

						bitfield temp;
						temp.assigned = l.grad->get_assigned(current_cell);
						temp.flag = isGlobalBoundary(k, y, z, l);

						if (l.grad->get_critical(current_cell) ) {
							temp.pair = 7;
						} else {
						temp.pair = 
							(unsigned char) ((mscRegularGrid3DImplicitMeshHandler*) l.mesh_handler)->pair_2_offset(
							l.grad->get_pair(current_cell) - current_cell);
						}
						temp.ldir = (unsigned char) l.grad->get_dim_asc_man(current_cell);
						dtype val = l.mesh_function->cell_value(current_cell);

					//   	if (current_cell % 1000 == 0)
					//printf("cell[%d]:(%d, %d, %d, %d)\n", (int) current_cell,
					//temp.pair, 
					//temp.ldir,
					//temp.flag,
					//temp.assigned);
	   
	   		
						fwrite(&temp, sizeof(bitfield), 1, fgrad);
						fwrite(&val, sizeof(dtype), 1, fdat);
						counter++;
						current_cell++;
						if (ferror(fdat)) printf("ERROR\n");
					}
					//printf("x seek- :%d",((globalSizes.sizes[0]*2-1) - (l.starts[0]*2 + (l.sizes[0]*2-1))));

					fseek(fgrad, sizeof(bitfield) * 
						((globalSizes.sizes[0]*2-1) - (l.starts[0]*2 + (l.sizes[0]*2-1))), SEEK_CUR);
					fseek(fdat, sizeof(dtype) * 
						((globalSizes.sizes[0]*2-1) - (l.starts[0]*2 + (l.sizes[0]*2-1))), SEEK_CUR);
					offsetcounter +=sizeof(dtype) * 
						((globalSizes.sizes[0]*2-1) - (l.starts[0]*2 + (l.sizes[0]*2-1)));
				//}
				//printf("y seek- :%d\n",(globalSizes.sizes[0]*2-1) *
				//	((globalSizes.sizes[1]*2-1)- (l.starts[1]*2 + (l.sizes[1]*2-1))));


			}
				fseek(fgrad, sizeof(bitfield) * (globalSizes.sizes[0]*2-1) *
					((globalSizes.sizes[1]*2-1)- (l.starts[1]*2 + (l.sizes[1]*2-1))), SEEK_CUR);

				fseek(fdat, sizeof(dtype) * (globalSizes.sizes[0]*2-1) *
					((globalSizes.sizes[1]*2-1)- (l.starts[1]*2 + (l.sizes[1]*2-1))), SEEK_CUR);


				offsetcounter += sizeof(dtype) * (globalSizes.sizes[0]*2-1) *
					((globalSizes.sizes[1]*2-1)- (l.starts[1]*2 + (l.sizes[1]*2-1)));		
		}
		printf("counter=%d    offsetcounter=%d\n", counter, offsetcounter);
		fclose(fdat);
		fclose(fgrad);
	}

	mscBasicMeshFunction<dtype>* getMeshFunction(INT_TYPE block_id) { 
		return blocks[block_id].mesh_function; 
	}
	mscBasicMeshHandler* getMeshHandler(INT_TYPE block_id) { 
		return blocks[block_id].mesh_handler; 
	}
	mscBasicGradientField* getGradientField(INT_TYPE block_id) { 
		return blocks[block_id].grad; 
	}
	mscBasicDataHandler<dtype>* getDataHandler(INT_TYPE block_id) { 
		return blocks[block_id].data_handler; 
	}
	void writeToDisk(char* filename) {}
};


#endif
