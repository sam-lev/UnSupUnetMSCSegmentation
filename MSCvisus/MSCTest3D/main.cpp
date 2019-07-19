#include <stdio.h>
#include "mscArrayFactory.h"
#include "mscDumbGradientField.h"
#include "mscRegularRawDataHandler.h"
#include "mscRegularGrid3DImplicitMeshHandler.h"
#include "mscDumbStoringMeshFunction.h"
#include "mscSimpleGradientBuilder.h"
#include "mscSimpleRandomGradientBuilder.h"
#include "mscSimpleGradientUsingAlgorithms.h"
#include "mscRegularGrid3DGradientField.h"
#include "mscRegularGrid3DMeshFunction.h"
#include "mscTwoWay3DGradientBuilder.h"
#include "mscConvergentGradientBuilder.h"
#include "mscNegatingMeshFunction.h"
#include "mscComplementMeshHandler.h"
#include "mscModifiedBoundaryMeshHandler.h"
#include "mscCombinationGradientField.h"
#include "mscNegatingDataHandler.h"
#include "mscShortcuttingConvergentGradientBuilder.h"
#include "mscBasicMSC.h"

#include "mscSegmentation.h"
#include "mscGenericComplexGenerator.h"
#include "mscGeomQuadSurface.h"
#include "mscGeometryGenerator.h"
#include "mscSimpleConstrainedGradientBuilder.h"

#include <map>

using namespace std;

class testme : public vector<node<float>*> {
};

struct xyz { float coords[3]; };

int main(int argc, char** argv) {

	// print command line args for sanity
	for (int i =0; i < argc; i++) {
		printf("%d=%s\n", i, argv[i]);
	}

	if (argc < 6) {
		printf("Usage: test filename x y z computemode\n\
			   \t 0 = randomized\n\
			   \t 1 = greedy\n\
			   \t 2 = 1 pass convergent (Morse complex)\n\
			   \t 3 = 2 pass convergent (inverse then Morse Complex)\n\
			   \t 4 = 2 pass convergent (Morse the inverse)\n");
		return 1;
	}

	// read in the command line arguments
	char* filename = argv[1];
	int X, Y, Z;
	sscanf(argv[2], "%d", &X);
	sscanf(argv[3], "%d", &Y);
	sscanf(argv[4], "%d", &Z);
	int gradienttype;
	sscanf(argv[5], "%d", &gradienttype);

	float cancelto = 5.0;
	if (argc > 6) 
		sscanf(argv[6], "%f", &cancelto);


	// declare array factory
	mscArrayFactory a(REGULAR_ARRAY);

	// load data of size X*Y*Z
	mscRegularRawDataHandler<float>* test_data;
	test_data = new mscRegularRawDataHandler<float>();
	if (! test_data->load_data(filename, X*Y*Z, &a)) printf("ERROR: DID NOT LOAD DATA\n");

	// some algorithms need the negative, so here it is
	mscNegatingDataHandler<float>* ndh = 
		new mscNegatingDataHandler<float>(test_data);

	// create mesh handler
	mscRegularGrid3DImplicitMeshHandler* bmsh = new mscRegularGrid3DImplicitMeshHandler(X,Y,Z);

	// create mesh function
	mscRegularGrid3DMeshFunction<float>* bmf = new mscRegularGrid3DMeshFunction<float>(test_data, bmsh, &a);
	bmf->initialize();

	//create data structure to store discrete gradient field
	mscRegularGrid3DGradientField* bgf = 
		new mscRegularGrid3DGradientField(bmsh, &a);

	char gradname[1024];
	sprintf(gradname,"%s.grad", filename);
	//if (! bgf->load_from_file(gradname)){
	
	
	// now pick the gradient construction algorithm and compute the discrete gradient field!
	printf("Computing discrete gradient: \n");
	if (gradienttype == 0) {
		// uses a randomized approach 
		printf("using randomized gradient\n");
		mscSimpleRandomGradientBuilder<float>* mscrb = 
			new mscSimpleRandomGradientBuilder<float>(bmf, bmsh, bgf, &a);
		mscrb->computeGradient();
	} else if (gradienttype == 1) {
		// constructs steepest descent discrete gradient
		printf("using greedy gradient\n");
		mscSimpleGradientBuilder<float>* mscb = 
			new mscSimpleGradientBuilder<float>(bmf, bmsh, bgf, &a);
		mscb->computeGradient();
	} else if (gradienttype == 2) {
		// integrates probabilities for accurate basin boundaries then accurate mountain boundaries
		printf("using convergent gradient\n");
		mscConvergentGradientBuilder<float>* msccb = 
			new mscConvergentGradientBuilder<float>(bmf, bmsh, bgf, &a);
		msccb->computeGradient();
		delete msccb;
		mscComplementMeshHandler* cmh = 
			new mscComplementMeshHandler(bmsh);
		mscNegatingMeshFunction<float>* nmf = 
			new mscNegatingMeshFunction<float>(bmf);
		mscModifiedBoundaryMeshHandler* mbmh = 
			new mscModifiedBoundaryMeshHandler(cmh, bgf);
		mscRegularGrid3DGradientField* bgf2 = 
			new mscRegularGrid3DGradientField(bmsh, &a);
		mscConvergentGradientBuilder<float>* msccb2 = 
			new mscConvergentGradientBuilder<float>(nmf, mbmh, bgf2, &a);
		msccb2->computeGradient();
		mscRegularGrid3DGradientField* tempgf = bgf;
		cellIterator it;
		iteratorOperator& all_cells = bmsh->all_cells_iterator(it);
		for (all_cells.begin(it); all_cells.valid(it); all_cells.advance(it)) {
			CELL_INDEX_TYPE cid = all_cells.value(it);
			bgf2->set_dim_asc_man(cid, bgf->get_dim_asc_man(cid));
		}
		bgf = bgf2;
		bgf->resetMeshHandler(bmsh);
	}else if (gradienttype == 5) {
		// convergent graidnet builder that uses a more efficient membership distribution
		printf("using NEW convergent gradient\n");
		mscShortcuttingConvergentGradientBuilder<float>* mscSHcb = 
			 new mscShortcuttingConvergentGradientBuilder<float>(bmf, bmsh, bgf, &a);
		mscSHcb->computeGradient();
		delete mscSHcb;
		mscComplementMeshHandler* cmh = 
			new mscComplementMeshHandler(bmsh);
		mscNegatingMeshFunction<float>* nmf = 
			new mscNegatingMeshFunction<float>(bmf);
		mscModifiedBoundaryMeshHandler* mbmh = 
			new mscModifiedBoundaryMeshHandler(cmh, bgf);
		mscRegularGrid3DGradientField* bgf2 = 
			new mscRegularGrid3DGradientField(bmsh, &a);
		mscShortcuttingConvergentGradientBuilder<float>* msccb2 = 
			new mscShortcuttingConvergentGradientBuilder<float>(nmf, mbmh, bgf2, &a);
		//	 msccb2->computeGradient();
		msccb2->computeGradient();
		mscRegularGrid3DGradientField* tempgf = bgf;
		cellIterator it;
		iteratorOperator& all_cells = bmsh->all_cells_iterator(it);
		for (all_cells.begin(it); all_cells.valid(it); all_cells.advance(it)) {
			CELL_INDEX_TYPE cid = all_cells.value(it);
			bgf2->set_dim_asc_man(cid, bgf->get_dim_asc_man(cid));
		}
		bgf = bgf2;
		bgf->resetMeshHandler(bmsh);
	} else if (gradienttype == 3){
		// first get accurate "mountains" then accurate basins
		printf("using convergent2 gradient - 2 pass\n");
		mscComplementMeshHandler* cmh = 
			new mscComplementMeshHandler(bmsh);
		mscNegatingMeshFunction<float>* nmf = 
			new mscNegatingMeshFunction<float>(bmf);
		mscRegularGrid3DGradientField* bgf2 = 
			new mscRegularGrid3DGradientField(bmsh, &a);	   
		mscConvergentGradientBuilder<float>* msccb2 = 
			new mscConvergentGradientBuilder<float>(nmf, cmh, bgf2, &a);
		msccb2->computeGradient();	   
		mscModifiedBoundaryMeshHandler* mbmh = 
			new mscModifiedBoundaryMeshHandler(bmsh, bgf2);

		mscConvergentGradientBuilder<float>* msccb3 = 
			new mscConvergentGradientBuilder<float>(bmf, mbmh, bgf, &a);
		msccb3->computeGradient();
	}  else if (gradienttype == 4){
		// only accurate basins
		printf("using convergent gradient -1 pass\n");
		mscConvergentGradientBuilder<float>* msccb = 
			new mscConvergentGradientBuilder<float>(bmf, bmsh, bgf, &a);
		msccb->computeGradient();
	}  else if (gradienttype == 6/* && argc > 6*/) {
		// READ IN SEGMENTATION
		//mscb->computeGradient();
		//USE_SEG = true;
		FILE* fseg = fopen("source_dest2.raw", "rb");
		int numV = X * Y * Z;
		int numC = (X-1) * (Y-1)* (Z-1);

		int* ascID = new int[numV];
		int* dscID = new int[numC];
		fread(ascID, sizeof(int), numV, fseg);
		fread(dscID, sizeof(int), numC, fseg);
		fclose(fseg);

		mscRegularGrid3DGradientField* G_mscg_TEMP = 
			new mscRegularGrid3DGradientField(bmsh, &a);

		
		mscPreClassifier* classes = new mscPreClassifier(ascID, dscID, bmsh);

		// THIS IS A MISNOMER all it actually does is set 
		mscSimpleConstrainedGradientBuilder<float>* scgb = 
			new mscSimpleConstrainedGradientBuilder<float>(classes, bmf,bmsh,G_mscg_TEMP, &a);
		scgb->computeGradient();


		mscModifiedBoundaryMeshHandler* mbmh = 
			new mscModifiedBoundaryMeshHandler(bmsh, G_mscg_TEMP);

		mscSimpleGradientBuilder<float>* mscb3 = 
			new mscSimpleGradientBuilder<float>(bmf, mbmh, bgf, &a);
		mscb3->computeGradient();

		//mscSimpleConstrainedGradientBuilder<float>* scgb2 = 
		//	new mscSimpleConstrainedGradientBuilder<float>(classes, bmf,mbmh,bgf, &a);
		//scgb2->computeGradient();



		bgf->output_to_renderer(filename);

		return 1;
		//G_mscmh = mbmh;

	} 





	//}
	// testing
	int cc[4]; for(int i=0;i<4;i++) cc[i]=0;
		int critcount = 0;
	printf("Testing:\n");
	for (CELL_INDEX_TYPE i = 0; i < bmsh->num_cells(); i++) {
		if (! bgf->get_assigned(i)) {
			printf("ERROR: cell %d is not assigned: %d\n", i, 
				bmsh->dimension(i));
		}

		if (bgf->get_critical(i)) {
			critcount++;
			cc[bmsh->dimension(i)]++;
		}
	}
	printf("critcount = %d\n", critcount);
	for (int i=0; i < 4; i++) printf("#index-%d = %d\n", i, cc[i]);

	bgf->output_to_renderer(filename);
	bmf->output_to_renderer(filename);

	return 1;

	BasicMSC<float>* msc = new BasicMSC<float>(bgf,bmsh,bmf);
	msc->ComputeFromGrad();

	// this is how we iterate through nodes
	// as an example we compute lowest min and highest max
	float lower_bound = (*msc->nodes.begin()).second->value; 
	float upper_bound = (*msc->nodes.begin()).second->value;
	for (map<CELL_INDEX_TYPE, node<float>* >::iterator it = msc->nodes.begin(); 
		it != msc->nodes.end(); it++){
			node<float>* n = (*it).second;
			if (n->value < lower_bound) lower_bound = n->value;
			if (n->value > upper_bound) upper_bound = n->value;
	}
	printf(" lower = %f, upper = %f\n", lower_bound, upper_bound);
	
	// now we will construct a simplification hierarchy up to 10% of function difference
	float tenpercent = (upper_bound - lower_bound) * 0.1;
	msc->ComputeHeirarchy(cancelto);

	// now set the selection persistence to 1% of function range
	msc->setAbsolutePersistence(cancelto  );
	//msc->setFirstNonzeroPersistence(0.00000001);


	printf("computing segmentation\n");
	mscComputeBasinSegmentation<float>* mCBS = 
		new mscComputeBasinSegmentation<float>(msc, bmsh);
	printf("gothere1\n");
	mCBS->Compute();
	printf("done\n");
	mscGenericComplexGenerator* mGCG = 
		new mscGenericComplexGenerator(bmsh, mCBS->GetSeg());
	mGCG->Generate();

	printf("now about to generate surfaces\n");
	mscGeomQuadSurface stuff;
	mscRegular3dGridSurfExtractor* surfgen = new mscRegular3dGridSurfExtractor(bmsh);
	
	//set<CELL_INDEX_TYPE> ALL_IDS;
	
	// compute centers of mass
	map<CELL_INDEX_TYPE, mscGenericCell*>& cells2 = mGCG->Complex()->GetCells(0);
	map<mscGenericCell*, xyz> centers;
	for(map<CELL_INDEX_TYPE, mscGenericCell*>::iterator it = cells2.begin(); it != cells2.end(); it++) {
		mscGenericCell* c = (*it).second;
		xyz center;
		for (int j=0;j<3;j++) center.coords[j]=0;
		set<CELL_INDEX_TYPE>& ids =c->ids;
		float count = (float) 1.0 / (float) ids.size();
		for (set<CELL_INDEX_TYPE>::iterator it2 = ids.begin(); it2 != ids.end(); it2++) {
			xyz other;
			bmsh->centroid(*it2, other.coords);
			for (int j=0;j<3;j++) {
				center.coords[j] += count * other.coords[j];
			}
		}

		centers[c]=center;
	}
	
	char sname[1024];
	map<CELL_INDEX_TYPE, mscGenericCell*>& cells = mGCG->Complex()->GetCells(1);
	for(map<CELL_INDEX_TYPE, mscGenericCell*>::iterator it = cells.begin(); it != cells.end(); it++) {
		mscGeomQuadSurface s;

		
		int countxz=0;
		int county=0;
		set<CELL_INDEX_TYPE>& ids = (*it).second->ids;
		float lowestvalue = bmf->cell_value(*ids.begin());
		for (set<CELL_INDEX_TYPE>::iterator it2 = ids.begin(); it2 != ids.end(); it2++) {
			if (bmsh->dimension(*it2) != 1) continue;
			CELL_INDEX_TYPE coords[3];
			bmsh->cellid_2_coords(*it2, coords);
			//printf("%d, %d, %d\n", coords[0], coords[1], coords[2]);
			if (coords[2]%2 == 1) { 
				countxz++;
			} else { 
				county++;
			}
			float tmp = bmf->cell_value(*it2);
			if (tmp < lowestvalue) lowestvalue = tmp;
		}
		float ratio;


		int pos = 0;
		mscGenericCell* negs[2];
		for (set<mscGenericCell*>::iterator it2 = (*it).second->neighbors.begin();
			it2 != (*it).second->neighbors.end(); it2++) {
				mscGenericCell* nb = *it2;
				if (nb->dim == 0) {
					negs[pos++]=nb;
				}
		}
		if (pos != 2) {
			//then this does not separate 2 asc mans so return
			continue;
		} 
		xyz& a = centers[negs[0]];
		xyz& b = centers[negs[1]];
		float mydist1 = sqrt(
			(a.coords[0] - b.coords[0]) *(a.coords[0] - b.coords[0]) +
			(a.coords[1] - b.coords[1]) *(a.coords[1] - b.coords[1]));
		float mydist2 = sqrt(
			(a.coords[2] - b.coords[2]) *(a.coords[2] - b.coords[2]) );

		float nratio;
		
		if (mydist2 == 0) {
			nratio = 9999;
		} else {
			nratio = mydist1/mydist2;
		}


		if (county == 0) {
			ratio = 9999;
		} else {
			ratio = (float) countxz / (float) county;
		}

		if (lowestvalue < 118) continue;
		if (ratio < 1.0) continue;
		//if (ratio > 1.0) {
		//	ALL_IDS.insert(ids.begin(), ids.end());
		//}
		surfgen->CreateSurface((*it).second->ids, s, 1);
		s.smooth(5);
		sprintf(sname,"surfs/surf_%f_%f_%f_%d.obj",mydist2, ratio, lowestvalue, countxz + county);
		s.dumpInObjWN(sname);
	}
	printf("donesurf\n");
	//mscGeomQuadSurface s2;
	//surfgen->CreateSurface(ALL_IDS, s2, 1);
	//s2.smooth(5);
	//s2.dumpInObjWN("golobalsurf.obj");

	//// now iterate through each node again, only doing things with "alive" ones
	//for (map<CELL_INDEX_TYPE, node<float>* >::iterator it = msc->nodes.begin(); 
	//	it != msc->nodes.end(); it++){
	//		node<float>* n = (*it).second;
	//		// if it's not alive in the current level of hierarchy skip
	//		if (! msc->isAlive(n)) continue;

	//		// NOW DO WHATEVER WITH n
	//		// (1) compute it's geometric embedding
	//		if (n->index == 1) {
	//			set<CELL_INDEX_TYPE> result; // space to store the result
	//			msc->fillGeometry(n, result); // compute the set of cells in the asc/dsc man of the node
	//			// now result contains all the cell indices in the ascending or descending 3 or 2 manifold of n
	//			// this now shows how that could be output as a file
	//			char fname[1024];
	//			sprintf(fname, "surfs/node_%d_surf.obj", n->cellid); // give each arc a file name
	//			mscGeomQuadSurface s;
	//			surfgen->CreateSurface(result, s);
	//			s.smooth(3);
	//			s.dumpInObjWN(fname);
	//		}
	//}

	//// iterating through arcs is the same
	//int counter = 0;
	//for (vector<arc<float>*>::iterator it = msc->arcs.begin(); 
	//	it != msc->arcs.end(); it++) {
	//		arc<float>* a = (*it);
	//		if (! msc->isAlive(a)) continue;
	//		// now do stuff with arc 

	//		//list of points in arc
	//		vector<CELL_INDEX_TYPE> result; // space to store the result
	//		msc->fillGeometry(a, result); // ordered list of cells from start to end of line
	//		// now result contains list of points in the geometric representation of the arc
	//		// this now shows how that could be output as a file
	//		char fname[1024];
	//		sprintf(fname, "arc_%d.xyz", counter++); // give each arc a file name
	//		FILE* fout = fopen(fname, "w");
	//		for (int i = 0; i < result.size(); i++) {
	//			float coords[3]; // space to store the coords of the centriod of a cell
	//			bmsh->centroid(result[i], coords); // compute centroid
	//			fprintf(fout, "%f %f %f\n", coords[0], coords[1], coords[2]); // write to file
	//		}
	//		fclose(fout); // done writing this arc
	//}

	return 1;


}
