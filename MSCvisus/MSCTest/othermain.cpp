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


#include <map>


int main2(int argc, char** argv) {

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

	char* filename = argv[1];
	int X, Y, Z;
	sscanf(argv[2], "%d", &X);
	sscanf(argv[3], "%d", &Y);
	sscanf(argv[4], "%d", &Z);

	int userand;
	sscanf(argv[5], "%d", &userand);

	mscArrayFactory a(REGULAR_ARRAY);

	// load data
	printf("loading data...\n");
   mscRegularRawDataHandler<float>* test_data;
   test_data = new mscRegularRawDataHandler<float>();
   test_data->load_data(filename, X*Y*Z, &a);

   mscNegatingDataHandler<float>* ndh = 
	   new mscNegatingDataHandler<float>(test_data);
   printf("Done!\n");
	printf("Computing discrete gradient: \n");


	if (userand == 5) {
		test_data->negate();
	}
   // create mesh handler
   mscRegularGrid3DImplicitMeshHandler* bmsh = new mscRegularGrid3DImplicitMeshHandler(X,Y,Z);

   // create mesh function
   mscRegularGrid3DMeshFunction<float>* bmf = new mscRegularGrid3DMeshFunction<float>(test_data, bmsh, &a);
  
   
   bmf->initialize();


   mscRegularGrid3DGradientField* bgf = 
	   new mscRegularGrid3DGradientField(bmsh, &a);

   // test the gradient builder!
   mscSimpleGradientBuilder<float>* mscb = 
	   new mscSimpleGradientBuilder<float>(bmf, bmsh, bgf, &a);

   mscSimpleRandomGradientBuilder<float>* mscrb = 
	   new mscSimpleRandomGradientBuilder<float>(bmf, bmsh, bgf, &a);

   mscConvergentGradientBuilder<float>* msccb = 
	   new mscConvergentGradientBuilder<float>(bmf, bmsh, bgf, &a);

   if (userand == 0) {
		printf("--using randomized gradient\n");
	   mscrb->computeGradient();
   } else if (userand == 1) {
	   printf("--using greedy gradient\n");
		mscb->computeGradient();
   } else if (userand == 2) {
	   printf("--using 2-pass convergent gradient, primal then complement\n");
	   msccb->computeGradient();
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

   } else if (userand == 3){
	   	   printf("--using 2-pass convergent gradient, complement then primal\n");


	   mscComplementMeshHandler* cmh = 
		   new mscComplementMeshHandler(bmsh);

	   mscNegatingMeshFunction<float>* nmf = 
		   new mscNegatingMeshFunction<float>(bmf);
	   mscRegularGrid3DGradientField* bgf2 = 
		   new mscRegularGrid3DGradientField(bmsh, &a);	   

	     mscConvergentGradientBuilder<float>* msccb2 = 
	   new mscConvergentGradientBuilder<float>(nmf, cmh, bgf2, &a);
		// msccb2->computeGradient();	   
	   
	   mscModifiedBoundaryMeshHandler* mbmh = 
		   new mscModifiedBoundaryMeshHandler(bmsh, bgf2);

	     mscConvergentGradientBuilder<float>* msccb3 = 
	   new mscConvergentGradientBuilder<float>(bmf, mbmh, bgf, &a);


	   //msccb3->computeGradient();

   }  else if (userand == 4){
	   	   printf("--using 1-pass convergent gradient on primal\n");


		  // msccb->computeGradient();

   } else if (userand == 5){
	   	   printf("--using 1-pass convergent gradient on primal on negated function\n");


		 //  msccb->computeGradient();
	 
   }
   printf("done!\n");

   // testing
   int cc[4]; for(int i=0;i<4;i++) cc[i]=0;

   int critcount = 0;
   printf("Testing:\n");
   for (CELL_INDEX_TYPE i = 0; i < bmsh->num_cells(); i++) {

	   if (bgf->get_critical(i) && bmsh->boundary_value(i) != 0) {
		   CELL_INDEX_TYPE coords[3];
		   bmsh->cellid_2_coords(i, coords);
		   //printf("(%d, %d, %d) = %d\n", 
			//   (int) coords[0], (int) coords[1], (int) coords[2], bmsh->boundary_value(i));
	   }


	   if (! bgf->get_assigned(i)) {
		   printf("ERROR: cell %d is not assigned: %d\n", i, 
			   bmsh->dimension(i));
	   }

	   if (bgf->get_critical(i)) {
		   critcount++;
		   cc[bmsh->dimension(i)]++;
	   }

	   if (bmsh->boundary_value(i) == 0) {
			
		   cellIterator it; 
		   iteratorOperator& fit = bmsh->facets(i, it);
			fit.begin(it);
		   int counter = 0;
		   while (fit.valid(it)) {
			   counter++;
			   fit.advance(it);
		   }
		   if (counter != bmsh->dimension(i) * 2) {
			   printf("ERROR, crappy facet iterator %d != %d\n",
				   counter, bmsh->dimension(i)*2);
		   }

		   cellIterator it2;
		   iteratorOperator& fit2 = bmsh->cofacets(i, it2);
		   fit2.begin(it2);
		   counter = 0;

		   while(fit2.valid(it2)) {
			   counter++;
			   fit.advance(it2);
		   }
		   if (counter != (3 - bmsh->dimension(i)) * 2) {
			   printf("ERROR, crappy cofacet iterator %d != %d\n",
				   counter, (3- bmsh->dimension(i))*2);
		   }

	   }

   }
	   printf("critcount = %d\n", critcount);

	   for (int i=0; i < 4; i++) printf("#index-%d = %d\n", i, cc[i]);
	printf("output result...\n");

   bgf->output_to_renderer(filename);
   bmf->output_to_renderer(filename);
   printf("ALL DONE!\n");

   return 1;


}
