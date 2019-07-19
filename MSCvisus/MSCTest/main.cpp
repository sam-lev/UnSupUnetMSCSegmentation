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

//#ifdef WIN32
//#include <windows.h>
//#endif
#include <map>

using namespace std;


int main(int argc, char** argv) {

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



	//FILE* fouttest = fopen("test.raw", "wb");
	//for (int i = 0; i < 3*3*3; i++) {
	//	float val = (float) test_func[i];
	//	fwrite(&val, sizeof(float), 1, fouttest);
	//}
	//fclose(fouttest);

   // declare array factory
	mscArrayFactory a(REGULAR_ARRAY);

	// load data
   mscRegularRawDataHandler<float>* test_data;
   test_data = new mscRegularRawDataHandler<float>();
   test_data->load_data(filename, X*Y*Z, &a);
  // test_data->logify();
   //test_data->negate();
   //test_data->hack_cut();
   mscNegatingDataHandler<float>* ndh = 
	   new mscNegatingDataHandler<float>(test_data);

//#ifdef WIN32      
   printf("defined win32\n");
		  //DWORD globalStart = GetTickCount();
//#endif

		  printf("Computing discrete gradient: \n");


   // create mesh handler
   mscRegularGrid3DImplicitMeshHandler* bmsh = new mscRegularGrid3DImplicitMeshHandler(X,Y,Z);
   // tests!!!

   printf("Gothere\n");
   cellIterator it;
   iteratorOperator& fit = bmsh->cofacets(1, it);
   fit.begin(it);

   while (fit.valid(it)) {
	   printf("neighbor=%llu\n", fit.value(it));
	   fit.advance(it);
   }



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



   //mscTwoWay3DGradientBuilder<float>* msctwb =
	  // new mscTwoWay3DGradientBuilder<float>(bmf, bmsh, bgf, &a);


   mscConvergentGradientBuilder<float>* msccb = 
	   new mscConvergentGradientBuilder<float>(bmf, bmsh, bgf, &a);
      mscShortcuttingConvergentGradientBuilder<float>* mscSHcb = 
	   new mscShortcuttingConvergentGradientBuilder<float>(bmf, bmsh, bgf, &a);

   if (userand == 0) {
		printf("using randomized gradient\n");
	   mscrb->computeGradient();
   } else if (userand == 1) {
	   printf("using greedy gradient\n");
		mscb->computeGradient();
   } else if (userand == 2) {
	   printf("using convergent gradient\n");
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

#ifdef FANCY_PROB_OUTPUT

		 char ufilename1[1024];
		 sprintf(ufilename1, "%s.asc", filename); 
		 test_data->dump_vals(ufilename1, X, Y, Z, msccb->getmaxvals());
		 sprintf(ufilename1, "%s.asc.prob", filename);

		 char ufilename2[1024];

		 sprintf(ufilename2, "%s.dsc", filename);
		 test_data->dump_vals(ufilename2, X, Y, Z, msccb2->getmaxvals());
   		 sprintf(ufilename2, "%s.dsc.prob", filename);

		 char ufilename3[1024];

		 sprintf(ufilename3, "%s.both", filename);
		 FILE* ff1 = fopen(ufilename1, "rb");
		 FILE* ff2 = fopen(ufilename2, "rb");
		 FILE* ffout = fopen(ufilename3, "wb");
		 while (! feof(ff1)) {
			 float v1, v2;
			 fread(&v1, sizeof(float), 1, ff1);
			 fread(&v2, sizeof(float), 1, ff2);
			 float res = (0.5-0.5f*v1) + 0.5f - (0.5-0.5f*v2);
			 fwrite(&res, sizeof(float), 1, ffout);
		 }
		 fclose(ff1);
		 fclose(ff2);
		 fclose(ffout);
#endif


   }else if (userand == 5) {
	   printf("using NEW convergent gradient\n");

#ifdef DEBUG_SIZE_OF_MD
	   mscSHcb->computeGradient(0);
#else
	   mscSHcb->computeGradient();
#endif
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
#ifdef DEBUG_SIZE_OF_MD
	   		 msccb2->computeGradient(1);

#else
	 		 msccb2->computeGradient();
#endif
		  mscRegularGrid3DGradientField* tempgf = bgf;
		  
		  cellIterator it;
		  iteratorOperator& all_cells = bmsh->all_cells_iterator(it);
		  for (all_cells.begin(it); all_cells.valid(it); all_cells.advance(it)) {
			  CELL_INDEX_TYPE cid = all_cells.value(it);
			  bgf2->set_dim_asc_man(cid, bgf->get_dim_asc_man(cid));
		  }
		  bgf = bgf2;

		  bgf->resetMeshHandler(bmsh);

#ifdef FANCY_PROB_OUTPUT

		 char ufilename1[1024];
		 sprintf(ufilename1, "%s.asc", filename); 
		 test_data->dump_vals(ufilename1, X, Y, Z, msccb->getmaxvals());
		 sprintf(ufilename1, "%s.asc.prob", filename);

		 char ufilename2[1024];

		 sprintf(ufilename2, "%s.dsc", filename);
		 test_data->dump_vals(ufilename2, X, Y, Z, msccb2->getmaxvals());
   		 sprintf(ufilename2, "%s.dsc.prob", filename);

		 char ufilename3[1024];

		 sprintf(ufilename3, "%s.both", filename);
		 FILE* ff1 = fopen(ufilename1, "rb");
		 FILE* ff2 = fopen(ufilename2, "rb");
		 FILE* ffout = fopen(ufilename3, "wb");
		 while (! feof(ff1)) {
			 float v1, v2;
			 fread(&v1, sizeof(float), 1, ff1);
			 fread(&v2, sizeof(float), 1, ff2);
			 float res = (0.5-0.5f*v1) + 0.5f - (0.5-0.5f*v2);
			 fwrite(&res, sizeof(float), 1, ffout);
		 }
		 fclose(ff1);
		 fclose(ff2);
		 fclose(ffout);
#endif

   } else if (userand == 3){
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
	
		  //mscRegularGrid3DGradientField* tempgf = bgf;
		  //
		  //cellIterator it;
		  //iteratorOperator& all_cells = bmsh->all_cells_iterator(it);
		  //for (all_cells.begin(it); all_cells.valid(it); all_cells.advance(it)) {
			 // CELL_INDEX_TYPE cid = all_cells.value(it);
			 // bgf2->set_dim_asc_man(cid, bgf->get_dim_asc_man(cid));
		  //}
		  //bgf = bgf2;

		  //bgf->resetMeshHandler(bmsh);



   }  else if (userand == 4){
	   	   printf("using convergent gradient -1 pass\n");


		   msccb->computeGradient();
	
		  //mscRegularGrid3DGradientField* tempgf = bgf;
		  //
		  //cellIterator it;
		  //iteratorOperator& all_cells = bmsh->all_cells_iterator(it);
		  //for (all_cells.begin(it); all_cells.valid(it); all_cells.advance(it)) {
			 // CELL_INDEX_TYPE cid = all_cells.value(it);
			 // bgf2->set_dim_asc_man(cid, bgf->get_dim_asc_man(cid));
		  //}
		  //bgf = bgf2;

		  //bgf->resetMeshHandler(bmsh);

#ifdef FANCY_PROB_OUTPUT
		   char ufilename1[1024];
		 sprintf(ufilename1, "%s.asc", filename); 
		 test_data->dump_vals(ufilename1, X, Y, Z, msccb->getmaxvals());
		 sprintf(ufilename1, "%s.asc.prob", filename);
#endif

   } /*else if (userand == 5){
	   	   printf("using convergent gradient -1 pass-reverse\n");
	   mscComplementMeshHandler* cmh = 
		   new mscComplementMeshHandler(bmsh);

	   mscNegatingMeshFunction<float>* nmf = 
		   new mscNegatingMeshFunction<float>(bmf);
	   mscRegularGrid3DGradientField* bgf2 = 
		   new mscRegularGrid3DGradientField(cmh, &a);	   

	     mscConvergentGradientBuilder<float>* msccb2 = 
	   new mscConvergentGradientBuilder<float>(nmf, cmh, bgf2, &a);
		 msccb2->computeGradient();	   

		 msccb = msccb2;
		 bgf = bgf2;



   }*/

//#ifdef WIN32        
		  //DWORD globalEnd = GetTickCount();
		  //printf(" --Computed discrete gradient in %.3f seconds\n", (globalEnd - globalStart)/1000.0);
//#endif


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


	//mscSimpleGradientUsingAlgorithms<float>* alg2 = 
	//		new mscSimpleGradientUsingAlgorithms<float>(bmf,bmsh, bgf, &a);
	//	alg2->count_critical_points(4);

   bgf->output_to_renderer(filename);
   bmf->output_to_renderer(filename);
   //for (CELL_INDEX_TYPE i = 0; i < bmsh->num_cells(); i++) {
	  // if (! bgf->get_assigned(i)) printf("ERROr %d not assigned\n", i);
	  // if (bgf->get_critical(i)) {
		 //  printf("critical %d-cell: id=%d, val=%f\n", bmsh->dimension(i), 
			//   i, bmf->cell_value(i));
	  // }
   //}

 //  map<CELL_INDEX_TYPE, int> pair_2_offset;
 //  pair_2_offset[1] = 0;
 //  pair_2_offset[-1] = 1;
 //  pair_2_offset[5] = 2;
 //  pair_2_offset[-5] = 3;
 //  pair_2_offset[25] = 4;
 //  pair_2_offset[-25] = 5;

 //  const char* dirs[7] = { "+X", "-X", "+Y", "-Y", "+Z", "-Z" , "CC" };
 //  int count = 0;
 //  for (int i = 0; i < 5; i++) {
	//   printf("\n");
	//   for (int j = 0; j < 5; j++) { 
	//	   printf("\n");
	//	   for (int k = 0; k  < 5; k++){
	//		   CELL_INDEX_TYPE asdf = bgf->get_pair(count) - count;
	//		   printf("(%d,%s) ", 
	//			   //bgf->get_assigned(count),
	//			   //bgf->get_mark(count),
	//			   bgf->get_dim_asc_man(count),
	//			   //bgf->get_pair(count),
	//			   //bgf->get_critical(count),
	//			   dirs[(bgf->get_critical(count) ? 6 : pair_2_offset[asdf])]);
	//		   count++;
	//	   }
	//	}
	//}
 //  printf("\n\n");

 /*     mscSimpleGradientUsingAlgorithms<float>* msca = 
	   new mscSimpleGradientUsingAlgorithms<float>(bmf, bmsh, bgf, &a);

	  for (int i = 0; i < 5*5*5; i++) {
		  if (bgf->get_critical(i)) {
			  printf("descMan of %d:\n", i);
			  vector<CELL_INDEX_TYPE> dsc_man;
			  vector<int> counts;
			  msca->trace_down_cells_restricted_counting(i, dsc_man, counts);
			  for (int j = 0; j < dsc_man.size(); j++) {
				  if(bgf->get_critical(dsc_man[j])) {
					   printf("(%d)%d ", dsc_man[j], counts[j]);
				  } else {
					   printf("%d.%d ", dsc_man[j], counts[j]);
				  }
			  }
			  printf("\n\n");
		  }
	  }*/


 //  cellIterator dit;
 //  iteratorOperator& dito = bmsh->d_cells_iterator(0, dit);
 //  dito.begin(dit);
	//CELL_INDEX_TYPE count = 0;
 //  while(dito.valid(dit)) {
	//   CELL_INDEX_TYPE i = dito.value(dit);
	//   printf("%d-cell[%d] = %f, count=%d\n", bmsh->dimension(i),
	//	   i, bmf->cell_value(i), count++);
	//   dito.advance(dit);
 //  }

 //  cellIterator it;
 //  //iteratorOperator& ito = bmsh->all_cells_iterator(it);
 //  //iteratorOperator& ito = bmsh->d_cells_iterator(2, it);
 //  
 //  CELL_INDEX_TYPE temp_id = (1) + (1) * 5 + (1) * 5 * 5;

 //  printf("Value of %d-cell %d = %f\n-->\n", 
	//   bmsh->dimension(temp_id), temp_id, bmf->cell_value(temp_id));

 //  iteratorOperator& ito = bmsh->facets(temp_id, it);
 //  ito.begin(it);
 //  while (ito.valid(it)) {
	//   CELL_INDEX_TYPE i = ito.value(it);
	//   printf("val of %d-cell %d = %f\n",
	//	   bmsh->dimension(i), i, bmf->cell_value(i));
	//   //CELL_INDEX_TYPE& val = *it;
	//   //printf("%d\n ", (int) val);
	//   ito.advance(it);
 //  }
 //  printf("\n");
   

   //A a(1, &num);
   //A b(2, &num);
   //A c= a;
   //a.printme();
   //b.printme();
   //c.printme();

   //
   //a.poop = 10;
   //a.fart = &num;
   //b.fart = &num;

   //a.printme();
   //b.printme();
   //c.printme();

   //printf("sizeofA = %d\n", sizeof(A));
//   mscRegularArray<int> stuff(num);
//  for (int i = 0; i < 10; i++) {
//    stuff[i*2] = i;
//}
//   printf("got here\n");
//
//   for (int i = 0; i < 20; i++) {
//      printf("val[%d] = %d\n", i, stuff[i]);
//   }


   return 1;


}
