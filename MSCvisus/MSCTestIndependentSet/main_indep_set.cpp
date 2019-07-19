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
#include "mscBasicMeshDecomposition.h"
#include "mscRegularGrid3DDecomposition.h"

#include <windows.h>
#include <map>

using namespace std;

void makeTest() {
	FILE* f = fopen("test.raw", "wb");
	for (int i; i < 100*100*100; i++) {
		float t = (float) i;
		fwrite(&t, sizeof(float), 1, f);
	}
	fclose(f);
}


int main(int argc, char** argv) {


	//makeTest();

	if (argc < 5) printf("usage: filename X Y Z nbX nbY nbZ nThreads\n");

	for (int i =0; i < argc; i++) {
		printf("%d=%s\n", i, argv[i]);
	}

	char* filename = argv[1];
	int X, Y, Z, nbX, nbY, nbZ;
	sscanf(argv[2], "%d", &X);
	sscanf(argv[3], "%d", &Y);
	sscanf(argv[4], "%d", &Z);
	sscanf(argv[5], "%d", &nbX);
	sscanf(argv[6], "%d", &nbY);
	sscanf(argv[7], "%d", &nbZ);




    // declare array factory
	mscArrayFactory a(REGULAR_ARRAY);
	// load data
   mscRegularRawDataHandler<float>* test_data;
   test_data = new mscRegularRawDataHandler<float>();
   test_data->load_data(filename, X*Y*Z, &a);
   // create mesh handler
   mscRegularGrid3DImplicitMeshHandler* bmsh = new mscRegularGrid3DImplicitMeshHandler(X,Y,Z);
   // create mesh function
   mscRegularGrid3DMeshFunction<float>* bmf = new mscRegularGrid3DMeshFunction<float>(test_data, bmsh, &a);
   bmf->initialize();


   

   // init fronteir with all minima?


   std::queue<CELL_INDEX_TYPE> fronteir;

	//char grad_name[1024];
	//sprintf(grad_name, "%s.grad", filename);
	//char dat_name[1024];
	//sprintf(dat_name, "%s.dat", filename);

	//// create a file of the right size
	//FILE* fgrad = fopen(grad_name, "wb");
	//FILE* fdat = fopen(dat_name, "wb");
	//CELL_INDEX_TYPE num_in_grad = 
	//	(X*2-1) * (Y*2-1) * (Z*2-1);
	//char zero = 0;
	//float val = 0.0f;
	//for (int i = 0; i < num_in_grad; i++) {
	//	fwrite(&zero, sizeof(char), 1, fgrad);
	//	fwrite(&val, sizeof(float), 1, fdat);
	//}
	//fclose(fgrad);
	//fclose(fdat);


	//////////mscArrayFactory a(REGULAR_ARRAY);

	//////////// load data
	//////////mscRegularRawDataHandler<float>* test_data;
	//////////test_data = new mscRegularRawDataHandler<float>();
	//////////test_data->load_data(filename, X*Y*Z, &a);
	//////////  //est_data->logify();
	//////////  // create mesh handler
	//////////  mscRegularGrid3DImplicitMeshHandler* bmsh = new mscRegularGrid3DImplicitMeshHandler(X,Y,Z);
	//////////  printf("POOH %d %d %d\n\n", X, Y, Z);

	//////////  // create mesh function
	//////////  mscRegularGrid3DMeshFunction<float>* bmf = new mscRegularGrid3DMeshFunction<float>(test_data, bmsh, &a);
	//////////  bmf->initialize();

	//////////   mscRegularGrid3DGradientField* bgf = 
	//////////   new mscRegularGrid3DGradientField(bmsh, &a);

	//////////  mscSimpleGradientBuilder<float>* mscb = 
	//////////   new mscSimpleGradientBuilder<float>(bmf, bmsh, bgf, &a);
	//////////  mscb->computeGradient();

	mscRegularGrid3DDecomposition<float>* mscD = 
		new mscRegularGrid3DDecomposition<float>(X, Y, Z, nbX, nbY, nbZ);







	//////////mscSimpleGradientUsingAlgorithms<float>* alg2 = 
	//////////		new mscSimpleGradientUsingAlgorithms<float>(bmf,bmsh, bgf, &a);
	//////////	alg2->count_critical_points(4);


	mscArrayFactory* array_factory = new mscArrayFactory();
	mscD->setInputFileName(filename);
	//mscD->setOutputFileNames(grad_name, dat_name);
	mscD->decompose();

	for (int i = 0; i < mscD->numBlocks(); i++) {

		printf("DOING BLOCK:\n");
		mscD->printBlockInfo(i);
		mscD->loadBlock(i);

		//		for (int j = 0; j < mscD->getMeshHandler(i)->num_cells(); j++) {
		//	if (bgf->get_pair(j) != mscD->getGradientField(i)->get_pair(j))
		//		printf("%d != %d \n", bgf->get_pair(j), mscD->getGradientField(i)->get_pair(j));
		//}

		mscSimpleGradientBuilder<float>* gradient_builder = 
			new mscSimpleGradientBuilder<float>(mscD->getMeshFunction(i), 
			mscD->getMeshHandler(i), mscD->getGradientField(i), array_factory);



		gradient_builder->computeGradient();


		//for (int j = 0; j < mscD->getMeshHandler(i)->num_cells(); j++) {
		//	if (bgf->get_pair(j) != mscD->getGradientField(i)->get_pair(j))
		//		printf("%d != %d \n", bgf->get_pair(j), mscD->getGradientField(i)->get_pair(j));
		//}
		mscSimpleGradientUsingAlgorithms<float>* alg = 
			new mscSimpleGradientUsingAlgorithms<float>(mscD->getMeshFunction(i), 
			mscD->getMeshHandler(i), mscD->getGradientField(i), array_factory);
		alg->count_critical_points(4);
		//delete gradient_builder;


		////bmsh = (mscRegularGrid3DImplicitMeshHandler*) mscD->getMeshHandler(i);
		////bgf = (mscRegularGrid3DGradientField*) mscD->getGradientField(i);

		//// // testing
		//// int critcount = 0;
		//// printf("Testing:\n");
		//// for (CELL_INDEX_TYPE i = 0; i < bmsh->num_cells(); i++) {
		////  
		////  
		////  

		////  
		////  
		////  
		////  if (! bgf->get_assigned(i)) {
		////   printf("ERROR: cell %d is not assigned: %d\n", i, 
		////	   bmsh->dimension(i));
		////  }

		////  if (bgf->get_critical(i)) critcount++;

		////  if (bmsh->boundary_value(i) == 0) {
		////	
		////   cellIterator it; 
		////   iteratorOperator& fit = bmsh->facets(i, it);
		////	fit.begin(it);
		////   int counter = 0;
		////   while (fit.valid(it)) {
		////	   counter++;
		////	   fit.advance(it);
		////   }
		////   if (counter != bmsh->dimension(i) * 2) {
		////	   printf("ERROR, crappy facet iterator %d != %d\n",
		////		   counter, bmsh->dimension(i)*2);
		////   }

		////   cellIterator it2;
		////   iteratorOperator& fit2 = bmsh->cofacets(i, it2);
		////   fit2.begin(it2);
		////   counter = 0;

		////   while(fit2.valid(it2)) {
		////	   counter++;
		////	   fit.advance(it2);
		////   }
		////   if (counter != (3 - bmsh->dimension(i)) * 2) {
		////	   printf("ERROR, crappy cofacet iterator %d != %d\n",
		////		   counter, (3- bmsh->dimension(i))*2);
		////   }

		////  }




		//printf("critcount = %d\n", critcount);




























		// now write the gradient
		printf("about to write\n");

		mscD->writeOutputs(i);
		mscD->unloadBlock(i);
		delete gradient_builder;
		delete alg;

	}

	for (int i= 0; i < mscD->numBlocks(); i++)
		mscD->printBlockInfo(i);
	//SEQUENCE OF OPERATIONS:

	// decompose domain -> produces a list of blocks


	// for each block in domain {

	//   load the block -> create mscBasicDataHandler, read in data
	//					-> create mscBasicMeshFunction
	//					-> create mscBasicMeshHandler
	//					-> create mscBasicGradientField

	//   compute local gradient usinb block->getFUnction,getMesh,getGrad

	//   write out gradient

	//   delete block?

	// }



	////	//FILE* fouttest = fopen("test.raw", "wb");
	////	//for (int i = 0; i < 3*3*3; i++) {
	////	//	float val = (float) test_func[i];
	////	//	fwrite(&val, sizeof(float), 1, fouttest);
	////	//}
	////	//fclose(fouttest);
	////
	////   // declare array factory
	////	mscArrayFactory a(REGULAR_ARRAY);
	////
	////	// load data
	////   mscRegularRawDataHandler<float>* test_data;
	////   test_data = new mscRegularRawDataHandler<float>();
	////   test_data->load_data(filename, X*Y*Z, &a);
	////   test_data->logify();
	////
	////#ifdef WIN32        
	////		  DWORD globalStart = GetTickCount();
	////#endif
	////
	////		  printf("Computing discrete gradient: \n");
	////
	////
	////   // create mesh handler
	////   mscRegularGrid3DImplicitMeshHandler* bmsh = new mscRegularGrid3DImplicitMeshHandler(X,Y,Z);
	////
	////   // create mesh function
	////   mscRegularGrid3DMeshFunction<float>* bmf = new mscRegularGrid3DMeshFunction<float>(test_data, bmsh, &a);
	////   bmf->initialize();
	////
	////
	////
	////   mscRegularGrid3DGradientField* bgf = 
	////	   new mscRegularGrid3DGradientField(bmsh, &a);
	////
	////   // test the gradient builder!
	////   mscSimpleGradientBuilder<float>* mscb = 
	////	   new mscSimpleGradientBuilder<float>(bmf, bmsh, bgf, &a);
	////
	////   mscSimpleRandomGradientBuilder<float>* mscrb = 
	////	   new mscSimpleRandomGradientBuilder<float>(bmf, bmsh, bgf, &a);
	////
	////   if (userand) {
	////		printf("using randomized gradient\n");
	////	   mscrb->computeGradient();
	////   } else{
	////	   printf("using greedy gradient\n");
	////		mscb->computeGradient();
	////   }
	////
	////#ifdef WIN32        
	////		  DWORD globalEnd = GetTickCount();
	////		  printf(" --Computed discrete gradient in %.3f seconds\n", (globalEnd - globalStart)/1000.0);
	////#endif
	////
	////   bgf->output_to_renderer(filename);
	////   bmf->output_to_renderer(filename);
	////


	return 1;


}
