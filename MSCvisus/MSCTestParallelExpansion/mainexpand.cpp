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

	if (argc < 5) printf("usage: filename X Y Z nThreads\n");

	for (int i =0; i < argc; i++) {
		printf("%d=%s\n", i, argv[i]);
	}

	char* filename = argv[1];
	int X, Y, Z, nbX, nbY, nbZ;
	sscanf(argv[2], "%d", &X);
	sscanf(argv[3], "%d", &Y);
	sscanf(argv[4], "%d", &Z);





	mscArrayFactory a(REGULAR_ARRAY);

	// load data
	mscRegularRawDataHandler<float>* test_data;
	test_data = new mscRegularRawDataHandler<float>();
	test_data->load_data(filename, X*Y*Z, &a);
	  //est_data->logify();
	  // create mesh handler
	  mscRegularGrid3DImplicitMeshHandler* bmsh = new mscRegularGrid3DImplicitMeshHandler(X,Y,Z);
	  printf("POOH %d %d %d\n\n", X, Y, Z);

	  // create mesh function
	  mscRegularGrid3DMeshFunction<float>* bmf = new mscRegularGrid3DMeshFunction<float>(test_data, bmsh, &a);
	  bmf->initialize();



	return 1;


}
