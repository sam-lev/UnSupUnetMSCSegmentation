


#include <iostream>
#include "dump_image.h"
#ifdef WIN32
#include <windows.h>
#endif


#include <stdio.h>
#include <stdlib.h>
#include <OpenGL/glu.h>  // openGL utilities
#include <OpenGL/gl.h>   // openGL declarations
#ifdef WIN32
#include "freeglut.h"
#else
#include "glut.h"
#endif

#include "Matrix4.h"
#ifdef WIN32
//#include "C:\local\CImg-1.3.2\CImg.h"
#include "C:\Users\jediati\Desktop\JEDIATI\libraries\CImg-1.5.7\CImg.h"
using namespace cimg_library;
#endif

#ifndef WIN32
#include <unistd.h>
#endif

#include <cmath>
#include <algorithm>



#include "mscArrayFactory.h"
#include "mscDumbGradientField.h"
#include "mscRegularRawDataHandler.h"
#include "mscRegularGrid3DImplicitMeshHandler.h"
#include "mscRegularGrid2DImplicitMeshHandler.h"
#include "mscDumbStoringMeshFunction.h"
#include "mscDumbStoringMinFunction.h"
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
#include "mscBasicMSC.h"
#include "mscNegatingDataHandler.h"
#include <map>
#include "mscTopoArray.h"
#include "mscAssistedGradientBuilder.h"
#include "mscSimpleConstrainedGradientBuilder.h"



#define EPSILON 0.001



bool USE_SEG = false;;
bool USE_SEG_LINES = false;
bool USE_SEG_AMAP = false;
bool USE_SEG_DMAP = false;

int gX, gY, gZ;
BasicMSC<float>* MSC;
mscBasicGradientField* G_mscg;
mscBasicGradientField* G_mscg_TEMP = NULL;
mscBasicMeshHandler* G_mscmh;
mscBasicMeshHandler* G_mscmh2;
mscBasicMeshFunction<float>* G_mscf;
mscConvergentGradientBuilder<float>* G_msccb;
mscConvergentGradientBuilder<float>* G_msccb2;

mscPreClassifier* classes;
float red_scale(float s) {
	Vector4 v; 
	v.vals[0] = s; //max(s, 1.0f-s);
	v.vals[1] = 1.0f-s;
	v.vals[2] =  max(s, 1.0f-s);
	Normalize3(&v);
	return v.vals[0];
	//return min(max(4.0*(0.5-s), 0.0), 1.0);
}

float green_scale(float s) {
	Vector4 v; 
	v.vals[0] = s;//max(s, 1.0f-s);
	v.vals[1] = 1.0f-s;
	v.vals[2] = max(s, 1.0f-s);
	Normalize3(&v);
	return v.vals[1];
	//return 0.2f;
	//return min(max(4.0*fabs(s-0.25)-1.0, 0.0), 1.0);
}

float blue_scale(float s) {
	Vector4 v; 
	v.vals[0] = s;//max(s, 1.0f-s);
	v.vals[1] = 1.0f-s;
	v.vals[2] = max(s, 1.0f-s);
	Normalize3(&v);
	return v.vals[2];
	// return  1.0-s;
	//return min(max(4.0*(0.75-s), 0.0), 1.0);
}


using namespace std;

void Initialize_GLUT(int argc, char **argv);


void draw_earth ( float x, float y, float z, float size )
{
	glEnable(GL_LIGHTING);
	glPushMatrix ( );
	glTranslatef ( x, y, z );
	//glScalef(size, size, size);
	GLUquadricObj* q = gluNewQuadric  ( );
	gluQuadricDrawStyle ( q, GLU_FILL   );
	gluQuadricNormals   ( q, GLU_SMOOTH );
	gluSphere ( q, size, 10, 10 );
	gluDeleteQuadric ( q );
	glPopMatrix ( );
	glDisable(GL_LIGHTING);
}


int seed;
int XMIN = 2;
int YMIN = 2;
int ZMIN = 2;


int XMAX = 2;
int YMAX = 2;
int ZMAX = 2;

float ballsize = 0.1f;
float percp = 0.0f;

template<class FType> struct gminmax {
	FType minval;
	FType maxval;
	float xmin;
	float xmax;
	float ymin;
	float ymax;

};

bool gusecutoffhack = false;
bool draw_flat = false;
bool draw_edges = true;
gminmax<float> fminmax;

bool clearcolor = 1;
bool g_draw_isolines = false;


bool draw_gradient = false;
float arrow_width = 0.1;
float line_width = 2.0;
float point_size = 4.0;
//template <unsigned char Dim, class FType>
//void centroid(UnstructuredSimplicialComplex<Dim, FType>* bcc, index_type cellid, float* verts) {
//
//    for (int i = 0; i < 3; i++) verts[i] = 0.0f;
//    Simplex<FType>& s = bcc->cells[cellid];
//    for (int i = 0; i < s.numberOfVerts; i++) {
//      BaseVertex<Dim>& v = bcc->verts[s.verts[i]];
//      verts[0] += v.position[0];
//      if (! draw_flat)
//	verts[1] += bcc->getValue(s.verts[i]);
//      else 
//	verts[1] += bcc->getValue(cellid);
//      verts[2] += v.position[1];
//    }
//    for (int i = 0; i < 3; i++) verts[i] /= s.numberOfVerts;
//
//};

void coordinates(CELL_INDEX_TYPE cellid, float* c){
	c[0] = (float) (cellid % XMAX);
	c[1] = (float) ((cellid / XMAX) % YMAX);
	c[2] = (float) (cellid / (XMAX*YMAX));
}

void coordinates2(int cellid, float* c){
	c[0] = (float) (cellid % gX) * 2;
	c[1] = (float) ((cellid / gX) % gY) * 2;
	c[2] = (float) (cellid / (gX*gY)) * 2;
}

void drawArrow( CELL_INDEX_TYPE cellid) {
	if (G_mscg->get_critical(cellid)) return;
	if (! G_mscg->get_assigned(cellid)) return;
	CELL_INDEX_TYPE pair = G_mscg->get_pair(cellid);
	if (G_mscmh->dimension(pair) > G_mscmh->dimension(cellid)) return;
	float start[3];
	float end[3];

	coordinates(pair, start);
	coordinates(cellid, end);

	//instead of a line
	//  glVertex3f(start[0], start[1]+ 2*EPSILON, start[2]);
	//  glVertex3f(end[0], end[1]+ 2*EPSILON, end[2]);


	//let's use code for computing the worlds most expensive arrow.
	double diff[3];
	diff[0] = end[0]-start[0];
	diff[1] = end[1]-start[1];
	diff[2] = end[2]-start[2];
	double len = diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2];
	len = sqrt(len);
	diff[0] /= len;
	diff[1] /= len;
	diff[2] /= len;

	const float RADDEG = 57.29578;

	float Q = atan2 (diff[1], diff[0]) * RADDEG;
	float P = acos  (diff[2]) * RADDEG;

	glColor3f(0.15,0.15,0.15);

	glPushMatrix();
	GLUquadricObj *quadObj = gluNewQuadric();
	GLUquadricObj *quadObj2 = gluNewQuadric();
	GLUquadricObj *quadObj3 = gluNewQuadric();
	GLUquadricObj *quadObj4 = gluNewQuadric();

	gluQuadricDrawStyle ( quadObj, GLU_FILL   );
	gluQuadricNormals   ( quadObj, GLU_SMOOTH );
	gluQuadricDrawStyle ( quadObj2, GLU_FILL   );
	gluQuadricNormals   ( quadObj2, GLU_SMOOTH );
	gluQuadricDrawStyle ( quadObj3, GLU_FILL   );
	gluQuadricNormals   ( quadObj3, GLU_SMOOTH );
	gluQuadricDrawStyle ( quadObj4, GLU_FILL   );
	gluQuadricNormals   ( quadObj4, GLU_SMOOTH );


	glTranslatef(start[0],start[1],start[2]);
	glRotatef (Q,0,0,1);
	glRotatef (P,0,1,0);
	gluCylinder(quadObj,0.4*arrow_width,0.4*arrow_width,len-2*arrow_width,10,10);
	gluDisk(quadObj2,0,0.4*arrow_width,10,10);

	glTranslatef(0,0,len-2.5*arrow_width);
	gluCylinder(quadObj3,arrow_width,0,2.5*arrow_width,10,10);
	gluDisk(quadObj4,0,arrow_width,10,10);

	gluDeleteQuadric(quadObj);
	gluDeleteQuadric(quadObj2);
	gluDeleteQuadric(quadObj3);
	gluDeleteQuadric(quadObj4);
	glPopMatrix();



};

bool DRAW_BACKGROUND = true;
bool DRAW_BACKGROUND2 = false;
bool DRAW_BACKGROUND3 = false;
bool DRAWASCLINES = false;
bool DRAWDSCLINES = false;
bool DRAWPOINTS = true;
bool DRAWDIMASCMAN = false;
bool DRAWGRAD = false;
bool DRAW_ISOLINES = true;
bool DRAW_PROBABILITIES1 = false;
bool DRAW_PROBABILITIES2 = false;
bool DRAW_GRID = false;
bool draw_mt = false;
bool flat_funct = false;
bool redrawstuff = true;
GLuint drawlist = -1;


bool DRAWNUMERICALDE = false;
bool DRAWNUMERICALDE_INIT = false;
int* NUMDE_sources;
int* NUMDE_dests;


struct DE {
	node<float>* nup;
	node<float>* ndown;
	int size;
	int interior;
	int boundary;
};

float distance(CELL_INDEX_TYPE id1, CELL_INDEX_TYPE id2) {
	float c1[3]; 
	coordinates(id1, c1);
	float c2[3]; 
	coordinates(id2, c2);
	return sqrt((c1[0]-c2[0])*(c1[0]-c2[0]) +(c1[1]-c2[1])*(c1[1]-c2[1])); 
};


void DrawDe() {
	CELL_INDEX_TYPE* aid= new CELL_INDEX_TYPE[G_mscmh->num_cells()]();
	CELL_INDEX_TYPE* did= new CELL_INDEX_TYPE[G_mscmh->num_cells()]();


	// get ascending man dimensions
	map< pair<CELL_INDEX_TYPE, CELL_INDEX_TYPE>, DE > pairmap;

	glBegin(GL_LINES);
	map<CELL_INDEX_TYPE, node<float>*>::iterator nit = MSC->nodes.begin();
	while (nit != MSC->nodes.end()) {
		node<float>* n = (*nit).second;
		nit++;
		if (! MSC->isAlive(n)) continue;
		if (n->index != 0) continue;

		set<CELL_INDEX_TYPE> g;
		MSC->fillGeometry(n, g);

		for (set<CELL_INDEX_TYPE>::iterator geom = g.begin(); geom != g.end(); geom++) {
			CELL_INDEX_TYPE id0 = *geom;
			cellIterator it1;
			iteratorOperator& cofacets1 = G_mscmh->cofacets(id0, it1);
			for (cofacets1.begin(it1); cofacets1.valid(it1); cofacets1.advance(it1)) {
				CELL_INDEX_TYPE id1 = cofacets1.value(it1);
				g.insert(id1);
				cellIterator it2;
				iteratorOperator& cofacets2 = G_mscmh->cofacets(id1, it2);
				for (cofacets2.begin(it2); cofacets2.valid(it2); cofacets2.advance(it2)) {
					CELL_INDEX_TYPE id2 = cofacets2.value(it2);
					g.insert(id2);
				}
			}
		}
		// now g has ids of all cells in thingy.
		for(set<CELL_INDEX_TYPE>::iterator sit = g.begin(); sit != g.end(); sit++) {
			aid[*sit] = n->cellid;
		}

		set<CELL_INDEX_TYPE> maxs;
		arc<float>* a = n->firstarc;
		while (a != NULL) {
			if (! MSC->isAlive(a)) {
				a = a->lower_next;
				continue;
			}
			node<float>* n1 = a->upper;
			arc<float>* a1 = n1->firstarc;
			while (a1 != NULL) {
				if (! MSC->isAlive(a1) || a1->lower != n1) {
					if (a1->lower == n1) {
						a1 = a1->lower_next;
					} else {
						a1 = a1->upper_next;
					}
					continue;
				}
				maxs.insert(a1->upper->cellid);
				a1 = a1->lower_next;
			}
			a = a->lower_next;
		}
		// now maxs has all maxs attached to this min
		float c0[3];
		coordinates(n->cellid, c0);
		for (set<CELL_INDEX_TYPE>::iterator sit = maxs.begin(); sit != maxs.end(); sit++) {
			float c1[3]; 
			coordinates((*sit), c1);
			glColor4f(.5,0,1,.5);
			glVertex3f(c0[0], c0[1], c0[2]);
			glColor4f(1,0,.5,.5);
			glVertex3f(c1[0], c1[1], c1[2]);

			//also insert into pairmap
			DE t; t.ndown = n; t.nup = MSC->nodes[(*sit)]; t.size = 0; t.interior = 0;
			t.boundary = 0;
			pairmap[pair<CELL_INDEX_TYPE, CELL_INDEX_TYPE>(n->cellid, (*sit))] = t;
		}
	}
	glEnd();

	nit = MSC->nodes.begin();
	while (nit != MSC->nodes.end()) {
		node<float>* n = (*nit).second;
		nit++;
		if (! MSC->isAlive(n)) continue;
		if (n->index != 2) continue;

		set<CELL_INDEX_TYPE> g;
		MSC->fillGeometry(n, g);
		for (set<CELL_INDEX_TYPE>::iterator geom = g.begin(); geom != g.end(); geom++) {
			CELL_INDEX_TYPE id0 = *geom;
			cellIterator it1;
			iteratorOperator& facets1 = G_mscmh->facets(id0, it1);
			for (facets1.begin(it1); facets1.valid(it1); facets1.advance(it1)) {
				CELL_INDEX_TYPE id1 = facets1.value(it1);
				g.insert(id1);
				cellIterator it2;
				iteratorOperator& facets2 = G_mscmh->facets(id1, it2);
				for (facets2.begin(it2); facets2.valid(it2); facets2.advance(it2)) {
					CELL_INDEX_TYPE id2 = facets2.value(it2);
					g.insert(id2);
				}
			}
		}
		// now g has ids of all cells in thingy.
		for(set<CELL_INDEX_TYPE>::iterator sit = g.begin(); sit != g.end(); sit++) {
			did[*sit] = n->cellid;
		}
	}
	//now aid and did have ids of ascending and descending manifolds

	for (CELL_INDEX_TYPE i = 0; i < G_mscmh->num_cells(); i++) {
		CELL_INDEX_TYPE id1 = aid[i];
		CELL_INDEX_TYPE id2 = did[i];
		//printf("%d, %d\n", (int) id1, (int) id2);
		pair<CELL_INDEX_TYPE, CELL_INDEX_TYPE> p(id1, id2);
		if (pairmap.count(p) != 0) {
			DE& d = pairmap[p];
			d.size++;

			bool isInterior = true;
			cellIterator fit;
			iteratorOperator& facets = G_mscmh->facets(i, fit);
			for (facets.begin(fit); facets.valid(fit); facets.advance(fit)) {
				CELL_INDEX_TYPE nid = facets.value(fit);
				if (aid[nid] != id1 || did[nid] != id2) isInterior = false;
			}
			iteratorOperator& cofacets = G_mscmh->cofacets(i, fit);
			for (cofacets.begin(fit); cofacets.valid(fit); cofacets.advance(fit)) {
				CELL_INDEX_TYPE nid = cofacets.value(fit);
				if (aid[nid] != id1 || did[nid] != id2) isInterior = false;
			}
			if(isInterior) {
				d.interior++;
			} else {
				d.boundary++;
			}
		}
	}

	//	FILE* fout = fopen("DISSIPATIONELEMENTS.txt", "w");
	//	for (map< pair<CELL_INDEX_TYPE, CELL_INDEX_TYPE>, DE>::iterator it = pairmap.begin();
	//		it != pairmap.end(); it++) {
	//			DE& d = (*it).second;
	//
	//			float c1[3]; 
	//			coordinates(d.ndown->cellid, c1);
	//			float c2[3]; 
	//			coordinates(d.nup->cellid, c2);
	//	
	//			fprintf(fout, "%f %f %f %d %f %d %d %f %f %f %f\n", d.ndown->value, d.nup->value,
	//				d.nup->value - d.ndown->value, d.size, distance(d.nup->cellid, d.ndown->cellid),
	//				d.interior, d.boundary, c1[0], c1[1], c2[0], c2[1]);
	//	}
	//	fclose(fout);

	delete[] aid;
	delete[] did;

}

struct TCOLOR {
	float r, g, b;
	int count;
};
map<int, map<int, TCOLOR> > sd2c;
vector<float> g_num_grad_pairs;

	struct MP {
		float c[2];
	};
	vector< vector<MP> > tlines;
	vector< MP> tlines2;

	void drawSEG() {

		
		if (USE_SEG_AMAP) {
			
			glBegin(GL_QUADS);
			cellIterator it;
			iteratorOperator& all = G_mscmh->d_cells_iterator(0, it);
			for (all.begin(it); all.valid(it); all.advance(it)) {
				//if (G_mscg->get_dim_asc_man(all.value(it)) != 3) 
				CELL_INDEX_TYPE cid = all.value(it);

				int v = classes->getId(cid); 

			float col[3];
			coordinates(v, col);
#define PI 3.14159265f

			int dd = 60;
			float hd = 60.0f;

			float v1 =  (sin(PI*((int) col[0]%dd)/((float) hd)))*0.8f+.2f;
			float v2 =  (sin((dd - (((int) col[1])%dd))/hd))*0.8f+.2f;
			float v3 = (sin(PI*((((int) col[0])%dd)-((int) col[1])%dd)/hd))*0.8f+.2f;

				float R =  0.3 + (((v ) * 54247) % 167) / 255.0;
				float G =  0.3 + (((v ) * 68321) % 167) / 255.0;
				float B = 0.3 + (((v ) * 5373) % 167) / 255.0;
				glColor3f(R,G,B);
					float c[3]; 
					coordinates(cid, c);
					glVertex3f(c[0]-1, c[1]-1, c[2]);
					glVertex3f(c[0]-1, c[1]+1, c[2]);
					glVertex3f(c[0]+1, c[1]+1, c[2]);
					glVertex3f(c[0]+1, c[1]-1, c[2]);
				}
			
			glEnd();
		}

		if (USE_SEG_DMAP) {
			
			glBegin(GL_QUADS);
			cellIterator it;
			iteratorOperator& all = G_mscmh->d_cells_iterator(2, it);
			for (all.begin(it); all.valid(it); all.advance(it)) {
				//if (G_mscg->get_dim_asc_man(all.value(it)) != 3) 
				CELL_INDEX_TYPE cid = all.value(it);

				int v = classes->getId(cid); 

			float col[3];
			coordinates(v, col);
#define PI 3.14159265f

			int dd = 60;
			float hd = 60.0f;

			//float v1 = (sin(PI*((int) col[0]%dd)/((float) hd)))*0.8f+.2f;
			//float v2 = (sin((dd - (((int) col[1])%dd))/hd))*0.8f+.2f;
			//float v3 = (sin(PI*((((int) col[0])%dd)-((int) col[1])%dd)/hd))*0.8f+.2f;
			float v1 =  (sin(PI*((int) col[0]%dd)/((float) hd)))*0.8f+.2f;
			float v2 =  (sin(PI*((int) col[1]%dd)/((float) hd)))*0.8f+.2f;
			float v3 =	(sin(PI*((int) col[0]%dd)/((float) hd)))*0.8f+.2f;

				float R =  0.3 + (((v ) * 54247) % 167) / 255.0;
				float G =  0.3 + (((v ) * 68321) % 167) / 255.0;
				float B = 0.3 + (((v ) * 5373) % 167) / 255.0;
				glColor3f(R,G,B);
					float c[3]; 
					coordinates(cid, c);
					glVertex3f(c[0]-1, c[1]-1, c[2]);
					glVertex3f(c[0]-1, c[1]+1, c[2]);
					glVertex3f(c[0]+1, c[1]+1, c[2]);
					glVertex3f(c[0]+1, c[1]-1, c[2]);
				}
			
			glEnd();
		}

	////	//if (g_num_grad_pairs.size() == 0) {
	////	FILE* fin = fopen("grad_dests.raw", "rb");
	////	g_num_grad_pairs.clear();
	////	while (! feof(fin)) {
	////		float s;
	////		fread(&s, sizeof(float), 1, fin);
	////		g_num_grad_pairs.push_back( s);
	////		//printf("%d %f\n", (g_num_grad_pairs.size() -1)%4, 
	////		//	g_num_grad_pairs[g_num_grad_pairs.size()-1]);
	////	}
	////	fclose(fin);
	////	printf("read %d thingies\n", g_num_grad_pairs.size()/4);
	////	//}

	////	glBegin(GL_LINES);
	////	glColor4f(0,0,0, 0.2);
	////	for (int i = 0; i < g_num_grad_pairs.size(); i+=4) {
	////		float x = g_num_grad_pairs[i]*2;
	////		float y = g_num_grad_pairs[i+1]*2;

	////		float dx = g_num_grad_pairs[i+2]*2;
	////		float dy = g_num_grad_pairs[i+3]*2;

	////		//float n = sqrt(dx*dx + dy*dy) * 0.5;

		////		glVertex3f(x,y, 0);					
		////		glVertex3f(dx,dy, 0);					
		////		//glVertex3f(x+dx/n, y+dy/n, 0);					
		////	}
		////	glEnd();

		////}

		if(USE_SEG_LINES) {
			printf("got here\n");

			if (tlines.size() == 0) {
				FILE* fin = fopen("lines_out.raw", "rb");
				while (! feof(fin)) {
					vector<MP> line;
					int num;
					fread(&num, sizeof(int), 1, fin);
					for (int i = 0; i < num; i++) {
						MP t;
						fread(t.c, sizeof(float), 2, fin);
						line.push_back(t);
					}
					tlines.push_back(line);
				}
				fclose(fin);
			}

			glLineWidth(3.0);
			glColor4f(0,0,0,0.2);
			for (int i = 0; i < tlines.size(); i++) {

				glBegin(GL_LINE_STRIP);
				for (int j = 0; j < tlines[i].size(); j++) {
					glVertex3f(tlines[i][j].c[0]*2, tlines[i][j].c[1]*2, 0);
				}
				glEnd();
			}

			//if (tlines2.size() == 0) {
			//	FILE* fin = fopen("lines_out2.raw", "rb");
			//	if (fin != NULL) {
			//		while (! feof(fin)) {
			//			int num;
			//			fread(&num, sizeof(int), 1, fin);
			//			for (int i = 0; i < num; i++) {
			//				MP t;
			//				fread(t.c, sizeof(float), 2, fin);
			//				tlines2.push_back(t);
			//			}
			//		}
			//		fclose(fin);
			//	}
			//}

			//glLineWidth(3.0);
			//glColor4f(0.1,0.4,0.2,0.7);
			//glBegin(GL_LINES);
			//for (int i = 0; i < tlines2.size(); i++) {
			//	if(i % 2 == 0)	glColor4f(0.1,0.4,0.2,0.2);
			//	if(i % 2 == 1)	glColor4f(0.1,0.4,0.2,0.7);
			//	glVertex3f(tlines2[i].c[0]*2, tlines2[i].c[1]*2, 0);
			//}
			//glEnd();
		}


		if (USE_SEG) {
			
			glBegin(GL_QUADS);
			cellIterator it;
			iteratorOperator& all = G_mscmh->all_cells_iterator(it);
			for (all.begin(it); all.valid(it); all.advance(it)) {
				//if (G_mscg->get_dim_asc_man(all.value(it)) != 3) 
				CELL_INDEX_TYPE cid = all.value(it);

				int v = G_mscg_TEMP->get_dim_asc_man(cid); 
				if (v > 0) {
					if (v == 1) {
						glColor3f(1,0,0);
					} else if (v == 2) {
						glColor3f(0,0,1);
					} else {
						glColor3f(0,1,0);
					}

					float c[3]; 
					coordinates(cid, c);
					glVertex3f(c[0]-0.5, c[1]-0.5, c[2]);
					glVertex3f(c[0]-0.5, c[1]+0.5, c[2]);
					glVertex3f(c[0]+0.5, c[1]+0.5, c[2]);
					glVertex3f(c[0]+0.5, c[1]-0.5, c[2]);
				}
			}
			glEnd();
		}
	
	}

int* global_dsc_ids = NULL;
int* global_asc_ids = NULL;

template <unsigned char Dim, class FType>
void DumpAscMan(/*UnstructuredSimplicialComplex<Dim, FType>* bcc*/) {
	printf("dumping ids\n");
	FILE* fdsc = fopen("dump_dsc.raw", "wb");
	
	typename map<CELL_INDEX_TYPE, node<float>* >::iterator nit = MSC->nodes.begin();
	if(global_dsc_ids == NULL) {
		global_dsc_ids = new int[(gX-1) * (gY-1)];
		for (int i = 0; i < (gX-1) * (gY-1); i++) global_dsc_ids[i] = 0;
	}

	while (nit != MSC->nodes.end()) {
		node<float>* n = (*nit).second;
		nit++;
	
		if (! MSC->isAlive(n)) continue;
		if (n->index != 2) continue;
		set<CELL_INDEX_TYPE> g;
		MSC->fillGeometry(n, g);

		unsigned int value = ((n->cellid + 52341) * 347) % 253;
		if (value == -1) printf("SD:FLKJSFL:KDJSF:LDKSJ %d, %d\n", value, n->cellid);
		for (set<CELL_INDEX_TYPE>::iterator sit = g.begin(); 
			sit != g.end(); sit++) {
			CELL_INDEX_TYPE currentcid = *sit;
			int x = (currentcid % XMAX) / 2;
			int y = (currentcid / XMAX) / 2;
			global_dsc_ids[x + y * (gX-1)] = value;
		}
	}
	fwrite(global_dsc_ids, (gX-1)*(gY-1), sizeof(int), fdsc);
	fclose(fdsc);

	for (int i = 0; i < (gX-1) * (gY-1); i++) {
		//if(global_dsc_ids[i] == -1) printf("WHnnnnnnOAthere -1\n");
	}

	FILE* fasc = fopen("dump_asc.raw", "wb");

	if(global_asc_ids == NULL) {
		global_asc_ids = new int[gX * gY];
		for (int i = 0; i < gX*gY; i++) global_asc_ids[i] = 0;
	}
	 nit = MSC->nodes.begin();
	while (nit != MSC->nodes.end()) {
		node<float>* n = (*nit).second;
		nit++;
	
		if (! MSC->isAlive(n)) continue;
		if (n->index != 0) continue;
		set<CELL_INDEX_TYPE> g;
		MSC->fillGeometry(n, g);
		
		unsigned int value = ((n->cellid + 52341) * 347) % 253;

		for (set<CELL_INDEX_TYPE>::iterator sit = g.begin(); 
			sit != g.end(); sit++) {
			CELL_INDEX_TYPE currentcid = *sit;
			int x = (currentcid % XMAX) / 2;
			int y = (currentcid / XMAX) / 2;
			global_asc_ids[x + y * (gX)] = value;
			//printf("[%d, %d] = %d, %d\n", x,y,n->cellid, gX);
		}
	}
	fwrite(global_asc_ids, (gX)*(gY), sizeof(int), fasc);
	fclose(fasc);
	for (int i = 0; i < gX*gY; i++) {
		//if (global_asc_ids[i] == -1) printf("wFFFFhaoascending -1\n");
	}

}



template <unsigned char Dim, class FType>
bool ArcTest(arc<FType>* a) {
	return true;
}


int* global_dsc_line_ids = NULL;
int* global_asc_line_ids = NULL;

template <unsigned char Dim, class FType>
void DumpLines(/*UnstructuredSimplicialComplex<Dim, FType>* bcc*/) {
	printf("dumping ids, %d, %d\n", (gX-1) ,(gY-1));
	FILE* fdsc = fopen("dump_asc_lines.raw", "wb");
	
	if(global_asc_line_ids == NULL) {
		global_asc_line_ids = new int[(gX-1) * (gY-1)];
		for (int i = 0; i < (gX-1) * (gY-1); i++) global_asc_line_ids[i] = 0;
	}

	for (int i = 0; i < MSC->arcs.size(); i++) {
		arc<float>* a = MSC->arcs[i];
		if (! MSC->isAlive(a)) continue;
		if (a->lower->index == 0 && ! DRAWDSCLINES) continue;
		if (a->lower->index == 1 && ! DRAWASCLINES) continue;

		if (! ArcTest<Dim, FType>(a)) continue;

		if (a->lower->index != 1) continue;

		vector<CELL_INDEX_TYPE> g;
		MSC->fillGeometry(a, g);

		for (int j = 0; j < g.size(); j++) {
			
			if (G_mscmh->dimension(g[j]) != 2) continue;
			int x = (g[j] % XMAX) / 2;
			int y = (g[j] / XMAX) / 2;
			global_asc_line_ids[x + y * (gX-1)] =128;

			//global_dsc_line_ids[g[j]] = 128;
		}
	}
	fwrite(global_asc_line_ids, (gX-1)*(gY-1), sizeof(int), fdsc);
	fclose(fdsc);
	fdsc = fopen("dump_dsc_lines.raw", "wb");
	
	if(global_dsc_line_ids == NULL) {
		global_dsc_line_ids = new int[gX * gY];
		for (int i = 0; i < (gX) * (gY); i++) global_dsc_line_ids[i] = 0;
	}

	for (int i = 0; i < MSC->arcs.size(); i++) {
		arc<float>* a = MSC->arcs[i];
		if (! MSC->isAlive(a)) continue;
		if (a->lower->index == 0 && ! DRAWDSCLINES) continue;
		if (a->lower->index == 1 && ! DRAWASCLINES) continue;

		if (! ArcTest<Dim, FType>(a)) continue;

		if (a->lower->index != 0) continue;

		vector<CELL_INDEX_TYPE> g;
		MSC->fillGeometry(a, g);

		for (int j = 0; j < g.size(); j++) {
			
			if (G_mscmh->dimension(g[j]) != 0) continue;
			int x = (g[j] % XMAX) / 2;
			int y = (g[j] / XMAX) / 2;
			global_dsc_line_ids[x + y * (gX)] =128;

			//global_dsc_line_ids[g[j]] = 128;
		}
	}
	fwrite(global_dsc_line_ids, (gX)*(gY), sizeof(int), fdsc);
	fclose(fdsc);
	//FILE* fasc = fopen("dump_asc_lines.raw", "wb");

	//if(global_asc_line_ids == NULL) {
	//	global_asc_line_ids = new int[gX * gY];
	//	for (int i = 0; i < gX*gY; i++) global_asc_line_ids[i] = 0;
	//}
	// nit = MSC->nodes.begin();
	//while (nit != MSC->nodes.end()) {
	//	node<float>* n = (*nit).second;
	//	nit++;
	//
	//	if (! MSC->isAlive(n)) continue;
	//	if (n->index != 0) continue;
	//	set<CELL_INDEX_TYPE> g;
	//	MSC->fillGeometry(n, g);

	//	for (set<CELL_INDEX_TYPE>::iterator sit = g.begin(); 
	//		sit != g.end(); sit++) {
	//		CELL_INDEX_TYPE currentcid = *sit;
	//		int x = (currentcid % XMAX) / 2;
	//		int y = (currentcid / XMAX) / 2;
	//		global_asc_ids[x + y * (gX)] = n->cellid;
	//		//printf("[%d, %d] = %d, %d\n", x,y,n->cellid, gX);
	//	}
	//}
	//fwrite(global_asc_ids, (gX)*(gY), sizeof(int), fasc);
	//fclose(fasc);


}
float myval = 1.0f;
bool g_use_test = false;

bool PassTest(arc<float>* a) {
	if (!g_use_test) return true;
	return a->lower->value < myval && a->upper->value < myval;
}
bool PassTest(node<float>* n) {
	if (!g_use_test) return true;
	return n->value < myval;
}
struct color {
	float r, g, b, a;
	color(float rr, float gg, float bb, float aa) : r(rr), g(gg), b(bb), a(aa){}
	color(){}
};

map<node<float>*, color> node_cols;

template <unsigned char Dim, class FType>
void drawStuff(/*UnstructuredSimplicialComplex<Dim, FType>* bcc*/) {
	if (! redrawstuff) {
		glCallList(drawlist);
		return;
	}

	if (drawlist != -1) 
		glDeleteLists(drawlist, 1);

	drawlist = glGenLists(1);
	glNewList(drawlist, GL_COMPILE_AND_EXECUTE);

	redrawstuff = false;


	//glBegin(GL_LINE_LOOP);

	//glColor3f(1,0,0);
	//glVertex3f(0,0,0);
	//glColor3f(0, 1, 0);
	//glVertex3f(0, YMAX-1, 0);

	//glColor3f(0,0,1);
	//glVertex3f(XMAX-1, YMAX-1, 0);
	//
	//glColor3f(1,1,0);
	//glVertex3f(XMAX-1, 0, 0);

	//glEnd();

	mscSimpleGradientUsingAlgorithms<float>* alg = new
		mscSimpleGradientUsingAlgorithms<float>(G_mscf, G_mscmh, G_mscg, NULL);

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0, 1.0);

	cellIterator it;
	if(DRAW_BACKGROUND) {

		glBegin(GL_QUADS);
		iteratorOperator& all = G_mscmh->all_cells_iterator(it);
		for (all.begin(it); all.valid(it); all.advance(it)) {
			//if (G_mscg->get_dim_asc_man(all.value(it)) != 3) 
			CELL_INDEX_TYPE cid = all.value(it);

			float sc = G_mscf->cell_value(cid);
			sc = (sc - fminmax.minval) / (fminmax.maxval - fminmax.minval);
			sc = sc*.8 + .2;
			glColor3f(sc+.01, sc, sc+.01);

			float c[3]; 
			coordinates(cid, c);
			glVertex3f(c[0]-0.5, c[1]-0.5, c[2]);
			glVertex3f(c[0]-0.5, c[1]+0.5, c[2]);
			glVertex3f(c[0]+0.5, c[1]+0.5, c[2]);
			glVertex3f(c[0]+0.5, c[1]-0.5, c[2]);
			//glVertex3f(c[0], c[1], c[2]);


		}
		glEnd();



	}







	glDisable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


	// draw grid
	if (DRAW_BACKGROUND2 || DRAW_BACKGROUND3) {
		DumpAscMan<Dim,FType>();
	}


	glBegin(GL_QUADS);
	typename map<CELL_INDEX_TYPE, node<float>* >::iterator nit = MSC->nodes.begin();
	if (DRAW_BACKGROUND2) {

		while (nit != MSC->nodes.end()) {
			node<float>* n = (*nit).second;
			nit++;
			if (! MSC->isAlive(n)) continue;
			if (n->index != 2) continue;

			set<CELL_INDEX_TYPE> g;
			MSC->fillGeometry(n, g);
			float col[3];
			coordinates(n->cellid, col);
#define PI 3.14159265f

			int dd = 60;
			float hd = 60.0f;

			//float v1 = (sin(PI*((int) col[0]%dd)/((float) hd)))*0.8f+.2f;
			//float v2 = (sin((dd - (((int) col[1])%dd))/hd))*0.8f+.2f;
			//float v3 = (sin(PI*((((int) col[0])%dd)-((int) col[1])%dd)/hd))*0.8f+.2f;
			if (node_cols.count(n) == 0) {
				float v1 = (rand() % 100) / 100.0f; // (sin(PI*((int) col[0]%dd)/((float) hd)))*0.8f+.2f;
				float v2 = (rand() % 100) / 100.0f; // (sin((dd - (((int) col[1])%dd))/hd))*0.8f+.2f;
				float v3 = (rand() % 100) / 100.0f; //(sin(PI*((((int) col[0])%dd)-((int) col[1])%dd)/hd))*0.8f+.2f;
				node_cols[n] = color(v1, v2, v3, 1.0);
			}
			color& mycolor = node_cols[n];
			glColor4f(mycolor.r, mycolor.g, mycolor.b, .6f);
			for (set<CELL_INDEX_TYPE>::iterator sit = g.begin(); 
				sit != g.end(); sit++) {
					CELL_INDEX_TYPE currentcid = *sit;
					float c[3]; 
					coordinates(currentcid, c);
					glVertex3f(c[0]-1.0, c[1]-1.0, c[2]);
					glVertex3f(c[0]-1.0, c[1]+1.0, c[2]);
					glVertex3f(c[0]+1.0, c[1]+1.0, c[2]);
					glVertex3f(c[0]+1.0, c[1]-1.0, c[2]);
			}

		}

	}
	if (DRAW_BACKGROUND3) {

		nit = MSC->nodes.begin();
		while (nit != MSC->nodes.end()) {
			node<float>* n = (*nit).second;
			nit++;
			if (! MSC->isAlive(n)) continue;
			if (n->index != 0) continue;

			set<CELL_INDEX_TYPE> g;
			MSC->fillGeometry(n, g);
			float col[3];
			coordinates(n->cellid, col);
			int dd = 60;
			float hd = 60.0f;
			//float v1 = (sin(PI*((int) col[0]%dd)/((float) hd)))*0.8f+.2f;
			//float v2 = (sin((dd - (((int) col[1])%dd))/hd))*0.8f+.2f;
			//float v3 = (sin(PI*((((int) col[0])%dd)-((int) col[1])%dd)/hd))*0.8f+.2f;
			//float v1 = ((rand()%100) * 0.01 + .3); //sin(PI*((int) col[0]%dd)/((float) hd)))/2.0+.4f;
			//float v2 =((rand()%100) * 0.01 + .3);// (sin((dd - (((int) col[1])%dd))/hd))/2.0+.4f;
			//float v3 = ((rand()%100) * 0.01 + .3);///*(sin(PI*((((int) col[0])%dd)-((int) col[1])%dd)/hd))/2.0+*/.4f;
			if (node_cols.count(n) == 0) {
				float v1 = (rand() % 100) / 100.0f; // (sin(PI*((int) col[0]%dd)/((float) hd)))*0.8f+.2f;
				float v2 = (rand() % 100) / 100.0f; // (sin((dd - (((int) col[1])%dd))/hd))*0.8f+.2f;
				float v3 = (rand() % 100) / 100.0f; //(sin(PI*((((int) col[0])%dd)-((int) col[1])%dd)/hd))*0.8f+.2f;
				node_cols[n] = color(v1, v2, v3, 1.0);
			}
			color& mycolor = node_cols[n];
			glColor4f(mycolor.r, mycolor.g, mycolor.b, .6f);
			//float v1 = (sin(PI*((int) col[0]%dd)/((float) hd)))/2.0+.4f;
			//float v2 = (sin((dd - (((int) col[1])%dd))/hd))/2.0+.4f;
			//float v3 = /*(sin(PI*((((int) col[0])%dd)-((int) col[1])%dd)/hd))/2.0+*/.4f;
			//
			//glColor4f(v2,v3,  v1, .6f);
			for (set<CELL_INDEX_TYPE>::iterator sit = g.begin(); 
				sit != g.end(); sit++) {
					CELL_INDEX_TYPE currentcid = *sit;
					float c[3]; 
					coordinates(currentcid, c);
					glVertex3f(c[0]-1.0, c[1]-1.0, c[2]);
					glVertex3f(c[0]-1.0, c[1]+1.0, c[2]);
					glVertex3f(c[0]+1.0, c[1]+1.0, c[2]);
					glVertex3f(c[0]+1.0, c[1]-1.0, c[2]);
			}

		}

	}

	glEnd();











	if (DRAW_GRID) {

		glLineWidth(line_width*2.1);
		glBegin(GL_LINES);
		glColor4f(0,0,0, .8);

		for (int i = 0; i < XMAX; i+=2) {

			glVertex3f(i, 0, 0);
			glVertex3f(i, YMAX, 0);
		}
		for (int i = 0; i < YMAX; i+=2) {

			glVertex3f(0, i, 0);
			glVertex3f(XMAX, i, 0);
		}

		glEnd();

		glLineWidth(line_width*1.1);













	}

	if (g_draw_isolines) {

		glBegin(GL_LINES);

		vector<float> isovals;
		isovals.push_back(1e-10);
		for (int i = 0; i < 50; i++) {
			//isovals.push_back(isovals[i] * 2.0);
			isovals.push_back(fminmax.minval + ((i+1)/ 51.0f) * (fminmax.maxval - fminmax.minval));
		}
		//isovals[0] = 0.00000001;
		//isovals[1] = 0.0000001;
		//isovals[2] = 0.000001;
		//isovals[3] = 0.00001;
		//isovals[4] = 0.0001;
		//isovals[5] = 0.001;
		//isovals[6] = 0.01;
		//isovals[7] = 0.1;
		//isovals[8] = 1.0;
		//isovals[9] = 10.0;
		//isovals[10] = 100.0;

		if (DRAW_ISOLINES) {
			glColor4f(1,1,1, .2);

			iteratorOperator& mall = G_mscmh->d_cells_iterator(2, it);
			for (mall.begin(it); mall.valid(it); mall.advance(it)) {
				//if (G_mscg->get_dim_asc_man(all.value(it)) != 3) 
				CELL_INDEX_TYPE cid = mall.value(it);

				for (int i = 0; i < isovals.size(); i++) {
					cellIterator fit;
					iteratorOperator& facets = G_mscmh->facets(cid, fit);
					for (facets.begin(fit); facets.valid(fit); facets.advance(fit)) {
						CELL_INDEX_TYPE eid = facets.value(fit);

						if (G_mscf->cell_value(eid) >= isovals[i] && 
							G_msccb->lowest_facet_value(eid) < isovals[i]) {

								float c[2][3]; 
								float vf[2];
								int vert = 0;
								cellIterator vit;
								iteratorOperator& viter = G_mscmh->facets(eid,vit);
								for (viter.begin(vit); viter.valid(vit); viter.advance(vit)) {
									CELL_INDEX_TYPE vid = viter.value(vit);
									vf[vert]= G_mscf->cell_value(vid);
									coordinates(vid, c[vert++]);
								}

								float t = (isovals[i] - vf[0]) / (vf[1] - vf[0]);
								glVertex3f(c[0][0] + t*(c[1][0] - c[0][0]),
									c[0][1] + t*(c[1][1] - c[0][1]),
									c[0][2] + t*(c[1][2] - c[0][2]));


						}
					}
				}
			}

		}

		glEnd();
	}
	//	cellIterator it2;
	//iteratorOperator& ait = G_mscmh->all_cells_iterator(it);
	//for (ait.begin(it); ait.valid(it); ait.advance(it)) {
	//	CELL_INDEX_TYPE cid = ait.value(it);
	//	
	//	//if (G_mscmh->dimension(cid) == 1) { //->get_dim_asc_man(cid) != 2)
	//		drawArrow(cid);
	//	//}

	//}

	glLineWidth(line_width*0.5);

	//if (g_draw_num_grad) {


	drawSEG();


	//}
	glLineWidth(line_width*1.1);



	glBegin(GL_QUADS);


	//if (DRAW_PROBABILITIES1) {
	//
	//vector<idfpair>& v = G_msccb->getmaxvals();
	//for (int i = 0; i < v.size(); i++) {
	//	idfpair p = v[i];
	//	if (G_mscmh->dimension(p.id) != 0) continue;
	//	if (p.prob < .98) {
	//		glColor4f(1, 0,0, 0.5*(1-p.prob));
	//	float c[3]; 
	//	coordinates(p.id, c);
	//	glVertex3f(c[0]-1.0, c[1]-1.0, c[2]);
	//	glVertex3f(c[0]-1.0, c[1]+1.0, c[2]);
	//	glVertex3f(c[0]+1.0, c[1]+1.0, c[2]);
	//	glVertex3f(c[0]+1.0, c[1]-1.0, c[2]);
	//	}

	//}
	//}

	//if (DRAW_PROBABILITIES2) {
	//
	//if (G_msccb != G_msccb2) {
	//	vector<idfpair>& v2 = G_msccb2->getmaxvals();
	//	for (int i = 0; i < v2.size(); i++) {
	//		idfpair p = v2[i];
	//	if (G_mscmh->dimension(p.id) != 2) continue;
	//		if (p.prob < .98) {
	//			glColor4f(0, 0,1, 0.5*(1-p.prob));
	//			float c[3]; 
	//			coordinates(p.id, c);
	//			glVertex3f(c[0]-1.0, c[1]-1.0, c[2]);
	//			glVertex3f(c[0]-1.0, c[1]+1.0, c[2]);
	//			glVertex3f(c[0]+1.0, c[1]+1.0, c[2]);
	//			glVertex3f(c[0]+1.0, c[1]-1.0, c[2]);
	//		}

	//	}


	//}
	//}
	glEnd();


	if(DRAWNUMERICALDE) {

		glLineWidth(0.5);
		glColor4f(0.1, 1.0, 0.1, 0.3);
		int numread = gX*gY;
		if (! DRAWNUMERICALDE_INIT) {
			FILE* fin = fopen("source_dest.raw", "rb");
			NUMDE_sources = new int[numread];
			NUMDE_dests = new int[numread];
			fread(NUMDE_sources, sizeof(int), numread, fin);
			fread(NUMDE_dests, sizeof(int), numread, fin);
			fclose(fin);
			DRAWNUMERICALDE_INIT = true;
			for (int i = 0; i < numread; i++) {
				if (sd2c.count(NUMDE_sources[i]) == 0 ||
					sd2c[NUMDE_sources[i]].count(NUMDE_dests[i]) == 0) {
						TCOLOR ccc; 
						ccc.r = (rand() % 100) * 0.01f;
						ccc.g = (rand() % 100) * 0.01f;
						ccc.b = (rand() % 100) * 0.01f;
						ccc.count = 1;
						sd2c[NUMDE_sources[i]][NUMDE_dests[i]] = ccc;
				} else {
					sd2c[NUMDE_sources[i]][NUMDE_dests[i]].count++;
				}
			}

			FILE* fout = fopen("DISSIPATIONELEMENTS_num.txt", "w");
			//for (map<int,  map<int, TCOLOR> >::iterator it = sd2c.begin();
			//	it != sd2c.end(); it++) {
			//		int src = (*it).first;
			//		int srcCID = (src % gX)*2 + (src / gX)*2 * XMAX;
			//		float srcF = G_mscf->cell_value(srcCID);
			//		for (map<int, TCOLOR>::iterator it2 = (*it).second.begin();
			//			it2 != (*it).second.end(); it2++) {
			//				int dst = (*it2).first;
			//				int dstCID = (dst % gX)*2 + (dst / gX)*2 * XMAX;
			//				float dstF = G_mscf->cell_value(dstCID);
			//				TCOLOR& c = (*it2).second;

			//				fprintf(fout, "%f %f %f %d %f %d %d %f %f %f %f\n", 
			//					dstF,srcF, srcF-dstF, 
			//					c.count, distance(srcCID, dstCID), 0, 0, 
			//					(dst % gX)*2.0f, (dst / gX)*2.0f, (src % gX)*2.0f, (src / gX)*2.0f);


			//		}



			//}
			fclose(fout);		


		}



		//		it != paiint, rmap.end(); it++) {
		//			DE& d = (*it).second;
		//
		//			float c1[3]; 
		//			coordinates(d.ndown->cellid, c1);
		//			float c2[3]; 
		//			coordinates(d.nup->cellid, c2);
		//	
		//			fprintf(fout, "%f %f %f %d %f %d %d %f %f %f %f\n", d.ndown->value, d.nup->value,
		//				d.nup->value - d.ndown->value, d.size, distance(d.nup->cellid, d.ndown->cellid),
		//				d.interior, d.boundary, c1[0], c1[1], c2[0], c2[1]);
		//	}


		//glBegin(GL_LINES);
		//for (int i = 0; i < numread; i++) {
		//	//glVertex3f((float) (i % gX) * 2, float (i / gX) * 2, 0);
		//	glVertex3f((float) (NUMDE_sources[i] % gX) * 2, float (NUMDE_sources[i] / gX) * 2, 0);

		//	//glVertex3f((float) (i % gX) * 2, float (i / gX) * 2, 0);
		//	glVertex3f((float) (NUMDE_dests[i] % gX) * 2, float (NUMDE_dests[i] / gX) * 2, 0);
		//}
		//glEnd();
		//glPointSize(2.0);
		glBegin(GL_QUADS);


		for (int i = 0; i < numread; i++) {
			TCOLOR& ttt = sd2c[NUMDE_sources[i]][NUMDE_dests[i]];
			glColor4f(ttt.r, ttt.g, ttt.b, 0.4);
		
			float x = (i % gX) * 2.0f;
			float y = (i / gX) * 2.0f;
			glVertex3f(x-1.0, y-1.0, 0);
			glVertex3f(x-1.0, y+1.0, 0);
			glVertex3f(x+1.0, y+1.0, 0);
			glVertex3f(x+1.0, y-1.0, 0);

			//glVertex3f((float) (i % gX) * 2, float (i / gX) * 2, 0);
		}
		glEnd();
		glLineWidth(line_width*2.5);
		glBegin(GL_LINES);
		glColor3f(1,0,0);
		for (map<int,  map<int, TCOLOR> >::iterator it = sd2c.begin();
			it != sd2c.end(); it++) {
			int src = (*it).first;
			float srcx = (src % gX)*2.0f;
			float srcy = (src / gX)*2.0f;
			for (map<int, TCOLOR>::iterator it2 = (*it).second.begin();
				it2 != (*it).second.end(); it2++) {
				int dst = (*it2).first;
				float dstx = (dst % gX)*2.0f;
				float dsty = (dst / gX)*2.0f;
				glVertex3f(srcx, srcy, 0.0);
				glVertex3f(dstx, dsty, 0.0);
			}
		}
		glEnd();
			
		glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
		//glPointSize(point_size);
		glPointSize(point_size);
		glBegin(GL_POINTS);
		
		for (map<int,  map<int, TCOLOR> >::iterator it = sd2c.begin();
			it != sd2c.end(); it++) {
			int src = (*it).first;
			float srcx = (src % gX)*2.0f;
			float srcy = (src / gX)*2.0f;
			for (map<int, TCOLOR>::iterator it2 = (*it).second.begin();
				it2 != (*it).second.end(); it2++) {
				int dst = (*it2).first;
				float dstx = (dst % gX)*2.0f;
				float dsty = (dst / gX)*2.0f;
				glColor3f( 0, 0,0.6);
				glVertex3f(dstx, dsty, 0.0);
				glColor3f(0.6, 0, 0);
				glVertex3f(srcx, srcy, 0.0);
			}
		}
		glEnd();
		glPointSize(point_size*.75);
		glBegin(GL_POINTS);
		for (map<int,  map<int, TCOLOR> >::iterator it = sd2c.begin();
			it != sd2c.end(); it++) {
			int src = (*it).first;
			float srcx = (src % gX)*2.0f;
			float srcy = (src / gX)*2.0f;
			for (map<int, TCOLOR>::iterator it2 = (*it).second.begin();
				it2 != (*it).second.end(); it2++) {
				int dst = (*it2).first;
				float dstx = (dst % gX)*2.0f;
				float dsty = (dst / gX)*2.0f;
				glColor3f( .3, .3,1);
				glVertex3f(dstx, dsty, 0.0);
				glColor3f(1, 0.3, 0.3);
				glVertex3f(srcx, srcy, 0.0);
			}
		}
		glEnd();

		////	glBegin(GL_QUADS);
		////	iteratorOperator& all2 = G_mscmh->d_cells_iterator(2,it);
		////for (all2.begin(it); all2.valid(it); all2.advance(it)) {
		////	//if (G_mscg->get_dim_asc_man(all.value(it)) != 3) 
		////	//drawArrow(all.value(it));
		////	CELL_INDEX_TYPE cid = all2.value(it);
		////	
		////	if (G_mscg->get_dim_asc_man(cid) == 2) continue;
		////	if (G_mscg->get_dim_asc_man(cid) == 1) glColor4f(.2, 1.0, .3,.6);
		////	if (G_mscg->get_dim_asc_man(cid) == 0) glColor4f(.4, 1.0, .7,.6); 
		////	
		////	float c[3]; 
		////	coordinates(cid, c);
		////	glVertex3f(c[0]-1.0, c[1]-1.0, c[2]);
		////	glVertex3f(c[0]-1.0, c[1]+1.0, c[2]);
		////	glVertex3f(c[0]+1.0, c[1]+1.0, c[2]);
		////	glVertex3f(c[0]+1.0, c[1]-1.0, c[2]);
		////	//glVertex3f(c[0], c[1], c[2]);


		////}
		////glEnd();
	}

	//glBegin(GL_POINTS);
	//cellIterator fuck;
	//iteratorOperator& allfuck = G_mscmh->all_cells_iterator(fuck);
	//for (allfuck.begin(fuck); allfuck.valid(fuck); allfuck.advance(fuck)) {
	//	CELL_INDEX_TYPE asdf = allfuck.value(fuck);
	//	DIM_TYPE ddd = G_mscmh->boundary_value(asdf); // G_mscg_TEMP->get_dim_asc_man(asdf);
	//	if (ddd != 0) {
	//		float c[3]; 
	//		coordinates(asdf, c);
	//		if (ddd == 1) { 
	//			glColor3f(1,0,0);
	//		} else if (ddd == 2) {
	//			glColor3f(0,0,1);
	//		} else if (ddd == 3) {
	//			glColor3f(0,1,0);
	//		} else {
	//			glColor3f(0,0,0);
	//		}

	//		glVertex3f(c[0], c[1], 0);
	//	}
	//}
	//glEnd();


	if(DRAWDIMASCMAN) {
		glBegin(GL_QUADS);
		iteratorOperator& all2 = G_mscmh->d_cells_iterator(2,it);
		for (all2.begin(it); all2.valid(it); all2.advance(it)) {
			//if (G_mscg->get_dim_asc_man(all.value(it)) != 3) 
			//drawArrow(all.value(it));
			CELL_INDEX_TYPE cid = all2.value(it);

			if (G_mscg->get_dim_asc_man(cid) == 2) continue;
			if (G_mscg->get_dim_asc_man(cid) == 1) glColor4f(.2, 1.0, .3,.6);
			if (G_mscg->get_dim_asc_man(cid) == 0) glColor4f(.4, 1.0, .7,.6); 

			float c[3]; 
			coordinates(cid, c);
			glVertex3f(c[0]-1.0, c[1]-1.0, c[2]);
			glVertex3f(c[0]-1.0, c[1]+1.0, c[2]);
			glVertex3f(c[0]+1.0, c[1]+1.0, c[2]);
			glVertex3f(c[0]+1.0, c[1]-1.0, c[2]);
			//glVertex3f(c[0], c[1], c[2]);


		}
		glEnd();
	}

	if(DRAWGRAD) {
		iteratorOperator& all3 = G_mscmh->all_cells_iterator(it);
		for (all3.begin(it); all3.valid(it); all3.advance(it)) {
			//if (G_mscg->get_dim_asc_man(all.value(it)) != 3) 
			drawArrow(all3.value(it));
		}


	}

	// white out stuff!
	if (clearcolor)
		glColor3f(1,1, 1);
	else 
		glColor3f(0,0,0);
	glBegin(GL_QUADS);
	glVertex3f(-2.0,-2.0, 0.0);
	glVertex3f(2.0+XMAX,-2.0, 0.0);
	glVertex3f(2.0+XMAX,0.0, 0.0);
	glVertex3f(-2.0,0.0, 0.0);

	glVertex3f(-2.0,YMAX+2.0, 0.0);
	glVertex3f(2.0+XMAX,YMAX+2.0, 0.0);
	glVertex3f(2.0+XMAX,YMAX-1.0, 0.0);
	glVertex3f(-2.0,YMAX-1.0, 0.0);

	glVertex3f(-2.0,-2.0, 0.0);
	glVertex3f(-2.0,YMAX+2.0, 0.0);
	glVertex3f(0.0,YMAX+2.0, 0.0);
	glVertex3f(0.0,-2.0, 0.0);

	glVertex3f(XMAX-1.0,-2.0, 0.0);
	glVertex3f(XMAX-1.0,YMAX+2.0, 0.0);
	glVertex3f(XMAX+2.0,YMAX+2.0, 0.0);
	glVertex3f(XMAX+2.0,-2.0, 0.0);

	glEnd();


	//
	//glBegin(GL_QUADS);
	//for (){
	//	

	//	float sc = G_mscf->cell_value(cid);
	//	sc = (sc - fminmax.minval) / (fminmax.maxval - fminmax.minval);

	//	glColor3f(sc, sc, sc);

	//	float c[3]; 
	//	coordinates(cid, c);
	//	glVertex3f(c[0]-0.5, c[1]-0.5, c[2]);
	//	glVertex3f(c[0]-0.5, c[1]+0.5, c[2]);
	//	glVertex3f(c[0]+0.5, c[1]+0.5, c[2]);
	//	glVertex3f(c[0]+0.5, c[1]-0.5, c[2]);
	//	//glVertex3f(c[0], c[1], c[2]);


	//}
	//glEnd();


	glDisable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(0.0, 0.0);




	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glLineWidth(line_width*2.5);
	//if (DRAWPOINTS) DrawDe();


	if (g_use_test) {

		glBegin(GL_LINES);


		glColor4f(0.6, 0, .8, 1.0);

		iteratorOperator& mall = G_mscmh->d_cells_iterator(2, it);
		for (mall.begin(it); mall.valid(it); mall.advance(it)) {
			//if (G_mscg->get_dim_asc_man(all.value(it)) != 3) 
			CELL_INDEX_TYPE cid = mall.value(it);
			cellIterator fit;
			int countintersects = 0;
			iteratorOperator& facets = G_mscmh->facets(cid, fit);
			for (facets.begin(fit); facets.valid(fit); facets.advance(fit)) {
				CELL_INDEX_TYPE eid = facets.value(fit);

//				if (G_mscf->cell_value(eid) >= myval &&
//					G_msccb->lowest_facet_value(eid) < myval) {

					float c[2][3];
					float vf[2];
					int vert = 0;
					cellIterator vit;
					iteratorOperator& viter = G_mscmh->facets(eid, vit);
					for (viter.begin(vit); viter.valid(vit); viter.advance(vit)) {
						CELL_INDEX_TYPE vid = viter.value(vit);
						vf[vert] = G_mscf->cell_value(vid);
						coordinates(vid, c[vert++]);
					}

					float t = (myval - vf[0]) / (vf[1] - vf[0]);
					if (t < 0 || t >= 1.0f) continue;
					glVertex3f(c[0][0] + t*(c[1][0] - c[0][0]),
						c[0][1] + t*(c[1][1] - c[0][1]),
						c[0][2] + t*(c[1][2] - c[0][2]));
					countintersects++;

			//	}
			}
			if (countintersects % 2 == 1) {
				printf("WHOATHERE uneven stuff\n");
				//glVertex3f(0, 0, 0);
			}



		}

		glEnd();
	}









	glLineWidth(line_width*3.1);

	for (int i = 0; i < MSC->arcs.size(); i++) {
		arc<float>* a = MSC->arcs[i];
		if (! MSC->isAlive(a) || ! PassTest(a)) continue;
		if (a->lower->index == 0 && ! DRAWDSCLINES) continue;
		if (a->lower->index == 1 && ! DRAWASCLINES) continue;

		if (! ArcTest<Dim, FType>(a)) continue;

		if (a->lower->index == 0) {
			glColor4f(.1, .2, 1,.8);
			//printf("ASDFL:KJSDFL:KJSDKL:FJSD\n");
		} else {
			glColor4f(1, .2, .1,.8);
		}

		vector<CELL_INDEX_TYPE> g;
		MSC->fillGeometry(a, g);

		glBegin(GL_LINE_STRIP);
		for (int j = 0; j < g.size(); j++) {
			float c[3]; 
			coordinates(g[j], c);
			glVertex3f(c[0], c[1], c[2]);
		}
		glEnd();
	}

	glLineWidth(line_width*2.5);
	for (int i = 0; i < MSC->arcs.size(); i++) {
		arc<float>* a = MSC->arcs[i];
		if (!MSC->isAlive(a) || !PassTest(a)) continue;

		if (a->lower->index == 0 && ! DRAWDSCLINES) continue;
		if (a->lower->index == 1 && ! DRAWASCLINES) continue;
		if (! ArcTest<Dim, FType>(a)) continue;

		if (a->lower->index == 0) {
			glColor4f(.6, .6, .9,.3);
			//printf("ASDFL:KJSDFL:KJSDKL:FJSD\n");
		} else {
			glColor4f(.9, .6, .6,.3);
		}


		//glColor4f(1,1,1,.5);
		vector<CELL_INDEX_TYPE> g;
		MSC->fillGeometry(a, g);

		glBegin(GL_LINE_STRIP);
		for (int j = 0; j < g.size(); j++) {
			float c[3]; 
			coordinates(g[j], c);
			glVertex3f(c[0], c[1], c[2]);
		}
		glEnd();
	}




	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glPointSize(point_size);

	if (DRAWPOINTS){
		glBegin(GL_POINTS);
		map<CELL_INDEX_TYPE, node<float>*>::iterator nit = MSC->nodes.begin();
		while (nit != MSC->nodes.end()) {
			node<float>* n = (*nit).second;
			nit++;
			if (! MSC->isAlive(n) || ! PassTest(n)) continue;
			switch(n->index) {
			case 0: glColor3f(0,0,.6); break;
			case 1: glColor3f(0,.6,0); break;
			case 2: glColor3f(.6,0,0); break;
			default: glColor3f(.3, .7, .7);
			}

			float c[3]; 
			coordinates(n->cellid, c);
			glVertex3f(c[0], c[1], c[2]);
			//glVertex3f(c[0]-0.6, c[1]-0.6, c[2]);
			//glVertex3f(c[0]-0.6, c[1]+0.6, c[2]);
			//glVertex3f(c[0]+0.6, c[1]+0.6, c[2]);
			//glVertex3f(c[0]+0.6, c[1]-0.6, c[2]);
		}
		glEnd();

		glPointSize(point_size*.75);
		glBegin(GL_POINTS);
		nit = MSC->nodes.begin();
		while (nit != MSC->nodes.end()) {
			node<float>* n = (*nit).second;
			nit++;
			if (!MSC->isAlive(n) || !PassTest(n)) continue;
			switch(n->index) {
			case 0: glColor3f(.3,.3,1); break;
			case 1: glColor3f(.3,1,.3); break;
			case 2: glColor3f(1,.3,.3); break;
			default: glColor3f(.3, .7, .7);
			}
			float c[3]; 
			coordinates(n->cellid, c);
			glVertex3f(c[0], c[1], c[2]);
			//glVertex3f(c[0]-0.6, c[1]-0.6, c[2]);
			//glVertex3f(c[0]-0.6, c[1]+0.6, c[2]);
			//glVertex3f(c[0]+0.6, c[1]+0.6, c[2]);
			//glVertex3f(c[0]+0.6, c[1]-0.6, c[2]);
		}

		glEnd();
	}



	glDisable(GL_BLEND);

	glEnable(GL_DEPTH_TEST);


	//if (draw_gradient) {
	//	//      glLineWidth(arrow_width);
	//	//      glColor3f(0, 0, 0);
	//	//      glBegin(GL_LINES);

	//	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	//	glEnable(GL_LIGHTING);

	//	for (int dd = 1; dd <= 2; dd++) {
	//		it = bcc->getCellIterator(dd);
	//		while (it.isValid()) {
	//			index_type cellid = *it.loc;
	//			drawArrow<Dim, FType>(bcc, cellid);
	//			it++;
	//		}
	//	}
	//	glDisable(GL_LIGHTING);
	//	//      glEnd();
	//}
	glEndList();




	glEndList();
};


//
//
//
//
//	if (! redrawstuff) {
//		glCallList(drawlist);
//		return;
//	}
//	
//	if (drawlist != -1) 
//		glDeleteLists(drawlist, 1);
//
//	drawlist = glGenLists(1);
//	glNewList(drawlist, GL_COMPILE_AND_EXECUTE);
//
//	redrawstuff = false;
//	
//	glEnable(GL_POLYGON_OFFSET_FILL);
//	glPolygonOffset(1.0, 1.0);
//
//
//// TRIANGLES	
//
//   //  glColor3f(.2, .2, .2); 
//   glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
//   if(! draw_flat)  glEnable(GL_LIGHTING);
//   glBegin(GL_TRIANGLES);
//   IndexIterator  it = bcc->getCellIterator(Dim);
//   while (it.isValid()) {
//      index_type cellid = *it.loc;
//      //total hack!
//      Simplex<FType>& s = bcc->cells[cellid];
//      if (s.numberOfVerts != 3) {
//         printf("SDJKFHDS:LJFSKL:D\n");
//      }
//
//      BaseVertex<Dim>& v1 = bcc->verts[s.verts[0] ]; 
//      BaseVertex<Dim>& v2 = bcc->verts[s.verts[1] ]; 
//      BaseVertex<Dim>& v3 = bcc->verts[s.verts[2] ];
//
//      bool is_assigned = bcc->getAssigned(cellid);
//	  //if (bcc->getNumUFacets(cellid) == 1) {
//	//	  glColor3f(0,1,0);
//	  //}
//      if (is_assigned) {
//        
//		  
//		  unsigned char asif = bcc->getDimAscMan(cellid);
//         if (asif == 2) {
//            glColor3f(0.8, 0.8, 0.0);
//         }	else if (asif == 1){ 
//            glColor3f(0.3, 0.3, 0.9);
//         } else {
//            glColor3f(1.0, 0, 0); 
//         }
//      }
//
//	  if (drawbasincolor) {
//		  float* c = bsn->livingCellColor(cellid, glivingcounter);
//		  glColor3f(c[0],c[1],c[2]);
//	  }
//
//      if (!draw_flat) {
//         // normal for lighting
//         Vector4 vv1; 
//         vv1.vals[0] = v1.position[0];
//         vv1.vals[1] = bcc->getValue(s.verts[0]);
//         vv1.vals[2] = v1.position[1];
//         vv1.vals[3] = 1.0f;
//         Vector4 vv2; 
//         vv2.vals[0] = v2.position[0];
//         vv2.vals[1] = bcc->getValue(s.verts[1]);
//         vv2.vals[2] = v2.position[1];
//         vv2.vals[3] = 1.0f;
//         Vector4 vv3; 
//         vv3.vals[0] = v3.position[0];
//         vv3.vals[1] = bcc->getValue(s.verts[2]);
//         vv3.vals[2] = v3.position[1];
//         vv3.vals[3] = 1.0f;
//         for (int i = 0; i < 4; i++) {
//            vv1.vals[i] = vv3.vals[i] - vv1.vals[i];
//            vv2.vals[i] = vv3.vals[i] - vv2.vals[i];
//         }
//         Normalize4(&vv1);
//         Normalize4(&vv2);
//         Vector4 n = Cross(vv1, vv2);
//         Normalize4(&n);
//         glNormal3f(n.vals[0], n.vals[1], n.vals[2]);
//      }
//
//      unsigned char asif = bcc->getDimAscMan(cellid);
//      if (! drawbasincolor && ! is_assigned /*&& bcc->getNumUFacets(cellid) == 0*/ && flat_funct) {
//	float nval = ((float) bcc->getValue(cellid) - (float) fminmax.minval)
//            /((float) fminmax.maxval - (float) fminmax.minval); 
//         float r = red_scale(nval);
//         float g = green_scale(nval);
//         float b = blue_scale(nval);
//         glColor3f(r,g,b);
//	
//
//      }
//      if (! drawbasincolor && ! is_assigned  /*&& bcc->getNumUFacets(cellid) == 0*/ && ! flat_funct){ 
//	float nval = ((float) bcc->getValue(s.verts[0]) - (float) fminmax.minval)
//	  /((float) fminmax.maxval - (float) fminmax.minval); 
//         float r = red_scale(nval);
//         float g = green_scale(nval);
//         float b = blue_scale(nval);
//         glColor3f(r,g,b);
//      }
//
//      if (! draw_flat) 
//         glVertex3f(v1.position[0],bcc->getValue(s.verts[0]),v1.position[1]);
//      else       
//         glVertex3f(v1.position[0],bcc->getValue(cellid),v1.position[1]);
//
//      if (! drawbasincolor && ! is_assigned  /*&& bcc->getNumUFacets(cellid) == 0 */&& ! flat_funct){ 
//         float nval = ((float) bcc->getValue(s.verts[1]) - (float) fminmax.minval)
//            /((float) fminmax.maxval - (float) fminmax.minval); 
//         float r = red_scale(nval);
//         float g = green_scale(nval);
//         float b = blue_scale(nval);
//         glColor3f(r,g,b);
//      } 
//      if (!draw_flat)
//         glVertex3f(v2.position[0],bcc->getValue(s.verts[1]),v2.position[1]);
//      else
//         glVertex3f(v2.position[0],bcc->getValue(cellid),v2.position[1]);
//
//      if (! drawbasincolor && ! is_assigned /*&& bcc->getNumUFacets(cellid) == 0*/ &&   ! flat_funct){ 
//         float nval = ((float) bcc->getValue(s.verts[2]) - (float) fminmax.minval)
//            /((float) fminmax.maxval - (float) fminmax.minval); 
//
//         float r = red_scale(nval);
//         float g = green_scale(nval);
//         float b = blue_scale(nval);
//         glColor3f(r,g,b);
//      }
//      
//      if (! draw_flat)
//         glVertex3f(v3.position[0],bcc->getValue(s.verts[2]),v3.position[1]);
//      else
//         glVertex3f(v3.position[0],bcc->getValue(cellid),v3.position[1]);
//
//
//
//      it++;
//   }
//   glEnd();
//   if (! draw_flat) glDisable(GL_LIGHTING);
//
//
//
//
//
//
//
//
//
//
//   // LINES
//   glDisable(GL_POLYGON_OFFSET_FILL);
//	glPolygonOffset(0.0, 0.0);
//
//	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
//	glEnable(GL_BLEND);
//   glLineWidth(line_width);
//
//
//
//   if (draw_mt) {
//		vector<index_type> vts;
//		mt->gl_fill_vertices2(vts);
//            
//		for (int i = 0; i < vts.size(); i+=2) {
//			BaseVertex<Dim>& v1 = bcc->verts[vts[i]];
//			index_type b = vts[i+1];
//			glColor3f(
//				(((float) ((b*321+93)%31))) /31.0f,
//				(((float) ((b*421+91)%51))) /51.0f,
//				(((float) ((b*221+92)%71))) /71.0f);
//			draw_earth(v1.position[0],bcc->getValue(vts[i]),v1.position[1], ballsize*3);
//		}
//		vts.clear();
//		mt->gl_fill_arcs(vts);
//		glBegin(GL_LINES);
//		for (int i = 0; i < vts.size(); i++) {
//			BaseVertex<Dim>& v1 = bcc->verts[vts[i]];
//			glVertex3f(v1.position[0],bcc->getValue(vts[i]) + EPSILON,v1.position[1]);
//		}
//		glEnd();
//		vts.clear();
//   }
//
//   if (draw_edges) {
//
//      glBegin(GL_LINES);
//      glColor4f(0.6, 0.5, 0.5, 1.0);
//      it = bcc->getCellIterator(1);
//      while (it.isValid()) {
//         index_type cellid = *it.loc;
//         //total hack!
//         Simplex<FType>& s = bcc->cells[cellid];
//         if (s.numberOfVerts != 2) {
//            printf("SDJKFHDS:LJFSKL:D\n");
//         }
//
//	 if (gusecutoffhack && bsn->livingCellValue(cellid, glivingcounter) < 0.01f) {
//	   it++;
//	   continue;
//	 }
//
//         BaseVertex<Dim>& v1 = bcc->verts[s.verts[0] ]; 
//         BaseVertex<Dim>& v2 = bcc->verts[s.verts[1] ]; 
//		 if (drawbasincolor) {
//			 float* c = bsn->livingCellColor(cellid, glivingcounter);
//			 glColor3f(c[0],c[1],c[2]);
//           
//			 glVertex3f(v1.position[0],bcc->getValue(s.verts[0]) + EPSILON,v1.position[1]);
//		     glVertex3f(v2.position[0],bcc->getValue(s.verts[1]) + EPSILON,v2.position[1]);		
//		 
//		 } else if (!draw_flat) {
//				 if (bcc->getNumUFacets(cellid) == 1) {
//
//					 //glColor4f(0,1, 0, 1.0);
//
//				  		
//				 if (bcc->getAssigned(cellid)) {
//
//
//					 if (bcc->getCritical(cellid)) {
//						 glColor4f(1,0,0, 1.0);
//					 } else {
//
//
//						 unsigned char asif = bcc->getDimAscMan(cellid);
//						 if (asif == 2) {
//							 glColor4f(.6, .6, 0, 1.0);
//						 } else if (asif == 1){ 
//							 glColor4f(0.3, 0.3, 0.9, 1.0);
//						 } else {
//							 glColor4f(1.0, 0, 0, 1.0); 
//						 }
//					 }
//				 }
//
//
//	     // float c[3];
//	     // centroid<Dim,FType>(bcc, cellid, c);
//	     // draw_earth(c[0],c[1],c[2],ballsize);
//            glVertex3f(v1.position[0],bcc->getValue(s.verts[0]) + EPSILON,v1.position[1]);
//            
//	 
//	     glVertex3f(v2.position[0],bcc->getValue(s.verts[1]) + EPSILON,v2.position[1]);
//	     
//				 } else {
//	     if (! flat_funct) {
//	       float nval = ((float) bcc->getValue(s.verts[0]) - (float) fminmax.minval)
//		 /((float) fminmax.maxval - (float) fminmax.minval); 
//	       float r = red_scale(nval);
//	       float g = green_scale(nval);
//	       float b = blue_scale(nval);
//	       float nval2 = ((float) bcc->getValue(s.verts[1]) - (float) fminmax.minval)
//		 /((float) fminmax.maxval - (float) fminmax.minval); 
//	       float r2 = red_scale(nval2);
//	       float g2 = green_scale(nval2);
//	       float b2 = blue_scale(nval2);
//	
//	       float a = 0.3f;
//	       glColor4f(r+a,g+a,b+a, 0.2);
//	       glVertex3f(v1.position[0],bcc->getValue(s.verts[0]) + EPSILON,v1.position[1]);
//	       
//	       glColor4f(r2+a, g2+a, b2+a, 0.2);
//	       glVertex3f(v2.position[0],bcc->getValue(s.verts[1]) + EPSILON,v2.position[1]);
//	     } else {
//	       float nval = ((float) bcc->getValue(cellid) - (float) fminmax.minval)
//		 /((float) fminmax.maxval - (float) fminmax.minval); 
//	       float r = red_scale(nval);
//	       float g = green_scale(nval);
//	       float b = blue_scale(nval);
//	       float a = 0.3f;
//	       glColor4f(r+a,g+a,b+a, 0.2);
//	       glVertex3f(v1.position[0],bcc->getValue(s.verts[0]) + EPSILON,v1.position[1]);
//	       glVertex3f(v2.position[0],bcc->getValue(s.verts[1]) + EPSILON,v2.position[1]);
//	     }
//	   }
//	   
//	   
//         } else {
//	   float nval = ((float) bcc->getValue(cellid) - (float) fminmax.minval)
//	       /((float) fminmax.maxval - (float) fminmax.minval); 
//	     float r = red_scale(nval);
//	     float g = green_scale(nval);
//	     float b = blue_scale(nval);
//
//	     float a = 0.3f;
//	     glColor4f(r+a,g+a,b+a, 0.2);
//            glVertex3f(v1.position[0],bcc->getValue(cellid) + EPSILON,v1.position[1]);
//            glVertex3f(v2.position[0],bcc->getValue(cellid) + EPSILON,v2.position[1]);
//         }
//
//
//         it++;
//      }
//      glEnd();
//
//   }
//   glDisable(GL_BLEND);
//
//
//   // VERTICES
//   it = bcc->getCellIterator(0);
//   index_type sizeshit = it.size;
//   while (it.isValid()) {
//      index_type cellid = *it.loc;
//      //total hack!
//      Simplex<FType>& s = bcc->cells[cellid];
//
//
//      BaseVertex<Dim>& v1 = bcc->verts[s.verts[0] ]; 
//
//      if (s.verts[0] != cellid) printf("ERORORORORORORO\n");
//
//
//      if (gusecutoffhack && bsn->livingCellValue(cellid, glivingcounter) < 0.01f) {
//	it++;
//	continue;
//      }
//      
//		 if (drawbasincolor && bcc->getCritical(cellid)) {
//			 float* c = bsn->livingCellColor(cellid, glivingcounter);			
//			 glColor3f(c[0],c[1],c[2]);             
//			 draw_earth(v1.position[0],bcc->getValue(s.verts[0]),v1.position[1], ballsize);		 		   
//			 
//		 } else if (bcc->getAssigned(cellid)) {
//		   
//         if (bcc->getCritical(cellid)) {
//            glColor3f(1,0,0);
//            draw_earth(v1.position[0],bcc->getValue(s.verts[0]),v1.position[1], ballsize);	
//         } else {
//            glColor3f(.6, .6, 0);
//            //draw_earth(v1.position[0],bcc->getValue(s.verts[0]),v1.position[1], ballsize/2.0);
//         }
//         //printf("reder min%d = %f %f %f\n", cellid,v1.position[0],bcc->getValue(s.verts[0]),v1.position[1] );
//
//		 } else {
//         float nval = ((float) bcc->getValue(s.verts[0]) - (float) fminmax.minval)
//            /((float) fminmax.maxval - (float) fminmax.minval); 
//         float r = red_scale(nval);
//         float g = green_scale(nval);
//         float b = blue_scale(nval);
//         glColor3f(r,g,b);
//         //draw_earth(v1.position[0],bcc->getValue(s.verts[0]),v1.position[1], ballsize/2.0);
//      }
//
//      it++;
//   }
//
//


//};
float gpers;

int main(int argc, char** argv) {

	if (argc < 6) {
		printf(
			"Usage: main_2d_viewer rawfilename X Y Z computemode [persistence]\n\
-- rawfilename is the name of a float32 binary file size X*Y*Z*4 bytes\n\
-- X Y Z are dimensions of the data (for the 2d viewer, Z=1\n\
-- computemod: 0=randomized, 1=steepest, 2-5=variations of accurate (hint: use 1)\n\
-- persistence: if a floating point number is supplied for persistence, the tool is COMMAND LINE ONLY and will dump nodes/arcs of the complex.\n");
		return 1;
	}

	for (int i =0; i < argc; i++) {
		printf("%d=%s\n", i, argv[i]);
	}

	char* filename = argv[1];
	int X, Y, Z;
	sscanf(argv[2], "%d", &X);
	sscanf(argv[3], "%d", &Y);
	sscanf(argv[4], "%d", &Z);

	gX = X;
	gY = Y;
	gZ = Z;
	XMIN = 0;
	XMAX = X*2-1;
	YMIN = 0; 
	YMAX = Y*2-1;
	ZMIN = 0; 
	ZMAX = 1;

	int userand;
	sscanf(argv[5], "%d", &userand);
	
	bool command_line_only = false;

	if (argc > 6) {
		sscanf(argv[6], "%f", &gpers);
		command_line_only = true;
	}





	mscSize stest(X, Y, Z);
	printf("stest=%d\n", (int) stest.count());



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
	//test_data->logify();
	//test_data->negate();
	mscNegatingDataHandler<float>* ndh = 
		new mscNegatingDataHandler<float>(test_data);
#ifdef WIN32        
	DWORD globalStart = GetTickCount();
#endif

	printf("Computing discrete gradient: \n");


	// create mesh handler
	mscRegularGrid2DImplicitMeshHandler* bmsh = new mscRegularGrid2DImplicitMeshHandler(X,Y);
	// tests!!!

	//printf("Go there\n");
	cellIterator it;
	iteratorOperator& fit = bmsh->cofacets(1, it);
	fit.begin(it);

	while (fit.valid(it)) {
		printf("neighbor=%llu\n", fit.value(it));
		fit.advance(it);
	}

	// ? ultimitely, what do I have to do with the data to use ttk
	// ? why use 3d mesh function
	// ?mesh handler assigns face/coface where id of sattle v critical done
	// with gradient field? once assigned every other approach for 2d used?
	// ? For TTK, can it compute the triangulation or do I assign the triang-
	// ulation and it traverses. (How to build triangulation, can it build)
	
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

	char gradname[2048];
	sprintf(gradname, "%s.grad", argv[1]);

	if (!bgf->load_from_file(gradname)) {


		mscConvergentGradientBuilder<float>* msccb =
			new mscConvergentGradientBuilder<float>(bmf, bmsh, bgf, &a);
		G_msccb2 = msccb;
		G_msccb = msccb;

		if (userand == 0) {
			printf("using randomized gradient\n");
			mscrb->computeGradient();
		}
		else if (userand == 1) {
			printf("using greedy gradient\n");
			mscb->computeGradient();
		}
		else if (userand == 2) {
			printf("using convergent gradient - 1 pass\n");
			msccb->computeGradient();
		}
		else if (userand == 3)  {
			printf("using convergent gradient - 2 pass\n");
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


			//char ufilename1[1024];
			//sprintf(ufilename1, "%s.asc", filename); 
			//test_data->dump_vals(ufilename1, X, Y, Z, msccb->getmaxvals());
			//sprintf(ufilename1, "%s.asc.prob", filename);

			//char ufilename2[1024];

			//sprintf(ufilename2, "%s.dsc", filename);
			//test_data->dump_vals(ufilename2, X, Y, Z, msccb2->getmaxvals());
			//		 sprintf(ufilename2, "%s.dsc.prob", filename);

			//char ufilename3[1024];

			//sprintf(ufilename3, "%s.both", filename);
			//FILE* ff1 = fopen(ufilename1, "rb");
			//FILE* ff2 = fopen(ufilename2, "rb");
			//FILE* ffout = fopen(ufilename3, "wb");
			//while (! feof(ff1)) {
			// float v1, v2;
			// fread(&v1, sizeof(float), 1, ff1);
			// fread(&v2, sizeof(float), 1, ff2);
			// float res = (0.5-0.5f*v1) + 0.5f - (0.5-0.5f*v2);
			// fwrite(&res, sizeof(float), 1, ffout);
			//}
			//fclose(ff1);
			//fclose(ff2);
			//fclose(ffout);

			G_msccb = msccb;
			G_msccb2 = msccb2;
		}
		else if (userand == 4) {
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
			G_msccb = msccb3;
			G_msccb2 = msccb2;
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



		}
		else if (userand == 5) {
			printf("using convergent3 gradient - 2 pass\n");


			mscComplementMeshHandler* cmh =
				new mscComplementMeshHandler(bmsh);
			//mscDumbStoringMinFunction<float>* dsminf = 
			// new mscDumbStoringMinFunction<float>(test_data, bmsh, &a);
			//dsminf->initialize();
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
			G_msccb = msccb3;
			G_msccb2 = msccb2;
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



		}
		else if (userand == 6/* && argc > 6*/) {
			// READ IN SEGMENTATION
			//mscb->computeGradient();
			USE_SEG = true;
			FILE* fseg = fopen("source_dest2.raw", "rb");
			int numV = X * Y;
			int numC = (X - 1) * (Y - 1);

			int* ascID = new int[numV];
			int* dscID = new int[numC];
			fread(ascID, sizeof(int), numV, fseg);
			fread(dscID, sizeof(int), numC, fseg);
			fclose(fseg);

			G_mscg_TEMP =
				new mscRegularGrid3DGradientField(bmsh, &a);


			classes = new mscPreClassifier(ascID, dscID, bmsh);

			// THIS IS A MISNOMER all it actually does is set 
			mscSimpleConstrainedGradientBuilder<float>* scgb =
				new mscSimpleConstrainedGradientBuilder<float>(classes, bmf, bmsh, G_mscg_TEMP, &a);
			scgb->computeGradient();


			mscModifiedBoundaryMeshHandler* mbmh =
				new mscModifiedBoundaryMeshHandler(bmsh, G_mscg_TEMP);

			mscSimpleGradientBuilder<float>* mscb3 =
				new mscSimpleGradientBuilder<float>(bmf, mbmh, bgf, &a);
			mscb3->computeGradient();

			//mscSimpleConstrainedGradientBuilder<float>* scgb2 = 
			//	new mscSimpleConstrainedGradientBuilder<float>(classes, bmf,mbmh,bgf, &a);
			//scgb2->computeGradient();


			//#define HACK

#ifdef HACK
			G_mscmh = mbmh;
#endif
		}
		else if (userand == 7) {
			printf("using assised gradient computation\n");
			mscAssistedGradientBuilder<float>* magb =
				new mscAssistedGradientBuilder<float>(bmf, bmsh, bgf, filename, &a);
			magb->computeGradient();




		}
		//else if (userand == 7)  {
		//	printf("using constrained convergent gradient - 2 pass\n");
		//	USE_SEG = true;
		//	FILE* fseg = fopen("source_dest2.raw", "rb");
		//	int numV = X * Y;
		//	int numC = (X-1) * (Y-1);

		//	int* ascID = new int[numV];
		//	int* dscID = new int[numC];
		//	fread(ascID, sizeof(int), numV, fseg);
		//	fread(dscID, sizeof(int), numC, fseg);
		//	fclose(fseg);

		//	G_mscg_TEMP = 
		//		new mscRegularGrid3DGradientField(bmsh, &a);

		//	
		//	classes = new mscPreClassifier(ascID, dscID, bmsh);

		//	// THIS IS A MISNOMER all it actually does is set 
		//	mscSimpleConstrainedGradientBuilder<float>* scgb = 
		//		new mscSimpleConstrainedGradientBuilder<float>(classes, bmf,bmsh,G_mscg_TEMP, &a);
		//	scgb->computeGradient();


		//	mscModifiedBoundaryMeshHandler* mbmh = 
		//		new mscModifiedBoundaryMeshHandler(bmsh, G_mscg_TEMP);

		//	mscConvergentGradientBuilder<float>* msccbX = 
		//		new mscConvergentGradientBuilder<float>(bmf, mbmh, bgf, &a);
		//	msccbX->computeGradient();


		//	mscComplementMeshHandler* cmh = 
		//		new mscComplementMeshHandler(mbmh);

		//	mscNegatingMeshFunction<float>* nmf = 
		//		new mscNegatingMeshFunction<float>(bmf);
		//	mscModifiedBoundaryMeshHandler* mbmh2 = 
		//		new mscModifiedBoundaryMeshHandler(cmh, bgf);
		//	mscRegularGrid3DGradientField* bgf2 = 
		//		new mscRegularGrid3DGradientField(bmsh, &a);
		//	mscConvergentGradientBuilder<float>* msccb2 = 
		//		new mscConvergentGradientBuilder<float>(nmf, mbmh2, bgf2, &a);
		//	msccb2->computeGradient();

		//	mscRegularGrid3DGradientField* tempgf = bgf;

		//	cellIterator it;
		//	iteratorOperator& all_cells = bmsh->all_cells_iterator(it);
		//	for (all_cells.begin(it); all_cells.valid(it); all_cells.advance(it)) {
		//		CELL_INDEX_TYPE cid = all_cells.value(it);
		//		bgf2->set_dim_asc_man(cid, bgf->get_dim_asc_man(cid));
		//	}
		//	bgf = bgf2;

		//	bgf->resetMeshHandler(bmsh);


		//	//char ufilename1[1024];
		//	//sprintf(ufilename1, "%s.asc", filename); 
		//	//test_data->dump_vals(ufilename1, X, Y, Z, msccb->getmaxvals());
		//	//sprintf(ufilename1, "%s.asc.prob", filename);

		//	//char ufilename2[1024];

		//	//sprintf(ufilename2, "%s.dsc", filename);
		//	//test_data->dump_vals(ufilename2, X, Y, Z, msccb2->getmaxvals());
		//	//		 sprintf(ufilename2, "%s.dsc.prob", filename);

		//	//char ufilename3[1024];

		//	//sprintf(ufilename3, "%s.both", filename);
		//	//FILE* ff1 = fopen(ufilename1, "rb");
		//	//FILE* ff2 = fopen(ufilename2, "rb");
		//	//FILE* ffout = fopen(ufilename3, "wb");
		//	//while (! feof(ff1)) {
		//	// float v1, v2;
		//	// fread(&v1, sizeof(float), 1, ff1);
		//	// fread(&v2, sizeof(float), 1, ff2);
		//	// float res = (0.5-0.5f*v1) + 0.5f - (0.5-0.5f*v2);
		//	// fwrite(&res, sizeof(float), 1, ffout);
		//	//}
		//	//fclose(ff1);
		//	//fclose(ff2);
		//	//fclose(ffout);

		//	G_msccb = msccb;
		//	G_msccb2 = msccb2;
		//	}
		bgf->output_to_file(gradname);

	}		
	G_mscg = bgf;

#ifdef HACK
		if(userand != 6) {
#endif
			G_mscmh = bmsh;
#ifdef HACK
		}
#endif
		G_mscf = bmf;



#ifdef WIN32        
		DWORD globalEnd = GetTickCount();
		printf(" --Computed discrete gradient in %.3f seconds\n", (globalEnd - globalStart) / 1000.0);
#endif




#ifndef HACK
	if (userand == 6) {
		MSC = new RestrictedMSC<float>(bgf, G_mscmh, bmf, G_mscg_TEMP);
		printf("using restricted\n");
	} else {
#endif
		MSC = new BasicMSC<float>(bgf, G_mscmh, bmf);
#ifndef HACK
	}
#endif


	MSC->ComputeFromGrad();
	MSC->ComputeHeirarchy();
	MSC->WriteCancelHistory();

	if (command_line_only) {

		MSC->setPercentPersistence(gpers);


		MSC->WriteComplex(argv[1]);

		return 1;

	}


	if (argc > 6)
		MSC->setAbsolutePersistence(gpers);
	else
		MSC->setPercentPersistence(5.0f);

	cellIterator ita;
	iteratorOperator& all = bmsh->all_cells_iterator( ita);
	//all.begin(ita); CELL_INDEX_TYPE cid = all.value(ita);
	fminmax.maxval = bmf->cell_value(0);
	fminmax.minval = bmf->cell_value(0);

	for (all.begin(ita); all.valid(ita); all.advance(ita)) {
		CELL_INDEX_TYPE cid = all.value(ita);
		float tmp = bmf->cell_value(cid);
		if (tmp > fminmax.maxval) fminmax.maxval = tmp;
		if (tmp < fminmax.minval) fminmax.minval = tmp;
	}

	printf("maxvalue = %f, minvalue = %f\n", fminmax.maxval, fminmax.minval);
	// const unsigned char purple[] = {255, 0, 255};
	// cimg_library::CImg<unsigned char>(640, 400, 1, 3, 0).draw_text(100, 100, "Hello World", purple).display("my first casasfd");


	/*
	int array[10];

	for (int i=0; i<10; i++) {
	array[i] = 2*i+5;
	}

	IndexIterator ii(array, 5);

	while (ii.isValid()) {
	cout << *ii.loc << endl;
	ii++;
	}*/



	Initialize_GLUT(argc, argv);

	return 0;
}


/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//////////                                       ////////////////
//////////    GLUT STUFF DOWN HERE               ////////////////
//////////                                       ////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

Matrix4 mat;
float translateoff[3];
bool g_mouseActive;
int mainWindow;

// Initial size of graphics window.
const int WIDTH  = 1200;
const int HEIGHT = 1200;

// Current size of window.
int width  = WIDTH;
int height = HEIGHT;

// Mouse positions, normalized to [0,1].
double xMouse = 0.0;
double yMouse = 0.0;

// Bounds of viewing frustum.
double nearPlane =  0.1;
double farPlane  = 30000;

// Viewing angle.
double fovy = 40.0;

// Variables.
double alpha = 0;                                  // Set by idle function.
double beta = 0;                                   // Set by mouse X.
double mdistance = - (farPlane - nearPlane) / 2;    // Set by mouse Y.
float dist;



GLfloat light_diffuse[] = {0.8, 0.8, 0.8, 1.0};  /* Red diffuse light. */
GLfloat light_position[] = {1.0, 1.0, 1.0, 0.0};  /* Infinite light location. */

float diffuseLight[] = {0.8f, 0.8f, 0.8f, 1.0f};
float specularLight[] = {1.0f, 1.0f, 1.0f, 1.0f};
float LightPosition[] = {1.1f, 0.0f, 8.0f, 1.0f};

struct Quat {
	float x, y, z, w;
};

Quat times(Quat a, Quat b) {
	Quat res;
	res.w = a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z;
	res.x = a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y;
	res.y = a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x;
	res.z = a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w;
	return res;
}

Quat rotation(Vector4 v, float angle) {
	Quat res;
	Normalize3(&v);
	res.w = cosf(angle/2.0);
	float tmp = sinf(angle/2.0);
	res.x = v.vals[0] * tmp;
	res.y = v.vals[1] * tmp;
	res.z = v.vals[2] * tmp;
	return res;
}

Matrix4 rotmatrix;
Quat rotationtotal;
void resetQuat(Quat &q) {
	q.x = 0;
	q.y = 0;
	q.z = 0;
	q.w = 1;
}

void setRotQuat(Matrix4 &m, Quat &q){
	{
		float x2 = 2.0f * q.x,  y2 = 2.0f * q.y,  z2 = 2.0f * q.z;

		float xy = x2 * q.y,  xz = x2 * q.z;
		float yy = y2 * q.y,  yw = y2 * q.w;
		float zw = z2 * q.w,  zz = z2 * q.z;

		m.vals[ 0] = 1.0f - ( yy + zz );
		m.vals[ 1] = ( xy - zw );
		m.vals[ 2] = ( xz + yw );
		m.vals[ 3] = 0.0f;

		float xx = x2 * q.x,  xw = x2 * q.w,  yz = y2 * q.z;

		m.vals[ 4] = ( xy +  zw );
		m.vals[ 5] = 1.0f - ( xx + zz );
		m.vals[ 6] = ( yz - xw );
		m.vals[ 7] = 0.0f;

		m.vals[ 8] = ( xz - yw );
		m.vals[ 9] = ( yz + xw );
		m.vals[10] = 1.0f - ( xx + yy );  
		m.vals[11] = 0.0f;  

		m.vals[12] = 0.0f;  
		m.vals[13] = 0.0f;   
		m.vals[14] = 0.0f;   
		m.vals[15] = 1.0f;
	}

}

void reset_view() {
	resetQuat(rotationtotal);
	mat = MIdentity();
	dist = 140*(XMAX-XMIN);
}

void render_view (){
	resetQuat(rotationtotal);
	mat = MIdentity();
	dist = 240*(ZMAX-ZMIN);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//glMultMatrixf(mat.vals);
	//glRotatef(-50, 0,1,0);
	Vector4 v; v.vals[0] = 0; v.vals[1] = 1; v.vals[2] = 0; 
	rotationtotal = rotation(v, -50*PI/180.0);
	setRotQuat(rotmatrix, rotationtotal);
	glMultMatrixf(rotmatrix.vals);

	//glMultMatrixf(mat.vals);

	//glGetFloatv(GL_MODELVIEW_MATRIX, mat.vals);
}

float read_p_from_file() {
	FILE* fin = fopen("input_persistence.txt", "r");
	float f; 
	fscanf(fin, "%f\n", &f);
	//printf("read %f\n", f);
	fclose(fin);
	return f;
}


float read_t_from_file() {
	FILE* fin = fopen("input_value.txt", "r");
	float f;
	fscanf(fin, "%f\n", &f);
	//printf("read %f\n", f);
	fclose(fin);

	return f * fminmax.maxval + fminmax.minval;
}

float oldf = -1.0;
void myidle() {
	//	cimg_library::cimg::sleep(10);
#ifdef WIN32
   Sleep(50);
#else
   usleep(50);
#endif
	float newf = read_p_from_file();
	if (newf != oldf) {
		oldf = newf;

		MSC->setPercentPersistence(oldf);
		redrawstuff = true;
		glutPostRedisplay();
	}
	if (g_use_test) {
		float ttestval = read_t_from_file();
		if (ttestval != myval) {
			myval = ttestval;
			redrawstuff = true;
			glutPostRedisplay();
		}
	}
}

bool dbbox = true;

void display()
{
	glPointSize(4.0);
	if (clearcolor) {
		glClearColor(1,1,1,1);
	} else {
		glClearColor(0,0,0,1);
	}
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// enable crystal ball transform here

	glTranslatef( 0,0,- dist / 100);
	glTranslatef(- translateoff[0], - translateoff[1], - translateoff[2]);

	setRotQuat(rotmatrix, rotationtotal);
	glMultMatrixf(rotmatrix.vals);


	glTranslatef(-(XMAX - XMIN)/2.0 - XMIN,
		-(YMAX - YMIN)/2.0 - YMIN,
		-(ZMAX - ZMIN)/2.0 - ZMIN);

	glEnable(GL_DEPTH_TEST);


	// draw the complex

	// if (dbbox) {
	// 	   if (clearcolor)
	// 		glColor3f(.05, .05, .05);
	// 	   else 
	// 		   glColor3f(.9,.9,.9);
	// 	 glBegin(GL_LINES);
	// 	 glVertex3f(0, 0, 0);
	// 	 glVertex3f(XDIM-1 , 0, 0);
	// 	 glVertex3f(0, YDIM-1, 0);
	// 	 glVertex3f(XDIM-1 , YDIM-1, 0);
	// 	 glVertex3f(0, 0, ZDIM-1);
	// 	 glVertex3f(XDIM-1 , 0, ZDIM-1);
	// 	 glVertex3f(0, YDIM-1, ZDIM-1);
	// 	 glVertex3f(XDIM-1 , YDIM-1, ZDIM-1);

	// 	glVertex3f(0, 0, 0);
	// 	glVertex3f(0, YDIM-1, 0);
	// 	glVertex3f(XDIM-1, 0, 0);
	// 	glVertex3f(XDIM-1, YDIM-1, 0);
	// 	glVertex3f(0, 0, ZDIM-1);
	// 	glVertex3f(0, YDIM-1, ZDIM-1);
	// 	glVertex3f(XDIM-1, 0, ZDIM-1);
	// 	glVertex3f(XDIM-1, YDIM-1, ZDIM-1);

	// 	 glVertex3f(0, 0, 0);
	// 	glVertex3f(0, 0, ZDIM-1);
	// 		 glVertex3f(XDIM-1, 0, 0);
	// 	glVertex3f(XDIM-1, 0, ZDIM-1);
	// 		 glVertex3f(0, YDIM-1, 0);
	// 	glVertex3f(0, YDIM-1, ZDIM-1);
	// 		 glVertex3f(XDIM-1, YDIM-1, 0);
	// 	glVertex3f(XDIM-1, YDIM-1, ZDIM-1);

	// 	glEnd();
	// }

	drawStuff<2, float>();
	glutSwapBuffers();


}



void *gl_copy_pixels = NULL;



void reshapeMainWindow (int newWidth, int newHeight)
{
	width = newWidth;
	height = newHeight;
	if (gl_copy_pixels != NULL) {
		delete(gl_copy_pixels);
		gl_copy_pixels = NULL;
	}
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(fovy, GLfloat(width) / GLfloat(height), nearPlane, farPlane);
	glutPostRedisplay();
}

int gcounter;

void outputimage() {
	unsigned char *rpixels = new unsigned char[width*height*3];
	unsigned char *mpixels = new unsigned char[width*height*3];
	char name[1024];
	sprintf(name, "im%d_%d.bmp", gcounter++, seed);	
	glReadPixels(0,0,width,height, GL_RGB, GL_UNSIGNED_BYTE, rpixels);
	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			int id1 = x + y*width;
			int id2 = x + (height - 1 - y)*width;
			mpixels[id1] = rpixels[id2*3];
			mpixels[id1 +width*height] = rpixels[id2*3+1];
			mpixels[id1 +2*width*height] = rpixels[id2*3+2];
		}
	}
#ifdef WIN32
	CImg< unsigned char > img(mpixels, width, height,1,  3, false);

	img.save_bmp(name);
#endif
	delete(rpixels);
	delete(mpixels);
	return;

}

void makeSequence() {

	//for (int i = 0; i < total_order.size(); i++) {
	//  usc->setAssigned(total_order[i], false);
	//}
	//for (int i = 0; i < total_pqlist.size(); i++) {
	//  usc->setNumUFacets(total_pqlist[i], 0);
	//}
	//int nc = total_counts[0];
	//for (int i = 0; i < nc; i++) {
	//  usc->setNumUFacets(total_pqlist[i], 1);
	//}

	//char fname[1024];

	//int counter = 0;
	//for (int i = 0; i < 700/*total_order.size()*/; i++) {
	// sprintf(fname, "out/full_%05d.png", counter++);
	// printf("outputtting: %s\n", fname);


	// if (usc->getCritical(total_order[i])) {
	//  for (int j = nc; j < nc+total_counts[i+1]; j++) {
	//	  usc->setNumUFacets(total_pqlist[j], 1);
	//  }
	//  nc += total_counts[i+1];

	//  usc->setAssigned(total_order[i], true);
	//  usc->setNumUFacets(total_order[i], 0);
	//  glutPostRedisplay();  
	//  display();
	//  glutSwapBuffers();
	// // outputimage(fname);

	// } else {
	//  for (int j = nc; j < nc+total_counts[i+1]; j++) {
	//	  usc->setNumUFacets(total_pqlist[j], 1);
	//  }
	//  nc += total_counts[i+1];

	//  usc->setAssigned(total_order[i++], true);

	//  for (int j = nc-1; j < nc+total_counts[i]; j++) {
	//	  usc->setNumUFacets(total_pqlist[j], 1);
	//  }
	//  nc += total_counts[i];
	//  usc->setAssigned(total_order[i], true);

	//  display();
	//  glutSwapBuffers();
	//  //glutPostRedisplay();  
	// // outputimage(fname);
	// }
	//}  
}



void graphicKeys (unsigned char key, int x, int y)
{

	switch (key)
	{

	case 'o':
		outputimage();
		break;
	case 't':
		draw_mt = ! draw_mt;
		redrawstuff = true;
		//bsn->dumpBasinCountPerPersistence("basins.txt");
		break;

	case 'y':
		char tmpname[1024];
		//sprintf(tmpname, "bsizes_%.6u.txt", glivingcounter);
		//bsn->dumpBasinSizes(tmpname, glivingcounter);
		break;

	case 'h':
		gusecutoffhack = ! gusecutoffhack;
		redrawstuff = true;
		break;

	case '1':
		MSC->setPercentPersistence(0);
		redrawstuff = true;
		break;
	case '2':
		MSC->setFirstNonzeroPersistence(2);
		redrawstuff = true;
		break;
	case '3':
		MSC->setPercentPersistence(2);

		redrawstuff = true;
		break;
	case '4':
		MSC->addPersistence(1);

		//	  glivingcounter = bsn->getDestrCount(0.01f);
		redrawstuff = true;
		break;
	case '5':
		MSC->addPersistence(-1);

		//	  glivingcounter = bsn->getDestrCount(0.05f);
		redrawstuff = true;
		break;
	case '6':
		MSC->setPercentPersistence(1);

		redrawstuff = true;
		break;
	case '7':
		MSC->setPercentPersistence(6);
		redrawstuff = true;
		break;
	case '8':
		MSC->setPercentPersistence(8);
		redrawstuff = true;
		break;
	case '9':
		MSC->setPercentPersistence(99);
		redrawstuff = true;
		break;

	case 'c':
		//		drawbasincolor = !drawbasincolor;
		redrawstuff = true;
		break;
	case 'j':
		//		if (glivingcounter > 0)
		//			glivingcounter--;	
		redrawstuff = true;
		break;	
	case 'k':
		//     if (glivingcounter < bsn->num_destroyed()) 
		//  		glivingcounter++;
		redrawstuff = true;
		break;
	case 'A':
		DRAW_BACKGROUND2 = !DRAW_BACKGROUND2;
		redrawstuff = true;
		break;
	case 'D':
		DRAW_BACKGROUND3 = !DRAW_BACKGROUND3;
		redrawstuff = true;
		break;
	case 'b':
		DRAW_BACKGROUND = !DRAW_BACKGROUND;
		redrawstuff = true;
		break;
	case 'B':
		clearcolor = !clearcolor;
		redrawstuff = true;
		break;

	case 'm':
		{
			/*         for (int i=0; i<1000; i++) {
			if (cutoff < total_order.size()) {
			if (usc->getCritical(total_order[cutoff])){
			usc->setAssigned(total_order[cutoff], true);
			if (cutoff < total_order.size()-1) cutoff++;
			} else {
			usc->setAssigned(total_order[cutoff], true);
			cutoff++;
			usc->setAssigned(total_order[cutoff], true);
			cutoff++;
			}
			}
			}
			*/     }
		redrawstuff = true;
		break;
	case 'n':
		DRAWPOINTS = !DRAWPOINTS;
		redrawstuff = true;
		break;	
	case 'M':
		{
			//   for (int i=0; i<1; i++) {
			//     if (cutoff < total_order.size()) {
			//if (usc->getCritical(total_order[cutoff])){
			//  usc->setAssigned(total_order[cutoff], true);
			//  if (cutoff < total_order.size()-1) cutoff++;
			//} else {
			//  usc->setAssigned(total_order[cutoff], true);
			//  cutoff++;
			//  usc->setAssigned(total_order[cutoff], true);
			//  cutoff++;
			//      }	       
			//     }
			//   }
		}
		redrawstuff = true;
		break;


	case 's':
		makeSequence();
		break;

	case 'i':
		g_draw_isolines = !g_draw_isolines;
		redrawstuff = true;
		break;
	case 'T':
		g_use_test = !g_use_test;
		redrawstuff = true;
		break;
	case 'p':
		DRAW_PROBABILITIES1 = !DRAW_PROBABILITIES1;
		redrawstuff = true;
		break;
	case 'P':
		DRAW_PROBABILITIES2 = !DRAW_PROBABILITIES2;
		redrawstuff = true;
		break;

	case '[':
		line_width *= 1.2f;
		redrawstuff = true;
		break;
	case ']': 
		line_width /= 1.2f;
		redrawstuff = true;
		break;
	case '{':
		point_size *= 1.2f;
		redrawstuff = true;
		break;
	case '}': 
		point_size /= 1.2f;
		redrawstuff = true;
		break;
	case 'e':
		draw_edges = !draw_edges;
		redrawstuff = true;
		break;
	case 'F':
		flat_funct = !flat_funct;
		redrawstuff = true;
		break;
	case 'f':
		draw_flat = !draw_flat;
		redrawstuff = true;
		break;
	case 'G':
		DRAW_GRID = !DRAW_GRID;
		redrawstuff = true;
		break;
	case 'K':
		DRAWNUMERICALDE = !DRAWNUMERICALDE;
		redrawstuff = true;
		break;
	case 'N':
		USE_SEG = !USE_SEG;
		redrawstuff = true;
		break;
	case '<':
		USE_SEG_LINES = !USE_SEG_LINES;
		redrawstuff = true;
		break;
	case '>':
		USE_SEG_AMAP = !USE_SEG_AMAP;
		redrawstuff = true;
		break;
	case '?':
		USE_SEG_DMAP = !USE_SEG_DMAP;
		redrawstuff = true;
		break;
	case 'g':
		draw_gradient = !draw_gradient;
		redrawstuff = true;
		break;

	case '(':
		arrow_width *= 1.2f;
		redrawstuff = true;
		break;

	case ')':
		arrow_width /= 1.2f;
		redrawstuff = true;
		break;

	case 'x':
		ballsize *= 1.2f;
		redrawstuff = true;
		break;

	case 'z':
		ballsize /= 1.2f;
		redrawstuff = true;
		break;
		//case 'w': 
		//  {
		//    char filename[1024];
		//    sprintf(filename, "out.bmp");
		//    outputimage(filename);

		//    break;
		//  }
	case 'L':
		DumpLines<2,float>();
		break;

	case 'a':
		DRAWASCLINES = ! DRAWASCLINES;
		redrawstuff = true;
		break;
	case 'd':
		DRAWDSCLINES = ! DRAWDSCLINES;
		redrawstuff = true;
		break;
	case 'r':
		reset_view();
		break;
	case 'q':
		exit(1);
		break;
	}
	glutPostRedisplay();
}

int xold;
int yold;
Vector4 vstart;
int buttonstate;
void mouseClicked (int button, int state, int x, int y)
{
	xold = x;
	yold = y;

	//printf("Mouse Clicked! button = %d, %d  x = %d   y = %d\n", state, button, x, y);
	buttonstate= button;
	if (state == GLUT_DOWN) {
		g_mouseActive = true;
	} else {
		g_mouseActive = false;
	}

	glutPostRedisplay();
}


void mouseMovement (int mx, int my)
{
	//printf("%d, %d\n", mx, my);
	if (! g_mouseActive) return;
	if (buttonstate == 2) {
		int dd = XMAX-XMIN;
		dist += dd*(my - yold);
		yold = my;

	}

	if (buttonstate == 1 || buttonstate == 0) {

		float ldim = (width > height ? width : height);


		vstart.vals[0] = xold - .5 * width ;
		vstart.vals[1] = .5 * height - yold;

		if ((.5 * ldim) * (.5* ldim) - vstart.vals[0]*vstart.vals[0] - vstart.vals[1]*vstart.vals[1] <= 0) return;
		vstart.vals[2] = sqrt((.5 * ldim) * (.5* ldim) - vstart.vals[0]*vstart.vals[0] - vstart.vals[1]*vstart.vals[1]);
		Normalize3(&vstart);

		xold = mx;
		yold = my;

		// calculate the vectors;
		Vector4 vend;


		vend.vals[0] = mx - .5 * width ;
		vend.vals[1] = .5 * height - my;

		if ((.5 * ldim) * (.5* ldim) - vend.vals[0]*vend.vals[0] - vend.vals[1]*vend.vals[1] <= 0) return;
		vend.vals[2] = sqrt((.5 * ldim) * (.5* ldim) - vend.vals[0]*vend.vals[0] - vend.vals[1]*vend.vals[1]);
		Normalize3(&vend);

		translateoff[0] += -(XMAX-XMIN)* (vend.vals[0] - vstart.vals[0]);
		translateoff[1] += -(XMAX-XMIN)* (vend.vals[1] - vstart.vals[1]);
		translateoff[2] += -0* (vend.vals[2] - vstart.vals[2]);
		glutPostRedisplay();
	}

	//if (buttonstate == 0) {
	//	float ldim = (width > height ? width : height);

	//	if (xold == mx && yold == my) return;

	//	vstart.vals[0] = xold - .5 * width ;
	//	vstart.vals[1] = .5 * height - yold;

	//	if ((.5 * ldim) * (.5* ldim) - vstart.vals[0]*vstart.vals[0] - vstart.vals[1]*vstart.vals[1] < 0) return;
	//	vstart.vals[2] = sqrt((.5 * ldim) * (.5* ldim) - vstart.vals[0]*vstart.vals[0] - vstart.vals[1]*vstart.vals[1]);
	//	if (Normalize3(&vstart)== 0) return;
	//	xold = mx;
	//	yold = my;

	//	// calculate the vectors;
	//	Vector4 vend;


	//	vend.vals[0] = mx - .5 * width ;
	//	vend.vals[1] = .5 * height - my;

	//	if ((.5 * ldim) * (.5* ldim) - vend.vals[0]*vend.vals[0] - vend.vals[1]*vend.vals[1] <= 0) return;
	//	vend.vals[2] = sqrt((.5 * ldim) * (.5* ldim) - vend.vals[0]*vend.vals[0] - vend.vals[1]*vend.vals[1]);
	//	if (Normalize3(&vend)== 0) return;



	//	float alpha = InteriorAngle(vstart, vend);
	//	if (alpha < 0.01) return;

	//	Vector4 cp = Cross(vstart, vend);
	//	if (Normalize3(&cp) == 0) return;

	//	//update the crystal ball matrix
	//	glMatrixMode(GL_MODELVIEW);
	//	glLoadIdentity();

	//	rotationtotal = times( rotationtotal,rotation(cp, alpha*PI/-180.0));

	//	setRotQuat(rotmatrix, rotationtotal);
	//	glMultMatrixf(rotmatrix.vals);

	//}


	glutPostRedisplay();

}

void initScene() {
	glDepthMask(true);
	if (clearcolor) {
		glClearColor(1,1,1,1);
	} else {
		glClearColor(0,0,0,1);
	}

	glDisable(GL_CULL_FACE);
	glDepthFunc(GL_LESS);
	glEnable(GL_DEPTH_TEST);

	glShadeModel(GL_SMOOTH);
	glLineWidth(1);
	glEnable(GL_LINE_SMOOTH);


	// Enable our light.

	glEnable(GL_LIGHT0);

	// Set up the material information for our objects.  Again this is just for show.
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specularLight);
	glMateriali(GL_FRONT_AND_BACK, GL_SHININESS, 128);


	// initialize crystal ball
	mat = MIdentity();
	dist = 300*(ZMAX-ZMIN);

}



void Initialize_GLUT(int argc, char **argv) {


	translateoff[0] =0.0;
	translateoff[1] =0.0;
	translateoff[2] =0.0;
	resetQuat(rotationtotal);
	g_mouseActive = false;

	// GLUT initialization.
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(width, height);
	glutInitWindowPosition(570, 80);
	mainWindow = glutCreateWindow("Morse3d Efficient Computation");

	//		glewInit();

	glutDisplayFunc(display);
	glutReshapeFunc(reshapeMainWindow);
	glutKeyboardFunc(graphicKeys);
	glutMouseFunc    (mouseClicked);
	//glutSpecialFunc(functionKeys);
	glutMotionFunc(mouseMovement);
	glutIdleFunc(myidle);

	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	initScene();
	//		setGlobalMinMax(usc);



	glutMainLoop();
}
