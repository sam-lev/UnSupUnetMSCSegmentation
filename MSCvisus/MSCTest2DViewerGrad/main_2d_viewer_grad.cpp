


#include <iostream>
#include "dump_image.h"
#include <windows.h>


#include <stdio.h>
	#include <GL/glu.h>  // openGL utilities
	#include <GL/gl.h>   // openGL declarations
#include "freeglut.h"
#include "Matrix4.h"
#include "C:\local\CImg-1.3.2\CImg.h"
using namespace cimg_library;

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
#include "mscContinuousGradient.h"

#include "interpolation.h"
using namespace alglib_impl;

#define EPSILON 0.001



mscContinuousGradient2D* cg_msc;


mscVectorToScalarGridFunc2D* v2d_msc;


BasicMSC<float>* MSC;
mscBasicGradientField* G_mscg;
mscBasicMeshHandler* G_mscmh;
mscBasicMeshFunction<float>* G_mscf;
mscConvergentGradientBuilder<float>* G_msccb;
mscConvergentGradientBuilder<float>* G_msccb2;

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



int XMIN = 2;
int YMIN = 2;
int ZMIN = 2;

int OLDX = 0;
int OLDY = 0;

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

unsigned char*	setup_LIC(double* VVV);

int LIC_X = 700;
int LIC_Y = 700;
unsigned char* GLOBAL_LIC_IMAGE;
  
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
bool DRAWASCLINES = true;
bool DRAWDSCLINES = true;
bool DRAWPOINTS = true;
bool DRAWDIMASCMAN = false;
bool DRAWGRAD = false;
bool DRAW_PROBABILITIES1 = false;
bool DRAW_PROBABILITIES2 = false;
bool DRAW_GRID = false;
bool draw_mt = false;
bool flat_funct = false;
bool redrawstuff = true;
GLuint drawlist = -1;


float GDEGREE;
float deg;
int oX, oY, oZ;
float* fbbb;

float MANALPHA = .3f;


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
	//if(DRAW_BACKGROUND) {

	//glBegin(GL_QUADS);
	//iteratorOperator& all = G_mscmh->all_cells_iterator(it);
	//for (all.begin(it); all.valid(it); all.advance(it)) {
	//	//if (G_mscg->get_dim_asc_man(all.value(it)) != 3) 
	//	CELL_INDEX_TYPE cid = all.value(it);

	//	float sc = G_mscf->cell_value(cid);
	//	sc = (sc - fminmax.minval) / (fminmax.maxval - fminmax.minval);
	//	 sc = sc*.8 + .2;
	//	glColor3f(sc+.01, sc, sc+.01);

	//	float c[3]; 
	//	coordinates(cid, c);
	//	glVertex3f(c[0]-0.5, c[1]-0.5, c[2]);
	//	glVertex3f(c[0]-0.5, c[1]+0.5, c[2]);
	//	glVertex3f(c[0]+0.5, c[1]+0.5, c[2]);
	//	glVertex3f(c[0]+0.5, c[1]-0.5, c[2]);
	//	//glVertex3f(c[0], c[1], c[2]);


	//}
	//glEnd();



	//}

	if(DRAW_BACKGROUND) {

		//glPushMatrix();
		//glScalef((float) OLDX/ (float) LIC_X, (float) OLDY / (float) LIC_Y, 1.0);
	glBegin(GL_QUADS);
	for (int i = 0; i < LIC_X; i++) {
		for (int j = 0; j < LIC_Y; j++) {

			unsigned char value = GLOBAL_LIC_IMAGE[i + j*LIC_X];
			float fval = (float) value / 255.0f; 
			glColor3f(fval, fval, fval);
			float c[2];
			float xfact = (float) OLDX/ (float) LIC_X;
			float yfact = (float) OLDY / (float) LIC_Y;
			c[0] = i*2 * xfact;
			c[1] = j*2 * yfact;

		glVertex3f(c[0]-xfact, c[1]-yfact,0);
		glVertex3f(c[0]-xfact, c[1]+yfact, 0);
		glVertex3f(c[0]+xfact, c[1]+yfact, 0);
		glVertex3f(c[0]+xfact, c[1]-yfact, 0);
		//glVertex3f(c[0], c[1], c[2]);
		}

	}
	glEnd();
	//glPopMatrix();


	}




	glDisable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	
	// draw grid


	if (DRAW_BACKGROUND) {
		glBegin(GL_QUADS);

		map<CELL_INDEX_TYPE, node<float>*>::iterator nit = MSC->nodes.begin();
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



			float centroid[2];
	centroid[0] = (float) oX/2.0f;
	centroid[1] = (float) oY/2.0f;
	float centroid2[2];
	centroid2[0] = (float) ((XMAX-1)/2)/2.0f;
	centroid2[1] = (float) ((YMAX-1)/2)/2.0f;
	
			float tx = (float) col[0]/2.0f - centroid2[0]; 
			float ty = (float) col[1]/2.0f - centroid2[1]; 

			float xnew = tx * cos(deg) - ty * sin(deg);
			float ynew = ty * cos(deg) + tx * sin(deg);

			xnew += centroid[0];
			ynew += centroid[1];



		int dd = 60;
		float hd = 60.0f;
		
		float v1 = (sin(PI*((int) xnew%dd)/((float) hd)))/2.0+.4f;
		float v2 = (sin((dd - (((int) ynew)%dd))/hd))/2.0+.4f;
		float v3 = /*(sin(PI*((((int) col[0])%dd)-((int) col[1])%dd)/hd))/2.0+*/.4f;
		
		glColor4f(v1, v2, v3, 0/*MANALPHA*/);
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


	nit = MSC->nodes.begin();
	while (nit != MSC->nodes.end()) {
		node<float>* n = (*nit).second;
		nit++;
		if (! MSC->isAlive(n)) continue;
		if (n->index != 0 ) continue;

		set<CELL_INDEX_TYPE> g;
		MSC->fillGeometry(n, g);
		float col[3];
		coordinates(n->cellid, col);

			float centroid[2];
	centroid[0] = (float) oX/2.0f;
	centroid[1] = (float) oY/2.0f;
	float centroid2[2];
	centroid2[0] = (float) ((XMAX-1)/2)/2.0f;
	centroid2[1] = (float) ((YMAX-1)/2)/2.0f;
	
			float tx = (float) col[0]/2.0f - centroid2[0]; 
			float ty = (float) col[1]/2.0f - centroid2[1]; 

			float xnew = tx * cos(deg) - ty * sin(deg);
			float ynew = ty * cos(deg) + tx * sin(deg);

			xnew += centroid[0];
			ynew += centroid[1];




		int dd = 10;
		float hd = 10.0f;
	
		float v1 = (sin(PI*((int) xnew%dd)/((float) hd)))/2.0+.4f;
		float v2 = (sin((dd - (((int) ynew)%dd))/hd))/2.0+.4f;
		float v3 = /*(sin(PI*((((int) col[0])%dd)-((int) col[1])%dd)/hd))/2.0+*/.4f;
		float alpha = MANALPHA;
		//if ( n->boundary) {
		//	v1 = 0.3;
		//	v2 = 0.3;
		//	v3 = 0.3;
		//	//alpha = 0.95;
		//}
		glColor4f(v3, v2, v1, alpha);
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



	glEnd();


	}








	if (DRAW_GRID) {

	glBegin(GL_LINES);
	glColor4f(0,0,0, .2);

	for (int i = 0; i < XMAX; i+=2) {

		glVertex3f(i, 0, 0);
		glVertex3f(i, YMAX, 0);
	}
	for (int i = 0; i < YMAX; i+=2) {

		glVertex3f(0, i, 0);
		glVertex3f(XMAX, i, 0);
	}

	glEnd();















	}
	
	
	//////glBegin(GL_LINES);
	//////
	//////vector<float> isovals;
	//////for (int i = 0; i < 20; i++) {
	//////	isovals.push_back(fminmax.minval + (i/ 20.0f) * (fminmax.maxval - fminmax.minval));
	//////}
	//////
	//////
	//////glColor4f(0,0,0, .2);

	//////iteratorOperator& mall = G_mscmh->d_cells_iterator(2, it);
	//////for (mall.begin(it); mall.valid(it); mall.advance(it)) {
	//////	//if (G_mscg->get_dim_asc_man(all.value(it)) != 3) 
	//////	CELL_INDEX_TYPE cid = mall.value(it);

	//////	for (int i = 0; i < isovals.size(); i++) {
	//////		cellIterator fit;
	//////		iteratorOperator& facets = G_mscmh->facets(cid, fit);
	//////		for (facets.begin(fit); facets.valid(fit); facets.advance(fit)) {
	//////			CELL_INDEX_TYPE eid = facets.value(fit);
	//////		
	//////			if (G_mscf->cell_value(eid) >= isovals[i] && 
	//////				G_msccb->lowest_facet_value(eid) < isovals[i]) {

	//////				float c[2][3]; 
	//////				float vf[2];
	//////				int vert = 0;
	//////				cellIterator vit;
	//////				iteratorOperator& viter = G_mscmh->facets(eid,vit);
	//////				for (viter.begin(vit); viter.valid(vit); viter.advance(vit)) {
	//////						CELL_INDEX_TYPE vid = viter.value(vit);
	//////						vf[vert]= G_mscf->cell_value(vid);
	//////						coordinates(vid, c[vert++]);
	//////				}
	//////				
	//////				float t = (isovals[i] - vf[0]) / (vf[1] - vf[0]);
	//////				glVertex3f(c[0][0] + t*(c[1][0] - c[0][0]),
	//////					c[0][1] + t*(c[1][1] - c[0][1]),
	//////					c[0][2] + t*(c[1][2] - c[0][2]));
	//////				

	//////			}
	//////		}
	//////	}
	//////}


	//////
	//////glEnd();
	
	//	cellIterator it2;
	//iteratorOperator& ait = G_mscmh->all_cells_iterator(it);
	//for (ait.begin(it); ait.valid(it); ait.advance(it)) {
	//	CELL_INDEX_TYPE cid = ait.value(it);
	//	
	//	//if (G_mscmh->dimension(cid) == 1) { //->get_dim_asc_man(cid) != 2)
	//		drawArrow(cid);
	//	//}

	//}
	
	
	
	
	glBegin(GL_QUADS);

#ifdef FANCY_PROB_OUTPUT
	if (DRAW_PROBABILITIES1) {
	
	vector<idfpair>& v = G_msccb->getmaxvals();
	for (int i = 0; i < v.size(); i++) {
		idfpair p = v[i];
		if (G_mscmh->dimension(p.id) != 0) continue;
		if (p.prob < .98) {
			glColor4f(1, 0,0, 0.5*(1-p.prob));
		float c[3]; 
		coordinates(p.id, c);
		glVertex3f(c[0]-1.0, c[1]-1.0, c[2]);
		glVertex3f(c[0]-1.0, c[1]+1.0, c[2]);
		glVertex3f(c[0]+1.0, c[1]+1.0, c[2]);
		glVertex3f(c[0]+1.0, c[1]-1.0, c[2]);
		}

	}
	}

	if (DRAW_PROBABILITIES2) {
	
	if (G_msccb != G_msccb2) {
		vector<idfpair>& v2 = G_msccb2->getmaxvals();
		for (int i = 0; i < v2.size(); i++) {
			idfpair p = v2[i];
		if (G_mscmh->dimension(p.id) != 2) continue;
			if (p.prob < .98) {
				glColor4f(0, 0,1, 0.5*(1-p.prob));
				float c[3]; 
				coordinates(p.id, c);
				glVertex3f(c[0]-1.0, c[1]-1.0, c[2]);
				glVertex3f(c[0]-1.0, c[1]+1.0, c[2]);
				glVertex3f(c[0]+1.0, c[1]+1.0, c[2]);
				glVertex3f(c[0]+1.0, c[1]-1.0, c[2]);
			}

		}


	}
	}
	glEnd();

#endif


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

	glLineWidth(line_width*1.1);

	for (int i = 0; i < MSC->arcs.size(); i++) {
		arc<float>* a = MSC->arcs[i];
		if (! MSC->isAlive(a)) continue;
		if (a->lower->index == 0 && ! DRAWDSCLINES) continue;
		if (a->lower->index == 1 && ! DRAWASCLINES) continue;

		if (a->lower->index == 0) {
			glColor4f(.1, .2, .5,.8);
			//printf("ASDFL:KJSDFL:KJSDKL:FJSD\n");
		} else {
			glColor4f(.5, .2, .1,.8);
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
	
	glLineWidth(line_width*0.5);
	for (int i = 0; i < MSC->arcs.size(); i++) {
		arc<float>* a = MSC->arcs[i];
		if (! MSC->isAlive(a)) continue;
		
		if (a->lower->index == 0 && ! DRAWDSCLINES) continue;
		if (a->lower->index == 1 && ! DRAWASCLINES) continue;
		
		if (a->lower->index == 0) {
			glColor4f(.6, .6, .9,.5);
			//printf("ASDFL:KJSDFL:KJSDKL:FJSD\n");
		} else {
			glColor4f(.9, .6, .6,.5);
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
	glBegin(GL_POINTS);

	if (DRAWPOINTS){
	map<CELL_INDEX_TYPE, node<float>*>::iterator nit = MSC->nodes.begin();
	while (nit != MSC->nodes.end()) {
		node<float>* n = (*nit).second;
		nit++;
		if (! MSC->isAlive(n)) continue;
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
		if (! MSC->isAlive(n)) continue;
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





	////glLineWidth(line_width*0.5);
	////glBegin(GL_LINES);
	////for (int i = 0; i < v2d_msc->barfedges.size(); i+= 2) {
	////	float c[3]; 
	////	
	////	glColor4f(1.0, 0, 0, .5);
	////	CELL_INDEX_TYPE tid = v2d_msc->barfedges[i];
	////	glVertex3f(2*(tid % OLDX), 2*((tid / OLDX) % OLDY), 0);
	////	
	////	glColor4f(0.0, 1.0, 0, .5);
	////	tid = v2d_msc->barfedges[i+1];
	////	glVertex3f(2*(tid % OLDX), 2*((tid / OLDX) % OLDY), 0);
	////	
	////}
	////glEnd();

	////glPointSize(point_size*2.75);
	////glBegin(GL_POINTS);
	////for (int i = 0; i < v2d_msc->barfproc.size(); i++) {
	////	glColor4f(1.0, 1.00, 0, .85);
	////	CELL_INDEX_TYPE tid = v2d_msc->barfproc[i];
	////	glVertex3f(2*(tid % OLDX), 2*((tid / OLDX) % OLDY), 0);
	////}
	////glEnd();


	////	glPointSize(point_size*4.75);
	////glBegin(GL_POINTS);
	////for (int i = 0; i < v2d_msc->barfpoop.size(); i+=2) {
	////	glColor4f(0.0, 1.00, 1.0, .85);
	////	CELL_INDEX_TYPE tid = v2d_msc->barfpoop[i];
	////	glVertex3f(2*(tid % OLDX), 2*((tid / OLDX) % OLDY), 0);
	////	glColor4f(1.0, 0.00, 1.0, .85);
	////	tid = v2d_msc->barfpoop[i+1];
	////	glVertex3f(2*(tid % OLDX), 2*((tid / OLDX) % OLDY), 0);
	////}
	////glEnd();










	glDisable(GL_BLEND);

		glEnable(GL_DEPTH_TEST);







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
//   if (draw_gradient) {
////      glLineWidth(arrow_width);
////      glColor3f(0, 0, 0);
////      glBegin(GL_LINES);
// 
//      glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
//      glEnable(GL_LIGHTING);
//
//      for (int dd = 1; dd <=2; dd++) {
//         it = bcc->getCellIterator(dd);
//         while (it.isValid()) {
//            index_type cellid = *it.loc;
//            drawArrow<Dim,FType>(bcc, cellid);
//            it++;
//         }
//      }
//      glDisable(GL_LIGHTING);
////      glEnd();
//   }
//	glEndList();
//
//};


float interpolate2d(float* in, float x, float y, int bX) {
	float xb = (float) ((int) x);
	float yb = (float) ((int) y);

	float v0 = in[ (int) (xb) + ((int) (yb)) * bX];
	float v1 = in[ (int) (xb+1) + ((int) (yb)) * bX];
	float v2 = in[ (int) (xb) + ((int) (yb+1)) * bX];
	float v3 = in[ (int) (xb+1) + ((int) (yb+1)) * bX];

	v0 = v0 * (1.0f - (x - xb)) + v1 * ((x - xb)) ;
	v2 = v2 * (1.0f - (x - xb)) + v3 * ((x - xb)) ;

	v0 = v0 * (1.0f - (y - yb)) + v2 * ((y - yb)) ;

	return v0;
}

float smoothed2d(float* in, float x, float y, int bX) {
	float v0 = 4 * interpolate2d(in, x, y, bX);
	float v1 = interpolate2d(in, x -1.0f, y, bX);
	float v2 = interpolate2d(in, x+1.0f, y, bX);
	float v3 = interpolate2d(in, x, y-1.0f, bX);
	float v4 = interpolate2d(in, x, y+1.0f, bX);

	return (v0 + v1 +v2 +v3 +v4 )/8.0f;
}



double* read_bin_vectors(char* filename){
    FILE* fp = fopen(filename, "rb");
    
    fflush(stdout);
    int nv = 0;
    fread(&nv,sizeof(int),1,fp);
    double* pData = new double[3*nv];
    
    // now read the data
    fread(pData,3*sizeof(double),nv,fp);
    
	double* out2d = new double[2*nv];
	for (int i =0; i < nv; i++) {
		out2d[2*i] = pData[3*i];
		out2d[2*i+1] = pData[3*i+1];
	}
    fclose(fp);
	delete[] pData;
	return out2d;
}


int main(int argc, char** argv) {


	mscTopoArray<float> flarray;
	flarray.append(1.2f);
	flarray.append(3.2f);
	flarray.append(3.4f);
	printf("size of flarray:%d\nfirst element: %f\n", flarray.size(), flarray[0]); 





		for (int i =0; i < argc; i++) {
		printf("%d=%s\n", i, argv[i]);
	}

	char* filename = argv[1];
	int X, Y, Z;
	//float deg;
	sscanf(argv[2], "%d", &oX);
	sscanf(argv[3], "%d", &oY);
	sscanf(argv[4], "%d", &oZ);
	sscanf(argv[5], "%d", &X);
	sscanf(argv[6], "%d", &Y);
	sscanf(argv[7], "%d", &Z);
	sscanf(argv[8], "%f", &deg);

	GDEGREE = -deg;

	deg = deg / 180.0f * 3.14159f;
	
	// input resample!
	FILE* fin = fopen(filename, "rb");
	fbbb = new float[oX*oY*oZ];
	fread(fbbb, sizeof(float), oX*oY*oZ, fin);
	fclose(fin);
	//centroid
	float centroid[2];
	centroid[0] = (float) oX/2.0f;
	centroid[1] = (float) oY/2.0f;
	float centroid2[2];
	centroid2[0] = (float) X/2.0f;
	centroid2[1] = (float) Y/2.0f;
	
	FILE* fout = fopen("tmp_rotated.raw", "wb");
	for (int j = 0; j < Y; j++) {
		for (int i = 0; i < X; i++) {
			float tx = (float) i - centroid2[0]; 
			float ty = (float) j - centroid2[1]; 

			float xnew = tx * cos(deg) - ty * sin(deg);
			float ynew = ty * cos(deg) + tx * sin(deg);

			xnew += centroid[0];
			ynew += centroid[1];

			float res;
			res = smoothed2d(fbbb, xnew, ynew, oX);
			fwrite(&res, sizeof(float), 1, fout);
		}
	}

	fclose(fout);

	filename = new char[1024];
	sprintf(filename, "tmp_rotated.raw");


	OLDX = X;
	OLDY = Y;

	XMIN = 0;
	XMAX = X*2-1;
	YMIN = 0; 
	YMAX = Y*2-1;
	ZMIN = 0; 
	ZMAX = 1;

	int userand;
	sscanf(argv[9], "%d", &userand);

		int seed;	
		if (argc > 10) {

		sscanf(argv[10], "%d", &seed);
		srand(seed);
	}


	mscSize stest(X, Y, Z);
	printf("stest=%d\n", (int) stest.count());



	//FILE* fouttest = fopen("test.raw", "wb");
	//for (int i = 0; i < 3*3*3; i++) {
	//	float val = (float) test_func[i];
	//	fwrite(&val, sizeof(float), 1, fouttest);
	//}
	//fclose(fouttest);

	cg_msc = new mscContinuousGradient2D();
	cg_msc->load_data(filename, X, Y);
	 

   // declare array factory
	mscArrayFactory a(REGULAR_ARRAY);

	// load data
   mscRegularRawDataHandler<float>* test_data;
   test_data = new mscRegularRawDataHandler<float>();
   test_data->load_data(filename, X*Y*Z, &a);





   //test_data->testfart();
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

      // TESTING
   // CENTRAL DIFFERENCES
   double* vectors = new double[X*Y*2];
   for (int i =0; i < X; i++) {
	   for (int j =0; j < Y; j++) {
		   int tmpi0 = (i == 0? 0: i - 1) + j * X;
		   int tmpi1 = (i == X-1? X-1: i + 1) + j * X;
		   vectors[(i + j*X)*2] = test_data->value(tmpi1) - 
			   test_data->value(tmpi0);
		   int tmpj0 = i + (j == 0? 0: j - 1) * X;
		   int tmpj1 = i + (j == Y-1? Y-1: j + 1) * X;
		   vectors[(i + j*X)*2+1] = test_data->value(tmpj1) - 
			   test_data->value(tmpj0);
		   
	   }
   }

   //double* vectors = new double[X*Y*2];
   //for (int i =0; i < X; i++) {
	  // for (int j =0; j < Y; j++) {
		 //  int tmpi0 = (i - 2 < 0? 0: i - 1) + j * X;
		 //  int tmpi1 = (i - 1 < 0? 0: i - 1) + j * X;
		 //  int tmpi2 = (i + 1 >= X-1? X-1: i + 1) + j * X;
		 //  int tmpi3 = (i + 2 >= X-1? X-1: i + 2) + j * X;
		 //  vectors[(i + j*X)*2] = (1.0/12.0)* (
			//   test_data->value(tmpi0) - 
			//   8 * test_data->value(tmpi1) +
			//   8 * test_data->value(tmpi2) -
			//   test_data->value(tmpi3));

		 //  int tmpj0 = i+(j - 2 < 0? 0: j - 1) * X;
		 //  int tmpj1 = i+(j - 1 < 0? 0: j - 1) * X;
		 //  int tmpj2 = i+(j + 1 >= Y-1? Y-1: j + 1) * X;
		 //  int tmpj3 = i+(j + 2 >= Y-1? Y-1: j + 2) * X;
		 //  vectors[(i + j*X)*2+1] = (1.0/12.0)* (
			//   test_data->value(tmpj0) - 
			//   8 * test_data->value(tmpj1) +
			//   8 * test_data->value(tmpj2) -
			//   test_data->value(tmpj3));

		 //  //int tmpj0 = i + (j == 0? 0: j - 1) * X;
		 //  //int tmpj1 = i + (j == Y-1? Y-1: j + 1) * X;
		 //  //vectors[(i + j*X)*2+1] = test_data->value(tmpj1) - 
			//  // test_data->value(tmpj0);
		 //  
	  // }
   //}

   if (argc > 10) {

	   char ffnname[1024];
	  
	   sprintf(ffnname, "ex%d_d.vbin", seed);
	  
	   printf("opening %s...\n", ffnname);
	   vectors = read_bin_vectors(ffnname);

   }

   printf("starting lic...\n");

   double* lic_vect = new double[LIC_X*LIC_Y*2];
   for (int i =0; i < LIC_X; i++) {
	   for (int j = 0; j < LIC_Y; j++) {
		   int ci = (i*OLDX) / LIC_X;
		   int cj = (j*OLDY) / LIC_Y;
		   lic_vect[(i + j*LIC_X)*2] = vectors[(ci + cj*OLDX)*2];
		   lic_vect[(i + j*LIC_X)*2+1] = vectors[(ci + cj*OLDX)*2+1];
	   }
   }


   GLOBAL_LIC_IMAGE = 	setup_LIC(lic_vect);
   printf("done with lic\n");

   v2d_msc = new mscVectorToScalarGridFunc2D(vectors, X, Y, bmsh);

   v2d_msc->init();

   mscVectorNegatorGridFunc2d* v2d_mscN = new mscVectorNegatorGridFunc2d(vectors, X, Y, bmsh);
   v2d_mscN->init();


   //printf("Gothere\n");
   //cellIterator it;
   //iteratorOperator& fit = bmsh->cofacets(1, it);
   //fit.begin(it);

   //while (fit.valid(it)) {
	  // printf("neighbor=%llu\n", fit.value(it));
	  // fit.advance(it);
   //}



   // create mesh function
   mscRegularGrid3DMeshFunction<float>* bmf = new mscRegularGrid3DMeshFunction<float>(v2d_msc, bmsh, &a);
   bmf->initialize();



   mscRegularGrid3DGradientField* bgf = 
	   new mscRegularGrid3DGradientField(bmsh, &a);

   // test the gradient builder!
   mscConvergentVFBuilder* cvfb =
	   new mscConvergentVFBuilder(v2d_msc, bmf, bmsh, bgf, &a);  
      G_msccb2 = cvfb;
	G_msccb = cvfb;
	cvfb->computeGradient();

		   mscComplementMeshHandler* cmh = 
		   new mscComplementMeshHandler(bmsh);

	   mscNegatingMeshFunction<float>* nmf = 
		   new mscNegatingMeshFunction<float>(bmf);
	   mscModifiedBoundaryMeshHandler* mbmh = 
		   new mscModifiedBoundaryMeshHandler(cmh, bgf);
	   mscRegularGrid3DGradientField* bgf2 = 
		   new mscRegularGrid3DGradientField(bmsh, &a);
	   mscConvergentVFBuilder* cvfb2 =
	   new mscConvergentVFBuilder(v2d_mscN, nmf, mbmh, bgf2, &a);  

		 cvfb2->computeGradient();

		  mscRegularGrid3DGradientField* tempgf = bgf;
		  
		  cellIterator it;
		  iteratorOperator& all_cells = bmsh->all_cells_iterator(it);
		  for (all_cells.begin(it); all_cells.valid(it); all_cells.advance(it)) {
			  CELL_INDEX_TYPE cid = all_cells.value(it);
			  bgf2->set_dim_asc_man(cid, bgf->get_dim_asc_man(cid));
		  }
		  bgf = bgf2;

		  bgf->resetMeshHandler(bmsh);


	/*	 char ufilename1[1024];
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
		 fclose(ffout);*/

		 G_msccb = cvfb;
		 G_msccb2 = cvfb2;






 ////  mscSimpleGradientBuilder<float>* mscb = 
	////   new mscSimpleGradientBuilder<float>(bmf, bmsh, bgf, &a);

 ////  mscSimpleRandomGradientBuilder<float>* mscrb = 
	////   new mscSimpleRandomGradientBuilder<float>(bmf, bmsh, bgf, &a);



 ////  //mscTwoWay3DGradientBuilder<float>* msctwb =
	////  // new mscTwoWay3DGradientBuilder<float>(bmf, bmsh, bgf, &a);


 ////  mscConvergentGradientBuilder<float>* msccb = 
	////   new mscConvergentGradientBuilder<float>(bmf, bmsh, bgf, &a);
 ////  G_msccb2 = msccb;
	////G_msccb = msccb;

 ////  if (userand == 0) {
	////	printf("using randomized gradient\n");
	////   mscrb->computeGradient();
 ////  } else if (userand == 1) {
	////   printf("using greedy gradient\n");
	////	mscb->computeGradient();
 ////  } else if (userand == 2) {
	////   printf("using convergent gradient - 1 pass\n");
	////   msccb->computeGradient();
 ////  } else if (userand == 3)  {
	////   printf("using convergent gradient - 2 pass\n");
	////   msccb->computeGradient();


	////   mscComplementMeshHandler* cmh = 
	////	   new mscComplementMeshHandler(bmsh);

	////   mscNegatingMeshFunction<float>* nmf = 
	////	   new mscNegatingMeshFunction<float>(bmf);
	////   mscModifiedBoundaryMeshHandler* mbmh = 
	////	   new mscModifiedBoundaryMeshHandler(cmh, bgf);
	////   mscRegularGrid3DGradientField* bgf2 = 
	////	   new mscRegularGrid3DGradientField(bmsh, &a);
	////     mscConvergentGradientBuilder<float>* msccb2 = 
	////   new mscConvergentGradientBuilder<float>(nmf, mbmh, bgf2, &a);
	////	 msccb2->computeGradient();

	////	  mscRegularGrid3DGradientField* tempgf = bgf;
	////	  
	////	  cellIterator it;
	////	  iteratorOperator& all_cells = bmsh->all_cells_iterator(it);
	////	  for (all_cells.begin(it); all_cells.valid(it); all_cells.advance(it)) {
	////		  CELL_INDEX_TYPE cid = all_cells.value(it);
	////		  bgf2->set_dim_asc_man(cid, bgf->get_dim_asc_man(cid));
	////	  }
	////	  bgf = bgf2;

	////	  bgf->resetMeshHandler(bmsh);


	////	 char ufilename1[1024];
	////	 sprintf(ufilename1, "%s.asc", filename); 
	////	 test_data->dump_vals(ufilename1, X, Y, Z, msccb->getmaxvals());
	////	 sprintf(ufilename1, "%s.asc.prob", filename);

	////	 char ufilename2[1024];

	////	 sprintf(ufilename2, "%s.dsc", filename);
	////	 test_data->dump_vals(ufilename2, X, Y, Z, msccb2->getmaxvals());
 ////  		 sprintf(ufilename2, "%s.dsc.prob", filename);

	////	 char ufilename3[1024];

	////	 sprintf(ufilename3, "%s.both", filename);
	////	 FILE* ff1 = fopen(ufilename1, "rb");
	////	 FILE* ff2 = fopen(ufilename2, "rb");
	////	 FILE* ffout = fopen(ufilename3, "wb");
	////	 while (! feof(ff1)) {
	////		 float v1, v2;
	////		 fread(&v1, sizeof(float), 1, ff1);
	////		 fread(&v2, sizeof(float), 1, ff2);
	////		 float res = (0.5-0.5f*v1) + 0.5f - (0.5-0.5f*v2);
	////		 fwrite(&res, sizeof(float), 1, ffout);
	////	 }
	////	 fclose(ff1);
	////	 fclose(ff2);
	////	 fclose(ffout);

	////	 G_msccb = msccb;
	////	 G_msccb2 = msccb2;
 ////  } else if (userand == 4) {
	////   	   printf("using convergent2 gradient - 2 pass\n");


	////   mscComplementMeshHandler* cmh = 
	////	   new mscComplementMeshHandler(bmsh);

	////   mscNegatingMeshFunction<float>* nmf = 
	////	   new mscNegatingMeshFunction<float>(bmf);
	////   mscRegularGrid3DGradientField* bgf2 = 
	////	   new mscRegularGrid3DGradientField(bmsh, &a);	   

	////     mscConvergentGradientBuilder<float>* msccb2 = 
	////   new mscConvergentGradientBuilder<float>(nmf, cmh, bgf2, &a);
	////	 msccb2->computeGradient();	   
	////   
	////   mscModifiedBoundaryMeshHandler* mbmh = 
	////	   new mscModifiedBoundaryMeshHandler(bmsh, bgf2);

	////     mscConvergentGradientBuilder<float>* msccb3 = 
	////   new mscConvergentGradientBuilder<float>(bmf, mbmh, bgf, &a);


	////   msccb3->computeGradient();
	////	 G_msccb = msccb3;
	////	 G_msccb2 = msccb2;
	////	  //mscRegularGrid3DGradientField* tempgf = bgf;
	////	  //
	////	  //cellIterator it;
	////	  //iteratorOperator& all_cells = bmsh->all_cells_iterator(it);
	////	  //for (all_cells.begin(it); all_cells.valid(it); all_cells.advance(it)) {
	////		 // CELL_INDEX_TYPE cid = all_cells.value(it);
	////		 // bgf2->set_dim_asc_man(cid, bgf->get_dim_asc_man(cid));
	////	  //}
	////	  //bgf = bgf2;

	////	  //bgf->resetMeshHandler(bmsh);



 ////  } else {
 ////  	   	   printf("using convergent3 gradient - 2 pass\n");


	////   mscComplementMeshHandler* cmh = 
	////	   new mscComplementMeshHandler(bmsh);
	////   //mscDumbStoringMinFunction<float>* dsminf = 
	////	  // new mscDumbStoringMinFunction<float>(test_data, bmsh, &a);
	////   //dsminf->initialize();
	////   mscNegatingMeshFunction<float>* nmf = 
	////	   new mscNegatingMeshFunction<float>(bmf);
	////   mscRegularGrid3DGradientField* bgf2 = 
	////	   new mscRegularGrid3DGradientField(bmsh, &a);	   

	////     mscConvergentGradientBuilder<float>* msccb2 = 
	////   new mscConvergentGradientBuilder<float>(nmf, cmh, bgf2, &a);
	////	 msccb2->computeGradient();	   
	////   
	////   mscModifiedBoundaryMeshHandler* mbmh = 
	////	   new mscModifiedBoundaryMeshHandler(bmsh, bgf2);

	////     mscConvergentGradientBuilder<float>* msccb3 = 
	////   new mscConvergentGradientBuilder<float>(bmf, mbmh, bgf, &a);


	////   msccb3->computeGradient();
	////	 G_msccb = msccb3;
	////	 G_msccb2 = msccb2;
	////	  //mscRegularGrid3DGradientField* tempgf = bgf;
	////	  //
	////	  //cellIterator it;
	////	  //iteratorOperator& all_cells = bmsh->all_cells_iterator(it);
	////	  //for (all_cells.begin(it); all_cells.valid(it); all_cells.advance(it)) {
	////		 // CELL_INDEX_TYPE cid = all_cells.value(it);
	////		 // bgf2->set_dim_asc_man(cid, bgf->get_dim_asc_man(cid));
	////	  //}
	////	  //bgf = bgf2;

	////	  //bgf->resetMeshHandler(bmsh);



 ////  }
   G_mscg = bgf;
   G_mscmh = bmsh;
   G_mscf = bmf;



#ifdef WIN32        
		  DWORD globalEnd = GetTickCount();
		  printf(" --Computed discrete gradient in %.3f seconds\n", (globalEnd - globalStart)/1000.0);
#endif




   MSC = new BasicMSC<float>(bgf, bmsh, bmf);
   MSC->ComputeFromGrad();
   MSC->ComputeHeirarchy();
   MSC->setPercentPersistence(0);

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
bool clearcolor = 1;
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

void myidle() {
	cimg_library::cimg::sleep(10);

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

	glRotatef(-GDEGREE, 0, 0, 1.0f);

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
		sprintf(name, "im%d_%d.bmp", gcounter++, (int) GDEGREE);	
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
		CImg< unsigned char > img(mpixels, width, height,1,  3, false);

		img.save_bmp(name);
		delete(rpixels);
		delete(mpixels);
		printf("output %s \n", name);
		return;

}


void outputASDF() {
	reset_view();
	MSC->setPercentPersistence(2);
	  redrawstuff = true;
	display();
	MSC->setPercentPersistence(2);
	  redrawstuff = true;
	display();

	outputimage();
	exit(1);
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
	case 'O':
		outputASDF();
		break;
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
		MSC->setPercentPersistence(2);
	  redrawstuff = true;
	  break;
	case '3':
		MSC->setPercentPersistence(15);

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
//	  glivingcounter = bsn->getDestrCount(0.1f); 
	  redrawstuff = true;
	  break;
	case '7':

	MSC->setPercentPersistence(0.015);

	  redrawstuff = true;


//	  glivingcounter = bsn->getDestrCount(0.2f);
	  redrawstuff = true;
	  break;
	case '8':
//	  glivingcounter = bsn->getDestrCount(0.5f);
	MSC->setPercentPersistence(0.15);

	  redrawstuff = true;
	  redrawstuff = true;
	  break;
	case '9':
	//  glivingcounter = bsn->getDestrCount(1.0f);
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
     {
 /*      for (int i=0; i<10; i++) {
	 if (cutoff > 0) {
	   if (usc->getCritical(total_order[cutoff])) {
	     usc->setAssigned(total_order[cutoff], false);
	     cutoff--;
	   } else {
	     usc->setAssigned(total_order[cutoff], false);
	     cutoff--;
	     usc->setAssigned(total_order[cutoff], false);
	     cutoff--;
	   }
	 }
       }
  */   }
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
	case 'N':
	  {
	 //   for (int i=0; i<1; i++) {
	 //     if (cutoff > 0) {
		//if (usc->getCritical(total_order[cutoff])) {
		//  usc->setAssigned(total_order[cutoff], false);
		//  cutoff--;
		//} else {
		//  usc->setAssigned(total_order[cutoff], false);
		//  cutoff--;
		//  usc->setAssigned(total_order[cutoff], false);
		//  cutoff--;
		//}
	 //     }
	 //   }
	  }
	  redrawstuff = true;
	  break;


	case 's':
	  makeSequence();
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
	case '+':
		MANALPHA += 0.1f;
		redrawstuff = true;
		break;
	case '-':
		MANALPHA -= 0.1f;
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

  if (buttonstate == 1) {

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

  if (buttonstate == 0) {
	float ldim = (width > height ? width : height);

	if (xold == mx && yold == my) return;

    vstart.vals[0] = xold - .5 * width ;
    vstart.vals[1] = .5 * height - yold;

    if ((.5 * ldim) * (.5* ldim) - vstart.vals[0]*vstart.vals[0] - vstart.vals[1]*vstart.vals[1] < 0) return;
    vstart.vals[2] = sqrt((.5 * ldim) * (.5* ldim) - vstart.vals[0]*vstart.vals[0] - vstart.vals[1]*vstart.vals[1]);
    if (Normalize3(&vstart)== 0) return;
    xold = mx;
    yold = my;

    // calculate the vectors;
    Vector4 vend;


    vend.vals[0] = mx - .5 * width ;
    vend.vals[1] = .5 * height - my;

    if ((.5 * ldim) * (.5* ldim) - vend.vals[0]*vend.vals[0] - vend.vals[1]*vend.vals[1] <= 0) return;
    vend.vals[2] = sqrt((.5 * ldim) * (.5* ldim) - vend.vals[0]*vend.vals[0] - vend.vals[1]*vend.vals[1]);
    if (Normalize3(&vend)== 0) return;

    

    float alpha = InteriorAngle(vstart, vend);
	if (alpha < 0.01) return;

    Vector4 cp = Cross(vstart, vend);
    if (Normalize3(&cp) == 0) return;
   
    //update the crystal ball matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
 
	rotationtotal = times( rotationtotal,rotation(cp, alpha*PI/-180.0));
	
	setRotQuat(rotmatrix, rotationtotal);
	glMultMatrixf(rotmatrix.vals);

  }
 

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






					   ////////////////////////////////////////////////////////////////////////////
					   ///		     Line Integral Convolution for Flow Visualization			///
					   ///							 Initial  Version							///
					   ///							   May 15, 1999								///
					   ///									by									///
					   ///                             Zhanping Liu								///
					   ///                      (zhanping@erc.msstate.edu)                      ///
					   ///								while with								///
					   ///						   Graphics  Laboratory							///
					   ///						    Peking  University							///
					   ///							   P. R.  China							    ///
					   ///----------------------------------------------------------------------///
					   ///							 Later  Condensed							///
					   ///                             May 4,  2002								///
					   ///          VAIL (Visualization Analysis & Imaging Laboratory)          ///
					   ///                  ERC  (Engineering Research Center)                  ///
					   ///                     Mississippi State University                     ///
					   ////////////////////////////////////////////////////////////////////////////

					
					   ////////////////////////////////////////////////////////////////////////////
					   ///This code was developed  based on  the original algorithm  proposed by///
					   ///Brian Cabral  and  Leith (Casey) Leedom  in the paper  "Imaging Vector///
					   ///Fields Using Line Integral Convolution", published  in  Proceedings of///
					   ///ACM  SigGraph 93,  Aug 2-6,  Anaheim,  California,  pp. 263-270, 1993.///
					   ///Permission to use, copy, modify, distribute and sell this code for any///
					   ///purpose is hereby granted without fee,  provided that the above notice///
					   ///appears in all copies  and  that both that notice and  this permission///
					   ///appear in supporting documentation. The  developer  of this code makes///
					   ///no  representations  about  the  suitability  of  this  code  for  any///
					   ///purpose.  It is provided "as is"  without express or implied warranty.///
					   //////////////////////////////////////////////////////////////////////////// 


#include "stdafx.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>



#define	 DISCRETE_FILTER_SIZE	1024
#define  LOWPASS_FILTR_LENGTH	10.00000f
#define	 LINE_SQUARE_CLIP_MAX	100000.0f
#define	 VECTOR_COMPONENT_MIN   0.050000f 



void	 NormalizVectrs(int  n_xres,  int     n_yres,  float*   pVectr);
void     GenBoxFiltrLUT(int  LUTsiz,  float*  p_LUT0,  float*   p_LUT1);
void     MakeWhiteNoise(int  n_xres,  int     n_yres,  unsigned char*  pNoise);
void	 FlowImagingLIC(int  n_xres,  int     n_yres,  float*   pVectr,   unsigned char*  pNoise,  
						unsigned char*  pImage,  float*  p_LUT0,  float*  p_LUT1,  float  krnlen);
//void 	 WriteImage2PPM(int  n_xres,  int     n_yres,  unsigned char*  pImage,     char*  f_name);


unsigned char*	setup_LIC(double* VVV)
{
		int				n_xres = LIC_X;
		int				n_yres = LIC_Y;
		float*			pVectr = (float*         ) malloc( sizeof(float        ) * n_xres * n_yres * 2 );
		float*			p_LUT0 = (float*		 ) malloc( sizeof(float        ) * DISCRETE_FILTER_SIZE);
		float*			p_LUT1 = (float*		 ) malloc( sizeof(float        ) * DISCRETE_FILTER_SIZE);
		unsigned char*	pNoise = (unsigned char* ) malloc( sizeof(unsigned char) * n_xres * n_yres     );
		unsigned char*	pImage = (unsigned char* ) malloc( sizeof(unsigned char) * n_xres * n_yres     );

		for (int i = 0; i < n_xres*n_yres*2; i++) pVectr[i] = (float) VVV[i];
		NormalizVectrs(n_xres, n_yres, pVectr);
		MakeWhiteNoise(n_xres, n_yres, pNoise);
		GenBoxFiltrLUT(DISCRETE_FILTER_SIZE, p_LUT0, p_LUT1);
		FlowImagingLIC(n_xres, n_yres, pVectr, pNoise, pImage, p_LUT0, p_LUT1, LOWPASS_FILTR_LENGTH);
		////WriteImage2PPM(n_xres, n_yres, pImage, "LIC.ppm");

		//free(pVectr);	pVectr = NULL;
		//free(p_LUT0);	p_LUT0 = NULL;
		//free(p_LUT1);	p_LUT1 = NULL;
		//free(pNoise);	pNoise = NULL;
		//free(pImage);	pImage = NULL;
		return pImage;
}





///		normalize the vector field     ///
void    NormalizVectrs(int  n_xres,  int  n_yres,  float*  pVectr)
{	
		for(int	 j = 0;  j < n_yres;  j ++)
    	for(int	 i = 0;  i < n_xres;  i ++)
    	{	
			int		index = (j * n_xres + i) << 1;
        	float	vcMag = float(  sqrt( double(pVectr[index] * pVectr[index] + pVectr[index + 1] * pVectr[index + 1]) )  );

			float	scale = (vcMag == 0.0f) ? 0.0f : 1.0f / vcMag;
			pVectr[index    ] *= scale;
            pVectr[index + 1] *= scale;
    	}
}


///		make white noise as the LIC input texture     ///
void	MakeWhiteNoise(int  n_xres,  int  n_yres,  unsigned char*  pNoise)
{		
		for(int  j = 0;   j < n_yres;  j ++)
    	for(int  i = 0;   i < n_xres;  i ++)
    	{	
			int  r = rand();
		    r = (  (r & 0xff) + ( (r & 0xff00) >> 8 )  ) & 0xff;
        	pNoise[j * n_xres + i] = (unsigned char) r;
    	}
}


///		generate box filter LUTs     ///
void    GenBoxFiltrLUT(int  LUTsiz,  float*  p_LUT0,  float*  p_LUT1)
{  		
   		for(int  i = 0;  i < LUTsiz;  i ++)  p_LUT0[i] = p_LUT1[i] = i;
}


///////		write the LIC image to a PPM file     ///
////void	WriteImage2PPM(int  n_xres,  int  n_yres,  unsigned char*  pImage,  char*  f_name)
////{	
////		FILE*	o_file;
////		if(   ( o_file = fopen(f_name, "w") )  ==  NULL   )  
////		{  
////			printf("Can't open output file\n");  
////			return;  
////		}
////
////  		fprintf(o_file, "P6\n%d %d\n255\n", n_xres, n_yres);
////
////  		for(int  j = 0;  j < n_yres;  j ++)
////   		for(int  i = 0;  i < n_xres;  i ++)
////		{
////			unsigned  char	unchar = pImage[j * n_xres + i];
////      		fprintf(o_file, "%c%c%c", unchar, unchar, unchar);
////		}
////
////  		fclose (o_file);	o_file = NULL;
////}


///		flow imaging (visualization) through Line Integral Convolution     ///
void	FlowImagingLIC(int     n_xres,  int     n_yres,  float*  pVectr,  unsigned char*  pNoise,  unsigned char*  pImage,  
					   float*  p_LUT0,  float*  p_LUT1,  float   krnlen)
{	
		int		vec_id;						///ID in the VECtor buffer (for the input flow field)
		int		advDir;						///ADVection DIRection (0: positive;  1: negative)
		int		advcts;						///number of ADVeCTion stepS per direction (a step counter)
		int		ADVCTS = int(krnlen * 3);	///MAXIMUM number of advection steps per direction to break dead loops	
		
		float	vctr_x;						///x-component  of the VeCToR at the forefront point
		float	vctr_y;						///y-component  of the VeCToR at the forefront point
		float	clp0_x;						///x-coordinate of CLiP point 0 (current)
		float	clp0_y;						///y-coordinate of CLiP point 0	(current)
		float	clp1_x;						///x-coordinate of CLiP point 1 (next   )
		float	clp1_y;						///y-coordinate of CLiP point 1 (next   )
		float	samp_x;						///x-coordinate of the SAMPle in the current pixel
		float	samp_y;						///y-coordinate of the SAMPle in the current pixel
		float	tmpLen;						///TeMPorary LENgth of a trial clipped-segment
		float	segLen;						///SEGment   LENgth
		float	curLen;						///CURrent   LENgth of the streamline
		float	prvLen;						///PReVious  LENgth of the streamline		
		float	W_ACUM;						///ACcuMulated Weight from the seed to the current streamline forefront
		float	texVal;						///TEXture VALue
		float	smpWgt;						///WeiGhT of the current SaMPle
		float	t_acum[2];					///two ACcUMulated composite Textures for the two directions, perspectively
		float	w_acum[2];					///two ACcUMulated Weighting values   for the two directions, perspectively
		float*	wgtLUT = NULL;				///WeiGhT Look Up Table pointing to the target filter LUT
  		float	len2ID = (DISCRETE_FILTER_SIZE - 1) / krnlen;	///map a curve LENgth TO an ID in the LUT

		///for each pixel in the 2D output LIC image///
		for(int  j = 0;	 j < n_yres;  j ++)
		for(int  i = 0;	 i < n_xres;  i ++)
		{	
			///init the composite texture accumulators and the weight accumulators///
			t_acum[0] = t_acum[1] = w_acum[0] = w_acum[1] = 0.0f;
		
			///for either advection direction///
        	for(advDir = 0;  advDir < 2;  advDir ++)
        	{	
				///init the step counter, curve-length measurer, and streamline seed///
				advcts = 0;
				curLen = 0.0f;
            	clp0_x = i + 0.5f;
				clp0_y = j + 0.5f;

				///access the target filter LUT///
				wgtLUT = (advDir == 0) ? p_LUT0 : p_LUT1;

				///until the streamline is advected long enough or a tightly  spiralling center / focus is encountered///
            	while( curLen < krnlen && advcts < ADVCTS ) 
	         	{
					///access the vector at the sample///
					vec_id = ( int(clp0_y) * n_xres + int(clp0_x) ) << 1;
					vctr_x = pVectr[vec_id    ];
					vctr_y = pVectr[vec_id + 1];

					///in case of a critical point///
					if( vctr_x == 0.0f && vctr_y == 0.0f )
					{	
						t_acum[advDir] = (advcts == 0) ? 0.0f : t_acum[advDir];		   ///this line is indeed unnecessary
						w_acum[advDir] = (advcts == 0) ? 1.0f : w_acum[advDir];
						break;
					}
					
					///negate the vector for the backward-advection case///
					vctr_x = (advDir == 0) ? vctr_x : -vctr_x;
					vctr_y = (advDir == 0) ? vctr_y : -vctr_y;

					///clip the segment against the pixel boundaries --- find the shorter from the two clipped segments///
					///replace  all  if-statements  whenever  possible  as  they  might  affect the computational speed///
					segLen = LINE_SQUARE_CLIP_MAX;
					segLen = (vctr_x < -VECTOR_COMPONENT_MIN) ? ( int(     clp0_x         ) - clp0_x ) / vctr_x : segLen;
					segLen = (vctr_x >  VECTOR_COMPONENT_MIN) ? ( int( int(clp0_x) + 1.5f ) - clp0_x ) / vctr_x : segLen;
					segLen = (vctr_y < -VECTOR_COMPONENT_MIN) ?	
							 (      (    (  tmpLen = ( int(     clp0_y)          - clp0_y ) / vctr_y  )  <  segLen    ) ? tmpLen : segLen      ) 
							: segLen;
					segLen = (vctr_y >  VECTOR_COMPONENT_MIN) ?
							 (      (    (  tmpLen = ( int( int(clp0_y) + 1.5f ) - clp0_y ) / vctr_y  )  <  segLen    ) ? tmpLen : segLen      ) 
							: segLen;
					
					///update the curve-length measurers///
					prvLen = curLen;
					curLen+= segLen;
					segLen+= 0.0004f;
			   
					///check if the filter has reached either end///
					segLen = (curLen > krnlen) ? ( (curLen = krnlen) - prvLen ) : segLen;

					///obtain the next clip point///
					clp1_x = clp0_x + vctr_x * segLen;
					clp1_y = clp0_y + vctr_y * segLen;

					///obtain the middle point of the segment as the texture-contributing sample///
					samp_x = (clp0_x + clp1_x) * 0.5f;
					samp_y = (clp0_y + clp1_y) * 0.5f;

					///obtain the texture value of the sample///
					texVal = pNoise[ int(samp_y) * n_xres + int(samp_x) ];

					///update the accumulated weight and the accumulated composite texture (texture x weight)///
					W_ACUM = wgtLUT[ int(curLen * len2ID) ];
					smpWgt = W_ACUM - w_acum[advDir];			
					w_acum[advDir]  = W_ACUM;								
					t_acum[advDir] += texVal * smpWgt;
				
					///update the step counter and the "current" clip point///
					advcts ++;
					clp0_x = clp1_x;
					clp0_y = clp1_y;

					///check if the streamline has gone beyond the flow field///
					if( clp0_x < 0.0f || clp0_x >= n_xres || clp0_y < 0.0f || clp0_y >= n_yres)  break;
				}	
         	}

			///normalize the accumulated composite texture///
         	texVal = (t_acum[0] + t_acum[1]) / (w_acum[0] + w_acum[1]);

			///clamp the texture value against the displayable intensity range [0, 255]
			texVal = (texVal <   0.0f) ?   0.0f : texVal;
			texVal = (texVal > 255.0f) ? 255.0f : texVal; 
			pImage[j * n_xres + i] = (unsigned char) texVal;
		} 	
}