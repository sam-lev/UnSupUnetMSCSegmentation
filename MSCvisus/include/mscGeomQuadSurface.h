#ifndef MSC_GEOM_QUAD_SURFACE
#define MSC_GEOM_QUAD_SURFACE

#include "mscIndexTypes.h"
#include <set>
#include <vector>
#include <map>

#ifndef WIN32
#include <cstdlib>
#endif

using namespace std;

struct geombitfield {
	unsigned char f1 : 1;
	unsigned char f2 : 1;
	unsigned char f3 : 1;
	unsigned char f4 : 5;
};

struct g_vertex {
	float coords[3];
	float norms[3];
	int id;
	geombitfield flags;
	set<int> adjacent;
	CELL_INDEX_TYPE insertID;

};

struct g_quad {
	geombitfield flags;
	int v1;
	int v2;
	int v3;
	int v4;
	float centroid[3];
};


struct g_edge {
	geombitfield flags;
	int v1;
	int v2;
	float centroid[3];
};

class mscGeomQuadSurface {
protected:


	int make_key(int v1, int v2) {
		int a = (int) min(v1, v2);
		int b = (int) max(v1, v2);
		int res = (int) a;
		res = (res << 16) + (int) b;
		return res;
	}

public:
	vector< g_quad > quads;
	vector< g_vertex > vertices;
	vector< g_edge > edges;

	map< CELL_INDEX_TYPE, int > vertex_map;	
	map< int, int > edge_map;

	mscGeomQuadSurface() {
	}

	// returns the id of the vertex added, or the id of the vertex
	// if it alredy exists
	int add_vertex(float* coords, CELL_INDEX_TYPE gid) {
		if(vertex_map.count(gid) > 0) {
			return vertex_map[gid];
		}
		g_vertex v;
		for (int i = 0; i < 3; i++) v.coords[i] = coords[i];
		v.flags.f4 = 0;
		v.insertID = gid;
		int position = (int) vertices.size();
		v.id = position;
		vertices.push_back(v);
		vertex_map[gid] = position;
		return position;
	}


	int add_edge(int v1, int v2, CELL_INDEX_TYPE keyval) {
		g_edge e;
		e.v1 = v1; e.v2 = v2; e.flags.f4 = 1;
		int pos = (int) edges.size();
		edges.push_back(e);
		edge_map[keyval] = pos;
		return pos;
	}

	void add_quad_edge(int v1, int v2) {
		int keyval = make_key(v1, v2);
		int edgepos;
		if (edge_map.count(keyval) > 0) {
			edgepos = edge_map[keyval];
			g_edge& e = edges[edgepos];
			e.flags.f4++;
		} else {
			edgepos = add_edge(v1, v2, keyval);
		}
		vertices[v1].adjacent.insert(edgepos);
		vertices[v2].adjacent.insert(edgepos);
	}

	int add_quad(int v1, int v2, int v3, int v4) {
		g_quad q;
		q.flags.f4 = 0;
		q.v1 = v1; q.v2 = v2; q.v3 = v3; q.v4 = v4;
		int pos = (int) quads.size();
		quads.push_back(q);		
		// now connect them
		// and increment the number of quads touching that vertex
		vertices[v1].flags.f4++;
		add_quad_edge(v1, v2);

		vertices[v2].flags.f4++;
		add_quad_edge(v2, v3);

		vertices[v3].flags.f4++;
		add_quad_edge(v3, v4);

		vertices[v4].flags.f4++;
		add_quad_edge(v4, v1);

		return pos;
	}

	void merge_surface( mscGeomQuadSurface* s) {
		for (int i = 0; i < s->vertices.size(); i++) {
			g_vertex& v = s->vertices[i];
			this->add_vertex(v.coords, v.insertID);
		}
		for (int i = 0; i < s->quads.size(); i++) {
			g_quad& q = s->quads[i];
			int v1 = add_vertex(NULL, s->vertices[q.v1].insertID);
			int v2 = add_vertex(NULL, s->vertices[q.v2].insertID);
			int v3 = add_vertex(NULL, s->vertices[q.v3].insertID);
			int v4 = add_vertex(NULL, s->vertices[q.v4].insertID);
			this->add_quad(v1, v2, v3, v4);
		}

	}

	void dumpInObjWN(char* filename) {
		FILE* obj = fopen(filename, "w");
		fprintf(obj, "# my object\n");

		for (int i = 0; i < vertices.size(); i++) {
			g_vertex &v = vertices[i];
			fprintf(obj, "v %.4f %.4f %.4f\n", v.coords[0], v.coords[1], v.coords[2]);
		}
		for (int i = 0; i < vertices.size(); i++) {
			g_vertex &v = vertices[i];
			fprintf(obj, "vn %.4f %.4f %.4f\n", v.norms[0], v.norms[1], v.norms[2]);
		}
		for (int i = 0; i < quads.size(); i++) {
			g_quad &q = quads[i];
			fprintf(obj, "f %d//%d %d//%d %d//%d %d//%d\n", q.v1 +1, q.v1 +1, 
				q.v2 +1 ,q.v2 +1 , q.v3 +1,q.v3 +1, q.v4 +1, q.v4 +1);
			//fprintf(obj, "f %d %d %d\n", q.v3 +1, q.v4 +1, q.v1 +1);
		}
		fclose(obj);
	}

	void smooth(int iterations) {

		for (int i = 0; i < vertices.size(); i++) {
			g_vertex &v = vertices[i];
			int minnum = 11110;
			set<int>::iterator it = v.adjacent.begin();
			while(it != v.adjacent.end()){
				g_edge &e = edges[*it];
				if (e.flags.f4 < minnum) minnum = e.flags.f4;
				it++;
			}
			v.flags.f4 = minnum;
		}
		while(iterations > 0) {
			iterations--;
		// compute centroids
		for (int i = 0; i < edges.size(); i++) {
			g_edge &e = edges[i];
			g_vertex& v1 = vertices[e.v1];
			for(int j=0;j<3;j++) e.centroid[j] = v1.coords[j];
			
			g_vertex& v2 = vertices[e.v2];
			for(int j=0;j<3;j++) e.centroid[j] += v2.coords[j];

			for(int j=0;j<3;j++) e.centroid[j] *= 0.5f;	
		}
		// now move vertices!
		for (int i = 0; i < vertices.size(); i++) {
			g_vertex& v = vertices[i];
			
			//printf(" v=%d %d \n", v.flags.f4, v.adjacent.size());
			//if (v.flags.f4 <= 1) continue;
			//clear the position
			for(int j=0;j<3;j++) v.coords[j] = 0.0f;
			float numadded = 0.0f;
			set<int>::iterator it = v.adjacent.begin();
			while(it != v.adjacent.end()){
				g_edge &e = edges[*it];
				if (v.flags.f4 >= e.flags.f4) {
					for(int j=0;j<3;j++) v.coords[j] += e.centroid[j];
					numadded += 1.0f;
				}
				it++;
			}
			float scale = 1.0f / numadded;
			for(int j=0;j<3;j++) v.coords[j] *= scale;			
		}
		}


	}
};


#endif
