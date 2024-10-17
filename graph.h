#pragma once
#ifndef __GRAPH_H
#define __GRAPH_H
#include "utility.h"

using namespace std;
using namespace std::chrono; 

struct EdgeInfo {
    bool exists;    // false if no edge, true if edge exists
    int count;  // Number of times the edge has been used
	// initialized to zero.
};

class Graph{
public:
	// Graph(std::string dir);
	Graph(string dir, bool is_Graph );
	Graph(const Graph& other); 
	Graph();
	~Graph() {
		// destroy();
	}
	// note that addEdgeRaw is processing raw ids.
	void addEdgeRaw(vid_t u, vid_t v);

	// addEdge is processing real ids. 
	void addEdge(vid_t u, vid_t v);
	// deleteEdge is processing real ids. 
	void deleteEdge(vid_t u, vid_t v);


public:

	void init(unsigned int num_v1, unsigned int num_v2);

	// for bipartite graphs:
	void loadGraph(std::string dir);

	// for general graphs
	void load_general_graph(std::string);

	void print_graph();

	void show_nodes();

	bool same_layer(vid_t u, vid_t v ){
		if(is_upper(u) && is_upper(v)) return true; 
		if(is_lower(u) && is_lower(v)) return true; 
		return false ; 
	}

	bool diff_layer(vid_t u, vid_t v ){
		return !same_layer(u,v);
	}

	bool is_upper(vid_t u){ 
		return (u < this->num_v1);
	}
	bool is_lower(vid_t u){
		return (u>= this->num_v1);
	}

	bool is_active(int vertex){ 
		return degree[vertex]>0;
	}
	int num_nodes(){ 
		if(is_bipartite){
			return this->num_v1+this->num_v2;
		}else{
			return num_vertices;
		}
		// return this->num_v1+this->num_v2;
	}

	// need to use neighborhash or global edge lookup to speed up this process.
	bool has(int i, int j); 

	bool has_computed(int u, int v){
		assert(u!=v);
		int smaller = (u < v) ? u : v;
		int larger = (u < v) ? v : u;
		if(memory_efficient){
			return __edge_vector[smaller].find(larger) != __edge_vector[smaller].end();
			// return amat[larger][smaller] >=0; 
		}else{
			return edge_vector[smaller].find(larger) != edge_vector[smaller].end();
		}
	}

	std::string dir;

	num_t num_v1;
	num_t num_v2;

	long num_edges, num_vertices;

	vector<vector<vid_t>> neighbor; 

	vector<int> degree; 

	vector<int> U, V;

	vector<int> nodes; 
public:
	bool memory_efficient = false;
	int v1_max_degree;
	int v2_max_degree;

	// the switch determing if the incoming graph is bipartite
	bool is_bipartite = false; 

	vector<unordered_map<int, bool>> edge_vector;

	vector<unordered_map<int, EdgeInfo>> __edge_vector;
};

#endif  /* __GRAPH_H */