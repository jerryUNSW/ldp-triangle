#include "graph.h"

using namespace std;

// a copy constructor: this only copies the skeleton of the graph, not the meaning members like neighbors, et.c 
Graph::Graph(const Graph& other) {

	// cout<<"cp constructor called "<<endl;

	this->is_bipartite = other.is_bipartite;
	if(this->is_bipartite){
		this->dir = other.dir;

		this->init(other.num_v1, other.num_v2);

		// this->num_edges = 0;
		this->v1_max_degree = 0;
		this->v2_max_degree = 0;
	}else{

		// cout<<"copying a general graph"<<endl;

		this->dir = other.dir;

		this-> num_vertices = other.num_vertices ; 

		neighbor.resize(this-> num_vertices);
		degree.resize(this-> num_vertices);
		fill_n(degree.begin(), this-> num_vertices, 0);

		// this->num_edges = 0;
		this->v1_max_degree = 0;
		// edge_map.clear();
		edge_vector.clear();
		edge_vector.resize(this-> num_vertices);
	}
}

// add a switch to this to distinguish between Graph and general graph
Graph::Graph(string dir, bool is_Graph )
{
	num_v1 = 0;
	num_v2 = 0;
	num_edges = 0;
	v1_max_degree=0; 
	v2_max_degree=0;
	neighbor.clear();
	degree.clear();
	this->dir = dir;

	// do this for bipartite graphs: 
	if(is_Graph){
		cout<<"processing bipartite graphs"<<endl;
		loadGraph(dir);
	}else{
		cout<<"processing general graphs"<<endl;
		this->is_bipartite =false;
		load_general_graph(dir);
	}
}

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2> &p) const {
        return std::hash<T1>()(p.first) ^ std::hash<T2>()(p.second);
    }
};

// default constructor 
Graph::Graph() {
	dir = "";
	num_v1 = 0;
	num_v2 = 0;
	num_edges = 0;
	v1_max_degree=0; 
	v2_max_degree=0;
	neighbor.clear();
	degree.clear();
}

void Graph::print_graph()
{
	print_dash(50);
	cout<<"print_graph() called\n";
	for(int i=0;i<degree.size();i++){
		vector<vid_t> NB = neighbor[i];
		cout<<i<<": ";
		vector_show<vid_t>(NB); 
	}
}
void Graph::show_nodes()
{
	cout<<"show_nodes() called\n";
	// cout<<"upper nodes: "; 
	for(int i=0;i<num_v1+num_v2;i++){
		if(degree[i]>0){
			cout<<i<<" ";
		}
	}
	cout<<endl;
}
// initialize the Graph with size requirements  
void Graph::init(unsigned int num1, unsigned int num2)
{
	num_v1 = num1;
	num_v2 = num2;
	num_edges = 0;
	neighbor.resize(num_v1+num_v2);

	degree.resize(num_v1+num_v2);

	// prio.resize(num_v1+num_v2);

	fill_n(degree.begin(), num_v1+num_v2, 0);
	// neighborHash.resize(num_v1+num_v2);
}
/*
void Graph::computePriority(){
	std::vector<std::pair<int, int>> vertexDegrees(num_nodes());
	for (int i = 0; i < num_nodes(); i++) vertexDegrees[i] = std::make_pair(i, degree[i]);
	// Sort the vertexDegrees based on degrees and IDs
	std::sort(vertexDegrees.begin(), vertexDegrees.end(), [](const auto& a, const auto& b) {
		// the vertex with higher priority has lower degree.
		if (a.second != b.second) {
			return a.second > b.second; 
		} else {
			return a.first > b.first; 
		}
	});
	for(int i=0; i<num_nodes() ;i++)
	{
		neighbor[i].shrink_to_fit();
		sort(neighbor[i].begin(), neighbor[i].end());
		// cout<<"rank = "<<i<<",  id = "<<vertexDegrees[i].first<<",  deg = "<<vertexDegrees[i].second<<endl;
		prio[vertexDegrees[i].first] = i; 
	}
}
*/
void Graph::loadGraph(string dir)
{
	unsigned int n1, n2;
	unsigned int edges = 0;
	int u, v, r;
	string metaFile = dir + ".meta";
	string edgeFile = dir + ".e";
	FILE * metaGraph = fopen(metaFile.c_str(), "r");
	FILE * edgeGraph = fopen(edgeFile.c_str(), "r");

	// bool include_amat = true ;
	// scan the meta file and read number of nodes 
	if (fscanf(metaGraph, "%d\n%d", &n1, &n2) != 2)
	{
		fprintf(stderr, "Bad file format: n1 n2 incorrect\n");
		exit(1);
	}
	fprintf(stdout, "n1: %d, n2: %d\n", n1, n2);

	n1 = n1 * vertex_ratio; 
	n2 = n2 * vertex_ratio; 

	fprintf(stdout, "sampled: n1: %d, n2: %d\n", n1, n2);

	init(n1, n2);

	// if(include_amat){a_mat = new map<int, int>[n1+n2+1];}

	while ((r = fscanf(edgeGraph, "%d %d", &u, &v)) != EOF)
	{
		//fprintf(stderr, "%d, %d\n", u, v);
		if (r != 2)
		{
			fprintf(stderr, "Bad file format: u v incorrect\n");
			exit(1);
		}
		if(u<n1 && v<n2)
			addEdgeRaw(u,v);
		// neighbor[u].push_back(v+num_v1);
		// neighbor[v+num_v1].push_back(u);

		// if(include_amat){a_mat[u][v+num_v1] = a_mat[v+num_v1][u] = 1;}
	}
	
	cout<<"|E| = "<<num_edges<<endl;

	cout<<"fill_rate = "<<num_edges*1.0/(num_v1 * num_v2)<<endl;

	fclose(metaGraph);
	fclose(edgeGraph);
}


bool Graph::has(int i, int j){
	// if(this->is_bipartite){
	return std::find(neighbor[i].begin(), neighbor[i].end(), j) != neighbor[i].end(); 
	// }
	// else{
	// 	// edge_vector based approach 
	// 	bool new_res = edge_vector[i].find(j) != edge_vector[i].end(); 
	// 	return new_res;
	// 	// it is possible that edge_vector[i][j] 
	// }
}


void Graph::load_general_graph(string dir)
{
	unsigned int n__, m__;

	int u, v, r;

	string edgeFile = dir; 

	FILE * edgeGraph = fopen(edgeFile.c_str(), "r");

	if (fscanf(edgeGraph, "%d\n%d", &n__, &m__) != 2)
	{
		fprintf(stderr, "Bad file format: n m incorrect\n");
		exit(1);
	}

	this-> num_vertices = n__ ; 
	this-> num_edges    = m__ ; 

	// cout<<"this-> num_vertices = "<< this-> num_vertices <<endl;
	// cout<<"this-> num_edges = "<< this-> num_edges <<endl;

	v1_max_degree = 0; // on general graphs, we use this to represent the maximum degree

	neighbor.resize(this-> num_vertices);
	degree.resize(this-> num_vertices);
	fill_n(degree.begin(), this-> num_vertices, 0);


	// edge_vector based approahc 
	// if(memory_efficient){
	// 	__edge_vector.resize(this-> num_vertices);
	// }else{
	// 	edge_vector.resize(this-> num_vertices);
	// }


	// initialize a_mat 
	// a_mat.resize(this-> num_vertices); 

	int edge_counter = 0;
	while ((r = fscanf(edgeGraph, "%d %d", &u, &v)) != EOF)
	{
		if (r != 2)
		{
			fprintf(stderr, "Bad file format: u v incorrect\n");
			exit(1);
		}
		// keep record of neighbors and degrees

		// no duplicate edges in this way
		// if(u>=v){
		// 	continue;
		// }

		neighbor[u].push_back(v);
		neighbor[v].push_back(u);		
		// do we really need this for the original graph G?
		// edge_vector[u][v] = true;
		// edge_vector[v][u] = true;

		// false indicate 0 entry 

		// counting the visited edges
		edge_counter++;
	}

	// computing the maximum degree.
	for(int u=0;u<num_nodes();u++){
		degree[u] = neighbor[u].size();
		v1_max_degree = v1_max_degree > degree[u] ? v1_max_degree : degree[u]; 
	}

	cout<<"|E| = "<<num_edges<<endl;

	cout<<"counted edges = "<<edge_counter <<endl;

	fclose(edgeGraph);
}

bool satisfy_bound(long double upper, long double lower, int u, int x, int v, int w, Graph& g){
	// return (g.degree[u] <= upper) & (g.degree[x] <= upper) & (g.degree[v] <= lower) & (g.degree[w] <= lower); 

	return (g.degree[u] >= upper) & (g.degree[x] >= upper) & (g.degree[v] >= lower) & (g.degree[w] >= lower); 
	// all heavy vertices. 
}

// u,v are raw ids, convert them into real ids, and then update neighbor and degree
void Graph::addEdgeRaw(vid_t u, vid_t v)
{
	neighbor[u].push_back(v+num_v1);

	degree[u]++; 

	neighbor[v+num_v1].push_back(u);

	degree[v+num_v1]++;

	num_edges++;
	v1_max_degree = v1_max_degree > degree[u] ? v1_max_degree : degree[u]; 
	v2_max_degree = v2_max_degree > degree[v+num_v1] ? v2_max_degree : degree[v+num_v1]; 
}

// maybe this function takes too long? 
void Graph::addEdge(vid_t u, vid_t v)
{	
	// it looks like only these two lines are necessary 
	neighbor[u].push_back(v);
	neighbor[v].push_back(u);

	degree[u]++; 
	degree[v]++;
	num_edges++;

	if(is_bipartite){
		v1_max_degree = v1_max_degree > degree[u] ? v1_max_degree : degree[u]; 
		v2_max_degree = v2_max_degree > degree[v] ? v2_max_degree : degree[v]; 
	}else{
		edge_vector[u][v] = true;
		edge_vector[v][u] = true;
	}
}

// u1, u2 are real ids
void Graph::deleteEdge(vid_t u1, vid_t u2)
{
	vid_t upper, lower; 
	if(same_layer(u1,u2)){
		fprintf(stderr, "Bad input for delete edge: u v on the same layer\n");
		exit(1);		
	}
	for (int i = 0; i < degree[u1]; ++i)
	{
		int vv = neighbor[u1][i];
		if (vv == u2)
		{	
			swap(neighbor[u1][i], neighbor[u1][degree[u1] - 1]);
			--degree[u1];
			neighbor[u1].pop_back();
			num_edges--;
			break;
		}
	}
	for (int i = 0; i < degree[u2]; ++i)
	{
		int uu = neighbor[u2][i];
		if (uu == u1)
		{
			swap(neighbor[u2][i], neighbor[u2][degree[u2] - 1]);
			--degree[u2];
			neighbor[u2].pop_back();
			break;
		}
	}
}
// this function tests the validity of the deleteEdge function.

