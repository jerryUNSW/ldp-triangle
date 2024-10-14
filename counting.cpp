#include "counting.h"
#include "mt19937ar.h"

using namespace std;

extern long int real_val;

extern bool budget_optimation; 

// per-ver-tri-cnt ground truth
extern vector<int> tri_cnt; 

extern vector<vector<double>> tri_cnt_estis;

// Single-ton based pruning
bool skip_single_ton = true; 

extern vector<int> high_degree_vertices;

// Sampling parameters
vector<int> upper_sample, lower_sample;
double sample_ratio = 1.0;
bool sampling_one_round = false;

// Privacy-related parameters
extern long double Eps, Eps0, Eps1, Eps2, p; 
extern vector<double> priv_deg; 
extern vector<long double> naive_estis;
extern int priv_dmax_1, priv_dmax_2; 
extern long double communication_cost; 
extern double gamma__;

extern vector<unordered_map<int, bool>> sup_vector; 

// Algorithm iteration and epsilon parameters
extern int iteration;
extern double d1, d2, epsilon;

extern int deg_threshold; 

// Random engine initialization (seeded with current time)
stats::rand_engine_t engine(std::time(0)); // Used to be 1776

// Clipping parameters
bool edge_clipping = true;

// int alpha = 10; // Alternative: int alpha = 150;

int alpha = 150; // Alternative: int alpha = 150;

bool clipping = false; // New switch for enabling/disabling clipping

// Vertex parameters
extern int q_vertex, x_vertex; 
extern long double RR_time, server_side_time, naive_server_side;
bool avg_switch___triangle  = true;
bool trim_vertices_triangle = false; // Not very effective

// Algorithm selection and sampling parameters
extern int algo_switch;


// evaluate A' usage:
bool noisy_matrix_usage = false ;

// noisy edge sampling
bool sampling_noisy_graph = false;
double p____ = 0.1; // Noisy edge sampling ratio
// VC is very sensitive to noisy edge sampling

bool sampling_vertex_pairs ;
double vertex_ratio = 0.1; // Vertex pair sampling ratio

// Notes on vertex ratio:
// - On wiki, when vertex_ratio = 0.01, does not outperform SOTA
// - Let's try 0.1
// - Tuning this is effective in changing efficiency, as vertex pairs dominate the cost

// Edge-centric vs. vertex-centric processing
// bool edge_centric = false; 

bool use_high_deg_vertices = false;
// true: per-edge count (too slow), maybe it is too slow on sparse graphs
// false: per-vertex count




long double solve_equation(long double epsilon, long double S1, long double S2) {
	double target = 0;
	auto target_func = [S1, S2, epsilon, &target](double x) {
		// f0 is the equation to be solved
		double f0 =  -(S1*pow(epsilon - x, 3)*pow(exp(x) - 1, 2) + 4*S1*pow(epsilon - x, 3)*exp(x) 
			- 4*S2*(epsilon - x)*(exp(x) - 1)*exp(x) 
			+ 8*S2*(epsilon - x)*exp(2*x) - 4*S2*(exp(x) - 1)*(exp(x) + 1)*exp(x))*exp(x)/(pow(epsilon - x, 3)*pow(exp(x) - 1, 3)*(exp(x) + 1)); 

		// f1 is the derivative of f0
		double f1 =  (1.0/4.0)*(2*S1*pow(epsilon - x, 4)*(exp(x) - 2)*pow(exp(x) - 1, 2)*exp(x) 
		- S1*pow(epsilon - x, 4)*pow(exp(x) - 1, 3)*(exp(x) + 1) 
		- 4*S1*pow(epsilon - x, 4)*(exp(x) - 1)*(exp(x) + 1)*exp(x) 
		+ 8*S1*pow(epsilon - x, 4)*(exp(x) - 1)*(2*exp(x) - 1)*exp(x) 
		+ 24*S1*pow(epsilon - x, 4)*exp(2*x) + 16*S2*pow(epsilon - x, 2)*(exp(x) - 2)*(exp(x) - 1)*exp(2*x) 
		+ 4*S2*pow(epsilon - x, 2)*pow(exp(x) - 1, 2)*(exp(x) + 1)*exp(x) 
		- 4*S2*pow(epsilon - x, 2)*pow(exp(x) - 1, 2)*(2*exp(x) - 1)*exp(x) 
		- 8*S2*pow(epsilon - x, 2)*(exp(x) - 1)*(exp(x) + 1)*exp(2*x) 
		+ 48*S2*pow(epsilon - x, 2)*exp(3*x) + 16*S2*(epsilon - x)*pow(exp(x) - 1, 2)*(exp(x) + 1)*exp(x) 
		- 32*S2*(epsilon - x)*(exp(x) - 1)*(exp(x) + 1)*exp(2*x) 
		+ 12*S2*pow(exp(x) - 1, 2)*pow(exp(x) + 1, 2)*exp(x))/(pow(epsilon - x, 4)*pow(exp(x) - 1, 4)*pow(cosh((1.0/2.0)*x), 2));
		return std::pair<double, double>(f0 - target, f1);
	};
	double guess = 0.5*epsilon;
	double min = 0.01 * epsilon;
	double max = 0.99 * epsilon;
	int sf = 10; // # of significant figures

	double root = boost::math::tools::newton_raphson_iterate(target_func, guess, min, max, sf);
    return root; // returns the best epsilon_1
}


// need to implement a function to combine both methods. 
// hybrid dont help a lot on email dataset, maybe still not big enough

// long double hybrid_triangle(BiGraph& g, unsigned long seed) {

// 	double t1 = omp_get_wtime();


// 	// Phase 1: use 0.05 budget to estimate vertex degree and divide high and low 
// 	double Eps01 = Eps * 0.2;
// 	cout << "epsilon 0 = " << Eps01 << endl;
// 	priv_deg.resize(g.num_nodes());
// 	for(int i=0;i<g.num_nodes();i++){
// 		priv_deg[i] = g.degree[i]+stats::rlaplace(0.0, 1/(Eps01), engine); 
// 		if(edge_clipping){
// 			priv_deg[i]+=alpha;
// 		}
// 	}
//     // Phase 2: noisy graph 
//     Eps1 = Eps * 0.4;
//     cout << "epsilon 1 = " << Eps1 << endl;
//     p = 1.0 / (exp(Eps1) + 1.0);
//     BiGraph g2(g);

    
//     // Phase 2: counting
//     Eps2 = Eps  - Eps1 - Eps01;
//     cout << "epsilon 2 = " << Eps2 << endl;
//     gamma__ = (1 - p) / (1 - 2 * p);
// 	// cout<<"Count Type I triangles among high vertices "<<endl;

// 	long double cnt_1 =0;


// 	int degreeThreshold = 10000;
// 	cout<<"deg threshold = "<<degreeThreshold <<endl;
// 	// cout<<"Divide vertices into high and low degree"<<endl;
//     std::vector<int> highDegreeVertices;
//     for (int i = 0; i < g.num_nodes(); ++i) {
//         if (priv_deg[i] >= degreeThreshold) {
//             highDegreeVertices.push_back(i);
//         } 
//     }
// 	cout<<"# hi = "<<highDegreeVertices.size()<<endl;
// 	cout<<"% hi = "<<highDegreeVertices.size()*1.0 /g.num_nodes()  <<endl;
// 	cout<<"% lo = "<< 1- highDegreeVertices.size()*1.0 /g.num_nodes()  <<endl;
// 	// this part is not concerned with H and L.

// 	/*
// 	double pair_sampling_thred = 1 ; 

// 	for (int i = 0; i < highDegreeVertices.size(); ++i) {
// 		for (int j = i + 1; j < highDegreeVertices.size(); ++j) {

// 			int u = highDegreeVertices[i]; 
// 			int v = highDegreeVertices[j]; 

// 			long double triangle_support = 0;
// 			if (!g2.has_computed(u, v)){
// 				randomized_response_single_bit(u, v, g, g2);
// 			}
// 			triangle_support = ((g2.edge_vector[min(u, v)][max(u, v)] ? 1 : 0) - p ) / (1 - 2 * p);

// 			// (2) estimate the (common neighbors of u, v) \cap high degree vertices. 
// 			// compute from neighborhood of u.
// 			long double w_u_v = 0;
// 			int Nx_cap_Nq_minus_x = 0, Nq_minus_Nx_minus_x = 0;
// 			for (auto nb : g.neighbor[u]) {
// 				// use priv degrees because otherwise the triangles dont add up.
// 				if(priv_deg[nb]< degreeThreshold)
// 					continue;
// 				if (nb != v) {
// 					if (!g2.has_computed(nb, v)){
// 						randomized_response_single_bit(nb, v, g, g2);
// 					}
// 					g2.edge_vector[min(nb, static_cast<unsigned int>(v))][max(nb, static_cast<unsigned int>(v))] 
// 						? Nx_cap_Nq_minus_x++ 
// 						: Nq_minus_Nx_minus_x++;

// 				}
// 			}
// 			w_u_v = Nx_cap_Nq_minus_x * gamma__ + Nq_minus_Nx_minus_x * (-p) / (1 - 2 * p);
// 			w_u_v += stats::rlaplace(0.0, (gamma__ / Eps2), engine);

// 			// compute from neighborhood of v.
// 			long double w_v_u = 0;
// 			int count_common_neighbors = 0, count_different_neighbors = 0;
// 			for (auto nb : g.neighbor[v]) {
// 				if(priv_deg[nb] < degreeThreshold)
// 					continue;
// 				if (nb != u) {
// 					if (!g2.has_computed(nb, u)){
// 						randomized_response_single_bit(nb, u, g, g2);
// 					}
// 					g2.edge_vector[min(nb, static_cast<unsigned int>(u))][max(nb, static_cast<unsigned int>(u))] 
// 						? count_common_neighbors++ 
// 						: count_different_neighbors++;
// 				}
// 			}
//             w_v_u = count_common_neighbors * gamma__ + count_different_neighbors * (-p) / (1 - 2 * p);
//             w_v_u += stats::rlaplace(0.0, (gamma__ / Eps2), engine);

// 			// take the avg 
// 			w_u_v = (w_u_v + w_v_u) / 2;

// 			// aggregate results 
// 			triangle_support *= w_u_v;
// 			cnt_1 += triangle_support;
// 		}
// 	}
	
// 	cnt_1 /= pair_sampling_thred;
// 	cnt_1 /= 3;
// 	*/


// 	// looks like just visiting the heavy vertices will be pretty good. 


// 	// this is incorrent. 
// 	// we actually need to 

// 	double real_type_one_vertices = 0 ; 
// 	/*
//     std::unordered_set<int> highDegreeSet(highDegreeVertices.begin(), highDegreeVertices.end());
// 	for (int t = 0; t < highDegreeVertices.size(); ++t) {
// 		int i = highDegreeVertices[t]; 
//         const auto& neighbors_i = g.neighbor[i];
// 		std::vector<int> intersection;
// 		for (int neighbor : neighbors_i) {
// 			if (highDegreeSet.find(neighbor) != highDegreeSet.end()) {
// 				intersection.push_back(neighbor);
// 			}
// 		}
// 		for (int j : intersection) {
//             if (i >= j)  continue;
// 			for (auto k: intersection) {
// 				if (j >= k) continue;
//                 const auto& neighbors_j = g.neighbor[j];
//                 if (std::find(neighbors_j.begin(), neighbors_j.end(), k) != neighbors_j.end()) {
//                     real_type_one_vertices++;
//                 }
//             }
//         }
//     }
// 	*/

// 	// can we grab the real type I triangle count? 


// 	// 
// 	long double sum___ = 0; 
// 	// sum___ += cnt_1; 

// 	// we still need to compute somthing from heavy vertices. 
// 	// consider the triangles containing the heavy vertices. 
// 		// 3 heavy: type 1
// 		// 2 heavy and one light: type 2
// 		// 1 heavy and two light: type 3



// 	// we use heavy_deg[i] and light_deg[i] to control which vertices to use

// 	// long double cnt_1_2 = 0 ;
// 	for(int u=0;u<g.num_nodes();u++){

// 		if (edge_clipping && priv_deg[u] <=1) {
// 			continue;
// 		}

// 		int degree_u = g.degree[u]; 
// 		if(edge_clipping && degree_u > priv_deg[u]){
// 			degree_u = priv_deg[u]; 
// 		}

// 		long double cnt_1_ =0, cnt_2 =0, cnt_3 =0, cnt_4 =0; 
		
// 		// for (int i = 0; i < clipped_nb.size(); ++i) {
// 			// int x = clipped_nb[i];

// 		if (priv_deg[u] >= degreeThreshold) {
// 			// heavy vertex 
// 			for (int i = 0; i < degree_u; ++i) {
// 				int x = g.neighbor[u][i];
// 				for (int j = i + 1; j < degree_u; ++j) {
// 					int w = g.neighbor[u][j];
// 					if (!g2.has_computed(x, w)){
// 						randomized_response_single_bit(x, w, g, g2);
// 					}
// 					double noisy_x_w = g2.edge_vector[min(x,w)][max(x,w)] ? 1.0 : 0.0; 
// 					noisy_x_w = (noisy_x_w - p * p____)/ (p____ * (1 - 2 * p));

// 					// Type 1 triangle
// 					if(priv_deg[x] >= degreeThreshold && priv_deg[w] >= degreeThreshold ){
// 						cnt_1_ += noisy_x_w; 
// 					}
// 					// Type 2 triangle
// 					if(priv_deg[x] >= degreeThreshold && priv_deg[w] < degreeThreshold ){
// 						cnt_2 += noisy_x_w; 
// 					}
// 					if(priv_deg[x] < degreeThreshold && priv_deg[w] >= degreeThreshold ){
// 						cnt_2 += noisy_x_w; 
// 					}
// 					// Type 3 triangle
// 					if(priv_deg[x] < degreeThreshold && priv_deg[w] < degreeThreshold ){
// 						cnt_3 += noisy_x_w; 
// 					}
// 				}
// 			}	

// 			// Laplace Mechanism
// 			double GS = 0; 
// 			int deg_up = edge_clipping ? priv_deg[u] : g.degree[u];

// 			GS = (1 - p ) * deg_up /  (1 - 2 * p); 

// 			// assert(GS>0);
// 			// sum___ += 2*cnt_1_/3 + cnt_2/2 + cnt_3 ;

// 			sum___ += cnt_1_/3 + cnt_2/4 + cnt_3/2 ;

// 			// sum___ += cnt_1_*4 + cnt_2*3 + cnt_3*6 ;

// 			GS/=2;
// 			sum___ += stats::rlaplace(0.0,  GS / Eps2, engine);


// 		}else{
// 			// light vertex
// 			for (int i = 0; i < degree_u; ++i) {
// 				int x = g.neighbor[u][i];
// 				for (int j = i + 1; j < degree_u; ++j) {
// 					int w = g.neighbor[u][j];

// 					if (!g2.has_computed(x, w)){
// 						randomized_response_single_bit(x, w, g, g2);
// 					}
// 					// double noisy_x_w = g2.edge_vector[min(x,w)][max(x,w)] ? 1.0 : 0.0; 
// 					double noisy_x_w; 
					
// 					if(g2.memory_efficient){
// 						noisy_x_w= g2.__edge_vector[min(x,w)][max(x,w)].exists ? 1.0 : 0.0;
// 						g2.__edge_vector[min(x,w)][max(x,w)].count ++;

// 						// this is slightly useful
// 						if(g2.__edge_vector[min(x,w)][max(x,w)].count == min(g.degree[x], g.degree[w])){
// 							g2.__edge_vector[min(x,w)].erase(max(x,w)); 
// 						}
// 					}else{
// 						noisy_x_w= g2.edge_vector[min(x,w)][max(x,w)] ? 1.0 : 0.0;
// 					}


// 					noisy_x_w = (noisy_x_w - p )/  (1 - 2 * p);

// 					// Type 2 = 0 
// 					if(priv_deg[x] >= degreeThreshold && priv_deg[w] >= degreeThreshold ){
// 						cnt_2 += noisy_x_w; 
// 					}
// 					// Type 3 = 0 
// 					if(priv_deg[x] >= degreeThreshold && priv_deg[w] < degreeThreshold ){
// 						cnt_3 += noisy_x_w; 
// 					}
// 					if(priv_deg[x] < degreeThreshold && priv_deg[w] >= degreeThreshold ){
// 						cnt_3 += noisy_x_w; 
// 					}
// 					// Type 4 = 0, this is basically all triangles
// 					if(priv_deg[x] < degreeThreshold && priv_deg[w] < degreeThreshold ){
// 						cnt_4 += noisy_x_w; 
// 					}
// 				}
// 			}	
// 			// Laplace Mechanism
// 			double GS = 0; 
// 			int deg_up = edge_clipping ? priv_deg[u] : g.degree[u];

// 			GS = (1 - p ) * deg_up /  (1 - 2 * p); 

// 			// assert(GS>0);
// 			// assert(cnt_2 ==0); 
// 			// assert(cnt_3 ==0); 

// 			sum___ += cnt_2/2 + cnt_3/4 + cnt_4/3 ;

// 			// sum___ += cnt_2*6 + cnt_3*3 + cnt_4*4 ;

// 			GS/=3;
// 			sum___ += stats::rlaplace(0.0,  GS / Eps2, engine);		
// 		}
// 	}

// 	// need to evaluate the relative error of all sub_estimators.
// 	// type 1 triangle:
// 	// cout<<"Real Type I = "<<real_type_one_vertices <<endl;

// 	// cout<<"PC: Type I error = " << abs(cnt_1 - real_type_one_vertices) * 1.0 / real_type_one_vertices <<endl;
// 	// cout<<"estimate = "<<cnt_1 <<endl;
// 	// // cout<<endl;
// 	// cout<<"VC: Type I error = "<< abs(cnt_1_2 - real_type_one_vertices) * 1.0 / real_type_one_vertices <<endl;
// 	// cout<<"estimate = "<<cnt_1_2 <<endl;
// 	// cout<<"\% of Type 1  = "<<real_type_one_vertices / real_val<<endl;




// 	double t2 = omp_get_wtime();
// 	cout << "running time: " << t2 - t1 << endl;

// 	// sum___ /=12 ;

// 	return sum___; 
// }


long double efficient_wedge_based_triangle(BiGraph& g, unsigned long seed) {
    cout << "vertex_ratio = " << vertex_ratio << endl;

    // Phase 0: Initialize
	double t1 = omp_get_wtime();

	if (algo_switch==1) {
		// PC algorithm
		Eps0 = Eps * 0.0;
		// Phase 1: noisy graph construction
		cout << "epsilon 1 = " << Eps1 << endl;
		p = 1.0 / (exp(Eps1) + 1.0);
		gamma__ = (1 - p) / (1 - 2 * p);
	}else{
		cout << "epsilon 0 = " << Eps0 << endl;
		cout<<"edge clipping = "<<alpha <<endl;
		private_estimate_of_degrees(g); 
	}

	BiGraph g2(g);


	if(g2.memory_efficient){
		g2.__edge_vector.resize(g2.num_nodes());
	}else{
		g2.edge_vector.resize(g2.num_nodes());
	}
	
	
    double t2 = omp_get_wtime();


    long double sum___ = 0;


	long double total_vertex_pairs = static_cast<long double>(g.num_nodes()) * 
										static_cast<long double>(g.num_nodes() - 1) / 2.0;
										
    long double num_computations = 0;
    int last_percent = 0;

    // Choose sampling method
    if (!sampling_noisy_graph) {
        p____ = 1;
    }else{
		cout<<"edge sampling = "<<p____ <<endl;
	}
	
	/*
    if (algo_switch==1) {
        cout << "Pair Centric Algorithm" << endl;


        cout << "vertex pair sampling ratio = " << vertex_ratio << endl;

		// need to recompute high_degree_vertices. 

		if(use_high_deg_vertices){
			cout << "Recomputing high_degree_vertices using private degrees." << endl;
			high_degree_vertices.clear();
			for (int i = 0; i < g.num_nodes(); i++) {
				if (priv_deg[i] >= deg_threshold) {
					high_degree_vertices.push_back(i);
				}
			}
			cout << "Num of high degree vertices = " << high_degree_vertices.size() << endl;
		}

		

		

		size_t N = (use_high_deg_vertices ? high_degree_vertices.size() : g.num_nodes()); 

		total_vertex_pairs = N*(N-1) * vertex_ratio/2; 


		if(total_vertex_pairs<1){
			cout<<"total_vertex_pairs ==0"<<endl;
			exit(1);
		}else{
			cout<<"total_vertex_pairs = "<<total_vertex_pairs<<endl;
			// exit(1);
		}

		int next_percentage = 1;  // Start at 10% for the first print


	#pragma omp parallel
	{
		#pragma omp for schedule(static)
		for (size_t i = 0; i < N; ++i) {

			for (size_t j = i + 1; j < N; ++j) {

				int u = use_high_deg_vertices ? high_degree_vertices[i] : i;
				int v = use_high_deg_vertices ? high_degree_vertices[j] : j;
				long double local_res = 0;


				// proceed with vertex pair with a fixed probability
				if (sampling_vertex_pairs && genrand_real2() >= vertex_ratio) {
					continue;
				}

				#pragma omp critical
				{
					if (!g2.has_computed(u, v)){
						randomized_response_single_bit(u, v, g, g2);
					}
				}
				
				local_res = ((g2.edge_vector[min(u, v)][max(u, v)] ? 1 : 0) - p * p____) / (p____ * (1 - 2 * p));

		        // Compute w_u_v
				long double w_u_v; 
				#pragma omp critical
				{
					w_u_v = compute_w_uv(u, v, g, g2, sampling_noisy_graph);
				}
		        local_res *= w_u_v;

		        sum___ += local_res;
		    }
		}
	}
    } 
	*/
	
	if (algo_switch==2) {
		cout << "vertex-centric algorithm" << endl;

		// after the first round, we can re-assign the epsilon 1 and epsilon2 .
		// this can often find almost the best budget allocation strategy 
		if(budget_optimation){
			long double S1 =0, S2 = 0; 
			for(int i=0;i<g.num_nodes();i++){
				long double deg_u = priv_deg[i]; 
				S1 += (deg_u*(deg_u-1)/2); 
				S2 += deg_u * deg_u;
			}
			// cout<<"remain = "<<"Eps - Eps0 = "<<Eps - Eps0 <<endl;
			// cout<<"Eps0 = "<<Eps0<<endl;
			cout<<"S1 = "<<S1 <<endl;
			cout<<"S2 = "<<S2 <<endl;
			Eps1 = solve_equation(Eps - Eps0, S1, S2); 
			Eps2 = Eps - Eps0 - Eps1; 
		}

		cout<<"Eps0 = "<<Eps0<<endl;
		cout<<"Eps1 = "<<Eps1<<endl;
		cout<<"Eps2 = "<<Eps2<<endl;
		
		p = 1.0 / (exp(Eps1) + 1.0);
		gamma__ = (1 - p) / (1 - 2 * p);
		
		cout<<"flipping = "<<p<<endl;
		cout<<"gamma = "<<gamma__<<endl;

        for (int u = 0; u < g.num_nodes(); ++u) {

			// proceed with some probability
			// cout<<"// proceed with some probability"<<endl;
            if (sampling_vertex_pairs && genrand_real2() >= vertex_ratio) {
                continue;
            }
            // priv_deg[u] = g.degree[u] + stats::rlaplace(0.0, 1 / (Eps0), engine);
            if (edge_clipping) {
                if (priv_deg[u] <= 0) {
                    continue;
                }
            }

            long double local_result = 0;
            int deg_up = edge_clipping ? priv_deg[u] : g.degree[u];

            // Compute local result
			// question: should I add some randomness in edge-clipping
			
			int degree_u = g.degree[u]; 
			if(edge_clipping && degree_u > priv_deg[u]){
				degree_u = priv_deg[u]; 
			}
			// cout<<"visiting neighborhood of u "<<endl;
            for (int i = 0; i < degree_u; ++i) {
                for (int j = i + 1; j < degree_u; ++j) {
					int x = g.neighbor[u][i];
                    int w = g.neighbor[u][j];

					if (!g2.has_computed(x, w)){
						// cout<<"computing noisy edge for "<<x<<"\t"<<w<<endl;
                        randomized_response_single_bit(x, w, g, g2);
						num_computations++;
						// need to record these u, w 
                    }

					double noisy_x_w; 
					if(g2.memory_efficient){
						// noisy_x_w = g2.amat[max(x, w)][min(x, w)]; // should be either 1 or 0
						noisy_x_w= g2.__edge_vector[min(x,w)][max(x,w)].exists ? 1.0 : 0.0;
						g2.__edge_vector[min(x,w)][max(x,w)].count ++;

						// this is slightly useful
						if(g2.__edge_vector[min(x,w)][max(x,w)].count == min(g.degree[x], g.degree[w])){
							g2.__edge_vector[min(x,w)].erase(max(x,w)); 
						}
					}else{
						noisy_x_w= g2.edge_vector[min(x,w)][max(x,w)] ? 1.0 : 0.0;
					}
					noisy_x_w = (noisy_x_w - p * p____)/ (p____ * (1 - 2 * p));

					local_result += noisy_x_w; 
                }
            }
		
            // Laplace Mechanism
			double GS = 0; 
			if(sampling_noisy_graph){
				GS = (1 - p * p____) * deg_up / (p____ * (1 - 2 * p)); 
				// GS = (1 - p * p____) * deg_up / ( 1 * (1 - 2 * p)); 
			}else{
				GS = (1 - p ) * deg_up /  (1 - 2 * p); 
			}
            long double noise = stats::rlaplace(0.0,  GS / Eps2, engine);

            local_result += noise;

            sum___ += local_result;
        }
	}
    double t3 = omp_get_wtime();
    cout << "running time: " << t3 - t1 << endl;

	// cout<<"usage of noisy graph  = "<<num_computations * 1.0 / total_vertex_pairs <<endl;
	if(noisy_matrix_usage){
		cout << scientific << setprecision(6) 
				<< "usage = " 
				<< num_computations * 1.0 / total_vertex_pairs 
				<< endl;
		exit(0);
	}


	if(sampling_vertex_pairs){
		sum___ /= vertex_ratio;
	}
    return sum___ / 3;
}

// Compute w_u_v helper function
long double compute_w_uv(int u, int v, BiGraph& g, BiGraph& g2, bool sampling_noisy_graph) {
    long double w_u_v = 0;
    int Nx_cap_Nq_minus_x = 0, Nq_minus_Nx_minus_x = 0;

    for (auto nb : g.neighbor[u]) {
        if (nb != v) {
            // if (g2.edge_vector[nb].find(v) == g2.edge_vector[nb].end()) { // if has not computed 
			if (!g2.has_computed(nb, v)){
                randomized_response_single_bit(nb, v, g, g2);
            }
			// g2.edge_vector[min(nb, v)][max(nb, v)] ? Nx_cap_Nq_minus_x++ : Nq_minus_Nx_minus_x++;
			g2.edge_vector[min(nb, static_cast<unsigned int>(v))][max(nb, static_cast<unsigned int>(v))] 
				? Nx_cap_Nq_minus_x++ 
				: Nq_minus_Nx_minus_x++;

        }
    }
    if (sampling_noisy_graph) {
        long double gamma_tmp = (1 - p * p____) / (p____ * (1 - 2 * p));
        w_u_v = Nx_cap_Nq_minus_x * gamma_tmp + Nq_minus_Nx_minus_x * (-p) / (1 - 2 * p);
        long double noise = stats::rlaplace(0.0, (gamma_tmp / Eps2), engine);
        w_u_v += noise;
    } else {
        w_u_v = Nx_cap_Nq_minus_x * gamma__ + Nq_minus_Nx_minus_x * (-p) / (1 - 2 * p);
        long double noise = stats::rlaplace(0.0, (gamma__ / Eps2), engine);
        w_u_v += noise;
    }

    if (avg_switch___triangle) {
        long double w_v_u = 0;
        int count_common_neighbors = 0, count_different_neighbors = 0;
        for (auto nb : g.neighbor[v]) {
            if (nb != u) {
				if (!g2.has_computed(nb, u)){
                    randomized_response_single_bit(nb, u, g, g2);
                }
				g2.edge_vector[min(nb, static_cast<unsigned int>(u))][max(nb, static_cast<unsigned int>(u))] 
					? count_common_neighbors++ 
					: count_different_neighbors++;
                // g2.edge_vector[min(nb, u)][max(nb, u)] ? count_common_neighbors++ : count_different_neighbors++;
            }
        }

        if (sampling_noisy_graph) {
            long double gamma_tmp = (1 - p * p____) / (p____ * (1 - 2 * p));
            w_v_u = count_common_neighbors * gamma_tmp + count_different_neighbors * (-p) / (1 - 2 * p);
            long double noise = stats::rlaplace(0.0, (gamma_tmp / Eps2), engine);
            w_v_u += noise;
        } else {
            w_v_u = count_common_neighbors * gamma__ + count_different_neighbors * (-p) / (1 - 2 * p);
            long double noise = stats::rlaplace(0.0, (gamma__ / Eps2), engine);
            w_v_u += noise;
        }

        w_u_v = (w_u_v + w_v_u) / 2;
    }
    return w_u_v;
}


// Function to perform rejection sampling to get approx (N choose 2) * vertex_ratio vertex pairs
std::vector<std::pair<int, int>> reject_sampling_pairs(int N) {
    std::set<std::pair<int, int>> sampled_pairs;
    std::vector<std::pair<int, int>> result;

    // Total number of possible pairs
	long double total_vertex_pairs = static_cast<long double>(N) * 
											static_cast<long double>(N - 1) / 2.0;
    // Target number of pairs to sample
    long double  target_pairs = (total_vertex_pairs * vertex_ratio);
    
	cout<<"target_pairs = "<<target_pairs <<endl;
    // Random number generator setup
    srand(static_cast<unsigned>(time(0)));

    int pairs_processed = 0; // Counter for processed pairs

    // Iterate over all possible vertex pairs
    for (int u = 0; u < N; ++u) {
		// if(skip_single_ton && degree_one_vertices.count(u)==0)
		// 	continue 
        for (int v = u + 1; v < N; ++v) {
            
			// if(skip_single_ton && degree_one_vertices.count(v)==0)
			// 	continue 
            // Print percentage of pairs processed
			// ++pairs_processed;
            // if (pairs_processed % (total_pairs / 10) == 0) { // Print every 10% of processing
            //     double percentage = (static_cast<double>(pairs_processed) / total_pairs) * 100;
            //     std::cout << "Processed " << percentage << "% of vertex pairs." << std::endl;
            // }

            // Decide whether to include the pair based on vertex_ratio
            if (genrand_real2() < vertex_ratio) {
                auto pair = std::make_pair(u, v);
                // Add the pair to the set if it's not already present
                if (sampled_pairs.find(pair) == sampled_pairs.end()) {
                    sampled_pairs.insert(pair);
                    result.push_back(pair);
                    // Stop if we've reached the target number of pairs
                    if (result.size() >= target_pairs) {
                        return result;
                    }
                }
            }
        }
    }

    return result;
}

// randomized responses do not use sampling. 
void randomized_response_single_bit(int u, int v, BiGraph& g, BiGraph& g2) {
	double keep_probability = g.has(u, v) ? 1 - p : p;
	if(sampling_noisy_graph){
		keep_probability *= p____;
	}
	// Randomly decide whether to keep the edge based on keep_probability
	assert(u!=v);

	if(g2.memory_efficient){
		// roll a dice
		// if(genrand_real2() < keep_probability){
		// 	g2.amat[max(u, v)][min(u, v)] = 1;
		// }else{
		// 	g2.amat[max(u, v)][min(u, v)] = 0;
		// }

		if(genrand_real2() < keep_probability){
			g2.__edge_vector[min(u, v)][max(u, v)].exists = true;
		}else{
			g2.__edge_vector[min(u, v)][max(u, v)].exists = false;
		}
		// initialize the count to be zero.
		g2.__edge_vector[min(u, v)][max(u, v)].count = 0;
	}else{
		g2.edge_vector[min(u, v)][max(u, v)] = (genrand_real2() < keep_probability);
	}
	// it is either 1 or 0
}

/*
// this version constructs the noisy graph first, making it much slower
long double wedge_based_triangle(BiGraph& g, unsigned long seed) {
	
    // Phase 0. deg_esti_time records the maximum degree perturbation time.
    Eps0 = Eps * 0 ;

    // Phase 1. RR
    double t1 = omp_get_wtime();
    cout << "construct_noisy_graph(g); " << endl;
    Eps1 = Eps * 0.7;  // maybe we need to vary this. 
	cout<<"epsilon 1 = "<<Eps1 <<endl;
	// need a better variance derived for triangle here to optimize the budget allocation 
    p = 1.0 / (exp(Eps1) + 1.0);
    BiGraph g2(g);
	construct_noisy_graph(g, g2, seed);  // upload noisy edges

	double t2 = omp_get_wtime();

	cout<<"noisy graph time: "<<t2-t1<<endl;

	// exit(1); // for now only try to make it faster

	// Phase 2: 
	Eps2 = Eps - Eps1 - Eps0;
	cout<<"epsilon 2 = "<<Eps2 <<endl;
	long double sum___ = 0;
	gamma__ = (1-p) / (1-2*p);

	int total_vertex_pairs = g.num_nodes()*(g.num_nodes()-1)/2; 

	// the time complexity is O(N^2) which is not very good. 
	// the common-neighbor-based approach the L2 loss is still very hard to evaluate. 
	// unless we impose the fact that for each (u < v), we only consider the common neighbors x such that u < v < x. 

	// num_needed_computation = 0; 

	// unordered_set<std::pair<int, int>> noisy_edge_set;
	unordered_set<std::pair<int, int>, PairHash> noisy_edge_set;
	int num_computations =0;

	#pragma omp parallel
	{
		#pragma omp for schedule(static)	
		for(int u=0; u<g.num_nodes();u++){
			for(int v=u+1; v<g.num_nodes();v++){

				// number of vertex pairs is (N choose 2) = O( N^2 )
				if(genrand_real2() >= vertex_ratio ) 
					continue; 
				
				long double local_res = 0; 
				
				// local_res  = ((g2.has(u, v) ? 1 : 0) - p)/(1-2*p); 
				if(sampling_noisy_graph){
					local_res  = ((g2.has(u, v) ? 1 : 0) - p*p____)/(p____*(1-2*p)); 

					if(u<v){
						noisy_edge_set.insert({u, v});
					}else{
						noisy_edge_set.insert({v, u});
					}
					num_computations++;
				}else{
					local_res  = ((g2.has(u, v) ? 1 : 0) - p)/(1-2*p); 
				}
				if(algo_switch==1){
					// local_res *= locally_compute_f_given_q_and_x(u, v, g, g2);
					// for each x in N(u)\v, push edge (x,v)
					long double w_u_v = 0;
					int Nx_cap_Nq_minus_x=0;
					int Nq_minus_Nx_minus_x=0;
					for(auto nb: g.neighbor[u]){
						if(nb != v){
							if(g2.has(nb, v)){
								Nx_cap_Nq_minus_x++;
							}else{
								Nq_minus_Nx_minus_x++;
							}
							if(nb<v){
								noisy_edge_set.insert({nb, v});
							}else{
								noisy_edge_set.insert({v, nb});
							}
							num_computations++;
						}
					}
					if(sampling_noisy_graph){
						long double gamma_tmp = (1 - p * p____) / (p____ * (1 - 2 * p));
						w_u_v = Nx_cap_Nq_minus_x * gamma_tmp + Nq_minus_Nx_minus_x * (-p) / (1 - 2 * p);
						long double noise = stats::rlaplace(0.0, (gamma_tmp/Eps2), engine); 
						w_u_v += noise; 
					}else{
						w_u_v = Nx_cap_Nq_minus_x * gamma__ + Nq_minus_Nx_minus_x * (-p)/(1-2*p); 
						long double noise = stats::rlaplace(0.0, (gamma__/Eps2), engine); 
						w_u_v += noise; 
					}					
					local_res *= w_u_v; 

				}
				if(algo_switch==2){// effectiveness slightly better but efficiency slower in general  
					local_res *= (locally_compute_f_given_q_and_x(u, v, g, g2)+ locally_compute_f_given_q_and_x(v, u, g, g2))/2 ; 
				}
				// (2) use just f1, but when considering common neighbors, just consider those x such that: u < v < x. 
				// try both 
				// local_res *= locally_compute_f_given_q_and_x(u, v, g, g2);
				// also need to estimate the number of common neighbors between u and v. 
				#pragma omp critical
				sum___ += local_res ; // this is getting all edges.
			}
		}
	}
	// cout<<"total pairs = "<<total_vertex_pairs <<endl;
	// cout<<"estimate = "<<sum___ <<endl;

	double t3 = omp_get_wtime();
	cout<<"counting time: "<<t3-t2<<endl;


	cout<<"# entries on A' = "<<g.num_nodes()*(g.num_nodes()-1)/2 <<endl;
	cout<<"# entries used in A' = "<<noisy_edge_set.size()<<endl;
	cout<<"\% entries used in A' = "<<noisy_edge_set.size()*1.0/(g.num_nodes()*(g.num_nodes()-1)/2)<<endl;

	// it is very weird that as I apply sampling, the number of needed entries increase. 

	// cout<<"num computations = "<<num_computations <<endl;

    // for (const auto& pair : noisy_edge_set) {
    //     std::cout << "(" << pair.first << ", " << pair.second << ")\n";
    // }

	if(sampling_vertex_pairs){
		sum___ /= vertex_ratio; 
	}
	// if (trim_vertices__triangle){
	// 	return sum___; 
	// }else{
	return sum___/3 ;
	// }
	
}
*/

void private_estimate_of_degrees(BiGraph& g){
	// Eps0 = 0.05;

	if(g.is_bipartite){
		// private estimate degrees. 
		priv_dmax_1 = 0 ;       
		priv_dmax_2 = 0 ; 
		priv_deg.resize(g.num_nodes());
		for(int i=0;i<g.num_nodes();i++){

			// priv_deg[i] = g.degree[i]+stats::rlaplace(0.0, 1/(Eps0), engine); 
			// need to re-run two-round algorithms with 1 --> 2 
			priv_deg[i] = g.degree[i]+stats::rlaplace(0.0, 1/(Eps0), engine); 
			
			if(edge_clipping){
				priv_deg[i]+=alpha;
			}
			// long double tmp = add_geometric_noise(g.degree[i], 1, Eps0);  
			// we found that geometric 
			// cout<<"deg = "<<g.degree[i]<<",\t"<<priv_deg[i]<<",\t"; cout<<tmp<<endl; 
			if(g.is_upper(i)){
				priv_dmax_1 = priv_dmax_1 > priv_deg[i] ? priv_dmax_1 : priv_deg[i] ; 
			}else{
				priv_dmax_2 = priv_dmax_2 > priv_deg[i] ? priv_dmax_2 : priv_deg[i] ; 
			}
			// if(choose_upper && g.is_upper(i))x1+=priv_deg[i];
			// if(!choose_upper && g.is_lower(i))x1+=priv_deg[i];
		}
	}else{
		// private estimate degrees. 
		priv_deg.resize(g.num_nodes());

		for(int i=0;i<g.num_nodes();i++){

			priv_deg[i] = g.degree[i]+stats::rlaplace(0.0, 1/(Eps0), engine); 
			
			if(edge_clipping){
				priv_deg[i]+=alpha;
			}
		}		
	}
}

double my_genrand_real2(){
	return genrand_real2(); 
}

void my_init_genrand(unsigned long seed){init_genrand(seed);}

void construct_noisy_graph(BiGraph& g, BiGraph& g2, unsigned long seed){

	// if(g.edge_vector[i][j]){ // looks like this is slower


	if(g.is_bipartite){
		// cout<<"is bigraph"<<endl;
		const int range_from  = g2.num_v1;
		const int range_to    = g2.num_nodes()-1;
		std::random_device                  rand_dev;
		std::mt19937                        generator(rand_dev());
		std::uniform_int_distribution<int>  distr(range_from, range_to);

		init_genrand(seed);
		int flip1 = 0; 
		int visited_vertices = 0;
		int total_vertices = g2.num_v1;
		int ten_percent = total_vertices / 5;
		long double max_time_per_user = -1;
		for(int i=0;i<g2.num_v1;i++){
			visited_vertices++; 
			// if (visited_vertices % ten_percent == 0) {
			// 	int progress = visited_vertices * 100 / total_vertices;
			// 	cout << "Processed " << progress << "% of vertices" << endl;
			// }
			double tx = omp_get_wtime();
			for(int j=g2.num_v1;j<g2.num_nodes();j++){
				if(std::find(g.neighbor[i].begin(), g.neighbor[i].end(), j) != g.neighbor[i].end() ){
					if(genrand_real2() >= p ){ // 1  --> 1 

						g2.addEdge(i,j);
						// #pragma omp atomic
						flip1++;
					}	
				}else{
					if(genrand_real2() < p){   // 0 --> 1

						g2.addEdge(i,j);
						// #pragma omp atomicss
						flip1++;	
					}
				}
			}
			double ty = omp_get_wtime();
			max_time_per_user = max_time_per_user > (ty-tx) ? max_time_per_user : (ty-tx);
		}
		cout<<"noisy edges = "<<flip1<<endl;

		communication_cost+= flip1*sizeof(int); 

		long double expected_E = g.num_edges*(1-p) +(g.num_v1*g.num_v2-g.num_edges)*p; 

		cout<<"expected E = "<< expected_E<<endl;

		// g2.computePriority();

	}else{
		
		// right now we still construct noisy graph by sampling edges. 
		// what if we sample vertices?? 
		if(!sampling_noisy_graph){
			cout<<"constructing a noisy graph for a GENERAL graph"<<endl;
		}else{
			cout<<"constructing a sampled noisy graph for a GENERAL graph with probability "<<p____<<endl;
		}
		// constructing a noisy graph for a general graph
		const int range_from  = 0;
		const int range_to    = g2.num_nodes()-1;
		std::random_device                  rand_dev;
		std::mt19937                        generator(rand_dev());
		std::uniform_int_distribution<int>  distr(range_from, range_to);
		init_genrand(seed);
		//
		int flip1 = 0; 
		int visited_vertices = 0;
		int total_vertices = g2.num_nodes();
		int ten_percent = total_vertices / 2;
		// long double max_time_per_user = -1;
		for(int i=0;i<g2.num_nodes();i++){
			visited_vertices++; 
			if (visited_vertices % ten_percent == 0) {
				int progress = visited_vertices * 100 / total_vertices;
				cout << "Processed " << progress << "% of vertices" << endl;
			}
			// double tx = omp_get_wtime();
			for(int j=0; j<i;j++){
				// the lower left triangle
				double keep_probability;
				if (std::find(g.neighbor[i].begin(), g.neighbor[i].end(), j) != g.neighbor[i].end()) {
					keep_probability = 1 - p;
				} else {
					keep_probability = p;
				}
				if (sampling_noisy_graph) {
					keep_probability *= p____;
				}
				if (genrand_real2() < keep_probability) {
					g2.addEdge(i, j); // O(1), this is used to add positive edges
					flip1++;
				}
			}
			double ty = omp_get_wtime();
			// max_time_per_user = max_time_per_user > (ty-tx) ? max_time_per_user : (ty-tx);
		}
		long double possible_edges = g.num_vertices*(g.num_vertices-1)/2; 

		// cout<<"# possible edges = "<<possible_edges <<endl;

		long double expected_E = g.num_edges*(1-p) +(possible_edges -g.num_edges) *p; 

		if(sampling_noisy_graph) expected_E*= p____;

		cout<<"expected E = "<< expected_E<<endl;

		cout<<"# noisy edges = "<<flip1<<endl;
	}
}

// if trim_vertices__triangle is true, then we only visit 
double locally_compute_f_given_q_and_x(int q, int x, BiGraph& g, BiGraph& g2) {

	// cout<<"q = "<<q <<endl; 
	// cout<<"x = "<<x <<endl; 	

	// cout<<"using Eps2 = "<<Eps2 <<endl;
	double res =-1;
	int start, end;
	if(g2.is_bipartite){
		start = ( g2.is_upper(q) ) ?   g2.num_v1 : 0 ;
		end =   ( g2.is_upper(q) ) ?   g2.num_nodes() : g2.num_v1;
	}else{
		start = 0;
		end = g2.num_nodes();
	}

	double Nx_cap_Nq_minus_x =0, Nq_minus_Nx_minus_x =0;
	// it looks like edge clipping makes it worse
	bool x_is_a_nb_of_q = false;


	// computing the set intersection and set difference 
	for(auto nb: g.neighbor[q]){
		if(nb == x){
			x_is_a_nb_of_q = true;// this cannot be true in bipartite case 
		}
		else{
			// at this point it is (Nq-x)
			if(g2.has(nb, x)){
				// computing Nx cap Nq -x 
				Nx_cap_Nq_minus_x++;
			}else{
				// computing Nq -Nx -x 
				Nq_minus_Nx_minus_x++;
			}
		}
	}
	
	// if x is not a neighbor of q, then 
	/*
	if(g2.is_bipartite){
		// in the case of bipartite graphs
		if(x_is_a_nb_of_q){
			cout<<"x cannot be a neighbor of q. They are from the same layer."<<endl;
			exit(1);
		}
	}else{
		// in the case of general graphs 
		// only check this when not trimming
		if(!trim_vertices__triangle){
			if(x_is_a_nb_of_q){
				assert(Nq_minus_Nx_minus_x == g.degree[q] - Nx_cap_Nq_minus_x-1);
			}else{
				assert(Nq_minus_Nx_minus_x == g.degree[q] - Nx_cap_Nq_minus_x); 
			}
		}
	}
	*/
	// locally the degree of q is known. 
	// res = Nx_cap_Nq_minus_x * gamma__ + Nq_minus_Nx_minus_x * (-p)/(1-2*p); 
	if(sampling_noisy_graph){
		long double gamma_tmp = (1 - p * p____) / (p____ * (1 - 2 * p));
		res = Nx_cap_Nq_minus_x * gamma_tmp + Nq_minus_Nx_minus_x * (-p) / (1 - 2 * p);
		// Laplacian noise
		long double noise = stats::rlaplace(0.0, (gamma_tmp/Eps2), engine); 
		res += noise; 
	}else{
		res = Nx_cap_Nq_minus_x * gamma__ + Nq_minus_Nx_minus_x * (-p)/(1-2*p); 
		// Laplacian noise
		long double noise = stats::rlaplace(0.0, (gamma__/Eps2), engine); 
		res += noise; 
	}

	return res;
}

long double get_laplace(long double parameter){
	return stats::rlaplace(0.0, parameter, engine); 
}

// let all u, w pairs share the same privacy budget allocations. 
// for each pairs of (u, w) in upper vertices. 
		// we individually allocated privacy budget to them.
		// in the end, we summarize the results. 
		// the running time will be much longer 


// we might also need to adopt some budget optimization strategy here. 
long double wedge_based_two_round_btf(BiGraph& g, unsigned long seed) {
	
    // Phase 0. deg_esti_time records the maximum degree perturbation time.
    // double t0 = omp_get_wtime();
    Eps0 = Eps * 0.1;
	vector<long double> deg_estis; 
	deg_estis.resize(g.num_nodes());
	for(int i=0;i<g.num_nodes();i++){
		deg_estis[i] = g.degree[i]+stats::rlaplace(0.0, 1/(Eps0), engine); 
	}
    // upload noisy degrees
    // if (eva_comm) communication_cost += g.num_nodes() * sizeof(int);

    // Phase 1. RR
    double t1 = omp_get_wtime();
    cout << "construct_noisy_graph(g); " << endl;
    Eps1 = Eps * 0.7;
    p = 1.0 / (exp(Eps1) + 1.0);
    BiGraph g2(g);
	construct_noisy_graph(g, g2, seed);  // upload noisy edges

    // Phase 2. local counting, 
	// each vertex u download the noisy graph. 
    double t2 = omp_get_wtime();
    cout << "local counting" << endl;
	Eps2 = Eps - Eps1 - Eps0;

	// cout<<"using Eps2 = "<<Eps2 <<endl;
	gamma__ = (1-p) / (1-2*p);
    // use eps2
    // long double global_sensitivity, sum = 0;
	long double res___ = 0; 

	bool vertex_pair_reduction = true; 
	bool averaging_f_estimates = true; 
	// what if we only consider upper vertices ?  --> better efficiency and effect  
	// #pragma omp parallel for reduction(+:res___)
	int start__, end__; 
	start__ = g.num_v1 < g.num_v2 ? 0 : g.num_v1; 
	end__ = g.num_v1 < g.num_v2 ? g.num_v1 : g.num_nodes(); 
	#pragma omp parallel
	{
	#pragma omp for schedule(static)
		for(int u =start__ ; u <end__ ; u++) {
			// if (omp_get_thread_num() == 0) {
			//     // This code runs only once in the parallel region to print the number of threads
			//     std::cout << "Number of threads: " << omp_get_num_threads() << std::endl;
			// }
			for(int w =start__ ; w <end__ ; w++) {
				if(u==w) 
					continue;
				if(vertex_pair_reduction && u<w) 
					continue; // only consider each pair once

				long double f_u_w = locally_compute_f_given_q_and_x(u, w, g, g2);
				if(averaging_f_estimates){
					// using the average of fu and fw.
					f_u_w += locally_compute_f_given_q_and_x(w, u, g, g2);
					f_u_w /= 2; 
				}
				long double local_res = f_u_w * f_u_w - f_u_w; 
				
				// compute the variance of tilde(f)
				long double esti_var_f = 2 * gamma__ * gamma__ / (Eps2 * Eps2); 
				esti_var_f += p * (1 - p) * deg_estis[u] / ((1 - 2 * p) * (1 - 2 * p));

				if(averaging_f_estimates){
					// computing the variance of f2, if the switch is on
					long double esti_var_f_2 = 2 * gamma__ * gamma__ / (Eps2 * Eps2); 
					esti_var_f_2 += p * (1 - p) * deg_estis[w] / ((1 - 2 * p) * (1 - 2 * p));	

					// take Var(\tilde(f)) = (var(f1) + var(f2))/4
					esti_var_f = (esti_var_f + esti_var_f_2)/4;
				}
				local_res -= esti_var_f; 
				// if the degree estimate of u is more accurate then it's better
				#pragma omp critical
				res___ += local_res; 
			}
		}
	}
	if (vertex_pair_reduction){
		return res___/2; 
	}else{
		return res___/4; 
	}

}

long double BFC_EVP(BiGraph& g) {
    long double BTF = 0;
#pragma omp parallel for reduction(+ : BTF)
    for (int u = 0; u < g.num_nodes(); u++) {
        if (g.degree[u] <= 1) continue;
        unordered_map<vid_t, int> count_wedge(0);
        for (auto v : g.neighbor[u]) {
            // u -> v
            for (auto w : g.neighbor[v]) {
                // u->v->w
                // this step is not invoked for g3.
                // if (g.com_p(w, v) & g.com_p(w, u)) {  // this is a lot faster.
				if ( w <v & w <u ) 
				{  // this is a lot faster.
                    count_wedge[w] = count_wedge[w] + 1;
                }
            }
        }
        // long double btf_u = 0;
        for (auto ele : count_wedge) {
            if (ele.second >= 2) {
                BTF += (ele.second - 1) * ele.second / 2;
            }
        }
    }
    return BTF;
}

bool count_per_ver = false;
bool count_per_edge = false;
unsigned long long int triangle_count(BiGraph& g) {
    unsigned long long int tri_num = 0;
	cout<<"Counting real triangle counts"<<endl;

	unsigned long long int edge_in_triangles = 0;

	#pragma omp parallel
	{
	#pragma omp for schedule(dynamic)
    for (int i = 0; i < g.num_nodes(); ++i) {
        const auto& neighbors_i = g.neighbor[i];

        // for (auto j_itr = neighbors_i.begin(); j_itr != neighbors_i.end(); ++j_itr) {
		for (int j : neighbors_i) {

            if (!count_per_ver && (i >= j)) 
				continue;

			for (auto k: neighbors_i) {

				if (j >= k) continue;

                const auto& neighbors_j = g.neighbor[j];
                if (std::find(neighbors_j.begin(), neighbors_j.end(), k) != neighbors_j.end()) {
					#pragma omp critical
                    tri_num++;
					if(count_per_ver){
						tri_cnt[i]++;
					}


					// i <j <k
					if(count_per_edge){
						if (!sup_vector[i][j]) {
							sup_vector[i][j] = true;
							edge_in_triangles++;
						}
						if (!sup_vector[i][k]) {
							sup_vector[i][k] = true;
							edge_in_triangles++;
						}
						if (!sup_vector[j][k]) {
							sup_vector[j][k] = true;
							edge_in_triangles++;
						}
					}

					// found triangle <i, j, k>
                }
            }
        }
    }
	}

	if(count_per_ver){
		tri_num/=3;
	}
	if(count_per_edge){
		cout<<"# of triangle edge "<<edge_in_triangles <<endl;
	}

	return tri_num;

}
// vector<unordered_map<int, bool>> edge_vector;

/*
long double weighted_sampling_triangle(BiGraph& g, unsigned long seed) {
	
	// conjecture: 
    // Define the threshold for high-degree vertices
    int degreeThreshold = 50;  // Example threshold, adjust as needed

    // Divide vertices into high and low degree based on the threshold
    std::vector<int> highDegreeVertices;
    std::vector<int> lowDegreeVertices;

    for (int i = 0; i < g.num_nodes(); ++i) {
        if (g.degree[i] >= degreeThreshold) {
            highDegreeVertices.push_back(i);
        } else {
            lowDegreeVertices.push_back(i);
        }
    }


	// maybe this is how we could combine VC and PC.
	// and 

    cout<<"// Count Type 1 Triangles: triangles among high-degree vertices"<<endl;
    // Count Type 1 Triangles: triangles among high-degree vertices
    // Iterate through each high-degree vertex
	int type1Triangles = 0;
    for (int u : highDegreeVertices) {
		for(auto v:g.neighbor[u]){
			if (std::find(highDegreeVertices.begin(), highDegreeVertices.end(), v) != highDegreeVertices.end() && v > u) {
                for (int w : g.neighbor[v]) {
					if (w > v && std::find(highDegreeVertices.begin(), highDegreeVertices.end(), w) != highDegreeVertices.end() && g.has(u, w)) {
                        ++type1Triangles;
                    }
                }
            }
        }
    }


    // Count Type 2 Triangles: triangles involving at least one low-degree vertex
    // Process each low-degree vertex
	int type2Triangles = 0;
	int type3Triangles = 0;
	int type4Triangles = 0;
    for (int i : lowDegreeVertices) {
        const auto& neighbors_i = g.neighbor[i];
        // Process pairs of neighbors of the low-degree vertex
        for (size_t j = 0; j < neighbors_i.size(); ++j) {
            int u = neighbors_i[j];
            // if (g.degree[u] >= degreeThreshold) continue; // Skip if neighbor is not a low-degree vertex

            for (size_t k = j + 1; k < neighbors_i.size(); ++k) {
                int v = neighbors_i[k];

				assert(u!=v);

                // Check if u and v are connected
                if (g.has(u, v)) {
                    // Type 2: Two high-degree, one low-degree
                    if (g.degree[u] >= degreeThreshold && g.degree[v] >= degreeThreshold) {
                        type2Triangles++;
                    }
                    // Type 3: One high-degree, two low-degree
                    else if (g.degree[u] >= degreeThreshold && g.degree[v] < degreeThreshold) {
						type3Triangles++;
                    }
                    // Type 3: One high-degree, two low-degree
                    else if (g.degree[u] < degreeThreshold && g.degree[v] >= degreeThreshold) {
						type3Triangles++;
                    }
                    // Type 4: All low-degree vertices
                    else {
                        type4Triangles++;
                    }
                }
            }
        }
    }

    // Output results
    std::cout << "Type 1 Triangles: " << type1Triangles << std::endl;
    std::cout << "Type 2 Triangles: " << type2Triangles << std::endl;
	std::cout << "Type 3 Triangles: " << type3Triangles << std::endl;
	std::cout << "Type 4 Triangles: " << type4Triangles << std::endl;
    std::cout << "Total Triangles: " << type1Triangles +  type2Triangles + type3Triangles/2 + type4Triangles/3 << std::endl;


	return -1;


    // Define thresholds for high and low degree vertices


	cout<<"weighted_sampling_triangle "<<endl;

	Eps0 = Eps * 0.1;
	cout << "epsilon 0 = " << Eps0 << endl;
	private_estimate_of_degrees(g); 

    cout << "epsilon 1 = " << Eps1 << endl;
    p = 1.0 / (exp(Eps1) + 1.0);
    BiGraph g2(g);

	g2.edge_vector.clear();

    g2.edge_vector.resize(g.num_v1 + g.num_v2);

    Eps2 = Eps - Eps1 - Eps0;
    cout << "epsilon 2 = " << Eps2 << endl;
    
    gamma__ = (1 - p) / (1 - 2 * p);

	// Now compare this direct sum with the computed S
	double degree_sqr_sum = 0.0;
	for (int i = 0; i < g.num_nodes(); ++i) {
		degree_sqr_sum += g.degree[i] * g.degree[i];  // Sum of squared degrees
	}
	double m = g.num_edges ; 
	double S = 4.0 * m * m - degree_sqr_sum;

	std::cout << "Computed S using 4m^2 - sum(deg^2): " << S << std::endl;


	// Calculate total degree sum using a for loop
	double total_degree_sum = 0.0;
	for (int i = 0; i < g.num_nodes(); ++i) {
		total_degree_sum += g.degree[i];
	}

	// Number of vertex pairs to sample
	int T = 1000000;
	std::map<int, std::map<int, int>> pair_frequency;

	// Sample T vertex pairs
	init_genrand(seed);
	for (int i = 0; i < T; ++i) {
		int u, v;
		do {
			// Sample u with probability proportional to its degree
			double rand_val_u = genrand_real2() * total_degree_sum;
			double cumulative_sum_u = 0.0;
			for (u = 0; u < g.num_nodes(); ++u) {
				cumulative_sum_u += g.degree[u];
				if (rand_val_u < cumulative_sum_u) break;
			}
			// Sample v with probability proportional to its degree
			double rand_val_v = genrand_real2() * total_degree_sum;
			double cumulative_sum_v = 0.0;
			for (v = 0; v < g.num_nodes(); ++v) {
				cumulative_sum_v += g.degree[v];
				if (rand_val_v < cumulative_sum_v) break;
			}
			
		} while (u == v);  // Reject if u == v

		pair_frequency[u][v]++;
	}

	// Output frequencies and compare with theoretical probabilities
	int cnt = 0;

	long double sum___ = 0;

	p____ =1 ;
	for (const auto& u_map : pair_frequency) {
		int u = u_map.first;
		for (const auto& v_pair : u_map.second) {
			int v = v_pair.first;
			int freq = v_pair.second;

			// here we have u, v freq.
			
			// number of triangles containing pair u, v
			
			long double local_res = 0;
			if (!g2.has_computed(u, v)){
				randomized_response_single_bit(u, v, g, g2);
			}
			local_res = ((g2.edge_vector[min(u, v)][max(u, v)] ? 1 : 0) - p ) /  (1 - 2 * p);

			long double w_u_v = compute_w_uv(u, v, g, g2, false);

			local_res *= w_u_v;
			
			sum___ += freq * local_res / (g.degree[u] * g.degree[v]*1.0);
			
			

			// long double local_res = 0;
			// if(g.has(u,v)){
			// 	// count the common neighbor of u and v
			// 	for (int xxx : g.neighbor[u]) {
			// 		if (find(g.neighbor[v].begin(), g.neighbor[v].end(), xxx) != g.neighbor[v].end()) {
			// 			local_res++;
			// 		}
			// 	}
			// }
			// sum___ += freq * local_res / (g.degree[u] * g.degree[v]*1.0);


			

			// sum += number of triangles containing u,v

			// double empirical_prob = static_cast<double>(freq) / T;
			// double theoretical_prob = (g.degree[u] * g.degree[v]*1.0) / S;
			// std::cout << "Pair (" << u << ", " << v << "): Emp. Prob = " << empirical_prob
			// 		<< ", Theo Prob = " << theoretical_prob;
			// cout<<"\t\tdu = "<<g.degree[u];
			// cout<<"\t dv = "<<g.degree[v]<<endl;
		}
	}
	sum___ *= S; 
	sum___ /= 6*T; // why is it not 3

	// if I combine this method with pair-wise method. 

	return sum___;
}
*/