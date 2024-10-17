#include "counting.h"
#include "mt19937ar.h"

using namespace std;

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

long double efficient_wedge_based_triangle(Graph& g, unsigned long seed) {
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

	Graph g2(g);


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
		if(budget_optimization){
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
long double compute_w_uv(int u, int v, Graph& g, Graph& g2, bool sampling_noisy_graph) {
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

    if (avg_switch_triangle) {
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
void randomized_response_single_bit(int u, int v, Graph& g, Graph& g2) {
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

void private_estimate_of_degrees(Graph& g){
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

void construct_noisy_graph(Graph& g, Graph& g2, unsigned long seed){

	// if(g.edge_vector[i][j]){ // looks like this is slower


	if(g.is_bipartite){
		// cout<<"is Graph"<<endl;
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
double locally_compute_f_given_q_and_x(int q, int x, Graph& g, Graph& g2) {

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

unsigned long long int triangle_count(Graph& g) {
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

            if (!count_per_vertex && (i >= j)) 
				continue;

			for (auto k: neighbors_i) {

				if (j >= k) continue;

                const auto& neighbors_j = g.neighbor[j];
                if (std::find(neighbors_j.begin(), neighbors_j.end(), k) != neighbors_j.end()) {
					#pragma omp critical
                    tri_num++;
					if(count_per_vertex){
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

	if(count_per_vertex){
		tri_num/=3;
	}
	if(count_per_edge){
		cout<<"# of triangle edge "<<edge_in_triangles <<endl;
	}

	return tri_num;

}