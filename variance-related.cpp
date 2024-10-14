	////////////////////////////////////////////////////////////
	/*
	cout<<" |V| = "<<g.num_nodes()<<endl;
	long double num_pairs = g.num_nodes() * (g.num_nodes() -1)/2;

	cout<<"num pairs = "<<num_pairs <<endl;

	long double num_pair_of_pairs = num_pairs * (num_pairs -1)/2;
	cout<<"num of pairs of pairs "<< num_pair_of_pairs <<endl;


    // Example graph with nodes and edges
    int num_nodes = g.num_nodes();
    std::vector<std::pair<int, int>> pairs;
    // Generate all pairs (u, v) where u < v
    for (int u = 0; u < num_nodes; ++u) {
        for (int v = u + 1; v < num_nodes; ++v) {
            pairs.push_back({u, v});
        }
    }
    // Vector to store all unique pairs of pairs
	long double cnt = 0; 
	long double cnt2 = 0; 

    // Generate all unique pairs of pairs
    for (size_t i = 0; i < pairs.size(); ++i) {
        for (size_t j = i + 1; j < pairs.size(); ++j) {
            const auto& p1 = pairs[i];
            const auto& p2 = pairs[j];
            // Ensure that the pairs are distinct
			if (p1 != p2) {
				int u1 = p1.first;
				int v1 = p1.second;
				int u2 = p2.first;
				int v2 = p2.second;

				// Store the vertices in a vector
				std::vector<int> cycle = {u1, v1, u2, v2};
				// think about how the 4-cycles are rotated 

				if(u1== u2){
					if(v1<v2){
						cnt++;
					}
					if(v1>v2){
						cnt++;
					}
				}
				if(u1 < u2){
					if(v1< u2 && u2 < v2){
						cnt++;
						// need to check u1, v1, u2, v2 is a four cycle? 
						if(g.has(u1, v1)  && g.has(v1, u2) && g.has(u2, v2) && g.has(v2, u1) ){
							cnt2++;
				
							// Print the sorted 4-cycle
							cout<<"cycle 1: ";
							for (int vertex : cycle) {
								std::cout << vertex << " ";
							}
							std::cout << std::endl;
						}
					}
					if(v1== u2 && u2 < v2){
						cnt++;
					}
					if(u2< v1 && v1 < v2){
						cnt++;
						if(g.has(u1, v1)  && g.has(v1, u2) && g.has(u2, v2) && g.has(v2, u1) ){
							cnt2++;
							cout<<"cycle 2: ";
							for (int vertex : cycle) {
								std::cout << vertex << " ";
							}
							std::cout << std::endl;
						}
					}
					if(u2 < v1 && v2 == v1){
						cnt++;

						// need to consider this as well. 
						// if(g.has(u1, v1) && g.has(u2, v1) ){
						// 	for(auto xxx : g.neighbor[u1]){
						// 		if(xxx==v1) continue;
						// 		if(g.has(u2,xxx)){
						// 			cnt2++;
						// 		}
						// 	}
						// }

					}
					if(u2 < v2 && v2 < v1){

						// cnt++;
						// if(g.has(u1, v1)  && g.has(v1, u2) && g.has(u2, v2) && g.has(v2, u1) ){
						// 	cnt2++;
						// 	cout<<"cycle 3: ";
						// 	for (int vertex : cycle) {
						// 		std::cout << vertex << " ";
						// 	}
						// 	std::cout << std::endl;
						// }
					}
				}
				// cnt++;
			}
        }
    }
	// it looks like different pairs of <u, v> may correspond to the same 4-cycle. 
    cout << "Number of unique pairs of pairs: " << cnt<< std::endl;

	cout<< "count 2 = " << cnt2 <<endl;

    int four_cycle_count = 0;
    // Iterate through all pairs of vertices (u, v) where u < v
    for (int u = 0; u < g.num_nodes(); ++u) {
        for (int v = u + 1; v < g.num_nodes(); ++v) {
            // Find common neighbors of u and v
            vector<int> common_neighbors;
			long double common_nb_size = 0;
            for (auto w : g.neighbor[v]) {
				if ( g.has(u,w) ) {
					common_nb_size++;
                }
            }
			four_cycle_count += common_nb_size * (common_nb_size-1)/2;
        }
    }
	four_cycle_count /=2; 
    cout << "Number of 4-cycles: " << four_cycle_count << endl;


	exit(1);



	cout<<"Lower bound of Variance computation: "<<endl;
	long double esp1__ = Eps * 0.7;
	long double p__ = 1.0 / (exp(esp1__) + 1.0);
	long double esp2__ = Eps - esp1__; 

	long double gamma____ = (p__ * (1 - p__)) / pow((1 - 2 * p__), 2);
	long double theta____ = pow((1 - p__), 2) / pow((1 - 2 * p__), 2);

	long double S1 = 2 * g.num_edges * theta____ ; 
	S1 += g.num_nodes() * (g.num_nodes() -1) * theta____ * gamma____ ; 
	S1 /= pow( esp2__ ,2); 
    for (int u = 0; u < g.num_nodes(); ++u) {
        long double du = g.degree[u] * 1.0;
        // Convert vector to unordered_set for faster operations
        std::unordered_set<int> neighbors_u(g.neighbor[u].begin(), g.neighbor[u].end());
        // Calculate |N(u) \ v| for each v
        for (int v = u + 1; v < g.num_nodes(); ++v) {
            // Calculate |N(u) \ {v}|, which is the size of N(u) minus 1 if v is a neighbor
            long double sizeN_u_minus_v = neighbors_u.size() - (neighbors_u.count(v) > 0 ? 1 : 0);

            S1 += sizeN_u_minus_v * pow(gamma____, 2);

            if (neighbors_u.count(v) > 0) { // Check if u is connected to v
                S1 += sizeN_u_minus_v * gamma____;
            }
        }
        S1 += gamma____ * du * (du - 1) / 2;
    }
	S1 /= 9; 
	// cout<<"Lower bound of variance = "<< S1 <<endl; 
	*/