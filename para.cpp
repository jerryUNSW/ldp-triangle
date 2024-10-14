		#pragma omp parallel for reduction(+:sum___) 
		for (size_t i = 0; i < N; ++i) {
		
			for (size_t j = i + 1; j < N; ++j) {

				int u = use_high_deg_vertices ? high_degree_vertices[i] : i;
				int v = use_high_deg_vertices ? high_degree_vertices[j] : j;
				long double local_res = 0;


				// proceed with vertex pair with a fixed probability
				if (sampling_vertex_pairs && genrand_real2() >= vertex_ratio) {
					continue;
				}

				// #pragma omp critical
				// {
					if (!g2.has_computed(u, v)){
						randomized_response_single_bit(u, v, g, g2);
					}
				// }
				local_res = ((g2.edge_vector[min(u, v)][max(u, v)] ? 1 : 0) - p * p____) / (p____ * (1 - 2 * p));

                // Compute w_u_v
				long double w_u_v; 
				// #pragma omp critical
				// {
					w_u_v = compute_w_uv(u, v, g, g2, sampling_noisy_graph);
				// }
                local_res *= w_u_v;


                #pragma omp critical
                sum___ += local_res;
            }
        }