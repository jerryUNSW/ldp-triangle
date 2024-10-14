#include "counting.h"

#define ONE_MILLION 1000000

using namespace std;
using namespace std::chrono; 

int num_rounds;

unsigned long long int real_val;

vector<long double> estis, relative_errors;

// vector<long double> cdp_estis, naive_estis, bs_estis, adv_estis;
// vector<long double> bs_estis2;

extern bool one_round, count_cate, naive_switch ; 

long double p; // this is the flipping probability, corresponding to Eps1.

double gamma__ ;

long double Eps, Eps0, Eps1, Eps2;  

vector<int> high_degree_vertices;

// private degrees: 
vector<double> priv_deg; 
int priv_dmax_1, priv_dmax_2; 

int iteration;

extern int alpha;

extern stats::rand_engine_t engine;  // Declare the engine as extern

// efficiency evaluation
long double RR_time, server_side_time, naive_server_side;

long double communication_cost=0;

int q_vertex = -1, x_vertex = -1; 

double d1, d2, epsilon; // these are global variables

int K2 = 5;

bool scalability = false; 


extern bool sampling_vertex_pairs;

extern double vertex_ratio;

int algo_switch ; 

int deg_threshold; 

extern bool edge_centric; 

vector<int> tri_cnt;

vector<vector<double>> tri_cnt_estis;
// just look at the variance of these estimates 
// make it a 2D array. 

bool budget_optimation = true; 

vector<unordered_map<int, bool>> sup_vector; 

int main(int argc, char *argv[]) {

	// input parser:
	Eps =stold(argv[1]);

	string dataset = argv[2]; 

	// the first paramer: 
	num_rounds = atoi(argv[3]); 

	// by default this should be 0 
	// int num_threads_param = atoi(argv[4]);// this is not even used

	vertex_ratio = stod(argv[4]);

	algo_switch = atoi(argv[5]);
	//

	// we use algo switch as the deg_threshold 

	bool eval_time = false, eval_com = false;

	// initialize time
	RR_time=0, server_side_time=0, naive_server_side=0;

	std::mt19937 rng(std::random_device{}()); // for seeding

	bool is_bipartite  = false ; 
	BiGraph g(dataset, is_bipartite);

    unsigned long seed__;
    long double esti_val_;
    long double sum_esti_btf___ = 0.0;
    const int num_runs = num_rounds;
    
    // Calculate the real_val BTF value
    // long double real_val = BFC_EVP(g);
	// omp_set_num_threads(num_threads_param);


	cout<<"vertex_ratio = "<<vertex_ratio <<endl;


	cout<<"n = "<<g.num_nodes()<<endl;

	std::unordered_map<std::string, long int> dataset_map = {
		// {"../graphs/test3.edges", 85334445},
		{"../graphs/email.edges", 727044},
		{"../graphs/youtube.edges", 3056386},
		{"../graphs/wiki.edges", 9203519},
		{"../graphs/gow.edges", 2273138},
		{"../graphs/dblp.edges", 2224385},
		{"../graphs/gplus.edges", 1073677742},
		{"../graphs/skitter.edges", 28769868},
		{"../graphs/imdb.edges",  3856982376}, 
		{"../graphs/lj.edges", 177820130},
		{"../graphs/orkut.edges", 627584181}
	};

	tri_cnt.resize(g.num_nodes()); 
	fill(tri_cnt.begin(), tri_cnt.end(),0); 

	auto it = dataset_map.find(dataset);
	if(dataset_map.find(dataset) != dataset_map.end()){
		real_val = it->second;
	}else{
		sup_vector.resize(g.num_nodes());
		cout<<"len = "<<sup_vector.size()<<endl;
		real_val = triangle_count(g);
	}
	printf("Triangle Count = %lld\n", real_val);






	map<int, int> degree_distribution;  // Using map to keep degrees sorted
	unordered_set<int> degree_one_vertices;
	int cnt = 0;

	int high_degree_count = 0; // Counter for vertices with degree >= threshold

	// cout<<"deg_threshold = " << deg_threshold <<endl;

	// Calculate the degree distribution

	double txxx = omp_get_wtime();

	// for ( vertex_ratio = 1.0; vertex_ratio >= 0.001; vertex_ratio /= 10.0) {
		// for each vertex sampling ratio, vary epsilon_1 

		// if not using such an extreme sampling ratio, it cannot be run within reasonable time 
		// if using such sampling, the effectiveness is undermined significantly. 

		// cout<<"vertex ratio = "<<vertex_ratio<<endl;
		if(vertex_ratio==1){
			sampling_vertex_pairs = false; 
			cout<<"no vertex sampling"<<endl;
		}else{
			sampling_vertex_pairs = true; 
			// cout << "total pairs = " << total_vertex_pairs << endl;
			// cout << "need computation = " << total_vertex_pairs * vertex_ratio << endl;
		}

		// fix this. 
		Eps0 = 0.1*Eps; 

		// for find plot.
		for (int kk = 1; kk<=9; kk++) {
			if(kk==9){
			// if(budget_optimation){
				budget_optimation = true;
				cout<<"compute epsilon1 ad hoc"<<endl;
			}else{
				budget_optimation = false; 
				Eps1 = (Eps - Eps0)*0.5; 
				Eps2 = (Eps - Eps0)*0.5; 
				cout<<"basic setting"<<endl;
				cout<<"Eps1 = "<<Eps1<<endl;
				cout<<"Eps2 = "<<Eps2<<endl;
			}
			relative_errors.clear();
			estis.clear();

			for (int i = 0; i < num_runs; ++i) {
				seed__ = rng();
				if (algo_switch==1){
					// PC algorithm
					esti_val_ = efficient_wedge_based_triangle(g, seed__);
				}
				if (algo_switch==2){
					// VC algorithm
					esti_val_ = efficient_wedge_based_triangle(g, seed__);
				}
				// testing the weighted sampling
				// esti_val_ = weighted_sampling_triangle(g, seed__);

				// sum_esti_btf___ += esti_val_;
				estis.push_back(esti_val_); 

				// Calculate relative error for current run
				cout << "estimate = " << esti_val_ << endl;
				long double relative_error = abs(esti_val_ - real_val) / real_val;
				cout<<"relative error = "<<relative_error <<endl;
				relative_errors.push_back(relative_error);
				cout << endl;
			}	
			// report mean relative errors:
			cout << "# mean rel err = " << calculateMean(relative_errors) << endl;

			long double loss = 0.0;
			// long double mean_estis = calculateMean(estis);
			for (size_t i = 0; i < estis.size(); ++i) {
				loss += static_cast<long double>((estis[i] - real_val) * (estis[i] - real_val));
			}
			loss /= estis.size(); 

			printf("# l2 loss = %.6Le\n", loss);
		}


	// }

	double tyyy = omp_get_wtime();
	double seconds = tyyy - txxx;


    // Output results
	// printf("# mean = %Lf\n", calculateMean(estis));

    // cout << "Mean rel err = " << calculateMean(relative_errors) << endl;

    cout << "real_val count = " << real_val << endl;

	printf("time:%f\n", seconds);
	

	
	/*
	// need to revisit this to better predict the L2 loss.
	long double estimated_variance = 0; 
	long double p__ = 1/(1 + exp(Eps1)); 
	long double gamma__ = p__ * (1-p__) / pow(1-2*p__,2); 
	long double lap_scale =  (1-p__) / (1-2*p__); 
	lap_scale /= Eps2; 
	
	long double S1 =0, S2 = 0; 
	for(int i=0;i<g.num_nodes();i++){

		long double deg_u = g.degree[i]; 
		S1 += (deg_u*(deg_u-1)/2);  // deg_u choose 2 
		S2 += deg_u * deg_u;		// deg_u^2
	}
	
	estimated_variance = gamma__ * S1 + 2 * pow(lap_scale,2) * S2; 
	estimated_variance /= 9; 
	*/

	
	
	// what we learned is that high degree vertices will result in higher standard deviations. 
	// how do we handle these vertices? 
		
	// variance is just the expected L2 loss. 
	// printf("# Variance = %Lf\n", calculateVariance(estis) );

	


	// printf("# Estimated Variance = %Lf\n", estimated_variance );

	// cout<<"Lower bound variance = " << S1 <<endl; 

	return 0;
}
