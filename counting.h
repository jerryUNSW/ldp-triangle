#pragma once
#ifndef COUNT_H
#define COUNT_H
#include "graph.h"
#include <sstream>


long double solve_equation(long double epsilon, long double S1, long double S2); 

// 
long double weighted_sampling_triangle(Graph& g, unsigned long seed); 

long double hybrid_triangle(Graph& g, unsigned long seed); 

void randomizedResponses(Graph& g, Graph& g2, 
    int query_vertex, int flipping_from, int flipping_to, double p);

// Define function fff(x, alpha)
double fff(double x, double alpha) ; 
// Define objective function to minimize
double objective(const gsl_vector *v, void *params) ;

std::vector<double> minimize(double d1__, double d2__, double epsilon__) ;

long double solve_equation(long double d, long double epsilon) ;
// Function to uniformly sample K vertices from a range
std::vector<int> uniformSampleKVertices(int K, bool useFirstRange, Graph& g); 

long double check_two_hop_path(int s, int v, int t, Graph& g2);

double locally_compute_f_given_q_and_x(int q, int x, Graph& g, Graph& g2) ; 

double compute_f_given_q_x(int q, int x, Graph& g2) ;

vector<double> locally_compute_num_two_paths_given_q(int q, Graph& g, Graph& g2);

vector<double> locally_compute_twoHop_given_q(int q, Graph& g, Graph& g2); 

long double locally_check_distance_two(int s, int t, Graph& g, Graph& g2) ; 


std::vector<std::pair<int, int>> reject_sampling_pairs(int N);


vector<double> compute_twoHop_given_q(int q, Graph& g2); 

long double one_round_even_partition(Graph& g, int k1, int k2, unsigned long seed);

long double one_round_biclique(Graph& g, unsigned long seed, int p__, int q__, int _switch);

std::vector<std::vector<int>> generateCombinations(const std::vector<int>& up_vs, int p) ;

void private_estimate_of_degrees(Graph& g); 

long double BFC_EVP(Graph& g);

bool satisfy_bound(long double upper, long double lower, int u, int x, int v, int w, Graph& g); 

double my_genrand_real2();

void my_init_genrand(unsigned long seed);

// layer based approach:
long double layer_based_two_round(Graph& g, unsigned long seed); 

// u>x && v>w && u >v , each btf is counted once.
long double reduce_redundancy_two_round(Graph& g, unsigned long seed); 

// with edge dependency.
long double reduce_redundancy_two_round_dependency(Graph& g, unsigned long seed);

// new method two round
long double new_method(Graph& g, unsigned long seed); 

// long double BFC_EVP_sample(Graph& g, long double sampling_prob); 

void compute_m3_m2(long double &m4, long double &m3, long double &m2, long double &m1, long double &m0, Graph& g2);

long double one_round_btf(Graph& g, unsigned long seed);

long double wedge_based_two_round_btf(Graph& g, unsigned long seed); 

// naive noisy motif counts
void get_noisy_naive(Graph& g, Graph& g2, long double& local_btfs, long double& local_cate, long double& res ); 
// process butterflies in batch 
void BFC_EVP_noisy(Graph& g, Graph& g2, long double& BTF, long double& cate, long double& res );

void approx_abcore(Graph& g, long double Eps0); 

void construct_noisy_graph(Graph& g, Graph& g2, unsigned long seed);

void construct_noisy_graph_2(int upper, int lower, Graph& g, Graph& g2, unsigned long seed);

void construct_noisy_graph_3(int upper, int lower, Graph& g, Graph& g2, unsigned long seed); 

long double get_wedges(Graph& g);

long double get_cate(Graph& g);

int findPercentile(std::vector<int>& data); 

void test_random_number(unsigned long seed); 

long double get_laplace(long double parameter);

void compute_m3_m2_2(long double &m4, long double &m3, long double &m2, long double &m1, long double &m0, Graph& g2);

// Custom hash function for std::vector<int>
struct VectorHasher {
    std::size_t operator()(const std::vector<int>& v) const {
        std::size_t seed = v.size();
        for (const auto& i : v) {
            seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

unsigned long long int triangle_count(Graph& g); 

long double wedge_based_triangle(Graph& g, unsigned long seed); 

long double efficient_wedge_based_triangle(Graph& g, unsigned long seed) ;


long double compute_w_uv(int u, int v, Graph& g, Graph& g2, bool sampling_noisy_graph); 



void randomized_response_single_bit(int u, int v, Graph& g, Graph& g2) ;

std::vector<std::vector<int>> generateDistinctMotifs(int p, int q, Graph& g, int num_motifs) ;
#endif