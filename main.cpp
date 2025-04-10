#include <igraph/igraph.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <cstdlib>
#include <ctime>

const int N = 100000;
int T;
int percent;

void initialize_network(igraph_t &graph, std::vector<bool> &nodes) {
    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, N, N * 10, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    nodes.resize(N);
    for (int i = 0; i < N; ++i) {
        nodes[i] = (std::rand() % 100 < percent) ? true : false;
    }
}

double compute_infected_fraction(const std::vector<bool> &nodes) {
    int infected = 0;
    for (const auto &node : nodes) {
        if(node) infected++;
    }
    return (double)infected / N;
}

bool find_edge(const igraph_vector_int_t* edges, int u, int v) {
    long int i;
    for (i = 0; i < igraph_vector_int_size(edges); i += 2) {
        if ((VECTOR(*edges)[i] == u && VECTOR(*edges)[i + 1] == v) ||
            (VECTOR(*edges)[i] == v && VECTOR(*edges)[i + 1] == u)) {
            return true;
        }
    }
    return false;
}

double simulate(igraph_t &graph, std::vector<bool> &nodes, double beta, double gamm, double omega) {
    double avg_f = 0.0;
    for (int t = 0; t < T; t++) {
        igraph_vector_int_t add_edges, del_edges;
        igraph_vector_int_init(&add_edges, 0);
        igraph_vector_int_init(&del_edges, 0);
        for (int i = 0; i < N; ++i) {
            if (!nodes[i]) {
                igraph_vector_int_t neighbors;
                igraph_vector_int_init(&neighbors, 0);
                igraph_neighbors(&graph, &neighbors, i, IGRAPH_ALL);
                for (long j = 0; j < igraph_vector_int_size(&neighbors); j++) {
                    int neighbor = VECTOR(neighbors)[j];
                    if (nodes[neighbor] && (std::rand() / double(RAND_MAX)) < beta) {
                        nodes[i] = true;
                        break;
                    }
                }
                igraph_vector_int_destroy(&neighbors);
            } 
            else {
                if ((std::rand() / double(RAND_MAX)) < gamm) {
                    nodes[i] = false;
                }
            }
        }
        for (int i = 0; i < N; ++i) {
            if (!nodes[i]) {
                igraph_vector_int_t neighbors;
                igraph_vector_int_init(&neighbors, 0);
                igraph_neighbors(&graph, &neighbors, i, IGRAPH_ALL);
                for (long j = 0; j < igraph_vector_int_size(&neighbors); ++j) {
                    int neighbor = VECTOR(neighbors)[j];
                    if (nodes[neighbor] && (std::rand() / double(RAND_MAX)) < omega) {
                        igraph_integer_t new_neighbor;
                        igraph_bool_t isConnected;
                        do {
                            new_neighbor = std::rand() % N;
                            igraph_are_adjacent(&graph, i, new_neighbor, &isConnected);
                        } while (new_neighbor == i || isConnected || nodes[new_neighbor] || find_edge(&add_edges, i, new_neighbor));
                        igraph_vector_int_push_back(&del_edges, i);
                        igraph_vector_int_push_back(&del_edges, neighbor);
                        igraph_vector_int_push_back(&add_edges, i);
                        igraph_vector_int_push_back(&add_edges, new_neighbor);
                        break;
                    }
                }
                igraph_vector_int_destroy(&neighbors);
            }
        }
        if (igraph_vector_int_size(&del_edges) > 0) {
            igraph_es_t es;
            igraph_es_pairs(&es, &del_edges, IGRAPH_UNDIRECTED);
            igraph_delete_edges(&graph, es);
            igraph_es_destroy(&es);
        }
        if (igraph_vector_int_size(&add_edges) > 0) {
            igraph_add_edges(&graph, &add_edges, 0);
        }
        igraph_vector_int_destroy(&add_edges);
        igraph_vector_int_destroy(&del_edges);
    }
    avg_f = compute_infected_fraction(nodes);
    std::cout << "Beta: " << beta << ", Infected rated: " << avg_f << '\n';
    return avg_f;
}
void knn_and_distribution(igraph_t &graph, std::vector<bool> &nodes, double omega, int initial, double beta) {
    std::string d_file = std::to_string(omega) + "_" + std::to_string(initial) + "_" + std::to_string(beta) + "_" + "degree_distribution.txt";
    std::string k_file = std::to_string(omega) + "_" + std::to_string(initial) + "_" + std::to_string(beta) + "_" + "mean_nearest-neighbor_degree.txt";
    std::ofstream outFile0(d_file.c_str());
    std::ofstream outFile1(k_file.c_str());
    igraph_vector_int_t degrees;
    igraph_vector_int_t neighbors;
    igraph_integer_t neighbors_0;
    igraph_vector_int_init(&degrees, 0);
    igraph_vector_int_init(&neighbors, 0);

    std::unordered_map<long, long> s_state, i_state, k, nn;
    long int num_s = 0;
    long int num_i = 0;
    long int max_s = 0;
    long int max_i = 0;
    long int max_k = 0;
    igraph_degree(&graph, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);
    for (long int i = 0; i < N; i++) {
        if (!nodes[i]) {
            s_state[VECTOR(degrees)[i]]++;
            num_s++;
            max_s = std::max(max_s, VECTOR(degrees)[i]);
        }
        else {
            i_state[VECTOR(degrees)[i]]++;
            num_i++;
            max_i = std::max(max_i, VECTOR(degrees)[i]);
        }
        igraph_neighbors(&graph, &neighbors, i, IGRAPH_ALL);
        long int avg = 0;
        long int n = igraph_vector_int_size(&neighbors);
        for (long int j = 0; j < n; ++j) {
            long int neighbor = VECTOR(neighbors)[j];
            igraph_degree_1(&graph, &neighbors_0, neighbor, IGRAPH_ALL, false);
            avg += neighbors_0;
        }
        if(n > 0) avg /= n;
        k[VECTOR(degrees)[i]] += avg;
        nn[VECTOR(degrees)[i]]++;
        max_k = std::max(max_s, max_i);
    }
    for(long int i = 0; i <= max_s; i++) {
        outFile0 << "Susceptibles degree: " << i << "->" << (double)s_state[i] / num_s << "\n";
    }
    for(long int i = 0; i <= max_i; i++) {
        outFile0 << "Infected degree: " << i << "->" << (double)i_state[i] / num_i << "\n";
    }
    for(long int i = 0; i <= max_k; i++) {
        outFile1 << "Knn: " << i << "->" << (nn[i] == 0 ? 0 : (double)k[i] / nn[i]) << "\n";
    }

    igraph_vector_int_destroy(&degrees);
    igraph_vector_int_destroy(&neighbors);
    outFile0.close();
    outFile1.close();
}
void output_results(const std::vector<bool> &nodes) {
    int susceptible = 0, infected = 0;
    for (const auto &node : nodes) {
        if(!node) ++susceptible;
        else ++infected;
    }
    std::cout << "Susceptible: " << susceptible << ", Infected: " << infected << std::endl;
}

int main(int argc, char** argv) {
    if(argc != 7) {
        std::cout << "Hi! Welcome to social network epidemic simulation model!" << std::endl;
        std::cout << "Please enter complete parameter to execute this program." << std::endl;
        std::cout << "./igraph_test [time] [omega] [beta] [gamma] [initial infection percentage] [output file]" << std::endl;
        return 1;
    }
    T = std::stoi(argv[1]);
    double omega = std::stof(argv[2]);
    double beta = std::stof(argv[3]);
    double gamm = std::stof(argv[4]);
    percent = std::stoi(argv[5]);
    std::string result = std::string(argv[6]);

    std::srand(std::time(nullptr));
    std::ofstream outFile2(result.c_str());
    double avg_f = 0.0;
    std::vector<bool> nodes;
    
    for(; beta <= 0.0081; beta += 0.0001) {
        std::cout << "########## Simulation Start ##########" << std::endl;
        igraph_t graph;
        initialize_network(graph, nodes);
        avg_f = simulate(graph, nodes, beta, gamm, omega);
        outFile2 << "Beta: " << beta << ", Infected rated: " << avg_f << '\n';
        output_results(nodes);
        if(getenv("KNN_AND_DEGREE_DISTRIBUTION")) knn_and_distribution(graph, nodes, omega, percent, beta);
        igraph_destroy(&graph);
        std::cout << "########## Simulation End   ##########" << std::endl;
        std::cout << std::endl;
    }
    outFile2.close();
    return 0;
}