#include <iostream>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/bipartite.hpp>
#include <string>

// TODO: is this the way this is done in C++?
#include "../include/graphloader.h"

using namespace boost;
using namespace std;

int main(int argc, char* argv[]) {

	cout << "Reading graph from stdin: " << endl;

	Graph g;
	loadGraph(g);

    int n = num_vertices(g);
	cout << "Read graph of size: " << n << endl;

    bool bipartite = is_bipartite(g);

	cout << "Graph is bipartite: " << bipartite << endl;

	////////////////////////
	// OPENMP HELLO WORLD //
	////////////////////////

	// This statement should only print once
	printf("Starting Program!\n");


#pragma omp parallel
	{
		// This statement will run on each thread.
		// If there are 4 threads, this will execute 4 times in total
		printf("Running on multiple threads\n");
	}

	// We're out of the parallelized secion.
	// Therefor, this should execute only once
	printf("Finished!\n");

	std::cout << "Press enter to exit...";
	while (std::cin.get() != '\n') {}

	return 0;
}