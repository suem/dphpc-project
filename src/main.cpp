#include <iostream>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <string>

using namespace boost;

int main(int argc, char* argv[]) {

	///////////////////////
	// GRAPH HELLO WORLD //
	///////////////////////

	// create a typedef for the Graph type
	typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;

	// Make convenient labels for the vertices
	enum { A, B, C, D, E, N };
	const int num_vertices = N;
	const char* name = "ABCDE";

	// writing out the edges in the graph
	typedef std::pair<int, int> Edge;
	Edge edge_array[] =
	{ Edge(A,B), Edge(A,D), Edge(C,A), Edge(D,C),
		Edge(C,E), Edge(B,D), Edge(D,E) };
	const int num_edges = sizeof(edge_array) / sizeof(edge_array[0]);

	// declare a graph object
	Graph g(num_vertices);

	// add the edges to the graph object
	for (int i = 0; i < num_edges; ++i)
		add_edge(edge_array[i].first, edge_array[i].second, g);

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

	std::cout << "Press a key to exit...";
	while (std::cin.get() != '\n') {} 

	return 0;
}