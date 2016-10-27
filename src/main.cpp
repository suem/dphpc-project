#include <iostream>
#include <string>
#include "graphtypes.h"
#include <boost/graph/bipartite.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

#include "graphloader.h"
#include "verifier.h"

using namespace boost;
using namespace std;

int main(int argc, char* argv[]) {

	cout << "Reading graph from stdin: " << endl;
	try {

		Graph g;
		loadGraph(g);
		verify_bipartite(g);

		vertex_size_t n = num_vertices(g);
		MateMap mates(n);

		boost::edmonds_maximum_cardinality_matching(g, &mates[0]);
		verify_matching(g, mates);

		vertex_size_t matchingSize = boost::matching_size(g, &mates[0]);

		cout << "Max Matching has cardinality: " << matchingSize << endl;
		cout << "Matchings: " << endl;
		VertexIterator start, end;
		for (tie(start, end) = vertices(g); start != end; start++) {
			Vertex u = *start;
			Vertex v = mates[u];
			if (v != g.null_vertex() && u < v) cout << "(" << u << " " << v << ")" << endl;
		}

	} catch (char const* error) {
		cerr << "Error: " << error << endl;
		return -1;
	}

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
