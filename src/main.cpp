#include <iostream>
#include <string>
#include "graphtypes.h"
#include <boost/graph/bipartite.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

#include "verifier.h"
#include "GraphHelper.h"
#include "pothen_fan.h"
#include "Timer.h" 

using namespace boost;
using namespace std; 

static const int NO_RUNS = 5;
//static const int NO_RUNS = 1;

void testGraphIO() {
    std::string inFile = "../test/graphs/small_graph_bi.txt";
    std::string outFile = "../test/out.txt";

    Graph g;
	Vertex first_right;
    GraphHelper::readGraphFromFile(g, first_right, inFile);
    GraphHelper::writeGraphToFile(outFile, g);
}

void testGraphGeneration() {
    Graph g = GraphHelper::generateRandomGraph(50, 1);
    for (EdgeIterator e = boost::edges(g).first; e != boost::edges(g).second; e++)
        std::cout << source(*e, g) << " " << target(*e, g) << std::endl;
    GraphHelper::writeGraphToFile("../test/out1.txt", g);
}

void runParallelPothenFan(const Graph& g, Vertex first_right, int n, vertex_size_t matching_size_solution, int numThreads) {
	std::cout << "parallel pothen fan with " << numThreads << std::endl;
	for (int i = 0; i < NO_RUNS; ++i) {

		VertexVector mates(n);

		Timer t = Timer();
		parallel_pothen_fan(g, first_right, mates, numThreads);
		double elapsed = t.elapsed();

		verify_matching(g, mates, matching_size_solution);
		volatile vertex_size_t matchingSize = boost::matching_size(g, &mates[0]);

		cout << matchingSize << "\t" <<  elapsed << endl;
	}
}

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
	if (argc < 2) {
		cerr << "invalid input, no file to read from" << endl;
		return -1;
	}

	try {

        Vertex first_right;
		Graph g;
		GraphHelper::readGraphFromFile(g, first_right, argv[1]);
		vertex_size_t n = num_vertices(g);

		verify_bipartite(g);

		VertexVector solution_mates(n);
		boost::edmonds_maximum_cardinality_matching(g, &solution_mates[0]);
		vertex_size_t matching_size_solution = boost::matching_size(g, &solution_mates[0]);


		std::cout << "pothen fan" << std::endl;
		for (int i = 0; i < NO_RUNS; ++i) {

			VertexVector mates(n);

			// run pf and measure time -------------
			Timer t = Timer();
			pothen_fan(g, first_right, mates);
			double elapsed = t.elapsed();
			// -------------------------------------

			verify_matching(g, mates, matching_size_solution);

			volatile vertex_size_t matchingSize = boost::matching_size(g, &mates[0]);

			cout << matchingSize << "\t" <<  elapsed << endl;
		}

        runParallelPothenFan(g, first_right, n, matching_size_solution, 10);
//		for (int i = 10; i < 251; i = i+30) {
//			runParallelPothenFan(g, first_right, n, matching_size_solution, i);
//		}

//		std::cout << "boost edmonds" << std::endl;
//		for (int i = 0; i < NO_RUNS; ++i) {
//
//			VertexVector mates(n);
//
//			Timer t = Timer();
//			boost::edmonds_maximum_cardinality_matching(g, &mates[0]);
//			double elapsed = t.elapsed();
//
//			volatile vertex_size_t matchingSize = boost::matching_size(g, &mates[0]);
//
//			cout << matchingSize << "\t" <<  elapsed << endl;
//		}


//		cout << "Max Matching has cardinality: " << matchingSize << endl;
//		cout << "Matchings: " << endl;
//		VertexIterator start, end;
//		for (tie(start, end) = vertices(g); start != end; start++) {
//			Vertex u = *start;
//			Vertex v = mates[u];
//			if (v != g.null_vertex() && u < v) cout << "(" << u << " " << v << ")" << endl;
//		}

	} catch (char const* error) {
		cerr << "Error: " << error << endl;
		return -1;
	}

	return 0;
}
