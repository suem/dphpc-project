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

void testGraphIO() {
    std::string inFile = "../test/small_graph.txt";
    std::string outFile = "../test/out.txt";

    Graph g = GraphHelper::readGraphFromFile(inFile);
    GraphHelper::writeGraphToFile(outFile, g);
}

void testGraphGeneration() {
    Graph g = GraphHelper::generateRandomGraph(50, 1);
    for (EdgeIterator e = boost::edges(g).first; e != boost::edges(g).second; e++)
        std::cout << source(*e, g) << " " << target(*e, g) << std::endl;
    GraphHelper::writeGraphToFile("../test/out1.txt", g);
}

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
	if (argc < 2) {
		cerr << "invalid input, no file to read from" << endl;
		return -1;
	}

	try {

		Graph g = GraphHelper::readGraphFromFile(argv[1]);
		verify_bipartite(g);

		vertex_size_t n = num_vertices(g);

		std::cout << "pothen fan" << std::endl;
		for (int i = 0; i < NO_RUNS; ++i) {

			VertexVector mates(n);

			Timer t = Timer();
			pothen_fan(g, mates);
			double elapsed = t.elapsed();

			verify_matching(g, mates);
			volatile vertex_size_t matchingSize = boost::matching_size(g, &mates[0]);

			cout << matchingSize << "\t" <<  elapsed << endl;
		}


		std::cout << "parallel pothen fan" << std::endl;
		for (int i = 0; i < NO_RUNS; ++i) {

			VertexVector mates(n);

			Timer t = Timer();
			parallel_pothen_fan(g, mates);
			double elapsed = t.elapsed();

			verify_matching(g, mates);
			volatile vertex_size_t matchingSize = boost::matching_size(g, &mates[0]);

			cout << matchingSize << "\t" <<  elapsed << endl;
		}


		std::cout << "boost edmonds" << std::endl;
		for (int i = 0; i < NO_RUNS; ++i) {

			VertexVector mates(n);

			Timer t = Timer();
			boost::edmonds_maximum_cardinality_matching(g, &mates[0]);
			double elapsed = t.elapsed();

			verify_matching(g, mates);
			volatile vertex_size_t matchingSize = boost::matching_size(g, &mates[0]);

			cout << matchingSize << "\t" <<  elapsed << endl;
		}





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
