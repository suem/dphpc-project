#include <iostream>
#include <string>
#include "graphtypes.h"
#include <boost/graph/bipartite.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

#include "graphloader.h"
#include "verifier.h"
#include "GraphHelper.h"
#include "pothen_fan.h"


using namespace boost;
using namespace std; 

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
        testGraphIO();
        //testGraphGeneration();    

    std::ios::sync_with_stdio(false);

	cout << "Reading graph from stdin: " << endl;
	try {

		Graph g;
		loadGraph(g);
		verify_bipartite(g);

		vertex_size_t n = num_vertices(g);
		VertexVector mates(n);

//		boost::edmonds_maximum_cardinality_matching(g, &mates[0]);
		parallel_pothen_fan(g, mates);

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

	return 0;
}
