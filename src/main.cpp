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
	std::string inFile = "../test/graphs/small_graph_bi.txt";
	std::string outFile = "../test/out.txt";

	Graph g;
	Vertex first_right;
	GraphHelper::readGraphFromFile(g, first_right, inFile);
	GraphHelper::writeGraphToFile(outFile, g);
}

void testGraphGeneration() {
	Graph g = GraphHelper::generateRandomGraph(40, 1);
	for (EdgeIterator e = boost::edges(g).first; e != boost::edges(g).second; e++)
		std::cout << source(*e, g) << " " << target(*e, g) << std::endl;
	GraphHelper::writeGraphToFile("../test/out1.txt", g);
}

void runParallelPothenFan(const Graph& g, Vertex first_right, size_t n, vertex_size_t matching_size_solution, VertexVector& initialMatching, int numThreads) {
	std::cout << "parallel pothen fan with " << numThreads << std::endl;
	for (int i = 0; i < NO_RUNS; ++i) {

		VertexVector mates = initialMatching;

		Timer t = Timer();
		parallel_pothen_fan(g, first_right, mates, numThreads);
		double elapsed = t.elapsed();

		verify_matching(g, mates, matching_size_solution);
		volatile vertex_size_t matchingSize = boost::matching_size(g, &mates[0]);

		cout << matchingSize << "\t" << elapsed << endl;
	}
}


void runPothenFan(const Graph& g, Vertex first_right, int n, vertex_size_t matching_size_solution, VertexVector& initialMatching) {
	std::cout << "pothen fan" << std::endl;
	for (int i = 0; i < NO_RUNS; ++i) {

		VertexVector mates = initialMatching;

		Timer t = Timer();
		pothen_fan(g, first_right, mates);
		double elapsed = t.elapsed();

		verify_matching(g, mates, matching_size_solution);
		volatile vertex_size_t matchingSize = boost::matching_size(g, &mates[0]);

		cout << matchingSize << "\t" << elapsed << endl;
	}
}

void testKarpSipser() {
	Graph g = GraphHelper::generateRandomGraph(5000, 0.5);
	VertexVector matching = GraphHelper::karpSipser(g);
	vertex_size_t matching_size = boost::matching_size(g, &matching[0]);

	VertexVector max_matching(boost::num_vertices(g));
	boost::edmonds_maximum_cardinality_matching(g, &max_matching[0]);
	vertex_size_t max_matching_size = boost::matching_size(g, &max_matching[0]);

	std::cout << "Matching Size:\t\t" << matching_size << std::endl;
	std::cout << "Max Matching Size:\t" << max_matching_size << std::endl;

	if (GraphHelper::isMaximumMatching(matching, g)) {
		std::cout << "Maximum matching found!" << std::endl;
	}
	std::cout << "Finished" << std::endl;
}

void compareInitialMatching(Graph& g) {
	std::cout << "Computing Karp-Sipser..." << std::endl;
	Timer t = Timer();
	VertexVector ksMatching = GraphHelper::karpSipser(g);
	double elapsed = t.elapsed();
	std::cout << "Computed Karp-Sipser in: " << elapsed << "s" << std::endl;
	std::cout << "Computing greedy matching..." << std::endl;
	t = Timer();
	VertexVector gMatching = GraphHelper::greedyMatching(g);
	elapsed = t.elapsed();
	std::cout << "Computed greedy matching in: " << elapsed << "s" << std::endl;
	std::cout << "Computing reference solution..." << std::endl;
	t = Timer();
	VertexVector maxMatching(boost::num_vertices(g));
	boost::edmonds_maximum_cardinality_matching(g, &maxMatching[0]);
	elapsed = t.elapsed();
	std::cout << "Computed reference solution in: " << elapsed << "s" << std::endl;

	vertex_size_t ksSize = boost::matching_size(g, &ksMatching[0]);
	vertex_size_t gSize = boost::matching_size(g, &gMatching[0]);
	vertex_size_t maxSize = boost::matching_size(g, &maxMatching[0]);

	std::cout << "greedy: " << gMatching.size() << std::endl;
	std::cout << "Greedy matching:\t" << (float)gSize / maxSize * 100 << "%" << std::endl;
	std::cout << "Karp Sipser:\t\t" << (float)ksSize / maxSize * 100 << "%" << std::endl;
}

int main(int argc, char* argv[]) {
	std::ios::sync_with_stdio(false);
	//testKarpSipser();
	//return 0;
	
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


		// Compute initial matching using karp-sister
		compareInitialMatching(g);
		return 0;
		//		VertexVector initialMatching = GraphHelper::karpSipser(g);
		VertexVector initialMatching = GraphHelper::greedyMatching(g);


		std::cout << "pothen fan" << std::endl;
		for (int i = 0; i < NO_RUNS; ++i) {

			// copy initial matching
			VertexVector mates = initialMatching;

			// run pf and measure time -------------
			Timer t = Timer();
			pothen_fan(g, first_right, mates);
			double elapsed = t.elapsed();
			// -------------------------------------

			verify_matching(g, mates, matching_size_solution);

			volatile vertex_size_t matchingSize = boost::matching_size(g, &mates[0]);

			cout << matchingSize << "\t" << elapsed << endl;
		}

		for (int i = 10; i < 251; i = i + 30) runParallelPothenFan(g, first_right, n, matching_size_solution, initialMatching, i);
		/*
				// Compute initial matching using karp-sister
				Timer t = Timer();
				VertexVector initialMatching = GraphHelper::karpSipser(g);
		//		VertexVector initialMatching = GraphHelper::greedyMatching(g);
				vertex_size_t matching_size_initial = boost::matching_size(g, &initialMatching[0]);
				double initial_matching_elapsed = t.elapsed();


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

				std::cout << "inital matching: " << matching_size_initial
						  << " of " << matching_size_solution << " (" << (((float) matching_size_initial / (float) matching_size_solution)) * 100 << "%)"
						  << " (" << initial_matching_elapsed << " s)" <<std::endl;

		*/
		runPothenFan(g, first_right, n, matching_size_solution, initialMatching);

		//		cout << "Max Matching has cardinality: " << matchingSize << endl;
		//		cout << "Matchings: " << endl;
		//		VertexIterator start, end;
		//		for (tie(start, end) = vertices(g); start != end; start++) {
		//			Vertex u = *start;
		//			Vertex v = mates[u];
		//			if (v != g.null_vertex() && u < v) cout << "(" << u << " " << v << ")" << endl;
		//		}

	}
	catch (char const* error) {
		cerr << "Error: " << error << endl;
		return -1;
	}

	return 0;
}
