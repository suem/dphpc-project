#include <iostream>
#include <string>
#include "graphtypes.h"
#include <boost/graph/bipartite.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

#include "verifier.h"
#include "GraphHelper.h"
#include "pothen_fan.h"
#include "unsync_pothen_fan.h"
#include "Timer.h" 

using namespace boost;
using namespace std;

static const int NO_RUNS = 5;
//static const int NO_RUNS = 2;

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

void runParallelPothenFan(const std::string& graphName, const Graph& g, Vertex first_right, size_t n, /*vertex_size_t matching_size_solution,*/ const VertexVector& initialMatching, int numThreads) {
	char buff[20];
	time_t now = time(NULL);
	strftime(buff, 20, "%Y-%m-%d %H:%M:%S", localtime(&now));
	BenchmarkResult result;
	result.algorithm = "parallel pothen fan";
	result.graphName = graphName;
	result.numEdges = num_edges(g);
	result.numVertices = num_vertices(g);
	result.numThreads = numThreads;
	result.timeStamp = std::string(buff);
	for (int i = 0; i < NO_RUNS; ++i) {

		VertexVector mates = initialMatching;

		Timer t = Timer();
		parallel_pothen_fan(g, first_right, mates, numThreads);
		double elapsed = t.elapsed();

//		verify_matching(g, mates, matching_size_solution);
	
		result.durations.push_back(elapsed);
	}

	GraphHelper::printOutput(result);
}


void runUnsyncParallelPothenFan(const std::string& graphName, const Graph& g, Vertex first_right, size_t n, /*vertex_size_t matching_size_solution,*/ const VertexVector& initialMatching, int numThreads) {
	char buff[20];
	time_t now = time(NULL);
	strftime(buff, 20, "%Y-%m-%d %H:%M:%S", localtime(&now));
	BenchmarkResult result;
	result.algorithm = "unsync parallel pothen fan";
	result.graphName = graphName;
	result.numEdges = num_edges(g);
	result.numVertices = num_vertices(g);
	result.numThreads = numThreads;
	result.timeStamp = std::string(buff);
	for (int i = 0; i < NO_RUNS; ++i) {

		VertexVector mates = initialMatching;

		Timer t = Timer();
		unsync_parallel_pothen_fan(g, first_right, mates, numThreads);
		double elapsed = t.elapsed();

//		verify_matching(g, mates, matching_size_solution);
	
		result.durations.push_back(elapsed);
	}

	GraphHelper::printOutput(result);
}

void runBoostEdmonds(const std::string& graphName, const Graph& g, const VertexVector& initialMatching) {
	char buff[20];
	time_t now = time(NULL);
	strftime(buff, 20, "%Y-%m-%d %H:%M:%S", localtime(&now));
	BenchmarkResult result;
	result.algorithm = "boost edmonds";
	result.graphName = graphName;
	result.numEdges = num_edges(g);
	result.numVertices = num_vertices(g);
	result.numThreads = 1;
	result.timeStamp = std::string(buff);
	for (int i = 0; i < NO_RUNS; ++i) {

		VertexVector mates = initialMatching;

		Timer t = Timer();
		boost::edmonds_maximum_cardinality_matching(g, &mates[0]);
		double elapsed = t.elapsed();

		result.durations.push_back(elapsed);
	}

	GraphHelper::printOutput(result);
}

void runPothenFan(const std::string& graphName, const Graph& g, Vertex first_right, int n, /*vertex_size_t matching_size_solution,*/ const VertexVector& initialMatching) {
	char buff[20];
	time_t now = time(NULL);
	strftime(buff, 20, "%Y-%m-%d %H:%M:%S", localtime(&now));
	BenchmarkResult result;
	result.algorithm = "pothen fan";
	result.graphName = graphName;
	result.numEdges = num_edges(g);	
	result.numVertices = num_vertices(g);
	result.numThreads = 1;
	result.timeStamp = std::string(buff);
	for (int i = 0; i < NO_RUNS; ++i) {

		VertexVector mates = initialMatching;

		Timer t = Timer();
		pothen_fan(g, first_right, mates);

		double elapsed = t.elapsed();

//		verify_matching(g, mates, matching_size_solution);
		
		result.durations.push_back(elapsed);
	}

	GraphHelper::printOutput(result);
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

void compareInitialMatching(const Graph& g, const Vertex first_right) {
	double elapsed;
	Timer t;

	std::cout << "Computing greedy matching..." << std::endl; t = Timer();
	VertexVector gMatching = GraphHelper::greedyMatching(g);
	elapsed = t.elapsed();
	std::cout << "Computed greedy matching in: " << elapsed << "s" << std::endl;

    /* Too slow for larger graphs
	std::cout << "Computing Karp-Sipser..." << std::endl;
	t = Timer();
	VertexVector ksMatchMatching = GraphHelper::karpSipser(g);
	elapsed = t.elapsed();
	std::cout << "Computed Karp-Sipser in: " << elapsed << "s" << std::endl;
    */

	std::cout << "Computing Karp-Sipser (no break solution)..." << std::endl;
	t = Timer();
	VertexVector ksFastMatching = GraphHelper::ks(g);
	elapsed = t.elapsed();
	std::cout << "Computed Karp-Sipser (no break solution) in: " << elapsed << "s" << std::endl;


	std::cout << "Computing Parallel Karp-Sipser..." << std::endl;
	t = Timer();
	VertexVector pksMatching = GraphHelper::parallelKarpSipser(g, first_right);
	elapsed = t.elapsed();
	std::cout << "Computed Parallel Karp-Sipser in: " << elapsed << "s" << std::endl;


	std::cout << "Computing reference solution (using pothen_fan)..." << std::endl;
	t = Timer();
	VertexVector maxMatching = ksFastMatching;
    pothen_fan(g, first_right, maxMatching);
	elapsed = t.elapsed();
	std::cout << "Computed reference solution in: " << elapsed << "s" << std::endl;

    /*
	vertex_size_t ksMatchingSize = boost::matching_size(g, &ksMatchMatching[0]);
	verify_matching(g, ksMatchMatching, ksMatchingSize);
    */

	vertex_size_t ksFastSize = boost::matching_size(g, &ksFastMatching[0]);
	verify_matching(g, ksFastMatching, ksFastSize);

	vertex_size_t pksMatchingSize = boost::matching_size(g, &pksMatching[0]);
	verify_matching(g, pksMatching, pksMatchingSize);

	vertex_size_t gSize = boost::matching_size(g, &gMatching[0]);
	verify_matching(g, gMatching, gSize);

	vertex_size_t maxSize = boost::matching_size(g, &maxMatching[0]);

	std::cout << "greedy: " << gSize << std::endl;
	//std::cout << "ks (matching): " << ksMatchingSize << std::endl;
	std::cout << "ks (fast): " << ksFastSize << std::endl;
	std::cout << "pks: " << pksMatchingSize << std::endl;

    std::cout << "Greedy matching:\t" << (float)gSize / maxSize * 100 << "%" << std::endl;
	//std::cout << "Karp Sipser (matching):\t\t" << (float)ksMatchingSize / maxSize * 100 << "%" << std::endl;
	std::cout << "Karp Sipser (fast):\t\t" << (float)ksFastSize / maxSize * 100 << "%" << std::endl;
	std::cout << "Parallel Karp Sipser:\t\t" << (float)pksMatchingSize / maxSize * 100 << "%" << std::endl;
}

void printMatchings(size_t matchingSize, const VertexVector& mates, const Graph& g) {
	cout << "Max Matching has cardinality: " << matchingSize << endl;
	cout << "Matchings: " << endl;
	VertexIterator start, end;
	for (tie(start, end) = vertices(g); start != end; start++) {
		Vertex u = *start;
		Vertex v = mates[u];
		if (v != g.null_vertex() && u < v) cout << "(" << u << " " << v << ")" << endl;
	}
}

int main(int argc, char* argv[]) {
	std::ios::sync_with_stdio(false);

	if (argc < 2) {
		cerr << "invalid input, no file to read from" << endl;
		return -1;
	}

	try {

		cout << "#Reading Graph" << endl;

		Vertex first_right;
		Graph g;
		GraphHelper::readGraphFromFile(g, first_right, argv[1]);
		vertex_size_t n = num_vertices(g);


		cout << "#Verifying if bipartite" << endl;
		verify_bipartite(g);

        /*
		VertexVector solution_mates(n);
		boost::edmonds_maximum_cardinality_matching(g, &solution_mates[0]);
		vertex_size_t matching_size_solution = boost::matching_size(g, &solution_mates[0]);
        */


//		VertexVector initialMatching = GraphHelper::greedyMatching(g);
        // Compute initial matching using karp-sipser
		cout << "#Run KarpSispser to get initial matching" << endl;
        Timer t = Timer();
		//VertexVector initialMatching = GraphHelper::karpSipser(g);
		VertexVector initialMatching = GraphHelper::ks(g);
        double el = t.elapsed();
		cout << "#KarpSispser took: " << el << endl;

		cout << "#Run pf" << endl;
		runPothenFan(argv[1], g, first_right, n,  /*matching_size_solution,*/ initialMatching);

		//cout << "#Run ppf" << endl;
		for (int i = 1; i < 251; i = i + 20) runParallelPothenFan(argv[1], g, first_right, n, /*matching_size_solution,*/ initialMatching, i);

		cout << "#Run unsync_ppf" << endl;
		for (int i = 1; i < 251; i = i + 20) runUnsyncParallelPothenFan(argv[1], g, first_right, n, /*matching_size_solution,*/ initialMatching, i);

//		runBoostEdmonds(argv[1], g, initialMatching);
	}
	catch (char const* error) {
		cerr << "Error: " << error << endl;
		return -1;
	}

	return 0;
}
