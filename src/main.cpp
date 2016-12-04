#include <iostream>
#include <string>
#include "graphtypes.h"
#include <boost/graph/bipartite.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

#include "verifier.h"
#include "GraphHelper.h"
#include "pothen_fan.h"
#include "karpsipser.h"
#include "greedy.h"
#include "pf.h"
#include "ppf1.h"
#include "ppf2.h"
#include "ppf3.h"
#include "ppf4.h"
#include "tree_grafting.h"
#include "unsync_pothen_fan.h"
#include "Timer.h" 

using namespace boost;
using namespace std;

//static const int NO_RUNS = 20;
static const int NO_RUNS = 10;

void testGraphGeneration(int numNodes, float density, char* path) {
	Vertex first_right;
	Graph g = GraphHelper::generateRandomGraph(numNodes, density, first_right);
	GraphHelper::writeGraphToFile(path, g, first_right);
}

void runParallelPothenFan(
		const std::string& graphName, const Graph& g, Vertex first_right, vertex_size_t matching_size_solution, const VertexVector& initialMatching,
		size_t noRuns,
		int numThreads) {

	cout << "#Running ppf " << noRuns << " times with " << numThreads << " threads:" << endl;

	std::vector<double> durations(noRuns);
    

	for (int i = 0; i < noRuns; ++i) {
		VertexVector mates = initialMatching;

		Timer t = Timer();
//		ppf1(g, first_right, mates, numThreads);
//		ppf2(g, first_right, mates, numThreads);
		ppf3(g, first_right, mates, numThreads); // best version so far

        /*
        // ppf4 ---
        vector<MateVisited> matchingVisited(initialMatching.size());
        for (Vertex v = 0; v < num_vertices(g); v++) {
            matchingVisited[v].iteration = 0;
            matchingVisited[v].mate = initialMatching[v];
        }
		ppf4(g, first_right, matchingVisited, numThreads);
        // --------
         */

		double elapsed = t.elapsed();
		durations[i] = elapsed;

        /*
        // ppf4 ------
        for (Vertex v = 0; v < num_vertices(g); v++) mates[v] = matchingVisited[v].mate;
        // -----------
         */

		verify_matching(g, mates, matching_size_solution);
	}

	double average_runtime = 0.0;
	for (int i = 1; i < noRuns; ++i) {
		average_runtime += durations[i];
	}
	average_runtime = average_runtime / (noRuns - 1);
	cout << "#ppf avg: " << average_runtime << endl;
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
	//result.numThreads = 1;
	result.timeStamp = std::string(buff);
	for (int i = 0; i < NO_RUNS; ++i) {

		VertexVector mates = initialMatching;

		Timer t = Timer();
		boost::edmonds_maximum_cardinality_matching(g, &mates[0]);
		double elapsed = t.elapsed();

		//result.durations.push_back(elapsed);
	}

	//GraphHelper::printOutput(result);
}

void runPothenFan(
		const std::string& graphName,
		const Graph& g, Vertex first_right, vertex_size_t matching_size_solution, const VertexVector& initialMatching,
		size_t noRuns
) {

	cout << "#Running pf " << noRuns << " times" << endl;

	std::vector<double> durations(noRuns);

	for (int i = 0; i < noRuns; ++i) {

		VertexVector mates = initialMatching;

		Timer t = Timer();
		pf(g, first_right, mates);
        durations[i] = t.elapsed();

        if (i == 0) verify_matching(g, mates, matching_size_solution);
	}

	double average_runtime = 0.0;
	for (int i = 1; i < noRuns; ++i) {
		average_runtime += durations[i];
	}
	average_runtime = average_runtime / (noRuns - 1);
	cout << "#pf avg: " << average_runtime << endl;
}

void runTreeGrafting(const std::string& graphName, const Graph& g, Vertex first_right, vertex_size_t n, vertex_size_t matching_size_solution, const VertexVector& initialMatching, int numThreads) {
	char buff[20];
	time_t now = time(NULL);
	strftime(buff, 20, "%Y-%m-%d %H:%M:%S", localtime(&now));
	BenchmarkResult result;
	result.algorithm = "tree grafting";
	result.graphName = graphName;
	result.numEdges = num_edges(g);
	result.numVertices = num_vertices(g);
	//result.numThreads = 1;
	result.timeStamp = std::string(buff);
	std::vector<double> times(NO_RUNS);
	for (int i = 0; i < NO_RUNS; ++i) {

		VertexVector mates = initialMatching;

		Timer t = Timer();
		TreeGrafting(g, first_right, mates, numThreads);

		double elapsed = t.elapsed();

		verify_matching(g, mates, matching_size_solution);

		times[i] = elapsed;
	}

	double avg = 0;
	for (auto i : times) {
		avg += i;
	}
	avg /= NO_RUNS;

	std::cout << avg << std::endl;
}

void testKarpSipser() {
	Vertex first_right;
	Graph g = GraphHelper::generateRandomGraph(10000, 0.5, first_right);
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
    /*
	double elapsed;
	Timer t;

	vertex_size_t n = num_vertices(g);

    // this used to be the messed up KS implementation
	std::cout << "Computing simple greedy matching..." << std::endl; t = Timer();
	VertexVector gMatching(n);
    simple_greedy_matching(g, gMatching);
	elapsed = t.elapsed();
	std::cout << "Computed simple greedy matching in: " << elapsed << "s" << std::endl;

	std::cout << "Computing Karp-Sipser..." << std::endl;
	t = Timer();
	VertexVector ksMatching(n);
	karp_sipser(g, first_right, ksMatching);
	elapsed = t.elapsed();
	std::cout << "Computed Karp-Sipser in: " << elapsed << "s" << std::endl;


	std::cout << "Computing Parallel Karp-Sipser..." << std::endl;
	t = Timer();
	VertexVector pksMatching(n);
	parallel_karp_sipser(g, first_right, pksMatching, 40);
	elapsed = t.elapsed();
	std::cout << "Computed Parallel Karp-Sipser in: " << elapsed << "s" << std::endl;


	std::cout << "Computing reference solution (using pf)..." << std::endl;
	t = Timer();
	VertexVector maxMatching = pksMatching;
    pf(g, first_right, maxMatching);
	elapsed = t.elapsed();
	std::cout << "Computed reference solution in: " << elapsed << "s" << std::endl;


    std::cout << "Verifying ks" << std::endl;
	vertex_size_t ksMatchingSize = boost::matching_size(g, &ksMatching[0]);
	verify_matching(g, ksMatching, ksMatchingSize);

    std::cout << "Verifying pks" << std::endl;
	vertex_size_t pksMatchingSize = boost::matching_size(g, &pksMatching[0]);
	verify_matching(g, pksMatching, pksMatchingSize);

    std::cout << "Verifying greedy" << std::endl;
	vertex_size_t gSize = boost::matching_size(g, &gMatching[0]);
	verify_matching(g, gMatching, gSize);

	vertex_size_t maxSize = boost::matching_size(g, &maxMatching[0]);

	std::cout << "simple-greedy: " << gSize << std::endl;
	std::cout << "ks: " << ksMatchingSize << std::endl;
	std::cout << "pks: " << pksMatchingSize << std::endl;

    std::cout << "Simple Greedy matching:\t" << (float)gSize / maxSize * 100 << "%" << std::endl;
	std::cout << "Karp Sipser :\t\t" << (float)ksMatchingSize / maxSize * 100 << "%" << std::endl;
	std::cout << "Parallel Karp Sipser:\t\t" << (float)pksMatchingSize / maxSize * 100 << "%" << std::endl;
     */
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
		//g = GraphHelper::generateRandomGraph(10000, 0.0001f, first_right);
		GraphHelper::readGraphFromFile(g, first_right, argv[1]);
		vertex_size_t n = num_vertices(g);
		vertex_size_t e = num_edges(g);
		cout << "vertices: " << n << " edges: " << e << endl;
	
		cout << "#Verifying if bipartite" << endl;
		verify_bipartite(g);

		cout << "#Run karp sipser to get initial matching" << endl;
        Timer t = Timer();
		//VertexVector initialMatching = GraphHelper::karpSipser(g);
		//VertexVector initialMatching = GraphHelper::greedyMatching(g);
		//VertexVector initialMatching = GraphHelper::ks(g);
        VertexVector initialMatchingKS(n);
        karp_sipser(g, first_right, initialMatchingKS);

        double el = t.elapsed();
		vertex_size_t matching_size_ks = boost::matching_size(g, &initialMatchingKS[0]);
		cout << "#karp sipser matching took: " << el << ", size = " << matching_size_ks << endl;


		cout << "#Run greedy to get initial matching" << endl;
		t = Timer();
		//VertexVector initialMatching = GraphHelper::karpSipser(g);
		//VertexVector initialMatching = GraphHelper::greedyMatching(g);
		//VertexVector initialMatching = GraphHelper::ks(g);
		VertexVector initialMatchingGreedy(n);
		simple_greedy_matching(g, initialMatchingGreedy);
		el = t.elapsed();

		vertex_size_t matching_size_greedy = boost::matching_size(g, &initialMatchingGreedy[0]);
		cout << "#greedy matching took: " << el << ", size = " << matching_size_greedy << endl;

		VertexVector& initialMatching = matching_size_greedy < matching_size_ks ? initialMatchingKS : initialMatchingGreedy;

		cout << "#Computing Solution with pf" << endl;
		VertexVector solution_mates = matching_size_greedy < matching_size_ks ? initialMatchingKS : initialMatchingGreedy;

		pf(g, first_right, solution_mates);
		vertex_size_t matching_size_solution = boost::matching_size(g, &solution_mates[0]);

		cout << "#Karp Sipser Initial Matching: " << (float) matching_size_ks / (float) matching_size_solution << "% optimal" << endl;
		cout << "#Greedy Initial Matching: " << (float) matching_size_greedy / (float) matching_size_solution << "% optimal" << endl;
		
		/*
		cout << "run tree grafting" << std::endl;
		for (int i = 1; i < 9; ++i) {
			cout << i << " threads" << std::endl;
			runTreeGrafting(GraphHelper::getGraphNameFromPath(argv[1]), g, first_right, n, matching_size_solution, initialMatching, i);
		}
		return 0;
		*/
		

		runPothenFan(argv[1], g, first_right, matching_size_solution, initialMatching, 2);

		for (int i = 10; i < 251; i = i + 20) runParallelPothenFan(argv[1], g, first_right, matching_size_solution, initialMatching, 10, i);



		/*
		char buff[20];
		time_t now = time(NULL);
		strftime(buff, 20, "%Y-%m-%d_%H%M%S", localtime(&now));

		BenchmarkResult result;
		result.algorithm = "parallel pothen fan";
		result.graphName = GraphHelper::getGraphNameFromPath(argv[1]);
		result.numEdges = num_edges(g);
		result.numVertices = num_vertices(g);
		result.timeStamp = std::string(buff);
		result.iter = NO_RUNS;

		cout << "#Run ppf" << endl;
		std::vector<int> numThreads;
		std::vector < std::vector<double>> durations;
		for (int i = 10; i < 251; i = i + 20) {
			numThreads.push_back(i);
			durations.push_back(runParallelPothenFan(argv[1], g, first_right, n, matching_size_solution, initialMatching, i));
		}

		result.numThreads = numThreads;
		result.durations = durations;

		GraphHelper::printOutput(result, "C:\\dphpc-test\\");

		 */

//		cout << "#Run unsync_ppf" << endl;
//		for (int i = 1; i < 251; i = i + 20) runUnsyncParallelPothenFan(argv[1], g, first_right, n, /*matching_size_solution,*/ initialMatching, i);
	}
	catch (const char* error) {
		cerr << "Error: " << error << endl;
		return -1;
	}

	return 0;
}
