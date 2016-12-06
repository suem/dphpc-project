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
#include "ppf5.h"
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
		//ppf1(g, first_right, mates, numThreads);
//		ppf2(g, first_right, mates, numThreads);
		//ppf3(g, first_right, mates, numThreads); // best version so far
		ppf5(g, first_right, mates, numThreads);

        /*
        // ppf4 --- not that good
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
	std::cout << "Max Matching has cardinality: " << matchingSize << endl;
	std::cout << "Matchings: " << endl;
	VertexIterator start, end;
	for (tie(start, end) = vertices(g); start != end; start++) {
		Vertex u = *start;
		Vertex v = mates[u];
		if (v != g.null_vertex() && u < v) std::cout << "(" << u << " " << v << ")" << endl;
	}
}

void runBenchmarks(const std::string& graphName) {
	std::cout << "#Reading Graph" << endl;

	Vertex first_right;
	Graph g;
	GraphHelper::readGraphFromFile(g, first_right, graphName);
	vertex_size_t n = num_vertices(g);
	vertex_size_t e = num_edges(g);
	std::cout << "#vertices: " << n << " edges: " << e << endl;

	std::cout << "#Verifying if bipartite" << endl;
	verify_bipartite(g);

	std::cout << "#Run karp sipser to get initial matching" << endl;
	Timer t = Timer();
	VertexVector initialMatchingKS(n);
	karp_sipser(g, first_right, initialMatchingKS);
	double el = t.elapsed();
	vertex_size_t matching_size_ks = boost::matching_size(g, &initialMatchingKS[0]);
	std::cout << "#karp sipser matching took: " << el << ", size = " << matching_size_ks << endl;


	std::cout << "#Run greedy to get initial matching" << endl;
	t = Timer();
	VertexVector initialMatchingGreedy(n);
	simple_greedy_matching(g, initialMatchingGreedy);
	el = t.elapsed();
	vertex_size_t matching_size_greedy = boost::matching_size(g, &initialMatchingGreedy[0]);
	std::cout << "#greedy matching took: " << el << ", size = " << matching_size_greedy << endl;

	std::cout << "#Computing Solution with pf" << endl;
	VertexVector solution_mates = matching_size_greedy < matching_size_ks ? initialMatchingKS : initialMatchingGreedy;
	pf(g, first_right, solution_mates);
	vertex_size_t matching_size_solution = boost::matching_size(g, &solution_mates[0]);

	std::cout << "#Karp Sipser Initial Matching: " << (float)matching_size_ks / (float)matching_size_solution << "% optimal" << endl;
	std::cout << "#Greedy Initial Matching: " << (float)matching_size_greedy / (float)matching_size_solution << "% optimal" << endl;

	int nOfRunsDefault = 101;
	int nOfMaxThreads = 31;

	char buff[20];
	time_t now;
	BenchmarkResult result;
	double elapsed;

	int actualIter;

	std::vector<int> numThreads;
	std::vector < std::vector<double>> durations;
	std::vector<double> durationsPerRun;

	VertexVector mates;

	std::cout << "#Run ppf3 with Karp Sipser" << endl;
	now = time(NULL);
	strftime(buff, 20, "%Y-%m-%d_%H%M%S", localtime(&now));

	result.algorithm = "ppf3_KS";
	result.graphName = GraphHelper::getGraphNameFromPath(graphName);
	result.numEdges = num_edges(g);
	result.numVertices = num_vertices(g);
	result.timeStamp = std::string(buff);
	result.iter = nOfRunsDefault;

	durationsPerRun.resize(result.iter);


	for (int nOfThreads = 1; nOfThreads < nOfMaxThreads; nOfThreads++) {

		//actualIter = nOfThreads * 10 + 1;
		//if (actualIter > result.iter) { actualIter = result.iter; }
		actualIter = result.iter;

		std::cout << "#Run " << actualIter << " times with " << nOfThreads << " threads" << std::endl;
			
		for (int run = 0; run < actualIter; run++) {
			mates = initialMatchingKS;
			t = Timer();
			ppf3(g, first_right, mates, nOfThreads);

			elapsed = t.elapsed();
			durationsPerRun[run] = elapsed;
		}

		numThreads.push_back(nOfThreads);
		durations.push_back(durationsPerRun);
	}

	result.numThreads = numThreads;
	result.durations = durations;

	GraphHelper::printOutput(result, "");

	numThreads.clear();
	durations.clear();
	durationsPerRun.clear();
	mates.clear();



	std::cout << "#Run ppf3 with Greedy" << endl;
	now = time(NULL);
	strftime(buff, 20, "%Y-%m-%d_%H%M%S", localtime(&now));

	result.algorithm = "ppf3_Greedy";
	result.graphName = GraphHelper::getGraphNameFromPath(graphName);
	result.numEdges = num_edges(g);
	result.numVertices = num_vertices(g);
	result.timeStamp = std::string(buff);
	result.iter = nOfRunsDefault;

	durationsPerRun.resize(result.iter);

	for (int nOfThreads = 1; nOfThreads < nOfMaxThreads; nOfThreads++) {

		//actualIter = nOfThreads * 10 + 1;
		//if (actualIter > result.iter) { actualIter = result.iter; }
		actualIter = result.iter;

		std::cout << "#Run " << actualIter << " times with " << nOfThreads << " threads" << std::endl;

		for (int run = 0; run < actualIter; run++) {
			mates = initialMatchingGreedy;
			t = Timer();
			ppf3(g, first_right, mates, nOfThreads);

			elapsed = t.elapsed();
			durationsPerRun[run] = elapsed;
		}

		numThreads.push_back(nOfThreads);
		durations.push_back(durationsPerRun);
	}

	result.numThreads = numThreads;
	result.durations = durations;

	GraphHelper::printOutput(result, "");	
	
	numThreads.clear();
	durations.clear();
	durationsPerRun.clear();
	mates.clear();




	std::cout << "#Run pf with KS" << endl;
	now = time(NULL);
	strftime(buff, 20, "%Y-%m-%d_%H%M%S", localtime(&now));

	result.algorithm = "pf_KS";
	result.graphName = GraphHelper::getGraphNameFromPath(graphName);
	result.numEdges = num_edges(g);
	result.numVertices = num_vertices(g);
	result.timeStamp = std::string(buff);
	result.iter = nOfRunsDefault;

	durationsPerRun.resize(result.iter);

	std::cout << "#Run " << result.iter << " times with " << 1 << " threads" << std::endl;
	for (int run = 0; run < result.iter; run++) {
		mates = initialMatchingKS;
		t = Timer();
		pf(g, first_right, mates);

		elapsed = t.elapsed();
		durationsPerRun[run] = elapsed;
	}

	numThreads.push_back(1);
	durations.push_back(durationsPerRun);

	result.numThreads = numThreads;
	result.durations = durations;

	GraphHelper::printOutput(result, "");

	numThreads.clear();
	durations.clear();
	durationsPerRun.clear();
	mates.clear();




	std::cout << "#Run pf with Greedy" << endl;
	now = time(NULL);
	strftime(buff, 20, "%Y-%m-%d_%H%M%S", localtime(&now));

	result.algorithm = "pf_Greedy";
	result.graphName = GraphHelper::getGraphNameFromPath(graphName);
	result.numEdges = num_edges(g);
	result.numVertices = num_vertices(g);
	result.timeStamp = std::string(buff);
	result.iter = nOfRunsDefault;

	durationsPerRun.resize(result.iter);


	std::cout << "#Run " << result.iter << " times with " << 1 << " threads" << std::endl;
	for (int run = 0; run < result.iter; run++) {
		mates = initialMatchingGreedy;
		t = Timer();
		pf(g, first_right, mates);

		elapsed = t.elapsed();
		durationsPerRun[run] = elapsed;
	}

	numThreads.push_back(1);
	durations.push_back(durationsPerRun);

	result.numThreads = numThreads;
	result.durations = durations;

	GraphHelper::printOutput(result, "");

}

int main(int argc, char* argv[]) {
	std::ios::sync_with_stdio(false);

	if (argc < 2) {
		cerr << "invalid input, no file to read from" << endl;
		return -1;
	}

	try {
		
		//runBenchmarks(argv[1]);
		//return 0;

		cout << "#Reading Graph" << endl;

		Vertex first_right;
		Graph g;
		GraphHelper::readGraphFromFile(g, first_right, argv[1]);
		vertex_size_t n = num_vertices(g);
		vertex_size_t e = num_edges(g);
		cout << "vertices: " << n << " edges: " << e << endl;
	
		cout << "#Verifying if bipartite" << endl;
		verify_bipartite(g);

		cout << "#Run karp sipser to get initial matching" << endl;
        Timer t = Timer();
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


		runPothenFan(argv[1], g, first_right, matching_size_solution, initialMatching, 2);

		runParallelPothenFan(argv[1], g, first_right, matching_size_solution, initialMatching, 10, 20);
		//for (int i = 10; i < 251; i = i + 20) runParallelPothenFan(argv[1], g, first_right, matching_size_solution, initialMatching, 10, i);

        /*
		cout << "run tree grafting" << std::endl;
		for (int i = 1; i < 9; ++i) {
			cout << i << " threads" << std::endl;
			runTreeGrafting(GraphHelper::getGraphNameFromPath(argv[1]), g, first_right, n, matching_size_solution, initialMatching, i);
		}
		for (int i = 10; i < 251; i+=10) {
		cout << i << " threads" << std::endl;
		runTreeGrafting(GraphHelper::getGraphNameFromPath(argv[1]), g, first_right, n, matching_size_solution, initialMatching, i);
		}
         */

	} catch (const char* error) {
		cerr << "Error: " << error << endl;
		return -1;
	}

	return 0;
}
