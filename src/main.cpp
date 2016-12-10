#include <iostream>
#include <string>
#include "graphtypes.h"
#include <boost/graph/bipartite.hpp>
#include <boost/graph/max_cardinality_matching.hpp>
#include <omp.h>

#include "verifier.h"
#include "GraphHelper.h"
#include "pothen_fan.h"
#include "InitialMatching.h"
#include "greedy.h"
#include "pf.h"
#include "ppf1.h"
#include "ppf2.h"
#include "ppf3.h"
#include "ppf4.h"
#include "tree_grafting.h"
#include "pothen_fan.h"
#include "Timer.h"

using namespace boost;
using namespace std;

typedef void (*run) (const Graph&, Vertex, VertexVector&, int);

void runPTG(const Graph& g, Vertex first_right, VertexVector& mate, int numThreads) {
	TreeGrafting(g, first_right, mate, numThreads);
}

void runBenchmarks(const std::string& graphName) {
	std::vector<std::pair<run, std::string>> functionArray;

	//functionArray.push_back(std::make_pair(&ppf1, "ppf1")); // dont use this, currently incorrect
	functionArray.push_back(std::make_pair(&parallel_pothen_fan, "ppf1"));
	functionArray.push_back(std::make_pair(&ppf2, "ppf2"));
	functionArray.push_back(std::make_pair(&ppf3, "ppf3"));
	functionArray.push_back(std::make_pair(&runPTG, "ptg"));

	bool run_pf = true;

	std::cout << "#Reading Graph" << endl;

	vertex_size_t matching_size_solution;
	Vertex first_right;
	Graph g;
	GraphHelper::readGraphFromFile(g, first_right, matching_size_solution, graphName);

	if (matching_size_solution == 0)
		throw "Invalid graph file! No Matching solution size given";

	vertex_size_t n = num_vertices(g);
	vertex_size_t e = num_edges(g);
	std::cout << "#vertices: " << n << " edges: " << e << endl;

	std::cout << "#Verifying if bipartite" << endl;
	verify_bipartite(g);

	std::cout << "#Run karp sipser to get initial matching" << endl;
	Timer t = Timer();
	VertexVector initialMatchingKS(n);
	InitialMatching::karp_sipser(g, first_right, initialMatchingKS);
	double el = t.elapsed();
	vertex_size_t matching_size_ks = boost::matching_size(g, &initialMatchingKS[0]);
	std::cout << "#karp sipser matching took: " << el << "s, size = " << matching_size_ks << endl;

	std::cout << "#Run greedy to get initial matching" << endl;
	t = Timer();
	VertexVector initialMatchingGreedy(n);
	simple_greedy_matching(g, initialMatchingGreedy);
	el = t.elapsed();
	vertex_size_t matching_size_greedy = boost::matching_size(g, &initialMatchingGreedy[0]);
	std::cout << "#greedy matching took: " << el << "s, size = " << matching_size_greedy << endl;

	std::cout << "#Maximum Matching Size: " << matching_size_solution << " edges" << endl;
	std::cout << "#Karp Sipser Initial Matching: " << 100.0 * (float)matching_size_ks / (float)matching_size_solution << "% optimal" << endl;
	std::cout << "#Greedy Initial Matching: " << 100.0 * (float)matching_size_greedy / (float)matching_size_solution << "% optimal" << endl;

	int nOfRunsDefault = 55;
	int nOfMinThreads = 1;
	int nOfMaxThreads = 251;
	int lengthStride = 1;

	char buff[20];
	time_t now;
	BenchmarkResult result;
	double elapsed;

	int actualIter;

	std::vector<int> numThreads;
	std::vector < std::vector<double>> durations;
	std::vector<double> durationsPerRun;

	VertexVector mates;


	std::vector<std::string> initialMatchingDesc;
	initialMatchingDesc.push_back("KS");
	initialMatchingDesc.push_back("Greedy");
	// loop over given algorithms
	for (auto descPair : functionArray) {
		run functionPointer = descPair.first;
		std::string algoName = descPair.second;

		// loop over initial matching
		for (std::string initialMatchingName : initialMatchingDesc) {

			std::cout << "#Run " << algoName << " with " << initialMatchingName << std::endl;
			now = time(NULL);
			strftime(buff, 20, "%Y-%m-%d_%H%M%S", localtime(&now));

			result.algorithm = algoName + "_" + initialMatchingName;
			result.graphName = GraphHelper::getGraphNameFromPath(graphName);
			result.numEdges = num_edges(g);
			result.numVertices = num_vertices(g);
			result.timeStamp = std::string(buff);
			result.iter = nOfRunsDefault;

			durationsPerRun.resize(result.iter);

			for (int nOfThreads = nOfMinThreads; nOfThreads < nOfMaxThreads; nOfThreads+=lengthStride) {

				//actualIter = nOfThreads * 10 + 1;
				//if (actualIter > result.iter) { actualIter = result.iter; }
				actualIter = result.iter;

				//std::cout << "#Run " << actualIter << " times with " << nOfThreads << " threads" << std::endl;

				std::cout << "# " << nOfThreads;

				for (int run = 0; run < actualIter; run++) {

					VertexVector mates;
					if (initialMatchingName == "KS") {
						mates = initialMatchingKS;
					} else if (initialMatchingName == "Greedy") {
						mates = initialMatchingGreedy;
					} else {
						throw "Invalid initial matching given";
					}

					t = Timer();
					functionPointer(g, first_right, mates, nOfThreads);

					elapsed = t.elapsed();

					verify_matching(g, mates, matching_size_solution);

					durationsPerRun[run] = elapsed;
					std::cout << "," << elapsed << std::flush;
				}
				std::cout << std::endl;

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
		}
	}

  // --------------------------------------------------------------------------------------------------------------------------
  // PF
  // --------------------------------------------------------------------------------------------------------------------------
	if (run_pf) {
		// loop over initial matching
		for (std::string initialMatchingName : initialMatchingDesc) {

			std::cout << "#Run pf with " << initialMatchingName << std::endl;
			now = time(NULL);
			strftime(buff, 20, "%Y-%m-%d_%H%M%S", localtime(&now));

			result.algorithm = "pf_" + initialMatchingName;
			result.graphName = GraphHelper::getGraphNameFromPath(graphName);
			result.numEdges = num_edges(g);
			result.numVertices = num_vertices(g);
			result.timeStamp = std::string(buff);
			result.iter = nOfRunsDefault;

			durationsPerRun.resize(result.iter);

			actualIter = result.iter;

			std::cout << "# 1";

			for (int run = 0; run < actualIter; run++) {
				VertexVector mates;
				if (initialMatchingName == "KS") {
					mates = initialMatchingKS;
				} else if (initialMatchingName == "Greedy") {
					mates = initialMatchingGreedy;
				} else {
					throw "Invalid initial matching given";
				}
				t = Timer();
				pf(g, first_right, mates);

				elapsed = t.elapsed();

				verify_matching(g, mates, matching_size_solution);

				durationsPerRun[run] = elapsed;
				std::cout << "," << elapsed << std::flush;
			}
			std::cout << std::endl;

			numThreads.push_back(1);
			durations.push_back(durationsPerRun);

			result.numThreads = numThreads;
			result.durations = durations;

			GraphHelper::printOutput(result, "");

			numThreads.clear();
			durations.clear();
			durationsPerRun.clear();
			mates.clear();
		}
	}
}

void getMatchingSizeSolution(const std::string& graphName) {
	std::cout << "#Reading Graph" << endl;

	vertex_size_t matching_size_solution_dummy;
	Vertex first_right;
	Graph g;
	GraphHelper::readGraphFromFile(g, first_right, matching_size_solution_dummy, graphName);
	std::cout << "#Current matching solution size read from file: " << matching_size_solution_dummy << " edges" << endl;

	vertex_size_t n = num_vertices(g);
	vertex_size_t e = num_edges(g);
	std::cout << "#vertices: " << n << " edges: " << e << endl;

	std::cout << "#Verifying if bipartite" << endl;
	verify_bipartite(g);

	std::cout << "#Run karp sipser to get initial matching" << endl;
	Timer t = Timer();
	VertexVector initialMatchingKS(n);
	InitialMatching::karp_sipser(g, first_right, initialMatchingKS);
	double el = t.elapsed();
	vertex_size_t matching_size_ks = boost::matching_size(g, &initialMatchingKS[0]);
	std::cout << "#karp sipser matching took: " << el << "s, size = " << matching_size_ks << endl;

	std::cout << "#Run greedy to get initial matching" << endl;
	t = Timer();
	VertexVector initialMatchingGreedy(n);
	simple_greedy_matching(g, initialMatchingGreedy);
	el = t.elapsed();
	vertex_size_t matching_size_greedy = boost::matching_size(g, &initialMatchingGreedy[0]);
	std::cout << "#greedy matching took: " << el << "s, size = " << matching_size_greedy << endl;

	std::cout << "#Computing Solution with pf" << endl;
	VertexVector solution_mates = matching_size_greedy < matching_size_ks ? initialMatchingKS : initialMatchingGreedy;
	pf(g, first_right, solution_mates);
	vertex_size_t matching_size_solution = boost::matching_size(g, &solution_mates[0]);

	std::cout << "#Maximum Matching Size: " << matching_size_solution << " edges" << endl;
	std::cout << "#Karp Sipser Initial Matching: " << 100.0 * (float)matching_size_ks / (float)matching_size_solution << "% optimal" << endl;
	std::cout << "#Greedy Initial Matching: " << 100.0 * (float)matching_size_greedy / (float)matching_size_solution << "% optimal" << endl;

	if (matching_size_solution == matching_size_solution_dummy) {
		std::cout << "#Matching solution size already part of the graph input file!" << endl;
	}
	else {
		std::cout << "#Enter matching solution size into graph file: " << matching_size_solution << endl;
	}
}

int main(int argc, char* argv[]) {
	std::ios::sync_with_stdio(false);

	if (argc < 2) {
		cerr << "invalid input, no file to read from" << endl;
		return -1;
	}

	try {
		runBenchmarks(argv[1]);
		return 0;

		getMatchingSizeSolution(argv[1]);
		return 0;

	} catch (const char* error) {
		cerr << "Error: " << error << endl;
		return -1;
	}

	return 0;
}
