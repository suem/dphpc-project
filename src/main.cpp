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
#include "unsync_pothen_fan.h"
#include "Timer.h" 

using namespace boost;
using namespace std;

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
	InitialMatching::karp_sipser(g, first_right, initialMatchingKS);
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
	int nOfMaxThreads = 251;

	char buff[20];
	time_t now;
	BenchmarkResult result;
	double elapsed;

	int actualIter;

	std::vector<int> numThreads;
	std::vector < std::vector<double>> durations;
	std::vector<double> durationsPerRun;

	VertexVector mates;

  // --------------------------------------------------------------------------------------------------------------------------
  // PPF 1
  // --------------------------------------------------------------------------------------------------------------------------
	std::cout << "#Run ppf1 with Karp Sipser" << endl;
	now = time(NULL);
	strftime(buff, 20, "%Y-%m-%d_%H%M%S", localtime(&now));

	result.algorithm = "ppf1_KS";
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
			ppf1(g, first_right, mates, nOfThreads);

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


	std::cout << "#Run ppf1 with Greedy" << endl;
	now = time(NULL);
	strftime(buff, 20, "%Y-%m-%d_%H%M%S", localtime(&now));

	result.algorithm = "ppf1_Greedy";
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
			ppf1(g, first_right, mates, nOfThreads);

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

  // --------------------------------------------------------------------------------------------------------------------------
  // PPF 2
  // --------------------------------------------------------------------------------------------------------------------------
	std::cout << "#Run ppf2 with Karp Sipser" << endl;
	now = time(NULL);
	strftime(buff, 20, "%Y-%m-%d_%H%M%S", localtime(&now));

	result.algorithm = "ppf2_KS";
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
			ppf2(g, first_right, mates, nOfThreads);

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

	std::cout << "#Run ppf2 with Greedy" << endl;
	now = time(NULL);
	strftime(buff, 20, "%Y-%m-%d_%H%M%S", localtime(&now));

	result.algorithm = "ppf2_Greedy";
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
			ppf2(g, first_right, mates, nOfThreads);

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
 
 
   // --------------------------------------------------------------------------------------------------------------------------
  // PPF 3
  // --------------------------------------------------------------------------------------------------------------------------
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
 
  // --------------------------------------------------------------------------------------------------------------------------
  // PPF 4
  // --------------------------------------------------------------------------------------------------------------------------
	std::cout << "#Run ppf4 with Karp Sipser" << endl;
	now = time(NULL);
	strftime(buff, 20, "%Y-%m-%d_%H%M%S", localtime(&now));

	result.algorithm = "ppf4_KS";
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
			vector<MateVisited> matchingVisited(initialMatchingKS.size());
      			for (Vertex v = 0; v < num_vertices(g); v++) {
				matchingVisited[v].iteration = 0;
			        matchingVisited[v].mate = initialMatchingKS[v];
			}	
			t = Timer();
			ppf4(g, first_right, matchingVisited, nOfThreads);

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

	std::cout << "#Run ppf4 with Greedy" << endl;
	now = time(NULL);
	strftime(buff, 20, "%Y-%m-%d_%H%M%S", localtime(&now));

	result.algorithm = "ppf4_Greedy";
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
			vector<MateVisited> matchingVisited(initialMatchingGreedy.size());
      			for (Vertex v = 0; v < num_vertices(g); v++) {
				matchingVisited[v].iteration = 0;
				matchingVisited[v].mate = initialMatchingGreedy[v];
			}
			t = Timer();
			ppf4(g, first_right, matchingVisited, nOfThreads);

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

  // --------------------------------------------------------------------------------------------------------------------------
  // PPF 5
  // --------------------------------------------------------------------------------------------------------------------------
  /*
	std::cout << "#Run ppf5 with Karp Sipser" << endl;
	now = time(NULL);
	strftime(buff, 20, "%Y-%m-%d_%H%M%S", localtime(&now));

	result.algorithm = "ppf5_KS";
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
			vector<MateVisited> matchingVisited(initialMatchingKS.size());
			for (Vertex v = 0; v < num_vertices(g); v++) {
				matchingVisited[v].iteration = 0;
				matchingVisited[v].mate = initialMatchingKS[v];
			}
			t = Timer();
			ppf5(g, first_right, mates, nOfThreads);

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

	std::cout << "#Run ppf5 with Greedy" << endl;
	now = time(NULL);
	strftime(buff, 20, "%Y-%m-%d_%H%M%S", localtime(&now));

	result.algorithm = "ppf5_Greedy";
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
			vector<MateVisited> matchingVisited(initialMatchingGreedy.size());
			for (Vertex v = 0; v < num_vertices(g); v++) {
				matchingVisited[v].iteration = 0;
				matchingVisited[v].mate = initialMatchingGreedy[v];
			}
			t = Timer();
			ppf5(g, first_right, mates, nOfThreads);

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

*/
  // --------------------------------------------------------------------------------------------------------------------------
  // PF
  // --------------------------------------------------------------------------------------------------------------------------

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
		
		runBenchmarks(argv[1]);
		return 0;

	} catch (const char* error) {
		cerr << "Error: " << error << endl;
		return -1;
	}

	return 0;
}
