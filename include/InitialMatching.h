#pragma once

#include <atomic>

#include "graphtypes.h"
#include "graphutils.h"

// static class to summarize initial matching methods
class InitialMatching {
public:
	// classical karp sipser
	static void karp_sipser(const Graph& g, const Vertex first_right, VertexVector& matching);
	// a faster karp sipser using another matching method
	static void karp_sipser_fast(const Graph& g, VertexVector& matching);
	// a very fast parallel karp sipser
	static void parallel_karp_sipser(const Graph& g, const Vertex first_right, VertexVector& matching, const int numThreads);
	// kind of karp sipser very fast but no break
	static void ks(const Graph& g, VertexVector& matching);
	// classical greedy matching
	static void greedy_matching(const Graph& g, VertexVector& matching);

private:
	// helper method for karp sipser
	static void ks_matchAndUpdate(
		const Vertex u,
		const Graph& g,
		VertexVector& matching,
		std::vector<size_t>& deg,
		std::vector<bool>& visited);

	// helper method for parallel karp sipser
	static void pks_matchAndUpdate(
		Vertex u,
		const Graph& g,
		VertexVector& matching,
		std::atomic_int* deg,
		std::atomic_flag* visited);

	// helper method for karp sipser
	static void markAdjacentEdges(Vertex v, const Graph& g, std::vector<size_t>& degree);
};