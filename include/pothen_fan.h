
#pragma once

#include <atomic>

#include "graphtypes.h"
#include "graphutils.h"

void parallel_pothen_fan(const Graph& g, Vertex first_right, VertexVector& mate, int numThreads);
void pothen_fan(const Graph& g, const Vertex first_right, VertexVector& mate);

////////////////////////////////////////////////////////
// find augmenting path and augmenting them functions //
////////////////////////////////////////////////////////

// atomic
bool find_path_atomic(const Vertex x0, const Graph& g, const Vertex first_right, VertexVector& mate, std::atomic_flag* visited);
bool find_path_recursive_atomic(const Vertex x0, const Graph& g, const Vertex first_right, VertexVector& mate, std::atomic_flag* visited);
bool find_path_la_recursive_atomic(
	const Vertex x0, const Graph& g, const Vertex first_right, VertexVector& mate, std::atomic_flag* visited, std::pair<AdjVertexIterator, AdjVertexIterator>* lookahead);

// non-atomic
bool find_path(const Vertex x0, const Graph& g, const Vertex first_right, VertexVector& mate, bool* visited, std::vector<FindPathElement>& stack);
bool find_path_recursive(const Vertex x0, const Graph& g, const Vertex first_right, VertexVector& mate, bool* visited);
bool find_path_la_recursive(const Vertex x0, const Graph& g, const Vertex first_right, VertexVector& mate, bool* visited, std::pair<AdjVertexIterator, AdjVertexIterator>* lookahead);

/////////////////////////
// lookahead functions //
/////////////////////////

bool lookahead_step(
	const Vertex x0,
	const Graph& g, const Vertex first_right, VertexVector& mate,
	bool* visited,
	std::pair<AdjVertexIterator, AdjVertexIterator>* lookahead);

bool lookahead_step_atomic(
	const Vertex x0,
	const Graph& g, const Vertex first_right, VertexVector& mate,
	//std::atomic_flag* visited,
	std::atomic<unsigned char>* visited,
	std::pair<AdjVertexIterator, AdjVertexIterator>* lookahead);

bool dfs_la(
	const Vertex v,
	const Graph& g, const Vertex first_right,
	VertexVector& mate,
	bool* visited,
	std::pair<AdjVertexIterator, AdjVertexIterator>* lookahead,
	std::vector<PathElement>& stack);

bool dfs_la_atomic(
	const Vertex v,
	const Graph& g, const Vertex first_right, VertexVector& mate,
	//std::atomic_flag* visited,
	std::atomic<unsigned char>* visited,
	std::pair<AdjVertexIterator, AdjVertexIterator>* lookahead,
	std::vector<PathElement>& stack);