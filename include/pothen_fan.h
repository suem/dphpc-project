
#pragma once

#include <atomic>

#include "graphtypes.h"
#include "graphutils.h"

struct FindPathElement {
	Vertex x0;
	Vertex y;
	AdjVertexIterator start;
	AdjVertexIterator end;
	bool found = false;
};

void parallel_pothen_fan(const Graph& g, Vertex first_right, VertexVector& mate, int numThreads);
bool find_path_atomic(const Vertex x0, const Graph& g, const Vertex first_right, VertexVector& mate, std::atomic_flag* visited);
bool find_path_la_atomic(
		const Vertex x0,
		const Graph& g,const Vertex first_right, VertexVector& mate,
//		std::atomic_flag* visited,
		std::atomic<unsigned int>* visited, unsigned int iteration,
		std::pair<AdjVertexIterator, AdjVertexIterator>* lookahead,
		std::vector<PathElement>& stack);

bool find_path_recursive_atomic(const Vertex x0, const Graph& g, const Vertex first_right, VertexVector& mate, std::atomic_flag* visited);
bool find_path_la_recursive_atomic(
		const Vertex x0, const Graph& g, const Vertex first_right, VertexVector& mate, std::atomic_flag* visited, std::pair<AdjVertexIterator, AdjVertexIterator>* lookahead);

void pothen_fan(const Graph& g, const Vertex first_right, VertexVector& mate);
bool find_path(const Vertex x0, const Graph& g, const Vertex first_right, VertexVector& mate, bool* visited, std::vector<FindPathElement>& stack);
bool find_path_la(
        const Vertex x0,
		const Graph& g,const Vertex first_right, VertexVector& mate,
		bool* visited,
		std::pair<AdjVertexIterator, AdjVertexIterator>* lookahead,
		std::vector<PathElement>& stack);


bool find_path_recursive(const Vertex x0, const Graph& g, const Vertex first_right, VertexVector& mate, bool* visited);
bool find_path_la_recursive(const Vertex x0, const Graph& g, const Vertex first_right, VertexVector& mate, bool* visited, std::pair<AdjVertexIterator, AdjVertexIterator>* lookahead);

