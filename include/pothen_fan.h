//
// Created by suem on 10/27/16.
//

#pragma once

#include <atomic>

#include "graphtypes.h"

struct FindPathElement {
	Vertex x0;
	Vertex y;
	AdjVertexIterator start;
	AdjVertexIterator end;
	bool found = false;
};

struct PathElement {
	Vertex x0;
	AdjVertexIterator yiter;
	AdjVertexIterator yiter_end;
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

inline bool claim_vertex(
        Vertex y, Vertex first_right,
		std::atomic<unsigned int>* visited, unsigned int iteration) {
	return std::atomic_exchange(&visited[y - first_right], iteration) != iteration;
}


inline bool is_right(const Vertex& v, const Vertex first_right) {
	return v >= first_right;
}

inline bool is_left(const Vertex& v, const Vertex first_right) {
	return v < first_right;
}

inline bool is_matched(const Vertex& v, const Graph& g, const VertexVector& mate) {
	return mate[v] != g.null_vertex();
}

inline bool is_unmatched(const Vertex& v, const Graph& g, const VertexVector& mate) {
	return mate[v] == g.null_vertex();
}

