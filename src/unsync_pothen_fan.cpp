//
// Created by suem on 10/27/16.
//
#include <iostream>
#include <stack>
#include <omp.h>
#include <atomic>

#include "graphutils.h"
#include "graphtypes.h"
#include "unsync_pothen_fan.h"

void unsync_parallel_pothen_fan(const Graph& g, Vertex first_right, VertexVector& mate, int numThreads) {

	const int NO_THREADS = numThreads; 

	const vertex_size_t n = boost::num_vertices(g);
	const vertex_size_t n_right = n - first_right;

	const int nt = std::min(static_cast<int>(n), NO_THREADS);

	volatile bool path_found;

//	std::atomic_flag* visited = new std::atomic_flag[n_right];

    std::atomic<unsigned int>* visited = new std::atomic<unsigned int>[n_right];
	memset(visited, 0, sizeof(std::atomic<unsigned int>) * n_right);
    unsigned int iteration = 0;


    // initialize lookahead
    std::pair<AdjVertexIterator, AdjVertexIterator>* lookahead = new std::pair<AdjVertexIterator, AdjVertexIterator>[first_right];
    for (Vertex i = 0; i < first_right; i++) lookahead[i] = boost::adjacent_vertices(i, g);

	do {
		path_found = false;

        // go to next iteration, equivalent to clearing all visited flags
        iteration++;

//		 // here we don't need atomic clears to reset the flags
//        memset(visited, 0, sizeof(std::atomic_flag) * n_right);

		std::vector<PathElement> stack;
#pragma omp parallel num_threads(nt) private(stack)
#pragma omp for
		for (Vertex v = 0; v < first_right; v++) {

			// skip if vertex is already matched
			if (is_matched(v, g, mate))  continue;

			bool path_found_v = dfs_la(v, g, first_right, mate, visited, iteration, lookahead, stack);
			if (path_found_v && !path_found) path_found = true;
		}

	} while (path_found);

	delete[] visited;
}



inline bool lookahead_step_atomic(
		const Vertex x0,
		const Graph& g,const Vertex first_right, VertexVector& mate,
//		std::atomic_flag* visited,
		std::atomic<unsigned int>* visited, unsigned int iteration,
		std::pair<AdjVertexIterator, AdjVertexIterator>* lookahead) {

	// lookahead phase
	AdjVertexIterator laStart, laEnd;
	for (std::tie(laStart, laEnd) = lookahead[x0]; laStart != laEnd; ++laStart) {
		Vertex y = *laStart;
		if (is_unmatched(y, g, mate)
			&& claim_vertex(y, first_right, visited, iteration)) {
//				!visited[y - first_right].test_and_set()) {

			// update matching
			mate[y] = x0;
			mate[x0] = y;

			lookahead[x0].first = laStart;
			return true;
		}
	}
	if (lookahead[x0].first != laStart) lookahead[x0].first = laStart;

	return false;
}

bool dfs_la(
		const Vertex v,
		const Graph& g,const Vertex first_right, VertexVector& mate,
//		std::atomic_flag* visited,
		std::atomic<unsigned int>* visited, unsigned int iteration,
		std::pair<AdjVertexIterator, AdjVertexIterator>* lookahead,
		std::vector<PathElement>& stack) {

	// do initial lookahead and return if successful ----------------------------------------------
	bool init_lookahead_success = lookahead_step_atomic(v, g, first_right, mate, visited, iteration, lookahead);
	if (init_lookahead_success) return true;
	// --------------------------------------------------------------------------------------------

	stack.clear();
	bool path_found = false;

	PathElement pe;
	std::tie(pe.yiter, pe.yiter_end) = boost::adjacent_vertices(v, g);
	pe.x0 = v;
	stack.push_back(pe);

	while(!stack.empty()) {
		PathElement& stack_top = stack.back();
		Vertex x0 = stack_top.x0;
		AdjVertexIterator& yiter = stack_top.yiter;
		AdjVertexIterator& yiter_end = stack_top.yiter_end;


		if (path_found) {
			// update matching
			const Vertex y = *yiter;
			mate[y] = x0;
			mate[x0] = y;
			// return e.g. pop stack
			stack.pop_back();
			continue;
		}


		// skip all visited vertices
//		while (yiter != yiter_end && visited[*yiter - first_right].test_and_set()) yiter++;
		while (yiter != yiter_end && !claim_vertex(*yiter, first_right, visited, iteration)) yiter++;

		// do dfs step on first unvisited neighbor
		if (yiter != yiter_end) { // if there are still neighbours to visit
			Vertex y = *yiter;
			Vertex x1 = mate[y];

			bool lookahead_success = lookahead_step_atomic(x1, g, first_right, mate, visited, iteration, lookahead);
			if (lookahead_success) {
				path_found = true;
				continue;
			}

			PathElement pe;
			std::tie(pe.yiter, pe.yiter_end) = boost::adjacent_vertices(x1, g);
			pe.x0 = x1;
			stack.push_back(pe);
			continue;
		} else {
			// pop stack otherwise, did not find any path and there are not more neighbors
			stack.pop_back();
		}

	}

	return path_found;
}

