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
#include "Timer.h"

union Flag {
    std::atomic_flag flag;
    bool flag_bool;
};

inline bool lookahead_step_atomic(
		const Vertex x0,
		const Graph& g, const Vertex first_right, VertexVector& mate,
		Flag* visited,
		std::vector<bool>& visited_local,
		std::pair<AdjVertexIterator, AdjVertexIterator>* lookahead) {

    bool found_unmatched = false;

	// lookahead phase
	AdjVertexIterator laStart, laEnd, laStartOld;
    std::tie(laStart, laEnd) = lookahead[x0];
    laStartOld = laStart;

	for (; laStart != laEnd; laStart++) {
		Vertex y = *laStart;
		if (is_unmatched(y, mate) && !visited[y - first_right].flag.test_and_set()) {

            visited_local[y - first_right] = false;
            // update matching
            mate[y] = x0;
            mate[x0] = y;
            found_unmatched = true;
            laStart++;
            visited[y - first_right].flag.clear();
            break;
            visited[y - first_right].flag.clear();
		} }
	lookahead[x0].first = laStart;

	return found_unmatched;
}

bool dfs_la(
		const Vertex v,
		const Graph& g, const Vertex first_right, VertexVector& mate,
		Flag* visited,
		std::vector<bool>& visited_local,
		std::pair<AdjVertexIterator, AdjVertexIterator>* lookahead,
		std::vector<PathElement>& stack) {

	// do initial lookahead and return if successful ----------------------------------------------
	bool init_lookahead_success = lookahead_step_atomic(v, g, first_right, mate, visited, visited_local, lookahead);
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
            visited[y - first_right].flag.clear();
			// return e.g. pop stack
			stack.pop_back();
			continue;
		} 

		// skip all visited vertices
		while (yiter != yiter_end) {
            const Vertex y = *yiter;
            // if y has not been visited in this dfs and if we can claim it, continue from there
            if (!visited_local[y - first_right] ) {
                visited_local[y - first_right] = true;
                if(!visited[y - first_right].flag.test_and_set()) break;
            }
            yiter++;
        }

		// do dfs step on first unvisited neighbor
		if (yiter != yiter_end) { // if there are still neighbours to visit
			Vertex y = *yiter;
			Vertex x1 = mate[y];

			bool lookahead_success = lookahead_step_atomic(x1, g, first_right, mate, visited, visited_local, lookahead);
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
            if(!stack.empty()) visited[*(stack.back().yiter++) - first_right].flag.clear();        
		}
	}

	return path_found;
}


void unsync_parallel_pothen_fan(const Graph& g, Vertex first_right, VertexVector& mate, int numThreads) {

	const int NO_THREADS = numThreads; 

	const vertex_size_t n = boost::num_vertices(g);
	const vertex_size_t n_right = n - first_right;

	const int nt = std::min(static_cast<int>(n), NO_THREADS);

	Flag* visited = new Flag[n_right];
	memset(visited, 0, sizeof(Flag) * n_right);

    // initialize lookahead
    std::pair<AdjVertexIterator, AdjVertexIterator>* lookahead = new std::pair<AdjVertexIterator, AdjVertexIterator>[first_right];
    for (Vertex i = 0; i < first_right; i++) lookahead[i] = boost::adjacent_vertices(i, g);

    std::vector<PathElement> stack;
    std::vector<bool> visited_local;
	volatile bool path_found;

	do {
		path_found = false;
#pragma omp parallel num_threads(nt) private(stack, visited_local)
        {
            visited_local.assign(n_right, false);
            bool found_local;
            found_local = false;
#pragma omp for
                for (int v = 0; v < first_right; v++) {

                    // skip if vertex is already matched
                    if (is_matched(v, mate)) {
                        continue;
                    }

                    bool path_found_v = dfs_la(v, g, first_right, mate, visited, visited_local, lookahead, stack);
                    if (path_found_v) {
                        //matching_size++;
                        if (!path_found) path_found = true;
                        if (!found_local) found_local = true;
                    }
                }
        }

	} while (path_found);

	delete[] visited;
}



