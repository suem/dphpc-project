// Sequential version of Pothen-Fan

#include <iostream>
#include <stack>
#include <omp.h>
#include <atomic>

#include "graphtypes.h"
#include "pf.h"

bool lookahead_step(
		const Vertex x0,
		const Graph& g, const Vertex first_right, VertexVector& mate,
		std::vector<bool>& visited,
		std::vector<Lookahead>& lookahead) {

	// lookahead phase
	AdjVertexIterator laStart, laEnd;
	for (std::tie(laStart, laEnd) = lookahead[x0]; laStart != laEnd; ++laStart) {
		Vertex y = *laStart;
		if (is_unmatched(y, mate) && !visited[y - first_right]) {
			visited[y - first_right] = true;

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
		const Graph& g, const Vertex first_right,
		VertexVector& mate,
		std::vector<bool>& visited,
		std::vector<Lookahead>& lookahead,
		std::vector<PathElement>& stack) {

	// do initial lookahead and return if successful ----------------------------------------------
	bool init_lookahead_success = lookahead_step(v, g, first_right, mate, visited, lookahead);
	if (init_lookahead_success) return true;
	// --------------------------------------------------------------------------------------------

	stack.clear();
	bool path_found = false;

	PathElement pe;
	std::tie(pe.yiter, pe.yiter_end) = boost::adjacent_vertices(v, g);
	pe.x0 = v;
	stack.push_back(pe);

	while (!stack.empty()) {
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
		while (yiter != yiter_end && visited[*yiter - first_right]) yiter++;

		// do dfs step on first unvisited neighbor
		if (yiter != yiter_end) { // if there are still neighbours to visit
			Vertex y = *yiter;
			Vertex x1 = mate[y];
			visited[y - first_right] = true; // set as visited

			bool lookahead_success = lookahead_step(x1, g, first_right, mate, visited, lookahead);
			if (lookahead_success) {
				path_found = true;
				continue;
			}

			PathElement pe;
			std::tie(pe.yiter, pe.yiter_end) = boost::adjacent_vertices(x1, g);
			pe.x0 = x1;
			stack.push_back(pe);
			continue;
		}
		else {
			// pop stack otherwise, did not find any path and there are not more neighbors
			stack.pop_back();
			if (!stack.empty()) stack.back().yiter++;
		}

	}

	return path_found;
}

void pf(const Graph& g, const Vertex first_right, VertexVector& mate) {

	const vertex_size_t  n = boost::num_vertices(g);
	const vertex_size_t  n_right = n - first_right;

	bool path_found;

	std::vector<bool> visited(n_right);

    // init stack
	std::vector<PathElement> stack;

    // initialize lookahead
	// collect all unmatched vertices

	std::vector<Lookahead> lookahead(first_right);
    std::vector<Vertex> unmatched;
	unmatched.reserve(first_right);
	for (Vertex v = 0; v < first_right; v++) {
        lookahead[v] = boost::adjacent_vertices(v, g);
		if (is_unmatched(v, mate)) unmatched.push_back(v);
	}
	size_t unmatched_size = unmatched.size();

    do {
		path_found = false; // reset path_found
        std::fill(visited.begin(), visited.end(), false);

        // iterate over all unmatched left vertices
        for (size_t i = 0; i < unmatched_size; i++) {
			Vertex x0 = unmatched[i];

			// skip if vertex is already matched
			if (is_matched(x0, mate)) continue;
			bool path_found_v = dfs_la(x0, g, first_right, mate, visited, lookahead, stack);
            if (path_found_v && !path_found) path_found = true;
		}
	} while(path_found);
}


