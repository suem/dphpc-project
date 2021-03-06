// PPF implementation based on the Azad. et. al. paper with Test-and-test-and-Set


#include <iostream>
#include <stack>
#include <omp.h>
#include <atomic>

#include "graphtypes.h"
#include "ppf3.h"

inline bool lookahead_step(
		const Vertex x0,
		const Graph& g,const Vertex first_right, VertexVector& mate,
		std::vector<std::atomic_size_t>& visited,
		size_t iteration,
		std::vector<Lookahead>& lookahead
) {

	// lookahead phase
	AdjVertexIterator laStart, laEnd;
	for (std::tie(laStart, laEnd) = lookahead[x0]; laStart != laEnd; ++laStart) {
		Vertex y = *laStart;
		if (is_unmatched(y, mate) && visited[y - first_right] < iteration) {
			if (std::atomic_exchange(&visited[y - first_right], iteration) < iteration) {
				// update matching
				mate[y] = x0;
				//mate[x0] = y; // this store can be done when we are done

				lookahead[x0].first = laStart;
				return true;
			}
		}
	}
	if (lookahead[x0].first != laStart) lookahead[x0].first = laStart;

	return false;
}

inline bool dfs_la_atomic(
		const Vertex v,
		const Graph& g,const Vertex first_right, VertexVector& mate,
        std::vector<std::atomic_size_t>& visited,
		size_t iteration,
		std::vector<Lookahead>& lookahead,
		std::vector<PathElement>& stack) {

	// do initial lookahead and return if successful ----------------------------------------------
	bool init_lookahead_success = lookahead_step(v, g, first_right, mate, visited, iteration, lookahead);
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
			//mate[x0] = y; // this write can be done later when we are done
			// return e.g. pop stack
			stack.pop_back();
			continue;
		}


		// skip all visited vertices
//		while (yiter != yiter_end && std::atomic_exchange(&visited[*yiter - first_right], iteration) >= iteration) {
//			yiter++;
//		}
		while (yiter != yiter_end) {
			const Vertex y_index = *yiter - first_right;
            if (visited[y_index] < iteration) {
				if (std::atomic_exchange(&visited[y_index], iteration) < iteration) {
					break;
				}
			}
			yiter++;
        }

		// do dfs step on first unvisited neighbor
		if (yiter != yiter_end) { // if there are still neighbours to visit
			Vertex y = *yiter;
			Vertex x1 = mate[y];

			bool lookahead_success = lookahead_step(x1, g, first_right, mate, visited, iteration, lookahead);
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
            // move yiter to next vertex for next iteration
            if (!stack.empty()) stack.back().yiter++;
		}

	}

	return path_found;
}

void ppf3(const Graph& g, Vertex first_right, VertexVector& mate, int numThreads) {

	const int NO_THREADS = numThreads; 

	const vertex_size_t n = boost::num_vertices(g);
	const vertex_size_t n_right = n - first_right;

	const int nt = std::min(static_cast<int>(n), NO_THREADS);

	volatile bool path_found;

	std::vector<std::atomic_size_t> visited(n_right);
	memset(&visited[0], 0, sizeof(visited[0]) * n_right);

    size_t iteration = 1;


    // initialize lookahead
    // collect all unmatched
	std::vector<Lookahead> lookahead(first_right);
    std::vector<UnmatchedVertex> unmatched;
	unmatched.reserve(first_right);

    for (Vertex v = 0; v < first_right; v++) {
        lookahead[v] = boost::adjacent_vertices(v, g);

        if (is_unmatched(v, mate)) {
            unmatched.push_back(UnmatchedVertex(v));
        }
    }
	size_t unmatched_size = unmatched.size();
    std::vector<std::vector<PathElement>> stacks(nt);

	do {
		path_found = false;

        iteration++;
        if (iteration == 0) { // integer overflow on iteration, reset visited flags
			memset(&visited[0], 0, sizeof(visited[0]) * n_right);
            iteration = 1;
        }

#pragma omp parallel num_threads(nt)
		{
			std::vector<PathElement> &stack = stacks[omp_get_thread_num()];
#pragma omp for
			for (int i = 0; i < unmatched_size; i++) {
				auto &urv = unmatched[i];

				// skip if vertex is already matched
				if (!urv.unmatched) continue;

				Vertex v = urv.x;

				bool path_found_v = dfs_la_atomic(v, g, first_right, mate, visited, iteration, lookahead, stack);
				if (path_found_v) {
					urv.unmatched = false;
					if (!path_found) path_found = true;
				}
			}
		}

	} while (path_found);


#pragma omp parallel num_threads(nt)
#pragma omp for
    // update matchings for left vertices
    for (int y = first_right; y < n; y++) if (is_matched(y, mate)) mate[mate[y]] = static_cast<int>(y);
}

