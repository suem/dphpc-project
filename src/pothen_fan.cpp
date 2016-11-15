//
// Created by suem on 10/27/16.
//
#include <iostream>
#include <stack>
#include <omp.h>
#include <atomic>

#include "graphtypes.h"
#include "pothen_fan.h"

void parallel_pothen_fan(const Graph& g, Vertex first_right, VertexVector& mate, int numThreads) {

	const int NO_THREADS = numThreads; 
//	const int NO_THREADS = omp_get_max_threads();

	const vertex_size_t n = boost::num_vertices(g);
	const vertex_size_t n_right = n - first_right;

	const int nt = std::min(static_cast<int>(n), NO_THREADS);

	volatile bool path_found;

	std::atomic_flag* visited = new std::atomic_flag[n_right];

    // initialize lookahead
    std::pair<AdjVertexIterator, AdjVertexIterator>* lookahead = new std::pair<AdjVertexIterator, AdjVertexIterator>[first_right];
    for (Vertex i = 0; i < first_right; i++) lookahead[i] = boost::adjacent_vertices(i, g);

	do {
		path_found = false;

		// here we don't need atomic clears to reset the flags
        memset(visited, 0, sizeof(std::atomic_flag) * n_right);

#pragma omp parallel num_threads(nt)
#pragma omp for
		for (Vertex v = 0; v < first_right; v++) {
			// if any thread  has found a path (including us) stop.
//				if  (path_found) v = first_right; // simulate break (not supported in omp for

			// skip if vertex is already matched
			if (is_matched(v, g, mate))  continue;

//			bool path_found_v = find_path_atomic(v, g, first_right, mate, visited);
//			bool path_found_v = find_path_recursive_atomic(v, g, first_right, mate, visited);
            bool path_found_v = find_path_la_recursive_atomic(v, g, first_right, mate, visited, lookahead);
			if (path_found_v && !path_found) path_found = true;
		}

	} while (path_found);

	delete[] visited;
}

bool find_path_atomic(const Vertex x0, const Graph& g, Vertex first_right, VertexVector& mate, std::atomic_flag* visited) {

	std::stack<FindPathElement> stack;
	FindPathElement e1;
	e1.x0 = x0;
	std::tie(e1.start, e1.end) = boost::adjacent_vertices(e1.x0, g);
	stack.push(e1);

	while (!stack.empty()) {
		FindPathElement& vars = stack.top();
		if (vars.found) {
			mate[vars.y] = vars.x0;
			mate[vars.x0] = vars.y;

			// handle return value
			stack.pop();
			if (stack.empty()) {
				return true;
			}

			stack.top().found = true;
			continue;
		}

		bool leaveWhile = false;

		for (; vars.start != vars.end; ++vars.start) {
			vars.y = *vars.start;

            // skip vertex if this was alredy visited by another thread
			if (visited[vars.y - first_right].test_and_set()) continue;

			leaveWhile = true;

			if (is_unmatched(vars.y, g, mate)) { // y is unmatched
				mate[vars.y] = vars.x0;
				mate[vars.x0] = vars.y;

				// handle return value
				stack.pop();
				if (stack.empty()) {
					return true;
				}

				stack.top().found = true;
				break;
			}

			// y is matched with x1
			FindPathElement e2;
			e2.x0 = mate[vars.y];
			std::tie(e2.start, e2.end) = boost::adjacent_vertices(e2.x0, g);
			stack.push(e2);

			break;
		}

		if (leaveWhile) {
			continue;
		}

		if (!stack.empty()) {
			stack.pop();
		}
	}

	return false;
}


bool find_path_recursive_atomic(const Vertex x0, const Graph& g, const Vertex first_right, VertexVector& mate, std::atomic_flag* visited) {

	AdjVertexIterator start, end;

	for (std::tie(start, end) = boost::adjacent_vertices(x0, g); start != end; ++start) {
		Vertex y = *start;

		// check if vertex y has already been visited in this DFS
		// NOTE: we subtract first_right because the visited array only stores the flag for all right vertices
		if (visited[y - first_right].test_and_set()) continue;

		if (is_unmatched(y, g, mate)) { // y is unmatched, we found an augmenting path

			// update matching while returning from DFS
			mate[y] = x0;
			mate[x0] = y;

			return true;
		}

		// y is matched with x1
		Vertex x1 = mate[y];

		bool path_found = find_path_recursive_atomic(x1, g, first_right, mate, visited);

		if (!path_found) {
			continue;
		}

		// we have found a path and can update the matching
		mate[y] = x0;
		mate[x0] = y;

		return true;
	}

	return false;
}


bool find_path_la_recursive_atomic(const Vertex x0, const Graph& g, const Vertex first_right, VertexVector& mate, std::atomic_flag* visited,
                                   std::pair<AdjVertexIterator, AdjVertexIterator>* lookahead) {


    // lookahead phase
    AdjVertexIterator laStart, laEnd;
    for (std::tie(laStart, laEnd) = lookahead[x0]; laStart != laEnd; ++laStart) {
        Vertex y = *laStart;
        if (is_unmatched(y, g, mate) && !visited[y - first_right].test_and_set()) {
            // update matching
            mate[y] = x0;
            mate[x0] = y;

            lookahead[x0].first = laStart;
            return true;
        }
    }

    if (lookahead[x0].first != laStart) lookahead[x0].first = laStart;

    // DFS phase
    AdjVertexIterator start, end;
    for (std::tie(start, end) = boost::adjacent_vertices(x0, g); start != end; ++start) {
        Vertex y = *start;

        // check if vertex y has already been visited in this DFS
        // NOTE: we subtract first_right because the visited array only stores the flag for all right vertices
        if (visited[y - first_right].test_and_set()) continue;

        // DFS
        bool path_found = find_path_la_recursive_atomic(mate[y], g, first_right, mate, visited, lookahead);

        if (!path_found) {
            continue;
        }

        // we have found a path and can update the matching
        mate[y] = x0;
        mate[x0] = y;

        return true;
    }

    return false;
}


void pothen_fan(const Graph& g, const Vertex first_right, VertexVector& mate) {

	const vertex_size_t  n = boost::num_vertices(g);
	const vertex_size_t  n_right = n - first_right;

	bool path_found = true;
	bool* visited = new bool[n_right];

    // initialize lookahead
    std::pair<AdjVertexIterator, AdjVertexIterator>* lookahead = new std::pair<AdjVertexIterator, AdjVertexIterator>[first_right];
    for (Vertex i = 0; i < first_right; i++) lookahead[i] = boost::adjacent_vertices(i, g);

	while (path_found) {
		path_found = false; // reset path_found

		memset(visited, 0, sizeof(bool) * n_right); // set all visited flags to 0
        // iterate over all left vertices
		for (Vertex v = 0; v < first_right; v++) {
			// skip if vertex is already matched
			if (is_matched(v, g, mate)) continue;


			// assert: v is on the left and unmatched
//			path_found = find_path(v, g, first_right, mate, visited) || path_found;
//			path_found = find_path_recursive(v, g, first_right, mate, visited) || path_found;
            path_found = find_path_la_recursive(v, g, first_right, mate, visited, lookahead) || path_found;
		}
	}

	delete[] visited;
}

//todo find out why recursive version is faster
bool find_path(const Vertex x0, const Graph& g, Vertex first_right, VertexVector& mate, bool* visited) {

	std::stack<FindPathElement> stack;
	FindPathElement e1;
	e1.x0 = x0;
	std::tie(e1.start, e1.end) = boost::adjacent_vertices(e1.x0, g);
	stack.push(e1);

	while (!stack.empty()) {
		FindPathElement& vars = stack.top();
		if (vars.found) {
			mate[vars.y] = vars.x0;
			mate[vars.x0] = vars.y;

			// handle return value
			stack.pop();
			if (stack.empty()) {
				return true;
			}

			stack.top().found = true;
			continue;
		}

		bool leaveWhile = false;

		for (; vars.start != vars.end; ++vars.start) {
			vars.y = *vars.start;

			if (visited[vars.y - first_right]) continue;
			visited[vars.y - first_right] = true;

			leaveWhile = true;

			if (is_unmatched(vars.y, g, mate)) { // y is unmatched
				mate[vars.y] = vars.x0;
				mate[vars.x0] = vars.y;

				// handle return value
				stack.pop();
				if (stack.empty()) {
					return true;
				}

				stack.top().found = true;
				break;
			}

			// y is matched with x1
			FindPathElement e2;
			e2.x0 = mate[vars.y];
			std::tie(e2.start, e2.end) = boost::adjacent_vertices(e2.x0, g);
			stack.push(e2);

			break;
		}

		if (leaveWhile) {
			continue;
		}

		if (!stack.empty()) {
			stack.pop();
		}
	}

	return false;
}

bool find_path_recursive(const Vertex x0, const Graph& g, const Vertex first_right, VertexVector& mate, bool* visited) {

	AdjVertexIterator start, end;

	for (std::tie(start, end) = boost::adjacent_vertices(x0, g); start != end; ++start) {
		Vertex y = *start;

        // check if vertex y has already been visited in this DFS
		// NOTE: we subtract first_right because the visited array only stores the flag for all right vertices
		if (visited[y - first_right]) continue;
		visited[y - first_right] = true;

		if (is_unmatched(y, g, mate)) { // y is unmatched, we found an augmenting path

			// update matching while returning from DFS
			mate[y] = x0;
			mate[x0] = y;

			return true;
		}

		// y is matched with x1
		Vertex x1 = mate[y];

		bool path_found = find_path_recursive(x1, g, first_right, mate, visited);

		if (!path_found) {
			continue;
		}

		// we have found a path and can update the matching
		mate[y] = x0;
		mate[x0] = y;

		return true;
	}

	return false;
}

bool find_path_la_recursive(const Vertex x0, const Graph& g, const Vertex first_right, VertexVector& mate, bool* visited, std::pair<AdjVertexIterator, AdjVertexIterator>* lookahead) {

    // lookahead phase
    AdjVertexIterator laStart, laEnd;
    for (std::tie(laStart, laEnd) = lookahead[x0]; laStart != laEnd; ++laStart) {
        Vertex y = *laStart;
        if (is_unmatched(y, g, mate) && !visited[y - first_right]) {
            visited[y - first_right] = true;

            // update matching
            mate[y] = x0;
            mate[x0] = y;

			lookahead[x0].first = laStart;
			return true;
		}
	}

	if (lookahead[x0].first != laStart) lookahead[x0].first = laStart;


    // DFS phase
    AdjVertexIterator start, end;
    for (std::tie(start, end) = boost::adjacent_vertices(x0, g); start != end; ++start) {
        Vertex y = *start;

        // check if vertex y has already been visited in this DFS
        // NOTE: we subtract first_right because the visited array only stores the flag for all right vertices
        if (visited[y - first_right]) continue;
        visited[y - first_right] = true;

        // recursive call
        bool path_found = find_path_la_recursive(mate[y], g, first_right, mate, visited, lookahead);

        if (!path_found) {
            continue;
        }

        // we have found a path and can update the matching
        mate[y] = x0;
        mate[x0] = y;

        return true;
    }

    return false;

}


