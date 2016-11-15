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

	bool path_found;

	vertex_size_t* visited = new vertex_size_t[n_right];
	memset(visited, 0, sizeof(bool) * n_right); // set all visited flags to 0
	vertex_size_t iteration = 0;

    // init parent array
    Vertex* parent = new Vertex[n_right];

    // init stack
	Vertex* stack = new Vertex[first_right];

    // initialize lookahead
    std::pair<AdjVertexIterator, AdjVertexIterator>* lookahead = new std::pair<AdjVertexIterator, AdjVertexIterator>[first_right];
    for (Vertex i = 0; i < first_right; i++) lookahead[i] = boost::adjacent_vertices(i, g);

	// collect all unmatched vertices
    std::vector<Vertex> unmatched;
	unmatched.reserve(first_right / 2);
	for (Vertex v = 0; v < first_right; v++) if (is_unmatched(v, g, mate)) unmatched.push_back(v);
	size_t unmatched_size = unmatched.size();

    do {
		path_found = false; // reset path_found

		iteration++;

        // iterate over all left vertices

        for (size_t i = 0; i < unmatched_size; i++) {
			Vertex x0 = unmatched[i];

			// skip if vertex is already matched
			if (is_matched(x0, g, mate)) continue;

            bool path_found_v = find_path_la(x0, g, first_right, mate, visited, iteration, parent, lookahead, stack);

            if (path_found_v && !path_found) path_found = true;
		}
	} while(path_found);

	delete[] visited;
	delete[] parent;
	delete[] stack;
    delete[] lookahead;
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

bool find_path_la(
		const Vertex x0,
		const Graph& g, Vertex first_right,
		VertexVector& mate,
		vertex_size_t* visited, vertex_size_t iteration,
		Vertex* parent,
		std::pair<AdjVertexIterator, AdjVertexIterator>* lookahead,
		Vertex* stack
) {

	int stack_pointer = 0;
	stack[stack_pointer] = x0;

	bool path_found = false;
	Vertex last_y;

	while(stack_pointer >= 0) {
		Vertex x = stack[stack_pointer];
        stack_pointer--;

		// lookahead phase
		AdjVertexIterator laStart, laEnd;
		for (std::tie(laStart, laEnd) = lookahead[x]; laStart != laEnd; ++laStart) {
			Vertex y = *laStart;
			if (visited[y - first_right] != iteration && is_unmatched(y, g, mate)) {
				visited[y - first_right] = iteration;
				lookahead[x].first = laStart;
                parent[y - first_right] = x;
                path_found = true;
				last_y = y;
                break;
			}
		}
		if (lookahead[x].first != laStart) lookahead[x].first = laStart;
		if (path_found) break;


		// DFS phase
		AdjVertexIterator start, end;
		for (std::tie(start, end) = boost::adjacent_vertices(x, g); start != end; ++start) {
			Vertex y = *start;

			if (visited[y - first_right] == iteration) continue;
			visited[y - first_right] = iteration;

			parent[y - first_right] = x;

            if (is_unmatched(y,g,mate)) {
				path_found = true;
				last_y = y;
				break;
			}

			stack_pointer++;
            stack[stack_pointer] = mate[y];
		}
		if (path_found) break;

	}

	if (path_found) {
        Vertex xp;
		do {
			xp = parent[last_y - first_right];
			Vertex next_y = mate[xp];
			mate[last_y] = xp;
			mate[xp] = last_y;
			last_y = next_y;
		} while (xp != x0);
	}

	return path_found;
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


