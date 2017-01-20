// NOTE: this was the first attempt of implementing Pothen-Fan is is only here for historical reference. See ppf3 for the currently best implementation
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

	//std::atomic_flag* visited = new std::atomic_flag[n_right];
	std::atomic<unsigned char>* visited = new std::atomic<unsigned char>[n_right];

    // initialize lookahead
    Lookahead* lookahead = new Lookahead[first_right];
#pragma omp parallel num_threads(nt)
#pragma omp for
    for (int i = 0; i < first_right; i++) lookahead[i] = boost::adjacent_vertices(i, g);

    // collect all unmatched
    std::vector<Vertex> unmatched;
	unmatched.reserve(first_right);
	for (Vertex v = 0; v < first_right; v++) if (is_unmatched(v, mate)) unmatched.push_back(v);
	size_t unmatched_size = unmatched.size();

	std::vector<std::vector<PathElement>> stacks(nt);

	do {
		path_found = false;

        memset(visited, 0, sizeof(std::atomic<unsigned char>) * n_right);

#pragma omp parallel num_threads(nt)
		{
			std::vector<PathElement>& stack = stacks[omp_get_thread_num()];
#pragma omp for
			for (int i = 0; i < unmatched_size; i++) {
				Vertex v = unmatched[i];

				// skip if vertex is already matched
				if (is_matched(v, mate)) continue;

				bool path_found_v = dfs_la_atomic(v, g, first_right, mate, visited, lookahead, stack);
				if (path_found_v && !path_found) path_found = true;
			}
		}

	} while (path_found);

	delete[] visited;
	delete[] lookahead;
}

void pothen_fan(const Graph& g, const Vertex first_right, VertexVector& mate) {

	const vertex_size_t  n = boost::num_vertices(g);
	const vertex_size_t  n_right = n - first_right;

	bool path_found;

	bool* visited = new bool[n_right];

    // init stack
	std::vector<PathElement> stack;

    // initialize lookahead
    Lookahead* lookahead = new Lookahead[first_right];
    for (Vertex i = 0; i < first_right; i++) lookahead[i] = boost::adjacent_vertices(i, g);

	// collect all unmatched vertices
    std::vector<Vertex> unmatched;
	unmatched.reserve(first_right);
	for (Vertex v = 0; v < first_right; v++) if (is_unmatched(v, mate)) unmatched.push_back(v);
	size_t unmatched_size = unmatched.size();

    do {
		path_found = false; // reset path_found
		memset(visited, 0, sizeof(bool) * n_right); // set all visited flags to 0

        // iterate over all unmatched left vertices
        for (size_t i = 0; i < unmatched_size; i++) {
			Vertex x0 = unmatched[i];

			// skip if vertex is already matched
			if (is_matched(x0, mate)) continue;
			bool path_found_v = dfs_la(x0, g, first_right, mate, visited, lookahead, stack);
            if (path_found_v && !path_found) path_found = true;
		}
	} while(path_found);

	delete[] visited;
    delete[] lookahead;
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

			if (is_unmatched(vars.y, mate)) { // y is unmatched
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

		if (is_unmatched(y, mate)) { // y is unmatched, we found an augmenting path

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
	Lookahead* lookahead) {


	// lookahead phase
	AdjVertexIterator laStart, laEnd;
	for (std::tie(laStart, laEnd) = lookahead[x0]; laStart != laEnd; ++laStart) {
		Vertex y = *laStart;
		if (is_unmatched(y, mate) && !visited[y - first_right].test_and_set()) {
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

bool find_path(const Vertex x0, const Graph& g, const Vertex first_right, VertexVector& mate, bool* visited, std::vector<FindPathElement>& stack) {

	stack.clear();

	FindPathElement e1;
	e1.x0 = x0;
	std::tie(e1.start, e1.end) = boost::adjacent_vertices(e1.x0, g);
	stack.push_back(e1);

	while (!stack.empty()) {
		FindPathElement& vars = stack.back();
		if (vars.found) {
			mate[vars.y] = vars.x0;
			mate[vars.x0] = vars.y;

			// handle return value
			stack.pop_back();
			if (stack.empty()) {
				return true;
			}

			stack.back().found = true;
			continue;
		}

		bool leaveWhile = false;

		for (; vars.start != vars.end; ++vars.start) {
			vars.y = *vars.start;

			if (visited[vars.y - first_right]) continue;
			visited[vars.y - first_right] = true;

			leaveWhile = true;

			if (is_unmatched(vars.y, mate)) { // y is unmatched
				mate[vars.y] = vars.x0;
				mate[vars.x0] = vars.y;

				// handle return value
				stack.pop_back();
				if (stack.empty()) {
					return true;
				}

				stack.back().found = true;
				break;
			}

			// y is matched with x1
			FindPathElement e2;
			e2.x0 = mate[vars.y];
			std::tie(e2.start, e2.end) = boost::adjacent_vertices(e2.x0, g);
			stack.push_back(e2);

			break;
		}

		if (leaveWhile) {
			continue;
		}

		if (!stack.empty()) {
			stack.pop_back();
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

		if (is_unmatched(y, mate)) { // y is unmatched, we found an augmenting path

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

bool find_path_la_recursive(const Vertex x0, const Graph& g, const Vertex first_right, VertexVector& mate, bool* visited, Lookahead* lookahead) {

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

bool lookahead_step_atomic(
	const Vertex x0,
	const Graph& g, const Vertex first_right, VertexVector& mate,
	//std::atomic_flag* visited,
	std::atomic<unsigned char>* visited,
	Lookahead* lookahead) {

	// lookahead phase
	AdjVertexIterator laStart, laEnd;
	for (std::tie(laStart, laEnd) = lookahead[x0]; laStart != laEnd; ++laStart) {
		Vertex y = *laStart;
		if (is_unmatched(y, mate) &&
			std::atomic_fetch_add(&visited[y - first_right], (unsigned char)1) == 0) {
			//!visited[y - first_right].test_and_set()) {

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

bool dfs_la_atomic(
	const Vertex v,
	const Graph& g, const Vertex first_right, VertexVector& mate,
	//std::atomic_flag* visited,
	std::atomic<unsigned char>* visited,
	Lookahead* lookahead,
	std::vector<PathElement>& stack) {

	// do initial lookahead and return if successful ----------------------------------------------
	bool init_lookahead_success = lookahead_step_atomic(v, g, first_right, mate, visited, lookahead);
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
		while (yiter != yiter_end
			//&& visited[*yiter - first_right].test_and_set()) {
			&& std::atomic_fetch_add(&visited[*yiter - first_right], (unsigned char)1) != 0) {
			yiter++;
		}

		// do dfs step on first unvisited neighbor
		if (yiter != yiter_end) { // if there are still neighbours to visit
			Vertex y = *yiter;
			Vertex x1 = mate[y];

			bool lookahead_success = lookahead_step_atomic(x1, g, first_right, mate, visited, lookahead);
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

bool lookahead_step(
	const Vertex x0,
	const Graph& g, const Vertex first_right, VertexVector& mate,
	bool* visited,
	Lookahead* lookahead) {

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
	bool* visited,
	Lookahead* lookahead,
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
