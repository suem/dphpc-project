#include <iostream>
#include <stack>
#include <omp.h>
#include <atomic>

#include "tree_grafting.h"

// value taken from the paper
const float TreeGrafting::ALPHA = 5.f;

// whether the Vertex v is in an active tree
inline bool is_active_tree(Vertex v, const VertexVector& root, const VertexVector& leaf) {
	return root[v] != Graph::null_vertex() && leaf[root[v]] == Graph::null_vertex();
}

// whether the Vertex v is in a renewable tree
inline bool is_renewable_tree(Vertex v, const VertexVector& root, const VertexVector& leaf) {
	return root[v] != Graph::null_vertex() && leaf[root[v]] != Graph::null_vertex();
}

TreeGrafting::TreeGrafting(const Graph& g, Vertex first_right, VertexVector& mate, int numThreads) :
	m_graph(g), m_mate(mate), m_firstRight(first_right), n(boost::num_vertices(g)), 
	activeX(n), activeY(n), renewableY(n) {
	n_right = n - first_right;
	m_numThreads = std::min(static_cast<int>(n), numThreads);

	init();
	ms_bfs_graft();
}

TreeGrafting::~TreeGrafting() {
	delete[] m_visited;
	delete m_augmentVisited;
}

void TreeGrafting::init() {
	m_visited = new bool[n];

	m_augmentVisited = new std::vector<std::atomic_size_t>(n_right);
	memset(&(*m_augmentVisited)[0], 0, sizeof((*m_augmentVisited)[0]) * n_right);

	// initialize the root vector
	m_root.resize(n);
	memset(&m_root[0], Graph::null_vertex(), sizeof(m_root[0]) * n);
	// initialize the parent vector
	m_parent.resize(n);
	memset(&m_parent[0], Graph::null_vertex(), sizeof(m_parent[0]) * n);
	// initialize the leaf vector
	m_leaf.resize(n);
	memset(&m_leaf[0], Graph::null_vertex(), sizeof(m_leaf[0]) * n);
	// initialize the visited array
	memset(m_visited, 0, sizeof(bool) * n);
}

void TreeGrafting::ms_bfs_graft() {
	volatile bool path_found;
	size_t iteration = 1;

	// F <- all unmatched X vertices
	// for all those unmatched vertices, set the root to itself
	TSVertexVector F(n);
	TSVertexVector R(n);
#pragma omp parallel num_threads(m_numThreads)
#pragma omp for
	for (int x = 0; x < m_firstRight; ++x) {
		if (is_unmatched(x, m_mate)) {
			F.push_back(x);
			m_root[x] = x;
		}
	}

	// initialize stack
	std::vector<PathElement> stack;
	// initialize lookahead
	Lookahead* lookahead = new Lookahead[m_firstRight];
#pragma omp parallel num_threads(m_numThreads)
#pragma omp for
	for (int i = 0; i < m_firstRight; i++) {
		lookahead[i] = boost::adjacent_vertices(i, m_graph);
	}

	do {
		path_found = false;

		// construct alternating BFS Forest
		while (!F.empty()) {
			// find numUnvisited
			size_t numUnvisitedY = 0;
			for (Vertex y = m_firstRight; y < n; ++y) {
				if (!m_visited[y]) {
					++numUnvisitedY;
				}
			}

			if (F.size() < numUnvisitedY / ALPHA) {
				top_down_bfs(F);
			}
			else {
				// fill unvisited Y vertices
				R.clear();
				for (Vertex y = m_firstRight; y < n; ++y) {
					if (!m_visited[y]) {
						R.push_back(y);
					}
				}
				bottom_up_bfs(R, F);
			}
		}

		iteration++;
		if (iteration == 0) { // integer overflow on iteration, reset visited flags
			memset(&(*m_augmentVisited)[0], 0, sizeof((*m_augmentVisited)[0]) * n_right);
			iteration = 1;
		}

		// step 2: augment matching. 
#pragma omp parallel num_threads(m_numThreads) private(stack)
#pragma omp for
		for (int x = 0; x < m_firstRight; ++x) {
			if (is_unmatched(x, m_mate)) {
				// if an augmenting path P from x is found then augment matching by P
				bool path_found_v = find_path_tg(x, iteration, lookahead, stack);
				if (path_found_v && !path_found) path_found = true;
			}
		}

		// step 3: construct frontier for the next phase
		if (path_found) {
			graft(F);
		}

	} while (path_found);
}

void TreeGrafting::top_down_bfs(TSVertexVector& F) {
	activeX.clear();

#pragma omp parallel num_threads(m_numThreads)
#pragma omp for
	for (int i = 0; i < F.size(); ++i) {
		Vertex x = F[i];
		if (is_active_tree(x, m_root, m_leaf)) {
			m_leaf[m_root[x]] = Graph::null_vertex();
			AdjVertexIterator start, end;
			for (std::tie(start, end) = boost::adjacent_vertices(x, m_graph); start != end; ++start) {
				Vertex y = *start;
				if (!m_visited[y]) {
					updatePointers(x, y, activeX);
				}
			}
		}
	}

	F = activeX;
}

void TreeGrafting::bottom_up_bfs(TSVertexVector& R, TSVertexVector& ret) {
	ret.clear();

#pragma omp parallel num_threads(m_numThreads)
#pragma omp for
	for (int i = 0; i < R.size(); ++i) {
		Vertex y = R[i];
		AdjVertexIterator start, end;
		for (std::tie(start, end) = boost::adjacent_vertices(y, m_graph); start != end; ++start) {
			Vertex x = *start;
			if (is_active_tree(x, m_root, m_leaf)) {
				m_leaf[m_root[x]] = Graph::null_vertex();
				updatePointers(x, y, ret);
				break;
			}
		}
	}
}

void TreeGrafting::graft(TSVertexVector& ret) {
	ret.clear();

	activeX.clear();
#pragma omp parallel num_threads(m_numThreads)
#pragma omp for
	for (int x = 0; x < m_firstRight; ++x) {
		if (is_active_tree(x, m_root, m_leaf)) {
			activeX.push_back(x);
		}
	}

	activeY.clear();
	renewableY.clear();
#pragma omp parallel num_threads(m_numThreads)
#pragma omp for
	for (int y = m_firstRight; y < n; ++y) {
		if (is_active_tree(y, m_root, m_leaf)) {
			activeY.push_back(y);
		}
		if (is_renewable_tree(y, m_root, m_leaf)) {
			renewableY.push_back(y);
		}
	}

#pragma omp parallel num_threads(m_numThreads)
#pragma omp for
	for (int i = 0; i < renewableY.size(); ++i) {
		Vertex y = renewableY[i];
		m_visited[y] = false;
		m_root[y] = Graph::null_vertex();
	}

	if (activeX.size() > renewableY.size() / ALPHA) {
		bottom_up_bfs(renewableY, ret);
	}
	else {
		// queue(ret) <- unmatched X vertices 
#pragma omp parallel num_threads(m_numThreads)
#pragma omp for
		for (int x = m_firstRight; x < n; ++x) {
			if (is_unmatched(x, m_mate)) {
				ret.push_back(x);
			}
		}

#pragma omp parallel num_threads(m_numThreads)
#pragma omp for
		for (int i = 0; i < activeY.size(); ++i) {
			Vertex y = activeY[i];
			m_visited[y] = false;
			m_root[y] = Graph::null_vertex();
		}

#pragma omp parallel num_threads(m_numThreads)
#pragma omp for
		for (int i = 0; i < activeX.size(); ++i) {
			// if x is in unmatchedX, exclude it
			Vertex x = activeX[i];
			if (is_unmatched(x, m_mate)) continue;
			m_root[x] = Graph::null_vertex();
		}
	}
}

void TreeGrafting::updatePointers(const Vertex x, const Vertex y, TSVertexVector& queue) {
	m_parent[y] = x;
	m_visited[y] = true;
	m_root[y] = m_root[x];
	if (is_matched(y, m_mate)) {
		queue.push_back(m_mate[y]);
		m_root[m_mate[y]] = m_root[y];
	}
	else {
		// an augmenting path is found, so we end the augmenting path
		m_leaf[m_root[x]] = y;
	}
}

bool TreeGrafting::find_path_tg(const Vertex v,
	size_t iteration,
	Lookahead* lookahead,
	std::vector<PathElement>& stack) {

	// do initial lookahead and return if successful ----------------------------------------------
	bool init_lookahead_success = lookahead_step(v, iteration, lookahead);
	if (init_lookahead_success) return true;
	// --------------------------------------------------------------------------------------------

	stack.clear();
	bool path_found = false;

	PathElement pe;
	std::tie(pe.yiter, pe.yiter_end) = boost::adjacent_vertices(v, m_graph);
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
			m_mate[y] = x0;
			m_mate[x0] = y;
			// return e.g. pop stack
			stack.pop_back();
			continue;
		}


		// skip all visited vertices
		while (yiter != yiter_end) {
			const Vertex y_index = *yiter - m_firstRight;
			if ((*m_augmentVisited)[y_index] < iteration) {
				if (std::atomic_exchange(&(*m_augmentVisited)[y_index], iteration) < iteration) {
					break;
				}
			}
			yiter++;
		}

		// do dfs step on first unvisited neighbor
		if (yiter != yiter_end) { // if there are still neighbours to visit
			Vertex y = *yiter;
			Vertex x1 = m_mate[y];

			bool lookahead_success = lookahead_step(x1, iteration, lookahead);
			if (lookahead_success) {
				path_found = true;
				continue;
			}

			PathElement pe;
			std::tie(pe.yiter, pe.yiter_end) = boost::adjacent_vertices(x1, m_graph);
			pe.x0 = x1;
			stack.push_back(pe);
			continue;
		}
		else {
			// pop stack otherwise, did not find any path and there are not more neighbors
			stack.pop_back();
			// move yiter to next vertex for next iteration
			if (!stack.empty()) stack.back().yiter++;
		}

	}

	return path_found;
}

bool TreeGrafting::lookahead_step(
	const Vertex x0,
	size_t iteration,
	Lookahead* lookahead) {

	// lookahead phase
	AdjVertexIterator laStart, laEnd;
	for (std::tie(laStart, laEnd) = lookahead[x0]; laStart != laEnd; ++laStart) {
		Vertex y = *laStart;
		if (is_unmatched(y, m_mate) && (*m_augmentVisited)[y - m_firstRight] < iteration) {
			if (std::atomic_exchange(&(*m_augmentVisited)[y - m_firstRight], iteration) < iteration) {
				// update matching
				m_mate[y] = x0;
				m_mate[x0] = y;
				lookahead[x0].first = laStart;
				return true;
			}
		}
	}
	if (lookahead[x0].first != laStart) lookahead[x0].first = laStart;

	return false;
}

