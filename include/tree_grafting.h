
#pragma once

#include <atomic>
#include "graphtypes.h"
#include "graphutils.h"

// a thread safe vertex vector for push back operations.
class TSVertexVector {
public:
	TSVertexVector(size_t capacity) : m_capacity(capacity), m_vector(capacity) { m_size = 0; }
	~TSVertexVector() { /*delete[] m_vector;*/ }

	const Vertex& operator[] (size_t n) { return m_vector[n]; }

	TSVertexVector& operator= (const TSVertexVector& obj) {
		assert(m_capacity == obj.capacity());
		m_size = obj.size();
		std::copy(obj.m_vector.begin(), obj.m_vector.begin() + m_size, m_vector.begin());
		return *this;
	}

	size_t size() const { return m_size; }
	size_t capacity() const { return m_capacity; } 
	bool empty() const { return m_size == 0; }
	void clear() { m_size = 0; }

	void push_back(Vertex v) {
		size_t old = std::atomic_fetch_add(&m_size, 1);
		m_vector[old] = v;
	}

private:
	//Vertex* m_vector;
	VertexVector m_vector;
	std::atomic_size_t m_size;
	const size_t m_capacity;
};

class TreeGrafting {
public:
	TreeGrafting(const Graph& g, Vertex first_right, VertexVector& mate, int numThreads);
	~TreeGrafting();

private:
	// initializes member arrays
	void init();
	// execute the tree grafting
	void ms_bfs_graft();

	// breadth first search from above
	void top_down_bfs(TSVertexVector& F);
	// breadth first search from below
	void bottom_up_bfs(TSVertexVector& R, TSVertexVector& ret);
	// rebuild frontier (queue) for the next phase
	void graft(TSVertexVector& ret);
	// update pointers in bfs traversals
	void updatePointers(const Vertex x, const Vertex y, TSVertexVector& queue);
	// find and augment a path
	bool find_path_tg(const Vertex x0, Lookahead* lookahead, std::vector<PathElement>& stack);	
	bool lookahead_step_atomic(const Vertex x0, Lookahead* lookahead);

private:
	int m_numThreads;

	const Graph& m_graph;
	VertexVector& m_mate;
	const Vertex m_firstRight;
	Vertex n_right;
	size_t n;

	VertexVector m_parent;
	VertexVector m_leaf;
	VertexVector m_root;

	TSVertexVector renewableY;
	TSVertexVector activeX;
	TSVertexVector activeY;

	bool* m_visited = nullptr;
	std::atomic_flag* m_augmentVisited = nullptr;

	static const float ALPHA;
};

