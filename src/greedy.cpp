#include "greedy.h"

#include "graphtypes.h"
#include "graphutils.h"

/**
 * Simple greedy matching that first tries to match all vertices with degree = 1
 * and then all other
 */
void simple_greedy_matching(const Graph& g, VertexVector& matching) {

    vertex_size_t n = boost::num_vertices(g);
    Vertex null_vertex = Graph::null_vertex();

    matching.assign(n, null_vertex);

    // match all vertices of degree 1
    for (Vertex v = 0; v < n; v++) {
        if (boost::degree(v, g) != 1 || is_matched(v, matching)) continue;
        AdjVertexIterator s, e;
        for (std::tie(s, e) = boost::adjacent_vertices(v, g); s != e; s++) {
            Vertex u = *s;
            if (is_matched(u, matching)) continue;
            matching[v] = u;
            matching[u] = v;
        }
    }

    // do greedy matching over all edges
    EdgeIterator start, end;
    for (std::tie(start, end) = boost::edges(g); start != end; ++start) {
        Edge e = *start;
        Vertex u = boost::source(e, g);
        Vertex v = boost::target(e, g);
        if (matching[u] == null_vertex && matching[v] == null_vertex) {
            matching[u] = v;
            matching[v] = u;
        }
    }

}
