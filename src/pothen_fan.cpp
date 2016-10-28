//
// Created by suem on 10/27/16.
//

#include "graphtypes.h"
#include "pothen_fan.h"
#include <iostream>
#include <stack>


void pothen_fan(const Graph& g, VertexVector& mate) {

    vertex_size_t  n = boost::num_vertices(g);

    // create initial greedy matching
    match_greedy(g, mate);

    Vertex null_vertex = g.null_vertex();

    bool path_found;
    bool* visited = new bool[n];

    do {
        path_found = false;
        std::fill_n(visited, n, false);

        VertexIterator start, end;
        for (std::tie(start, end) = boost::vertices(g); start != end; start++) {
            const Vertex v = *start;
            if (is_right(v) || mate[v] != null_vertex) {
                // skip if vertex is on the right or if vertex is already matched
                continue;
            }
            // assert: v is on the left and unmatched
            path_found = find_path(v, g, mate, visited);
        }
    } while (path_found);

    delete[] visited;

}

bool find_path(const Vertex& x0, const Graph& g, VertexVector& mate, bool* visited) {

    AdjVertexIterator start, end;

    for (std::tie(start, end) = boost::adjacent_vertices(x0, g); start != end; start++) {
        Vertex y = *start;

        if (visited[y]) continue;
        visited[y] = true;

        if (mate[y] == g.null_vertex()) { // y is unmatched
            mate[y] = x0;
            mate[x0] = y;

            return true;
        }

        // y is matched with x1

        Vertex x1 = mate[y];

        bool path_found = find_path(x1, g, mate, visited);

        if (!path_found) {
            return false;
        }

        mate[y] = x0;
        mate[x0] = y;

        return true;
    }

    return false;

}



void match_greedy(const Graph& g, VertexVector& mate) {

    Vertex null_vertex = g.null_vertex();

    // set all mates to null vector
    std::fill(mate.begin(), mate.end(), null_vertex);

    // do greedy matching over all edges
    EdgeIterator start, end;
    for (std::tie(start, end) = boost::edges(g); start != end; start++) {
        Edge e = *start;
        Vertex u = boost::source(e, g);
        Vertex v = boost::target(e, g);
        if (mate[u] == null_vertex && mate[v] == null_vertex) {
            mate[u] = v;
            mate[v] = u;
        }
    }

}
