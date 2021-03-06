//
// Created by suem on 10/27/16.
//

#pragma once

#include <string>
#include "graphtypes.h"
#include <boost/graph/bipartite.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

inline void verify_bipartite(const Graph& g) {
    bool bp = boost::is_bipartite(g);
    if (!bp) throw "Graph not bipartite";
}

inline void verify_matching(const Graph& g, const VertexVector& mate, unsigned long max_matching_size) {
    unsigned long matching_size = boost::matching_size(g, &mate[0]);

    VertexIterator start, end;
    for (std::tie(start, end) = boost::vertices(g); start != end; start++) {
        Vertex v = *start;
        if (mate[v] != g.null_vertex() && !boost::edge(v, mate[v], g).second) {
            std::cout << v << ", " << mate[v] << std::endl;
            throw "matching incorrect, invalid edge";
        }
    }

    if (matching_size < max_matching_size) {
        std::cout <<  " Matching too small  " << std::endl;
		std::cout << matching_size  << std::endl;
        std::cout << " < "<< std::endl;
        std::cout << max_matching_size << std::endl;
		throw "matching not maximal";
    } else if (matching_size > max_matching_size) {
        std::cout <<  " Matching too large  " << std::endl;
		std::cout << matching_size << std::endl;
        std::cout << " > " << std::endl;
        std::cout << max_matching_size << std::endl;
        throw "matching incorrect, too large";
    }
}