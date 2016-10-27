//
// Created by suem on 10/25/16.
//
#pragma once

#include <iostream>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "graphtypes.h"

inline void loadGraph(Graph& graph) {
	unsigned int u, v;
	while (std::cin >> u && std::cin >> v) {
		boost::add_edge(u, v, graph);
	}
}
