//
// Created by suem on 10/25/16.
//

#ifndef MAXIMUMMATCHING_GRAPHLOADER_H
#define MAXIMUMMATCHING_GRAPHLOADER_H

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS> Graph;

void loadGraph(Graph& graph);


#endif //MAXIMUMMATCHING_GRAPHLOADER_H
