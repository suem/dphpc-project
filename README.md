# DPHPC 2016 Project

*Isabelle Roesch, Samuel Ueltschi, Conradin Roffler, Thomas Meier*

## Parallel Maximum Cardinality Matching on Bipartite Graphs

We implement and optimize two existing parallel algorithms that solve the problem of maximum cardinality matching in bipartite graphs: Parallel Pothen Fan (PPF) and Parallel Tree Grafting (PTG). Our work is focused mainly on the analysis and possible optimizations of the PPF algorithm. We reason about the theoretical performance of PPF using the PRAM model and give a detailed description of the additional optimizations we implemented. We benchmark and test our implementations in a highly multi-threaded environment using a Xeon Phi multicore processor. The results confirm that our optimizations for PPF improve the original algorithm. 

