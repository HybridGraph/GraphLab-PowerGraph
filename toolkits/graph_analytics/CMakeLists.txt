project(GraphProcessing)

add_graphlab_executable(pagerank pagerank.cpp)
add_graphlab_executable(lpa label_propagation.cpp)
add_graphlab_executable(sc semicluster.cpp)
add_graphlab_executable(sa simulate_advertisement.cpp)
add_graphlab_executable(sssp sssp.cpp)

# add_graphlab_executable(warp_pagerank warp_pagerank.cpp)
# add_graphlab_executable(warp_pagerank2 warp_pagerank2.cpp)
# add_graphlab_executable(warp_coloring warp_coloring.cpp)
# add_graphlab_executable(warp_bond_percolation warp_bond_percolation.cpp)
# add_graphlab_executable(warp_pagerank_sweeps warp_pagerank_sweeps.cpp)

#add_graphlab_executable(TSC TSC.cpp)
#requires_eigen(TSC)
