# Disk-resident GraphLab-PowerGraph
A modified GraphLab PowerGraph. Vertices and edges reside on local disks.

##1. Introduction
The disk-resident GraphLab PowerGraph is based on the latest codebase of GraphLab PowerGraph (a distributed graph computation framework written in C++). The disk version supports accessing vertices and edges on disk. Specially, all edges are spilled onto disk since they are read sequentially. And for vertices, the disk version stores them in blocks and manages the blocks using the LRU replacing strategy. 

Note that the disk version cannot improve the data processing capacity, compared with the original GraphLab PowerGraph. This is because the modified version still employs the vertex-cut graph partitioning mechanism in the original system. Currently, the vertex-cut method works in memory, which limits the data processing capacity. We just use the disk-version to validate the I/O-inefficiency of existing pull-based techniques (GraphLab PowerGraph is a well-known memory-resident pull-based system).

##2. Building
The building, installation and tutorial of PowerLyra fully follow that of GraphLab PowerGraph.

##3. Running a job
###3.1 PageRank
_cd $GRAPHLAB_HOME/release/toolkits/graph_analytics_  
_hadoop dfs -rmr graphlab/output_  
_mpiexec -n 5 env CLASSPATH=\`hadoop classpath\` pagerank --format=adj --graph_opts="ingress=oblivious" --graph=hdfs:///user/root/graphlab/livej/ --saveprefix=hdfs:///user/root/graphlab/output/pagerank --iterations=5 --ver_block_size=10000 --ver_buf_block_num=$5 --edge_block_size=$6 --useVerDisk=1 --useEdgeDisk=1_

About input arguments:  
-n:  the number of MPI instances to run this job  
--graph:  the input directory on HDFS  
--saveprefix: the output directory on HDFS  
--iterations: the maximal number of iterations  
--ver_block_size: the number of vertices in each vertex block (used in disk version)  
--ver_buf_block_num:  the number of vertex blocks which are kept in memory  
--edge_block_size:  the number of vertices associated with edges in each edge block  
--useVerDisk:  1(storing vertices on disk):0(keeping all vertices in memory)  
--useEdgeDisk:  1(storing edges on disk):0(keeping all edges in memory)

###3.2 Single Source Shortest Path (SSSP)  
_cd $GRAPHLAB_HOME/release/toolkits/graph_analytics_  
_hadoop dfs -rmr graphlab/output_  
_mpiexec -n 5 env CLASSPATH=\`hadoop classpath\` sssp --format=adj --graph_opts="ingress=oblivious" --graph=hdfs:///user/root/graphlab/livej/ --saveprefix=hdfs:///user/root/graphlab/output/pagerank --iterations=100 --ver_block_size=10000 --ver_buf_block_num=$5 --edge_block_size=$6 --useVerDisk=1 --useEdgeDisk=1 --source=2_

About input arguments:  
Compared with PageRank, the additional argument is "--source", the source vertex id.
