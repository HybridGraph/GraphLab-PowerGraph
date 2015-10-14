/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */

#include <vector>
#include <cstdlib>
#include <ctime>
#include <string>
#include <fstream>
#include <graphlab.hpp>


/**
 * \brief The type used to measure distances in the graph.
 */
typedef double distance_type;
size_t ITERATIONS;
bool USE_DELTA_CACHE = false;


/**
 * \brief The current distance of the vertex.
 * Moreover prior to each call to init the vertex program's previous state is cleared.
 * Therefore any persistent state must be saved into the vertex data.
 */
struct vertex_data : graphlab::IS_POD_TYPE {
  distance_type dist; // initialized by the maximum of FLOAT, disk.
  vertex_data(distance_type dist = std::numeric_limits<distance_type>::max()) :
    dist(dist) { }
}; // end of vertex data

/**
 * \brief This class is used as the gather type.
 */
struct min_distance_type : graphlab::IS_POD_TYPE {
  distance_type dist;
  min_distance_type(distance_type dist = 
                    std::numeric_limits<distance_type>::max()) : dist(dist) { }
  min_distance_type& operator+=(const min_distance_type& other) {
    dist = std::min(dist, other.dist);
    return *this;
  }
};

typedef vertex_data vertex_data_type;
typedef graphlab::empty edge_data_type;
typedef min_distance_type gather_data_type;

/**
 * \brief The graph type encodes the distances between vertices and
 * edges
 */
typedef graphlab::distributed_graph<vertex_data_type, edge_data_type> graph_type;

/**
 * \brief The single source shortest path vertex program.
 *
   *\li graph_type: the type of graph used to store the data for this
   * vertex program.  This currently always the distributed_graph.
   *
   * \li gather_type: the type used in the gather phase and must
   * implement the operator+= function.
   *
   * \li message_type: The type used for signaling and is typically
   * empty.  However if a message type is desired it must implement
   * the operator+= to allow message merging across the network.  In
   * addition the message type may also implement the priority()
   * function which returns a double assigning a priority to the
   * reception of the message (used by the asynchronous engines). We
   * provide a basic set of simple prioritized messages in
   * \ref graphlab::signals.
   *
   *
   * template<typename Graph,
           typename GatherType,
           typename MessageType = graphlab::empty>
  class ivertex_program {
 */
class sssp :
  public graphlab::ivertex_program<graph_type, gather_data_type> {
  bool changed;
  //distance_type lastchange;
public:

  edge_dir_type gather_edges(icontext_type& context, 
                             const vertex_type& vertex) const { 
    return graphlab::IN_EDGES;
  }

  gather_data_type gather(icontext_type& context, const vertex_type& vertex,
      		edge_type& edge) const {
	gather_data_type msg;
	msg.dist = edge.source().data().dist
			+ (std::rand()%10000)*0.0001; // the edge weight is 0-1, double, disk.
	/*std::cout << "vid=" << vertex.id()
			  << ", edge=(" << edge.source().id() << "," << edge.target().id() << ")\n";*/
    return msg;
  }

  /**
   * \brief If the distance is smaller then update
   */
  void apply(icontext_type& context, vertex_type& vertex,
             const gather_type& total) {
    changed = false;
    if(context.iteration() == 0) {
    	changed = true;
    	vertex.data().dist = 0;
    	context.setUpdateFlag(changed);
    	return;
    	//lastchange = vertex.data().dist;
    }
    if(vertex.data().dist > total.dist) {
      changed = true;
      vertex.data().dist = total.dist;
      //lastchange = vertex.data().dist;
    }

    context.setUpdateFlag(changed);
    //std::cout << "vid=" << vertex.id() << ", val=" << vertex.data().dist << "\n";
  }

  /**
   * \brief Determine if SSSP should run on all edges or just in edges
   */
  edge_dir_type scatter_edges(icontext_type& context, 
                             const vertex_type& vertex) const {
    if(changed) {
    	return graphlab::OUT_EDGES;
    } else {
    	return graphlab::NO_EDGES;
    }
  }

  /**
   * \brief The scatter function just signal adjacent pages 
   */
  void scatter(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const {
    const vertex_type other = edge.target();
    /*if (USE_DELTA_CACHE) {
    	gather_data_type delta;
    	delta.dist = lastchange;
    	context.post_delta(other, delta);
    }*/

    context.signal(other);
  }

  void save(graphlab::oarchive& oarc) const {
	  /*if (USE_DELTA_CACHE) {
		  oarc << lastchange;
	  }*/
	  oarc << changed; // otherwise, the value of changed cannot be available for different functions.
  }
  void load(graphlab::iarchive& iarc) {
	  /*if (USE_DELTA_CACHE) {
		  iarc >> lastchange;
	  }*/
	  iarc >> changed;
  }
}; // end of shortest path vertex program

/**
 * \brief We want to save the final graph so we define a write which will be
 * used in graph.save("path/prefix", pagerank_writer()) to save the graph.
 */
struct shortest_path_writer {
  std::string save_vertex(const graph_type::vertex_type& vtx) {
    std::stringstream strm;
    strm << vtx.id() << "\t" << vtx.data().dist << "\n";
    return strm.str();
  }
  std::string save_edge(graph_type::edge_type e) { return ""; }
}; // end of shortest_path_writer


int main(int argc, char** argv) {
  std::cout << "Single Source Shortest Path Algorithm\n";

  // Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  global_logger().set_log_level(LOG_INFO);

  // Parse command line options -----------------------------------------------
  graphlab::command_line_options 
    clopts("Single Source Shortest Path Algorithm.");
  std::string graph_dir;
  std::string saveprefix;
  std::string format = "adj";
  std::string exec_type = "synchronous";
  size_t powerlaw = 0;
  unsigned int source = std::numeric_limits<unsigned int>::max();

  //disk
  size_t ver_block_size = 1;
  size_t ver_buf_block_num = 1;
  size_t edge_block_size = 1;
  size_t useVerDisk = 0;
  size_t useEdgeDisk = 0;

  clopts.attach_option("graph", graph_dir,
                       "The graph file.  If none is provided "
                       "then a toy graph will be created");
  clopts.add_positional("graph");
  clopts.attach_option("format", format,
                       "graph format");
  clopts.attach_option("source", source,
                       "The source vertex");
  clopts.add_positional("source");

  clopts.attach_option("engine", exec_type, 
                       "The engine type synchronous or asynchronous");
  
  clopts.attach_option("powerlaw", powerlaw,
                       "Generate a synthetic powerlaw out-degree graph. ");
  clopts.attach_option("saveprefix", saveprefix,
                       "If set, will save the resultant pagerank to a "
                       "sequence of files with prefix saveprefix");
  clopts.attach_option("iterations", ITERATIONS,
                           "If set, will force the use of the synchronous engine"
                           "overriding any engine option set by the --engine parameter. "
                           "Runs complete (non-dynamic) label propagation for a fixed "
                           "number of iterations. Also overrides the iterations "
                           "option in the engine");
  /*clopts.attach_option("use_delta", USE_DELTA_CACHE,
                         "Use the delta cache to reduce time in gather.");*/

  clopts.attach_option("ver_block_size", ver_block_size,
  		  "the number of vertices per vertex block");
  clopts.add_positional("ver_block_size");

  clopts.attach_option("ver_buf_block_num", ver_buf_block_num,
		  "the number of blocks in vertex buffer");
  clopts.add_positional("ver_buf_block_num");

  clopts.attach_option("edge_block_size", edge_block_size,
		  "the size of each edge block, i.e. the number of source/target vertices");
  clopts.add_positional("edge_block_size");

  clopts.attach_option("useVerDisk", useVerDisk,
  		  "useVerDisk or not? =0:false, =1:true");
  clopts.add_positional("useVerDisk");

  clopts.attach_option("useEdgeDisk", useEdgeDisk,
  		  "useEdgeDisk or not? =0:false, =1:true");
  clopts.add_positional("useEdgeDisk");

  if(!clopts.parse(argc, argv)) {
    dc.cout() << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
  if(source == std::numeric_limits<unsigned int>::max()) {
  	dc.cout() << "Source vertex not specified. Cannot continue";
  	return EXIT_FAILURE;
  } else {
  	dc.cout() << "Using source=" << source << " to sssp\n" << std::endl;
  }

  // Enable gather caching in the engine
  /*if (USE_DELTA_CACHE) {
	  clopts.get_engine_args().set_option("use_cache", USE_DELTA_CACHE);
  }*/
  clopts.get_engine_args().set_option("use_cache", false);

  if (ITERATIONS) {
      // make sure this is the synchronous engine
      dc.cout() << "--iterations set. Forcing Synchronous engine, and running "
                << "for " << ITERATIONS << " iterations." << std::endl;
      clopts.get_engine_args().set_option("type", "synchronous");
      clopts.get_engine_args().set_option("max_iterations", ITERATIONS);
      clopts.get_engine_args().set_option("sched_allv", false);
  }

  clopts.setVerBlockSize(ver_block_size);
  clopts.setVerBufBlockNum(ver_buf_block_num);
  clopts.setEdgeBlockSize(edge_block_size);
  clopts.setUseVerDisk(useVerDisk);
  clopts.setUseEdgeDisk(useEdgeDisk);

  // Build the graph ----------------------------------------------------------
  graph_type graph(dc, clopts);
  if(powerlaw > 0) { // make a synthetic graph
    dc.cout() << "Loading synthetic Powerlaw graph." << std::endl;
    graph.load_synthetic_powerlaw(powerlaw, false, 2, 100000000);
  } else if (graph_dir.length() > 0) { // Load the graph from a file
    dc.cout() << "Loading graph in format: "<< format << std::endl;
    graph.load_format(graph_dir, format);
  } else {
    dc.cout() << "graph or powerlaw option must be specified" << std::endl;
    clopts.print_description();
    return EXIT_FAILURE;
  }
  // must call finalize before querying the graph
  graph.finalize();
  dc.cout() << "#vertices:  " << graph.num_vertices() << std::endl
            << "#edges:     " << graph.num_edges() << std::endl;

  // Running The Engine -------------------------------------------------------
  graphlab::omni_engine<sssp> engine(dc, graph, exec_type, clopts);

  engine.signal(source);
  engine.start();
  const float runtime = engine.elapsed_seconds();
  dc.cout() << "Finished Running engine in " << runtime
            << " seconds." << std::endl;

  // Save the final graph -----------------------------------------------------
  if (saveprefix != "") {
    graph.save(saveprefix, shortest_path_writer(),
    		   false,  // do not gzip
    		   true,   // save vertices
    		   false,  // do not save edges
    		   1);     // #_of_files per machine
  }

  // Tear-down communication layer and quit -----------------------------------
  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // End of main


// We render this entire program in the documentation



