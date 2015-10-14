/*
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
 *
 */

/**
 * Simulate advertisement, coded by xxxx.
 */

#include <string>
#include <vector>
#include <fstream>
#include <graphlab.hpp>

struct vertex_data {
  int adverid;
  int advernum; // how many in-neighbors like this advertisement
  vertex_data() :
      adverid(0), advernum(0) {
  }

  void save(graphlab::oarchive& oarc) const {
    oarc << adverid << advernum;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> adverid >> advernum;
  }
};

struct gather_data {
  std::map<int, int> adver_count; // adverid, counter
  gather_data() {
  }

  ~gather_data() {
	  adver_count.clear();
	  std::map<int, int>().swap(adver_count);
  }

  gather_data& operator+=(const gather_data& other) {
    for ( std::map<int, int>::const_iterator iter = other.adver_count.begin();
              iter != other.adver_count.end(); ++iter ) {
    	adver_count[iter->first] += iter->second;
    }

    return *this;
  }

  void save(graphlab::oarchive& oarc) const {
	oarc << adver_count;
  }

  void load(graphlab::iarchive& iarc) {
	iarc >> adver_count;
  }
};

size_t ITERATIONS = 0;

typedef vertex_data vertex_data_type;
typedef graphlab::empty edge_data_type;
typedef gather_data gather_data_type;

// The graph type is determined by the vertex and edge data types
typedef graphlab::distributed_graph<vertex_data_type, edge_data_type> graph_type;

//set advertisement id at vertex id
void initialize_vertex(graph_type::vertex_type& v) {
  v.data().adverid = v.id();
  v.data().advernum = 0;
}

class simulate_advertisement :
  public graphlab::ivertex_program<graph_type, gather_data_type>,
  public graphlab::IS_POD_TYPE {
    bool changed;

  public:
    edge_dir_type gather_edges(icontext_type& context, const vertex_type& vertex) const {
      if (context.iteration() == 0) {
    	  return graphlab::NO_EDGES;
      } else {
    	  return graphlab::IN_EDGES;
      }
    }

    /**
     * IN_EDGES, target vertex id == vertex.id().
     */
    gather_data_type gather(icontext_type& context, const vertex_type& vertex,
    		edge_type& edge) const {
      vertex_data_type neighbor_adver = edge.source().data();
      /*std::cout << "vid=" << vertex.id() << ", in_edge: ("
    		    << edge.source().id() << "," << edge.target().id() << ")\n";*/

      // make a adver_counter and place the neighbor data in it
      gather_data counter;
      counter.adver_count[neighbor_adver.adverid] = 1;

      //  += will add neighbor counts to the adver_count map.
      return counter;
    }

    void apply(icontext_type& context, vertex_type& vertex, const gather_type& total) {
    	changed = false;
    	if (context.iteration() == 0) {
    		changed = true;
    		context.setUpdateFlag(changed);
    		return;
    	}

    	int maxCount = 0;
    	std::vector<int> ids;

    	// Figure out which advertisement of the vertex's neighbors' ads is most common
    	for (std::map<int, int>::const_iterator iter = total.adver_count.begin();
    			iter != total.adver_count.end(); ++iter) {
    		if (iter->second > maxCount) {
    			maxCount = iter->second;
    			ids.clear();
    			std::vector<int>().swap(ids);
    			ids.push_back(iter->first);
    		} else if (iter->second == maxCount) {
    			ids.push_back(iter->first);
    		}
    	}

    	vertex_data_type maxAdver;
    	int idx = 0;
    	if (ids.size() > 2) {
    		idx = (std::rand()%ids.size());
    	}
    	maxAdver.adverid = ids[idx];
    	maxAdver.advernum = maxCount;

    	// if maxAdver or maxCount varies, mark vertex as changed and update its data.
    	if (vertex.data().adverid != maxAdver.adverid
    			|| vertex.data().advernum < maxAdver.advernum) {
    		changed = true;
    		vertex.data().adverid = maxAdver.adverid;
    		vertex.data().advernum = maxCount;
    	}
    	context.setUpdateFlag(changed);
    }

    edge_dir_type scatter_edges(icontext_type& context, const vertex_type& vertex) const {
      // if vertex data changes, scatter to all out-edges.
      if (changed) {
        return graphlab::OUT_EDGES;
      } else {
        return graphlab::NO_EDGES;
      }
    }

    void scatter(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
      /*std::cout << "prepare to send new_msg from " << vertex.id()
    		    << " to " << edge.target().id() << "\n";*/

      context.signal(edge.target());
    }
  };

struct simulateadvertisement_writer {
  std::string save_vertex(graph_type::vertex_type v) {
    std::stringstream strm;
    strm << v.id() << "\t" << v.data().adverid << "," << v.data().advernum << "\n";
    return strm.str();
  }
  std::string save_edge (graph_type::edge_type e) { return ""; }
};


int main(int argc, char** argv) {
  std::cout << "Simulate Advertisement\n";

  // Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  global_logger().set_log_level(LOG_INFO);

  // Parse command line options -----------------------------------------------
  graphlab::command_line_options clopts("Simulate Advertisement algorithm.");
  std::string graph_dir;
  std::string saveprefix;
  std::string format = "adj";
  std::string execution_type = "synchronous";
  unsigned int source = std::numeric_limits<unsigned int>::max();

  //disk
  size_t ver_block_size = 1;
  size_t ver_buf_block_num = 1;
  size_t edge_block_size = 1;
  size_t useVerDisk = 0;
  size_t useEdgeDisk = 0;

  clopts.attach_option("graph", graph_dir, "The graph file. Required ");
  clopts.add_positional("graph");
  clopts.attach_option("format", format,
                         "The graph file format");
  clopts.attach_option("saveprefix", saveprefix,
                       "If set, will save the resultant label propagation to a "
                       "sequence of files with prefix saveprefix");
  clopts.attach_option("iterations", ITERATIONS,
                         "If set, will force the use of the synchronous engine"
                         "overriding any engine option set by the --engine parameter. "
                         "Runs complete (non-dynamic) label propagation for a fixed "
                         "number of iterations. Also overrides the iterations "
                         "option in the engine");
  clopts.attach_option("source", source,
                         "The source vertex");
  clopts.add_positional("source");

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
  if (graph_dir == "") {
    dc.cout() << "Graph not specified. Cannot continue";
    return EXIT_FAILURE;
  }
  if(source == std::numeric_limits<unsigned int>::max()) {
	dc.cout() << "Source vertex not specified. Cannot continue";
	return EXIT_FAILURE;
  } else {
	dc.cout() << "Using source=" << source << " to simulate advertisement\n" << std::endl;
  }

  // Enable gather caching in the engine
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
  dc.cout() << "Loading graph in format: "<< format << std::endl;
  graphlab::timer load_timer;
  graph.load_format(graph_dir, format);
  // must call finalize before querying the graph
  graph.finalize();
  graph.transform_vertices(initialize_vertex);
  dc.cout() << "Finished loading in " << load_timer.current_time() << std::endl;
  dc.cout() << "#vertices: " << graph.num_vertices()
		    << "\n#edges:" << graph.num_edges() << std::endl;

  // Run the engine
  graphlab::omni_engine<simulate_advertisement> engine(dc, graph, execution_type, clopts);
  engine.signal(source);
  engine.start();
  const float runtime = engine.elapsed_seconds();
  dc.cout() << "Finished running engine in " << runtime << " seconds." << std::endl;

  if (saveprefix != "") {
    graph.save(saveprefix, simulateadvertisement_writer(),
    		false,  // do not gzip
    		true,   // save vertices
    		false,  // do not save edges
    		1);     // #_of_files per machine
  }

  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
}
