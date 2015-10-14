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
 * Source code is located: /usr/local/termite/graphlab/apps/label_propagation.
 * Modified by zhigang wang, in order to be consistent with other systems.
 */

#include <string>
#include <fstream>
#include <graphlab.hpp>

struct vertex_data {
  int labelid;
  vertex_data() :
      labelid(0) {
  }

  void save(graphlab::oarchive& oarc) const {
    oarc << labelid;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> labelid;
  }
};

struct label_counter {
  std::map<int, int> label_count;

  label_counter() {
  }

  ~label_counter() {
	  label_count.clear();
	  std::map<int, int>().swap(label_count);
  }

  label_counter& operator+=(const label_counter& other) {
    for ( std::map<int, int>::const_iterator iter = other.label_count.begin();
              iter != other.label_count.end(); ++iter ) {
            label_count[iter->first] += iter->second;
    }

    return *this;
  }

  void save(graphlab::oarchive& oarc) const {
    oarc << label_count;
  }

  void load(graphlab::iarchive& iarc) {
    iarc >> label_count;
  }
};

size_t ITERATIONS = 0;

// The vertex data is its label
typedef vertex_data vertex_data_type;
typedef graphlab::empty edge_data_type;
typedef label_counter gather_data_type;

// The graph type is determined by the vertex and edge data types
typedef graphlab::distributed_graph<vertex_data_type, edge_data_type> graph_type;

//set label id at vertex id
void initialize_vertex(graph_type::vertex_type& v) {
  v.data().labelid = v.id();
}

class labelpropagation :
  public graphlab::ivertex_program<graph_type, gather_data_type> {

  public:
    edge_dir_type gather_edges(icontext_type& context, const vertex_type& vertex) const {
      //return graphlab::ALL_EDGES;
      return graphlab::IN_EDGES; //disk
    }

    gather_data_type gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
      // figure out which data to get from the edge.
      //bool isEdgeSource = (vertex.id() == edge.source().id());
      //vertex_data_type neighbor_label = isEdgeSource ? edge.target().data() : edge.source().data();
      vertex_data_type neighbor_label = edge.source().data();

      // make a label_counter and place the neighbor data in it
      gather_data_type counter;
      counter.label_count[neighbor_label.labelid] = 1;

      // gather_type is a label counter, so += will add neighbor counts to the
      // label_count map.
      return counter;
    }

    void apply(icontext_type& context, vertex_type& vertex, const gather_type& total) {
    	bool changed = true;

      int maxCount = 0;
      vertex_data_type maxLabel;

      // Figure out which label of the vertex's neighbors' labels is most common
      for ( std::map<int, int>::const_iterator iter = total.label_count.begin();
                iter != total.label_count.end(); ++iter ) {
              if (iter->second > maxCount) {
                maxCount = iter->second;
                maxLabel.labelid = iter->first;
              }
      }


      // if maxLabel differs to vertex data, mark vertex as changed and update
      // its data.
      if (vertex.data().labelid != maxLabel.labelid) {
        vertex.data().labelid = maxLabel.labelid;
      }

      context.setUpdateFlag(changed);
      context.signal(vertex);
    }

    edge_dir_type scatter_edges(icontext_type& context, const vertex_type& vertex) const {
      return graphlab::NO_EDGES;
      // if vertex data changes, scatter to all edges.
      /*if (changed) {
        return graphlab::OUT_EDGES;
      } else {
        return graphlab::NO_EDGES;
      }*/
      //return graphlab::OUT_EDGES; //disk
    }

    void scatter(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
      //bool isEdgeSource = (vertex.id() == edge.source().id());
      //context.signal(isEdgeSource ? edge.target() : edge.source());
      //context.signal(edge.target()); //disk
    }

    void save(graphlab::oarchive& oarc) const {
    }

    void load(graphlab::iarchive& iarc) {
    }
  };

struct labelpropagation_writer {
  std::string save_vertex(graph_type::vertex_type v) {
    std::stringstream strm;
    strm << v.id() << "\t" << v.data().labelid << "\n";
    return strm.str();
  }
  std::string save_edge (graph_type::edge_type e) { return ""; }
};


int main(int argc, char** argv) {
  std::cout << "Label Propagation\n";

  // Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  global_logger().set_log_level(LOG_INFO);

  // Parse command line options -----------------------------------------------
  graphlab::command_line_options clopts("Label Propagation algorithm.");
  std::string graph_dir;
  std::string saveprefix;
  std::string format = "adj";
  std::string execution_type = "synchronous";

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
  clopts.add_positional("useVerDisk");

  if(!clopts.parse(argc, argv)) {
    dc.cout() << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
  if (graph_dir == "") {
    dc.cout() << "Graph not specified. Cannot continue";
    return EXIT_FAILURE;
  }

  // Enable gather caching in the engine
  clopts.get_engine_args().set_option("use_cache", false);

  if (ITERATIONS) {
      // make sure this is the synchronous engine
      dc.cout() << "--iterations set. Forcing Synchronous engine, and running "
                << "for " << ITERATIONS << " iterations." << std::endl;
      clopts.get_engine_args().set_option("type", "synchronous");
      clopts.get_engine_args().set_option("max_iterations", ITERATIONS);
      clopts.get_engine_args().set_option("sched_allv", true);
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
  graphlab::omni_engine<labelpropagation> engine(dc, graph, execution_type, clopts);
  engine.signal_all();
  engine.start();
  const float runtime = engine.elapsed_seconds();
  dc.cout() << "Finished running engine in " << runtime << " seconds." << std::endl;

  if (saveprefix != "") {
    graph.save(saveprefix, labelpropagation_writer(),
    		false,  // do not gzip
    		true,   // save vertices
    		false,  // do not save edges
    		1);     // #_of_files per machine
  }

  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
}
