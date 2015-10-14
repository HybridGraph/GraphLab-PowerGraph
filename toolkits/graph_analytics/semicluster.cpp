/**
 * SemiCluster Algorithm, written by zhigang wang.
 */

#include <string>
#include <vector>
#include <set>
#include <fstream>
#include <algorithm>
#include <graphlab.hpp>

#include <boost/range/iterator_range.hpp>

struct cluster {
	std::set<int> vertices; // set of vertices, no multiple elements
	double score; // score of current semi_cluster
	double innerScore; // inner score
	double boundaryScore; // boundary score

	cluster() {
		score = 1.0;
		innerScore = 0.0;
		boundaryScore = 0.0;
	}

	cluster(const cluster& c) {
		for (std::set<int>::const_iterator iter = c.vertices.begin();
				iter != c.vertices.end(); ++iter) {
			vertices.insert(*iter);
		}
		score = c.score;
		innerScore = c.innerScore;
		boundaryScore = c.boundaryScore;
	}

	~cluster() {
		vertices.clear();
		std::set<int>().swap(vertices);
	}

	void addVertex(const int sid,
			const std::map<int, double>& edges, const double scoreFactor) {
		if (vertices.insert(sid).second) {
			if (size() == 1) {
				for (std::map<int, double>::const_iterator iter = edges.begin();
						iter != edges.end(); ++iter) {
					boundaryScore += iter->second;
				}
				score = 0.0;
			} else {
				for (std::map<int, double>::const_iterator iter = edges.begin();
						iter != edges.end(); ++iter) {
					if (vertices.find(iter->first) != vertices.end()) {
						innerScore += iter->second;
						boundaryScore -= iter->second;
					} else {
						boundaryScore += iter->second;
					}
				}
				score = (innerScore-scoreFactor*boundaryScore)/(size()*(size()-1)/2);
			}
		} // insert successfully, i.e. new vertex id
	}

	int size() {
		return vertices.size();
	}

	void save(graphlab::oarchive& oarc) const {
	    oarc << vertices << score << innerScore << boundaryScore;
	}

	void load(graphlab::iarchive& iarc) {
	    iarc >> vertices >> score >> innerScore >> boundaryScore;
	}

	// override the < operator, used for less order
	bool operator<(const cluster& other) const {
		return (score<other.score);
	}

	// override the > operator, used for greater order
	bool operator>(const cluster& other) const {
		return (score>other.score);
	}
};

struct vertex_data {
	std::vector<cluster> cluster_set;

	vertex_data() {
	}

	/** only maintain the top-maxClusters clusters based on scores */
	void update(int maxClusters) {
		std::sort(cluster_set.begin(), cluster_set.end(), std::greater<cluster>()); // sort, greater
		int delNum = cluster_set.size() - maxClusters;
		if (delNum > 0) {
			for (int i = 0; i < delNum; ++i) {
				cluster_set.pop_back(); // delete the cluster with smaller scores
			}
		}
	}

	/** return the string of vertex data */
	std::string toString() {
		std::stringstream strm;
		for (std::vector<cluster>::const_iterator cIter = cluster_set.begin();
				cIter != cluster_set.end(); ++cIter) {
			cluster c = *cIter;
			strm << "  (";
			for (std::set<int>::const_iterator vIter = c.vertices.begin();
					vIter != c.vertices.end(); ++vIter) {
				strm << *vIter << ",";
			}
			strm << c.score << ")";
		}
		return strm.str();
	}

	void save(graphlab::oarchive& oarc) const {
		oarc << cluster_set;
	}
	void load(graphlab::iarchive& iarc) {
		iarc >> cluster_set;
	}
};

struct gather_data {
	std::map<int, double> edges;
	std::vector<vertex_data> cluster_sets;

	gather_data() {
	}

	~gather_data() {
		edges.clear();
		std::map<int, double>().swap(edges);
		cluster_sets.clear();
		std::vector<vertex_data>().swap(cluster_sets);
	}

	gather_data& operator+=(const gather_data& other) {
		for (std::map<int, double>::const_iterator iter = other.edges.begin();
				iter != other.edges.end(); ++iter) {
			edges[iter->first] = iter->second;
		}

		for (std::vector<vertex_data>::const_iterator iter = other.cluster_sets.begin();
				iter != other.cluster_sets.end(); ++iter) {
			cluster_sets.push_back(*iter);
		}

    return *this;
  }

  void save(graphlab::oarchive& oarc) const {
    oarc << edges << cluster_sets;
  }

  void load(graphlab::iarchive& iarc) {
    iarc >> edges >> cluster_sets;
  }
};

size_t ITERATIONS = 0;
int MaxClusters = 2;
int MaxCapacity = 2;
double ScoreFactor = 0.5;
double EdgeWeight = 0.01;

typedef vertex_data vertex_data_type;
typedef graphlab::empty edge_data_type;
typedef gather_data gather_data_type;

// The graph type is determined by the vertex and edge data types
typedef graphlab::distributed_graph<vertex_data_type, edge_data_type> graph_type;

void initialize_vertex(graph_type::vertex_type& v) {
	cluster c;
	std::map<int, double> edges;
	c.addVertex(v.id(), edges, ScoreFactor);
	v.data().cluster_set.push_back(c);
} // now, we don't have the edge info, just use the null edges

class semicluster :
	public graphlab::ivertex_program<graph_type, gather_data_type>,
	public graphlab::IS_POD_TYPE {

	public:
		edge_dir_type gather_edges(icontext_type& context,
				const vertex_type& vertex) const {
			return graphlab::ALL_EDGES;
		}

		gather_data_type gather(icontext_type& context,
				const vertex_type& vertex, edge_type& edge) const {
			// figure out which data to get from the edge.
			bool isEdgeSource = (vertex.id()==edge.source().id());
			gather_data_type result;
			if (isEdgeSource) {
				std::pair<int, double> p(edge.target().id(), EdgeWeight);
				result.edges.insert(p);
				result.cluster_sets.push_back(edge.target().data());
			} else {
				std::pair<int, double> p(edge.source().id(), EdgeWeight);
				result.edges.insert(p);
				result.cluster_sets.push_back(edge.source().data());
			}

			return result;
		}

		void apply(icontext_type& context, vertex_type& vertex,
				const gather_type& total) {
			bool changed = true;

			for (std::vector<vertex_data>::const_iterator iter = total.cluster_sets.begin();
					iter != total.cluster_sets.end(); ++iter) {
				vertex_data neighbor = *iter;
				int len = neighbor.cluster_set.size();
				for (int i = 0; i < len; ++i) {
					cluster c = neighbor.cluster_set[i];
					bool isContained = (c.vertices.find(vertex.id()) != c.vertices.end());
					if (!isContained && c.size()<MaxCapacity) {
						cluster newC(c);
						newC.addVertex(vertex.id(), total.edges, ScoreFactor);
						vertex.data().cluster_set.push_back(newC);
					} else if (isContained) {
						vertex.data().cluster_set.push_back(c);
					}
				}
			}
			vertex.data().update(MaxClusters);

			context.setUpdateFlag(changed);
			context.signal(vertex);
			//vertex.num_in_edges();
			//std::cout << vertex.id() << "\t" << vertex.data().toString() << "\n";
			//std::cout << vertex.in_edges() << std::endl; // interrupted, no implemented!
		}

		edge_dir_type scatter_edges(icontext_type& context,
				const vertex_type& vertex) const {
			//return graphlab::ALL_EDGES;
			return graphlab::NO_EDGES;
		}

		void scatter(icontext_type& context,
				const vertex_type& vertex, edge_type& edge) const {
			bool isEdgeSource = (vertex.id()==edge.source().id());
			context.signal(isEdgeSource ? edge.target() : edge.source());
		}
	};

struct semicluster_writer {
  std::string save_vertex(graph_type::vertex_type v) {
    std::stringstream strm;
    strm << v.id() << "\t" << v.data().toString() << "\n";
    return strm.str();
  }
  std::string save_edge (graph_type::edge_type e) { return ""; }
};


int main(int argc, char** argv) {
  std::cout << "SemiCluster\n";

  // Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  global_logger().set_log_level(LOG_INFO);

  // Parse command line options -----------------------------------------------
  graphlab::command_line_options clopts("SemiCluster algorithm.");
  std::string graph_dir;
  std::string saveprefix;
  std::string format = "adj";
  std::string execution_type = "synchronous";

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

  //! All input parameters cannot be used before .parse() is invoked.
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
  dc.cout() << "global #vertices: " << graph.num_vertices()
		    << "\nglobal #edges:" << graph.num_edges() << std::endl;

  // Run the engine
  graphlab::omni_engine<semicluster> engine(dc, graph, execution_type, clopts);
  engine.signal_all();
  engine.start();
  const float runtime = engine.elapsed_seconds();
  dc.cout() << "Finished running engine in " << runtime << " seconds." << std::endl;

  if (saveprefix != "") {
    graph.save(saveprefix, semicluster_writer(),
       false,  // do not gzip
       true,   // save vertices
       false,  // do not save edges
       1);     // #_of_files per machine
  }

  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
}
