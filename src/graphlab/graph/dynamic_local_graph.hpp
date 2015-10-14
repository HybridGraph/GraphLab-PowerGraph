#ifndef GRAPHLAB_DYNAMIC_LOCAL_GRAPH_HPP
#define GRAPHLAB_DYNAMIC_LOCAL_GRAPH_HPP

#include <cmath>
#include <string>
#include <sstream>
#include <list>
#include <vector>
#include <set>
#include <map>

#include <queue>
#include <algorithm>
#include <functional>
#include <fstream>

#include <boost/bind.hpp>
#include <boost/unordered_set.hpp>
#include <boost/type_traits.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/range/iterator_range.hpp>

#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/graph/local_edge_buffer.hpp>
#include <graphlab/util/random.hpp>
#include <graphlab/util/generics/shuffle.hpp>
#include <graphlab/util/generics/counting_sort.hpp>
#include <graphlab/util/generics/dynamic_csr_storage_disk.hpp>
#include <graphlab/parallel/atomic.hpp>

#include <graphlab/logger/logger.hpp>
#include <graphlab/logger/assertions.hpp>

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>

#include <graphlab/util/random.hpp>
#include <graphlab/macros_def.hpp>

#include <sys/stat.h>


namespace graphlab {
  template<typename VertexData, typename EdgeData>
  class dynamic_local_graph {
  public:

    /** The type of the vertex data stored in the local_graph. */
    typedef VertexData vertex_data_type;

    /** The type of the edge data stored in the local_graph. */
    typedef EdgeData edge_data_type;

    typedef graphlab::vertex_id_type vertex_id_type;
    typedef graphlab::edge_id_type edge_id_type;

    typedef std::vector<VertexData> vertex_block_unit;

  private:
    class edge_iterator;
    class edge_iterator_disk;

  public:
    typedef boost::iterator_range<edge_iterator> edge_list_type;
    typedef boost::iterator_range<edge_iterator_disk> edge_list_disk_type;

    /** Vertex object which provides access to the vertex data
     * and information about it.
     */
    class vertex_type;

      /** Edge object which provides access to the edge data
     * and information about it.
     */
    class edge_type;

  public:

    // CONSTRUCTORS ============================================================>
    /** Create an empty local_graph. */
    dynamic_local_graph()  {
    	vDataDir = "data/vertices/";
    	eDataDirCSR = "data/edges/csr/";
    	eDataDirCSC = "data/edges/csc/";

    	clcdir_Prefix = "rm -rf ";
    	mkdir_Prefix = "mkdir -p ";
    	system((clcdir_Prefix+vDataDir).c_str());
    	system((mkdir_Prefix+vDataDir).c_str());
    	system((clcdir_Prefix+eDataDirCSR).c_str());
    	system((mkdir_Prefix+eDataDirCSR).c_str());
    	system((clcdir_Prefix+eDataDirCSC).c_str());
    	system((mkdir_Prefix+eDataDirCSC).c_str());

    	vDataFile_Prefix = "block-";
    	eDataFile_Prefix = "block-";

    	vertex_block_size = 2; //each block includes 2 vertices
    	/**
    	 * vertex_buffer includes 2 vertex_block_unit,
    	 * total vertices in memory buffer = vertex_block_size*vertex_buffer_block_num.
    	 */
    	vertex_buffer_block_num = 2;
    	vertex_buffer_block_counter = 0;
    	edge_block_size = 2; //each block includes 2 vertices' edges
    	tmp_counter = 0;
    	total_block_num = 0;
    	total_vertex_num = 0;

    	useVerDisk = false;
    	useEdgeDisk = false; //edges are resident on disk?
    	closeLocalReadWriteModel();
    }

    /** Create a local_graph with nverts vertices. */
    dynamic_local_graph(size_t nverts) :
    	vertices(nverts)  {
    	vDataDir = "data/vertices/";
    	eDataDirCSR = "data/edges/csr/";
    	eDataDirCSC = "data/edges/csc/";

    	clcdir_Prefix = "rm -rf ";
    	mkdir_Prefix = "mkdir -p ";
    	system((clcdir_Prefix+vDataDir).c_str());
    	system((mkdir_Prefix+vDataDir).c_str());
    	system((clcdir_Prefix+eDataDirCSR).c_str());
    	system((mkdir_Prefix+eDataDirCSR).c_str());
    	system((clcdir_Prefix+eDataDirCSC).c_str());
    	system((mkdir_Prefix+eDataDirCSC).c_str());

    	vDataFile_Prefix = "block-";
    	eDataFile_Prefix = "block-";

    	vertex_block_size = 2; //each block includes 2 vertices
    	/**
    	 * vertex_buffer includes 2 vertex_block_unit,
    	 * total vertices in memory buffer = vertex_block_size*vertex_buffer_block_num.
    	 */
    	vertex_buffer_block_num = 2;
    	vertex_buffer_block_counter = 0;
    	edge_block_size = 2; //each block includes 2 vertices' edges
    	tmp_counter = 0;
    	total_block_num = 0;
    	total_vertex_num = 0;

    	useVerDisk = false;
    	useEdgeDisk = false;  //edges are resident on disk?
    	closeLocalReadWriteModel();
    }

    std::string getVertexDataFileName() {
    	return (vDataDir+vDataFile_Prefix);
    }

    std::string getVertexDataFileName(int block_id) {
    	std::stringstream ss;
    	ss << block_id;
    	return (vDataDir+vDataFile_Prefix+ss.str());
    }

    // METHODS =================================================================>

    static bool is_dynamic() {
      return true;
    }

    /**
     * \brief Resets the local_graph state.
     */
    void clear() {
    	clear_vertex_buffer();
    	edges.clear();
    	_csc_storage.clear();
    	_csr_storage.clear();
    	std::vector<vertex_block_unit>().swap(vertex_buffer);
    	std::vector<EdgeData>().swap(edges);
    	edge_buffer.clear();
    }

    /** \brief Get the number of vertices */
    size_t num_vertices() const {
      return (useVerDisk)? total_vertex_num:vertices.size();
    } // end of num vertices

    /** \brief Get the number of edges */
    size_t num_edges() const {
        return edges.size();
    } // end of num edges


    /**
     * \brief Creates a vertex containing the vertex data and returns the id
     * of the new vertex id. Vertex ids are assigned in increasing order with
     * the first vertex having id 0.
     */
    void add_vertex(lvid_type vid, const VertexData& vdata = VertexData()) {
    	if(vid >= vertices.size()) {
    		// Enable capacity doubling if resizing beyond capacity
    		if(vid >= vertices.capacity()) {
    			const size_t new_size = std::max(2 * vertices.capacity(),
    					size_t(vid));
    			vertices.reserve(new_size);
    		}
    		vertices.resize(vid+1);
    	}
    	vertices[vid] = vdata;
    	//logstream(LOG_INFO) << "add vertex, lvid=" << vid << std::endl;
    } // End of add vertex;

    void reserve(size_t num_vertices) {
      logstream(LOG_FATAL) << "this function has been deleted" << std::endl;
      //ASSERT_GE(num_vertices, vertices.size());
      //vertices.reserve(num_vertices);
    }

    /**
     * \brief Add additional vertices up to provided num_vertices.  This will
     * fail if resizing down.
     */
    void resize(size_t num_vertices) {
      ASSERT_GE(num_vertices, vertices.size());
      vertices.resize(num_vertices);
    } // End of resize

    void reserve_edge_space(size_t n) {
      edge_buffer.reserve_edge_space(n);
    }
    /**
     * \brief Creates an edge connecting vertex source to vertex target.  Any
     * existing data will be cleared. Should not be called after finalization.
     */
    edge_id_type add_edge(lvid_type source, lvid_type target,
                          const EdgeData& edata = EdgeData()) {
      if(source == target) {
        logstream(LOG_FATAL)
          << "Attempting to add self edge (" << source << " -> " << target <<  ").  "
          << "This operation is not permitted in GraphLab!" << std::endl;
        ASSERT_MSG(source != target, "Attempting to add self edge!");
      }

      //logstream(LOG_INFO) << "add edge, (" << source << "," << target << ")\n";
      if(source >= num_vertices() || target >= num_vertices())
        add_vertex(std::max(source, target));

      // Add the edge to the set of edge data (this copies the edata)
      edge_buffer.add_edge(source, target, edata);

      size_t tmp = std::max(source, target);
      tmp_counter = std::max(tmp, tmp_counter);

      // This is not the final edge_id, so we always return 0.
      return 0;
    } // End of add edge

    /**
     * \brief Add edges in block.
     */
    void add_edges(const std::vector<lvid_type>& src_arr,
                   const std::vector<lvid_type>& dst_arr,
                   const std::vector<EdgeData>& edata_arr) {
      ASSERT_TRUE((src_arr.size() == dst_arr.size())
                  && (src_arr.size() == edata_arr.size()));

      for (size_t i = 0; i < src_arr.size(); ++i) {
        lvid_type source = src_arr[i];
        lvid_type target = dst_arr[i];
        if ( source >= num_vertices()
             || target >= num_vertices() ) {
          logstream(LOG_FATAL)
            << "Attempting add_edge (" << source
            << " -> " << target
            << ") when there are only " << num_vertices()
            << " vertices" << std::endl;
          ASSERT_MSG(source < num_vertices(), "Invalid source vertex!");
          ASSERT_MSG(target < num_vertices(), "Invalid target vertex!");
        }

        if(source == target) {
          logstream(LOG_FATAL)
            << "Attempting to add self edge (" << source << " -> " << target <<  ").  "
            << "This operation is not permitted in GraphLab!" << std::endl;
          ASSERT_MSG(source != target, "Attempting to add self edge!");
        }
      }
      edge_buffer.add_block_edges(src_arr, dst_arr, edata_arr);
    } // End of add block edges


    /** \brief Returns a vertex of given ID. */
    vertex_type vertex(lvid_type vid) {
      ASSERT_LT(vid, num_vertices());
      return vertex_type(*this, vid);
    }

    /** \brief Returns a vertex of given ID. */
    const vertex_type vertex(lvid_type vid) const {
    	ASSERT_LT(vid, num_vertices());
    	return vertex_type(*this, vid);
    }

    /** \brief Returns a reference to the data stored on the vertex v. */
    VertexData& vertex_data(lvid_type v) {
    	ASSERT_LT(v, num_vertices());
    	if (useVerDisk) {
    		int block_id = getBlockId((int)v);
    		if (!isMemBlock(block_id)) {
    			if (isOverFlow()) {
    				int remove_block_id = getLRUBid();
    				if (isReadWriteModel()) {
    					save_vertex_block(remove_block_id);
    				} else {
    					clear_vertex_block(remove_block_id);
    				}
    				vertex_buffer_block_counter--;
    			}

    			load_vertex_block(block_id);
    			vertex_buffer_block_counter++;
    		} else {
    			incHits(block_id);
    		}

    		int localIdxInBlock = getLocalIdxInBlock((int)v);
    		return vertex_buffer[block_id][localIdxInBlock];
    	} else {
    		return vertices[v];
    	}
    } // end of data(v)

    /** \brief Returns a constant reference to the data stored on the vertex v. */
    const VertexData& vertex_data(lvid_type v) const {
    	ASSERT_LT(v, num_vertices());
    	if (useVerDisk) {
    		int block_id = getBlockId((int)v);
    		if (!isMemBlock(block_id)) {
    			if (isOverFlow()) {
    				int remove_block_id = getLRUBid();
    				if (isReadWriteModel()) {
    					save_vertex_block(remove_block_id);
    				} else {
    					clear_vertex_block(remove_block_id);
    				}
    				decVerBufBlockCounter();
    			}

    			load_vertex_block(block_id);
    			incVerBufBlockCounter();
    		} else {
    			incHits(block_id);
    		}

    		int localIdxInBlock = getLocalIdxInBlock((int)v);
    		return vertex_buffer[block_id][localIdxInBlock];
    	} else {
    		return vertices[v];
    	}
    } // end of data(v)

    /**
     * \brief Finalize the local_graph data structure by
     * sorting edges to maximize the efficiency of graphlab.
     * This function takes O(|V|log(degree)) time and will
     * fail if there are any duplicate edges.
     * Detail implementation depends on the type of graph_storage.
     * This is also automatically invoked by the engine at start.
     */
    void finalize() {

      graphlab::timer mytimer; mytimer.start();
#ifdef DEBUG_GRAPH
      logstream(LOG_DEBUG) << "Graph2 finalize starts." << std::endl;
#endif
      std::vector<edge_id_type> src_permute;
      std::vector<edge_id_type> dest_permute;
      std::vector<edge_id_type> src_counting_prefix_sum;
      std::vector<edge_id_type> dest_counting_prefix_sum;
      tmp_counter++;

#ifdef DEBUG_GRAPH
      logstream(LOG_DEBUG) << "Graph2 finalize: Sort by source vertex" << std::endl;
#endif
      counting_sort(edge_buffer.source_arr, dest_permute, &src_counting_prefix_sum);
      std::vector< std::pair<lvid_type, edge_id_type> >  csr_values;
      csr_values.reserve(dest_permute.size());
      edge_id_type begineid = edges.size();
      for (size_t i = 0; i < dest_permute.size(); ++i) {
    	  csr_values.push_back(std::pair<lvid_type, edge_id_type> (edge_buffer.target_arr[dest_permute[i]],
    			  begineid + dest_permute[i]));
      }

#ifdef DEBUG_GRAPH
      logstream(LOG_DEBUG) << "Graph2 finalize: Sort by dest id" << std::endl;
#endif
      counting_sort(edge_buffer.target_arr, src_permute, &dest_counting_prefix_sum);
      std::vector< std::pair<lvid_type, edge_id_type> >  csc_values;
      csc_values.reserve(src_permute.size());
      for (size_t i = 0; i < src_permute.size(); ++i) {
        csc_values.push_back(std::pair<lvid_type, edge_id_type> (edge_buffer.source_arr[src_permute[i]],
                                                                 begineid + src_permute[i]));
      }

      ASSERT_EQ(csc_values.size(), csr_values.size());

      // fast path with first time insertion.
      // note: now, we do not support inserting edges during computing. disk
      if (edges.size() == 0) {
    	  // this is invoked when loading and partitioning the original graph data.
        edges.swap(edge_buffer.data);
        edge_buffer.clear();

        _csr_storage.wrap(src_counting_prefix_sum, csr_values);
        dest_permute.clear();
        std::vector<edge_id_type>().swap(dest_permute);

        _csc_storage.wrap(dest_counting_prefix_sum, csc_values);
        src_permute.clear();
        std::vector<edge_id_type>().swap(src_permute);
      } else {
    	  //inserting new edges
    	logstream(LOG_FATAL) << "this function has been deleted!" << std::endl;
      }
      //ASSERT_EQ(_csr_storage.num_values(), _csc_storage.num_values());
      //ASSERT_EQ(_csr_storage.num_values(), edges.size());

      //

#ifdef DEBUG_GRAPH
      logstream(LOG_DEBUG) << "End of finalize." << std::endl;
#endif
      logstream(LOG_INFO) << "Graph finalized in " << mytimer.current_time()
                          << " secs" << std::endl;

#ifdef DEBUG_GRAPH
      _csr_storage.meminfo(std::cerr);
      _csc_storage.meminfo(std::cerr);
#endif
    } // End of finalize

    // return the location in one block.
    int getLocalIdxInBlock(int lvid) {
    	return (lvid%vertex_block_size);
    }

    // return the block_id this vertex belongs to
    int getBlockId(int lvid) {
    	return (lvid/vertex_block_size);
    }

    bool isMemBlock(int block_id) {
    	return bid2InMem[block_id];
    }

    void incHits(int bid) {
    	bid2Hits[bid]++;
    }

    void incVerBufBlockCounter() {
    	vertex_buffer_block_counter++;
    }

    void decVerBufBlockCounter() {
    	vertex_buffer_block_counter--;
    }

    bool isOverFlow() {
    	return (vertex_buffer_block_counter >= vertex_buffer_block_num);
    }

    int getLRUBid() {
    	long count = std::numeric_limits<long>::max();
    	int bid = -1;
    	for (int i = 0; i < (int)bid2InMem.size(); i++) {
    		if (bid2InMem[i] && count>bid2Hits[i]) {
    			count = bid2Hits[i];
    			bid = i;
    		}
    	}
    	return bid;
    }

    void load_vertex_block(int block_id) {
    	std::string filename = getVertexDataFileName(block_id);
    	//logstream(LOG_INFO) << "open vertex block file=" << filename << std::endl;
    	std::ifstream fin(filename.c_str());
    	iarchive iarc(fin);
    	ASSERT_TRUE(fin.good());
    	int loopLen = 0;
    	iarc >> loopLen;

    	vertex_block_unit vblock;
    	vblock.resize(loopLen);
    	for (int i = 0; i < loopLen; i++) {
    		iarc >> vblock[i];
    	}
    	vertex_buffer[block_id] = vblock;
    	bid2InMem[block_id] = true;
    	bid2Hits[block_id] = 0;

    	vblock.clear();
    	vertex_block_unit().swap(vblock);

    	fin.close();
    	if (useVerDisk) {
    		struct stat info;
    		stat(filename.c_str(), &info);
    		read_bytes = read_bytes + info.st_size;
    	}

    	/*logstream(LOG_INFO) << "disk debug. load " << loopLen << " vertices from block-"
    	    			<< block_id << std::endl;*/
    }

    /**
     * disk
     * Block file structure
     * num_of_vertices vertices.....
     */
    void save_vertex_block(int block_id) {
    	std::string filename = getVertexDataFileName(block_id);
    	remove(filename.c_str());
    	std::ofstream fout(filename.c_str(), std::ios::binary);
    	oarchive oarc(fout);
    	ASSERT_TRUE(fout.good());
    	int loopLen = vertex_buffer[block_id].size();

    	oarc << loopLen;
    	for (int i = 0; i < loopLen; i++) {
    		oarc << vertex_buffer[block_id][i];
    	}

    	fout.close();
    	//logstream(LOG_INFO) << "disk debug. save " << loopLen << " vertices to block-" << block_id << std::endl;
    	if (useVerDisk) {
    		struct stat info;
    		stat(filename.c_str(), &info);
    		write_bytes = write_bytes + info.st_size;
    	}
    	clear_vertex_block(block_id);
    }

    void clear_vertex_block(int block_id) {
    	vertex_buffer[block_id].clear();
    	vertex_block_unit().swap(vertex_buffer[block_id]);

    	bid2InMem[block_id] = false;
    	bid2Hits[block_id] = 0;
    }

    void clear_vertex_buffer() {
    	for (int i = 0; i < (int)bid2InMem.size(); i++) {
    		if (bid2InMem[i]) {
    			clear_vertex_block(i);
    			vertex_buffer_block_counter--;
    		}
    	}
    }

    //disk
    void save_vertices(bool _useVerDisk) {
    	useVerDisk = _useVerDisk;
    	if (!useVerDisk) {
    		return;
    	}

    	total_vertex_num = vertices.size();
    	int localIdxInBlock = 0;
    	vertex_block_unit vblock;
    	for (size_t lvid = 0; lvid < total_vertex_num; lvid++) {
    		localIdxInBlock = getLocalIdxInBlock((int)lvid); //offset/location in this block
    		if (localIdxInBlock >= (int)vblock.size()) {
    			// Enable capacity doubling if resizing beyond capacity
    			if (localIdxInBlock >= (int)vblock.capacity()) {
    				size_t new_size_tmp = std::max(2 * vblock.capacity(),
    						(size_t)localIdxInBlock);
    				size_t new_size_final = std::min(new_size_tmp, (size_t)vertex_block_size);
    				vblock.reserve(new_size_final);
    			}
    			vblock.resize(localIdxInBlock+1);
    		}
    		vblock[localIdxInBlock] = vertices[lvid];

    		if ((int)vblock.size() == vertex_block_size) {
    			int block_id = total_block_num;
    			if (block_id >= (int)vertex_buffer.size()) {
    				if (block_id >= (int)vertex_buffer.capacity()) {
    					size_t new_size = std::max(2 * vertex_buffer.capacity(),
    							(size_t)block_id);
    					vertex_buffer.reserve(new_size);
    					bid2InMem.reserve(new_size);
    					bid2Hits.reserve(new_size);
    				}
    				vertex_buffer.resize(block_id+1);
    				bid2InMem.resize(block_id+1);
    				bid2Hits.resize(block_id+1);
    			}

    			//buffer is not overflow, put it into buffer, otherwise, write it onto disk.
    			/*if (vertex_buffer_block_counter < vertex_buffer_block_num) {
    				vertex_buffer[block_id] = vblock;
    				bid2InMem[block_id] = true;
    				bid2Hits[block_id] = 0;
    				vertex_buffer_block_counter++;
    			} else {
    				vertex_buffer[block_id] = vblock;
    				save_vertex_block(block_id);
    			}*/

    			vertex_buffer[block_id] = vblock;
    			save_vertex_block(block_id);

    			vblock.clear();
    			vertex_block_unit().swap(vblock);
    			total_block_num++;
    		}
    	}
    	vertices.clear();
    	std::vector<VertexData>().swap(vertices);

    	// process the last block, i.e., current buffer
    	if (vblock.size() > 0) {
    		int block_id = total_block_num;
    		if (block_id >= (int)vertex_buffer.size()) {
    			if (block_id >= (int)vertex_buffer.capacity()) {
    				size_t new_size = std::max(2 * vertex_buffer.capacity(),
    						(size_t)block_id);
    				vertex_buffer.reserve(new_size);
    				bid2InMem.reserve(new_size);
    				bid2Hits.reserve(new_size);
    			}
    			vertex_buffer.resize(block_id+1);
    			bid2InMem.resize(block_id+1);
    			bid2Hits.resize(block_id+1);
    		}

    		//buffer is not overflow, put it into buffer, otherwise, write it onto disk.
    		/*if (vertex_buffer_block_counter < vertex_buffer_block_num) {
    			vertex_buffer[block_id] = vblock;
    			bid2InMem[block_id] = true;
    			bid2Hits[block_id] = 0;
    			vertex_buffer_block_counter++;
    		} else {
    			vertex_buffer[block_id] = vblock;
    			save_vertex_block(block_id);
    		}*/
    		vertex_buffer[block_id] = vblock;
    		save_vertex_block(block_id);

    		vblock.clear();
    		vertex_block_unit().swap(vblock);
    		total_block_num++;
    	}

    	logstream(LOG_INFO) << "save " << total_vertex_num
    			<< " vertices to " << total_block_num << " blocks" << std::endl;
    }

    void save_edges(bool _useEdgeDisk) {
    	useEdgeDisk = _useEdgeDisk;
    	if (useEdgeDisk) {
    		_csr_storage.init(eDataDirCSR+eDataFile_Prefix, edge_block_size, (int)num_vertices());
    		save_csr();
    		_csr_storage.clear();

    		_csc_storage.init(eDataDirCSC+eDataFile_Prefix, edge_block_size, (int)num_vertices());
    		save_csc();
    		_csc_storage.clear();
    	}
    }

    //disk
    void save_csr() {
    	if (!useEdgeDisk) {
    		return;
    	}

      int vid, len, eid, eCounter = 0, current_block_id = 0;
      int loopLen = (int)num_vertices();

      std::string filename = _csr_storage.getFileName(current_block_id);
      remove(filename.c_str());
      std::ofstream fout(filename.c_str(), std::ios::binary);
      ASSERT_TRUE(fout.good());
      for (vid = 0; vid < loopLen; ++vid) {
    	  int bid = _csr_storage.getBlockId(vid);
    	  if (current_block_id != bid) {
    		  fout.close();
    		  current_block_id = bid;
    		  filename = _csr_storage.getFileName(current_block_id);
    		  remove(filename.c_str());
    		  fout.open(filename.c_str(), std::ios::binary);
    		  ASSERT_TRUE(fout.good());
    	  }

    	  edge_list_type ls = out_edges((lvid_type)vid);
    	  len = (int)ls.size();
    	  fout.write((char*)(&vid), sizeof(vid));
    	  fout.write((char*)(&len), sizeof(len));
    	  _csr_storage.set_num_edges(vid, len);

    	  foreach(edge_type e, ls) {
    		  eid = (int)(e.target().id());
    		  fout.write((char*)(&eid), sizeof(eid));
    		  eCounter++;
    	  }
      }

      if (fout.is_open()) {
    	  fout.close();
      }
      logstream(LOG_INFO) << "save " << eCounter << " edges of "
    		  << loopLen  << " vertices to " << (current_block_id+1) << " blocks\n";
    }

    void save_csc() {
    	if (!useEdgeDisk) {
    		return;
    	}

    	int vid, len, eid, eCounter = 0, current_block_id = 0;
    	int loopLen = (int)num_vertices();

    	std::string filename = _csc_storage.getFileName(current_block_id);
    	remove(filename.c_str());
    	std::ofstream fout(filename.c_str(), std::ios::binary);
    	ASSERT_TRUE(fout.good());
    	for (vid = 0; vid < loopLen; ++vid) {
    		int bid = _csc_storage.getBlockId(vid);
    		if (current_block_id != bid) {
    			fout.close();
    			current_block_id = bid;
    			filename = _csc_storage.getFileName(current_block_id);
    			remove(filename.c_str());
    			fout.open(filename.c_str(), std::ios::binary);
    			ASSERT_TRUE(fout.good());
    		}

    		edge_list_type ls = in_edges((lvid_type)vid);
            len = (int)ls.size();
            fout.write((char*)(&vid), sizeof(vid));
            fout.write((char*)(&len), sizeof(len));
            _csc_storage.set_num_edges(vid, len);

            foreach(edge_type e, ls) {
            	eid = (int)(e.source().id());
            	fout.write((char*)(&eid), sizeof(eid));
            	eCounter++;
            }
    	}

    	if (fout.is_open()) {
    		fout.close();
    	}
    	logstream(LOG_INFO) << "save " << eCounter << " edges of "
    			<< loopLen  << " vertices to " << (current_block_id+1) << " blocks\n";
    }

    //disk.
    //note: _userVerDisk and _useEdgeDisk are not used now!!! we initialize them at save_vertices and save_edges.
    //Because, when loading graph data, we should guarantee that all data are resident in memory.
    void initLocalDiskArgs(bool _useVerDisk, bool _useEdgeDisk, int _ver_block_size,
    		int _ver_buf_block_num, int _edge_block_size) {
    	vertex_block_size = _ver_block_size;
    	vertex_buffer_block_num = _ver_buf_block_num;
    	edge_block_size = _edge_block_size; //each block includes 2 vertices' edges

    	vertex_buffer.resize(vertex_buffer_block_num);
    	bid2InMem.resize(vertex_buffer_block_num);
    	bid2Hits.resize(vertex_buffer_block_num);
    }

    //disk.
    void resetLocalBefIte() {
    	read_write = false;
    	read_bytes = 0;
    	write_bytes = 0;
    	if (useEdgeDisk) {
    		_csr_storage.resetLocalBefIte();
    		_csc_storage.resetLocalBefIte();
    	}
    }

    //disk.
    long getLocalIOBytes_Read() {
    	return (read_bytes + _csr_storage.getReadBytes() + _csc_storage.getReadBytes());
    }

    //disk.
    long getLocalIOBytes_Write() {
    	return write_bytes;
    }

    /** If the current block is in memory, use it, and save it if it should be replaced another block. */
    void openLocalReadWriteModel() {
    	read_write = true;
    }

    void closeLocalReadWriteModel() {
    	if (total_block_num > vertex_buffer_block_num) {
    		for (int i = 0; i < (int)bid2InMem.size(); i++) {
    			if (bid2InMem[i]) {
    				save_vertex_block(i);
    				vertex_buffer_block_counter--;
    			}
    		}
    	}

    	read_write = false;
    }

    bool isReadWriteModel() {
    	return read_write;
    }

    //disk.
    void clearLocalTmpFiles() {
    	if (useEdgeDisk) {
    		system((clcdir_Prefix+vDataDir).c_str());
    		system((clcdir_Prefix+eDataDirCSR).c_str());
    		system((clcdir_Prefix+eDataDirCSC).c_str());
    	}
    }

    /** \brief Load the local_graph from an archive */
    void load(iarchive& arc) {
      clear();
      // read the vertices
      arc >> vertex_buffer
          >> edges
          >> _csr_storage
          >> _csc_storage;
    } // end of load

    /** \brief Save the local_graph to an archive */
    void save(oarchive& arc) const {
      // Write the number of edges and vertices
      arc << vertex_buffer
          << edges
          << _csr_storage
          << _csc_storage;
    } // end of save

    /** swap two graphs */
    void swap(dynamic_local_graph& other) {
    	std::swap(vertices, other.vertices);
      std::swap(vertex_buffer, other.vertex_buffer);
      std::swap(edges, other.edges);
      std::swap(_csr_storage, other._csr_storage);
      std::swap(_csc_storage, other._csc_storage);
    } // end of swap


    /** \brief Load the local_graph from a file */
    void load(const std::string& filename) {
      std::ifstream fin(filename.c_str());
      iarchive iarc(fin);
      iarc >> *this;
      fin.close();
    } // end of load

    /**
     * \brief save the local_graph to the file given by the filename
     */
    void save(const std::string& filename) const {
      std::ofstream fout(filename.c_str());
      oarchive oarc(fout);
      oarc << *this;
      fout.close();
    } // end of save

    /**
     * \brief save the adjacency structure to a text file.
     *
     * Save the adjacency structure as a text file in:
     *    src_Id, dest_Id \n
     *    src_Id, dest_Id \n
     * format.
     */
    void save_adjacency(const std::string& filename) const {
    	remove(filename.c_str());
      std::ofstream fout(filename.c_str());
      ASSERT_TRUE(fout.good());

      for (size_t i = 0; i < num_vertices(); ++i) {
        vertex_type v(i);
        edge_list_type ls = v.out_edges();
        foreach(edge_type e, ls) {
          fout << (lvid_type)i << ", " << e.target().id() << "\n";
          ASSERT_TRUE(fout.good());
        }
      }
      fout.close();
    }

 /****************************************************************************
 *                       Internal Functions                                 *
 *                     ----------------------                               *
 * These functions functions and types provide internal access to the       *
 * underlying local_graph representation. They should not be used unless you      *
 * *really* know what you are doing.                                        *
 ****************************************************************************/
    /**
     * \internal
     * \brief Returns the number of in edges of the vertex with the given id. */
    size_t num_in_edges(const lvid_type v) const {
    	//logstream(LOG_INFO) << "request num_in_edges, lvid=" << v << std::endl;
    	return useEdgeDisk? _csc_storage.num_edges(v):_csc_storage.begin(v).pdistance_to(_csc_storage.end(v));
    }

    /**
     * \internal
     * \brief Returns the number of in edges of the vertex with the given id. */
    size_t num_out_edges(const lvid_type v) const {
    	//logstream(LOG_INFO) << "request num_out_edges, lvid=" << v << std::endl;
    	return useEdgeDisk? _csr_storage.num_edges(v):_csr_storage.begin(v).pdistance_to(_csr_storage.end(v));
    }

    /**
     * \internal
     * \brief Returns a list of in edges of the vertex with the given id. */
    edge_list_type in_edges(lvid_type v) {
      edge_iterator begin = edge_iterator(*this, edge_iterator::CSC,
                                          _csc_storage.begin(v), v);
      edge_iterator end = edge_iterator(*this, edge_iterator::CSC,
                                        _csc_storage.end(v), v);
      return boost::make_iterator_range(begin, end);
    }

    //disk, only used for disk version in iteration computation.
    edge_list_disk_type in_edges_disk(lvid_type v) {
    	edge_iterator_disk begin = edge_iterator_disk(*this, edge_iterator_disk::CSC,
    			_csc_storage.begin_disk(v), v);
    	edge_iterator_disk end = edge_iterator_disk(*this, edge_iterator_disk::CSC,
    			_csc_storage.end_disk(v), v);
    	return boost::make_iterator_range(begin, end);
    }

    //disk, only used for disk version in iteration computation.
    edge_list_disk_type out_edges_disk(lvid_type v) {
    	edge_iterator_disk begin = edge_iterator_disk(*this, edge_iterator_disk::CSR,
    			_csr_storage.begin_disk(v), v);
    	edge_iterator_disk end = edge_iterator_disk(*this, edge_iterator_disk::CSR,
    			_csr_storage.end_disk(v), v);
    	return boost::make_iterator_range(begin, end);
    }

    /**
     * \internal
     * \brief Returns a list of out edges of the vertex with the given id. */
    edge_list_type out_edges(lvid_type v) {
      edge_iterator begin = edge_iterator(*this, edge_iterator::CSR,
                                          _csr_storage.begin(v), v);
      edge_iterator end = edge_iterator(*this, edge_iterator::CSR,
                                        _csr_storage.end(v), v);
      return boost::make_iterator_range(begin, end);
    }

    /**
     * \internal
     * \brief Returns edge data of edge_type e
     * */
    EdgeData& edge_data(edge_id_type eid) {
      ASSERT_LT(eid, num_edges());
      return edges[eid];
    }
    /**
     * \internal
     * \brief Returns const edge data of edge_type e
     * */
    const EdgeData& edge_data(edge_id_type eid) const {
      ASSERT_LT(eid, num_edges());
      return edges[eid];
    }

    /**
     * \internal
     * \brief Returns the estimated memory footprint of the local_graph. */
    size_t estimate_sizeof() const {
      size_t vlist_size = sizeof(vertex_buffer) +
        sizeof(VertexData) * vertex_buffer.capacity()
        + sizeof(vertices) + sizeof(VertexData)*vertices.capacity();
      size_t elist_size = _csr_storage.estimate_sizeof()
          + _csc_storage.estimate_sizeof()
          + sizeof(edges) + sizeof(EdgeData)*edges.capacity();
      size_t ebuffer_size = edge_buffer.estimate_sizeof();
      return vlist_size + elist_size + ebuffer_size;
    }

    /** \internal
     * \brief For debug purpose, returns the largest vertex id in the edge_buffer
     */
    const lvid_type maxlvid() const {
      if (edge_buffer.size()) {
        lvid_type max(0);
        foreach(lvid_type i, edge_buffer.source_arr)
         max = std::max(max, i);
        foreach(lvid_type i, edge_buffer.target_arr)
         max = std::max(max, i);
        return max;
      } else {
        return lvid_type(-1);
      }
    }

  private:
    /**
     * \internal
     * CSR/CSC storage types
     */
    typedef dynamic_csr_storage_disk<std::pair<lvid_type, edge_id_type>, edge_id_type> csr_type;

    typedef typename csr_type::iterator csr_edge_iterator;

    // PRIVATE DATA MEMBERS ===================================================>
    //
    /** Disk file path and name */
    std::string vDataDir;
    std::string eDataDirCSR;
    std::string eDataDirCSC;
    std::string vDataFile_Prefix;
    std::string eDataFile_Prefix;
    std::string clcdir_Prefix;
    std::string mkdir_Prefix;
    bool read_write;

    /** The vertex buffer with a user-specified size */
    std::vector<VertexData> vertices;
    std::vector<vertex_block_unit> vertex_buffer;
    bool useVerDisk;
    bool useEdgeDisk;

    int vertex_block_size;
    int vertex_buffer_block_num;
    int vertex_buffer_block_counter;
    std::vector<bool> bid2InMem; //[i]=true: block_i is in memory.
    std::vector<long> bid2Hits; //[i]=20: hit block_i 20 times.

    int total_block_num;
    size_t total_vertex_num;
    int edge_block_size;
    size_t tmp_counter;

    long read_bytes;
    long write_bytes;

    /** Stores the edge data and edge relationships. */
    csr_type _csr_storage;
    csr_type _csc_storage;
    std::vector<EdgeData> edges;

    /** The edge data is a vector of edges where each edge stores its
        source, destination, and data. Used for temporary storage. The
        data is transferred into CSR+CSC representation in
        Finalize. This will be cleared after finalized.*/
    local_edge_buffer<VertexData, EdgeData> edge_buffer;

    /**************************************************************************/
    /*                                                                        */
    /*                            declare friends                             */
    /*                                                                        */
    /**************************************************************************/
    friend class local_graph_test;
  }; // End of class dynamic_local_graph


  template<typename VertexData, typename EdgeData>
  std::ostream& operator<<(std::ostream& out,
                           const dynamic_local_graph<VertexData, EdgeData>& local_graph) {
    for(lvid_type vid = 0; vid < local_graph.num_vertices(); ++vid) {
      foreach(edge_id_type eid, local_graph.out_edge_ids(vid))
        out << vid << ", " << local_graph.target(eid) << '\n';
    }
    return out;
  }
} // end of namespace graphlab


/////////////////////// Implementation of Helper Class ////////////////////////////

namespace graphlab {
  template<typename VertexData, typename EdgeData>
  class dynamic_local_graph<VertexData, EdgeData>::vertex_type {
     public:
       vertex_type(dynamic_local_graph& lgraph_ref, lvid_type vid):lgraph_ref(lgraph_ref),vid(vid) { }

       /// \brief Returns a constant reference to the data on the vertex.
       const vertex_data_type& data() const {
         return lgraph_ref.vertex_data(vid);
       }
       /// \brief Returns a reference to the data on the vertex.
       vertex_data_type& data() {
         return lgraph_ref.vertex_data(vid);
       }
       /// \brief Returns the number of in edges of the vertex.
       size_t num_in_edges() const {
         return lgraph_ref.num_in_edges(vid);
       }
       /// \brief Returns the number of out edges of the vertex.
       size_t num_out_edges() const {
         return lgraph_ref.num_out_edges(vid);
       }
       /// \brief Returns the ID of the vertex.
       lvid_type id() const {
         return vid;
       }
       /// \brief Returns a list of in edges.
       edge_list_type in_edges() {
         return lgraph_ref.in_edges(vid);
       }
       /// \brief Returns a list of out edges.
       edge_list_type out_edges() {
         return lgraph_ref.out_edges(vid);
       }
     private:
       dynamic_local_graph& lgraph_ref;
       lvid_type vid;
    };

    template<typename VertexData, typename EdgeData>
    class dynamic_local_graph<VertexData, EdgeData>::edge_type {
     public:
      edge_type(dynamic_local_graph& lgraph_ref, lvid_type _source, lvid_type _target, edge_id_type _eid) :
        lgraph_ref(lgraph_ref), _source(_source), _target(_target), _eid(_eid) { }

      /// \brief Returns a constant reference to the data on the edge.
      const edge_data_type& data() const {
        return lgraph_ref.edge_data(_eid);
      }
      /// \brief Returns a reference to the data on the edge.
      edge_data_type& data() {
        return lgraph_ref.edge_data(_eid);
      }
      /// \brief Returns the source vertex of the edge.
      vertex_type source() const {
        return vertex_type(lgraph_ref, _source);
      }

      lvid_type source_lvid() const {
    	  return _source;
      }

      /// \brief Returns the target vertex of the edge.
      vertex_type target() const {
        return vertex_type(lgraph_ref, _target);
      }

      lvid_type target_lvid() const {
    	  return _target;
      }

      /// \brief Returns the internal ID of this edge
      edge_id_type id() const { return _eid; }

     private:
      dynamic_local_graph& lgraph_ref;
      lvid_type _source;
      lvid_type _target;
      edge_id_type _eid;
    };

    template<typename VertexData, typename EdgeData>
    class dynamic_local_graph<VertexData, EdgeData>::edge_iterator :
        public boost::iterator_facade < edge_iterator,
                                        edge_type,
                                        boost::random_access_traversal_tag,
                                        edge_type> {
         public:
           enum list_type {CSR, CSC};

           edge_iterator(dynamic_local_graph& lgraph_ref, list_type _type,
                         csr_edge_iterator _iter, lvid_type _vid)
               : lgraph_ref(lgraph_ref), _type(_type), _iter(_iter), _vid(_vid) {}

         private:
           friend class boost::iterator_core_access;

           void increment() {
             ++_iter;
           }
           bool equal(const edge_iterator& other) const
           {
             ASSERT_EQ(_type, other._type);
             return _iter == other._iter;
           }
           edge_type dereference() const {
             return make_value();
           }
           void advance(int n) {
             _iter += n;
           }
           ptrdiff_t distance_to(const edge_iterator& other) const {
             return (other._iter - _iter);
           }
         private:
           edge_type make_value() const {
            typename csr_edge_iterator::reference ref = *_iter;
             switch (_type) {
              case CSC: {
                return edge_type(lgraph_ref, ref.first, _vid, ref.second);
              }
              case CSR: {
                return edge_type(lgraph_ref, _vid, ref.first, ref.second);
              }
              default: return edge_type(lgraph_ref, -1, -1, -1);
             }
           }
           dynamic_local_graph& lgraph_ref;
           const list_type _type;
           csr_edge_iterator _iter;
           const lvid_type _vid;
        }; // end of edge_iterator

    //disk, only used for disk version.
    template<typename VertexData, typename EdgeData>
        class dynamic_local_graph<VertexData, EdgeData>::edge_iterator_disk :
            public boost::iterator_facade < edge_iterator_disk,
                                            edge_type,
                                            boost::random_access_traversal_tag,
                                            edge_type> {
             public:
               enum list_type {CSR, CSC};

               edge_iterator_disk(dynamic_local_graph& lgraph_ref, list_type _type,
                             std::vector<int>::const_iterator _iter, lvid_type _vid)
                   : lgraph_ref(lgraph_ref), _type(_type), _iter(_iter), _vid(_vid) {}

             private:
               friend class boost::iterator_core_access;

               void increment() {
                 ++_iter;
               }
               bool equal(const edge_iterator_disk& other) const
               {
                 ASSERT_EQ(_type, other._type);
                 return _iter == other._iter;
               }
               edge_type dereference() const {
                 return make_value();
               }
               void advance(int n) {
                 _iter += n;
               }
               ptrdiff_t distance_to(const edge_iterator_disk& other) const {
                 return (other._iter - _iter);
               }
             private:
               edge_type make_value() const {
                 switch (_type) {
                  case CSC: {
                    return edge_type(lgraph_ref, (lvid_type)*_iter, _vid, -1); //don't support edge data
                  }
                  case CSR: {
                    return edge_type(lgraph_ref, _vid, (lvid_type)*_iter, -1);
                  }
                  default: return edge_type(lgraph_ref, -1, -1, -1);
                 }
               }
               dynamic_local_graph& lgraph_ref;
               const list_type _type;
               std::vector<int>::const_iterator _iter;
               const lvid_type _vid;
            }; // end of edge_iterator

} // end of namespace


namespace std {
  /**
   * Swap two graphs
   */
  template<typename VertexData, typename EdgeData>
  inline void swap(graphlab::dynamic_local_graph<VertexData,EdgeData>& a,
                   graphlab::dynamic_local_graph<VertexData,EdgeData>& b) {
    a.swap(b);
  } // end of swap
}; // end of namespace std


#include <graphlab/macros_undef.hpp>
#endif
