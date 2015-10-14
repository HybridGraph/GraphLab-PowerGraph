#ifndef GRAPHLAB_DYNAMIC_CSR_STORAGE
#define GRAPHLAB_DYNAMIC_CSR_STORAGE

#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>

#include <graphlab/graph/dynamic_local_graph.hpp>

#include <graphlab/util/generics/counting_sort.hpp>
#include <graphlab/util/generics/block_linked_list.hpp>

#include <graphlab/logger/logger.hpp>
#include <graphlab/logger/assertions.hpp>

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>

#include <boost/iterator/permutation_iterator.hpp>
#include <sys/stat.h>

namespace graphlab {
  /**
   * A compact key-value(s) data structure using Compressed Sparse Row format.
   * The key has type size_t and can be assolicated with multiple values of valuetype.
   *
   * The core operation of is querying the list of values associated with the query key
   * and returns the begin and end iterators via <code>begin(id)</code> and <code>end(id)</code>.
   *
   * Also, this class supports insert (and batch insert) values associated with any key. 
   */

  template<typename valuetype, typename sizetype=size_t,
           uint32_t blocksize=(4096-20)/(4*sizeof(valuetype))> // the block size makes the block fit in a memory page
  class dynamic_csr_storage_disk {
   public:
     typedef block_linked_list<valuetype, blocksize> block_linked_list_t;
     typedef typename block_linked_list_t::iterator iterator;
     typedef typename block_linked_list_t::const_iterator const_iterator;
     typedef typename block_linked_list_t::blocktype blocktype;
     typedef valuetype value_type;

   private:
        std::vector<iterator> value_ptrs;
        block_linked_list_t values;

        size_t ver_num;
        std::vector<int> degree;
        std::vector<int> edges;
        std::string fileName_Prefix;
        int idx;
        std::ifstream fin;

        int block_size;
        int block_id;
        int block_num;

        long read_bytes;

   public:
     dynamic_csr_storage_disk() {
    	 idx = -1;
    	 block_size = 1;
    	 block_id = -1;
    	 block_num = 1;
    	 read_bytes = 0;
     }

     /** Should be invoked before computing */
     void init(const std::string _fileName_Prefix, int _block_size, int _ver_num) {
    	 block_size = _block_size;
    	 fileName_Prefix = _fileName_Prefix;
    	 degree.resize(_ver_num);
    	 ver_num = (size_t)_ver_num;
    	 //logstream(LOG_INFO) << "ver_num=" << ver_num << std::endl;
     }

     int getBlockId(int lvid) {
    	 /*logstream(LOG_INFO) << "lvid=" << lvid << ", block_size=" << block_size
    			 << ", bid=" << (lvid/block_size) << std::endl;*/
    	 return (lvid/block_size);
     }

     std::string getFileName(int _block_id) {
    	 std::stringstream ss;
    	 ss << _block_id;
    	 return (fileName_Prefix+ss.str());
     }

     /** disk, only used by iteration computation for disk version. */
     inline std::vector<int>::iterator begin_disk(size_t id) {
    	 //logstream(LOG_INFO) << "disk debug. vid in storage=" << id << std::endl;
    	 bool yes = read((int)id);
    	 ASSERT_TRUE(yes);

    	 return edges.begin();
     }

     // disk, only used by disk version when computing.
     inline std::vector<int>::iterator end_disk(size_t id) {
    	 ASSERT_EQ((size_t)idx, id);
    	 return edges.end();
     }

     inline void set_num_edges(int id, int num) {
    	 degree[id] = num;
    	 //logstream(LOG_INFO) << "id=" << id << ", num=" << num << std::endl;
     }

     inline int num_edges(size_t id) const {
    	 //ASSERT_EQ((size_t)idx, id);
    	 //return edges.size();
    	 //ASSERT_LT(id, ver_num);
    	 //logstream(LOG_INFO) << " request num_edges, id=" << id << std::endl;
    	 if (id < ver_num) {
    		 return degree[id];
    	 } else {
    		 return 0;
    	 }
     }

     void resetLocalBefIte() {
    	 ASSERT_TRUE(closeEdgeFile(0));
    	 idx = -1;
    	 block_id = -1;
    	 read_bytes = 0;
     }

     long getReadBytes() {
    	 return read_bytes;
     }

     bool openEdgeFile(int iter_counter) {
    	 if (block_id < 0) {
    		 return false;
    	 }
    	 std::string fileName = getFileName(block_id);
    	 if (fin.is_open()) {
    		 fin.close();
    	 }
    	 fin.open(fileName.c_str(), std::ios::binary);
    	 if (!fin) {
    		 return false;
    	 } else {
    		 //logstream(LOG_INFO) << "open edge block-" << block_id << std::endl;
    		 struct stat info;
    		 stat(fileName.c_str(), &info);
    		 read_bytes = read_bytes + info.st_size;
    		 return true;
    	 }
     }

     bool closeEdgeFile(int iter_counter) {
    	 if (fin.is_open()) {
    		 fin.close();
    		 //logstream(LOG_INFO) << "close edge block-" << block_id << std::endl;
    	 }
    	 return true;
     }

     bool read(int gid) {
    	 //logstream(LOG_INFO) << "gid=" << gid << ", idx=" << idx << std::endl;
    	 if (gid > idx) {
    		 idx = gid;
    	 }

    	 int bid = getBlockId(idx);
    	 if (bid != block_id) {
    		 ASSERT_TRUE(closeEdgeFile(0));

    		 block_id = bid;
    		 ASSERT_TRUE(openEdgeFile(0));
    	 }

    	 int vid, len, eid;
    	 bool find = false;
    	 std::vector<int>().swap(edges);
    	 //logstream(LOG_INFO) << "num_of_edges = " << degree[idx] << ", total_ver_num=" << ver_num << std::endl;
    	 while (fin.peek()!=EOF) {
    		 //logstream(LOG_INFO) << "begin to read file..." << std::endl;
    		 fin.read((char*)(&vid), sizeof(vid));
    		 fin.read((char*)(&len), sizeof(len));
    		 //logstream(LOG_INFO) << "read vid=" << vid << ", current idx=" << idx <<std::endl;
    		 if (vid > idx) {
    			 //logstream(LOG_INFO) << "judge, read vid=" << vid << ", current idx=" << idx <<std::endl;
    			 break;
    		 } else if (vid < idx) {
    			 for (int i = 0; i < len; i++) {
    				 fin.read((char*)(&eid), sizeof(eid));
    				 //logstream(LOG_INFO) << "skip eid=" << eid << std::endl;
    			 } // skip no-used edges
    			 continue;
    		 } else {
    			 for (int i = 0; i < len; i++) {
    				 fin.read((char*)(&eid), sizeof(eid));
    				 edges.push_back(eid);
    				 //logstream(LOG_INFO) << "eid=" << eid << std::endl;
    			 }
    			 find = true;
    			 break;
    		 }
    	 }

    	 if (find) {
    		 //logstream(LOG_INFO) << "find" << std::endl;
    		 return true;
    	 } else {
    		 //logstream(LOG_INFO) << "not find" << std::endl;
    		 return false;
    	 }
     }

     /// Number of keys in the storage.
     inline size_t num_keys() const { return value_ptrs.size(); }

     /// Number of values in the storage.
     inline size_t num_values() const { return values.size(); }

     /// Return iterator to the begining value with key == id
     inline iterator begin(size_t id) {
       return id < num_keys() ? value_ptrs[id] : values.end();
     }

     /// Return iterator to the ending+1 value with key == id 
     inline iterator end(size_t id) {
    	 //logstream(LOG_INFO) << "disk debug. id=" << id << std::endl;
       return (id+1) < num_keys() ? value_ptrs[id+1] : values.end();
     }

     /// Return iterator to the begining value with key == id 
     inline const_iterator begin(size_t id) const {
    	 //logstream(LOG_INFO) << "disk debug. vid in mem=" << id << std::endl;
       return id < num_keys() ? value_ptrs[id] : values.end();
     } 

     /// Return iterator to the ending+1 value with key == id 
     inline const_iterator end(size_t id) const {
       return (id+1) < num_keys() ? value_ptrs[id+1] : values.end();
     }

     /**
      * Wrap the index vector and value vector into csr_storage.
      * Check the property of the input vector.
      * Clean up the input on finish.
      */
     void wrap(std::vector<sizetype>& valueptr_vec,
    		 std::vector<valuetype>& value_vec) {
    	 //logstream(LOG_INFO) << "!!!!!!!!!!!!wrap() is invoked" << std::endl;
    	 for (ssize_t i = 1; i < (ssize_t)valueptr_vec.size(); ++i) {
    		 ASSERT_LE(valueptr_vec[i-1], valueptr_vec[i]);
    		 ASSERT_LT(valueptr_vec[i], value_vec.size());
    	 }

    	 values.assign(value_vec.begin(), value_vec.end());
    	 sizevec2ptrvec(valueptr_vec, value_ptrs);

    	 std::vector<value_type>().swap(value_vec);
    	 std::vector<sizetype>().swap(valueptr_vec);
     }

     void swap(dynamic_csr_storage_disk<valuetype, sizetype>& other) {
    	 logstream(LOG_FATAL) << "this function has been deleted!" << std::endl;
     }

     void clear() {
    	 //logstream(LOG_INFO) << "!!!!!!!!!!!!clear() is invoked" << std::endl;
    	 value_ptrs.clear();
       std::vector<iterator>().swap(value_ptrs);
       values.clear();
       edges.clear();
       std::vector<int>().swap(edges);
     }

     void load(iarchive& iarc) { 
    	 logstream(LOG_FATAL) << "this function has been deleted!" << std::endl;
     }

     void save(oarchive& oarc) const { 
    	 logstream(LOG_FATAL) << "this function has been deleted!" << std::endl;
     }

     ////////////////////// Internal APIs /////////////////
   public:
     /**
      * \internal
      */

     size_t estimate_sizeof() const {
    	 logstream(LOG_FATAL) << "this function has been deleted!" << std::endl;
    	 return (size_t)0;
     }

     void meminfo(std::ostream& out) {
    	 logstream(LOG_FATAL) << "this function has been deleted!" << std::endl;
     }

     ///////////////////// Helper Functions /////////////
   private:

     // Convert integer pointers into block_linked_list::value_iterator
     // Assuming all blocks are fully packed.
     void sizevec2ptrvec (const std::vector<sizetype>& ptrs,
                          std::vector<iterator>& out) {
       ASSERT_EQ(out.size(), 0);
       out.reserve(ptrs.size());

       // for efficiency, we advance pointers based on the previous value
       // because block_linked_list is mostly forward_traversal.
       iterator it = values.begin();
       sizetype prev = 0;
       for (size_t i = 0; i < ptrs.size(); ++i) {
         sizetype cur = ptrs[i];
         it += (cur-prev);
         out.push_back(it);
         prev = cur;
       }
     }
  }; // end of class
} // end of graphlab 
#endif
