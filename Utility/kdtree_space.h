
//part of the source code is from nanoflann lib at https://github.com/jlblancoc/nanoflann

#ifndef  KD_TREE_SPACE_PARTITION_HPP_
#define  KD_TREE_SPACE_PARTITION_HPP_

#include <vector>
#include <cassert>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <utility> 

namespace KDTreeSpace
{
	/** @addtogroup memalloc_grp Memory allocation
	  * @{ */

	/**
	 * Allocates (using C's malloc) a generic type T.
	 *
	 * Params:
	 *     count = number of instances to allocate.
	 * Returns: pointer (of type T*) to memory buffer
	 */
	template <typename T>
	inline T* allocate(size_t count = 1)
	{
		T* mem = static_cast<T*>( ::malloc(sizeof(T)*count));
		return mem;
	}


	/**
	 * Pooled storage allocator
	 *
	 * The following routines allow for the efficient allocation of storage in
	 * small chunks from a specified pool.  Rather than allowing each structure
	 * to be freed individually, an entire pool of storage is freed at once.
	 * This method has two advantages over just using malloc() and free().  First,
	 * it is far more efficient for allocating small objects, as there is
	 * no overhead for remembering all the information needed to free each
	 * object or consolidating fragmented memory.  Second, the decision about
	 * how long to keep an object is made at the time of allocation, and there
	 * is no need to track down all the objects to free them.
	 *
	 */

	const size_t     WORDSIZE=16;
	const size_t     BLOCKSIZE=8192;

	class PooledAllocator
	{
		/* We maintain memory alignment to word boundaries by requiring that all
		    allocations be in multiples of the machine wordsize.  */
		/* Size of machine word in bytes.  Must be power of 2. */
		/* Minimum number of bytes requested at a time from	the system.  Must be multiple of WORDSIZE. */


		size_t  remaining;  /* Number of bytes left in current block of storage. */
		void*   base;     /* Pointer to base of current block of storage. */
		void*   loc;      /* Current location in block to next allocate memory. */

		void internal_init()
		{
			remaining = 0;
			base = NULL;
			usedMemory = 0;
			wastedMemory = 0;
		}

	public:
		size_t  usedMemory;
		size_t  wastedMemory;

		/**
		    Default constructor. Initializes a new pool.
		 */
		PooledAllocator() {
			internal_init();
		}

		/**
		 * Destructor. Frees all the memory allocated in this pool.
		 */
		~PooledAllocator() {
			free_all();
		}

		/** Frees all allocated memory chunks */
		void free_all()
		{
			while (base != NULL) {
				void *prev = *(static_cast<void**>( base)); /* Get pointer to prev block. */
				::free(base);
				base = prev;
			}
			internal_init();
		}

		/**
		 * Returns a pointer to a piece of new memory of the given size in bytes
		 * allocated from the pool.
		 */
		void* malloc(const size_t req_size)
		{
			/* Round size up to a multiple of wordsize.  The following expression
			    only works for WORDSIZE that is a power of 2, by masking last bits of
			    incremented size to zero.
			 */
			const size_t size = (req_size + (WORDSIZE - 1)) & ~(WORDSIZE - 1);

			/* Check whether a new block must be allocated.  Note that the first word
			    of a block is reserved for a pointer to the previous block.
			 */
			if (size > remaining) {

				wastedMemory += remaining;

				/* Allocate new storage. */
				const size_t blocksize = (size + sizeof(void*) + (WORDSIZE-1) > BLOCKSIZE) ?
							size + sizeof(void*) + (WORDSIZE-1) : BLOCKSIZE;

				// use the standard C malloc to allocate memory
				void* m = ::malloc(blocksize);
				if (!m) {
					fprintf(stderr,"Failed to allocate memory.\n");
					return NULL;
				}

				/* Fill first word of new block with pointer to previous block. */
				static_cast<void**>(m)[0] = base;
				base = m;

				size_t shift = 0;
				//int size_t = (WORDSIZE - ( (((size_t)m) + sizeof(void*)) & (WORDSIZE-1))) & (WORDSIZE-1);

				remaining = blocksize - sizeof(void*) - shift;
				loc = (static_cast<char*>(m) + sizeof(void*) + shift);
			}
			void* rloc = loc;
			loc = static_cast<char*>(loc) + size;
			remaining -= size;

			usedMemory += size;

			return rloc;
		}

		/**
		 * Allocates (using this pool) a generic type T.
		 *
		 * Params:
		 *     count = number of instances to allocate.
		 * Returns: pointer (of type T*) to memory buffer
		 */
		template <typename T>
		T* allocate(const size_t count = 1)
		{
			T* mem = static_cast<T*>(this->malloc(sizeof(T)*count));
			return mem;
		}

	};
	/** @} */

	template <typename ElementType,  typename IndexType = size_t>
	class PartitioningKDTree
	{
	private:
		/** Hidden copy constructor, to disallow copying indices (Not implemented) */
		PartitioningKDTree(const PartitioningKDTree<ElementType, IndexType>&)=delete;
	protected:

		/**
		 *  Array of indices to vectors in the dataset.
		 */
		std::vector<IndexType> m_vind;

		/**
		 * The dataset used by this index
		 */
		const std::vector<std::vector<ElementType>> &m_dataset; //!< The source of our data

		size_t m_size; //!< Number of current points in the dataset
		int m_dim;  //!< Dimensionality of each data point

		/*--------------------- Internal Data Structures --------------------------*/
		struct Node
		{
			/** Union used because a node can be either a LEAF node or a non-leaf node, so both data fields are never used simultaneously */
			union{			
				struct {
					int    idx_region;
				}lr;
				struct{
					int          divfeat; //!< Dimension used for subdivision.
					ElementType pivot; // pivot value for division
					IndexType    idx_sample;		// index
					ElementType low, high;	//boundary for the box to be cutted.
				}sub;
			};			
			Node* child1, *child2;  //!< Child nodes (both=NULL mean its a leaf node)
		};
		typedef Node* NodePtr;

		/** The KD-tree used to find regions */
		NodePtr m_root;

	 	typedef  struct {
			std::vector<std::pair<ElementType, ElementType>> box;
			double rat=1.0;
		}BoundingBox;
		BoundingBox m_rootBbox;

		/**
		 * Pooled memory allocator.
		 *
		 * Using a pooled memory allocator is more efficient
		 * than allocating memory directly when there is a large
		 * number small of memory allocations.
		 */
		PooledAllocator pool;
		double m_lrat = 0, m_srat = 1;
		int m_lbox=0, m_sbox=0;
	public:
		std::vector<BoundingBox> region;

		PartitioningKDTree(const int dimensionality, const std::vector<std::vector<ElementType>>& inputData, const std::vector<std::pair<ElementType, ElementType>> &initBBox) :
			m_dataset(inputData), m_root(NULL), m_rootBbox(initBBox)
		{
			m_size = inputData.size();
			m_dim = dimensionality;
			
			// Create a permutable array of indices to the input vectors.
			init_vind();
		}

		PartitioningKDTree(const int dimensionality, const std::vector<std::vector<ElementType>>& inputData) :
			m_dataset(inputData), m_root(NULL)
		{
			m_size = inputData.size();
			m_dim = dimensionality;

			// Create a permutable array of indices to the input vectors.
			init_vind();
		}
		void setInitBox(const std::vector<std::pair<ElementType, ElementType>> &initBBox){
			m_rootBbox.box = initBBox;
		}
		/** Standard destructor */
		~PartitioningKDTree() { }

		/** Frees the previously-built index. Automatically called within buildIndex(). */
		void freeIndex()
		{
			pool.free_all();
			m_root=NULL;
			region.clear();
		}

		/**
		 * Builds the index
		 */
		void buildIndex(int mode=1)
		{
			init_vind();
			freeIndex();
			if (mode == 1) // randomly constr
				m_root = divideTree(0, m_size, m_rootBbox);   // construct the tree
			else
				m_root = equalDivideTree(0, m_size, m_rootBbox);
		}

		/** Returns number of points in dataset  */
		size_t size() const { return m_size; }

		/** Returns the length of each point in the dataset */
		size_t veclen() const {
			return m_dim;
		}

		/**
		 * Computes the inde memory usage
		 * Returns: memory used by the index
		 */
		size_t usedMemory() const
		{
			return pool.usedMemory+pool.wastedMemory+m_dataset.size()*sizeof(IndexType);  // pool memory and vind array memory
		}

		size_t get_regionIdx(const std::vector<ElementType> & p){
			return enqury(p, m_root);
		}
		void get_leafParent(const IndexType idx, IndexType &sidx, IndexType &cutfeat, ElementType& low, ElementType&high){
		
			NodePtr result=NULL;
			leafParent(idx, m_root, NULL, result);
			sidx = result->sub.idx_sample;
			cutfeat = result->sub.divfeat;
			low = result->sub.low;
			high = result->sub.high;
		}
		
		const BoundingBox & get_rootBox(){
			return m_rootBbox;
		}
		int smallestBox(){ return m_sbox; }
		int largestBox(){ return m_lbox; }
	private:
		/** Make sure the auxiliary list \a vind has the same size than the current dataset, and re-generate if size has changed. */
		void init_vind()
		{
			// Create a permutable array of indices to the input vectors.
			m_size = m_dataset.size();
			if (m_vind.size()!=m_size) m_vind.resize(m_size);
			size_t k = 0;
			for (auto &i : m_vind) i = k++;
		}

		/// Helper accessor to the dataset points:
		inline ElementType dataset_get(size_t idx, int component) const {
			return m_dataset[idx][component];
		}


		/**
		 * Create a tree node that subdivides the list of vecs from vind[first]
		 * to vind[last].  The routine is called recursively on each sublist.
		 *
		 * @param left index of the first vector
		 * @param right index of the last vector
		 */
		NodePtr divideTree(const IndexType left, const IndexType right, BoundingBox& bbox, int depth =0)
		{
			NodePtr node = pool.allocate<Node>(); // allocate memory

			/*a leaf node,create a sub-region. */
			if ( (right-left) <= 0) {
				node->child1 = node->child2 = NULL;    /* Mark as leaf node. */	
				node->lr.idx_region = region.size();
				region.push_back(bbox);
				boxRatio(region.back(), region.size()-1);
			}
			else {
				IndexType idx;
				int cutfeat;
				ElementType cutval;
				middleSplit_(&m_vind[0]+left, right-left, idx, cutfeat, cutval, bbox,depth);

				node->sub.idx_sample = m_vind[left + idx];
				node->sub.divfeat = cutfeat;
				node->sub.low = bbox.box[cutfeat].first;
				node->sub.high = bbox.box[cutfeat].second;
				BoundingBox left_bbox(bbox);
				left_bbox.box[cutfeat].second = cutval;
				node->child1 = divideTree(left, left + idx, left_bbox,  depth + 1);

				BoundingBox right_bbox(bbox);
				right_bbox.box[cutfeat].first = cutval;
				node->child2 = divideTree(left + idx+1, right, right_bbox, depth + 1);
				node->sub.pivot = cutval;				
			}
			return node;
		}

		NodePtr equalDivideTree(const IndexType left, const IndexType right, BoundingBox& bbox, int depth = 0)
		{
			NodePtr node = pool.allocate<Node>(); // allocate memory

			/*a leaf node,create a sub-region. */
			if ((right - left) <= 0) {
				node->child1 = node->child2 = NULL;    /* Mark as leaf node. */
				node->lr.idx_region = region.size();
				region.push_back(bbox);
				boxRatio(region.back(), region.size()-1);
			}
			else {
				IndexType idx=0;
				int cutfeat = depth%m_dim;
				ElementType cutval = dataset_get(left,cutfeat);
				
				node->sub.idx_sample = m_vind[left + idx];
				node->sub.divfeat = cutfeat;
				node->sub.low = bbox.box[cutfeat].first;
				node->sub.high = bbox.box[cutfeat].second;
				BoundingBox left_bbox(bbox);
				left_bbox.box[cutfeat].second = cutval;
				node->child1 = equalDivideTree(left, left + idx, left_bbox, depth + 1);

				BoundingBox right_bbox(bbox);
				right_bbox.box[cutfeat].first = cutval;
				node->child2 = equalDivideTree(left + idx + 1, right, right_bbox, depth + 1);
				node->sub.pivot = cutval;
			}
			return node;
		}


		void computeMinMax(IndexType* ind, IndexType count, int element, ElementType& min_elem, ElementType& max_elem)
		{
			min_elem = dataset_get(ind[0],element);
			max_elem = dataset_get(ind[0],element);
			for (IndexType i=1; i<count; ++i) {
				ElementType val = dataset_get(ind[i],element);
				if (val<min_elem) min_elem = val;
				if (val>max_elem) max_elem = val;
			}
		}

		void middleSplit_(IndexType* ind, IndexType count, IndexType& index, int& cutfeat, ElementType& cutval, const BoundingBox& bbox, int depth)
		{
			
			cutfeat = depth%m_dim;
			// for a balanced kd-tree, split in the median value
			std::vector<IndexType> cur_idx(count);
			for (IndexType i = 0; i < count; ++i) {
				cur_idx[i] = ind[i];
			}
			std::nth_element(cur_idx.begin(), cur_idx.begin() + cur_idx.size() / 2, cur_idx.end(), [this, &cutfeat](const IndexType a, const IndexType b){
				return this->dataset_get(a, cutfeat) < this->dataset_get(b, cutfeat);
			});
			ElementType split_val = dataset_get(cur_idx[cur_idx.size() / 2], cutfeat);
			//.....

			ElementType min_elem, max_elem;
			computeMinMax(ind, count, cutfeat, min_elem, max_elem);

			if (split_val<min_elem) cutval = min_elem;
			else if (split_val>max_elem) cutval = max_elem;
			else cutval = split_val;

			IndexType lim1, lim2;
			planeSplit(ind, count, cutfeat, cutval, lim1, lim2);

			if (lim1>count/2) index = lim1;
			else if (lim2<count/2) index = lim2;
			else index = count/2;
		}


		/**
		 *  Subdivide the list of points by a plane perpendicular on axe corresponding
		 *  to the 'cutfeat' dimension at 'cutval' position.
		 *
		 *  On return:
		 *  dataset[ind[0..lim1-1]][cutfeat]<cutval
		 *  dataset[ind[lim1..lim2-1]][cutfeat]==cutval
		 *  dataset[ind[lim2..count]][cutfeat]>cutval
		 */
		void planeSplit(IndexType* ind, const IndexType count, int cutfeat, ElementType cutval, IndexType& lim1, IndexType& lim2)
		{
			/* Move vector indices for left subtree to front of list. */
			IndexType left = 0;
			IndexType right = count-1;
			for (;; ) {
				while (left<=right && dataset_get(ind[left],cutfeat)<cutval) ++left;
				while (right && left<=right && dataset_get(ind[right],cutfeat)>=cutval) --right;
				if (left>right || !right) break;  // "!right" was added to support unsigned Index types
				std::swap(ind[left], ind[right]);
				++left;
				--right;
			}
			/* If either list is empty, it means that all remaining features
			 * are identical. Split in the middle to maintain a balanced tree.
			 */
			lim1 = left;
			right = count-1;
			for (;; ) {
				while (left<=right && dataset_get(ind[left],cutfeat)<=cutval) ++left;
				while (right && left<=right && dataset_get(ind[right],cutfeat)>cutval) --right;
				if (left>right || !right) break;  // "!right" was added to support unsigned Index types
				std::swap(ind[left], ind[right]);
				++left;
				--right;
			}
			lim2 = left;
		}

		size_t enqury(const std::vector<ElementType> & p, NodePtr node){
			if (node->child1 == NULL&&node->child2 == NULL){
				return node->lr.idx_region;
			}
			if (p[node->sub.divfeat] < dataset_get(node->sub.idx_sample, node->sub.divfeat)){
				return enqury(p, node->child1);
			}
			else{
				return enqury(p, node->child2);
			}

		}
		void leafParent(IndexType idx_region, NodePtr node, NodePtr parent, NodePtr &result){
			if (node->child1 == NULL&&node->child2 == NULL){
				if (node->lr.idx_region == idx_region){ 
					if(node!=m_root) result = parent; 
					else{
						result = m_root;
					}
				}
				return;
			}
			if (node->child1 != NULL&&result == NULL)  leafParent(idx_region, node->child1, node, result);
			if (node->child2 != NULL&&result == NULL)  leafParent(idx_region, node->child2, node, result);
		}

		void boxRatio(BoundingBox &it,unsigned idx){
			it.rat = 1;
			for (int i = 0; i < m_dim; ++i){
				it.rat *= (it.box[i].second - it.box[i].first) / (m_rootBbox.box[i].second - m_rootBbox.box[i].first);
			}
			if (it.rat > m_lrat){
				m_lrat = it.rat;
				m_lbox = idx;
			}
			if (it.rat < m_srat){
				m_srat = it.rat;
				m_sbox = idx;
			}
		}
	}; 
/** @} */ // end of grouping
} // end of NS


#endif /* kdtree_space_HPP_ */
