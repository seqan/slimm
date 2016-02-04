/**
 * @file 
 * @brief A single decision tree
 * The implementation is a port of the Fortran implementation. 
*
 */

#ifndef _TREE_H_
#define _TREE_H_
#include "librf/types.h"
#include "librf/tree_node.h"
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <math.h>
using namespace std;

namespace librf {

class InstanceSet;
class weight_list;
class DiscreteDist;

/**
 * @brief
 * A single decision tree
 *
 *
 *
 * Strategy:
 *  - Binary tree in a fixed size array
 *    - #leaves is bounded by #instances!
 *  - Sorted instances kept in an special matrix such that:
 *    - Each node has a contiguous group of rows in the matrix
 *    - Each column gives the sorted instance order with respect
 *          to a feature
 *
 * Trees can only be created in two ways:
 *  -# load from a saved model
 *  -# grown from a certain bagging of a dataset 
 
*/
class Tree {
    public:
        /// Construct a new tree by loading it from a file
        Tree(istream& in);
        /// Construct a new tree by training
        Tree(const InstanceSet& set, weight_list* weights,
             int K, int min_size = 1,
             float min_gain = 0, unsigned int seed =0);
         ~Tree();  // clean up 
        /// predict an instance from a set
        int predict(const InstanceSet& set, int instance_no, int *terminal = NULL) const;
        int predict(const InstanceSet& set, int instance_no, vector<pair<int, float> >*) const;
        int terminal_node(const InstanceSet& set, int i) const;

        int predict_skew(const InstanceSet& set, int instance_no, float* skew, int *terminal = NULL) const;
        void compute_skewed_proximity(const InstanceSet& set,
                               vector<vector<float> >* prox,
                               bool oob,
                               int limit) const;
        void compute_proximity(const InstanceSet& set,
                               vector<vector<float> >* prox,
                               bool oob = false, int limit = -1) const;
        // void write_dot(const string& s) const;
        /// Return the accuracy for the training set
        float training_accuracy() const;
        // predict all the instances in testset and return the accuracy
        float testing_accuracy(const InstanceSet& testset) const;
        float oob_accuracy() const;
        void oob_predictions(vector<DiscreteDist> *) const;

        void variable_importance(vector<float>* scores, unsigned int* seed) const;
        void print() const;
        // do all the work -- separated this from constructor to 
        // facilitate threading
        void grow();
        void write(ostream& o) const;
        void read(istream& i);

        bool oob(int instance_no) const;
    private:
        void copy_instances();
        void move_data(tree_node* n, uint16 split_attr, uint16 split_idx);
        void find_best_split(tree_node* n,
                             const vector<int>& attrs,
                             int* split_attr, int* split_idx,
                             float* split_point, float* split_gain);
        void find_best_split_for_attr(tree_node* n,
                                      int attr,
                                      float prior,
                                      int* split_idx,
                                      float *split_point,
                                      float* best_gain);

        // Node marking
        void add_node(uint16 start, uint16 size, uchar depth);
        void mark_terminal(tree_node* n);
        void mark_split(tree_node* n, uint16 split_attr, float split_point);

        void build_tree(int min_size);
        void build_node(uint16 node_num, uint16 min_size);
        void print_node(int n) const;

        void permuteOOB(int m, double *x);
        vector<tree_node> nodes_;
        set<uint16> vars_used_;
        uint16 terminal_nodes_;
        uint16 split_nodes_;
        // get sorted indices
        // Const reference to instance set -- we don't get to delete it
        const InstanceSet& set_;
        // array of instance nums sorted by attributes 
        // this is the block array that stores which instances belong to
        // which node
        // ex: sorted_inum_[attr*stride + start]
        // uint16 * sorted_inum_;
        // uint16 stride_;
        //
        // Turns out there is not much of a gain in batch allocating the
        // 2d array (perhaps because, we only access a single column at a time)
        uint16** sorted_inum_;
        // label population 
        // uchar * sorted_labels_; necessary?
        // A single weight list for all of the instances 
        weight_list* weight_list_;
        // Depth of current tree
        // uint16 max_depth_; DEPRECATE?
        uint16 K_;
        uint16 min_size_;
        float min_gain_;
        // scratch space
        int* temp;
        uchar* move_left;
        uint16 num_instances_;
        uint16 num_attributes_;
        unsigned int rand_seed_;
        // Constants
        static const int kLeft;
        static const int kRight;
};

} // namespace
#endif
