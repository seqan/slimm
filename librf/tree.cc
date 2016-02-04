/**
 * @file
 * @brief Tree implementation
 * TODO: fix nodes to be contiguous
 * use ncur idea
 * if current node is split
 * last_node +1 (left)
 * last_node +2 (right)
 * ncur+=2
 */
#include "librf/tree.h"
#include "librf/instance_set.h"
#include "librf/discrete_dist.h"
#include "librf/weights.h"
#include "librf/utils.h"
#include <float.h>
#include <deque>
#include <set>
#include <map>
// binary tree implcit in array
// rows in sorted_inum matrix are arranged in a similar style

namespace librf {

const int Tree::kLeft = 0;
const int Tree::kRight = 1;

Tree::Tree(istream& in):
              // if we load the tree from disk, there is no training data set
              set_(InstanceSet()),
              // also there is no list of weights
              weight_list_(NULL),
              sorted_inum_(NULL),
              temp(NULL),
              move_left(NULL)
{
  read(in);
}
/**
 * Constructor
 * @param set training set
 * @param weights weights --- presumably from bagging process
 * @param K number of vars to try per split
 * @param min_size minimum number of instances in a node
 * @param min_gain minimum information gain for making a split
 * @param seed random seed
 */
Tree::Tree(const InstanceSet& set,
           weight_list* weights,
           int K,
           int min_size,
           float min_gain,
           unsigned int seed
           ) :
                             set_(set),
                             weight_list_(weights),
                             K_(K),
                             min_size_(min_size),
                             min_gain_(min_gain),
                             num_attributes_(set.num_attributes()),
                             num_instances_(set.size()),
                             // stride_(set.size()),
                             split_nodes_(0), terminal_nodes_(0),
                             rand_seed_(seed)
{
}
/***
 * Copy the sorted indices from the training set
 * into our special matrix (sorted_inum)
 */
void Tree::copy_instances() {
  sorted_inum_ = new uint16*[num_attributes_];
  for (int i = 0; i <num_attributes_; ++i) {
    const vector<int>& sorted = set_.get_sorted_indices(i);
    sorted_inum_[i] = new uint16[num_instances_];
    for (int j = 0; j < num_instances_; ++j) {
      sorted_inum_[i][j] = sorted[j];
    }
  }
  temp = new int[num_instances_];
  move_left = new uchar[num_instances_];
}


/**
 * Do the work of growing the tree
 * - Copy the data into special matrix
 * - Build the tree
 * - Delete special matrix
 */
void Tree::grow() {
  copy_instances();
  build_tree(min_size_);
  // delete sorted_inum
  if (sorted_inum_ != NULL) {
    for (int i = 0; i <num_attributes_; ++i) {
      delete [] sorted_inum_[i];
    }
    delete [] sorted_inum_;
  }
  delete [] temp;
  delete [] move_left;
}


Tree::~Tree() {
  delete  weight_list_;
}


/** Save a tree to disk
 * Important things to record:
 *    - Nodes (all active nodes need to be written)
 *    - WeightList ? - this seems irrelevant without the instance set
 *    - Statistics? 
 */
void Tree::write(ostream& o) const {
  o << "Tree: " << nodes_.size() << endl;
  // Loop through active nodes
  for (int i = 0; i < nodes_.size(); ++i) {
    // Write the node number
    o << i << " ";
    nodes_[i].write(o);
  }
}

/**
 * Read the tree from disk
 */
void Tree::read(istream& in) {
  string spacer;
  int num_nodes;
  in >> spacer >> num_nodes;
  nodes_.resize( num_nodes);
  for (int i = 0; i < num_nodes; ++i) {
    int cur_node;
    in >> cur_node;
    nodes_[cur_node].read(in);
  }
}


void Tree::build_tree(int min_size) {
  int built_nodes = 0;
  // set up ROOT NODE (constains all instances)
  add_node(0, num_instances_, 0);
  do {
    build_node(built_nodes, min_size_);
    built_nodes++;
  } while (built_nodes < nodes_.size());
}

void Tree::mark_terminal(tree_node* n) {
  n->status = TERMINAL;
  terminal_nodes_++;
}

void Tree::mark_split(tree_node* n, uint16 split_attr, float split_point) {
  n->status = SPLIT;
  n->attr = split_attr;
  n->split_point = split_point;
  n->left = nodes_.size(); // last_node + 1 (due to zero indexing)
  n->right = nodes_.size() + 1; // last_node +2 
  split_nodes_++;
  vars_used_.insert(split_attr);
}


void Tree::add_node(uint16 start, uint16 size, uchar depth) {
  tree_node n;
  n.status = BUILD_ME;
  n.start = start;
  n.size = size;
  n.depth = depth;
  nodes_.push_back(n);
}

void Tree::build_node(uint16 node_num, uint16 min_size) {
  // cout << "building node " << node_num <<endl;
  uint16 nodes_created = 0;
  assert(node_num < nodes_.size());
  tree_node* n = &nodes_[node_num];
  // Calculate starting entropy
  DiscreteDist d;
  uint16 nstart = n->start;
  uint16 nend = n->start + n->size;
  for (uint16 i = n->start; i < nend; ++i) {
    uint16 instance = sorted_inum_[0][i];
    d.add(set_.label(instance), (*weight_list_)[instance]);
  }
  n->entropy = d.entropy_over_classes();
  // cout << "entropy: " << n-> entropy << endl;
  n->label = d.mode();

  // Min_size or completely pure check
  if (n->size <= min_size || n->entropy == 0) { //|| n->depth >= (max_depth_ -1)) {
    mark_terminal(n);
    /* if (n->size <= min_size) {
       cout << "terminal due to size of " << n->size << endl;
    } else if (n->entropy==0) {
      cout << "terminal due to zero entropy" << endl;
    } else  {
      cout << "terminal due to depth: " << int(n->depth) << endl;
    }*/
    return;
  }

  int split_attr, split_idx;
  float split_point, split_gain;
  vector<int> attrs;
  random_sample(num_attributes_, K_, &attrs, &rand_seed_);
  find_best_split(n, attrs, &split_attr, &split_idx, &split_point, &split_gain);
  if (split_gain > min_gain_) {
    mark_split(n, split_attr, split_point);
    move_data(n, split_attr, split_idx);
    uint16 left_size = split_idx - n->start + 1;
    uint16 right_size = n->size - left_size;
    add_node(n->start, left_size, n->depth + 1);
    add_node(split_idx + 1, right_size, n->depth + 1);
   } else {
    // cout << "couldn't find a split" << endl;
    mark_terminal(n);
   }
}

void Tree::move_data(tree_node* n, uint16 split_attr, uint16 split_idx) {
// PRE-CONDITION
// the same number of distinct case numbers are found in
// sorted_inum_[m][nstart-nend] for all m

// Step 1:
// Create an indicator bit set for moving left
// ex. move_left[instance number]
  // vector<uchar> move_left(set_.size(), 0);
  uint16 nstart = n->start;
  uint16 nend = nstart + n->size;
  for (int i =0; i < num_instances_; ++i) {
    move_left[i] = 0;
  }
  // int split_stride = split_attr * stride_;
  for (uint16 i = nstart; i <=split_idx; ++i) {
    //int instance_num = sorted_inum_[split_stride + i];
    int instance_num = sorted_inum_[split_attr][i];
    move_left[instance_num] = 1;
  }

// Step 2:
// For every attribute
//    fill a temporary vector -- move left | move right
//    write this back to the sorted_inum_
  // vector<int> temp(n->size);
  for (int attr = 0; attr < num_attributes_; ++attr) {
    int left = n->start;
    int right = split_idx + 1;
    // int attr_stride = attr * stride_;
    // Move instance numbers Left and right
    for (uint16 i = nstart; i < nend; ++i) {
      // int instance_num = sorted_inum_[attr_stride + i];
       int instance_num = sorted_inum_[attr][i];
      //cout << i << " ";
      //cout << instance_num << " ";
      //cout << int(move_left[instance_num]) << endl;
      if (move_left[instance_num]) {
        assert(left < num_instances_);
        temp[left++] = instance_num;
      } else {
        assert(right < num_instances_);
        temp[right++] = instance_num;
      }
    }
    // Write temp back to sorted_inum
    for (uint16 i = nstart; i < nend; ++i) {
      //sorted_inum_[attr_stride + i] = temp[i];
      sorted_inum_[attr][i] = temp[i];
    }
   }
// POST-CONDITION
// sorted_instances_[m][nstart-split] are consistent
// sorted_instances_[m][split-nend] are consistent
}
void Tree::find_best_split(tree_node* n, const vector<int>& attrs,
                           int* split_attr, int* split_idx,
                           float* split_point, float* split_gain) {
  float best_gain = -DBL_MAX;
	int best_attr = -1;
  int best_split_idx = -1;
	float best_split_point = -DBL_MAX;
	for (int i = 0; i < attrs.size(); ++i) {
    int attr =attrs[i];
    // cout << "investigating attr #" << attr <<endl;
		int curr_split_idx = -999;
    float curr_split_point = -999;
		float curr_gain = -DBL_MAX;
		find_best_split_for_attr(n, attr, n->entropy, &curr_split_idx,
                             &curr_split_point, &curr_gain);
    // cout << attr << ":" << curr_split_point << "->" << curr_gain <<endl;
		if (curr_gain > best_gain) {
				best_gain = curr_gain;
				best_split_idx = curr_split_idx;
        best_split_point = curr_split_point;
				best_attr = attr;
        assert(best_split_idx >=0);
        assert(best_split_idx < num_instances_);
		}
	}
  // get the split point
	*split_point = best_split_point;
	*split_attr = best_attr;
  *split_idx = best_split_idx;
	*split_gain = best_gain;
}



/*
void Tree::find_best_split_for_attr_gini(tree_node* n,
                                    int attr,
                                    float prior_entropy,
                                    int* split_idx,
                                    float* split_point,
                                    float* best_gain) {
  int nstart = n->start;
  int nend = n->start + n->size;
  DiscreteDist split_dist[2];
  // int attr_stride = attr * stride_;
  // Move all the instances into the right split at first
  for (int i = nstart; i < nend; ++i) {
    //int inst_no = sorted_inum_[attr_stride + i];
    int inst_no = sorted_inum_[attr][i];
    split_dist[kRight].add(set_.label(inst_no),
                           (*weight_list_)[inst_no]);
  }
  //cout << "Right Dist: " << endl;
  // split_dist[kRight].print();
  // cout << "debug sorted_inum" << endl;
  //for (int i = nstart; i < nend; ++i) {
  //  cout << int(set_.label(sorted_inum_[attr][i])) << endl;
  //}
  // set up initial values
  *best_gain = -DBL_MAX;
  //int next = sorted_inum_[attr_stride + nstart];
  int next = sorted_inum_[attr][nstart];
  float next_value = set_.get_attribute(next, attr);
  // Look for splits
  for (int i = nstart; i < nend - 1; ++i) {
    int cur = next;
    //next = sorted_inum_[attr_stride +  i + 1];
    next = sorted_inum_[attr][ i + 1];
    // cout << "cur is : " << cur <<endl;
    int label = int(set_.label(cur));
    int weight = (*weight_list_)[cur];
    // cout << "moving " << label << " weight " << weight << endl;
    split_dist[kRight].remove(label, weight);
    split_dist[kLeft].add(label, weight);
    float cur_value =  next_value;
    next_value = set_.get_attribute(next, attr);
    if (cur_value < next_value) {
      // Calculate gain (can be sped up with incremental calculation!
      float split_entropy = DiscreteDist::entropy_conditioned(split_dist, 2);
      float curr_gain = prior_entropy - split_entropy;
      // cout << "split point: " << (cur_value + next_value)/2.0 << " gain: " << curr_gain << endl;
      if (curr_gain > *best_gain) {
        *best_gain = curr_gain;
        *split_idx = i;
        *split_point = (cur_value + next_value)/2.0;
      }
    }
  }
}
*/

void Tree::find_best_split_for_attr(tree_node* n,
                                    int attr,
                                    float prior_entropy,
                                    int* split_idx,
                                    float* split_point,
                                    float* best_gain) {
  int nstart = n->start;
  int nend = n->start + n->size;
  DiscreteDist split_dist[2];
  // int attr_stride = attr * stride_;
  // Move all the instances into the right split at first
  for (int i = nstart; i < nend; ++i) {
    //int inst_no = sorted_inum_[attr_stride + i];
    int inst_no = sorted_inum_[attr][i];
    split_dist[kRight].add(set_.label(inst_no),
                           (*weight_list_)[inst_no]);
  }
  //cout << "Right Dist: " << endl;
  // split_dist[kRight].print();
  // cout << "debug sorted_inum" << endl;
  //for (int i = nstart; i < nend; ++i) {
  //  cout << int(set_.label(sorted_inum_[attr][i])) << endl;
  //}
  // set up initial values
  *best_gain = -DBL_MAX;
  //int next = sorted_inum_[attr_stride + nstart];
  int next = sorted_inum_[attr][nstart];
  float next_value = set_.get_attribute(next, attr);
  //uint16 next_value = set_.get_rank(next, attr);
  // Look for splits
  for (int i = nstart; i < nend - 1; ++i) {
    int cur = next;
    //next = sorted_inum_[attr_stride +  i + 1];
    next = sorted_inum_[attr][ i + 1];
    // cout << "cur is : " << cur <<endl;
    int label = int(set_.label(cur));
    int weight = (*weight_list_)[cur];
    // cout << "moving " << label << " weight " << weight << endl;
    split_dist[kRight].remove(label, weight);
    split_dist[kLeft].add(label, weight);
    float cur_value =  next_value;
    //uint16 cur_value =  next_value;
    next_value = set_.get_attribute(next, attr);
    // next_value = set_.get_rank(next, attr);
    if (cur_value < next_value) {
      // Calculate gain (can be sped up with incremental calculation!
      float split_entropy = DiscreteDist::entropy_conditioned(split_dist, 2);
      float curr_gain = prior_entropy - split_entropy;
      // cout << "split point: " << (cur_value + next_value)/2.0 << " gain: " << curr_gain << endl;
      if (curr_gain > *best_gain) {
        *best_gain = curr_gain;
        *split_idx = i;
        *split_point = (cur_value + next_value)/2.0;
      }
    }
  }
}
/*
void Tree::write_dot(ostream& out) {
  // if root node 
  //   write header information
  // Recursive call to children
  subtrees_[kLeft]->write_dot(out);
  subtrees_[kRight]->write_dot(out);
  // if root node
  //   write footer information
}
*/
/*
int Tree::predict(const Instance& c) const {
  // base case
  if (terminal_ ==  true) {
    return label_;
  }
  assert(subtrees_[kLeft] != NULL);
  assert(subtrees_[kRight] != NULL);
  // otherwise we will have to do the attribute test
  if (c.features_[split_attr_] < split_point_) {
    return subtrees_[kLeft]->predict(c);
  } else {
    return subtrees_[kRight]->predict(c);
  }
}
*/
int Tree::predict(const InstanceSet& set, int instance_no, int *terminal) const {
  //base case
  bool result = false;
  int cur_node = 0;
  int label = 0;
  while (!result) {
    const tree_node* n = &nodes_[cur_node];
    assert(n->status == TERMINAL || n->status == SPLIT);
    if (n->status == TERMINAL) {
      result = true;
      label = n->label;
    } else {
      if (set.get_attribute(instance_no, n->attr) < n->split_point) {
        cur_node = n->left;
      } else {
        cur_node = n->right;
      }
    }
  }
  if (terminal != NULL) {
    *terminal = cur_node;
  }
  return label;
}

int Tree::predict_skew(const InstanceSet& set, int instance_no, float* skew,
                       int* terminal) const {
  //base case
  bool result = false;
  int cur_node = 0;
  int label = 0;
  (*skew) = 0;
  int len = 0;
  while (!result) {
    const tree_node* n = &nodes_[cur_node];
    assert(n->status == TERMINAL || n->status == SPLIT);
    if (n->status == TERMINAL) {
      result = true;
      label = n->label;
    } else {
      len++;
      if (set.get_attribute(instance_no, n->attr) < n->split_point) {
        cur_node = n->left;
      } else {
        (*skew)++;
        cur_node = n->right;
      }
    }
  }
   (*skew) = (*skew) / len;

  if (terminal != NULL) {
    *terminal = cur_node;
  }
  return label;
}




/***
 * Special logging predict method
 * Will track which split variables and split points were used
 */
int Tree::predict(const InstanceSet& set, int instance_no,
                  vector<pair<int, float> > * nodes) const {
  //base case
  bool result = false;
  int cur_node = 0;
  int label = 0;
  while (!result) {
    const tree_node* n = &nodes_[cur_node];
    assert(n->status == TERMINAL || n->status == SPLIT);
    if (n->status == TERMINAL) {
      result = true;
      label = n->label;
    } else {
      nodes->push_back(make_pair(n->attr, n->split_point));
      if (set.get_attribute(instance_no, n->attr) < n->split_point) {
        cur_node = n->left;
      } else {
        cur_node = n->right;
      }
    }
  }
  return label;
}





float Tree::oob_accuracy () const{
  int correct = 0;
  int total = 0;
  // Loop through training set looking for instances with weight 0
  for (int i =0; i < num_instances_; ++i) {
    if ((*weight_list_)[i] ==0) {
      if (predict(set_, i) == set_.label(i))
        correct++;
      total++;
    }
  }
  return float(correct) / total;
}

void Tree::oob_predictions(vector<DiscreteDist>* predicts) const{
  // Loop through training set looking for instances with weight 0
  for (int i =0; i < num_instances_; ++i) {
    if ((*weight_list_)[i] ==0) {
      (*predicts)[i].add(predict(set_,i));
    }
  }
}

float Tree::training_accuracy() const {
  int correct = 0;
  for (int i =0; i < num_instances_; ++i) {
    if (predict(set_, i) == set_.label(i))
      correct++;
  }
  return float(correct) / num_instances_;
}

float Tree::testing_accuracy(const InstanceSet& set) const {
  int correct = 0;
  for (int i =0; i < set.size(); ++i) {
    if (predict(set, i) == set.label(i))
      correct++;
  }
  return float(correct) / set.size();
}


void Tree::print() const {
  int cur_node = 0;
  print_node(cur_node);
  cout << "Training acc: " << training_accuracy() << endl;
  int nonzero = 0;
  for (int i = 0; i < set_.size(); ++i) {
    if ((*weight_list_)[i] != 0)
        nonzero++;
  }
  cout << "nonzero instances: " << nonzero << endl;
}

void Tree::print_node(int n) const{
  const tree_node& node = nodes_[n];
  cout << "Tree with " << nodes_.size() << " nodes " << endl;
  cout << "Split nodes: " << split_nodes_ <<endl;
  cout << "Terminal nodes: " << terminal_nodes_ << endl;
}

//generate scores for all variables
void Tree::variable_importance(vector<float>* score,
                                unsigned int* seed) const{
  // build subset
  InstanceSet* subset = InstanceSet::create_subset(set_, *weight_list_);
  // get the oob accuracy before we start
  int correct = 0;
  for (int i = 0; i < subset->size(); ++i) {
    if (predict(*subset,i) == subset->label(i)) {
      correct++;
    }
  }
  float oob_acc = oob_accuracy();
  score->resize(set_.num_attributes());
  for (int i = 0; i < set_.num_attributes(); ++i) {
    if (vars_used_.find(i) != vars_used_.end()) {
      // make a backup copy of the variable 
      vector<float> backup;
      subset->save_var(i, &backup);
      // shuffle the values in this variable around
      subset->permute(i, seed);
      int permuted = 0;
      for (int j = 0; j < subset->size(); ++j) {
        if (predict(*subset,j) == subset->label(j)) {
          permuted++;
        }
      }
      // decrease in accuracy!
      (*score)[i] = (correct - permuted);
      // restore the proper stuff
      subset->load_var(i, backup);
    } else {
      (*score)[i] = 0.0;
    }
  }
  delete subset;
}
bool Tree::oob(int instance_no) const {
  return ((*weight_list_)[instance_no] == 0);
}


void Tree::compute_proximity(const InstanceSet& set,
                               vector<vector<float> >* prox,
                               bool oob,
                               int limit) const{
  // Limit is useful for unsupervised case
  // No need to calculate proximity for synthetic test set
  if (limit == -1) { // default sentinel value
    limit = set.size();
  }
  vector<int> nodes;
  for (int i = 0; i < limit; ++i) {
    int node_i;
    predict(set, i, &node_i);
    nodes.push_back(node_i);
  }
  for (int i = 0; i < limit; ++i) {
    int node_i = nodes[i];
    for (int j = i + 1; j < limit; ++j) {
      if (node_i == nodes[j]) {
        (*prox)[i][j] += 1.0;
        (*prox)[j][i] += 1.0;
      }
    }
  }
}



void Tree::compute_skewed_proximity(const InstanceSet& set,
                               vector<vector<float> >* prox,
                               bool oob,
                               int limit) const{
  // Limit is useful for unsupervised case
  // No need to calculate proximity for synthetic test set
  if (limit == -1) { // default sentinel value
    limit = set.size();
  }
  vector<int> nodes;
  vector<float> skews;
  for (int i = 0; i < limit; ++i) {
      int node_i;
      float skew;
      predict_skew(set, i, &skew, &node_i);
      nodes.push_back(node_i);
      skews.push_back(skew);
  }
  for (int i = 0; i < limit; ++i) {
    int node_i = nodes[i];
    for (int j = i + 1; j < limit; ++j) {
      int node_j = nodes[j];
      if (node_i == node_j) {
        (*prox)[i][j] += skews[i];
        (*prox)[j][i] += skews[i];
      }
    }
  }
}

} //namespace
