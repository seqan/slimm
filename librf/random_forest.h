/** 
 * @file 
 * @brief Randomforest interface
 * This is the interface to manage a random forest
 */
#ifndef _RANDOM_FOREST_H_
#define _RANDOM_FOREST_H_

#include <vector>

using namespace std;

namespace librf {
class DiscreteDist;
class InstanceSet;
class Tree;
/**
 * @brief
 * RandomForest class.  Interface for growing random forests from training
 * data or loading a random forest from disk.
 */
class RandomForest {
  public:
    /// Empty constructor
    RandomForest();
    /// Constructor. (Build from training data)
    RandomForest(const InstanceSet& set,
                 int num_trees,
                 int K,
                 const vector<int>& weights = vector<int>());
    ~RandomForest();
     /// Method to predict the label
     // int predict(const Instance& c) const;
     /// Method that returns the class probability 
     // float predict_prob(const Instance& c) const;
     /// Method to predict the label
     int predict(const InstanceSet& set, int instance_no) const;

     /// Special logging method to predict the label
     int predict(const InstanceSet& set, int instance_no, vector<pair<int, float> >*) const;
     int oob_predict(int instance_no, vector<pair<int, float> >* ) const;
     float oob_predict_prob( int instance_no,
                            int label) const;
     /// Predict probability of given label
     float predict_prob(const InstanceSet& set, int instance_no, int label) const;
     /// Returns test accuracy of a labeled test set
     float testing_accuracy(const InstanceSet& testset) const;
     /// Returns training accuracy 
     float training_accuracy() const;
     void oob_votes(const InstanceSet& set, int instance_no,
                    DiscreteDist* ) const;
     void oob_predictions(vector<DiscreteDist>* predicts) const;
     void confusion_matrix(int num_labels,
                           const vector<int>&,
                           const vector<int>&) const;
     /// Returns OOB accuracy (unbiased estimate of test accuracy)
     float oob_accuracy() const;
     void oob_confusion() const;
     void test_confusion(const InstanceSet& set) const;
     /// Variable importance ranking of features
     void variable_importance(vector< pair<float, int> >* ranking,
                              unsigned int* seed) const;
     void variable_importance2(vector< pair<float, int> >* ranking,
                              unsigned int* seed) const;

     void reliability_diagram(int bins,
                              vector<pair<float, float> >*,
                              vector<int>*,
                              int label = 1) const;
     void reliability_diagram(const InstanceSet& set,
                              int bins,
                              vector<pair<float, float> >*,
                              vector<int>*,
                              int label = 1) const;
     void compute_proximity(const InstanceSet& set,
                            vector<vector<float> >* prox,
                            int limit = -1) const;
    void compute_skewed_proximity(const InstanceSet& set,
                            vector<vector<float> >* prox,
                            int limit = -1) const;


     void compute_outliers(const InstanceSet& set, int label,
                           const vector<vector<float> >& mat,
                           vector<pair< float, int> >* ranking) const;
     /// Load random forest
     void read(istream& i);
     /// Save random forest
     void write(ostream& o);
     /// Debug output
     void print() const;
  private:
    const InstanceSet& set_;  // training data set
    vector<Tree*> trees_;     // component trees in the forest
    // int max_depth_;           // maximum depth of trees (DEPRECATED)
    int K_;                   // random vars to try per split
    vector< pair<float, int> > var_ranking_; // cached var_ranking
    vector<int> class_weights_;
};

} // namespace
#endif
