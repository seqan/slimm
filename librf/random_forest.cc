/**
 * @file
 * @brief random forest implementation
 */
#include "librf/random_forest.h"
#include "librf/tree.h"
#include "librf/instance_set.h"
#include "librf/weights.h"
#include <fstream>
#include <algorithm>
namespace librf {

RandomForest::RandomForest() : set_(InstanceSet()) {}
/**
 * @param set training data
 * @param num_trees #trees to train
 * @param K #random vars to consider at each split
 * number of instances)
 */
RandomForest::RandomForest(const InstanceSet& set,
                           int num_trees,
                           int K,
                           const vector<int>& weights) :set_(set), K_(K) {
  if (weights.size() == 0) {
    class_weights_.resize(2, 1); //HARDCODE
  } else {
    class_weights_ = weights;
  }
  // cout << "RandomForest Constructor " << num_trees << endl;
  for (int i = 0; i < num_trees; ++i) {
    weight_list* w = new weight_list(set.size(), set.size());
    // sample with replacement
    for (int j = 0; j < set.size(); ++j) {
      int instance = rand() % set.size();
      w->add(instance, class_weights_[set.label(instance)]);
    }
    Tree* tree = new Tree(set, w,  K, 1, 0, rand());
    tree->grow();
    cout << "Grew tree " << i << endl;
    trees_.push_back(tree);
  }
}

RandomForest::~RandomForest() {
  for (int i = 0; i < trees_.size(); ++i) {
    delete trees_[i];
  }
}

void RandomForest::print() const {
  for (int i = 0; i < trees_.size(); ++i) {
     trees_[i]->print();
  }
}

void RandomForest::write(ostream& o) {
  o << trees_.size() << " " << K_ << endl;
  for (int i = 0; i < trees_.size(); ++i) {
    trees_[i]->write(o);
  }
}

void RandomForest::read(istream& in) {
  int num_trees, K;
  in >> num_trees >> K;
  for (int i = 0; i < num_trees; ++i) {
    trees_.push_back(new Tree(in));
  }
}

int RandomForest::predict(const InstanceSet& set, int instance_no) const {
  // Gather the votes from each tree
  DiscreteDist votes;
  for (int i = 0; i < trees_.size(); ++i) {
    int predict = trees_[i]->predict(set, instance_no);
    votes.add(predict);
  }
  return votes.mode();
}

int RandomForest::predict(const InstanceSet& set, int instance_no,
                          vector<pair<int, float> >*nodes) const {
  // Gather the votes from each tree
  DiscreteDist votes;
  for (int i = 0; i < trees_.size(); ++i) {
    int predict = trees_[i]->predict(set, instance_no, nodes);
    votes.add(predict);
  }
  return votes.mode();
}

void RandomForest::reliability_diagram(const InstanceSet&set,
                                       int bins,
                                       vector<pair<float, float> >*out,
                                       vector<int>* count, int label) const {
  float increment = 1.0 / bins;
  float half = increment / 2.0;
  vector<DiscreteDist> bin_dists(bins);
  count->resize(bins, 0);
  for (int i = 0; i < set.size(); ++i) {
    float prob = predict_prob(set, i, label);
    int bin_no = int(floor(prob/increment));
    if (bin_no == bins) {
      bin_no = bins - 1;
    }
    bin_dists[bin_no].add(set.label(i));
    (*count)[bin_no]++;
  }
  float x = half;
  for (int i = 0; i < bins; ++i) {
    float y = bin_dists[i].percentage(label);
    out->push_back(make_pair(x,y));
    x += increment;
  }
}



void RandomForest::reliability_diagram(int bins,
                                       vector<pair<float, float> >*out,
                                       vector<int>* count, int label) const {
  float increment = 1.0 / bins;
  float half = increment / 2.0;
  vector<DiscreteDist> bin_dists(bins);
  count->resize(bins, 0);
  for (int i = 0; i < set_.size(); ++i) {
    float prob = oob_predict_prob(i, label);
    int bin_no = int(floor(prob/increment));
    if (bin_no == bins) {
      bin_no = bins - 1;
    }
    bin_dists[bin_no].add(set_.label(i));
    (*count)[bin_no]++;
  }
  float x = half;
  for (int i = 0; i < bins; ++i) {
    float y = bin_dists[i].percentage(label);
    out->push_back(make_pair(x,y));
    x += increment;
  }
}

void RandomForest::compute_skewed_proximity(const InstanceSet& set,
                            vector<vector<float> >* prox,
                            int limit) const {
  for (int i = 0; i < trees_.size(); ++i) {
    trees_[i]->compute_skewed_proximity(set, prox, false, limit);
  }
  for (int i = 0; i < prox->size(); ++i) {
    for ( int j = 0; j < prox->size(); ++j) {
        (*prox)[i][j]/= trees_.size();
    }
  }
}



void RandomForest::compute_proximity(const InstanceSet& set,
                            vector<vector<float> >* prox,
                            int limit) const {
  for (int i = 0; i < trees_.size(); ++i) {
    trees_[i]->compute_proximity(set, prox, false, limit);
  }
  for (int i = 0; i < prox->size(); ++i) {
    for ( int j = 0; j < prox->size(); ++j) {
        (*prox)[i][j]/= trees_.size();
    }
  }
}

void RandomForest::compute_outliers(const InstanceSet& set, int label,
                                   const vector<vector<float> >& mat,
                                   vector< pair< float, int> >*ranking) const {

  vector<float> average_proximity(set.size(), 0.0);
  for (int i = 0; i < set.size(); ++i) {
    if (set.label(i) == label) {
      for (int j = 0; j < set.size(); ++j) {
        if ((i != j) && (set.label(j) == label)) {
          float prox = mat[i][j];
          average_proximity[i] += prox*prox;
        }
      }
    }
  }
  for (int i = 0; i <set.size(); ++i) {
    if (set.label(i) ==label) {
      float out_score = float(set.size())/average_proximity[i];
      ranking->push_back(make_pair(out_score, i));
    }
  }
  sort(ranking->begin(), ranking->end(), greater<pair<float,int> >());
}



float RandomForest::oob_predict_prob(int instance_no,
                                     int label) const {
   // Gather the votes from each tree
  DiscreteDist votes;
  int total = 0;
  for (int i = 0; i < trees_.size(); ++i) {
    if (trees_[i]->oob(instance_no)) {
      int predict = trees_[i]->predict(set_, instance_no);
      votes.add(predict);
      total++;
    }
  }
  return votes.percentage(label);
}

int RandomForest::oob_predict(int instance_no,
                          vector<pair<int, float> >*nodes) const {
  // Gather the votes from each tree
  DiscreteDist votes;
  for (int i = 0; i < trees_.size(); ++i) {
    if (trees_[i]->oob(instance_no)) {
      int predict = trees_[i]->predict(set_, instance_no, nodes);
      votes.add(predict);
    }
  }
  return votes.mode();
}




float RandomForest::predict_prob(const InstanceSet& set, int instance_no, int label) const {
  // Gather the votes from each tree
  DiscreteDist votes;
  for (int i = 0; i < trees_.size(); ++i) {
    int predict = trees_[i]->predict(set, instance_no);
    votes.add(predict);
  }
  return votes.percentage(label);
//  int count = votes.weight(label);
//  return float(count) / trees_.size();
}


void RandomForest::oob_predictions(vector<DiscreteDist>* predicts) const{
  predicts->resize(set_.size(),DiscreteDist(2));
  for (int i = 0; i < trees_.size(); ++i) {
    trees_[i]->oob_predictions(predicts);
  }
}


float RandomForest::oob_accuracy() const {
  vector<DiscreteDist> predicts;
  oob_predictions(&predicts);
  int total = 0;
  for (int i = 0; i < set_.size(); ++i) {
    if (predicts[i].mode() == set_.label(i)) {
      total++;
    }
  }
  return float(total)/set_.size();
}


void RandomForest::test_confusion(const InstanceSet& set) const {
  vector<int> predictions;
  vector<int> labels;
  for (int i = 0; i < set.size(); ++i) {
    labels.push_back(set.label(i));
    predictions.push_back(predict(set, i));
  }
  // HARDCODED
  // TODO: Fix me!
  confusion_matrix(2, predictions, labels);
}





void RandomForest::confusion_matrix(int num_labels,
                                    const vector<int>& predicts,
                                    const vector<int>& actual) const {
  assert(actual.size() == predicts.size());
  vector< vector<int> > matrix;
  matrix.resize(num_labels);
  for (int i = 0; i < num_labels; ++i) {
    matrix[i].resize(num_labels, 0);
  }
  for (int i = 0; i < actual.size(); ++i) {
    matrix[actual[i]][predicts[i]]++;
  }
  for (int i = 0; i < num_labels; ++i) {
    for (int j = 0; j < num_labels; ++j) {
      cout << matrix[i][j] << " ";
    }
    cout <<endl;
  }
}



void RandomForest::oob_confusion() const {
  vector<DiscreteDist> predicts;
  oob_predictions(&predicts);
  vector<int> prediction;
  vector<int> labels;
  for (int i = 0; i < predicts.size(); ++i) {
    prediction.push_back(predicts[i].mode());
    labels.push_back(set_.label(i));
  }
  //HARDCODED
  //TODO: fix me!
  confusion_matrix(2, prediction, labels);
}

/*
int RandomForest::predict(const Instance& c) const {
  // Gather the votes from each tree
  DiscreteDist votes;
  for (int i = 0; i < trees_.size(); ++i) {
    int predict = trees_[i]->predict(c);
    votes.add(predict);
  }
  return votes.mode();
}*/

float RandomForest::training_accuracy() const {
  int correct = 0;
  for (int i =0; i < set_.size(); ++i) {
    if (predict(set_, i) == set_.label(i))
      correct++;
  }
  return float(correct) / set_.size();
}

float RandomForest::testing_accuracy(const InstanceSet& set) const {
  int correct = 0;
  for (int i = 0; i < set.size(); ++i) {
    if (predict(set, i) == set.label(i))
      correct++;
  }
  return float(correct) / set.size();
}



void RandomForest::variable_importance(vector< pair< float, int> >*ranking,
                                       unsigned int* seed) const {
  vector<float> importances(set_.num_attributes(),
                            0.00);
  // Zero-out importances
  for (int i = 0; i < trees_.size(); ++i) {
    vector<float> tree_importance;
    trees_[i]->variable_importance(&tree_importance, seed);
    //aggregate
    for (int j = 0; j < tree_importance.size(); ++j) {
      importances[j] +=  tree_importance[j];
    }
  }
  // Get the mean of scores
  vector<float> raw_scores;
  float sum = 0;
  float sum_of_squares = 0;
  for (int i = 0; i < importances.size(); ++i) {
    float avg = importances[i] / trees_.size();
    assert(avg == avg);
    raw_scores.push_back(avg);
    sum += avg;
    sum_of_squares += (avg * avg);
  }
  float mean = sum / importances.size();
  assert(mean == mean);
  float std = sqrt(sum_of_squares/importances.size() - mean*mean);
  assert(std == std);

  // Write the z-scores
  for (int i = 0; i < raw_scores.size(); ++i) {
    float raw = raw_scores[i];
    float zscore = 0;
    if (std != 0) {
      zscore = (raw - mean) / std;
    }
    assert(zscore == zscore);
    ranking->push_back(make_pair(zscore, i));
  }
  // Sort
  sort(ranking->begin(), ranking->end(), greater<pair<float,int> >());
}


/*
float RandomForest::predict_prob(const Instance& c) const {
  // Gather the votes from each tree
  DiscreteDist votes;
  for (int i = 0; i < trees_.size(); ++i) {
    int predict = trees_[i]->predict(c);
    votes.add(predict);
  }
  int count = votes.weight(1);
  return float(count) / trees_.size();
}*/

} // namespace
