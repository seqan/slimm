#include "librf/instance_set.h"
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <float.h>
#include "librf/weights.h"
#include "librf/types.h"
#include "librf/stringutils.h"

using namespace std;

namespace librf {
InstanceSet::InstanceSet(){}
/***
 * Named constructor for loading from a csv file and a label file
 * Makes simpler to have a separate label file.
 * @param csv_data CSV filename
 * @param header whether there is a header with var. names
 * @param delim CSV delimiter - defaults to ','
 */
InstanceSet* InstanceSet::load_csv_and_labels(const string& csv_data,
                                      const string& label_file,
                                      bool header,
                                      const string& delim) {
  return new InstanceSet(csv_data, label_file, header, delim);
}

/***
 * Named constructor for loading unsupervised
 * @param csv_data CSV filename
 * @param header whether there is a header with var. names
 * @param delim CSV delimiter - defaults to ','
 */
InstanceSet* InstanceSet::load_unsupervised(const string& csv_data,
                                              unsigned int * seed,
                                      bool header,
                                      const string& delim) {
  return new InstanceSet(csv_data, seed, header, delim);
}




/**
 * Named constructor for creating a subset from an existing set
 * 
 * @param set existing set 
 * @param wl non-zeroes weights
 */
InstanceSet* InstanceSet::create_subset(const InstanceSet& set,
                                     const weight_list& wl) {
  return new InstanceSet(set, wl);
}
/***
 * Unnamed private constructor for loading csv file/label file
 */
InstanceSet::InstanceSet(const string& csv_data,
            const string& label_file,
            bool header, const string& delim) {
  ifstream data(csv_data.c_str());
  load_csv(data, header, delim);

  ifstream labels(label_file.c_str());
  load_labels(labels);
  create_sorted_indices();
  assert(attributes_.size() > 0);
  assert(attributes_[0].size() == labels_.size());
}

/***
 * Unnamed private constructor for loading CSV for unsupervised
 */
InstanceSet::InstanceSet(const string& csv_data, unsigned int* seed,
                         bool header, const string& delim) {
  ifstream data(csv_data.c_str());
  load_csv(data, header, delim);
  // organic set gets 0 label
  assert(attributes_.size() > 0);
  labels_.resize(attributes_[0].size(), 0);
  int original_set_size = attributes_[0].size();
  for (int i = 0; i < original_set_size; ++i) {
    for (int j = 0; j < attributes_.size(); ++j) {
      // uniformly sample from attributes_[j][0-original_set_size-1]
      int select = rand_r(seed) % original_set_size;
      float sample = attributes_[j][select];
      attributes_[j].push_back(sample);
    }
    // synthetic gets 1 label
    labels_.push_back(1);
  }
  assert(attributes_[0].size() == labels_.size());
  create_sorted_indices();
}

/**
 * Named constructor for feature selection 
 * @param set existing InstanceSet
 * @param attrs attributes/features wanted 
 */
InstanceSet* InstanceSet::feature_select(const InstanceSet& set,
                                      const vector<int>& attrs) {
  return new InstanceSet(set, attrs);
}

/**
 * Private constructor for feature selection
 */
InstanceSet::InstanceSet(const InstanceSet& set,
                         const vector<int>& attrs) {
  // Copy labels
  labels_ = set.labels_;
  // Only copy given attrs
  attributes_.resize(attrs.size());
  sorted_indices_.resize(attrs.size());
  var_names_.resize(attrs.size());
  for (int i = 0; i < attrs.size(); ++i) {
    attributes_[i] = set.attributes_[attrs[i]];
    sorted_indices_[i] = set.sorted_indices_[attrs[i]];
    var_names_[i] = set.var_names_[attrs[i]];
  }
}
/***
 * Load labels from an istream
 *
 * LIMITATION: assumes binary labels... this needs to change
 */
void InstanceSet::load_labels(istream&in) {
  float label;
  int true_label;
  while (in >> label) {
    if (label == -1.0) {
      true_label =0;
    } else if (label ==0.0) {
      true_label =0;
    } else if (label ==1.0) {
      true_label =1;
    } else {
      cerr << "Incorrect label (only +1, 0, -1 supported)" << endl;
      assert(false);
    }
    labels_.push_back(true_label);
    distribution_.add(true_label);
  }
}

/**
 *
 */
void InstanceSet::write_csv(ostream& out, bool header,
                            const string& delim) {
  if (header) {
    for (int i = 0; i < num_attributes() - 1; ++i) {
      out <<  var_names_[i] << delim;
    }
    out << var_names_[num_attributes() - 1] << endl;
  }
  assert(attributes_.size() > 0);
  assert(attributes_[0].size() > 0);
  for (int i = 0; i < attributes_[0].size(); ++i) {
    for (int j = 0; j < num_attributes() - 1; ++j) {
      out << attributes_[j][i] << delim;
    }
    out << attributes_[num_attributes() - 1][i] << endl;
  }
}

void InstanceSet::write_transposed_csv(ostream& out,
                            const string& delim) {
  assert(attributes_.size() > 0);
  assert(attributes_[0].size() > 0);
  for (int j = 0; j < num_attributes(); ++j) {
    for (int i = 0; i < attributes_[0].size() - 1; ++i) {
      out << attributes_[j][i] << delim;
    }
    out << attributes_[j][attributes_[0].size() - 1] << endl;
  }
}




/***
 * Load csv from an istream
 */
void InstanceSet::load_csv(istream&in, bool header, const string& delim) {
  // read variable names
  int num_features = -1;
  string buffer;
  getline(in, buffer);
  if (header) {
    StringUtils::split(buffer, &var_names_, delim);
    num_features = var_names_.size();
  } else {
    vector<string> ary;
    StringUtils::split(buffer, &ary, delim);
    num_features = ary.size();
    // reset istream
    in.seekg (0, ios::beg);
    // create dummy var names
  }
  create_dummy_var_names(num_features);
  attributes_.resize(num_features);
  // get all data
  while(getline(in, buffer)) {
    vector<string> ary;
    StringUtils::split(buffer, &ary, delim);
    assert(ary.size() == num_features);
    for (int i = 0; i < ary.size(); ++i) {
      stringstream ss(ary[i]);
      //convert to float
      float val;
      ss >> val;
      attributes_[i].push_back(val);
    }
  }
}

void InstanceSet::create_dummy_var_names(int n) {
  for (int i = 0; i <  n; ++i) {
    stringstream ss;
    ss << i;
    var_names_.push_back(ss.str());
  }
}
void InstanceSet::create_sorted_indices() {
    // allocate sorted_indices_
    sorted_indices_.resize(attributes_.size());
    // sort 
    for (int i = 0; i < attributes_.size(); ++i) {
        sort_attribute(attributes_[i], &sorted_indices_[i]);
    }
}

void InstanceSet::sort_attribute(const vector<float>& attribute,
                                 vector<int>*indices) {
    vector<pair<float, int> > pairs;
    for (int i = 0; i < attribute.size(); ++i) {
        pairs.push_back(make_pair(attribute[i],i));
    }
    sort(pairs.begin(), pairs.end());
    for (int i = 0; i < pairs.size(); ++i) {
        indices->push_back(pairs[i].second);
    }
}

// Grab a subset of the instance (for getting OOB data
InstanceSet::InstanceSet(const InstanceSet& set,
                         const weight_list& weights) : attributes_(set.num_attributes()){
  // Calculate the number of OOB cases
  //cout << "creating OOB subset for weight list of size "
  //     << weights.size() << endl;
  for (int i = 0; i < weights.size(); ++i) {
    if (weights[i] == 0) {
      // append instance
      for (int j = 0; j < set.num_attributes(); ++j) {
          attributes_[j].push_back(set.get_attribute(i, j));
      }
      labels_.push_back(set.label(i));
    }
  }
}


/**
 * Permute method
 * Used for variable importance
 */
void InstanceSet::permute(int var, unsigned int *seed) {
  vector<float>& attr = attributes_[var];
  for (int i = 0; i < attr.size(); ++i) {
    int idx = rand_r(seed) % labels_.size(); // randomly select an index
    float tmp = attr[i];  // swap last value with random index value
    attr[i] = attr[idx];
    attr[idx] =  tmp;
  }
}

void InstanceSet::load_var(int var, const vector<float>& source) {
  // use the STL built-in copy/assignment
  attributes_[var] = source;
}

void InstanceSet::save_var(int var, vector<float>* target) {
  const vector<float>& attr = attributes_[var];
  // use the STL built-in copy/assignment
  *target = attr;
}

} // namespace
