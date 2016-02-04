/**
 * @file
 * @brief Instance Set
 * This is the abstraction for a data set
 * -- Currently libSVM, CSV
 * -- want to support ARFF
 */
#ifndef _INSTANCE_SET_H_
#define _INSTANCE_SET_H_

#include <string>
#include <vector>
#include <fstream>
#include "discrete_dist.h"

using namespace std;

namespace librf {

class weight_list;
/**
 * @brief
 * InstanceSet class. Interface for loading training/testing data.
 */
class InstanceSet {
    public:
        /// Empty constructor
        InstanceSet();
        /// Named constructor - create a subset from an existing instanceset
        static InstanceSet* create_subset(const InstanceSet&, const weight_list&);
        /// Named constructor - no labels
        static InstanceSet* load_unsupervised(const string& data,
                                              unsigned int* seed,
                                              bool header = false,
                                              const string& delim =",");
        /// Named constructor - feature selection
        static InstanceSet* feature_select(const InstanceSet&, const vector<int>&);
        /// Named constructor - load from csv file and a label file
        static InstanceSet* load_csv_and_labels(const string& data,
                                                    const string& labels,
                                                    bool header = false,
                                                    const string& delim =",");
        /// copy a variable array out 
        void save_var(int var, vector<float> *target);
        /// load a variable array in
        void load_var(int var, const vector<float>&);
        /// permute a variable's instances (shuffle)
        void permute(int var, unsigned int * seed);
        /// sort the variables 
        void create_sorted_indices();
        /// Sorted indices (available after create_sorted_indices)
        const vector<int>& get_sorted_indices(int attribute) const{
            return sorted_indices_[attribute];
        }
        /// Most common label
        int mode_label() const {
          return distribution_.mode();
        }
        /// Get a particular instance's label
        unsigned char label(int i) const{
          return labels_[i];
        }
        /// Number of instances
        unsigned int size() const {
            return labels_.size();
        }
        /// Number of attributes
        unsigned int num_attributes() const {
          return attributes_.size();
        }
        /// Get a particular instance's attribute
        float get_attribute(int i, int attr) const {
          return attributes_[attr][i];
        }
        /// Get a variable name (useful if there is a header with var
        //names)
        string get_varname(int i) const {
          return var_names_[i];
        }
        void write_csv(ostream& out, bool header, const string& delim);
        void write_transposed_csv(ostream& out, const string& delim);
        //float class_entropy() const{
        //  return distribution_.entropy_over_classes();
        //}
    private:
        /// Load from csv file and labels
        InstanceSet(const string& csv_data, const string& labels,
                    bool header=false, const string& delim=",");
        // Load csv for unsupervised
        InstanceSet(const string& csv_data, unsigned int* seed,
                    bool header=false, const string& delim=",");
        /// Load from libsvm format
        InstanceSet(const string& filename, int num);
        /// Get a subset of an existing instance set
        InstanceSet(const InstanceSet&, const weight_list&);
        /// Feature select from existing instance set
        InstanceSet(const InstanceSet&, const vector<int>&);
        void load_labels(istream& in);
        void load_csv(istream& in, bool header, const string& delim);
        void load_svm(istream& in);
        void create_dummy_var_names(int n);
        void sort_attribute(const vector<float>&attribute, vector<int>*indices);
        DiscreteDist distribution_;
        // List of Attribute Lists
        // Thus access is attributes_ [attribute] [ instance]
        vector< vector<float> > attributes_;
        // List of true labels
        // access is labels_ [instance]
        vector<unsigned char> labels_;
        vector<string> var_names_;
        vector< vector<int> > sorted_indices_;
};

}  // namespace
#endif
