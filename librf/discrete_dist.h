/*******************
Discrete distribution
********************/
#ifndef _DISCRETE_DIST_H_
#define _DISCRETE_DIST_H_
#include <assert.h>
#include <vector>
#include <iostream>
#include <math.h>
#include "types.h"

using namespace std;

namespace librf {

class DiscreteDist {
  public:
    DiscreteDist(int size = 2) : sum_(0), size_(size),counter_(size, 0)
    {
      //counter_ = new unsigned int[size];
      //for (int i =0; i < size; ++i) {
      //  counter_[i] = 0;
      //}
    }
    /*~DiscreteDist() {
      //delete [] counter_;
    }*/
    void add(int value, unsigned int weight=1) {
      counter_[value] += weight;
      sum_ += weight;
    }
    void remove(int value, unsigned int weight=1) {
      counter_[value] -= weight;
      sum_ -= weight;
    }
    unsigned int sum() const {
      return sum_;
    }
    int mode() const {
      int max = -1;
      int mode = -10;
      for (int i = 0; i< size_; ++i) {
        int val = counter_[i];
        if (val > max) {
          max = val;
          mode = i;
        }
      }
      return mode;
    }
    void print() {
      for (int i = 0; i < size_; ++i) {
        cout << i << ":" << int(counter_[i]) << endl;
      }
    }
		unsigned int num_labels() const {
			return size_;
		}
		unsigned int weight(int i) const {
			return counter_[i];
		}
    float percentage(int i) const {
      return float(weight(i)) / sum();
    }
		static const double kLog2;
    static float entropy_conditioned_naive(const DiscreteDist* sets,
                                           int num_dists) {
      float H = 0;
      // H(Y |X) = Sum Prob(X=x) H(Y | x = x)
      float total = 0;
      for (int i = 0; i < num_dists; ++i) {
        float split_entropy = 0;
        float split_total = 0;
        for (int j = 0; j< sets[i].num_labels(); ++j) {
          float weight = sets[i].weight(j);
          split_entropy -= lnFunc(weight);
          split_total += weight;
          total += weight;
          cerr << j << ":" << weight <<endl;
        }
        if (split_total == 0) {
          split_entropy = 0;
        } else {
        split_entropy = (split_entropy + lnFunc(split_total) ) /
            (split_total *kLog2);
        }
        cerr << "Split " << i << ":" << split_entropy <<endl
          ;
        H += split_total * split_entropy;
      }
      return H / (total);
    }
		static float entropy_conditioned(const DiscreteDist* sets, int num_dists) {
			float returnValue = 0;
			float total = 0;
			float sumForSet;

			for (int i = 0; i < num_dists; ++i ) {
				sumForSet = 0;
				for (int j = 0; j < sets[i].num_labels(); ++j) {
					float weight = sets[i].weight(j);
					returnValue += lnFunc(weight);
					sumForSet += weight;
				}
				returnValue -= lnFunc(sumForSet);
				total += sumForSet;
			}
			if (total == 0){
				return 0;
			}
      returnValue = -returnValue /(total *kLog2);
      assert (returnValue == returnValue);
			return returnValue;
		}

		// Adapted from ContingencyTables.java: entropyOverColumns
		float entropy_over_classes() const{
			float returnValue = 0;
			float total = 0;
			for (int i = 0; i < size_; ++i) {
				returnValue -= lnFunc(counter_[i]);
				total += counter_[i];
			}
			if (total == 0) {
				return 0;
			}
			return (returnValue + lnFunc(total)) / (total * kLog2);
		}
  private:
    unsigned int sum_;
    unsigned int size_;
  static float lnFunc(float num) {
			if (num  < 1e-6) {
				return 0;
			} else {
				return num * log(num);
			}
		}
    vector<unsigned int> counter_;
    //unsigned int* counter_;
};
} // namespace
#endif
