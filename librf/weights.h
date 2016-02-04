/* weights.h
 * Instance Weights:
 *
 * this is trying to be a memory "smart" container that will decide
 * whether or not to use sparse/dense representation yet
 * still provide a uniform interface.
 *
 * It is possible to do this because - we know how many training 
 * instances there are apriori - and each time we split a node
 * we will know how many instances are going to the left/right
 * subtrees
 * 
 * approx mem usage:
 *  sparse: density * 5 bytes + map overhead
 *  dense: num_instances bytes
 */
#ifndef _WEIGHTS_H_
#define _WEIGHTS_H_

#include <vector>
#include <iostream>
#include "librf/types.h"

using namespace std;

namespace librf {

/*
using ::__gnu_cxx::hash_map;
class weight_interface {
  public:
    weight_interface(int max_size) {
      sum_ = 0;
    }
    unsigned int sum() {
       return sum_;
    }
    void add_instance(int instance, byte num =1) {
      sum_ += num;
      add_instance_impl(instance, num);
    }
    byte get_weight(int instance) const{
      return get_weight_impl(instance);
    }
  private:
    unsigned int sum_;
    virtual void add_instance_impl(int instance, byte num) = 0;
    virtual byte get_weight_impl(int instance) const = 0;
};


// need a smart iterator
class sparse_weight : public weight_interface {
 public:
 sparse_weight(int size) : weight_interface(size) {
 }
 private:
    virtual void add_instance_impl(int i, byte num) {
      byte_map::iterator it = store_.find(i);
      if (it != store_.end()) {
        it->second += num;
      } else {
        store_[i] = num;
      }
    }
    virtual byte get_weight_impl(int i) const {
      byte_map::const_iterator it = store_.find(i);
      if (it != store_.end()) {
        return it->second;
      }
      return 0;
    }
    typedef google::dense_hash_map<int, byte> byte_map;
    byte_map store_;
};

class dense_weight : public weight_interface {
  public:
  dense_weight(int max_size) : weight_interface(max_size),
                               store_(max_size,0) {
  }
  private:
    // virtual init_impl(int max_size) {
    //  store_.resize(max_size, 0);
    // }
    virtual void add_instance_impl(int i, byte num) {
      // cout << int(this) <<" dense add " << i << "," << int(num) << endl;
      store_[i] += num;
    }
    virtual byte get_weight_impl(int i) const {
      return  store_[i];
    }
    vector<byte> store_;
};

/*
class weight_list {
  public:
    weight_list(int n, int density) : num_instances_(n){
      int sparse_cost = density * 5;
      if (sparse_cost < num_instances_) {
        // cout << "creating sparse" << endl;
        impl = new sparse_weight(num_instances_);
      } else {
        // cout << "creating dense " << density << "/" << num_instances_ << endl;
        impl = new dense_weight(num_instances_);
      }
    }
    byte operator[](int i) const{
      return impl->get_weight(i);
    }
    void add(int i, byte num = 1) {
      impl->add_instance(i, num);
    }
    int size() const {
      return num_instances_;
    }
    int sum() const {
      return impl->sum();
    }
    ~weight_list() {
       delete impl;
    }
  private:
    int num_instances_;
    weight_interface *impl;
};*/

class weight_list {
  public:
   weight_list(int n, int density) : num_instances_(n), sum_(0){
      array_ = new byte[n];
      for(int i =0; i < n; ++i) {
        array_[i] = 0;
      }
   }
   byte operator[](int i) const {
    return array_[i];
   }
   void add(int i, byte num =1){
    array_[i]+=num;
    sum_+=num;
   }
   int size() const { return num_instances_;}
   int sum() const {
      return sum_;
   }
   ~weight_list(){
     delete[] array_;
   }
  private:
   int sum_;
   int num_instances_;
   byte* array_;
};

} // namespace
#endif
