/**
 * tree_node.h
 * @file
 * @brief structure for a single node in a decision tree 
 */

#ifndef _TREE_NODE_H_
#include "librf/types.h"
#include <fstream>
using namespace std;

namespace librf {

typedef enum {EMPTY, BUILD_ME, TERMINAL, SPLIT} NodeStatusType;

struct tree_node {
/// default constructor initializes to garbage known state
  tree_node() :status(EMPTY),
               label(99),
               attr(99),
               start(99),
               size(99),
               split_point(-999.0),
               entropy(-9999.0),
               left(0),
               right(0){}
  NodeStatusType status;
  uchar label;
  uint16 attr;
  uint16 start;
  uint16 size;
  uint16 left;
  uint16 right;
  float entropy;
  float split_point;
  uchar depth;

  void write(ostream& o) const;
  void read(istream& i);
};
} //namespace
#endif
