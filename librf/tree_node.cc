/* tree_node.cc
 * 
 * The only implementation here is in reading/writing
 *
 */
#include "librf/tree_node.h"
namespace librf {
void tree_node::write(ostream& o) const {
  // we shouldn't be saving any other kind of node
  assert(status == TERMINAL || status == SPLIT);
  o << int(status);
  switch(status) {
    case TERMINAL:
      o << " " << int(label) << endl;
    break;
    case SPLIT:
      o << " " << left << " " << right << " " <<  attr << " " << split_point << endl;
    break;
  }
}

void tree_node::read(istream& i) {
  int status_int;
  i >> status_int;
  status = NodeStatusType(status_int);
  assert(status != EMPTY);
  switch(status) {
    case TERMINAL:
      int label_int;
      i >> label_int;
      label = uchar(label_int);
    break;
    case SPLIT:
      i >> left >> right >> attr >> split_point;
    break;
  }
}

} // namespace
