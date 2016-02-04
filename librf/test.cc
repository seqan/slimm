#include <vector>
#include <iostream>

using namespace std;

struct tree_node {
  int node;
};
int main(int argc, char*argv[]) {
  vector<tree_node> v;
  tree_node a;
  v.push_back(a);
  tree_node &ref = v[0];
  ref.node = 23;
  cout << v[0].node <<endl;
}
