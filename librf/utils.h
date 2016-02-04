
#ifndef _UTILS_H_
#define _UTILS_H_

#include <vector>
#include <stdlib.h>

namespace librf {
/// Sample without replacment
// Fix provided by naoki.tanaka@naoki.tanaka@jp.astellas.com
void random_sample(int n, int K, vector<int>*v, unsigned int* seed) {
  if (K < n) {
  int pop = n;
  v->reserve(K);
  for (int i = K; i > 0; --i) {
    float cumprob = 1.0;
    float x;
    do {
      x = float(rand_r(seed))/RAND_MAX;
    } while (x == 1);
    for (; x < cumprob; pop--) {
      cumprob -= cumprob * i /pop;
    }
    v->push_back(n - pop - 1);
  }
  } else {
    for (int i =0; i < n; i++) {
      v->push_back(i);
    }
  }
}
// slow and stupid median
/*
float median(const vector<float>& list) {
  vector<float>  sorted = list;
  sort(sorted.begin(), sorted.end());
  float half = sorted.size() / 2;
  if ( sorted.size()&1 ==0) { //even case
    return (sorted[half] + sorted[half - 1])/2.0;
  } else {
    return (sorted[int(half)]);
  }
}
*/
} // namespace
#endif
