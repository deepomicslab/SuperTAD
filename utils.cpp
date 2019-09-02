//
// Created by wang mengbo on 2019-09-02.
//

#include "utils.h"

namespace utils {
  // including low but excluding high
  int randInt (int low, int high) {
    double d = (double) rand () / (double) RAND_MAX * (high - low);
    double intpart;
    if (modf (d, &intpart) > 0.5) {
      int r = ceil (d);
      if (r >= high)
        return high - 1;
      return r;
    }
    else
      return floor (d);
  }
  
  
  double randDouble (double low, double high) {
    return (double) rand () / (double) RAND_MAX * (high - low);
  }
}