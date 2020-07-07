//
// Created by wang mengbo on 2019-09-01.
//

#ifndef PROGRAM_PARAMS_H
#define PROGRAM_PARAMS_H

#include <string>


extern bool _BINARY_;
extern bool _MULTI_;
extern bool _MULTI_H_;

extern int _K_;
extern bool _DETERMINE_K_;
extern int _N_;
extern int _H_;

extern float _THRESHOLD_;  //avoid the precision issue
extern std::string _WORK_DIR_;
extern std::string _INPUT_;
extern std::string _OUTPUT_;
extern bool _FILTERING_;

extern bool _VERBOSE_;

extern std::string _TMP_PATH_;

extern bool _FAST_;
extern int _PENALTY_;

extern bool _DEBUG_;

// for matrix input
extern std::string _CHROM1_;
extern std::string _CHROM2_;
extern int64_t _CHROM1_START_;
extern int64_t _CHROM2_START_;
extern int _RESOLUTION_;

// for input in short bed format

// for output in bedpe
extern bool _BEDPE_;

#endif //PROGRAM_PARAMS_H
