//
// Created by wang mengbo on 2019-09-01.
//

#ifndef PROGRAM_PARAMS_H
#define PROGRAM_PARAMS_H

#include <string>


namespace SuperTAD
{
    extern std::string _WORK_DIR_;
    extern std::string _INPUT_;
    extern std::string _OUTPUT_;
    extern std::string _RESULT_;
    extern std::string _RESULT_1_;
    extern std::string _RESULT_2_;
    extern bool _BINARY_;
    extern bool _MULTI_;
    extern bool _MULTI_H_;
    extern bool _FILTER_;
    extern bool _COMPARE_;
    extern int _K_;
    extern int _MinSize_;
    extern int _optimalK_;
    extern bool _DETERMINE_K_;
    extern int _N_;
    extern int _H_;
    // for multi_2d mode
    extern int _HU_;
    extern int _HD_;
    extern std::string _PRE_;
    extern float _THRESHOLD_;
    extern bool _PRE_LOG_;
    extern bool _FILTERING_;
    extern bool _FAST_;
    extern int _PENALTY_;
    extern bool _VERBOSE_;
    // for matrix input
    extern std::string _CHROM1_;
    extern std::string _CHROM2_;
    extern int64_t _CHROM1_START_;
    extern int64_t _CHROM2_START_;
    extern int _RESOLUTION_;
    // output format
    extern bool _BEDPE_;
    extern bool _SHORT_;
    extern bool _BIN_LIST_;
    // debug
    extern bool _DEBUG_;
    extern std::string _TMP_PATH_;
}

//avoid the precision issue

// for input in short bed format

#endif //PROGRAM_PARAMS_H
