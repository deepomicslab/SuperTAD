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
    extern std::string _PRE_RESULT_;    // for filter and multi_2d mode
    extern std::string _RESULT_1_;  // for compare mode
    extern std::string _RESULT_2_;
    extern std::string _CHROM1_;
    extern std::string _CHROM2_;

    extern bool _BINARY_;
    extern bool _MULTI_;
    extern bool _MULTI_H_;
    extern bool _FILTER_;
    extern bool _COMPARE_;
    extern bool _DEEPBINARY_;
    extern bool _DETERMINE_K_;
    extern bool _PRUNE_;    // optimal, nut subcommand
    extern bool _FAST_; // for SuperTAD-Fast
    extern bool _VERBOSE_;
    extern bool _SPARSE_;   // for Bayesian approach
    extern bool _BEDPE_;
    extern bool _SHORT_;
    extern bool _BIN_LIST_;
    extern bool _TURBO_PRUNE_;  // unknown command
    extern bool _V1_;   // run version1 where K is kept

    extern int _HU_;
    extern int _HD_;
    extern int _PENALTY_;   // under testing ###
    extern int _PRUNE_METHOD_;  // for deepbinary mode
    extern int _K_;
    extern int _MAX_SIZE_;  // the max size of domain, default: 10 Mb
    extern int _OPTIMAL_K_; // for version1.0
    extern int _N_;
    extern int _H_;
    extern int _RESOLUTION_;
    extern int _WINDOW_;
    extern int _STEP_;
    extern double _BAYESFACTOR_;

    extern double _THRESHOLD_;

    extern int64_t _CHROM1_START_;
    extern int64_t _CHROM2_START_;

    // debug
    extern bool _DEBUG_;
    extern std::string _TMP_PATH_;
}

#endif //PROGRAM_PARAMS_H
