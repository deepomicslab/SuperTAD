//
// Created by wang mengbo on 2019-09-01.
//

#include "params.h"


namespace SuperTAD
{
    std::string _WORK_DIR_;
    std::string _INPUT_;
    std::string _OUTPUT_;
    std::string _PRE_RESULT_;    // for filter and multi_2d mode
    std::string _RESULT_1_;
    std::string _RESULT_2_;
    std::string _CHROM1_ = "chr1";
    std::string _CHROM2_ = "chr1";

    bool _BINARY_ = false;
    bool _MULTI_ = false;
    bool _MULTI_H_ = false;
    bool _FILTER_ = false;
    bool _COMPARE_ = false;
    bool _DEEPBINARY_ = false;
    bool _PRUNE_ = true;
    bool _FAST_ = true; // if multi mode, default run SuperTAD-Fast
    bool _VERBOSE_ = false;
    bool _SPARSE_ = false;
    bool _BEDPE_ = false;
    bool _SHORT_ = false;
    bool _BIN_LIST_ = false;
    bool _DETERMINE_K_ = true;  // if version1.0, default automatically select K
    bool _TURBO_PRUNE_ = false;
    bool _V1_ = false;

    int _HU_ = 1;
    int _HD_ = 2;
    int _PENALTY_ = -1;
    int _PRUNE_METHOD_ = 0;
    int _K_ = -999;
    int _MAX_SIZE_ = -999;
    int _OPTIMAL_K_ = -999;
    int _N_ = -999;
    int _H_ = 2;
    int _RESOLUTION_ = 1;
    double _THRESHOLD_ = 1e-6;
    int _STEP_= -999;
    int _WINDOW_ = 5;
    double _BAYESFACTOR_ = 1;   // for bayesian corrected version, default 1

    int64_t _CHROM1_START_ = 0;
    int64_t _CHROM2_START_ = 0;

// debug
    bool _DEBUG_ = false;
    std::string _TMP_PATH_;

}