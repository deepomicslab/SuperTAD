//
// Created by wang mengbo on 2019-09-01.
//

#include "params.h"


namespace SuperTAD
{
    std::string _WORK_DIR_;
    std::string _INPUT_;
    std::string _OUTPUT_;
    std::string _RESULT_;
    std::string _RESULT_1_;
    std::string _RESULT_2_;
    std::string _SE_RESULT_PATH_;
    std::string _CHROM1_ = "chr1";
    std::string _CHROM2_ = "chr1";
    std::string _PRE_;

    bool _BINARY_ = false;
    bool _MULTI_ = false;
    bool _MULTI_H_ = false;
    bool _FILTER_ = false;
    bool _COMPARE_ = false;
    bool _DEEPBINARY_ = false;
    bool _PRUNE_ = false;
    bool _PRE_LOG_ = true;
    bool _FILTERING_ = true;
    bool _FAST_ = true;
    bool _VERBOSE_ = false;
    bool _SPARSE_ = false;
    bool _BEDPE_ = false;
    bool _SHORT_ = false;
    bool _BIN_LIST_ = false;
    bool _NO_OUTPUT_ = false;
    bool _DETERMINE_K_ = true;
    bool _APPEND_RESULT_ = false;
    bool _TURBO_PRUNE_ = false;

    int _K_ = -999;
    int _MinSize_ = 2;
    int _optimalK_ = -999;
    int _N_ = -999;
    int _H_ = 2;
    int _HU_ = 1;
    int _HD_ = 2;
    int _PENALTY_ = -1;
    int _RESOLUTION_ = 1;
    int _PRUNE_METHOD_ = 0;
    double _THRESHOLD_ = 1e-13;

    int64_t _CHROM1_START_ = 0;
    int64_t _CHROM2_START_ = 0;

// debug
    bool _DEBUG_ = false;
    std::string _TMP_PATH_;

}