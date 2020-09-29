//
// Created by wang mengbo on 2019-09-01.
//

#include "params.h"


std::string _WORK_DIR_ = "";
std::string _INPUT_;
std::string _OUTPUT_;
std::string _RESULT_;
std::string _RESULT_1_;
std::string _RESULT_2_;

bool _BINARY_ = false;
bool _MULTI_ = false;
bool _MULTI_H_ = false;
bool _FILTER_ = false;
bool _COMPARE_ = false;

int _K_ = -999;
int _optimalK_ = -999;
bool _DETERMINE_K_ = true;
int _N_ = -999;
int _H_ = 2;

// for multi_2d mode
int _HU_ = 1;
int _HD_ = 2;
std::string _PRE_ = "";

float _THRESHOLD_ = 1e-6;
bool _PRE_LOG_ = true;

bool _FILTERING_ = true;

bool _FAST_ = true;
int _PENALTY_ = -1;

bool _VERBOSE_ = false;

// for matrix input
std::string _CHROM1_ = "chr1";
std::string _CHROM2_ = "chr1";
int64_t _CHROM1_START_ = 0;
int64_t _CHROM2_START_ = 0;
int _RESOLUTION_ = 1;

// output format
bool _BEDPE_ = false;
bool _SHORT_ = false;
bool _BIN_LIST_ = false;

// debug
bool _DEBUG_ = false;
std::string _TMP_PATH_;
