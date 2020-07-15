//
// Created by wang mengbo on 2019-09-01.
//

#include "params.h"


std::string _WORK_DIR_ = "";
std::string _INPUT_;
std::string _OUTPUT_;

bool _BINARY_ = true;
bool _MULTI_ = false;
bool _MULTI_H_ = false;

int _K_ = -999;
bool _DETERMINE_K_ = true;
int _N_ = -999;
int _H_ = 2;

float _THRESHOLD_ = 1e-6;
bool _PRE_LOG_ = true;

bool _FILTERING_ = true;

bool _FAST_ = true;
int _PENALTY_ = -1;

bool _VERBOSE_ = false;

// for matrix input
std::string _CHROM1_ = "1";
std::string _CHROM2_ = "1";
int64_t _CHROM1_START_ = 0;
int64_t _CHROM2_START_ = 0;
int _RESOLUTION_ = 1;

bool _BEDPE_ = false;

// debug
bool _DEBUG_ = false;
bool _TEST_LOG2_TIME_ = false;
bool _TEST_FAST_ = false;
std::string _TMP_PATH_;
