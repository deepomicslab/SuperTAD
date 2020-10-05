//
// Created by wang mengbo on 2019-09-01.
//

#include "params.h"


std::string SuperTAD::_WORK_DIR_ = "";
std::string SuperTAD::_INPUT_;
std::string SuperTAD::_OUTPUT_;
std::string SuperTAD::_RESULT_;
std::string SuperTAD::_RESULT_1_;
std::string SuperTAD::_RESULT_2_;

bool SuperTAD::_BINARY_ = false;
bool SuperTAD::_MULTI_ = false;
bool SuperTAD::_MULTI_H_ = false;
bool SuperTAD::_FILTER_ = false;
bool SuperTAD::_COMPARE_ = false;

int SuperTAD::_K_ = -999;
int SuperTAD::_MinSize_ = 2;
int SuperTAD::_optimalK_ = -999;
bool SuperTAD:: _DETERMINE_K_ = true;
int SuperTAD::_N_ = -999;
int SuperTAD::_H_ = 2;

// for multi_2d mode
int SuperTAD::_HU_ = 1;
int SuperTAD::_HD_ = 2;
std::string SuperTAD::_PRE_ = "";

float SuperTAD::_THRESHOLD_ = 1e-6;
bool SuperTAD::_PRE_LOG_ = true;

bool SuperTAD::_FILTERING_ = true;

bool SuperTAD::_FAST_ = true;
int SuperTAD::_PENALTY_ = -1;

bool SuperTAD::_VERBOSE_ = false;

// for matrix input
std::string SuperTAD::_CHROM1_ = "chr1";
std::string SuperTAD::_CHROM2_ = "chr1";
int64_t SuperTAD::_CHROM1_START_ = 0;
int64_t SuperTAD::_CHROM2_START_ = 0;
int SuperTAD::_RESOLUTION_ = 1;

// output format
bool SuperTAD::_BEDPE_ = false;
bool SuperTAD::_SHORT_ = false;
bool SuperTAD::_BIN_LIST_ = false;
bool SuperTAD::_NO_OUTPUT_ = false;

// debug
bool SuperTAD::_DEBUG_ = false;
std::string SuperTAD::_TMP_PATH_;
