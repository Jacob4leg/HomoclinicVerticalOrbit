#include<iostream>
#include<tuple>
#include "capd/capdlib.h"
#include "cr3bp.h"
#include "linalg_helper.h"

#ifndef _PROVE_FIXED_POINT_H_
#define _PROVE_FIXED_POINT_H_

bool prove_fixed_point(capd::LDVector v_E, capd::Interval X, capd::Interval Y, capd::Interval E, const capd::LDMatrix &T);
// capd::Interval cone_coeff(long double eps_x, long double eps_y, capd::Interval E, const capd::LDMatrix &T);
void eval_rectangle(long double eps_x, long double eps_y, capd::Interval E);
#endif