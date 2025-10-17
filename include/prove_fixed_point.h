#include<iostream>
#include "capd/capdlib.h"
#include "cr3bp.h"
#include "linalg_helper.h"

#ifndef _PROVE_FIXED_POINT_H_
#define _PROVE_FIXED_POINT_H_

bool prove_fixed_point(long double eps_x, long double eps_y, long double E_eps);
capd::Interval cone_coeff(long double eps_x, long double eps_y, capd::Interval E);
void eval_rectangle(long double eps_x, long double eps_y, capd::Interval E);
void test(capd::Interval E);

#endif