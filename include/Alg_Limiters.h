#ifndef ALG_LIMITERS_H
#define ALG_LIMITERS_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include "Data_Type.h"
#include "LowDim_Alg.h"
using namespace std;

#ifndef K_VANLEER_MAX
#define K_VANLEER_MAX 10
#endif

#ifndef K_VANALBADA_MAX
#define K_VANALBADA_MAX 10
#endif

/* This library calculate the 1D SLOPE limiter rather than the 
 * mathematically equivalent FLUX limiter. There are two reasons.
 * 
 * The first reason is: 
 * On non-uniform 1D grids, the flux limiters are calculated at the cell faces,
 * and it will depend on the local characteristic velocity direction (upwinding).
 * However, the slope limiters are calculated at the cell centers, and it turns
 * out that the slope limiter has unified expressions for both left going
 * and right going characteristics. 
 * Hence though equivalent, slope limiter gives simpler expressions than flux
 * limiters.
 *
 * The second reason is:
 * For FVM, data is always assumed to be stored at the center of the cell, 
 * hence the data slopes involved in linear reconstruction in each cell should 
 * also be a cell-information rather than nodal (cell-faces) information.
 */

/* For each limiter, there are three versions: scalar, VEC3 and VEC4.
 * They are built from the same principal, but with different interfaces
 * for convenience.
 */

/* This library only deals with algebraic form of limiters, no fluid 
 * is considered.
 * But it does assume that Finite Volumn Schemes are used.
 */

/* For the following limiters, both uniform-mesh and non-uniform-mesh versions
 * are available:
 *
 * Minmod
 * Superbee
 * MC
 *
 * For the following limiters, only uniform-mesh version is implemented:
 *
 * van Leer
 * van Albada
 */

/* Remark:
 * Typically, Minmod and Superbee are not recommended, since the L_inf error
 * may not converge in second order accuracy when they are used.
 * This is mainly due to the fact that when these two limiters are used, 
 * at the corner, when mesh is refined, not all the data are converged in the 
 * same direction (i.e. "< convergent value" or "> convergent value")
 */

/* General rules:
 * Each function is assumed to calculate the limited slope of data for a particular
 * cell, i.e. U_j (the index j is not given to the function).
 *
 * In order to calculate the second-order limiter (i.e. linear reconstruction),
 * U_{j-1}, U_j, U_{j+1} are needed, which bear the names:
 * U_l, U_c, U_r in all function calls.
 *
 * For limiters on non-uniform grid, the three cell sizes:
 * dx_{j-1}, dx_j, dx_{j+1} are also required, which have the names:
 * dx_l, dx_c, dx_r
 *
 * CAUTION: They always appear in the input list as:
 * center, left, right sequence!
 *
 * Typically the smoothness monitor theta_j for cell j is calculated as:
 *			 U_j-U_{j-1}
 * theta_j = ___________
 * 			 U_{j+1}-U_j
 * However, to avoid the most smoothest possible data (constant data), the value
 * of limiter function phi_j(theta_j) is not explicitly calculated.
 * Instead, the limited slope is calculated and returned:
 * 			  1
 * sigma_j = ____ phi_j(theta_j)(U_{j+1}-U_j)
 * 			 dx_j
 *
 * For limiters implemented for both uniform and non-uniform meshes, the uniform case 
 * is considered as default.
 */

// No limiter: linear reconstruction with no limiter
// three coefficients: coefficients used to calculate
// the slope from left / right / or central difference
// On uniform grids, the optimized configuration for order
// is _L = 1.0; _R = 1.0; _C = -1.0
#define ALG_LIMSLOPE_LINEAR_L 0.0
#define ALG_LIMSLOPE_LINEAR_R 1.0
#define ALG_LIMSLOPE_LINEAR_C (1.0-ALG_LIMSLOPE_LINEAR_L-ALG_LIMSLOPE_LINEAR_R)
DType Alg_LimSlope_Linear(DType U_c, DType U_l, DType U_r, DType dx_c, DType dx_l=-1, DType dx_r=-1);
VEC3 Alg_LimSlope_Linear(VEC3 U_c, VEC3 U_l, VEC3 U_r, DType dx_c, DType dx_l=-1, DType dx_r=-1);
VEC4 Alg_LimSlope_Linear(VEC4 U_c, VEC4 U_l, VEC4 U_r, DType dx_c, DType dx_l=-1, DType dx_r=-1);

// Minmod limiter: non-uniform mesh implemented
DType Alg_LimSlope_Minmod(DType U_c, DType U_l, DType U_r, DType dx_c, DType dx_l=-1, DType dx_r=-1, DType coef_a=-1, DType coef_b=-1);
VEC3 Alg_LimSlope_Minmod(VEC3 U_c, VEC3 U_l, VEC3 U_r, DType dx_c, DType dx_l=-1, DType dx_r=-1, DType coef_a=-1, DType coef_b=-1);
VEC4 Alg_LimSlope_Minmod(VEC4 U_c, VEC4 U_l, VEC4 U_r, DType dx_c, DType dx_l=-1, DType dx_r=-1, DType coef_a=-1, DType coef_b=-1);

// Superbee limiter: non-uniform mesh implemented
DType Alg_LimSlope_Superbee(DType U_c, DType U_l, DType U_r, DType dx_c, DType dx_l=-1, DType dx_r=-1, DType coef_a=-1, DType coef_b=-1);
VEC3 Alg_LimSlope_Superbee(VEC3 U_c, VEC3 U_l, VEC3 U_r, DType dx_c, DType dx_l=-1, DType dx_r=-1, DType coef_a=-1, DType coef_b=-1);
VEC4 Alg_LimSlope_Superbee(VEC4 U_c, VEC4 U_l, VEC4 U_r, DType dx_c, DType dx_l=-1, DType dx_r=-1, DType coef_a=-1, DType coef_b=-1);

// MC limiter: non-uniform mesh implemented
DType Alg_LimSlope_MC(DType U_c, DType U_l, DType U_r, DType dx_c, DType dx_l=-1, DType dx_r=-1, DType coef_a=-1, DType coef_b=-1);
VEC3 Alg_LimSlope_MC(VEC3 U_c, VEC3 U_l, VEC3 U_r, DType dx_c, DType dx_l=-1, DType dx_r=-1, DType coef_a=-1, DType coef_b=-1);
VEC4 Alg_LimSlope_MC(VEC4 U_c, VEC4 U_l, VEC4 U_r, DType dx_c, DType dx_l=-1, DType dx_r=-1, DType coef_a=-1, DType coef_b=-1);

// van Leer limiter
DType Alg_LimSlope_vanLeer(DType U_c, DType U_l, DType U_r, DType dx_c, DType dx_l=-1, DType dx_r=-1, int k_vl=-1, DType denom=-1, DType nomtr=-1, DType coef_a=-1, DType coef_b=-1);
VEC3 Alg_LimSlope_vanLeer(VEC3 U_c, VEC3 U_l, VEC3 U_r, DType dx_c, DType dx_l=-1, DType dx_r=-1, int k_vl=-1, DType denom=-1, DType nomtr=-1, DType coef_a=-1, DType coef_b=-1);
VEC4 Alg_LimSlope_vanLeer(VEC4 U_c, VEC4 U_l, VEC4 U_r, DType dx_c, DType dx_l=-1, DType dx_r=-1, int k_vl=-1, DType denom=-1, DType nomtr=-1, DType coef_a=-1, DType coef_b=-1);

// calculate the parameter k_vl in generalized van Leer limiter
int Alg_LimSlope_vanLeer_calc_k(DType a_in, DType b_in, DType& denom_ot, DType& nomtr_ot);

// van Albada limiter
DType Alg_LimSlope_vanAlbada(DType U_c, DType U_l, DType U_r, DType dx_c, DType dx_l=-1, DType dx_r=-1, int k_va=-1, DType coef_a=-1, DType coef_b=-1);
VEC3 Alg_LimSlope_vanAlbada(VEC3 U_c, VEC3 U_l, VEC3 U_r, DType dx_c, DType dx_l=-1, DType dx_r=-1, int k_va=-1, DType coef_a=-1, DType coef_b=-1);
VEC4 Alg_LimSlope_vanAlbada(VEC4 U_c, VEC4 U_l, VEC4 U_r, DType dx_c, DType dx_l=-1, DType dx_r=-1, int k_va=-1, DType coef_a=-1, DType coef_b=-1);

// calculate the parameter k_va in generalized van Albada limiter
int Alg_LimSlope_vanAlbada_calc_k(DType a_in, DType b_in);

// Berger's first limiter
DType Alg_LimSlope_Berger1(DType U_c, DType U_l, DType U_r, DType dx_c, DType dx_l=-1, DType dx_r=-1);
VEC3 Alg_LimSlope_Berger1(VEC3 U_c, VEC3 U_l, VEC3 U_r, DType dx_c, DType dx_l=-1, DType dx_r=-1);
VEC4 Alg_LimSlope_Berger1(VEC4 U_c, VEC4 U_l, VEC4 U_r, DType dx_c, DType dx_l=-1, DType dx_r=-1);

// Berger's second limiter
DType Alg_LimSlope_Berger2(DType U_c, DType U_l, DType U_r, DType dx_c, DType dx_l=-1, DType dx_r=-1);
VEC3 Alg_LimSlope_Berger2(VEC3 U_c, VEC3 U_l, VEC3 U_r, DType dx_c, DType dx_l=-1, DType dx_r=-1);
VEC4 Alg_LimSlope_Berger2(VEC4 U_c, VEC4 U_l, VEC4 U_r, DType dx_c, DType dx_l=-1, DType dx_r=-1);



#endif
