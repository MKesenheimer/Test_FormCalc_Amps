#if 0
* vars.h
* variable declarations
* generated by FormCalc 9.4 (7 Jun 2016) on 10-Mar-2017 15:46
#endif

#ifndef vars_h
#define vars_h

#define SQUAREDME
#define LEGS 4

#include "decl.h"

#else

#include "decl.h"

	RealType S, T, U
	common /varXa/ S, T, U

	HelType F1, F2, F3, F4, F5, F6, F7, F8, Pair1, Pair2, Pair3
	HelType Pair4, Sub1, Sub2
	common /varXh/ F1, F2, F3, F4, F5, F6, F7, F8, Pair1, Pair2
	common /varXh/ Pair3, Pair4, Sub1, Sub2

	HelType Ctree(HelDim(1))
	ComplexType MatSUN(1,1)
	common /formfactors/ Ctree, MatSUN

#endif
