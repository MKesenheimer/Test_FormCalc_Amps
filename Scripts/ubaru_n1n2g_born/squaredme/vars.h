#if 0
* vars.h
* variable declarations
* generated by FormCalc 9.4 (7 Jun 2016) on 24-Mar-2017 11:09
#endif

#ifndef vars_h
#define vars_h

#define SQUAREDME
#define LEGS 5

#include "decl.h"

#else

#include "decl.h"

	ComplexType Opt1, Sub1(2), Sub3, Sub4(2), Sub5(2), Sub7(2)
	ComplexType Sub8(2), Sub9, Sub10(2), Sub12(2), Sub14(2)
	ComplexType Sub15(2), Sub16(2), Sub20, Sub22
	common /varXs/ Opt1, Sub1, Sub3, Sub4, Sub5, Sub7, Sub8, Sub9
	common /varXs/ Sub10, Sub12, Sub14, Sub15, Sub16, Sub20
	common /varXs/ Sub22

	ComplexType Sub31, Sub32
	RealType S, T, T14, U, T24, S34
	common /varXa/ Sub31, Sub32, S, T, T14, U, T24, S34

	HelType F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12
	HelType F13, F14, F15, F16, F17, F18, F19, F20, F21, F22
	HelType F23, F24, F25, F26, F27, F28, F29, F30, F31, F32
	HelType F33, F34, F35, F36, F37, F38, F39, F40, F41, F42
	HelType F43, F44, F45, F46, F47, F48, F49, F50, F51, F52
	HelType F53, F54, Pair1, Pair2, Pair3, Pair4
	HelType Sub17(HelDim(2)), Sub18(HelDim(2)), Sub19, Sub21
	HelType Sub23, Sub24, Sub25, Sub26, Sub27, Sub28, Sub29
	HelType Sub30, Sub34(HelDim(2)), Abb1, Abb2, Sub2, Sub6
	HelType Sub11, Sub13, Sub33(HelDim(2)), Sub35(HelDim(2))
	common /varXh/ F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11
	common /varXh/ F12, F13, F14, F15, F16, F17, F18, F19, F20
	common /varXh/ F21, F22, F23, F24, F25, F26, F27, F28, F29
	common /varXh/ F30, F31, F32, F33, F34, F35, F36, F37, F38
	common /varXh/ F39, F40, F41, F42, F43, F44, F45, F46, F47
	common /varXh/ F48, F49, F50, F51, F52, F53, F54, Pair1
	common /varXh/ Pair2, Pair3, Pair4, Sub17, Sub18, Sub19
	common /varXh/ Sub21, Sub23, Sub24, Sub25, Sub26, Sub27
	common /varXh/ Sub28, Sub29, Sub30, Sub34, Abb1, Abb2, Sub2
	common /varXh/ Sub6, Sub11, Sub13, Sub33, Sub35

	integer Sfe6
	common /indices/ Sfe6

	HelType Ctree(HelDim(1))
	ComplexType MatSUN(1,1)
	common /formfactors/ Ctree, MatSUN

#endif
