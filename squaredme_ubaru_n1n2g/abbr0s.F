#if 0
* abbr0s.F
* abbreviations for {-F[3, {1, Col1}], F[3, {1, Col2}]} -> {F[11, {1}], F[11, {2}], V[5, {Glu5}]}
* generated by FormCalc 9.4 (7 Jun 2016) on 24-Mar-2017 11:09
#endif

#include "vars.h"
#include "inline.h"
#include "contains.h"


	subroutine abbr0s
	implicit none

#include "vars.h"
#include "inline.h"

	MatSUN(1,1) = 4

	Opt1 = 3*MU2 - S + MNeu2(1) + MNeu2(2)

	Sub3 = 3*CW*ZNeu(1,2) + SW*ZNeu(1,1)

	Sub9 = 3*CW*ZNeu(2,2) + SW*ZNeu(2,1)

	Sub20 = ZNeu(2,3)*ZNeuC(1,3) - ZNeu(2,4)*ZNeuC(1,4)

	Sub22 = ZNeu(1,3)*ZNeuC(2,3) - ZNeu(1,4)*ZNeuC(2,4)

	LOOP(Sfe6, 1,2,1)

        Sub1(Sfe6) = 4*MW*(SB*SW*USf(Sfe6,2,3,1)*ZNeu(2,1)) - 
     &    3*CW*MU*(USf(Sfe6,1,3,1)*ZNeu(2,4))

        Sub4(Sfe6) = MW*(SB*Sub3*USfC(Sfe6,1,3,1)) + 
     &    3*CW*MU*(USfC(Sfe6,2,3,1)*ZNeu(1,4))

        Sub5(Sfe6) = 4*MW*
     &     (SB*SW*USfC(Sfe6,2,3,1)*ZNeuC(1,1)) - 
     &    3*CW*MU*(USfC(Sfe6,1,3,1)*ZNeuC(1,4))

        Sub7(Sfe6) = MW*(SB*USf(Sfe6,1,3,1)*ZNeuC(2,2)) + 
     &    MU*(USf(Sfe6,2,3,1)*ZNeuC(2,4))

        Sub8(Sfe6) = 3*CW*Sub7(Sfe6) + 
     &    MW*(SB*SW*USf(Sfe6,1,3,1)*ZNeuC(2,1))

        Sub10(Sfe6) = MW*(SB*Sub9*USfC(Sfe6,1,3,1)) + 
     &    3*CW*MU*(USfC(Sfe6,2,3,1)*ZNeu(2,4))

        Sub12(Sfe6) = 4*MW*(SB*SW*USf(Sfe6,2,3,1)*ZNeu(1,1)) - 
     &    3*CW*MU*(USf(Sfe6,1,3,1)*ZNeu(1,4))

        Sub14(Sfe6) = MW*(SB*USf(Sfe6,1,3,1)*ZNeuC(1,2)) + 
     &    MU*(USf(Sfe6,2,3,1)*ZNeuC(1,4))

        Sub15(Sfe6) = 3*CW*Sub14(Sfe6) + 
     &    MW*(SB*SW*USf(Sfe6,1,3,1)*ZNeuC(1,1))

        Sub16(Sfe6) = 4*MW*
     &     (SB*SW*USfC(Sfe6,2,3,1)*ZNeuC(2,1)) - 
     &    3*CW*MU*(USfC(Sfe6,1,3,1)*ZNeuC(2,4))

	ENDLOOP(Sfe6)

#include "contains.h"

	end

