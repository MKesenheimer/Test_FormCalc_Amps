#if 0
* abbr0h.F
* abbreviations for {-F[3, {1, Col1}], F[3, {1, Col2}]} -> {V[1], V[5, {Glu4}]}
* generated by FormCalc 9.4 (7 Jun 2016) on 10-Mar-2017 15:46
#endif

#include "vars.h"
#include "inline.h"
#include "contains.h"


	subroutine abbr0h
	implicit none

#include "vars.h"
#include "inline.h"

	F1 = ChainV1(Spinor(1,-1,2),0,ec(3),0,Spinor(2,1,1))

        F2 = ChainV3(Spinor(1,-1,2),0,ec(3),ec(4),k(3),0,
     &    Spinor(2,1,1))

	F3 = ChainV1(Spinor(1,-1,2),0,ec(4),0,Spinor(2,1,1))

	F4 = -ChainB1(Spinor(1,-1,2),0,ec(3),0,Spinor(2,1,1))

        F5 = -ChainB3(Spinor(1,-1,2),0,ec(3),ec(4),k(3),0,
     &     Spinor(2,1,1))

	F6 = -ChainB1(Spinor(1,-1,2),0,ec(4),0,Spinor(2,1,1))

	F7 = ChainV1(Spinor(1,-1,2),0,k(3),0,Spinor(2,1,1))

	F8 = -ChainB1(Spinor(1,-1,2),0,k(3),0,Spinor(2,1,1))

	Pair1 = Pair(ec(4),k(3))

	Pair2 = Pair(ec(3),k(2))

	Pair3 = Pair(ec(3),k(1))

	Pair4 = Pair(ec(3),ec(4))

	Sub1 = F2 + F5 - 2*(F1 + F4)*Pair1 + 2*(F3 + F6)*Pair3

	Sub2 = (F3 + F6)*Pair2 + (F7 + F8)*Pair4

#include "contains.h"

	end

