	subroutine full_ubaru_gamg(p,born)
	implicit none
#include "PhysPars.h"
	double precision pi
	parameter (pi = 4.D0*datan(1.D0))
	double precision p(0:3,4)
	integer i, j
	double precision born
	double precision k1(0:3)
	double precision k2(0:3)
	double precision k3(0:3)
	double precision k4(0:3)
	double precision S, T, U
	integer Placeholder

	double precision Epsilonk, DotP, Den, Kronecker
	double precision momsq, momsum2sq, momsum3sq
	double complex cDotP
	external Epsilonk, DotP, Den, Kronecker
	external momsq, momsum2sq, momsum3sq
	external cDotP

	born = 0D0
	S   = momsum2sq(p(:,1), p(:,2))
	T   = momsum2sq(p(:,1),-p(:,3))
	U   = momsum2sq(p(:,2),-p(:,3))

	do i=0,3
	if(i.eq.0) then
	k1(i) = p(i,1)
	k2(i) = p(i,2)
	k3(i) = p(i,3)
	k4(i) = p(i,4)
	else
	k1(i) = -p(i,1)
	k2(i) = -p(i,2)
	k3(i) = -p(i,3)
	k4(i) = -p(i,4)
	endif
	enddo

	born = born + (        -(Alfa*Alfas*Pi**2*((-512*S**2*Den(T,MU2)*Den(U,MU2))/9. + 
     -      (512*(T**2*Den(T,MU2)**2 + U**2*Den(U,MU2)**2))/9.)))

	end
