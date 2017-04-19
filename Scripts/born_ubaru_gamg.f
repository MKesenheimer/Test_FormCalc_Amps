	subroutine born_ubaru_gamg(p,born)
	implicit none
#include "PhysPars.h"
	double precision pi
	parameter (pi = 4.D0*datan(1.D0))
	double precision p(0:3,4)
	integer hel
	double complex e4(0:3), ec4(0:3)
	integer i, j
	double precision born
	double complex k1(0:3)
	double complex k2(0:3)
	double complex k3(0:3)
	double complex k4(0:3)
	double precision S, T, U
	integer Placeholder
	double complex Pair10,Pair13,Pair9,Pair5,Pair6,Pair1,Abb57
	double complex Abb58,Abb59,Abb55,Abb53,Abb54,Abb56,Sub1

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
	k1(i) = dcmplx(p(i,1))
	k2(i) = dcmplx(p(i,2))
	k3(i) = dcmplx(p(i,3))
	k4(i) = dcmplx(p(i,4))
	enddo
	do hel=-1,1,2

	call polvector(dreal(k4), hel, e4)
	ec4(:) = dconjg(e4(:))

      Pair10 = cDotP(e4,k1)
      Pair13 = cDotP(e4,k2)
      Pair9 = cDotP(e4,k3)
      Pair5 = cDotP(ec4,k1)
      Pair6 = cDotP(ec4,k2)
      Pair1 = cDotP(ec4,k3)
      Abb57 = (Pair10 - Pair13)*(Pair5 - Pair6)
      Abb58 = 2*Pair10*Pair5 - Pair1*Pair9
      Abb59 = 2*Pair13*Pair6 - Pair1*Pair9
      Abb55 = Pair10*(4*Pair5 + Pair6) + Pair13*(Pair5 + 2*Pair6)
      Abb53 = Pair10*(2*Pair5 + Pair6) + Pair13*(Pair5 + 4*Pair6)
      Abb54 = Pair1*Pair10 + Pair5*Pair9
      Abb56 = Pair1*Pair13 + Pair6*Pair9
      Sub1 = S**2 + 2*Abb59*(S + U) - 2*(Abb57*S + Abb58*U)

	born = born + dreal(        -(Alfa*Alfas*Pi**2*((-256*Sub1*Den(T,MU2)*Den(U,MU2))/9. + 
     -      (256*(T*(Abb53 - Abb54 + T)*Den(T,MU2)**2 + 
     -           U*(Abb55 - Abb56 + U)*Den(U,MU2)**2))/9.)))
	enddo

	end
