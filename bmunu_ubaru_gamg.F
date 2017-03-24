      subroutine bmunu_ubaru_gamg(p,bmunu)
      implicit none
#include "PhysPars.h"
      double precision pi
      parameter (pi = 4.D0*datan(1.D0))
      double precision p(0:3,4)
      integer hels, helt
      double complex es4(0:3), ecs4(0:3)
      double complex et4(0:3), ect4(0:3)
      integer alind, beind, i, j
      double precision bmunu(0:3,0:3,4)
      double complex k1(0:3)
      double complex k2(0:3)
      double complex k3(0:3)
      double complex k4(0:3)
      double precision S, T, U
      integer Placeholder
      double complex Pair18,Pair16,Pair14,Pair15,Pair21,Pair20,Pair17
      double complex Pair19,Abb61,Abb63,Abb60,Abb62,Sub2,Sub4,Sub3

      double precision Epsilon, DotP, Den, Kronecker
      double precision momsq, momsum2sq, momsum3sq
      double complex cDotP
      external Epsilon, DotP, Den, Kronecker
      external momsq, momsum2sq, momsum3sq
      external cDotP

      bmunu(:,:,:) = 0D0
      S   = momsum2sq(p(:,1), p(:,2))
      T   = momsum2sq(p(:,1),-p(:,3))
      U   = momsum2sq(p(:,2),-p(:,3))

      do i=0,3
      k1(i) = dcmplx(p(i,1))
      k2(i) = dcmplx(p(i,2))
      k3(i) = dcmplx(p(i,3))
      k4(i) = dcmplx(p(i,4))
      enddo

      do alind=0,3
      do beind=0,3
      do hels=-1,1,2
      do helt=-1,1,2

      call polvector(dreal(k4), hels, es4)
      call polvector(dreal(k4), helt, et4)
      ecs4(:) = dconjg(es4(:))
      ect4(:) = dconjg(et4(:))

      Pair18 = cDotP(k1,ect4)
      Pair16 = cDotP(k1,es4)
      Pair14 = cDotP(k2,ect4)
      Pair15 = cDotP(k2,es4)
      Pair21 = cDotP(k3,ect4)
      Pair20 = cDotP(k3,es4)
      Pair17 = cDotP(k4,ect4)
      Pair19 = cDotP(k4,es4)
      Abb61 = -(Pair14*Pair15) + Pair16*Pair18
      Abb63 = Pair15*Pair17 + 4*Pair16*Pair18 + Pair14*Pair19
      Abb60 = 4*Pair14*Pair15 + Pair16*Pair17 + Pair18*Pair19
      Abb62 = 2*Pair16*Pair18 - Pair17*Pair20 - Pair19*(Pair17 + Pair21)
      Sub2 = 2*Abb60*T*es4(alind) + T**2*et4(alind)
      Sub4 = 2*Abb63*U*es4(alind) + U**2*et4(alind)
      Sub3 = -4*(Abb62*S + 2*Abb61*U)*es4(alind) + S**2*et4(alind)

      bmunu(alind,beind,4) = bmunu(alind,beind,4) +
     -   dreal(-(Alfa*Alfas*Pi**2*((-128*Sub3*Den(T,MU2)*Den(U,MU2))/9. +
     -      (128*(Sub2*Den(T,MU2)**2 + Sub4*Den(U,MU2)**2))/9.)*
     -    ect4(beind)))

      enddo
      enddo
      enddo
      enddo

      end
