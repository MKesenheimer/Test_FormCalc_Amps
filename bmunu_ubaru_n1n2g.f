      subroutine bmunu_ubaru_n1n2g(p,bmunu)
      implicit none
#include "PhysPars.h"
      double precision pi
      parameter (pi = 4.D0*datan(1.D0))
      double precision p(0:3,5)
      integer hels, helt
      double complex es5(0:3), ecs5(0:3)
      double complex et5(0:3), ect5(0:3)
      integer alind, beind, i, j
      double precision bmunu(0:3,0:3,5)
      double complex k1(0:3)
      double complex k2(0:3)
      double complex k3(0:3)
      double complex k4(0:3)
      double complex k5(0:3)
      double precision S, T, U, S34, T14, T24
      integer Placeholder
      double complex Sub102,Sub2,Sub1,Sub5,Pair15,Pair16,Eps29,Eps32
      double complex Eps30,Sub174,Sub22,Sub181,Sub110,Pair11,Pair9
      double complex Pair10,Pair12,Pair13,Pair14,Eps37,Eps31,Eps38
      double complex Abb467,Abb460,Abb457,Abb456,Abb455,Abb491,Abb472
      double complex Abb499,Abb470,Abb480,Abb502,Abb505,Sub99,Sub97
      double complex Sub131,Sub143,Sub135,Sub100,Sub130,Sub11,Sub167
      double complex Abb459,Abb452,Abb453,Abb469,Abb468,Sub177,Sub137
      double complex Sub179,Sub113,Sub175,Sub104,Sub152,Sub106,Sub182
      double complex Sub105,Sub176,Sub136,Sub170,Sub107,Sub180,Sub108
      double complex Sub119,Sub109,Sub111,Eps33,Eps36,Eps34,Eps35
      double complex Eps40,Eps39,Abb545,Abb540,Abb541,Abb465,Abb466
      double complex Abb500,Abb486,Abb498,Abb487,Abb504,Abb463,Abb461
      double complex Abb482,Abb510,Abb503,Abb492,Abb462,Abb549,Abb497
      double complex Abb548,Abb546,Abb511,Abb550,Abb551,Abb528,Abb516
      double complex Abb543,Abb481,Abb488,Abb493,Abb494,Abb464,Sub133
      double complex Sub120,Sub92,Sub93,Sub9,Sub132,Sub98,Sub96
      double complex Abb506,Abb483,Abb547,Abb542,Abb454,Sub153,Sub148
      double complex Sub183,Sub165,Sub169,Sub134,Sub101,Sub147,Sub114
      double complex Sub168,Abb507,Abb458,Sub172,Sub112,Sub171,Sub138
      double complex Sub178,Sub103,Sub173,Sub139,Eps6,Eps45,Eps43
      double complex Abb485,Abb534,Abb519,Abb512,Abb526,Abb523,Abb517
      double complex Abb532,Abb552,Abb538,Abb518,Abb471,Abb531,Abb474
      double complex Abb520,Sub163,Sub67,Sub6,Abb490,Abb501,Abb539
      double complex Abb473,Abb513,Sub155,Sub122,Sub164,Sub118,Sub117
      double complex Sub187,Sub184,Sub149,Sub121,Sub150,Sub156,Eps42
      double complex Eps41,Eps50,Eps49,Abb536,Abb484,Abb509,Abb535
      double complex Abb508,Abb533,Abb529,Abb530,Abb514,Abb475,Abb524
      double complex Abb476,Abb489,Abb537,Abb544,Abb495,Sub94,Sub161
      double complex Sub128,Sub166,Sub160,Abb515,Abb525,Abb496,Sub144
      double complex Sub140,Sub189,Sub95,Sub129,Sub162,Sub123,Sub188
      double complex Sub141,Sub190,Sub145,Sub115,Sub154,Sub146,Sub142
      double complex Sub116,Sub151,Sub185,Sub186,Sub157,Sub158,Sub159
      double complex Eps54,Eps53,Eps52,Eps51,Eps48,Eps47,Eps46,Eps44
      double complex Opt150,Opt151,Abb521,Abb522,Abb527,Sub191,Sub192
      double complex Sub193,Sub194,Abb477,Abb478,Abb479,Sub124,Sub125
      double complex Sub126,Sub127

      double precision Epsilonk, DotP, Den, Kronecker
      double precision momsq, momsum2sq, momsum3sq
      double complex cDotP
      external Epsilonk, DotP, Den, Kronecker
      external momsq, momsum2sq, momsum3sq
      external cDotP

      bmunu(:,:,:) = 0D0
      S   = momsum2sq(p(:,1), p(:,2))
      T   = momsum2sq(p(:,1),-p(:,3))
      U   = momsum2sq(p(:,2),-p(:,3))
      S34 = momsum2sq(p(:,3), p(:,4))
      T14 = momsum2sq(p(:,1),-p(:,4))
      T24 = momsum2sq(p(:,2),-p(:,4))

      do i=0,3
      k1(i) = dcmplx(p(i,1))
      k2(i) = dcmplx(p(i,2))
      k3(i) = dcmplx(p(i,3))
      k4(i) = dcmplx(p(i,4))
      k5(i) = dcmplx(p(i,5))
      enddo

      do alind=0,3
      do beind=0,3
      do hels=-1,1,2
      do helt=-1,1,2

      call polvector(dreal(k5), hels, es5)
      call polvector(dreal(k5), helt, et5)
      ecs5(:) = dconjg(es5(:))
      ect5(:) = dconjg(et5(:))

      Sub102 = S*S34 - (T + T14)*(T24 + U)
      Sub2 = ZNeu(2,3)*ZNeuC(1,3) - ZNeu(2,4)*ZNeuC(1,4)
      Sub1 = ZNeu(1,3)*ZNeuC(2,3) - ZNeu(1,4)*ZNeuC(2,4)
      Sub5 = -1 + 8*CW2 - 16*(CW2**2 + SW2**2)
      Pair15 = cDotP(k4,ect5)
      Pair16 = cDotP(k4,es5)
      Eps29 = epsilonk(k1,k3,es5,ect5)
      Eps32 = epsilonk(k2,k3,es5,ect5)
      Eps30 = epsilonk(k2,k4,es5,ect5)
      Sub174 = (T + U)*(2*S34 + T + U) - 4*S34*MNeu2(1)
      Sub22 = Sub1*Sub2*(MNeu2(1) - MNeu2(2))**2*(MNeu2(1) + MNeu2(2)) + 
     -  (Sub1**2 + Sub2**2)*MNeu(1)*MNeu(2)*(MNeu2(1) + MNeu2(2))**2
      Sub181 = Sub5*(T + U)*et5(alind)*(3*MNeu2(1) + MNeu2(2)) + 
     -  6*Eps29*(1 - 4*CW2 + 4*SW2)*es5(alind)*(5*MNeu2(1) + MNeu2(2))
      Sub110 = Sub5*et5(alind)*(MNeu2(1) - MNeu2(2))**2 + 
     -  12*(1 - 4*CW2 + 4*SW2)*es5(alind)*
     -   (Eps30*MNeu2(1) + Eps32*MNeu2(2))
      Pair11 = cDotP(k1,ect5)
      Pair9 = cDotP(k1,es5)
      Pair10 = cDotP(k2,ect5)
      Pair12 = cDotP(k2,es5)
      Pair13 = cDotP(k3,ect5)
      Pair14 = cDotP(k3,es5)
      Eps37 = epsilonk(k1,k2,es5,ect5)
      Eps31 = epsilonk(k1,k4,es5,ect5)
      Eps38 = epsilonk(k1,k5,es5,ect5)
      Abb467 = (-Pair10 + Pair11)*(-Pair12 + Pair9)
      Abb460 = Eps29 - Eps32
      Abb457 = -Eps30 + Eps37
      Abb456 = -Eps32 + Eps37
      Abb455 = Eps37 - Eps38
      Abb491 = Pair14*Pair15 + Pair13*Pair16
      Abb472 = Pair11*Pair16 + Pair15*Pair9
      Abb499 = -2*Pair15*Pair16*S + Abb472*T
      Abb470 = -(Eps32*S34) + Eps29*(-T + T14)
      Abb480 = Eps30*S34 + Eps31*(-T + T14)
      Abb502 = (Eps30 + Eps31)*S34 + Eps31*(-T + T14 + 3*T24)
      Abb505 = (Eps29 + Eps32)*S34 + Eps29*(T - T14 + 3*U)
      Sub99 = T24*(S - 2*T14 + T24) - S34*(3*T14 + T24)
      Sub97 = U*(S - 2*T + U) - S34*(3*T + U)
      Sub131 = -(T*(S + T - 2*U)) + S34*(T + 3*U)
      Sub143 = 4*Pair15*Pair16*MNeu2(1)**2 + Abb472*(S - S34)*(2*MNeu2(1) - MNeu2(2))
      Sub135 = 2*S*S34 - T*(T14 + T24) - (T14 - 3*T24)*U
      Sub100 = 2*S*S34 + T*(3*T14 - T24) - (T14 + T24)*U
      Sub130 = S34*(T14*T24 + T*U) - (T24 - U)*(T*T24 - T14*U) + 
     -  S*(2*T24*U + S34*(T24 + U))
      Sub11 = MNeu2(1)**2 + MNeu2(1)*MNeu2(2) + MNeu2(2)**2
      Sub167 = 2*S**2 + T*(T14 - 3*T24) - (3*T14 - T24)*U + 
     -  S*(6*S34 + T + T14 + T24 + U)
      Abb459 = -Eps30 + Eps31
      Abb452 = Pair11*Pair12 + Pair10*Pair9
      Abb453 = (2*Pair10 + Pair11)*Pair12 + Pair10*Pair9
      Abb469 = Pair10*Pair9 + Pair11*(Pair12 + 2*Pair9)
      Abb468 = Eps37*S34 + Abb455*(T24 + U)
      Sub177 = 4*Abb467*(S + 2*S34)*es5(alind) - Sub167*et5(alind)
      Sub137 = Abb469*Sub5 - 3*Abb455*(1 - 4*CW2 + 4*SW2)
      Sub179 = 2*Abb467*Sub5 + 3*Eps37*(1 - 4*CW2 + 4*SW2)
      Sub113 = -24*Abb470*(1 - 4*CW2 + 4*SW2)*es5(alind) + Sub5*Sub99*et5(alind)
      Sub175 = -3*Abb460*(Sub1**2 + Sub2**2)*(1 - 4*CW2 + 4*SW2)*MNeu(1)*
     -   MNeu(2) + 4*Abb491*Sub1*Sub2*Sub5*MNeu2(1)
      Sub104 = -3*(Eps31 + Eps32)*(Sub1**2 + Sub2**2)*(1 - 4*CW2 + 4*SW2)*
     -   MNeu(1)*MNeu(2) + 5*Pair10*Pair12*Sub1*Sub2*Sub5*MNeu2(1)
      Sub152 = Sub1*Sub130*Sub2 + (Sub1**2 + Sub2**2)*(T + T14 + T24 + U)*
     -   MNeu(1)*MNeu(2)*(MNeu2(1) + MNeu2(2))
      Sub106 = Abb472*S*(2*MNeu2(1) + MNeu2(2)) + Abb499*(3*MNeu2(1) + MNeu2(2))
      Sub182 = Sub5*(T14 + T24)*et5(alind)*(MNeu2(1) + 3*MNeu2(2)) + 
     -  6*Eps31*(1 - 4*CW2 + 4*SW2)*es5(alind)*(MNeu2(1) + 5*MNeu2(2))
      Sub105 = -3*(Eps29 + Eps30)*(Sub1**2 + Sub2**2)*(1 - 4*CW2 + 4*SW2)*
     -   MNeu(1)*MNeu(2)*MNeu2(1) + 
     -  Pair15*Pair16*Sub1*Sub2*Sub5*MNeu2(1)**2 + Sub104*MNeu2(2)
      Sub176 = -4*Sub175*es5(alind) + Sub1*Sub174*Sub2*Sub5*et5(alind)
      Sub136 = -3*Abb468*(1 - 4*CW2 + 4*SW2) + Sub5*(Abb452*S34 + Abb469*(T24 + U))
      Sub170 = Sub1*Sub2*Sub5*(T14 + T24)*(2*S34 + T14 + T24)*et5(alind) + 
     -  12*Abb459*(Sub1**2 + Sub2**2)*(1 - 4*CW2 + 4*SW2)*es5(alind)*
     -   MNeu(1)*MNeu(2)
      Sub107 = Sub1*Sub106*Sub2*es5(alind) - 
     -  (Sub1**2 + Sub2**2)*et5(alind)*MNeu(1)*MNeu(2)*
     -   (MNeu2(1) + MNeu2(2))**2
      Sub180 = Sub179*(Sub1**2 + Sub2**2)*MNeu(1)*MNeu(2) + 
     -  3*Sub1*Sub2*(1 - 4*CW2 + 4*SW2)*
     -   (Eps30*MNeu2(1) + Eps32*MNeu2(2))
      Sub108 = Abb453*Sub5 - 9*Eps37*(1 - 4*CW2 + 4*SW2)
      Sub119 = Abb452*Sub5 + 2*Eps37*(1 - 4*CW2 + 4*SW2)
      Sub109 = 4*Sub108*es5(alind) - Sub5*(T + T14 + T24 + U)*et5(alind)
      Sub111 = Sub1*Sub110*Sub2 + Sub109*(Sub1**2 + Sub2**2)*MNeu(1)*MNeu(2)
      Eps33 = epsilonk(k1,k2,k3,ect5)
      Eps36 = epsilonk(k1,k2,k3,es5)
      Eps34 = epsilonk(k1,k2,k4,ect5)
      Eps35 = epsilonk(k1,k2,k4,es5)
      Eps40 = epsilonk(k1,k2,k5,ect5)
      Eps39 = epsilonk(k1,k2,k5,es5)
      Abb545 = Eps37*(S34 + T24 + U)
      Abb540 = -Eps29 + Eps37
      Abb541 = -Eps31 + Eps37
      Abb465 = -Eps29 + 2*Eps37
      Abb466 = -Eps31 + 2*Eps37
      Abb500 = Eps35*Pair13 - Eps34*Pair14
      Abb486 = Pair12*Pair13 + Pair10*Pair14
      Abb498 = Eps36*Pair15 - Eps33*Pair16
      Abb487 = Pair12*Pair15 + Pair10*Pair16
      Abb504 = 3*Pair12*Pair15 + (3*Pair10 + 4*Pair15)*Pair16
      Abb463 = (-2*Pair10 + Pair11)*Pair12 + Pair10*Pair9
      Abb461 = Pair10*Pair12 + Pair11*Pair9
      Abb482 = Pair11*Pair14 + Pair13*Pair9
      Abb510 = Pair11*(Pair12 + Pair14) + (Pair10 + Pair13)*Pair9
      Abb503 = (3*Pair11 + 4*Pair15)*Pair16 + 3*Pair15*Pair9
      Abb492 = Pair11*(2*Pair12 + Pair9) + Pair10*(-Pair12 + 2*Pair9)
      Abb462 = -(Pair10*Pair9) + Pair11*(-Pair12 + 2*Pair9)
      Abb549 = Pair11*(-2*Pair12 + Pair9) - Pair10*(Pair12 + 2*Pair9)
      Abb497 = -2*Pair13*Pair14*S + Abb482*T14
      Abb548 = -2*Pair13*Pair14*S + Abb486*T24
      Abb546 = Abb469 + Pair15*Pair16
      Abb511 = Abb452 + Abb486
      Abb550 = -2*Pair10*Pair12 - (4*Pair10 + 3*Pair13)*Pair9 + 
     -  Pair11*(-4*Pair12 - 3*Pair14 + 2*Pair9)
      Abb551 = -2*Pair10*Pair12 - (4*Pair10 + 3*Pair15)*Pair9 + 
     -  Pair11*(-4*Pair12 - 3*Pair16 + 2*Pair9)
      Abb528 = Pair10*Pair12*T + Abb491*(T14 + T24) + Pair11*Pair9*U
      Abb516 = Pair10*Pair12*T14 + Pair11*Pair9*T24 + Abb491*(T + U)
      Abb543 = -Abb491 + Pair12*(Pair13 + Pair15) + Pair10*(Pair14 + Pair16)
      Abb481 = Abb452 + 2*Pair10*Pair12 + Pair15*Pair16
      Abb488 = -Abb491 + Pair11*(Pair14 + Pair16) + (Pair13 + Pair15)*Pair9
      Abb493 = -2*Pair10*Pair12 + 3*Pair12*Pair13 + 3*Pair10*Pair14 + 
     -  4*Pair10*Pair9 + 2*Pair11*(2*Pair12 + Pair9)
      Abb494 = -2*Pair10*Pair12 + 3*Pair12*Pair15 + 3*Pair10*Pair16 + 
     -  4*Pair10*Pair9 + 2*Pair11*(2*Pair12 + Pair9)
      Abb464 = -(Eps39*Pair10) + Eps40*Pair12 + 
     -  Eps36*(Pair10 + Pair11 - 2*Pair13) - Eps39*Pair13 + 
     -  Eps40*Pair14 - Eps33*(Pair12 - 2*Pair14 + Pair9)
      Sub133 = -(T14*(S + T14 - 2*T24)) + S34*(T14 + 3*T24)
      Sub120 = Abb486*(MNeu2(1) - 2*MNeu2(2)) - Abb487*(2*MNeu2(1) - MNeu2(2))
      Sub92 = MNeu2(1)**2 - 5*MNeu2(2)**2
      Sub93 = 5*MNeu2(1)**2 - MNeu2(2)**2
      Sub9 = 2*MNeu2(1)*MNeu2(2) + 3*(MNeu2(1)**2 + MNeu2(2)**2)
      Sub132 = 2*Abb472*T - 2*Abb482*T14 - Abb487*T24 + Abb486*U
      Sub98 = Abb482*T - Abb472*T14 - 2*Abb486*T24 + 2*Abb487*U
      Sub96 = S*(2*T*T14 + S34*(T + T14)) - T**2*T24 + 
     -  T14*(S34*T24 - T14*U) + T*(S34*U + T14*(T24 + U))
      Abb506 = Eps39*(3*Pair10 + Pair11) - Eps40*(3*Pair12 + Pair9)
      Abb483 = -(Abb482*S34) + 2*Pair15*Pair16*(3*T - T14 + T24) + 2*Abb481*U
      Abb547 = -(Abb486*S34) + 2*Abb546*T + 2*Pair15*Pair16*(T14 - T24 + 3*U)
      Abb542 = Abb541*S34*T24 + Eps29*T24**2 + Abb540*S34*U + 
     -  (Eps37 + Eps38)*T24*U + Eps31*U**2
      Abb454 = -(Eps34*Pair12) + Eps35*(Pair10 + 2*Pair11 - Pair13) + 
     -  Eps36*(Pair11 + Pair13) - Eps33*Pair14 + Eps34*Pair14 - 
     -  Eps33*Pair9 - 2*Eps34*Pair9
      Sub153 = 2*Abb543*es5(alind) - (S + S34)*et5(alind)
      Sub148 = Abb482*(S - S34)*(MNeu2(1) - 2*MNeu2(2)) + 
     -  (Abb472*T24 + Abb491*U)*(3*MNeu2(1) - MNeu2(2))
      Sub183 = -(Abb461*Sub9) - Pair13*Pair14*Sub92 + Pair15*Pair16*Sub93
      Sub165 = Abb504*T - 2*(Abb482*T14 + Abb486*T24) + Abb503*U
      Sub169 = 2*Abb461*S - Abb463*(T + T14) + Abb462*(T24 + U)
      Sub134 = -2*Abb549*S34 + Abb491*(4*S + T + T14) - Abb550*T24 - Abb551*U
      Sub101 = 2*Abb492*S34 + Abb494*T + Abb493*T14 + Abb491*(4*S + T24 + U)
      Sub147 = -4*U*MNeu2(1)**2 + U**2*(MNeu2(1) - 2*MNeu2(2)) - 
     -  T24**2*(2*MNeu2(1) - MNeu2(2)) + Sub133*MNeu2(2)
      Sub114 = 2*Sub96 + 4*T*MNeu2(1)**2 - T**2*(MNeu2(1) - 2*MNeu2(2)) + 
     -  T14**2*(2*MNeu2(1) - MNeu2(2)) - 4*(T24 + U)*MNeu2(1)*MNeu2(2)
      Sub168 = 2*Abb464 - (Eps37 + Eps38)*S - (Eps30 + Eps37)*T - 
     -  (Eps32 + Eps37)*T14 - Abb465*T24 - Abb466*U
      Abb507 = Abb506 + Eps37*S + Eps38*S + Eps30*T + Eps32*T14 + 
     -  Eps37*(S34 + T24 + U)
      Abb458 = 2*Abb454 - Abb455*S - Eps37*S34 - Abb457*T - Abb456*T14 + 
     -  Eps31*T24 + Eps29*U
      Sub172 = Sub165*Sub5 - 6*Abb502*(1 - 4*CW2 + 4*SW2)
      Sub112 = 24*Abb480*(1 - 4*CW2 + 4*SW2)*es5(alind) + 
     -  Sub5*(-2*Abb483*es5(alind) + Sub97*et5(alind))
      Sub171 = Sub169*Sub5 - 3*Sub168*(1 - 4*CW2 + 4*SW2)
      Sub138 = Abb542 + Eps37*(MNeu2(1)**2 + MNeu2(2)**2)
      Sub178 = Sub177*Sub5 + 12*Abb507*(1 - 4*CW2 + 4*SW2)*es5(alind)
      Sub103 = 3*Abb458*(1 - 4*CW2 + 4*SW2) + Sub5*(Abb452*S34 + Abb453*(T + T14))
      Sub173 = -2*Sub171*(Sub1**2 + Sub2**2)*MNeu(1)*MNeu(2) + 
     -  Sub1*Sub172*Sub2*MNeu2(1)
      Sub139 = 3*Sub1*Sub138*Sub2*(1 - 4*CW2 + 4*SW2) + 
     -  Sub137*(Sub1**2 + Sub2**2)*MNeu(1)*MNeu(2)*(MNeu2(1) + MNeu2(2))
      Eps6 = epsilonk(k1,k2,k3,k4)
      Eps45 = epsilonk(k2,k3,k4,ect5)
      Eps43 = epsilonk(k2,k3,k4,es5)
      Abb485 = -(Pair12*Pair15) + Pair10*(2*Pair12 - Pair16)
      Abb534 = -(Pair14*Pair15) + Pair13*(2*Pair14 - Pair16)
      Abb519 = Pair14*Pair15 - Pair13*Pair16
      Abb512 = Pair13*Pair14 - Pair15*Pair16
      Abb526 = Eps43*(2*Pair13 + Pair15) - Eps45*(2*Pair14 + Pair16)
      Abb523 = Eps43*(Pair13 + 2*Pair15) - Eps45*(Pair14 + 2*Pair16)
      Abb517 = Pair11*(Pair12 + Pair16) + (Pair10 + Pair15)*Pair9
      Abb532 = -(Pair13*Pair9) + Pair11*(-Pair14 + 2*Pair9)
      Abb552 = -2*Pair15*Pair16*S + Abb487*U
      Abb538 = Abb469 + Pair13*Pair14
      Abb518 = Abb452 + Abb487
      Abb471 = Abb452 + 2*Pair10*Pair12 + Pair13*Pair14
      Abb531 = -Abb487 + 4*Pair15*Pair16 - (Pair10 + 3*Pair15)*Pair9 + 
     -  Pair11*(-Pair12 - 3*Pair16 + 2*Pair9)
      Abb474 = -(Eps34*Pair12) + Eps35*(Pair10 + 3*Pair11 - Pair13) + 
     -  Eps36*(2*Pair11 + Pair13) - Eps33*Pair14 + Eps34*Pair14 - 
     -  2*Eps33*Pair9 - 3*Eps34*Pair9
      Abb520 = -(Eps33*Pair12) + 3*Eps40*Pair12 + 
     -  Eps36*(Pair10 + Pair11 - 2*Pair13) - 
     -  Eps39*(3*Pair10 + Pair11 + Pair13) + 2*Eps33*Pair14 + 
     -  Eps40*Pair14 - Eps33*Pair9 + Eps40*Pair9
      Sub163 = S**2 - T*T24 - T14*U
      Sub67 = 2*S34**2 + S34*(T + T14 + T24 + U) - 2*(T*T24 + T14*U)
      Sub6 = MNeu2(1)**2 + 4*MNeu2(1)*MNeu2(2) + MNeu2(2)**2
      Abb490 = Eps43*Pair13 - Eps45*Pair14
      Abb501 = Eps43*Pair15 - Eps45*Pair16
      Abb539 = Abb487*S34 - 2*Abb538*T14 - 2*Pair13*Pair14*(T + 3*T24 - U)
      Abb473 = -(Abb472*S34) + 2*Abb471*T24 + 2*Pair13*Pair14*(-T + 3*T14 + U)
      Abb513 = Abb512*(2*S + S34) + Pair11*Pair9*(T - T14) + Pair10*Pair12*(-T24 + U)
      Sub155 = Abb487*S*(2*MNeu2(1) + MNeu2(2)) + 
     -  Abb552*(3*MNeu2(1) + MNeu2(2)) + 
     -  Abb486*S*(MNeu2(1) + 2*MNeu2(2))
      Sub122 = Eps38*Sub6 + Eps31*T24*(3*MNeu2(1) + MNeu2(2)) + 
     -  Eps29*U*(MNeu2(1) + 3*MNeu2(2))
      Sub164 = 2*S34*Sub163 + S*Sub67 + (-T**2 + T*T14)*T24 + 
     -  (-T14**2 + (T + T14)*T24)*U - T14*U**2 + T*(-T24**2 + T14*U)
      Sub118 = 6*Abb501*(1 - 4*CW2 + 4*SW2) + Sub5*(Abb491*T + Abb487*T14)
      Sub117 = -6*Abb490*(1 - 4*CW2 + 4*SW2) + Sub5*(Abb486*T + Abb491*T14)
      Sub187 = -6*Abb490*(1 - 4*CW2 + 4*SW2) + Sub5*(Abb511*T + Abb510*U)
      Sub184 = Sub164 + 8*MNeu2(1)*MNeu2(2)*(MNeu2(1) + MNeu2(2))
      Sub149 = (Abb491*T24 + Abb482*U)*(MNeu2(1) - 3*MNeu2(2)) + 
     -  Sub132*(MNeu2(1) - MNeu2(2)) + Abb539*MNeu2(2)
      Sub121 = (S - S34)*Sub120 + Abb473*MNeu2(2) - 
     -  Abb482*S*(MNeu2(1) + 2*MNeu2(2)) + 
     -  2*Abb488*(MNeu2(1)**2 + MNeu2(2)**2)
      Sub150 = -2*(Sub148 - Sub149)*es5(alind) - Sub147*et5(alind) + 
     -  (-2*Abb547*es5(alind) - Sub131*et5(alind))*MNeu2(1) + 
     -  4*(2*Pair13*Pair14*es5(alind) + T24*et5(alind))*MNeu2(2)**2 + 
     -  (-2*Sub134*es5(alind) - Sub135*et5(alind))*(MNeu2(1) + MNeu2(2))
      Sub156 = Sub155 + 6*Abb452*(MNeu2(1) + MNeu2(2))**2
      Eps42 = epsilonk(k1,k3,k4,ect5)
      Eps41 = epsilonk(k1,k3,k4,es5)
      Eps50 = epsilonk(k3,k4,k5,ect5)
      Eps49 = epsilonk(k3,k4,k5,es5)
      Abb536 = Eps29 - Eps31
      Abb484 = -(Pair12*Pair13) + Pair10*(2*Pair12 - Pair14)
      Abb509 = 3*Pair12*Pair13 + (3*Pair10 + 4*Pair13)*Pair14
      Abb535 = Pair14*Pair15 + (Pair13 - 2*Pair15)*Pair16
      Abb508 = (3*Pair11 + 4*Pair13)*Pair14 + 3*Pair13*Pair9
      Abb533 = -(Pair15*Pair9) + Pair11*(-Pair16 + 2*Pair9)
      Abb529 = Pair10*(Pair12 - Pair14) - Pair13*(Pair12 - 2*Pair14 + Pair9) + 
     -  Pair11*(-Pair14 + Pair9)
      Abb530 = -Abb486 + 4*Pair13*Pair14 - (Pair10 + 3*Pair13)*Pair9 + 
     -  Pair11*(-Pair12 - 3*Pair14 + 2*Pair9)
      Abb514 = Eps41*Pair13 - Eps42*Pair14
      Abb475 = Eps41*(Pair13 - 2*Pair15) - Eps42*(Pair14 - 2*Pair16)
      Abb524 = Eps41*Pair15 - Eps42*Pair16
      Abb476 = Eps41*(2*Pair13 - Pair15) + Eps42*(-2*Pair14 + Pair16)
      Abb489 = Eps49*(Pair13 + Pair15) - Eps50*(Pair14 + Pair16)
      Abb537 = -(Eps31*S34) + Abb536*T24
      Abb544 = Eps29*S34 + Abb536*U
      Abb495 = -(Eps33*Pair12) - Eps34*Pair12 + 3*Eps41*Pair13 + 
     -  Eps36*(Pair10 + 5*Pair11 + Pair13) - Eps33*Pair14 - 
     -  3*Eps42*Pair14 - 3*Eps41*Pair15 + 
     -  Eps35*(Pair10 + 5*Pair11 + Pair15) - Eps34*Pair16 + 
     -  3*Eps42*Pair16 - 5*Eps33*Pair9 - 5*Eps34*Pair9
      Sub94 = -(Abb484*S34) + Abb463*T + Abb486*T14
      Sub161 = Abb533*S34 + Abb535*T - Abb530*T14
      Sub128 = Abb533*S34 + Abb462*T24 - Abb472*U
      Sub166 = Abb509*T14 + Abb508*T24 - 2*(Abb472*T + Abb487*U)
      Sub160 = 2*Abb529*S + Abb485*T + Abb484*T14 + Abb532*T24
      Abb515 = Abb514 + Eps29*T24
      Abb525 = Abb524 - Eps31*U
      Abb496 = Abb495 - 2*Abb455*S - 2*(Eps37 + Eps38)*S34 + 2*Eps30*T + 
     -  2*Eps32*T14 + 2*Eps29*T24 + 2*Eps31*U
      Sub144 = Abb544 - Eps29*(MNeu2(1) - MNeu2(2))
      Sub140 = Abb537 - Eps31*(MNeu2(1) - MNeu2(2))
      Sub189 = (Abb518*T14 + Abb517*T24)*(2*MNeu2(1) - MNeu2(2)) + Sub166*MNeu2(2)
      Sub95 = Abb485*S34*T - Abb487*T**2 - Sub94*T14
      Sub129 = -(Abb532*S34*T24) + Abb482*T24**2 - Sub128*U
      Sub162 = -(S34*Sub160) + (Abb531*T + Abb534*T14)*T24 - Sub161*U
      Sub123 = Abb489*(MNeu2(1) - MNeu2(2)) + Abb500*(5*MNeu2(1) + MNeu2(2)) + 
     -  Abb498*(MNeu2(1) + 5*MNeu2(2))
      Sub188 = Abb516*Sub5 + 3*Abb515*(1 - 4*CW2 + 4*SW2)
      Sub141 = 3*Sub140*(1 - 4*CW2 + 4*SW2) + 5*Pair11*Pair9*Sub5*MNeu2(1)
      Sub190 = -Sub189 + 2*Abb513*(MNeu2(1) - MNeu2(2))
      Sub145 = -(Sub144*MNeu2(1)) + Abb545*(MNeu2(1) + MNeu2(2))
      Sub115 = 4*Pair11*Pair9*Sub11 - 4*Sub95 + Sub98*(MNeu2(1) - MNeu2(2)) + 
     -  Abb497*(MNeu2(1) + 3*MNeu2(2))
      Sub154 = 4*(Pair10*Pair12*Sub11 + Sub129) + Abb548*(MNeu2(1) + 3*MNeu2(2))
      Sub146 = Sub143*Sub5 + 12*Sub145*(1 - 4*CW2 + 4*SW2)
      Sub142 = 2*Sub141*es5(alind) - Sub5*(T + T14)*et5(alind)*MNeu2(1)
      Sub116 = 2*Sub115*es5(alind) + Sub114*et5(alind) + 
     -  4*(2*Pair13*Pair14*es5(alind) + T14*et5(alind))*MNeu2(2)**2 + 
     -  (-2*Sub101*es5(alind) - Sub100*et5(alind))*
     -   (MNeu2(1) + MNeu2(2)) - 
     -  S*T*et5(alind)*(3*MNeu2(1) + 2*MNeu2(2)) - 
     -  S*T14*et5(alind)*(2*MNeu2(1) + 3*MNeu2(2)) + 
     -  2*(S + S34)*et5(alind)*(MNeu2(1)**2 + MNeu2(2)**2)
      Sub151 = Sub150*Sub5 + 2*Sub146*es5(alind) + 4*Sub142*MNeu2(2)
      Sub185 = Sub162 + Abb452*(MNeu2(1) + MNeu2(2))**2 + 
     -  Abb528*(2*MNeu2(1) + MNeu2(2))
      Sub186 = 4*Sub183*es5(alind) + Sub184*et5(alind) + 
     -  4*(Sub185*es5(alind) + S*Sub11*et5(alind))
      Sub157 = 2*Sub154*es5(alind) + 2*Sub156*es5(alind) - 
     -  S*U*et5(alind)*(3*MNeu2(1) + 2*MNeu2(2)) - 
     -  S*T24*et5(alind)*(2*MNeu2(1) + 3*MNeu2(2)) - 
     -  2*Sub153*(MNeu2(1)**2 + MNeu2(2)**2)
      Sub158 = Sub1*Sub157*Sub2 + 2*(Sub152 - Sub22)*et5(alind)
      Sub159 = -(Sub1*Sub151*Sub2) - Sub158*Sub5 + 8*Sub139*es5(alind) - 
     -  2*(Sub1**2 + Sub2**2)*
     -   (4*Sub136*es5(alind) + Sub102*Sub5*et5(alind))*MNeu(1)*MNeu(2)
      Eps54 = epsilonk(k1,k3,k5,ect5)
      Eps53 = epsilonk(k1,k3,k5,es5)
      Eps52 = epsilonk(k1,k4,k5,ect5)
      Eps51 = epsilonk(k1,k4,k5,es5)
      Eps48 = epsilonk(k2,k3,k5,ect5)
      Eps47 = epsilonk(k2,k3,k5,es5)
      Eps46 = epsilonk(k2,k4,k5,ect5)
      Eps44 = epsilonk(k2,k4,k5,es5)
      Opt150 = -2*Eps44*Pair13 + 2*Eps46*Pair14 + 2*Eps47*Pair15 - 2*Eps48*Pair16
      Opt151 = -2*Eps51*Pair13 + 2*Eps52*Pair14 + 2*Eps53*Pair15 - 2*Eps54*Pair16
      Abb521 = Opt151 + 3*Eps41*Pair13 - 3*Eps42*Pair14
      Abb522 = Opt151 + 3*Eps41*Pair15 - 3*Eps42*Pair16
      Abb527 = 8*Abb519*Eps6 - Abb520*S34 + Eps37*S*S34 + Eps38*S*S34 - 
     -  Abb523*T + Eps30*S34*T + Abb526*T14 + Eps32*S34*T14 - 
     -  Abb522*T24 + Eps29*S34*T24 + Eps37*S34*T24 - Abb465*T*T24 + 
     -  Eps29*T24*(-T14 + T24) + Abb521*U + Eps31*S34*U + Eps37*S34*U - 
     -  Eps31*T*U - Abb466*T14*U - Abb455*T24*U + Eps31*U**2
      Sub191 = Abb527 + Eps37*(MNeu2(1) - MNeu2(2))**2 + 
     -  Abb501*(2*MNeu2(1) - MNeu2(2)) + 
     -  Abb525*(2*MNeu2(1) + MNeu2(2)) + 
     -  Abb500*(3*MNeu2(1) + MNeu2(2)) + Abb498*(MNeu2(1) + 3*MNeu2(2))
      Sub192 = Sub190*Sub5 + Sub187*(MNeu2(1) - 2*MNeu2(2)) + 
     -  2*Sub188*(MNeu2(1) + 2*MNeu2(2)) + 
     -  6*(1 - 4*CW2 + 4*SW2)*(-Sub191 + Abb505*MNeu2(2))
      Sub193 = Sub186*Sub5 + 2*Sub192*es5(alind) - 2*Sub182*MNeu2(1) - 
     -  2*Sub181*MNeu2(2)
      Sub194 = Sub1*Sub193*Sub2 - 2*Sub173*es5(alind) + Sub170*MNeu2(1) + 
     -  Sub176*MNeu2(2) + (Sub1*Sub178*Sub2 - 4*Sub180*es5(alind))*
     -   (MNeu2(1) + MNeu2(2))
      Abb477 = Abb490 + Opt150
      Abb478 = Abb501 + Opt150
      Abb479 = Abb474*S34 - Abb455*S*S34 + Abb477*T + Eps30*S34*T - 
     -  Abb478*T14 + Eps32*S34*T14 + 2*Eps37*T*T14 + Abb476*T24 + 
     -  Abb455*S34*T24 + Eps29*(T - T14)*T24 + Abb475*U + 
     -  Abb455*S34*U - Eps31*T*U + Eps31*T14*U
      Sub124 = 2*Abb479 - 2*Sub122 + Sub123 - Abb496*(MNeu2(1) + MNeu2(2))
      Sub125 = Sub121*Sub5 + 6*Sub124*(1 - 4*CW2 + 4*SW2) - 
     -  Sub117*(MNeu2(1) - 3*MNeu2(2)) + 
     -  Sub118*(3*MNeu2(1) - MNeu2(2)) - 
     -  6*Sub119*(MNeu2(1) + MNeu2(2))**2
      Sub126 = Sub116*Sub5 - 2*Sub125*es5(alind) + Sub112*MNeu2(1) + Sub113*MNeu2(2)
      Sub127 = Sub1*Sub126*Sub2 + 2*Sub107*Sub5 + 8*Sub105*es5(alind) + 
     -  2*(Sub1**2 + Sub2**2)*
     -   (4*Sub103*es5(alind) + Sub102*Sub5*et5(alind))*MNeu(1)*MNeu(2)\
     -   - 2*Sub111*(MNeu2(1) + MNeu2(2))

      bmunu(alind,beind,5) = bmunu(alind,beind,5) + dreal((Alfa2*Alfas*Pi**3*Den(S34,MZ2 - ii*MZ*WZ)*
     -    Den(S34,MZ2 + ii*MZ*WZ)*
     -    ((Sub159*Den(3*MU2 - S - T - T14 + MNeu2(1) + MNeu2(2),MU2)**
     -          2)/9. - (2*Sub194*
     -         Den(3*MU2 - S - T - T14 + MNeu2(1) + MNeu2(2),MU2)*
     -         Den(3*MU2 - S - T24 - U + MNeu2(1) + MNeu2(2),MU2))/9. - 
     -      (Sub127*Den(3*MU2 - S - T24 - U + MNeu2(1) + MNeu2(2),MU2)**
     -          2)/9.)*ect5(beind))/(CW2**2*SW2**2))

      enddo
      enddo
      enddo
      enddo

      end
