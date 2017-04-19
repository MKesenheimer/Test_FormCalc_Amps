c############### Main.F ################################################
      program main
        implicit none
        integer nlegborn,ndiminteg
        parameter (nlegborn=5)
        parameter (ndiminteg=(nlegborn-1)*3-4+2-1)
        double precision p(0:3,nlegborn)
        double precision xborn(ndiminteg-3)
        integer flav(nlegborn), i, j
        double precision amp1, amp2
        double precision bmunu(0:3,0:3,nlegborn)
        
        ! generate random number
        call srand(1234)
        
        ! process parameters
        !flav(1) = -2
        !flav(2) = 2
        !flav(3) = -2
        !flav(4) = 2
        flav(1) = -2
        flav(2) = 2
        flav(3) = 1000022
        flav(4) = 1000023
        flav(5) = 0

        call init_couplings
        
        do i=1,5
          ! choose a point
          do j=1,ndiminteg-3
            xborn(j) = rand()
          enddo
        
          if(nlegborn.eq.4) call born_phsp4(xborn,flav,p)
          if(nlegborn.eq.5) call born_phsp5(xborn,flav,p)
        
#if DEBUG >= 0
          print*
          print*,"p(:,:) = "
          do j=1,nlegborn
            print*,p(:,j)
          enddo
#endif
          call formCalcAmp(p,flav,amp1)

          print*,"born = ", amp1
          
          ! test the spin correlated amplitude
          !call bmunu_ubaru_gamg(p,bmunu)
          call bmunu_ubaru_n1n2g(p,bmunu)
          bmunu(:,:,:) = bmunu(:,:,:) / 9.
          print*,"bmunu(0) = ", bmunu(0,:,5)
          print*,"bmunu(1) = ", bmunu(1,:,5)
          print*,"bmunu(2) = ", bmunu(2,:,5)
          print*,"bmunu(3) = ", bmunu(3,:,5)
          print*,"gmunu*bmunu = ",bmunu(0,0,5)-bmunu(1,1,5)-bmunu(2,2,5)-bmunu(3,3,5)
        enddo  
      end
c############### end Main ##############################################

c############### subroutine formCalcVecSet #############################
c calculate unit vectors and set the vectors for FormCalc
      subroutine formCalcVecSet(nlegs,p,flav)
        implicit none
        integer nlegs
        double precision p(0:3,nlegs)
        integer  flav(nlegs)
        double precision mass(nlegs)
        double precision ex, ey, ez, pNorm
        integer i
        ! function to calculate the squared momentum sum
        double precision momsum3sq, momsum2sq, dotp
        external momsum2sq, momsum3sq, dtop
        ! load process dependent masses
        call set_process(nlegs,flav(1:nlegs),mass(1:nlegs))
        do i=1,nlegs
          pNorm=dsqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2)
          ex=p(1,i)/pNorm
          ey=p(2,i)/pNorm
          ez=p(3,i)/pNorm
          call VecSet(i, mass(i), pNorm, ex,ey,ez)
        enddo
      end
c############### end subroutine formCalcVecSet #########################

c############### subroutine getAmps ####################################
      subroutine formCalcAmp(p,flav,amp)
        implicit none
#include "PhysPars.h"
        integer nlegborn,ndiminteg
        parameter (nlegborn=5)
        parameter (ndiminteg=(nlegborn-1)*3-4+2-1)
        double precision p(0:3,nlegborn)
        ! helicity and flags parameter for formcalc       
        integer*8 helicities
        integer flags
        integer  flav(nlegborn)
        ! formcalc results
        double precision fc_result(2), amp, symfac

        ! reset
        fc_result(:) = 0D0
        amp = 0D0
        
        call set_ebe_couplings
        
        ! process parameters
        ! unpolarized fermion: B01010  FC8Guid.pdf p.49
        ! unpolarized gluon: B01110
        !helicities = 338382
        helicities = 10824014
        flags = 2 ! = 010d, Bit0: set mass, Bit1: reset, Bit2: loop
        symfac = 36
        ! set the formCalc vectors
        call formCalcVecSet(nlegborn,p,flav)
        
        !print*,MW,MZ,WZ,WW
        !print*,Alfa,Alfa2,AlfaS,AlfaS2
        
        call squaredME(fc_result,helicities,flags)
        amp = fc_result(1)/symfac
      end
c############### subroutine getAmps ####################################
