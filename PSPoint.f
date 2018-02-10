c Copyright (C) Matthias Kesenheimer - All Rights Reserved
c Written by Matthias Kesenheimer <m.kesenheimer@gmx.net>, 2017

!#######################################################################

!Dieses Programm wertet einen zufälligen Phasenraumpunkt von squaredME.a
!und renconst.a aus und gibt das Ergebnis aus.

!kompilieren (Mac OS X): "gfortran -O3 -g -ffixed-line-length-none -march=native -Wa,-q -ff2c -cpp -dM -L/Users/kesenheimer/Applications/LoopTools/x86_64-Darwin/lib -looptools renconst.a squaredme.a testLib.o PSPoint.f -o PSPoint"

!kompilieren (Linux): "gfortran -O3 -g -ffixed-line-length-none -march=native -Wa,-q -ff2c -cpp -dM -I/usr/local/LoopTools/x86_64-Linux/include/ PSPoint.f squaredme.a renconst.a /usr/local/LoopTools/x86_64-Linux/lib64/libooptools.a -o PSPoint"

!oder CMakeLists.txt anpassen und cmake . && make ausführen

!#####################Program Start#####################################

      program readPS
	
        implicit none

#include "looptools.h"
#include "const.h"
#include "model_sm.h"
#include "model_mssm.h"
#include "Functions.h"
#include "phaseSpace.h"

!#####################Define Variables##################################	

!       Mandelstam variables
        double precision s, t, u
!       phase space limits
        double precision tp, tm

!       variables for squaredME
        double precision result(2)
        integer*8 helicities
        integer flags
        
!       auxiliary variables
        integer i       !loop counter
        double precision x

!       console variables
        character*512 argv(5)
        integer argc

!       write in file
        integer io_error
        integer n

!       DEBUG variables

!       Define Tab for output
        character*1 tab
        tab = char(9)
        
        
!#####################Calculate Model Parameters########################

        write(*,*) "###################################################"
        write(*,*) "#                                                 #"
        write(*,*) "#    Calculate SquaredME at a single PS Point     #"
        write(*,*) "#                                                 #"
        write(*,*) "###################################################"
        write(*,*) ""

        call setSMParameters()
        call setMSSMParameters()        

!#####################Console Magic#####################################
        
        argc = iargc()
        
        if(argc < 1) then
          print*,"To evaluate the matrix element in a single PS point:"
          print*,"-> ./PSPoint s [GeV**2], t [GeV**2]"
          print*,"To evaluate the matrix element in a PS range:"
          print*,"-> ./PSPoint s [GeV**2]"
          stop
        end if

!       read CMS energy
        call getarg(1,argv(1))
        read(argv(1),*) s      !read s
        
!       get phase space limits
        tm = getLimitTm(s, MU, MU, MNeu(1), MNeu(1))
        tp = getLimitTp(s, MU, MU, MNeu(1), MNeu(1))
        
!       read t
        if(argc .eq. 2) then
          call getarg(2,argv(2))
          read(argv(2),*) t
          if(.not. (t <= tp .and. t >= tm)) then
            write(*,*) "t should be in [tm:tp] = [", tm, "," , tp, "]"
            stop
          end if
        end if
        
        call printModelParameters()

!#####################Main##############################################

!       init of LoopTools
        call ltini
        call setversionkey(0)

!       unpolarized particles: B01010 01010 01010 01010 = D338250
        helicities = 338250.D0
!       flags: Bit0 (reset) = 1, Bit1 (loop) = 1 -> B11 = D3
        flags = 3

        if(argc .eq. 2) then
!         set up phase space
          call setUpPhaseSpace(s, t, u, MU, MU, MNeu(1), MNeu(1))
          print*, "bla"

!         calculate renormalization consts and matrix element
          call clearcache
          call calcRenConst
          result(1) = 0.D0
          result(2) = 0.D0
          call squaredME(result, helicities, flags)

!         output result
          print*
          print*, "== Results =="
          print*, "result =", result(1), result(2)
          print*

!         test for finiteness, lambda operates also as IR-cutoff
          call setMuDim(1D100)
          ! call setLambda(1D0) !Test IR finiteness
          call setDelta(1D7)    !Test UV finiteness

!         calculate renormalization consts and matrix element
          call clearcache
          call calcRenConst
          result(1) = 0.D0
          result(2) = 0.D0
          call squaredME(result, helicities, flags)

!         output result
          print*
          print*, "== UV Check =="
          print*, "result =", result(1), result(2)
          print*

        end if

        if(argc .eq. 1) then

          print*
          print*, "== Results =="
          print*
        
          open(unit=20,file="results.txt",status="new",action="write",
     &         iostat=io_error)
          write(20,*) "Sqrt(s) = ", dsqrt(s)
          write(20,*) "s       = ", s
          write(20,*)
          write(20,"(A,5a,A,5a)") "t",tab,"u",tab,"tree",tab,"loop"
          write(20,*)

          x = 0.00001     !regulate

          do i=1,100

            t = tm - x*(tm - tp)
            call setUpPhaseSpace(s, t, u, MU, MU, MNeu(1), MNeu(1))
            
            call clearcache
            call calcRenConst
            
            result(1) = 0.D0
            result(2) = 0.D0
            call squaredME(result, helicities, flags)

            print*, "result =", result(1), result(2)
            if(io_error .eq. 0) then
              write(20,*) t, u, result(1), result(2)
            end if
          
            x = x + 0.01

          end do

          print*
          close(unit=20)

        end if

      end program readPS

!#######################################################################

