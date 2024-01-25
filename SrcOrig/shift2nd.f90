!***********************************************************************
      module m_shift2nd
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/05/19
!     Modification: 2003/07/15, 2003/11/05, 2003/11/28, 2003/12/12,
!                   2004/03/05, 2004/05/31, 2004/08/20, 2004/10/12,
!                   2005/02/10, 2005/10/05, 2006/01/10, 2006/02/13,
!                   2006/04/03, 2006/07/21, 2006/09/21, 2006/11/06,
!                   2006/12/04, 2007/01/05, 2007/01/20, 2007/01/31,
!                   2007/11/26, 2008/05/02, 2008/08/25, 2009/01/30,
!                   2009/02/27, 2011/08/18, 2011/09/22, 2011/11/10,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     exchange the value horizontally in the case the 4th order
!     calculation is performed.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bc4news
      use m_bcycle
      use m_combuf
      use m_comindx
      use m_getbufgx
      use m_getbufgy
      use m_getbufsx
      use m_getbufsy
      use m_getiname
      use m_putbufgx
      use m_putbufgy
      use m_putbufsx
      use m_putbufsy
      use m_shiftgx
      use m_shiftgy
      use m_shiftsx
      use m_shiftsy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: shift2nd, s_shift2nd

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface shift2nd

        module procedure s_shift2nd

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_shift2nd(fpwbc,fpebc,fpsbc,fpnbc,                    &
     &                      fpadvopt,fpsmtopt,fpcphopt,fphaiopt,        &
     &                      fpqcgopt,fpaslopt,fptrkopt,fptubopt,        &
     &                      fmois,fproc,ni,nj,nk,nqw,nnw,nqi,nni,nqa,   &
     &                      u,v,w,pp,ptp,qv,qwtr,nwtr,qice,nice,        &
     &                      qcwtr,qcice,qasl,qt,tke)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      character(len=10), intent(in) :: fproc
                       ! Control flag of processing type

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpsbc
                       ! Formal parameter of unique index of sbc

      integer, intent(in) :: fpnbc
                       ! Formal parameter of unique index of nbc

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fpsmtopt
                       ! Formal parameter of unique index of smtopt

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fphaiopt
                       ! Formal parameter of unique index of haiopt

      integer, intent(in) :: fpqcgopt
                       ! Formal parameter of unique index of qcgopt

      integer, intent(in) :: fpaslopt
                       ! Formal parameter of unique index of aslopt

      integer, intent(in) :: fptrkopt
                       ! Formal parameter of unique index of trkopt

      integer, intent(in) :: fptubopt
                       ! Formal parameter of unique index of tubopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

      integer, intent(in) :: nnw
                       ! Number of categories of water concentrations

      integer, intent(in) :: nqi
                       ! Number of categories of ice hydrometeor

      integer, intent(in) :: nni
                       ! Number of categories of ice concentrations

      integer, intent(in) :: nqa(0:4)
                       ! Number of types of aerosol

! Input and output variables

      real, intent(inout) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(inout) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(inout) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

      real, intent(inout) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation

      real, intent(inout) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

      real, intent(inout) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

      real, intent(inout) :: qwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor

      real, intent(inout) :: nwtr(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations

      real, intent(inout) :: qice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor

      real, intent(inout) :: nice(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations

      real, intent(inout) :: qcwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water

      real, intent(inout) :: qcice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice

      real, intent(inout) :: qasl(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio

      real, intent(inout) :: qt(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio

      real, intent(inout) :: tke(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions
      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

      integer advopt   ! Option for advection scheme
      integer smtopt   ! Option for numerical smoothing
      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes
      integer qcgopt   ! Option for charging distribution
      integer aslopt   ! Option for aerosol processes
      integer trkopt   ! Option for mixing ratio tracking
      integer tubopt   ! Option for turbulent mixing

      integer ib       ! Exchanging variables number
      integer nb       ! Number of exchanging variables

      integer n        ! Array index in 4th direction

      integer ni_sub   ! Substitute for ni
      integer nj_sub   ! Substitute for nj

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpsbc,sbc)
      call getiname(fpnbc,nbc)
      call getiname(fpadvopt,advopt)
      call getiname(fpsmtopt,smtopt)
      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)
      call getiname(fpqcgopt,qcgopt)
      call getiname(fpaslopt,aslopt)
      call getiname(fptrkopt,trkopt)
      call getiname(fptubopt,tubopt)

! -----

! Set the substituted variables.

      ni_sub=ni
      nj_sub=nj

! -----

!!!!! Exchange the value horizontally or set the periodic boundary
!!!!! conditions.

      if(mod(wbc,10).eq.6.or.mod(ebc,10).eq.6.or.                       &
     &   mod(sbc,10).eq.6.or.mod(nbc,10).eq.6.or.advopt.ge.2.or.        &
     &   tubopt.ge.1.or.mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

!!!! Exchange the value horizontally.

! Count the number of variables.

        nb=0

        if(fproc(1:1).eq.'o') then
          nb=nb+1
        end if

        if(fproc(2:2).eq.'o') then
          nb=nb+1
        end if

        if(fproc(3:3).eq.'o') then
          nb=nb+1
        end if

        if(fproc(4:4).eq.'o') then
          nb=nb+1
        end if

        if(fproc(5:5).eq.'o') then
          nb=nb+1
        end if

        if(fmois(1:5).eq.'moist'.and.fproc(6:6).eq.'o') then
          nb=nb+1
        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'o') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then
              nb=nb+nqw
            end if

            if(abs(cphopt).ge.2) then

              if(abs(cphopt).eq.2) then
                nb=nb+nqi+1
              else
                nb=nb+nqi+nni
              end if

            end if

            if(abs(cphopt).eq.4) then
              nb=nb+nnw
            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then
                nb=nb+nqw
              end if

              nb=nb+nqi

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then
              nb=nb+nqw+nnw
            end if

            if(abs(cphopt).eq.12) then
              nb=nb+nqi+nni
            end if

          end if

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'-') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then
              nb=nb+nqw-1
            end if

            if(abs(cphopt).ge.2) then
              nb=nb+nqi+nni-2
            end if

            if(abs(cphopt).eq.4) then
              nb=nb+nnw-1
            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then
                nb=nb+nqw-1
              end if

              nb=nb+nqi-1

            end if

          end if

        end if

        if(fproc(8:8).eq.'o'.and.aslopt.ge.1) then
          nb=nb+nqa(0)
        end if

        if(fproc(9:9).eq.'o'.and.trkopt.ge.1) then
          nb=nb+1
        end if

        if(fproc(10:10).eq.'o'.and.tubopt.ge.2) then
          nb=nb+1
        end if

! -----

!!! Exchange the value horizontally between sub domain.

!! In x direction.

! Put the exchanging buffer.

        ib=0

        if(fproc(1:1).eq.'o') then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'all',4,ni-3,ni,nj,nk,u,          &
     &                    ib,nb,sbuf)

        end if

        if(fproc(2:2).eq.'o') then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,v,          &
     &                    ib,nb,sbuf)

        end if

        if(fproc(3:3).eq.'o') then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,w,          &
     &                    ib,nb,sbuf)

        end if

        if(fproc(4:4).eq.'o') then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,pp,         &
     &                    ib,nb,sbuf)

        end if

        if(fproc(5:5).eq.'o') then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,ptp,        &
     &                    ib,nb,sbuf)

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(6:6).eq.'o') then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,qv,         &
     &                    ib,nb,sbuf)

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'o') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qwtr(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qwtr(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nwtr(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nwtr(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qice(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(abs(cphopt).eq.2) then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nice(0,0,1,1),ib,nb,sbuf)

            else if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nice(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          nice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qcwtr(0,0,1,1),ib,nb,sbuf)

                ib=ib+1

                call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qcwtr(0,0,1,2),ib,nb,sbuf)

              end if

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qcice(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qcice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qcice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qcice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                ib=ib+1

                call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qwtr(0,0,1,n),ib,nb,sbuf)

              end do

              do n=1,nnw

                ib=ib+1

                call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          nwtr(0,0,1,n),ib,nb,sbuf)

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                ib=ib+1

                call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qice(0,0,1,n),ib,nb,sbuf)

              end do

              do n=1,nni

                ib=ib+1

                call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          nice(0,0,1,n),ib,nb,sbuf)

              end do

            end if

          end if

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'-') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qwtr(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                       nwtr(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          nice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qcwtr(0,0,1,2),ib,nb,sbuf)

              end if

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qcice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qcice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qcice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

          end if

        end if

        if(fproc(8:8).eq.'o'.and.aslopt.ge.1) then

          do n=1,nqa(0)

            ib=ib+1

            call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,          &
     &                      qasl(0,0,1,n),ib,nb,sbuf)

          end do

        end if

        if(fproc(9:9).eq.'o'.and.trkopt.ge.1) then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,qt,         &
     &                    ib,nb,sbuf)

        end if

        if(fproc(10:10).eq.'o'.and.tubopt.ge.2) then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,tke,        &
     &                    ib,nb,sbuf)

        end if

! -----

! Call the exchanger.

        call s_shiftsx(idwbc,idebc,'all',nj,nk,nb,sbuf,rbuf)

! -----

! Get the exchanging buffer.

        ib=0

        if(fproc(1:1).eq.'o') then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'all',0,ni+1,ni,nj,nk,u,          &
     &                    ib,nb,rbuf)

        end if

        if(fproc(2:2).eq.'o') then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,v,        &
     &                    ib,nb,rbuf)

        end if

        if(fproc(3:3).eq.'o') then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,w,        &
     &                    ib,nb,rbuf)

        end if

        if(fproc(4:4).eq.'o') then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,pp,       &
     &                    ib,nb,rbuf)

        end if

        if(fproc(5:5).eq.'o') then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,ptp,      &
     &                    ib,nb,rbuf)

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(6:6).eq.'o') then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,qv,       &
     &                    ib,nb,rbuf)

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'o') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qwtr(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qwtr(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nwtr(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nwtr(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qice(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(abs(cphopt).eq.2) then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nice(0,0,1,1),ib,nb,rbuf)

            else if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nice(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          nice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qcwtr(0,0,1,1),ib,nb,rbuf)

                ib=ib+1

                call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qcwtr(0,0,1,2),ib,nb,rbuf)

              end if

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qcice(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qcice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qcice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qcice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                ib=ib+1

                call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qwtr(0,0,1,n),ib,nb,rbuf)

              end do

              do n=1,nnw

                ib=ib+1

                call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          nwtr(0,0,1,n),ib,nb,rbuf)

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                ib=ib+1

                call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qice(0,0,1,n),ib,nb,rbuf)

              end do

              do n=1,nni

                ib=ib+1

                call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          nice(0,0,1,n),ib,nb,rbuf)

              end do

            end if

          end if

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'-') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qwtr(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nwtr(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          nice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qcwtr(0,0,1,2),ib,nb,rbuf)

              end if

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qcice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qcice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qcice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

          end if

        end if

        if(fproc(8:8).eq.'o'.and.aslopt.ge.1) then

          do n=1,nqa(0)

            ib=ib+1

            call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,        &
     &                      qasl(0,0,1,n),ib,nb,rbuf)

          end do

        end if

        if(fproc(9:9).eq.'o'.and.trkopt.ge.1) then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,qt,       &
     &                    ib,nb,rbuf)

        end if

        if(fproc(10:10).eq.'o'.and.tubopt.ge.2) then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,tke,      &
     &                    ib,nb,rbuf)

        end if

! -----

!! -----

!! In y direction.

! Put the exchanging buffer.

        ib=0

        if(fproc(1:1).eq.'o') then

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,u,          &
     &                    ib,nb,sbuf)

        end if

        if(fproc(2:2).eq.'o') then

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'all',4,nj-3,ni,nj,nk,v,          &
     &                    ib,nb,sbuf)

        end if

        if(fproc(3:3).eq.'o') then

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,w,          &
     &                    ib,nb,sbuf)

        end if

        if(fproc(4:4).eq.'o') then

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,pp,         &
     &                    ib,nb,sbuf)

        end if

        if(fproc(5:5).eq.'o') then

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,ptp,        &
     &                    ib,nb,sbuf)

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(6:6).eq.'o') then

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,qv,         &
     &                    ib,nb,sbuf)

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'o') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qwtr(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qwtr(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        nwtr(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        nwtr(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qice(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          qice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(abs(cphopt).eq.2) then

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        nice(0,0,1,1),ib,nb,sbuf)

            else if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        nice(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        nice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        nice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          nice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          qcwtr(0,0,1,1),ib,nb,sbuf)

                ib=ib+1

                call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          qcwtr(0,0,1,2),ib,nb,sbuf)

              end if

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qcice(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qcice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qcice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          qcice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                ib=ib+1

                call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          qwtr(0,0,1,n),ib,nb,sbuf)

              end do

              do n=1,nnw

                ib=ib+1

                call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          nwtr(0,0,1,n),ib,nb,sbuf)

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                ib=ib+1

                call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          qice(0,0,1,n),ib,nb,sbuf)

              end do

              do n=1,nni

                ib=ib+1

                call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          nice(0,0,1,n),ib,nb,sbuf)

              end do

            end if

          end if

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'-') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qwtr(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        nwtr(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          qice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        nice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        nice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          nice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          qcwtr(0,0,1,2),ib,nb,sbuf)

              end if

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qcice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qcice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          qcice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

          end if

        end if

        if(fproc(8:8).eq.'o'.and.aslopt.ge.1) then

          do n=1,nqa(0)

            ib=ib+1

            call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,          &
     &                      qasl(0,0,1,n),ib,nb,sbuf)

          end do

        end if

        if(fproc(9:9).eq.'o'.and.trkopt.ge.1) then

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,qt,         &
     &                    ib,nb,sbuf)

        end if

        if(fproc(10:10).eq.'o'.and.tubopt.ge.2) then

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,tke,        &
     &                    ib,nb,sbuf)

        end if

! -----

! Call the exchanger.

        call s_shiftsy(idsbc,idnbc,'all',ni,nk,nb,sbuf,rbuf)

! -----

! Get the exchanging buffer.

        ib=0

        if(fproc(1:1).eq.'o') then

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,u,        &
     &                    ib,nb,rbuf)

        end if

        if(fproc(2:2).eq.'o') then

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'all',0,nj+1,ni,nj,nk,v,          &
     &                    ib,nb,rbuf)

        end if

        if(fproc(3:3).eq.'o') then

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,w,        &
     &                    ib,nb,rbuf)

        end if

        if(fproc(4:4).eq.'o') then

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,pp,       &
     &                    ib,nb,rbuf)

        end if

        if(fproc(5:5).eq.'o') then

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,ptp,      &
     &                    ib,nb,rbuf)

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(6:6).eq.'o') then

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,qv,       &
     &                    ib,nb,rbuf)

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'o') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qwtr(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qwtr(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        nwtr(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        nwtr(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qice(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          qice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(abs(cphopt).eq.2) then

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        nice(0,0,1,1),ib,nb,rbuf)

            else if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        nice(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        nice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        nice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          nice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          qcwtr(0,0,1,1),ib,nb,rbuf)

                ib=ib+1

                call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          qcwtr(0,0,1,2),ib,nb,rbuf)

              end if

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qcice(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qcice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qcice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          qcice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                ib=ib+1

                call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          qwtr(0,0,1,n),ib,nb,rbuf)

              end do

              do n=1,nnw

                ib=ib+1

                call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          nwtr(0,0,1,n),ib,nb,rbuf)

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                ib=ib+1

                call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          qice(0,0,1,n),ib,nb,rbuf)

              end do

              do n=1,nni

                ib=ib+1

                call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          nice(0,0,1,n),ib,nb,rbuf)

              end do

            end if

          end if

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'-') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qwtr(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        nwtr(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          qice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        nice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        nice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          nice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          qcwtr(0,0,1,2),ib,nb,rbuf)

              end if

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qcice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qcice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          qcice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

          end if

        end if

        if(fproc(8:8).eq.'o'.and.aslopt.ge.1) then

          do n=1,nqa(0)

            ib=ib+1

            call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,        &
     &                      qasl(0,0,1,n),ib,nb,rbuf)

          end do

        end if

        if(fproc(9:9).eq.'o'.and.trkopt.ge.1) then

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,qt,       &
     &                    ib,nb,rbuf)

        end if

        if(fproc(10:10).eq.'o'.and.tubopt.ge.2) then

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,tke,      &
     &                    ib,nb,rbuf)

        end if

! -----

!! -----

!!! -----

!!! Exchange the value horizontally between group domain.

!! In x direction.

! Put the exchanging buffer.

        ib=0

        if(fproc(1:1).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'all',4,ni-3,ni,nj,nk,u,          &
     &                    ib,nb,sbuf)

        end if

        if(fproc(2:2).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,v,          &
     &                    ib,nb,sbuf)

        end if

        if(fproc(3:3).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,w,          &
     &                    ib,nb,sbuf)

        end if

        if(fproc(4:4).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,pp,         &
     &                    ib,nb,sbuf)

        end if

        if(fproc(5:5).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,ptp,        &
     &                    ib,nb,sbuf)

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(6:6).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,qv,         &
     &                    ib,nb,sbuf)

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'o') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qwtr(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qwtr(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nwtr(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nwtr(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qice(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(abs(cphopt).eq.2) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nice(0,0,1,1),ib,nb,sbuf)

            else if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nice(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          nice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qcwtr(0,0,1,1),ib,nb,sbuf)

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qcwtr(0,0,1,2),ib,nb,sbuf)

              end if

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qcice(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qcice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qcice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qcice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qwtr(0,0,1,n),ib,nb,sbuf)

              end do

              do n=1,nnw

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          nwtr(0,0,1,n),ib,nb,sbuf)

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qice(0,0,1,n),ib,nb,sbuf)

              end do

              do n=1,nni

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          nice(0,0,1,n),ib,nb,sbuf)

              end do

            end if

          end if

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'-') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qwtr(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nwtr(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          nice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qcwtr(0,0,1,2),ib,nb,sbuf)

              end if

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qcice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qcice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qcice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

          end if

        end if

        if(fproc(8:8).eq.'o'.and.aslopt.ge.1) then

          do n=1,nqa(0)

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,          &
     &                      qasl(0,0,1,n),ib,nb,sbuf)

          end do

        end if

        if(fproc(9:9).eq.'o'.and.trkopt.ge.1) then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,qt,         &
     &                    ib,nb,sbuf)

        end if

        if(fproc(10:10).eq.'o'.and.tubopt.ge.2) then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,tke,        &
     &                    ib,nb,sbuf)

        end if

! -----

! Call the exchanger.

        call s_shiftgx(idwbc,idebc,'all',nj,nk,nb,sbuf,rbuf)

! -----

! Get the exchanging buffer.

        ib=0

        if(fproc(1:1).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'all',0,ni+1,ni,nj,nk,u,          &
     &                    ib,nb,rbuf)

        end if

        if(fproc(2:2).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,v,        &
     &                    ib,nb,rbuf)

        end if

        if(fproc(3:3).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,w,        &
     &                    ib,nb,rbuf)

        end if

        if(fproc(4:4).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,pp,       &
     &                    ib,nb,rbuf)

        end if

        if(fproc(5:5).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,ptp,      &
     &                    ib,nb,rbuf)

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(6:6).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,qv,       &
     &                    ib,nb,rbuf)

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'o') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qwtr(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qwtr(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nwtr(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nwtr(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qice(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(abs(cphopt).eq.2) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nice(0,0,1,1),ib,nb,rbuf)

            else if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nice(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          nice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qcwtr(0,0,1,1),ib,nb,rbuf)

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qcwtr(0,0,1,2),ib,nb,rbuf)

              end if

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qcice(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qcice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qcice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qcice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qwtr(0,0,1,n),ib,nb,rbuf)

              end do

              do n=1,nnw

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          nwtr(0,0,1,n),ib,nb,rbuf)

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qice(0,0,1,n),ib,nb,rbuf)

              end do

              do n=1,nni

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          nice(0,0,1,n),ib,nb,rbuf)

              end do

            end if

          end if

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'-') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qwtr(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nwtr(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          nice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qcwtr(0,0,1,2),ib,nb,rbuf)

              end if

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qcice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qcice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qcice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

          end if

        end if

        if(fproc(8:8).eq.'o'.and.aslopt.ge.1) then

          do n=1,nqa(0)

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,        &
     &                      qasl(0,0,1,n),ib,nb,rbuf)

          end do

        end if

        if(fproc(9:9).eq.'o'.and.trkopt.ge.1) then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,qt,       &
     &                    ib,nb,rbuf)

        end if

        if(fproc(10:10).eq.'o'.and.tubopt.ge.2) then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,tke,      &
     &                    ib,nb,rbuf)

        end if

! -----

!! -----

!! In y direction.

! Put the exchanging buffer.

        ib=0

        if(fproc(1:1).eq.'o') then

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,u,          &
     &                    ib,nb,sbuf)

        end if

        if(fproc(2:2).eq.'o') then

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'all',4,nj-3,ni,nj,nk,v,          &
     &                    ib,nb,sbuf)

        end if

        if(fproc(3:3).eq.'o') then

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,w,          &
     &                    ib,nb,sbuf)

        end if

        if(fproc(4:4).eq.'o') then

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,pp,         &
     &                    ib,nb,sbuf)

        end if

        if(fproc(5:5).eq.'o') then

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,ptp,        &
     &                    ib,nb,sbuf)

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(6:6).eq.'o') then

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,qv,         &
     &                    ib,nb,sbuf)

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'o') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qwtr(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qwtr(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        nwtr(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        nwtr(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qice(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          qice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(abs(cphopt).eq.2) then

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        nice(0,0,1,1),ib,nb,sbuf)

            else if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        nice(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        nice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        nice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          nice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          qcwtr(0,0,1,1),ib,nb,sbuf)

                ib=ib+1

                call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          qcwtr(0,0,1,2),ib,nb,sbuf)

              end if

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qcice(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qcice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qcice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          qcice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                ib=ib+1

                call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          qwtr(0,0,1,n),ib,nb,sbuf)

              end do

              do n=1,nnw

                ib=ib+1

                call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          nwtr(0,0,1,n),ib,nb,sbuf)

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                ib=ib+1

                call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          qice(0,0,1,n),ib,nb,sbuf)

              end do

              do n=1,nni

                ib=ib+1

                call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          nice(0,0,1,n),ib,nb,sbuf)

              end do

            end if

          end if

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'-') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qwtr(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        nwtr(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          qice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        nice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        nice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          nice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          qcwtr(0,0,1,2),ib,nb,sbuf)

              end if

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qcice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,        &
     &                        qcice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,      &
     &                          qcice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

          end if

        end if

        if(fproc(8:8).eq.'o'.and.aslopt.ge.1) then

          do n=1,nqa(0)

            ib=ib+1

            call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,          &
     &                      qasl(0,0,1,n),ib,nb,sbuf)

          end do

        end if

        if(fproc(9:9).eq.'o'.and.trkopt.ge.1) then

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,qt,         &
     &                    ib,nb,sbuf)

        end if

        if(fproc(10:10).eq.'o'.and.tubopt.ge.2) then

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,tke,        &
     &                    ib,nb,sbuf)

        end if

! -----

! Call the exchanger.

        call s_shiftgy(idsbc,idnbc,'all',ni,nk,nb,sbuf,rbuf)

! -----

! Get the exchanging buffer.

        ib=0

        if(fproc(1:1).eq.'o') then

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,u,        &
     &                    ib,nb,rbuf)

        end if

        if(fproc(2:2).eq.'o') then

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'all',0,nj+1,ni,nj,nk,v,          &
     &                    ib,nb,rbuf)

        end if

        if(fproc(3:3).eq.'o') then

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,w,        &
     &                    ib,nb,rbuf)

        end if

        if(fproc(4:4).eq.'o') then

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,pp,       &
     &                    ib,nb,rbuf)

        end if

        if(fproc(5:5).eq.'o') then

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,ptp,      &
     &                    ib,nb,rbuf)

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(6:6).eq.'o') then

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,qv,       &
     &                    ib,nb,rbuf)

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'o') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qwtr(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qwtr(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        nwtr(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        nwtr(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qice(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          qice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(abs(cphopt).eq.2) then

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        nice(0,0,1,1),ib,nb,rbuf)

            else if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        nice(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        nice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        nice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          nice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          qcwtr(0,0,1,1),ib,nb,rbuf)

                ib=ib+1

                call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          qcwtr(0,0,1,2),ib,nb,rbuf)

              end if

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qcice(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qcice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qcice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          qcice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                ib=ib+1

                call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          qwtr(0,0,1,n),ib,nb,rbuf)

              end do

              do n=1,nnw

                ib=ib+1

                call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          nwtr(0,0,1,n),ib,nb,rbuf)

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                ib=ib+1

                call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          qice(0,0,1,n),ib,nb,rbuf)

              end do

              do n=1,nni

                ib=ib+1

                call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          nice(0,0,1,n),ib,nb,rbuf)

              end do

            end if

          end if

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'-') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qwtr(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        nwtr(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          qice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        nice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        nice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          nice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          qcwtr(0,0,1,2),ib,nb,rbuf)

              end if

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qcice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,      &
     &                        qcice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,    &
     &                          qcice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

          end if

        end if

        if(fproc(8:8).eq.'o'.and.aslopt.ge.1) then

          do n=1,nqa(0)

            ib=ib+1

            call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,        &
     &                      qasl(0,0,1,n),ib,nb,rbuf)

          end do

        end if

        if(fproc(9:9).eq.'o'.and.trkopt.ge.1) then

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,qt,       &
     &                    ib,nb,rbuf)

        end if

        if(fproc(10:10).eq.'o'.and.tubopt.ge.2) then

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,tke,      &
     &                    ib,nb,rbuf)

        end if

! -----

!! -----

!! In x direction again.

! Put the exchanging buffer.

        ib=0

        if(fproc(1:1).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'all',4,ni-3,ni,nj,nk,u,          &
     &                    ib,nb,sbuf)

        end if

        if(fproc(2:2).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,v,          &
     &                    ib,nb,sbuf)

        end if

        if(fproc(3:3).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,w,          &
     &                    ib,nb,sbuf)

        end if

        if(fproc(4:4).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,pp,         &
     &                    ib,nb,sbuf)

        end if

        if(fproc(5:5).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,ptp,        &
     &                    ib,nb,sbuf)

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(6:6).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,qv,         &
     &                    ib,nb,sbuf)

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'o') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qwtr(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qwtr(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nwtr(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nwtr(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qice(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(abs(cphopt).eq.2) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nice(0,0,1,1),ib,nb,sbuf)

            else if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nice(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          nice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qcwtr(0,0,1,1),ib,nb,sbuf)

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qcwtr(0,0,1,2),ib,nb,sbuf)

              end if

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qcice(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qcice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qcice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qcice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qwtr(0,0,1,n),ib,nb,sbuf)

              end do

              do n=1,nnw

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          nwtr(0,0,1,n),ib,nb,sbuf)

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qice(0,0,1,n),ib,nb,sbuf)

              end do

              do n=1,nni

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          nice(0,0,1,n),ib,nb,sbuf)

              end do

            end if

          end if

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'-') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qwtr(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nwtr(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        nice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          nice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qcwtr(0,0,1,2),ib,nb,sbuf)

              end if

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qcice(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,        &
     &                        qcice(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,      &
     &                          qcice(0,0,1,4),ib,nb,sbuf)

              end if

            end if

          end if

        end if

        if(fproc(8:8).eq.'o'.and.aslopt.ge.1) then

          do n=1,nqa(0)

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,          &
     &                      qasl(0,0,1,n),ib,nb,sbuf)

          end do

        end if

        if(fproc(9:9).eq.'o'.and.trkopt.ge.1) then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,qt,         &
     &                    ib,nb,sbuf)

        end if

        if(fproc(10:10).eq.'o'.and.tubopt.ge.2) then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,tke,        &
     &                    ib,nb,sbuf)

        end if

! -----

! Call the exchanger.

        call s_shiftgx(idwbc,idebc,'all',nj,nk,nb,sbuf,rbuf)

! -----

! Get the exchanging buffer.

        ib=0

        if(fproc(1:1).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'all',0,ni+1,ni,nj,nk,u,          &
     &                    ib,nb,rbuf)

        end if

        if(fproc(2:2).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,v,        &
     &                    ib,nb,rbuf)

        end if

        if(fproc(3:3).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,w,        &
     &                    ib,nb,rbuf)

        end if

        if(fproc(4:4).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,pp,       &
     &                    ib,nb,rbuf)

        end if

        if(fproc(5:5).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,ptp,      &
     &                    ib,nb,rbuf)

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(6:6).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,qv,       &
     &                    ib,nb,rbuf)

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'o') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qwtr(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qwtr(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nwtr(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nwtr(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qice(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(abs(cphopt).eq.2) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nice(0,0,1,1),ib,nb,rbuf)

            else if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nice(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          nice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qcwtr(0,0,1,1),ib,nb,rbuf)

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qcwtr(0,0,1,2),ib,nb,rbuf)

              end if

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qcice(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qcice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qcice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qcice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qwtr(0,0,1,n),ib,nb,rbuf)

              end do

              do n=1,nnw

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          nwtr(0,0,1,n),ib,nb,rbuf)

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qice(0,0,1,n),ib,nb,rbuf)

              end do

              do n=1,nni

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          nice(0,0,1,n),ib,nb,rbuf)

              end do

            end if

          end if

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'-') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qwtr(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nwtr(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        nice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          nice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qcwtr(0,0,1,2),ib,nb,rbuf)

              end if

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qcice(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,      &
     &                        qcice(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,    &
     &                          qcice(0,0,1,4),ib,nb,rbuf)

              end if

            end if

          end if

        end if

        if(fproc(8:8).eq.'o'.and.aslopt.ge.1) then

          do n=1,nqa(0)

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,        &
     &                      qasl(0,0,1,n),ib,nb,rbuf)

          end do

        end if

        if(fproc(9:9).eq.'o'.and.trkopt.ge.1) then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,qt,       &
     &                    ib,nb,rbuf)

        end if

        if(fproc(10:10).eq.'o'.and.tubopt.ge.2) then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,tke,      &
     &                    ib,nb,rbuf)

        end if

! -----

!! -----

!!! -----

!!!! -----

! Set the periodic boundary conditions.

        if(fproc(1:1).eq.'o') then

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                4,0,ni-3,ni+1,3,0,nj-3,nj_sub,ni,nj,nk,u)

        end if

        if(fproc(2:2).eq.'o') then

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                3,0,ni-3,ni_sub,4,0,nj-3,nj+1,ni,nj,nk,v)

        end if

        if(fproc(3:3).eq.'o') then

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,ni,nj,nk,w)

        end if

        if(fproc(4:4).eq.'o') then

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,ni,nj,nk,pp)

        end if

        if(fproc(5:5).eq.'o') then

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,ni,nj,nk,ptp)

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(6:6).eq.'o') then

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,ni,nj,nk,qv)

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'o') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,            &
     &                      ni,nj,nk,qwtr(0,0,1,1))

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,            &
     &                      ni,nj,nk,qwtr(0,0,1,2))

            end if

            if(abs(cphopt).eq.4) then

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,            &
     &                      ni,nj,nk,nwtr(0,0,1,1))

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,            &
     &                      ni,nj,nk,nwtr(0,0,1,2))

            end if

            if(abs(cphopt).ge.2) then

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,            &
     &                      ni,nj,nk,qice(0,0,1,1))

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,            &
     &                      ni,nj,nk,qice(0,0,1,2))

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,            &
     &                      ni,nj,nk,qice(0,0,1,3))

              if(haiopt.eq.1) then

                call s_bcycle(idwbc,idebc,idsbc,idnbc,                  &
     &                        3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,          &
     &                        ni,nj,nk,qice(0,0,1,4))

              end if

            end if

            if(abs(cphopt).eq.2) then

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,            &
     &                      ni,nj,nk,nice(0,0,1,1))

            else if(abs(cphopt).ge.3) then

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,            &
     &                      ni,nj,nk,nice(0,0,1,1))

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,            &
     &                      ni,nj,nk,nice(0,0,1,2))

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,            &
     &                      ni,nj,nk,nice(0,0,1,3))

              if(haiopt.eq.1) then

                call s_bcycle(idwbc,idebc,idsbc,idnbc,                  &
     &                        3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,          &
     &                        ni,nj,nk,nice(0,0,1,4))

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                call s_bcycle(idwbc,idebc,idsbc,idnbc,                  &
     &                        3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,          &
     &                        ni,nj,nk,qcwtr(0,0,1,1))

                call s_bcycle(idwbc,idebc,idsbc,idnbc,                  &
     &                        3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,          &
     &                        ni,nj,nk,qcwtr(0,0,1,2))

              end if

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,            &
     &                      ni,nj,nk,qcice(0,0,1,1))

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,            &
     &                      ni,nj,nk,qcice(0,0,1,2))

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,            &
     &                      ni,nj,nk,qcice(0,0,1,3))

              if(haiopt.eq.1) then

                call s_bcycle(idwbc,idebc,idsbc,idnbc,                  &
     &                        3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,          &
     &                        ni,nj,nk,qcice(0,0,1,4))

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                call s_bcycle(idwbc,idebc,idsbc,idnbc,                  &
     &                        3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,          &
     &                        ni,nj,nk,qwtr(0,0,1,n))

              end do

              do n=1,nnw

                call s_bcycle(idwbc,idebc,idsbc,idnbc,                  &
     &                        3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,          &
     &                        ni,nj,nk,nwtr(0,0,1,n))

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                call s_bcycle(idwbc,idebc,idsbc,idnbc,                  &
     &                        3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,          &
     &                        ni,nj,nk,qice(0,0,1,n))

              end do

              do n=1,nni

                call s_bcycle(idwbc,idebc,idsbc,idnbc,                  &
     &                        3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,          &
     &                        ni,nj,nk,nice(0,0,1,n))

              end do

            end if

          end if

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'-') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,            &
     &                      ni,nj,nk,qwtr(0,0,1,2))

            end if

            if(abs(cphopt).eq.4) then

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,            &
     &                      ni,nj,nk,nwtr(0,0,1,2))

            end if

            if(abs(cphopt).ge.2) then

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,            &
     &                      ni,nj,nk,qice(0,0,1,2))

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,            &
     &                      ni,nj,nk,qice(0,0,1,3))

              if(haiopt.eq.1) then

                call s_bcycle(idwbc,idebc,idsbc,idnbc,                  &
     &                        3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,          &
     &                        ni,nj,nk,qice(0,0,1,4))

              end if

            end if

            if(abs(cphopt).ge.3) then

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,            &
     &                      ni,nj,nk,nice(0,0,1,2))

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,            &
     &                      ni,nj,nk,nice(0,0,1,3))

              if(haiopt.eq.1) then

                call s_bcycle(idwbc,idebc,idsbc,idnbc,                  &
     &                        3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,          &
     &                        ni,nj,nk,nice(0,0,1,4))

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                call s_bcycle(idwbc,idebc,idsbc,idnbc,                  &
     &                        3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,          &
     &                        ni,nj,nk,qcwtr(0,0,1,2))

              end if

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,            &
     &                      ni,nj,nk,qcice(0,0,1,2))

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,            &
     &                      ni,nj,nk,qcice(0,0,1,3))

              if(haiopt.eq.1) then

                call s_bcycle(idwbc,idebc,idsbc,idnbc,                  &
     &                        3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,          &
     &                        ni,nj,nk,qcice(0,0,1,4))

              end if

            end if

          end if

        end if

        if(fproc(8:8).eq.'o'.and.aslopt.ge.1) then

          do n=1,nqa(0)

            call s_bcycle(idwbc,idebc,idsbc,idnbc,                      &
     &                    3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,              &
     &                    ni,nj,nk,qasl(0,0,1,n))

          end do

        end if

        if(fproc(9:9).eq.'o'.and.trkopt.ge.1) then

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,ni,nj,nk,qt)

        end if

        if(fproc(10:10).eq.'o'.and.tubopt.ge.2) then

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,ni,nj,nk,tke)

        end if

! -----

! Set the boundary conditions at the four corners.

        if(fproc(1:1).eq.'o') then

          call bc4news(idwbc,idebc,idsbc,idnbc,0,ni+1,0,nj_sub,         &
     &                 ni,nj,nk,u)

        end if

        if(fproc(2:2).eq.'o') then

          call bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj+1,         &
     &                 ni,nj,nk,v)

        end if

        if(fproc(3:3).eq.'o') then

          call bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub,       &
     &                 ni,nj,nk,w)

        end if

        if(fproc(4:4).eq.'o') then

          call bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub,       &
     &                 ni,nj,nk,pp)

        end if

        if(fproc(5:5).eq.'o') then

          call bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub,       &
     &                 ni,nj,nk,ptp)

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(6:6).eq.'o') then

          call bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub,       &
     &                 ni,nj,nk,qv)

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'o') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              call s_bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub, &
     &                       ni,nj,nk,qwtr(0,0,1,1))

              call s_bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub, &
     &                       ni,nj,nk,qwtr(0,0,1,2))

            end if

            if(abs(cphopt).eq.4) then

              call s_bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub, &
     &                       ni,nj,nk,nwtr(0,0,1,1))

              call s_bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub, &
     &                       ni,nj,nk,nwtr(0,0,1,2))

            end if

            if(abs(cphopt).ge.2) then

              call s_bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub, &
     &                       ni,nj,nk,qice(0,0,1,1))

              call s_bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub, &
     &                       ni,nj,nk,qice(0,0,1,2))

              call s_bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub, &
     &                       ni,nj,nk,qice(0,0,1,3))

              if(haiopt.eq.1) then

                call s_bc4news(idwbc,idebc,idsbc,idnbc,                 &
     &                         0,ni_sub,0,nj_sub,ni,nj,nk,qice(0,0,1,4))

              end if

            end if

            if(abs(cphopt).eq.2) then

              call s_bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub, &
     &                       ni,nj,nk,nice(0,0,1,1))

            else if(abs(cphopt).ge.3) then

              call s_bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub, &
     &                       ni,nj,nk,nice(0,0,1,1))

              call s_bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub, &
     &                       ni,nj,nk,nice(0,0,1,2))

              call s_bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub, &
     &                       ni,nj,nk,nice(0,0,1,3))

              if(haiopt.eq.1) then

                call s_bc4news(idwbc,idebc,idsbc,idnbc,                 &
     &                         0,ni_sub,0,nj_sub,ni,nj,nk,nice(0,0,1,4))

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                call s_bc4news(idwbc,idebc,idsbc,idnbc,                 &
     &                        0,ni_sub,0,nj_sub,ni,nj,nk,qcwtr(0,0,1,1))

                call s_bc4news(idwbc,idebc,idsbc,idnbc,                 &
     &                        0,ni_sub,0,nj_sub,ni,nj,nk,qcwtr(0,0,1,2))

              end if

              call s_bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub, &
     &                       ni,nj,nk,qcice(0,0,1,1))

              call s_bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub, &
     &                       ni,nj,nk,qcice(0,0,1,2))

              call s_bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub, &
     &                       ni,nj,nk,qcice(0,0,1,3))

              if(haiopt.eq.1) then

                call s_bc4news(idwbc,idebc,idsbc,idnbc,                 &
     &                        0,ni_sub,0,nj_sub,ni,nj,nk,qcice(0,0,1,4))

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                call s_bc4news(idwbc,idebc,idsbc,idnbc,                 &
     &                         0,ni_sub,0,nj_sub,ni,nj,nk,qwtr(0,0,1,n))

              end do

              do n=1,nnw

                call s_bc4news(idwbc,idebc,idsbc,idnbc,                 &
     &                         0,ni_sub,0,nj_sub,ni,nj,nk,nwtr(0,0,1,n))

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                call s_bc4news(idwbc,idebc,idsbc,idnbc,                 &
     &                         0,ni_sub,0,nj_sub,ni,nj,nk,qice(0,0,1,n))

              end do

              do n=1,nni

                call s_bc4news(idwbc,idebc,idsbc,idnbc,                 &
     &                         0,ni_sub,0,nj_sub,ni,nj,nk,nice(0,0,1,n))

              end do

            end if

          end if

        end if

        if(fmois(1:5).eq.'moist'.and.fproc(7:7).eq.'-') then

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              call s_bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub, &
     &                       ni,nj,nk,qwtr(0,0,1,2))

            end if

            if(abs(cphopt).eq.4) then

              call s_bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub, &
     &                       ni,nj,nk,nwtr(0,0,1,2))

            end if

            if(abs(cphopt).ge.2) then

              call s_bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub, &
     &                       ni,nj,nk,qice(0,0,1,2))

              call s_bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub, &
     &                       ni,nj,nk,qice(0,0,1,3))

              if(haiopt.eq.1) then

                call s_bc4news(idwbc,idebc,idsbc,idnbc,                 &
     &                         0,ni_sub,0,nj_sub,ni,nj,nk,qice(0,0,1,4))

              end if

            end if

            if(abs(cphopt).ge.3) then

              call s_bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub, &
     &                       ni,nj,nk,nice(0,0,1,2))

              call s_bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub, &
     &                       ni,nj,nk,nice(0,0,1,3))

              if(haiopt.eq.1) then

                call s_bc4news(idwbc,idebc,idsbc,idnbc,                 &
     &                         0,ni_sub,0,nj_sub,ni,nj,nk,nice(0,0,1,4))

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                call s_bc4news(idwbc,idebc,idsbc,idnbc,                 &
     &                        0,ni_sub,0,nj_sub,ni,nj,nk,qcwtr(0,0,1,2))

              end if

              call s_bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub, &
     &                       ni,nj,nk,qcice(0,0,1,2))

              call s_bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub, &
     &                       ni,nj,nk,qcice(0,0,1,3))

              if(haiopt.eq.1) then

                call s_bc4news(idwbc,idebc,idsbc,idnbc,                 &
     &                        0,ni_sub,0,nj_sub,ni,nj,nk,qcice(0,0,1,4))

              end if

            end if

          end if

        end if

        if(fproc(8:8).eq.'o'.and.aslopt.ge.1) then

          do n=1,nqa(0)

            call s_bc4news(idwbc,idebc,idsbc,idnbc,                     &
     &                     0,ni_sub,0,nj_sub,ni,nj,nk,qasl(0,0,1,n))

          end do

        end if

        if(fproc(9:9).eq.'o'.and.trkopt.ge.1) then

          call bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub,       &
     &                 ni,nj,nk,qt)

        end if

        if(fproc(10:10).eq.'o'.and.tubopt.ge.2) then

          call bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub,       &
     &                 ni,nj,nk,tke)

        end if

! -----

      end if

!!!!! -----

      end subroutine s_shift2nd

!-----7--------------------------------------------------------------7--

      end module m_shift2nd
