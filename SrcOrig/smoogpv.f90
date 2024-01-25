!***********************************************************************
      module m_smoogpv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/03/05
!     Modification: 2004/04/15, 2004/09/10, 2006/01/10, 2006/09/21,
!                   2007/05/07, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2009/03/23, 2011/08/18, 2011/09/22, 2013/02/13,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     smooth the interpolated GPV data.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_getcname
      use m_getiname
      use m_gsmoos
      use m_gsmoou
      use m_gsmoov
      use m_gsmoow
      use m_inichar

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: smoogpv, s_smoogpv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface smoogpv

        module procedure s_smoogpv

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_smoogpv(fpgpvvar,fpcphopt,fphaiopt,fpgsmcnt,         &
     &                     ni,nj,nk,nqw,nqi,ugpv,vgpv,wgpv,ppgpv,       &
     &                     ptpgpv,qvgpv,qwgpv,qigpv,tmp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpgpvvar
                       ! Formal parameter of unique index of gpvvar

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fphaiopt
                       ! Formal parameter of unique index of haiopt

      integer, intent(in) :: fpgsmcnt
                       ! Formal parameter of unique index of gsmcnt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

      integer, intent(in) :: nqi
                       ! Number of categories of ice hydrometeor

! Input and output variables

      real, intent(inout) :: ugpv(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity of GPV data
                       ! at marked time

      real, intent(inout) :: vgpv(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity of GPV data
                       ! at marked time

      real, intent(inout) :: wgpv(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity of GPV data
                       ! at marked time

      real, intent(inout) :: ppgpv(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation of GPV data
                       ! at marked time

      real, intent(inout) :: ptpgpv(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation of GPV data
                       ! at marked time

      real, intent(inout) :: qvgpv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio of GPV data
                       ! at marked time

      real, intent(inout) :: qwgpv(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor of GPV data at marked time

      real, intent(inout) :: qigpv(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor of GPV data at marked time

! Internal shared variables

      character(len=108) gpvvar
                       ! Control flag of input GPV data variables

      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes

      integer gsmcnt   ! Iteration count

      integer iit      ! Index of iteration

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(gpvvar)

! -----

! Get the required namelist variables.

      call getcname(fpgpvvar,gpvvar)
      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)
      call getiname(fpgsmcnt,gsmcnt)

! -----

! Smooth the interpolated GPV data.

      do iit=1,gsmcnt

        call gsmoou(idgsmcoe,ni,nj,nk,ugpv,tmp1)
        call gsmoov(idgsmcoe,ni,nj,nk,vgpv,tmp1)
        call gsmoos(idgsmcoe,ni,nj,nk,ppgpv,tmp1)
        call gsmoos(idgsmcoe,ni,nj,nk,ptpgpv,tmp1)

        if(gpvvar(1:1).eq.'o') then

          call gsmoow(idgsmcoe,ni,nj,nk,wgpv,tmp1)

        end if

        if(gpvvar(2:2).eq.'o') then

          call gsmoos(idgsmcoe,ni,nj,nk,qvgpv,tmp1)

        end if

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            if(gpvvar(3:3).eq.'o') then

              call s_gsmoos(idgsmcoe,ni,nj,nk,qwgpv(0,0,1,1),tmp1)

            end if

            if(gpvvar(4:4).eq.'o') then

              call s_gsmoos(idgsmcoe,ni,nj,nk,qwgpv(0,0,1,2),tmp1)

            end if

          end if

          if(abs(cphopt).ge.2) then

            if(gpvvar(5:5).eq.'o') then

              call s_gsmoos(idgsmcoe,ni,nj,nk,qigpv(0,0,1,1),tmp1)

            end if

            if(gpvvar(6:6).eq.'o') then

              call s_gsmoos(idgsmcoe,ni,nj,nk,qigpv(0,0,1,2),tmp1)

            end if

            if(gpvvar(7:7).eq.'o'.or.gpvvar(8:8).eq.'o') then

              call s_gsmoos(idgsmcoe,ni,nj,nk,qigpv(0,0,1,3),tmp1)

            end if

            if(gpvvar(8:8).eq.'o') then

              if(haiopt.eq.1) then

                call s_gsmoos(idgsmcoe,ni,nj,nk,qigpv(0,0,1,4),tmp1)

              end if

            end if

          end if

        end if

      end do

! -----

      end subroutine s_smoogpv

!-----7--------------------------------------------------------------7--

      end module m_smoogpv
