!***********************************************************************
      module m_satadjst
!***********************************************************************

!     Author      : Sakakibara Atsushi, Naito Daisuke
!     Date        : 2004/12/17
!     Modification: 2005/04/04, 2005/10/05, 2006/01/10, 2006/02/13,
!                   2006/09/30, 2007/11/26, 2008/05/02, 2008/07/01,
!                   2008/08/25, 2009/02/27, 2011/03/29, 2011/08/18,
!                   2013/02/13

!     Author      : Satoki Tsujino
!     Modification: 2024/12/25

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the inferior procedures for the saturation adjustment.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_dmpcld
      use m_getiname
      use m_getexner
      use m_siadjst
      use m_swadjst

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: satadjst, s_satadjst

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface satadjst

        module procedure s_satadjst

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
      subroutine s_satadjst(fpcphopt,ni,nj,nk,nqw,nnw,nqi,nni,pbr,ptbr, &
     &                      w,pp,ptp,qv,qwtr,nwtr,qice,nice,pi,p)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

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

      real, intent(in) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

      real, intent(in) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Presure perturbation

! Input and output variables

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

! Internal shared variables

      integer cphopt   ! Option for cloud micro physics

      real, intent(inout) :: pi(0:ni+1,0:nj+1,1:nk)
                       ! Exnar function

      real, intent(inout) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure

! Internal private variables

      real :: tmpdmp(0:ni+1,0:nj+1,1:nk)

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpcphopt,cphopt)

! -----

! Setting ptp to vdvc (by satoki)

      call s_dmpcld( ni, nj, nk, ptp, 'vdvc ', 'o', 1.0e0 )

! -----

!! Perform the saturation adjustment.

      if(abs(cphopt).ge.1) then

        if(abs(cphopt).lt.10) then

! Calculate the total pressure variable and Exner function.

          call getexner(ni,nj,nk,pbr,pp,pi,p)

! -----

! For the bulk warm rain cloud physics.

          if(abs(cphopt).eq.1) then

            call s_swadjst(idcphopt,idthresq,ni,nj,nk,ptbr,pi,w,p,      &
     &                     ptp,qv,qwtr(0,0,1,1),nwtr(0,0,1,1))

! -----

! For the bulk cold rain cloud physics.

          else if(abs(cphopt).ge.2) then

            call s_siadjst(idthresq,ni,nj,nk,ptbr,pi,p,ptp,qv,          &
     &                     qice(0,0,1,1),nice(0,0,1,1))

            call s_swadjst(idcphopt,idthresq,ni,nj,nk,ptbr,pi,w,p,      &
     &                     ptp,qv,qwtr(0,0,1,1),nwtr(0,0,1,1))

          end if

! -----

        end if

      end if

!! -----

! Setting ptp to vdvc and adding vdvc to dqdt and dqadt (by satoki)

      call s_dmpcld( ni, nj, nk, ptp, 'vdvc ', 'm', 1.0e0 )
      call s_dmpcld( ni, nj, nk, tmpdmp, 'vdvc ', 'r', 1.0e0 )
      call s_dmpcld( ni, nj, nk, tmpdmp, 'dqdt ', 'p', 1.0e0 )
      call s_dmpcld( ni, nj, nk, tmpdmp, 'dqadt', 'p', 1.0e0 )

! -----

      end subroutine s_satadjst

!-----7--------------------------------------------------------------7--

      end module m_satadjst
