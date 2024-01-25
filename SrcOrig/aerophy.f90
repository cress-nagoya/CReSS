!***********************************************************************
      module m_aerophy
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/08/18
!     Modification: 2011/09/22, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the inferior procedures for the aerosol processes.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comaero
      use m_comphy
      use m_getiname

      use m_copy3d
      use m_copy4d
      use m_setcst3d

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: aerophy, s_aerophy

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface aerophy

        module procedure s_aerophy

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

!     none

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_aerophy(fpadvopt,fpcphopt,dtb,ni,nj,nk,              &
     &                     nqw,nnw,nqi,nni,nqa,jcb,pbr,ptbr,rbr,rst,    &
     &                     ppp,ptpp,qvp,qwtrp,nwtrp,qicep,nicep,qaslp,  &
     &                     qaslf,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

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

      integer, intent(in) :: nqa(0:4)
                       ! Number of types of aerosol

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: ppp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at past

      real, intent(in) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(in) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(in) :: qwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at past

      real, intent(in) :: nwtrp(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at past

      real, intent(in) :: qicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at past

      real, intent(in) :: nicep(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at past

      real, intent(in) :: qaslp(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at past

! Input and output variable

      real, intent(inout) :: qaslf(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at future

! Internal shared variables

      integer advopt   ! Option for advection scheme
      integer cphopt   ! Option for cloud micro physics

      integer n        ! Array index in 4th direction

      real dtb_sub     ! Substitute for dtb

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp4(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp5(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp6(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp7(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpadvopt,advopt)
      call getiname(fpcphopt,cphopt)

! -----

! Reset the large time steps interval.

      if(advopt.le.3) then
        dtb_sub=2.e0*dtb
      else
        dtb_sub=dtb
      end if

! -----

!! Perform the aerosol processes.

! Now underconstruction, set array with zero or use copy procedure
! to avoid warning.

      call copy3d(0,ni+1,0,nj+1,1,nk,jcb,tmp1)

      call copy3d(0,ni+1,0,nj+1,1,nk,pbr,tmp1)
      call copy3d(0,ni+1,0,nj+1,1,nk,ptbr,tmp1)
      call copy3d(0,ni+1,0,nj+1,1,nk,rbr,tmp1)

      call copy3d(0,ni+1,0,nj+1,1,nk,rst,tmp1)

      call copy3d(0,ni+1,0,nj+1,1,nk,ppp,tmp1)
      call copy3d(0,ni+1,0,nj+1,1,nk,ptpp,tmp1)
      call copy3d(0,ni+1,0,nj+1,1,nk,qvp,tmp1)

      do n=1,nqw

        call s_copy3d(0,ni+1,0,nj+1,1,nk,qwtrp(0,0,1,n),tmp1)

      end do

      do n=1,nnw

        call s_copy3d(0,ni+1,0,nj+1,1,nk,nwtrp(0,0,1,n),tmp1)

      end do

      do n=1,nqi

        call s_copy3d(0,ni+1,0,nj+1,1,nk,qicep(0,0,1,n),tmp1)

      end do

      do n=1,nni

        call s_copy3d(0,ni+1,0,nj+1,1,nk,nicep(0,0,1,n),tmp1)

      end do

      call copy4d(0,ni+1,0,nj+1,1,nk,1,nqa(0),qaslp,qaslf)

      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp1)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp2)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp3)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp4)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp5)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp6)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp7)

! -----

!! -----

      end subroutine s_aerophy

!-----7--------------------------------------------------------------7--

      end module m_aerophy
