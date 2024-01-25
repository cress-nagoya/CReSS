!***********************************************************************
      module m_inibinbw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/09/30
!     Modification: 2007/03/29, 2007/05/07, 2008/01/11, 2008/05/02,
!                   2008/07/01, 2008/08/25, 2009/01/30, 2009/02/27,
!                   2011/08/18

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the inferior procedures for the initial bin distribution
!     of total water.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_inidisbw
      use m_sadjstbw
      use m_sumbin

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: inibinbw, s_inibinbw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface inibinbw

        module procedure s_inibinbw

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
      subroutine s_inibinbw(dtb,ni,nj,nk,nqw,nnw,jcb8w,ptbr,rbv,        &
     &                      pi,p,w,t,brw,rbrw,ptp,qv,mwbin,nwbin,qwtr,  &
     &                      ptptmp,qvtmp,qwtmp,tmp1,tmp2,tmp3,tmp4,     &
     &                      tmp5,tmp6,tmp7,tmp8,tmp9,tmp10)
!***********************************************************************

! Input variables

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

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: rbv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of base state density [cm^3/g]

      real, intent(in) :: pi(0:ni+1,0:nj+1,1:nk)
                       ! Exnar function

      real, intent(in) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure [dyn/cm^2]

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

      real, intent(in) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

      real, intent(in) :: brw(1:nqw+1)
                       ! Radius of water bin boundary [cm]

      real, intent(in) :: rbrw(1:nqw+1,1:8)
                       ! Related parameters of brw

! Input and output variables

      real, intent(inout) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

      real, intent(inout) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio [g/g]

      real, intent(inout) :: mwbin(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Total water mass [g/cm^3]

      real, intent(inout) :: nwbin(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations [1/cm^3]

! Internal shared variables

      real, intent(inout) :: qwtr(0:ni+1,0:nj+1,1:nk)
                       ! Total water mixing ratio

      real, intent(inout) :: ptptmp(0:ni+1,0:nj+1,1:nk)
                       ! Temporary alternative
                       ! potential temperature perturbation

      real, intent(inout) :: qvtmp(0:ni+1,0:nj+1,1:nk)
                       ! Temporary alternative
                       ! water vapor mixing ratio

      real, intent(inout) :: qwtmp(0:ni+1,0:nj+1,1:nk)
                       ! Temporary alternative
                       ! total water mixing ratio

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp4(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp5(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp6(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp7(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp8(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp9(0:ni+1,0:nj+1,1:nqw+1)
                       ! Temporary array

      real, intent(inout) :: tmp10(0:ni+1,0:nj+1,1:nqw+1)
                       ! Temporary array

!-----7--------------------------------------------------------------7--

! Get the total mixing ratio.

      call sumbin(ni,nj,nk,nqw,rbv,mwbin,qwtr)

! -----

! Perform the saturation adjustment for water in the case the generated
! total water is more than critical value.

      call sadjstbw(idthresq,ni,nj,nk,ptbr,pi,p,w,t,ptp,qv,qwtr,        &
     &              ptptmp,qvtmp,qwtmp)

! -----

! Set the initial bin distribution of total water.

      call inidisbw(iddziv,dtb,ni,nj,nk,nqw,nnw,jcb8w,rbv,pi,w,t,       &
     &              brw,rbrw,ptptmp,qvtmp,qwtmp,ptp,qv,mwbin,nwbin,     &
     &              tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10)

! -----

      end subroutine s_inibinbw

!-----7--------------------------------------------------------------7--

      end module m_inibinbw
