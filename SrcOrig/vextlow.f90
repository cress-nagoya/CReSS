!***********************************************************************
      module m_vextlow
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/01/09
!     Modification: 2004/04/15, 2004/07/01, 2004/08/20, 2004/09/10,
!                   2005/08/05, 2006/09/21, 2007/01/31, 2007/03/23,
!                   2007/06/27, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2009/11/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     extrapolate the base state variables in lowest layer.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: vextlow, s_vextlow

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vextlow

        module procedure s_vextlow

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic log
      intrinsic min

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_vextlow(etrvar,z1d,u1d,v1d,p1d,pt1d,qv1d)
!***********************************************************************

! Input variable

      character(len=6), intent(in) :: etrvar
                       ! Control flag of extrapolating method

! Input and output variables

      real, intent(inout) :: z1d(0:2)
                       ! Horizontally averaged z coordinates

      real, intent(inout) :: u1d(0:2)
                       ! Horizontally averaged x components of velocity

      real, intent(inout) :: v1d(0:2)
                       ! Horizontally averaged y components of velocity

      real, intent(inout) :: p1d(0:2)
                       ! Horizontally averaged pressure

      real, intent(inout) :: pt1d(0:2)
                       ! Horizontally averaged potential temperature

      real, intent(inout) :: qv1d(0:2)
                       ! Horizontally averaged water vapor mixing ratio

! Internal shared variables

      real gdvcp2      ! 2.0 x g / cp

      real cpdvrd      ! cp / rd
      real rddvcp      ! rd / cp

      real p0iv        ! 1.0 / p0

      real pi0         ! Exner function in lowest layer

      real ptv0        ! Virtual potential temperature in lowest layer
      real ptv1        ! Virtual potential temperature in second layer

!-----7--------------------------------------------------------------7--

! Set the physical parameters.

      gdvcp2=2.e0*g/cp

      cpdvrd=cp/rd
      rddvcp=rd/cp

      p0iv=1.e0/p0

! -----

! Perform extrapolation in the lowest layer.

      z1d(0)=min(2.e0*z1d(1)-z1d(2),0.e0)

      if(etrvar(1:1).eq.'o') then

        u1d(0)=((z1d(2)-z1d(0))*u1d(1)-(z1d(1)-z1d(0))*u1d(2))          &
     &    /(z1d(2)-z1d(1))

      else

        u1d(0)=u1d(1)

      end if

      if(etrvar(2:2).eq.'o') then

        v1d(0)=((z1d(2)-z1d(0))*v1d(1)-(z1d(1)-z1d(0))*v1d(2))          &
     &    /(z1d(2)-z1d(1))

      else

        v1d(0)=v1d(1)

      end if

      if(etrvar(5:5).eq.'o') then

        pt1d(0)=((z1d(2)-z1d(0))*pt1d(1)-(z1d(1)-z1d(0))*pt1d(2))       &
     &    /(z1d(2)-z1d(1))

      else

        pt1d(0)=pt1d(1)

      end if

      if(etrvar(6:6).eq.'o') then

        qv1d(0)=((z1d(2)-z1d(0))*qv1d(1)-(z1d(1)-z1d(0))*qv1d(2))       &
     &    /(z1d(2)-z1d(1))

      else

        qv1d(0)=qv1d(1)

      end if

      ptv0=pt1d(0)*(1.e0+qv1d(0))/(1.e0+epsav*qv1d(0))
      ptv1=pt1d(1)*(1.e0+qv1d(1))/(1.e0+epsav*qv1d(1))

      pi0=exp(rddvcp*log(p0iv*p1d(1)))                                  &
     &  +gdvcp2*(z1d(1)-z1d(0))/(ptv0+ptv1)

      p1d(0)=p0*exp(cpdvrd*log(pi0))

! -----

      end subroutine s_vextlow

!-----7--------------------------------------------------------------7--

      end module m_vextlow
