!***********************************************************************
      module m_getbase
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/06/29
!     Modification: 2002/04/02, 2002/06/18, 2002/07/15, 2002/08/15,
!                   2002/09/09, 2003/04/30, 2003/05/19, 2003/12/12,
!                   2004/01/09, 2004/07/01, 2005/04/04, 2006/09/21,
!                   2007/05/07, 2008/05/02, 2008/07/01, 2008/08/25,
!                   2009/02/27, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     extract the base state variables for the model grid from the
!     interpolated GPV variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getcname
      use m_inichar
      use m_var8w8s
      use m_var8w8u
      use m_var8w8v
      use m_vint13

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getbase, s_getbase

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getbase

        module procedure s_getbase

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
      subroutine s_getbase(fpgpvvar,ni,nj,nk,zph,z1d,u1d,v1d,           &
     &                     p1d,pt1d,qv1d,ubr,vbr,pbr,ptbr,qvbr,z8uvs)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpgpvvar
                       ! Formal parameter of unique index of gpvvar

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(in) :: z1d(0:4*nk-3)
                       ! Horizontally averaged z coordinates

      real, intent(in) :: u1d(0:4*nk-3)
                       ! Horizontally averaged x components of velocity

      real, intent(in) :: v1d(0:4*nk-3)
                       ! Horizontally averaged y components of velocity

      real, intent(in) :: p1d(0:4*nk-3)
                       ! Horizontally averaged pressure

      real, intent(in) :: pt1d(0:4*nk-3)
                       ! Horizontally averaged potential temperature

      real, intent(in) :: qv1d(0:4*nk-3)
                       ! Horizontally averaged water vapor mixing ratio

! Output variables

      real, intent(out) :: ubr(0:ni+1,0:nj+1,1:nk)
                       ! Base state x components of velocity

      real, intent(out) :: vbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state y components of velocity

      real, intent(out) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(out) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(out) :: qvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state water vapor mixing ratio

! Internal shared variables

      character(len=108) gpvvar
                       ! Control flag of input GPV data variables

      real, intent(inout) :: z8uvs(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates at u, v or slcar points

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(gpvvar)

! -----

! Get the required namelist variable.

      call getcname(fpgpvvar,gpvvar)

! -----

! Calculate the z physical coordinates at u points.

      call var8w8u(ni,nj,nk,zph,z8uvs)

! -----

! Get the base state x components of velocity.

      call s_vint13('ox',0,ni+1,0,nj+1,2,nk-2,z8uvs(0,0,2),ubr(0,0,2),  &
     &              4*nk-3,z1d,u1d)

! -----

! Calculate the z physical coordinates at v points.

      call var8w8v(ni,nj,nk,zph,z8uvs)

! -----

! Get the base state y components of velocity.

      call s_vint13('xo',0,ni+1,0,nj+1,2,nk-2,z8uvs(0,0,2),vbr(0,0,2),  &
     &              4*nk-3,z1d,v1d)

! -----

! Calculate the z physical coordinates at scalar points.

      call var8w8s(ni,nj,nk,zph,z8uvs)

! -----

! Get the base state pressure.

      call s_vint13('xx',0,ni+1,0,nj+1,2,nk-2,z8uvs(0,0,2),pbr(0,0,2),  &
     &              4*nk-3,z1d,p1d)

! -----

! Get the base state potential temperature.

      call s_vint13('xx',0,ni+1,0,nj+1,2,nk-2,z8uvs(0,0,2),ptbr(0,0,2), &
     &              4*nk-3,z1d,pt1d)

! -----

! Get the base state water vapor mixing ratio.

      if(gpvvar(2:2).eq.'o') then

       call s_vint13('xx',0,ni+1,0,nj+1,2,nk-2,z8uvs(0,0,2),qvbr(0,0,2),&
     &               4*nk-3,z1d,qv1d)

      end if

! -----

      end subroutine s_getbase

!-----7--------------------------------------------------------------7--

      end module m_getbase
