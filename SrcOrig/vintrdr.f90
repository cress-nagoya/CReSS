!***********************************************************************
      module m_vintrdr
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/09/09
!     Modification: 2003/04/30, 2003/05/19, 2006/01/10, 2006/09/21,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2013/02/13,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the interior procedures for vertical interpolating the
!     radar data to the flat plane.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_copy3d
      use m_getcname
      use m_inichar
      use m_vint31r

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: vintrdr, s_vintrdr

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vintrdr

        module procedure s_vintrdr

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
      subroutine s_vintrdr(fprdrvar,nid,njd,nkd,zdat,nk,z,varef,        &
     &                     km,udat,vdat,wdat,qpdat)
!***********************************************************************

! Input variables

      integer, intent(in) :: fprdrvar
                       ! Formal parameter of unique index of rdrvar

      integer, intent(in) :: nid
                       ! Data dimension in x direction

      integer, intent(in) :: njd
                       ! Data dimension in y direction

      integer, intent(in) :: nkd
                       ! Data dimension in z direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: km
                       ! Maximum dimension in z direction

      real, intent(in) :: zdat(1:nid,1:njd,1:nkd)
                       ! z physical coordinates in data

      real, intent(in) :: z(1:nk)
                       ! zeta coordinates

! Input and output variables

      real, intent(inout) :: udat(1:nid,1:njd,1:km)
                       ! x components of velocity in data

      real, intent(inout) :: vdat(1:nid,1:njd,1:km)
                       ! y components of velocity in data

      real, intent(inout) :: wdat(1:nid,1:njd,1:km)
                       ! z components of velocity in data

      real, intent(inout) :: qpdat(1:nid,1:njd,1:km)
                       ! Precipitarion mixing ratio in data

! Internal shared variables

      character(len=108) rdrvar
                       ! Control flag of input radar data variables

      real, intent(inout) :: varef(1:nid,1:njd,1:nk)
                       ! Optional vertically interpolated variable

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(rdrvar)

! -----

! Get the required namelist variable.

      call getcname(fprdrvar,rdrvar)

! -----

! Interpolate the x components of velocity.

      if(rdrvar(1:1).eq.'o') then

        call vint31r(nid,njd,nkd,zdat,udat,nk,z,varef)

        call copy3d(1,nid,1,njd,1,nk,varef,udat)

      end if

! -----

! Interpolate the y components of velocity.

      if(rdrvar(2:2).eq.'o') then

        call vint31r(nid,njd,nkd,zdat,vdat,nk,z,varef)

        call copy3d(1,nid,1,njd,1,nk,varef,vdat)

      end if

! -----

! Interpolate the z components of velocity.

      if(rdrvar(3:3).eq.'o') then

        call vint31r(nid,njd,nkd,zdat,wdat,nk,z,varef)

        call copy3d(1,nid,1,njd,1,nk,varef,wdat)

      end if

! -----

! Interpolate the precipitation mixing ratio.

      if(rdrvar(4:4).eq.'o') then

        call vint31r(nid,njd,nkd,zdat,qpdat,nk,z,varef)

        call copy3d(1,nid,1,njd,1,nk,varef,qpdat)

      end if

! -----

      end subroutine s_vintrdr

!-----7--------------------------------------------------------------7--

      end module m_vintrdr
