!***********************************************************************
      module m_vintbase
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/07/01
!     Modification: 2004/08/01, 2004/08/20, 2004/09/10, 2006/01/10,
!                   2006/02/03, 2006/09/21, 2007/03/23, 2008/05/02,
!                   2008/07/01, 2008/08/25, 2008/12/11, 2009/02/27,
!                   2011/12/17, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the interior procedure for vertical interpolating the
!     GPV data to the flat plane to create the base state variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chksat
      use m_comindx
      use m_copy2d
      use m_copy3d
      use m_getcname
      use m_getiname
      use m_getkref
      use m_inichar
      use m_vint31g
      use m_vint31uv

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: vintbase, s_vintbase

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vintbase

        module procedure s_vintbase

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
      subroutine s_vintbase(fpgpvvar,fprefsfc_gpv,                      &
     &                      idstr,idend,jdstr,jdend,kref,               &
     &                      nid,njd,nkd,zdat,zlow,nk,z,km,pdat,         &
     &                      ptdat,udat,vdat,qvdat,pbdat,ptbdat,dtmp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpgpvvar
                       ! Formal parameter of unique index of gpvvar

      integer, intent(in) :: fprefsfc_gpv
                       ! Formal parameter of unique index of refsfc_gpv

      integer, intent(in) :: idstr
                       ! Minimum index of model grid in data region
                       ! in x direction

      integer, intent(in) :: idend
                       ! Maximum index of model grid in data region
                       ! in x direction

      integer, intent(in) :: jdstr
                       ! Minimum index of model grid in data region
                       ! in y direction

      integer, intent(in) :: jdend
                       ! Maximum index of model grid in data region
                       ! in y direction

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

      real, intent(in) :: pdat(1:nid,1:njd,1:km)
                       ! Pressure in data

      real, intent(in) :: ptdat(1:nid,1:njd,1:km)
                       ! Potential temperature in data

! Input and output variables

      real, intent(inout) :: udat(1:nid,1:njd,1:km)
                       ! x components of velocity in data

      real, intent(inout) :: vdat(1:nid,1:njd,1:km)
                       ! y components of velocity in data

      real, intent(inout) :: qvdat(1:nid,1:njd,1:km)
                       ! Water vapor mixing ratio in data

! Output variables

      integer, intent(out) :: kref
                       ! Reference index

      real, intent(out) :: pbdat(1:nid,1:njd,1:km)
                       ! Base state pressure in data

      real, intent(out) :: ptbdat(1:nid,1:njd,1:km)
                       ! Base state potential temperature in data

! Internal shared variables

      character(len=108) gpvvar
                       ! Control flag of input GPV data variables

      integer refsfc_gpv
                       ! Option for
                       ! surface data reference in interpolating

      real, intent(inout) :: zlow(1:nid,1:njd)
                       ! Lowest z physical coordinates in data

      real, intent(inout) :: dtmp1(1:nid,1:njd,1:km)
                       ! Temporary array

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(gpvvar)

! -----

! Get the required namelist variables.

      call getcname(fpgpvvar,gpvvar)
      call getiname(fprefsfc_gpv,refsfc_gpv)

! -----

! Get the terrain height.

      if(refsfc_gpv.eq.1) then

        call s_copy2d(1,nid,1,njd,zdat(1,1,1),zlow)

      end if

! -----

! Get the index of the base state reference pressure.

      call s_getkref(idrefsfc_gpv,idstr,idend,jdstr,jdend,kref,         &
     &               nid,njd,nkd,zlow,zdat,dtmp1,nk,z)

! -----

! Interpolate the x components of velocity.

      call vint31uv(idetrvar_gpv,idrefsfc_gpv,1,nid,njd,nkd,zlow,zdat,  &
     &              udat,nk,z,dtmp1)

      call copy3d(1,nid,1,njd,1,nk,dtmp1,udat)

! -----

! Interpolate the y components of velocity.

      call vint31uv(idetrvar_gpv,idrefsfc_gpv,2,nid,njd,nkd,zlow,zdat,  &
     &              vdat,nk,z,dtmp1)

      call copy3d(1,nid,1,njd,1,nk,dtmp1,vdat)

! -----

! Interpolate the pressure.

      call vint31g(idetrvar_gpv,idrefsfc_gpv,4,nid,njd,nkd,zlow,zdat,   &
     &             pdat,nk,z,pbdat)

! -----

! Interpolate the potential temperature.

      call vint31g(idetrvar_gpv,idrefsfc_gpv,5,nid,njd,nkd,zlow,zdat,   &
     &             ptdat,nk,z,ptbdat)

! -----

! Interpolate the water vapor mixing ratio.

      if(gpvvar(2:2).eq.'o') then

        call vint31g(idetrvar_gpv,idrefsfc_gpv,6,nid,njd,nkd,zlow,zdat, &
     &               qvdat,nk,z,dtmp1)

        call copy3d(1,nid,1,njd,1,nk,dtmp1,qvdat)

      end if

! -----

! Check and avoid the super saturation mixing ratio.

      if(gpvvar(2:2).eq.'o') then

        call chksat('bar  ','ooo',1,nid,1,njd,1,nk,pbdat,ptbdat,        &
     &              pdat,ptdat,qvdat)

      end if

! -----

      end subroutine s_vintbase

!-----7--------------------------------------------------------------7--

      end module m_vintbase
