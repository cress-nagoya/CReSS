!***********************************************************************
      module m_vintgpv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/07/03
!     Modification: 2002/07/15, 2002/09/02, 2003/04/30, 2003/05/19,
!                   2003/12/12, 2004/01/09, 2004/03/05, 2004/05/07,
!                   2004/08/20, 2006/01/10, 2006/02/03, 2006/09/21,
!                   2007/03/23, 2007/07/30, 2008/05/02, 2008/07/01,
!                   2008/08/25, 2008/12/11, 2009/02/27, 2011/09/22,
!                   2011/12/17, 2013/01/28, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the interior procedures for vertical interpolating the
!     GPV data to the flat plane.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chksat
      use m_comindx
      use m_comkind
      use m_copy2d
      use m_copy3d
      use m_getcname
      use m_getiname
      use m_inichar
      use m_vint13
      use m_vint31g
      use m_vint31uv

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: vintgpv, s_vintgpv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vintgpv

        module procedure s_vintgpv

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
      subroutine s_vintgpv(fpgpvvar,fprefsfc_gpv,it,nid,njd,nkd,zdat,   &
     &                     zlow,nk,z,z1d,p1d,pt1d,varef,pbdat,ptbdat,   &
     &                     km,udat,vdat,wdat,ppdat,ptpdat,qvdat,        &
     &                     qcdat,qrdat,qidat,qsdat,qgdat,qhdat)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpgpvvar
                       ! Formal parameter of unique index of gpvvar

      integer, intent(in) :: fprefsfc_gpv
                       ! Formal parameter of unique index of refsfc_gpv

      integer(kind=i8), intent(in) :: it
                       ! Index of main do loop in upper procedure

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

      real, intent(in) :: z1d(0:4*nk-3)
                       ! Horizontally averaged z coordinates

      real, intent(in) :: p1d(0:4*nk-3)
                       ! Horizontally averaged pressure

      real, intent(in) :: pt1d(0:4*nk-3)
                       ! Horizontally averaged potential temperature

! Input and output variables

      real, intent(inout) :: udat(1:nid,1:njd,1:km)
                       ! x components of velocity in data

      real, intent(inout) :: vdat(1:nid,1:njd,1:km)
                       ! y components of velocity in data

      real, intent(inout) :: wdat(1:nid,1:njd,1:km)
                       ! z components of velocity in data

      real, intent(inout) :: ppdat(1:nid,1:njd,1:km)
                       ! Pressure perturbation in data

      real, intent(inout) :: ptpdat(1:nid,1:njd,1:km)
                       ! Potential temperature perturbation in data

      real, intent(inout) :: qvdat(1:nid,1:njd,1:km)
                       ! Water vapor mixing ratio in data

      real, intent(inout) :: qcdat(1:nid,1:njd,1:km)
                       ! Cloud water mixing ratio in data

      real, intent(inout) :: qrdat(1:nid,1:njd,1:km)
                       ! Rain water mixing ratio in data

      real, intent(inout) :: qidat(1:nid,1:njd,1:km)
                       ! Cloud ice mixing ratio in data

      real, intent(inout) :: qsdat(1:nid,1:njd,1:km)
                       ! Snow mixing ratio in data

      real, intent(inout) :: qgdat(1:nid,1:njd,1:km)
                       ! Graupel mixing ratio in data

      real, intent(inout) :: qhdat(1:nid,1:njd,1:km)
                       ! Hail mixing ratio in data

! Internal shared variables

      character(len=108) gpvvar
                       ! Control flag of input GPV data variables

      integer refsfc_gpv
                       ! Option for
                       ! surface data reference in interpolating

      real, intent(inout) :: zlow(1:nid,1:njd)
                       ! Lowest z physical coordinates in data

      real, intent(inout) :: varef(1:nid,1:njd,1:nk)
                       ! Optional vertically interpolated variable

      real, intent(inout) :: pbdat(1:nid,1:njd,1:nk)
                       ! Base state pressure in data

      real, intent(inout) :: ptbdat(1:nid,1:njd,1:nk)
                       ! Base state potential temperature in data

! Internal private variables

      integer id       ! Array index in x direction
      integer jd       ! Array index in y direction

      integer k        ! Array index in z direction

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

! Interpolate the x components of velocity.

      if(it.ne.1_i8) then

        call vint31uv(idetrvar_gpv,idrefsfc_gpv,1,nid,njd,nkd,zlow,zdat,&
     &                udat,nk,z,varef)

        call copy3d(1,nid,1,njd,1,nk,varef,udat)

      end if

! -----

! Interpolate the y components of velocity.

      if(it.ne.1_i8) then

        call vint31uv(idetrvar_gpv,idrefsfc_gpv,2,nid,njd,nkd,zlow,zdat,&
     &                vdat,nk,z,varef)

        call copy3d(1,nid,1,njd,1,nk,varef,vdat)

      end if

! -----

! Interpolate the pressure.

      call vint31g(idetrvar_gpv,idrefsfc_gpv,4,nid,njd,nkd,zlow,zdat,   &
     &             ppdat,nk,z,varef)

      call copy3d(1,nid,1,njd,1,nk,varef,ppdat)

! -----

! Interpolate the potential temperature.

      call vint31g(idetrvar_gpv,idrefsfc_gpv,5,nid,njd,nkd,zlow,zdat,   &
     &             ptpdat,nk,z,varef)

      call copy3d(1,nid,1,njd,1,nk,varef,ptpdat)

! -----

! Interpolate the z components of velocity.

      if(gpvvar(1:1).eq.'o') then

        call vint31g(idetrvar_gpv,idrefsfc_gpv,3,nid,njd,nkd,zlow,zdat, &
     &               wdat,nk,z,varef)

        call copy3d(1,nid,1,njd,1,nk,varef,wdat)

      end if

! -----

! Interpolate the water vapor mixing ratio.

      if(gpvvar(2:2).eq.'o'.and.it.ne.1_i8) then

        call vint31g(idetrvar_gpv,idrefsfc_gpv,6,nid,njd,nkd,zlow,zdat, &
     &               qvdat,nk,z,varef)

        call copy3d(1,nid,1,njd,1,nk,varef,qvdat)

      end if

! -----

! Interpolate the cloud water mixing ratio.

      if(gpvvar(3:3).eq.'o') then

        call vint31g(idetrvar_gpv,idrefsfc_gpv,7,nid,njd,nkd,zlow,zdat, &
     &               qcdat,nk,z,varef)

        call copy3d(1,nid,1,njd,1,nk,varef,qcdat)

      end if

! -----

! Interpolate the rain water mixing ratio.

      if(gpvvar(4:4).eq.'o') then

        call vint31g(idetrvar_gpv,idrefsfc_gpv,7,nid,njd,nkd,zlow,zdat, &
     &               qrdat,nk,z,varef)

        call copy3d(1,nid,1,njd,1,nk,varef,qrdat)

      end if

! -----

! Interpolate the cloud ice mixing ratio.

      if(gpvvar(5:5).eq.'o') then

        call vint31g(idetrvar_gpv,idrefsfc_gpv,7,nid,njd,nkd,zlow,zdat, &
     &               qidat,nk,z,varef)

        call copy3d(1,nid,1,njd,1,nk,varef,qidat)

      end if

! -----

! Interpolate the snow mixing ratio.

      if(gpvvar(6:6).eq.'o') then

        call vint31g(idetrvar_gpv,idrefsfc_gpv,7,nid,njd,nkd,zlow,zdat, &
     &               qsdat,nk,z,varef)

        call copy3d(1,nid,1,njd,1,nk,varef,qsdat)

      end if

! -----

! Interpolate the graupel mixing ratio.

      if(gpvvar(7:7).eq.'o') then

        call vint31g(idetrvar_gpv,idrefsfc_gpv,7,nid,njd,nkd,zlow,zdat, &
     &               qgdat,nk,z,varef)

        call copy3d(1,nid,1,njd,1,nk,varef,qgdat)

      end if

! -----

! Interpolate the hail mixing ratio.

      if(gpvvar(8:8).eq.'o') then

        call vint31g(idetrvar_gpv,idrefsfc_gpv,7,nid,njd,nkd,zlow,zdat, &
     &               qhdat,nk,z,varef)

        call copy3d(1,nid,1,njd,1,nk,varef,qhdat)

      end if

! -----

!! Check and avoid the super saturation mixing ratio.

! Get the 3 dimensional constant z physical coordinates.

      if(gpvvar(2:2).eq.'o') then

!$omp parallel default(shared) private(k)

        do k=1,nk

!$omp do schedule(runtime) private(id,jd)

          do jd=1,njd
          do id=1,nid
            varef(id,jd,k)=z(k)
          end do
          end do

!$omp end do

        end do

!$omp end parallel

! -----

! Get the base state variables.

        call vint13('oo',1,nid,1,njd,1,nk,varef,pbdat,4*nk-3,z1d,p1d)

        call vint13('oo',1,nid,1,njd,1,nk,varef,ptbdat,4*nk-3,z1d,pt1d)

! -----

! Perform checking.

        call chksat('total','ooo',1,nid,1,njd,1,nk,pbdat,ptbdat,        &
     &              ppdat,ptpdat,qvdat)

! -----

      end if

!! -----

      end subroutine s_vintgpv

!-----7--------------------------------------------------------------7--

      end module m_vintgpv
