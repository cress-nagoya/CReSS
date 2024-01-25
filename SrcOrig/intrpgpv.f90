!***********************************************************************
      module m_intrpgpv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/07/03
!     Modification: 2002/07/15, 2003/04/30, 2003/05/19, 2003/12/12,
!                   2004/01/09, 2004/03/05, 2004/05/07, 2004/08/20,
!                   2004/09/01, 2005/02/10, 2005/04/04, 2006/02/03,
!                   2006/09/21, 2006/09/30, 2006/11/06, 2007/01/20,
!                   2007/03/23, 2007/05/07, 2008/05/02, 2008/07/01,
!                   2008/08/25, 2008/12/11, 2009/02/27, 2011/09/22,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the interior procedures for interpolating the GPV data to
!     the model grid points.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_getcname
      use m_getrij
      use m_hint3d
      use m_inichar
      use m_var8w8s
      use m_var8w8u
      use m_var8w8v
      use m_vint133g

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: intrpgpv, s_intrpgpv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface intrpgpv

        module procedure s_intrpgpv

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
      subroutine s_intrpgpv(fpgpvvar,x0,y0,cpj,x0gpv,y0gpv,cpjgpv,      &
     &                      ni,nj,nk,z,zph,u,v,w,pp,ptp,qv,qc,qr,       &
     &                      qi,qs,qg,qh,ri,rj,di,dj,z8uvs,tmp1,         &
     &                      nid,njd,udat,vdat,wdat,ppdat,ptpdat,        &
     &                      qvdat,qcdat,qrdat,qidat,qsdat,qgdat,qhdat)
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

      integer, intent(in) :: nid
                       ! Data dimension in x direction

      integer, intent(in) :: njd
                       ! Data dimension in y direction

      real, intent(in) :: x0
                       ! x origin of model grid

      real, intent(in) :: y0
                       ! y origin of model grid

      real, intent(in) :: cpj(1:7)
                       ! Map projection parameters of model grid

      real, intent(in) :: x0gpv
                       ! x origin of GPV data grid

      real, intent(in) :: y0gpv
                       ! y origin of GPV data grid

      real, intent(in) :: cpjgpv(1:7)
                       ! Map projection parameters of GPV data grid

      real, intent(in) :: z(1:nk)
                       ! zeta coordinates

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(in) :: udat(1:nid,1:njd,1:nk)
                       ! Reference x components of velocity

      real, intent(in) :: vdat(1:nid,1:njd,1:nk)
                       ! Reference y components of velocity

      real, intent(in) :: wdat(1:nid,1:njd,1:nk)
                       ! Reference z components of velocity

      real, intent(in) :: ppdat(1:nid,1:njd,1:nk)
                       ! Reference pressure perturbation

      real, intent(in) :: ptpdat(1:nid,1:njd,1:nk)
                       ! Reference potential temperature perturbation

      real, intent(in) :: qvdat(1:nid,1:njd,1:nk)
                       ! Reference water vapor mixing ratio

      real, intent(in) :: qcdat(1:nid,1:njd,1:nk)
                       ! Reference cloud water mixing ratio

      real, intent(in) :: qrdat(1:nid,1:njd,1:nk)
                       ! Reference rain water mixing ratio

      real, intent(in) :: qidat(1:nid,1:njd,1:nk)
                       ! Reference cloud ice mixing ratio

      real, intent(in) :: qsdat(1:nid,1:njd,1:nk)
                       ! Reference snow mixing ratio

      real, intent(in) :: qgdat(1:nid,1:njd,1:nk)
                       ! Reference graupel mixing ratio

      real, intent(in) :: qhdat(1:nid,1:njd,1:nk)
                       ! Reference hail mixing ratio

! Output variables

      real, intent(out) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(out) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(out) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

      real, intent(out) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation

      real, intent(out) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

      real, intent(out) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

      real, intent(out) :: qc(0:ni+1,0:nj+1,1:nk)
                       ! Cloud water mixing ratio

      real, intent(out) :: qr(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixing ratio

      real, intent(out) :: qi(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixing ratio

      real, intent(out) :: qs(0:ni+1,0:nj+1,1:nk)
                       ! Snow mixing ratio

      real, intent(out) :: qg(0:ni+1,0:nj+1,1:nk)
                       ! Graupel mixing ratio

      real, intent(out) :: qh(0:ni+1,0:nj+1,1:nk)
                       ! Hail mixing ratio

! Internal shared variables

      character(len=108) gpvvar
                       ! Control flag of input GPV data variables

      integer nk_sub   ! Substitute for nk

      real, intent(inout) :: ri(0:ni+1,0:nj+1)
                       ! Real indices in data region in x direction

      real, intent(inout) :: rj(0:ni+1,0:nj+1)
                       ! Real indices in data region in y direction

      real, intent(inout) :: di(0:ni+1,0:nj+1)
                       ! Distance between data and model grid points
                       ! in x direction

      real, intent(inout) :: dj(0:ni+1,0:nj+1)
                       ! Distance between data and model grid points
                       ! in y direction

      real, intent(inout) :: z8uvs(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates at u, v or slcar points

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Remark

!     di,dj: These variables are also temporary.

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(gpvvar)

! -----

! Get the required namelist variable.

      call getcname(fpgpvvar,gpvvar)

! -----

! Set the substituted variable.

      nk_sub=nk

! -----

!! Interpolate the x components of velocity.

! Calculate the z physical coordinates at u points.

      call var8w8u(ni,nj,nk,zph,z8uvs)

! -----

! Calculate the real indices at the model grid points in the data
! region.

      call s_getrij(idmpopt_gpv,idnspol_gpv,idtlon_gpv,                 &
     &              iddxiv_gpv,iddyiv_gpv,'gridata ',7,'ox',            &
     &              x0,y0,cpj,x0gpv,y0gpv,cpjgpv,ni,nj,ri,rj,di,dj,     &
     &              tmp1(0,0,1),tmp1(0,0,2))

! -----

! Interpolate horizontally.

      call hint3d(idmpopt_gpv,idintopt_gpv,'ox','cal ',ni,nj,nk,ri,rj,  &
     &            di,dj,tmp1,nid,njd,udat)

! -----

! Interpolate vertically.

      call vint133g('oxx',ni,nj,nk,z8uvs,u,nk_sub,z,tmp1)

! -----

!! -----

!! Interpolate the y components of velocity to the model grid points.

! Calculate the z physical coordinates at v points.

      call var8w8v(ni,nj,nk,zph,z8uvs)

! -----

! Calculate the real indices at the model grid points in the data
! region.

      call s_getrij(idmpopt_gpv,idnspol_gpv,idtlon_gpv,                 &
     &              iddxiv_gpv,iddyiv_gpv,'gridata ',7,'xo',            &
     &              x0,y0,cpj,x0gpv,y0gpv,cpjgpv,ni,nj,ri,rj,di,dj,     &
     &              tmp1(0,0,1),tmp1(0,0,2))

! -----

! Interpolate horizontally.

      call hint3d(idmpopt_gpv,idintopt_gpv,'xo','cal ',ni,nj,nk,ri,rj,  &
     &            di,dj,tmp1,nid,njd,vdat)

! -----

! Interpolate vertically.

      call vint133g('xox',ni,nj,nk,z8uvs,v,nk_sub,z,tmp1)

! -----

!! -----

! Calculate the z physical coordinates at scalar points.

      call var8w8s(ni,nj,nk,zph,z8uvs)

! -----

! Calculate the real indices at the model grid points in the data
! region.

      call s_getrij(idmpopt_gpv,idnspol_gpv,idtlon_gpv,                 &
     &              iddxiv_gpv,iddyiv_gpv,'gridata ',7,'xx',            &
     &              x0,y0,cpj,x0gpv,y0gpv,cpjgpv,ni,nj,ri,rj,di,dj,     &
     &              tmp1(0,0,1),tmp1(0,0,2))

! -----

!! Interpolate the pressure.

! Interpolate horizontally.

      call hint3d(idmpopt_gpv,idintopt_gpv,'xx','cal ',ni,nj,nk,ri,rj,  &
     &            di,dj,tmp1,nid,njd,ppdat)

! -----

! Interpolate vertically.

      call vint133g('xxx',ni,nj,nk,z8uvs,pp,nk_sub,z,tmp1)

! -----

!! -----

!! Interpolate the potential temperature.

! Interpolate horizontally.

      call hint3d(idmpopt_gpv,idintopt_gpv,'xx','skip',ni,nj,nk,ri,rj,  &
     &            di,dj,tmp1,nid,njd,ptpdat)

! -----

! Interpolate vertically.

      call vint133g('xxx',ni,nj,nk,z8uvs,ptp,nk_sub,z,tmp1)

! -----

!! -----

!! Interpolate the z components of velocity.

      if(gpvvar(1:1).eq.'o') then

! Interpolate horizontally.

        call hint3d(idmpopt_gpv,idintopt_gpv,'xx','skip',ni,nj,nk,ri,rj,&
     &              di,dj,tmp1,nid,njd,wdat)

! -----

! Interpolate vertically.

        call vint133g('xxo',ni,nj,nk,zph,w,nk_sub,z,tmp1)

! -----

      end if

!! -----

!! Interpolate the water vapor mixing ratio.

      if(gpvvar(2:2).eq.'o') then

! Interpolate horizontally.

        call hint3d(idmpopt_gpv,idintopt_gpv,'xx','skip',ni,nj,nk,ri,rj,&
     &              di,dj,tmp1,nid,njd,qvdat)

! -----

! Interpolate vertically.

        call vint133g('xxx',ni,nj,nk,z8uvs,qv,nk_sub,z,tmp1)

! -----

      end if

!! -----

!! Interpolate the cloud water mixing ratio.

      if(gpvvar(3:3).eq.'o') then

! Interpolate horizontally.

        call hint3d(idmpopt_gpv,idintopt_gpv,'xx','skip',ni,nj,nk,ri,rj,&
     &              di,dj,tmp1,nid,njd,qcdat)

! -----

! Interpolate vertically.

        call vint133g('xxx',ni,nj,nk,z8uvs,qc,nk_sub,z,tmp1)

! -----

      end if

!! -----

!! Interpolate the rain water mixing ratio.

      if(gpvvar(4:4).eq.'o') then

! Interpolate horizontally.

        call hint3d(idmpopt_gpv,idintopt_gpv,'xx','skip',ni,nj,nk,ri,rj,&
     &              di,dj,tmp1,nid,njd,qrdat)

! -----

! Interpolate vertically.

        call vint133g('xxx',ni,nj,nk,z8uvs,qr,nk_sub,z,tmp1)

! -----

      end if

!! -----

!! Interpolate the cloud ice mixing ratio.

      if(gpvvar(5:5).eq.'o') then

! Interpolate horizontally.

        call hint3d(idmpopt_gpv,idintopt_gpv,'xx','skip',ni,nj,nk,ri,rj,&
     &              di,dj,tmp1,nid,njd,qidat)

! -----

! Interpolate vertically.

        call vint133g('xxx',ni,nj,nk,z8uvs,qi,nk_sub,z,tmp1)

! -----

      end if

!! -----

!! Interpolate the snow mixing ratio.

      if(gpvvar(6:6).eq.'o') then

! Interpolate horizontally.

        call hint3d(idmpopt_gpv,idintopt_gpv,'xx','skip',ni,nj,nk,ri,rj,&
     &              di,dj,tmp1,nid,njd,qsdat)

! -----

! Interpolate vertically.

        call vint133g('xxx',ni,nj,nk,z8uvs,qs,nk_sub,z,tmp1)

! -----

      end if

!! -----

!! Interpolate the graupel mixing ratio.

      if(gpvvar(7:7).eq.'o') then

! Interpolate horizontally.

        call hint3d(idmpopt_gpv,idintopt_gpv,'xx','skip',ni,nj,nk,ri,rj,&
     &              di,dj,tmp1,nid,njd,qgdat)

! -----

! Interpolate vertically.

        call vint133g('xxx',ni,nj,nk,z8uvs,qg,nk_sub,z,tmp1)

! -----

      end if

!! -----

!! Interpolate the hail mixing ratio.

      if(gpvvar(8:8).eq.'o') then

! Interpolate horizontally.

        call hint3d(idmpopt_gpv,idintopt_gpv,'xx','skip',ni,nj,nk,ri,rj,&
     &              di,dj,tmp1,nid,njd,qhdat)

! -----

! Interpolate vertically.

        call vint133g('xxx',ni,nj,nk,z8uvs,qh,nk_sub,z,tmp1)

! -----

      end if

!! -----

      end subroutine s_intrpgpv

!-----7--------------------------------------------------------------7--

      end module m_intrpgpv
