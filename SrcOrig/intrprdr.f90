!***********************************************************************
      module m_intrprdr
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/09/09
!     Modification: 2003/04/30, 2003/05/19, 2003/12/12, 2004/08/20,
!                   2004/09/01, 2005/02/10, 2005/04/04, 2005/08/05,
!                   2006/09/21, 2006/09/30, 2006/11/06, 2007/01/20,
!                   2007/03/23, 2007/05/07, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the interior procedures for interpolating the radar data
!     to the model grid points.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_getcname
      use m_getrij
      use m_hintrdr
      use m_inichar
      use m_var8w8s
      use m_var8w8u
      use m_var8w8v
      use m_vint133r

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: intrprdr, s_intrprdr

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface intrprdr

        module procedure s_intrprdr

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
      subroutine s_intrprdr(fprdrvar,x0,y0,cpj,x0rdr,y0rdr,cpjrdr,      &
     &                      ni,nj,nk,z,zph,u,v,w,qp,ri,rj,di,dj,z8uvs,  &
     &                      tmp1,nid,njd,udat,vdat,wdat,qpdat)
!***********************************************************************

! Input variables

      integer, intent(in) :: fprdrvar
                       ! Formal parameter of unique index of rdrvar

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

      real, intent(in) :: x0rdr
                       ! x origin of radar data grid

      real, intent(in) :: y0rdr
                       ! y origin of radar data grid

      real, intent(in) :: cpjrdr(1:7)
                       ! Map projection parameters of radar data grid

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

      real, intent(in) :: qpdat(1:nid,1:njd,1:nk)
                       ! Reference precipitation mixing ratio

! Output variables

      real, intent(out) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(out) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(out) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

      real, intent(out) :: qp(0:ni+1,0:nj+1,1:nk)
                       ! Precipitation mixing ratio

! Internal shared variables

      character(len=108) rdrvar
                       ! Control flag of input radar data variables

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

      call inichar(rdrvar)

! -----

! Get the required namelist variable.

      call getcname(fprdrvar,rdrvar)

! -----

! Set the substituted variable.

      nk_sub=nk

! -----

!! Interpolate the x components of velocity.

      if(rdrvar(1:1).eq.'o') then

! Calculate the z physical coordinates at u points.

        call var8w8u(ni,nj,nk,zph,z8uvs)

! -----

! Calculate the real indices at the model grid points in the data
! region.

        call s_getrij(idmpopt_rdr,idnspol_rdr,idtlon_rdr,               &
     &                iddxiv_rdr,iddyiv_rdr,'radata  ',6,'ox',          &
     &                x0,y0,cpj,x0rdr,y0rdr,cpjrdr,ni,nj,ri,rj,di,dj,   &
     &                tmp1(0,0,1),tmp1(0,0,2))

! -----

! Interpolate horizontally.

        call hintrdr(idmpopt_rdr,'ox','cal ',ni,nj,nk,ri,rj,di,dj,tmp1, &
     &               nid,njd,udat)

! -----

! Interpolate vertically.

        call vint133r('oxx',ni,nj,nk,z8uvs,u,nk_sub,z,tmp1)

! -----

      end if

!! -----

!! Interpolate the y components of velocity to the model grid points.

      if(rdrvar(2:2).eq.'o') then

! Calculate the z physical coordinates at v points.

        call var8w8v(ni,nj,nk,zph,z8uvs)

! -----

! Calculate the real indices at the model grid points in the data
! region.

        call s_getrij(idmpopt_rdr,idnspol_rdr,idtlon_rdr,               &
     &                iddxiv_rdr,iddyiv_rdr,'radata  ',6,'xo',          &
     &                x0,y0,cpj,x0rdr,y0rdr,cpjrdr,ni,nj,ri,rj,di,dj,   &
     &                tmp1(0,0,1),tmp1(0,0,2))

! -----

! Interpolate horizontally.

        call hintrdr(idmpopt_rdr,'xo','cal ',ni,nj,nk,ri,rj,di,dj,tmp1, &
     &               nid,njd,vdat)

! -----

! Interpolate vertically.

        call vint133r('xox',ni,nj,nk,z8uvs,v,nk_sub,z,tmp1)

! -----

      end if

!! -----

! Calculate the real indices at the model grid points in the data
! region.

      if(rdrvar(3:3).eq.'o'.or.rdrvar(4:4).eq.'o') then

        call s_getrij(idmpopt_rdr,idnspol_rdr,idtlon_rdr,               &
     &                iddxiv_rdr,iddyiv_rdr,'radata  ',6,'xx',          &
     &                x0,y0,cpj,x0rdr,y0rdr,cpjrdr,ni,nj,ri,rj,di,dj,   &
     &                tmp1(0,0,1),tmp1(0,0,2))

      end if

! -----

!! Interpolate the z components of velocity.

      if(rdrvar(3:3).eq.'o') then

! Interpolate horizontally.

        call hintrdr(idmpopt_rdr,'xx','cal ',ni,nj,nk,ri,rj,di,dj,tmp1, &
     &               nid,njd,wdat)

! -----

! Interpolate vertically.

        call vint133r('xxo',ni,nj,nk,zph,w,nk_sub,z,tmp1)

! -----

      end if

!! -----

!! Interpolate the precipitation mixing ratio.

      if(rdrvar(4:4).eq.'o') then

! Calculate the z physical coordinates at scalar points.

        call var8w8s(ni,nj,nk,zph,z8uvs)

! -----

! Interpolate horizontally.

        if(rdrvar(3:3).eq.'o') then

         call hintrdr(idmpopt_rdr,'xx','skip',ni,nj,nk,ri,rj,di,dj,tmp1,&
     &                nid,njd,qpdat)

        else

         call hintrdr(idmpopt_rdr,'xx','cal ',ni,nj,nk,ri,rj,di,dj,tmp1,&
     &                nid,njd,qpdat)

        end if

! -----

! Interpolate vertically.

        call vint133r('xxx',ni,nj,nk,z8uvs,qp,nk_sub,z,tmp1)

! -----

      end if

!! -----

      end subroutine s_intrprdr

!-----7--------------------------------------------------------------7--

      end module m_intrprdr
