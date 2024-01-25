!***********************************************************************
      module m_intrpasl
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/09/22
!     Modification:

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the interior procedures for interpolating the aerosol data
!     to the model grid points.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_getrij
      use m_hint3d
      use m_var8w8s
      use m_vint133a

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: intrpasl, s_intrpasl

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface intrpasl

        module procedure s_intrpasl

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
      subroutine s_intrpasl(x0,y0,cpj,x0asl,y0asl,cpjasl,ni,nj,nk,nqa,  &
     &                      z,zph,qasl,ri,rj,di,dj,zph8s,tmp1,          &
     &                      nid,njd,qadat)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqa(0:4)
                       ! Number of types of aerosol

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

      real, intent(in) :: x0asl
                       ! x origin of aerosol data grid

      real, intent(in) :: y0asl
                       ! y origin of aerosol data grid

      real, intent(in) :: cpjasl(1:7)
                       ! Map projection parameters of aerosol data grid

      real, intent(in) :: z(1:nk)
                       ! zeta coordinates

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(in) :: qadat(1:nid,1:njd,1:nk,1:nqa(0))
                       ! Reference aerosol mixing ratio

! Output variable

      real, intent(out) :: qasl(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio

! Internal shared variables

      integer nk_sub   ! Substitute for nk

      integer n        ! Array index in 4th direction

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

      real, intent(inout) :: zph8s(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates at slcar points

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Remark

!     di,dj: These variables are also temporary.

!-----7--------------------------------------------------------------7--

! Set the substituted variable.

      nk_sub=nk

! -----

! Calculate the z physical coordinates at scalar points.

      call var8w8s(ni,nj,nk,zph,zph8s)

! -----

! Calculate the real indices at the model grid points in the data
! region.

      call s_getrij(idmpopt_asl,idnspol_asl,idtlon_asl,                 &
     &              iddxiv_asl,iddyiv_asl,'asldata ',7,'xx',            &
     &              x0,y0,cpj,x0asl,y0asl,cpjasl,ni,nj,ri,rj,di,dj,     &
     &              tmp1(0,0,1),tmp1(0,0,2))

! -----

!! Interpolate the all categories of aerosol.

      do n=1,nqa(0)

! Interpolate horizontally.

        if(n.eq.1) then

          call s_hint3d(idmpopt_asl,idintopt_asl,'xx','cal ',ni,nj,nk,  &
     &                  ri,rj,di,dj,tmp1,nid,njd,qadat(1,1,1,n))

        else

          call s_hint3d(idmpopt_asl,idintopt_asl,'xx','skip',ni,nj,nk,  &
     &                  ri,rj,di,dj,tmp1,nid,njd,qadat(1,1,1,n))

        end if

! -----

! Interpolate vertically.

        call s_vint133a(ni,nj,nk,zph8s,qasl(0,0,1,n),nk_sub,z,tmp1)

! -----

      end do

!! -----

      end subroutine s_intrpasl

!-----7--------------------------------------------------------------7--

      end module m_intrpasl
