!***********************************************************************
      module m_intrpsnd
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/05/10,
!                   1999/05/20, 1999/06/28, 1999/07/05, 1999/08/03,
!                   1999/08/18, 1999/08/23, 1999/09/30, 1999/10/12,
!                   1999/11/01, 2000/01/17, 2000/03/08, 2000/07/05,
!                   2001/01/15, 2001/04/15, 2001/05/29, 2001/06/29,
!                   2001/10/17, 2002/04/02, 2002/06/06, 2002/06/18,
!                   2002/08/15, 2002/09/09, 2003/04/30, 2003/05/19,
!                   2003/11/05, 2003/12/12, 2004/01/09, 2004/09/10,
!                   2005/02/10, 2005/04/04, 2006/09/21, 2006/09/30,
!                   2006/12/04, 2007/01/05, 2007/01/20, 2007/05/07,
!                   2008/05/02, 2008/07/01, 2008/08/25, 2009/02/27,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     extract the base state variables from interpolated sounding data.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bc8u
      use m_bc8v
      use m_bcyclex
      use m_bcycley
      use m_combuf
      use m_comindx
      use m_getbufgx
      use m_getbufgy
      use m_getbufsx
      use m_getbufsy
      use m_getcname
      use m_inichar
      use m_putbufgx
      use m_putbufgy
      use m_putbufsx
      use m_putbufsy
      use m_shiftgx
      use m_shiftgy
      use m_shiftsx
      use m_shiftsy
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

      public :: intrpsnd, s_intrpsnd

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface intrpsnd

        module procedure s_intrpsnd

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
      subroutine s_intrpsnd(fpsndtyp,ni,nj,nk,zph,ubr,vbr,pbr,ptbr,qvbr,&
     &                      zph8s,zph8u,zph8v,nlev,z1d,u1d,v1d,p1d,     &
     &                      pt1d,qv1d)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpsndtyp
                       ! Formal parameter of unique index of sndtyp

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nlev
                       ! Horizontally averaged vertical dimension

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(in) :: z1d(0:nlev)
                       ! Horizontally averaged z physical coordinates

      real, intent(in) :: u1d(0:nlev)
                       ! Horizontally averaged x components of velocity

      real, intent(in) :: v1d(0:nlev)
                       ! Horizontally averaged y components of velocity

      real, intent(in) :: p1d(0:nlev)
                       ! Horizontally averaged pressure

      real, intent(in) :: pt1d(0:nlev)
                       ! Horizontally averaged potential temrerature

      real, intent(in) :: qv1d(0:nlev)
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

      character(len=108) sndtyp
                       ! Control flag of sounding data type

      integer nlev25   ! nlev / 4

      integer ni_sub   ! Substitute for ni
      integer nj_sub   ! Substitute for nj

      real, intent(inout) :: zph8s(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates at scalar points

      real, intent(inout) :: zph8u(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates at u points

      real, intent(inout) :: zph8v(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates at v points

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(sndtyp)

! -----

! Get the required namelist variable.

      call getcname(fpsndtyp,sndtyp)

! -----

! Set the common used variable.

      nlev25=nlev/4

! -----

! Set the substituted variables.

      ni_sub=ni
      nj_sub=nj

! -----

! Calculate the z physical coordinates at the scalar, u and v points.

      call var8w8s(ni,nj,nk,zph,zph8s)
      call var8w8u(ni,nj,nk,zph,zph8u)
      call var8w8v(ni,nj,nk,zph,zph8v)

! -----

!! Set the lateral boundary conditions.

! Set the west and east bondary conditions.

      call s_putbufsx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,zph8u,1,1,sbuf)

      call s_shiftsx(idwbc,idebc,'bnd',nj,nk,1,sbuf,rbuf)

      call s_getbufsx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,zph8u,1,1,    &
     &                rbuf)

      call s_putbufgx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,zph8u,1,1,sbuf)

      call s_shiftgx(idwbc,idebc,'bnd',nj,nk,1,sbuf,rbuf)

      call s_getbufgx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,zph8u,1,1,    &
     &                rbuf)

      call bcyclex(idwbc,idebc,3,1,ni-2,ni_sub,ni,nj,nk,zph8u)

      call bc8u(idwbc,idebc,ni,nj,nk,zph8u)

! -----

! Set the south and north bondary conditions.

      call s_putbufsy(idsbc,idnbc,'bnd',3,nj-2,ni,nj,nk,zph8v,1,1,sbuf)

      call s_shiftsy(idsbc,idnbc,'bnd',ni,nk,1,sbuf,rbuf)

      call s_getbufsy(idsbc,idnbc,'bnd',1,nj_sub,ni,nj,nk,zph8v,1,1,    &
     &                rbuf)

      call s_putbufgy(idsbc,idnbc,'bnd',3,nj-2,ni,nj,nk,zph8v,1,1,sbuf)

      call s_shiftgy(idsbc,idnbc,'bnd',ni,nk,1,sbuf,rbuf)

      call s_getbufgy(idsbc,idnbc,'bnd',1,nj_sub,ni,nj,nk,zph8v,1,1,    &
     &                rbuf)

      call bcycley(idsbc,idnbc,3,1,nj-2,nj_sub,ni,nj,nk,zph8v)

      call bc8v(idsbc,idnbc,ni,nj,nk,zph8v)

! -----

!! -----

! Extract the base state variables from interpolated sounding data.

      if(sndtyp(1:1).eq.'z') then

       call s_vint13('ox',0,ni+1,0,nj+1,2,nk-2,zph8u(0,0,2),ubr(0,0,2), &
     &               nlev,z1d,u1d)

       call s_vint13('xo',0,ni+1,0,nj+1,2,nk-2,zph8v(0,0,2),vbr(0,0,2), &
     &               nlev,z1d,v1d)

       call s_vint13('xx',0,ni+1,0,nj+1,2,nk-2,zph8s(0,0,2),pbr(0,0,2), &
     &               nlev,z1d,p1d)

       call s_vint13('xx',0,ni+1,0,nj+1,2,nk-2,zph8s(0,0,2),ptbr(0,0,2),&
     &               nlev,z1d,pt1d)

       call s_vint13('xx',0,ni+1,0,nj+1,2,nk-2,zph8s(0,0,2),qvbr(0,0,2),&
     &               nlev,z1d,qv1d)

      else if(sndtyp(1:1).eq.'p') then

       call s_vint13('ox',0,ni+1,0,nj+1,2,nk-2,zph8u(0,0,2),ubr(0,0,2), &
     &               nlev25,z1d,u1d)

       call s_vint13('xo',0,ni+1,0,nj+1,2,nk-2,zph8v(0,0,2),vbr(0,0,2), &
     &               nlev25,z1d,v1d)

       call s_vint13('xx',0,ni+1,0,nj+1,2,nk-2,zph8s(0,0,2),pbr(0,0,2), &
     &               nlev25,z1d,p1d)

       call s_vint13('xx',0,ni+1,0,nj+1,2,nk-2,zph8s(0,0,2),ptbr(0,0,2),&
     &               nlev25,z1d,pt1d)

       call s_vint13('xx',0,ni+1,0,nj+1,2,nk-2,zph8s(0,0,2),qvbr(0,0,2),&
     &               nlev25,z1d,qv1d)

      end if

! -----

      end subroutine s_intrpsnd

!-----7--------------------------------------------------------------7--

      end module m_intrpsnd
