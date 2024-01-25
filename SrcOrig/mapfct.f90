!***********************************************************************
      module m_mapfct
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/03/17, 1999/03/25, 1999/04/06, 1999/06/07,
!                   1999/08/23, 1999/09/30, 1999/10/12, 1999/11/01,
!                   2000/01/17, 2000/07/05, 2001/02/13, 2001/05/29,
!                   2001/11/20, 2002/04/02, 2003/01/04, 2003/03/21,
!                   2003/04/30, 2003/05/19, 2003/10/31, 2003/11/05,
!                   2003/12/12, 2004/05/31, 2004/06/10, 2004/08/20,
!                   2004/09/01, 2005/01/31, 2005/02/10, 2006/04/03,
!                   2006/11/06, 2006/12/04, 2007/01/05, 2007/01/20,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2011/08/09, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the map scale factors.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bcyclex
      use m_bcycley
      use m_combuf
      use m_comindx
      use m_commath
      use m_getbufgx
      use m_getbufgy
      use m_getbufsx
      use m_getbufsy
      use m_getiname
      use m_getrname
      use m_putbufgx
      use m_putbufgy
      use m_putbufsx
      use m_putbufsy
      use m_shiftgx
      use m_shiftgy
      use m_shiftsx
      use m_shiftsy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: mapfct, s_mapfct

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface mapfct

        module procedure s_mapfct

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic cos
      intrinsic sin
      intrinsic tan
      intrinsic exp
      intrinsic log
      intrinsic real
      intrinsic sqrt

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_mapfct(fpmpopt,fpnspol,fpadvopt,fptubopt,            &
     &                    fpdisr,fpdxiv,fpdyiv,cpj,ni,nj,x,lat,         &
     &                    mf,mf8u,mf8v,rmf,rmf8u,rmf8v,tmp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpnspol
                       ! Formal parameter of unique index of nspol

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fptubopt
                       ! Formal parameter of unique index of tubopt

      integer, intent(in) :: fpdisr
                       ! Formal parameter of unique index of disr

      integer, intent(in) :: fpdxiv
                       ! Formal parameter of unique index of dxiv

      integer, intent(in) :: fpdyiv
                       ! Formal parameter of unique index of dyiv

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      real, intent(in) :: cpj(1:7)
                       ! Map projection parameters

      real, intent(in) :: x(0:ni+1)
                       ! x coordinates at scalar points

      real, intent(in) :: lat(0:ni+1,0:nj+1)
                       ! Latitude

! Output variables

      real, intent(out) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(out) :: mf8u(0:ni+1,0:nj+1)
                       ! Map scale factors at u points

      real, intent(out) :: mf8v(0:ni+1,0:nj+1)
                       ! Map scale factors at v points

      real, intent(out) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

      real, intent(out) :: rmf8u(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at u points

      real, intent(out) :: rmf8v(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at v points

! Internal shared variables

      integer mpopt    ! Option for map projection
      integer nspol    ! Option for projected region

      integer advopt   ! Option for advection scheme
      integer tubopt   ! Option for turbulent mixing

      real disr        ! Distance from center of circular cylinder
                       ! to origin of calculation domain

      real dxiv        ! Inverse of dx
      real dyiv        ! Inverse of dy

      real dxv625      ! 0.0625 x dxiv
      real dyv625      ! 0.0625 x dyiv

      real pol05       ! nspol / 2.0
      real pold2r      ! real(nspol) x d2r

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpmpopt,mpopt)
      call getiname(fpnspol,nspol)
      call getiname(fpadvopt,advopt)
      call getiname(fptubopt,tubopt)
      call getrname(fpdisr,disr)
      call getrname(fpdxiv,dxiv)
      call getrname(fpdyiv,dyiv)

! -----

! Set the common used variables.

      pol05=.5e0*real(nspol)
      pold2r=real(nspol)*d2r

      dxv625=.0625e0*dxiv
      dyv625=.0625e0*dyiv

! -----

!!! Calculate the map scale factors.

!$omp parallel default(shared)

!! Calculate the map scale factor at scalar points.

! Calculate the map scale factor with the spherical coordinates system.

      if(mpopt.eq.0.or.mpopt.eq.10) then

!$omp do schedule(runtime) private(i,j)

        do j=0,nj
        do i=0,ni
          mf(i,j)=1.e0/(cos(lat(i,j)*d2r)+eps)
        end do
        end do

!$omp end do

! -----

! Calculate the map scale factor with the Polar Stereographic
! projection method.

      else if(mpopt.eq.1) then

!$omp do schedule(runtime) private(i,j)

        do j=0,nj
        do i=0,ni
          mf(i,j)=cpj(1)/(1.e0+sin(pold2r*lat(i,j))+eps)
        end do
        end do

!$omp end do

! -----

! Calculate the map scale factor with the Lambert Conformal Conic
! projection method.

      else if(mpopt.eq.2) then

!$omp do schedule(runtime) private(i,j)

        do j=0,nj
        do i=0,ni
          mf(i,j)=(cpj(1)/(cos(lat(i,j)*d2r)+eps))                      &
     &      *exp(cpj(4)*log(tan((45.e0-pol05*lat(i,j))*d2r)*cpj(3)+eps))
        end do
        end do

!$omp end do

! -----

! Calculate the map scale factor with the Mercator projection method.

      else if(mpopt.eq.3.or.mpopt.eq.13) then

!$omp do schedule(runtime) private(i,j)

        do j=0,nj
        do i=0,ni
          mf(i,j)=cpj(1)/(cos(lat(i,j)*d2r)+eps)
        end do
        end do

!$omp end do

! -----

! Calculate the map scale factor without any projection method.

      else if(mpopt.eq.4) then

!$omp do schedule(runtime) private(i,j)

        do j=0,nj
        do i=0,ni
          mf(i,j)=1.e0
        end do
        end do

!$omp end do

! -----

! Calculate the map scale factor with the circular cylinder coordinates
! system.

      else if(mpopt.eq.5) then

!$omp do schedule(runtime) private(i,j)

        do j=0,nj
        do i=0,ni
          mf(i,j)=1.e0/(x(i)+disr)
        end do
        end do

!$omp end do

      end if

! -----

!! -----

! Calculate the map scale factors at the u and v points.

!$omp do schedule(runtime) private(i,j)

      do j=0,nj
      do i=1,ni
        mf8u(i,j)=.5e0*(mf(i-1,j)+mf(i,j))
      end do
      end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

      do j=1,nj
      do i=0,ni
        mf8v(i,j)=.5e0*(mf(i,j-1)+mf(i,j))
      end do
      end do

!$omp end do

! -----

! Calculate the map scale factor squared and the inverse of map scale
! factor.

!$omp do schedule(runtime) private(i,j)

      do j=0,nj
      do i=0,ni
        rmf(i,j,1)=mf(i,j)*mf(i,j)
        rmf(i,j,2)=1.e0/mf(i,j)
        rmf(i,j,3)=rmf(i,j,2)*rmf(i,j,2)
        rmf(i,j,4)=sqrt(rmf(i,j,2))
      end do
      end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

      do j=0,nj
      do i=1,ni
        rmf8u(i,j,1)=mf8u(i,j)*mf8u(i,j)
        rmf8u(i,j,2)=1.e0/mf8u(i,j)
      end do
      end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

      do j=1,nj
      do i=0,ni
        rmf8v(i,j,1)=mf8v(i,j)*mf8v(i,j)
        rmf8v(i,j,2)=1.e0/mf8v(i,j)
      end do
      end do

!$omp end do

! -----

! Calculate the differential of map scale factor x 0.0625.

!$omp do schedule(runtime) private(i,j)

      do j=0,nj
      do i=1,ni
        tmp1(i,j)=(mf(i,j)-mf(i-1,j))*dxv625
      end do
      end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

      do j=0,nj
      do i=1,ni-1
        rmf8u(i,j,3)=tmp1(i,j)+tmp1(i+1,j)
      end do
      end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

      do j=1,nj
      do i=0,ni
        tmp1(i,j)=(mf(i,j)-mf(i,j-1))*dyv625
      end do
      end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=0,ni
        rmf8v(i,j,3)=tmp1(i,j)+tmp1(i,j+1)
      end do
      end do

!$omp end do

! -----

!$omp end parallel

!!! -----

!! Exchange the value in the case the 4th order calculation is
!! performed.

      if(advopt.ge.2.or.tubopt.ge.1) then

! Exchange the value in x direction.

        call s_putbufsx(idwbc,idebc,'all',4,ni-3,ni,nj,1,mf8u,1,1,sbuf)

        call s_shiftsx(idwbc,idebc,'all',nj,1,1,sbuf,rbuf)

        call s_getbufsx(idwbc,idebc,'all',0,ni+1,ni,nj,1,mf8u,1,1,rbuf)

        call s_putbufgx(idwbc,idebc,'all',4,ni-3,ni,nj,1,mf8u,1,1,sbuf)

        call s_shiftgx(idwbc,idebc,'all',nj,1,1,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'all',0,ni+1,ni,nj,1,mf8u,1,1,rbuf)

        call s_bcyclex(idwbc,idebc,4,0,ni-3,ni+1,ni,nj,1,mf8u)

! -----

! Exchange the value in y direction.

        call s_putbufsy(idsbc,idnbc,'all',4,nj-3,ni,nj,1,mf8v,1,1,sbuf)

        call s_shiftsy(idsbc,idnbc,'all',ni,1,1,sbuf,rbuf)

        call s_getbufsy(idsbc,idnbc,'all',0,nj+1,ni,nj,1,mf8v,1,1,rbuf)

        call s_putbufgy(idsbc,idnbc,'all',4,nj-3,ni,nj,1,mf8v,1,1,sbuf)

        call s_shiftgy(idsbc,idnbc,'all',ni,1,1,sbuf,rbuf)

        call s_getbufgy(idsbc,idnbc,'all',0,nj+1,ni,nj,1,mf8v,1,1,rbuf)

        call s_bcycley(idsbc,idnbc,4,0,nj-3,nj+1,ni,nj,1,mf8v)

! -----

      end if

!! -----

      end subroutine s_mapfct

!-----7--------------------------------------------------------------7--

      end module m_mapfct
