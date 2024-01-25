!***********************************************************************
      module m_gsmoou
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/03/05
!     Modification: 2004/04/15, 2004/08/01, 2004/09/10, 2004/09/25,
!                   2005/02/10, 2006/09/30, 2004/12/04, 2007/01/05,
!                   2007/01/20, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2011/09/22, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     smooth the x components of velocity of interpolated GPV data.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bc4news
      use m_bcgsmu
      use m_bcycle
      use m_combuf
      use m_comindx
      use m_commath
      use m_getbufgx
      use m_getbufgy
      use m_getbufsx
      use m_getbufsy
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

      public :: gsmoou, s_gsmoou

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface gsmoou

        module procedure s_gsmoou

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
      subroutine s_gsmoou(fpgsmcoe,ni,nj,nk,ugpv,dfu)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpgsmcoe
                       ! Formal parameter of unique index of gsmcoe

      integer, intent(in) :: ni
                       ! Model dimensionn in x direction

      integer, intent(in) :: nj
                       ! Model dimensionn in y direction

      integer, intent(in) :: nk
                       ! Model dimensionn in z direction

! Input and output variable

      real, intent(inout) :: ugpv(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity of GPV data

! Internal shared variables

      integer ni_sub   ! Substitute for ni

      real gsmcoe      ! Smoothing coefficient

      real dtcoe       ! Time interval of iteration x gsmcoe

      real, intent(inout) :: dfu(0:ni+1,0:nj+1,1:nk)
                       ! Diffusion term
                       ! of x components of velocity of GPV

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in y direction

      real a           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getrname(fpgsmcoe,gsmcoe)

! -----

! Set the common used variable.

      dtcoe=twod15*gsmcoe

! -----

! Set the substituted variable.

      ni_sub=ni

! -----

!! Perform the 2nd order smoothing.

! Calculate the diffusion term.

!$omp parallel default(shared) private(k)

      do k=3,nk-3

!$omp do schedule(runtime) private(i,j,a)

        do j=2,nj-2
        do i=2,ni-1
          a=2.e0*ugpv(i,j,k)

          dfu(i,j,k)=(ugpv(i,j,k+1)+ugpv(i,j,k-1)-a)                    &
     &      +((ugpv(i+1,j,k)+ugpv(i-1,j,k)-a)                           &
     &      +(ugpv(i,j+1,k)+ugpv(i,j-1,k)-a))

        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

! Exchange the value horizontally.

      call s_putbufsx(idwbc,idebc,'all',3,ni-2,ni,nj,nk,dfu,1,1,sbuf)

      call s_shiftsx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

      call s_getbufsx(idwbc,idebc,'all',1,ni_sub,ni,nj,nk,dfu,1,1,rbuf)

      call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,dfu,1,1,sbuf)

      call s_shiftsy(idsbc,idnbc,'all',ni,nk,1,sbuf,rbuf)

      call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,dfu,1,1,rbuf)

      call s_putbufgx(idwbc,idebc,'all',3,ni-2,ni,nj,nk,dfu,1,1,sbuf)

      call s_shiftgx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

      call s_getbufgx(idwbc,idebc,'all',1,ni_sub,ni,nj,nk,dfu,1,1,rbuf)

      call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,dfu,1,1,sbuf)

      call s_shiftgy(idsbc,idnbc,'all',ni,nk,1,sbuf,rbuf)

      call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,dfu,1,1,rbuf)

      call s_putbufgx(idwbc,idebc,'all',3,ni-2,ni,nj,nk,dfu,1,1,sbuf)

      call s_shiftgx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

      call s_getbufgx(idwbc,idebc,'all',1,ni_sub,ni,nj,nk,dfu,1,1,rbuf)

! -----

! Set the periodic boundary conditions.

      call bcycle(idwbc,idebc,idsbc,idnbc,                              &
     &            3,1,ni-2,ni_sub,2,1,nj-2,nj-1,ni,nj,nk,dfu)

! -----

! Set the boundary conditions at the four corners.

      call bc4news(idwbc,idebc,idsbc,idnbc,1,ni_sub,1,nj-1,ni,nj,nk,dfu)

! -----

! Set the lateral boundary conditions.

      call bcgsmu(idwbc,idebc,idsbc,idnbc,ni,nj,nk,dfu)

! -----

! Update the GPV data.

!$omp parallel default(shared) private(k)

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni
          ugpv(i,j,k)=ugpv(i,j,k)+dtcoe*dfu(i,j,k)
        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

!! -----

      end subroutine s_gsmoou

!-----7--------------------------------------------------------------7--

      end module m_gsmoou
