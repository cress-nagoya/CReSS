!***********************************************************************
      module m_trilat
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/03/25, 1999/04/06, 1999/05/10, 1999/08/23,
!                   1999/09/30, 1999/10/12, 1999/11/01, 2000/01/17,
!                   2000/07/05, 2000/12/18, 2001/03/13, 2002/04/02,
!                   2003/01/04, 2003/04/30, 2003/05/19, 2003/10/31,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the Coriolis parameters x 0.25.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comphy
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: trilat, s_trilat

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface trilat

        module procedure s_trilat

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic sin
      intrinsic sqrt

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_trilat(fpcoropt,ni,nj,lat,fc)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpcoropt
                       ! Formal parameter of unique index of coropt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      real, intent(in) :: lat(0:ni+1,0:nj+1)
                       ! Latitude

! Output variables

      real, intent(out) :: fc(0:ni+1,0:nj+1,1:2)
                       ! 0.25 x Coriolis parameters

! Internal shared variables

      integer coropt   ! Option for Coriolis force

      real omega5      ! 0.5 x omega

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      real sinlat      ! sin(latitude)

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpcoropt,coropt)

! -----

! Calculate 0.5 omega.

      omega5=.5e0*omega

! -----

! Calculate the Coriolis parameters x 0.25.

!$omp parallel default(shared)

      if(coropt.eq.1) then

!$omp do schedule(runtime) private(i,j,sinlat)

        do j=1,nj-1
        do i=1,ni-1
          sinlat=sin(lat(i,j)*d2r)

          fc(i,j,1)=omega5*sinlat
          fc(i,j,2)=0.e0

        end do
        end do

!$omp end do

      else if(coropt.eq.2) then

!$omp do schedule(runtime) private(i,j,sinlat)

        do j=1,nj-1
        do i=1,ni-1
          sinlat=sin(lat(i,j)*d2r)

          fc(i,j,1)=omega5*sinlat
          fc(i,j,2)=omega5*sqrt(1.e0-sinlat*sinlat)

        end do
        end do

!$omp end do

      end if

!$omp end parallel

! -----

      end subroutine s_trilat

!-----7--------------------------------------------------------------7--

      end module m_trilat
