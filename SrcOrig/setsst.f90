!***********************************************************************
      module m_setsst
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/11/10
!     Modification: 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the interpolated sea surface temperature.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: setsst, s_setsst

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface setsst

        module procedure s_setsst

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
      subroutine s_setsst(fpsstitv,ird,ni,nj,sst,sstd)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpsstitv
                       ! Formal parameter of unique index of sstitv

      integer, intent(in) :: ird
                       ! Index of count to read out in rdsstnxt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

! Input and output variables

      real, intent(inout) :: sst(0:ni+1,0:nj+1)
                       ! Sea surface temperature of external data
                       ! at marked time

      real, intent(inout) :: sstd(0:ni+1,0:nj+1)
                       ! Time tendency of
                       ! sea surface temperature of external data

! Internal shared variables

      real sstitv      ! Time interval of sea surface temperature data

      real sstiv       ! Inverse of sstitv

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getrname(fpsstitv,sstitv)

! -----

! Set the common used variable.

      sstiv=1.e0/sstitv

! -----

!! Set the interpolated sea surface temperature.

!$omp parallel default(shared)

! Set the time tendency of sea surface temperature at current marked
! time.

      if(ird.eq.1) then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj
        do i=1,ni
          sstd(i,j)=(sstd(i,j)-sst(i,j))*sstiv
        end do
        end do

!$omp end do

      end if

! -----

! Set the sea surface temperature at current marked time.

      if(ird.eq.2) then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj
        do i=1,ni
          sst(i,j)=sstd(i,j)
        end do
        end do

!$omp end do

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_setsst

!-----7--------------------------------------------------------------7--

      end module m_setsst
