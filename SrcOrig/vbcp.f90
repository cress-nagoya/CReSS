!***********************************************************************
      module m_vbcp
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/08/18,
!                   1999/08/23, 1999/09/06, 1999/10/12, 1999/12/06,
!                   2000/01/17, 2000/03/17, 2001/04/15, 2001/06/29,
!                   2001/07/13, 2001/08/07, 2001/12/11, 2002/04/02,
!                   2002/07/23, 2002/08/15, 2003/04/30, 2003/05/19,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the vertical boundary conditions for the pressure.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: vbcp, s_vbcp

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vbcp

        module procedure s_vbcp

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
      subroutine s_vbcp(fpbbc,ni,nj,nk,ppf)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpbbc
                       ! Formal parameter of unique index of bbc

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

! Input and output variable

      real, intent(inout) :: ppf(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at future

! Internal shared variables

      integer bbc      ! Option for bottom boundary conditions

      integer nkm1     ! nk - 1
      integer nkm2     ! nk - 2

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpbbc,bbc)

! -----

! Set the common used variables.

      nkm1=nk-1
      nkm2=nk-2

! -----

!! Set the bottom and top boundary conditions.

!$omp parallel default(shared)

! Set the bottom boundary conditions.

      if(bbc.eq.2) then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          ppf(i,j,1)=2.e0*ppf(i,j,2)-ppf(i,j,3)
        end do
        end do

!$omp end do

      else if(bbc.eq.3) then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          ppf(i,j,1)=ppf(i,j,2)
        end do
        end do

!$omp end do

      end if

! -----

! Set the top boundary conditions.

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=1,ni-1
        ppf(i,j,nkm1)=ppf(i,j,nkm2)
      end do
      end do

!$omp end do

! -----

!$omp end parallel

!! -----

      end subroutine s_vbcp

!-----7--------------------------------------------------------------7--

      end module m_vbcp
