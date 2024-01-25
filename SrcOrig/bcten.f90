!***********************************************************************
      module m_bcten
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/07/05
!     Modification: 1999/08/23, 1999/09/30, 1999/10/12, 1999/11/01,
!                   1999/11/19, 2000/01/17, 2001/06/06, 2001/12/11,
!                   2002/04/02, 2002/07/23, 2003/04/30, 2003/05/19,
!                   2003/12/12, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the boundary conditions for optional tensor.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: bcten, s_bcten

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface bcten

        module procedure s_bcten

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
      subroutine s_bcten(fpbbc,fptbc,ni,nj,nk,ten)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpbbc
                       ! Formal parameter of unique index of bbc

      integer, intent(in) :: fptbc
                       ! Formal parameter of unique index of tbc

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

! Input and output variable

      real, intent(inout) :: ten(0:ni+1,0:nj+1,1:nk)
                       ! Optional tensor

! Internal shared variables

      integer bbc      ! Option for bottom boundary conditions
      integer tbc      ! Option for top boundary conditions

      integer nkm2     ! nk - 2

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpbbc,bbc)
      call getiname(fptbc,tbc)

! -----

! Set the common used variable.

      nkm2=nk-2

! -----

!! Set the bottom and top boundary conditions.

!$omp parallel default(shared)

! Set the bottom boundary conditions.

      if(bbc.eq.2) then

!$omp do schedule(runtime) private(i,j)

        do j=0,nj+1
        do i=0,ni+1
          ten(i,j,1)=-ten(i,j,3)
        end do
        end do

!$omp end do

      else if(bbc.ge.3) then

!$omp do schedule(runtime) private(i,j)

        do j=0,nj+1
        do i=0,ni+1
          ten(i,j,1)=ten(i,j,3)
        end do
        end do

!$omp end do

      end if

! -----

! Set the top boundary conditions.

      if(tbc.eq.2) then

!$omp do schedule(runtime) private(i,j)

        do j=0,nj+1
        do i=0,ni+1
          ten(i,j,nk)=-ten(i,j,nkm2)
        end do
        end do

!$omp end do

      else if(tbc.ge.3) then

!$omp do schedule(runtime) private(i,j)

        do j=0,nj+1
        do i=0,ni+1
          ten(i,j,nk)=ten(i,j,nkm2)
        end do
        end do

!$omp end do

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_bcten

!-----7--------------------------------------------------------------7--

      end module m_bcten
