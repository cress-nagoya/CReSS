!***********************************************************************
      module m_totalq
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/01/10
!     Modification: 2006/09/11, 2006/09/21, 2006/09/30, 2007/05/14,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/01/30,
!                   2009/02/27, 2011/08/18, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the total mixing ratio or concentrations of same material.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: totalq, s_totalq

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface totalq

        module procedure s_totalq

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
      subroutine s_totalq(istr,iend,ni,nj,nk,nctg,q,qall)
!***********************************************************************

! Input variables

      integer, intent(in) :: istr
                       ! Start index of categories or types

      integer, intent(in) :: iend
                       ! End index of categories or types

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nctg
                       ! Number of categories or types
                       ! of optional material

      real, intent(in) :: q(0:ni+1,0:nj+1,1:nk,1:nctg)
                       ! Mixing ratio or concentrations
                       ! of optional material

! Output variable

      real, intent(out) :: qall(0:ni+1,0:nj+1,1:nk)
                       ! Total mixing ratio

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer n        ! Array index in 4th direction

!-----7--------------------------------------------------------------7--

! Get the total mixing ratio or concentrations of same material.

!$omp parallel default(shared) private(k,n)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          qall(i,j,k)=q(i,j,k,istr)
        end do
        end do

!$omp end do

      end do

      if(istr.lt.iend) then

        do n=istr+1,iend

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              qall(i,j,k)=qall(i,j,k)+q(i,j,k,n)
            end do
            end do

!$omp end do

          end do

        end do

      end if

!$omp end parallel

! -----

      end subroutine s_totalq

!-----7--------------------------------------------------------------7--

      end module m_totalq
