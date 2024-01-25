!***********************************************************************
      module m_repsit3d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2005/07/01
!     Modification: 2005/08/05, 2006/09/21, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     reposition the restructed variable form original restart variable.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: repsit3d, s_repsit3d

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface repsit3d

        module procedure s_repsit3d

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
      subroutine s_repsit3d(xo,istr,iend,jstr,jend,                     &
     &                      di,dj,ni,nj,nk,ni_rst,nj_rst,var,var_rst)
!***********************************************************************

! Input variables

      character(len=2), intent(in) :: xo
                       ! Control flag of variable arrangement

      integer, intent(in) :: istr
                       ! Minimum do loops index in x direction

      integer, intent(in) :: iend
                       ! Maximum do loops index in x direction

      integer, intent(in) :: jstr
                       ! Minimum do loops index in y direction

      integer, intent(in) :: jend
                       ! Maximum do loops index in y direction

      integer, intent(in) :: di
                       ! Differential index to istr

      integer, intent(in) :: dj
                       ! Differential index to jstr

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: ni_rst
                       ! Restructed files dimension in x direction

      integer, intent(in) :: nj_rst
                       ! Restructed files dimension in y direction

      real, intent(in) :: var(0:ni+1,0:nj+1,1:nk)
                       ! Optional variable

! Input and output variable

      real, intent(inout) :: var_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! var in restructed domain

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Reposition the restructed variable form original restart variable.

!$omp parallel default(shared) private(k)

      if(xo(1:2).eq.'ox') then

        do k=1,nk

!$omp do schedule(runtime) private(i,j)

          do j=jstr,jend
          do i=istr,iend+1
            var_rst(di+i,dj+j,k)=var(i,j,k)
          end do
          end do

!$omp end do

        end do

      else if(xo(1:2).eq.'xo') then

        do k=1,nk

!$omp do schedule(runtime) private(i,j)

          do j=jstr,jend+1
          do i=istr,iend
            var_rst(di+i,dj+j,k)=var(i,j,k)
          end do
          end do

!$omp end do

        end do

      else

        do k=1,nk

!$omp do schedule(runtime) private(i,j)

          do j=jstr,jend
          do i=istr,iend
            var_rst(di+i,dj+j,k)=var(i,j,k)
          end do
          end do

!$omp end do

        end do

      end if

!$omp end parallel

! -----

      end subroutine s_repsit3d

!-----7--------------------------------------------------------------7--

      end module m_repsit3d
