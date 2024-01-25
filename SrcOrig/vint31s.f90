!***********************************************************************
      module m_vint31s
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/04/15
!     Modification: 2001/05/29, 2002/04/02, 2002/06/18, 2003/04/30,
!                   2003/05/19, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     interpolate the 3 dimensional input variable to the 1 dimensional
!     flat plane.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: vint31s, s_vint31s

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vint31s

        module procedure s_vint31s

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
      subroutine s_vint31s(ni,nj,nk,z1d,zph8s,invar,outvar)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: z1d(1:nk)
                       ! Output z coordinates at scalar points

      real, intent(in) :: zph8s(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates at scalar points

      real, intent(in) :: invar(0:ni+1,0:nj+1,1:nk)
                       ! Optional input variable

! Output variable

      real, intent(out) :: outvar(0:ni+1,0:nj+1,1:nk)
                       ! Optional interpolated variable

! Internal shared variable

      integer nkm2     ! nk - 2

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer ki       ! Input array index in z direction

      real dk          ! Distance in z direction
                       ! between model and data points

!-----7--------------------------------------------------------------7--

! Set the common used variable.

      nkm2=nk-2

! -----

!! Interpolate the 3 dimensional input variable to the 1 dimensional
!! flat plane.

!$omp parallel default(shared) private(k,ki)

! Fill in the undifined value outside of the flat plane.

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2

          if(zph8s(i,j,2).gt.z1d(k).or.zph8s(i,j,nkm2).le.z1d(k)) then

            outvar(i,j,k)=lim35n

          end if

        end do
        end do

!$omp end do

      end do

! -----

! Interpolate the 3 dimensional input variable.

      do ki=2,nk-3

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j,dk)

          do j=2,nj-2
          do i=2,ni-2

            if(zph8s(i,j,ki).le.z1d(k)                                  &
     &        .and.zph8s(i,j,ki+1).gt.z1d(k)) then

              dk=(z1d(k)-zph8s(i,j,ki))/(zph8s(i,j,ki+1)-zph8s(i,j,ki))

              outvar(i,j,k)=(1.e0-dk)*invar(i,j,ki)+dk*invar(i,j,ki+1)

            end if

          end do
          end do

!$omp end do

        end do

      end do

! -----

!$omp end parallel

!! -----

      end subroutine s_vint31s

!-----7--------------------------------------------------------------7--

      end module m_vint31s
