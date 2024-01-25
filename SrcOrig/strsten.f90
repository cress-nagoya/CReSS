!***********************************************************************
      module m_strsten
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/07/05
!     Modification: 1999/07/19, 1999/07/21, 1999/08/18, 1999/10/12,
!                   1999/11/01, 2000/01/17, 2000/07/05, 2000/12/19,
!                   2001/05/29, 2001/06/06, 2001/11/20, 2002/04/02,
!                   2003/01/04, 2003/04/30, 2003/05/19, 2003/12/12,
!                   2004/02/01, 2006/02/13, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the stress tensor.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bcten
      use m_comindx
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: strsten, s_strsten

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface strsten

        module procedure s_strsten

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
      subroutine s_strsten(fpsfcopt,ni,nj,nk,ufrc,vfrc,rkh,rkv,         &
     &                     t11,t22,t33,t12,t13,t23,t31,t32)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpsfcopt
                       ! Formal parameter of unique index of sfcopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: ufrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in u equation

      real, intent(in) :: vfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in v equation

      real, intent(in) :: rkh(0:ni+1,0:nj+1,1:nk)
                       ! rbr x horizontal eddy viscosity

      real, intent(in) :: rkv(0:ni+1,0:nj+1,1:nk)
                       ! rbr x vertical eddy viscosity

! Input and output variables

      real, intent(inout) :: t11(0:ni+1,0:nj+1,1:nk)
                       ! x-x components of stress tensor

      real, intent(inout) :: t22(0:ni+1,0:nj+1,1:nk)
                       ! y-y components of stress tensor

      real, intent(inout) :: t33(0:ni+1,0:nj+1,1:nk)
                       ! z-z components of stress tensor

      real, intent(inout) :: t12(0:ni+1,0:nj+1,1:nk)
                       ! x-y components of stress tensor

      real, intent(inout) :: t13(0:ni+1,0:nj+1,1:nk)
                       ! x-z components of stress tensor

      real, intent(inout) :: t23(0:ni+1,0:nj+1,1:nk)
                       ! y-z components of stress tensor

      real, intent(inout) :: t31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of stress tensor

      real, intent(inout) :: t32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of stress tensor

! Internal shared variable

      integer sfcopt   ! Option for surface physics

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpsfcopt,sfcopt)

! -----

!! Calculate the stress tensor.

!$omp parallel default(shared) private(k)

! Calculate the diagonal and the x-y components of the stress tensor.

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          t11(i,j,k)=rkh(i,j,k)*t11(i,j,k)
          t22(i,j,k)=rkh(i,j,k)*t22(i,j,k)
          t33(i,j,k)=rkv(i,j,k)*t33(i,j,k)
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-1
        do i=2,ni-1
          t12(i,j,k)=.25e0*t12(i,j,k)                                   &
     &      *((rkh(i-1,j-1,k)+rkh(i,j,k))+(rkh(i-1,j,k)+rkh(i,j-1,k)))
        end do
        end do

!$omp end do

      end do

! -----

! Calculate the x-z and the y-z components of the stress tensor.

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-1
          t13(i,j,k)=.25e0*t13(i,j,k)                                   &
     &      *((rkv(i-1,j,k-1)+rkv(i,j,k))+(rkv(i-1,j,k)+rkv(i,j,k-1)))
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-1
        do i=2,ni-2
          t23(i,j,k)=.25e0*t23(i,j,k)                                   &
     &      *((rkv(i,j-1,k-1)+rkv(i,j,k))+(rkv(i,j-1,k)+rkv(i,j,k-1)))
        end do
        end do

!$omp end do

      end do

      if(sfcopt.ge.1) then

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-1
          t13(i,j,2)=ufrc(i,j,1)
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-1
        do i=2,ni-2
          t23(i,j,2)=vfrc(i,j,1)
        end do
        end do

!$omp end do

      end if

! -----

! Calculate the z-x and the z-y components of the stress tensor.

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-1
          t31(i,j,k)=.25e0*t31(i,j,k)                                   &
     &      *((rkh(i-1,j,k-1)+rkh(i,j,k))+(rkh(i-1,j,k)+rkh(i,j,k-1)))
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-1
        do i=2,ni-2
          t32(i,j,k)=.25e0*t32(i,j,k)                                   &
     &      *((rkh(i,j-1,k-1)+rkh(i,j,k))+(rkh(i,j-1,k)+rkh(i,j,k-1)))
        end do
        end do

!$omp end do

      end do

! -----

!$omp end parallel

!! -----

! Set the boundary conditions.

      call bcten(idbbc,idtbc,ni,nj,nk,t31)
      call bcten(idbbc,idtbc,ni,nj,nk,t32)

! -----

      end subroutine s_strsten

!-----7--------------------------------------------------------------7--

      end module m_strsten
