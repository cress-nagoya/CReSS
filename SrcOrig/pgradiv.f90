!***********************************************************************
      module m_pgradiv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/12/20
!     Modification: 2000/01/17, 2001/06/29, 2002/04/02, 2003/04/30,
!                   2003/05/19, 2003/12/12, 2003/12/26, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the pressure gradient force vertically with the
!     horizontally explicit and vertically implicit method.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: pgradiv, s_pgradiv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface pgradiv

        module procedure s_pgradiv

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
      subroutine s_pgradiv(fpdziv,fpweicoe,dts,ni,nj,nk,jcb,wfrc,fp,fw, &
     &                     fpdvj)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpdziv
                       ! Formal parameter of unique index of dziv

      integer, intent(in) :: fpweicoe
                       ! Formal parameter of unique index of weicoe

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dts
                       ! Small time steps interval

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: wfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in w equation in large time steps

      real, intent(in) :: fp(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in pressure equation

! Input and output variable

      real, intent(inout) :: fw(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in w equation

! Internal shared variables

      real dziv        ! Inverce of dz

      real weicoe      ! Weighting coefficient for implicit method

      real dtdzw       ! dts x dziv x weicoe

      real, intent(inout) :: fpdvj(0:ni+1,0:nj+1,1:nk)
                       ! fp / jcb

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getrname(fpdziv,dziv)
      call getrname(fpweicoe,weicoe)

! -----

! Set the common used variable.

      dtdzw=dts*dziv*weicoe

! -----

! Calculate the pressure gradient force vertically.

!$omp parallel default(shared) private(k)

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          fpdvj(i,j,k)=fp(i,j,k)/jcb(i,j,k)
        end do
        end do

!$omp end do

      end do

      do k=3,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          fw(i,j,k)                                                     &
     &      =wfrc(i,j,k)+fw(i,j,k)+(fpdvj(i,j,k-1)-fpdvj(i,j,k))*dtdzw
        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_pgradiv

!-----7--------------------------------------------------------------7--

      end module m_pgradiv
