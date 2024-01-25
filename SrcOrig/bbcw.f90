!***********************************************************************
      module m_bbcw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/07/01
!     Modification: 2006/11/06, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2011/08/09, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the bottom boundary conditions for the z components of
!     velocity.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: bbcw, s_bbcw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface bbcw

        module procedure s_bbcw

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
      subroutine s_bbcw(fpmpopt,fpmfcopt,ni,nj,nk,j31,j32,mf,aa,        &
     &                  uf,vf,wf,j31u2,j32v2)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: j31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of Jacobian

      real, intent(in) :: j32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of Jacobian

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: aa(0:ni+1,0:nj+1,1:nk)
                       ! Coefficient matrix

      real, intent(in) :: uf(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at future

      real, intent(in) :: vf(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at future

! Input and output variable

      real, intent(inout) :: wf(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at future

! Internal shared variables

      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor

      real, intent(inout) :: j31u2(0:ni+1,0:nj+1)
                       ! 2.0 x j31 x u

      real, intent(inout) :: j32v2(0:ni+1,0:nj+1)
                       ! 2.0 x j32 x v

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)

! -----

! Set the bottom boundary conditions.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=1,ni
        j31u2(i,j)=(uf(i,j,1)+uf(i,j,2))*j31(i,j,2)
      end do
      end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

      do j=1,nj
      do i=1,ni-1
        j32v2(i,j)=(vf(i,j,1)+vf(i,j,2))*j32(i,j,2)
      end do
      end do

!$omp end do

      if(mfcopt.eq.0) then

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          wf(i,j,3)=wf(i,j,3)+.25e0*aa(i,j,3)                           &
     &      *((j31u2(i,j)+j31u2(i+1,j))+(j32v2(i,j)+j32v2(i,j+1)))
        end do
        end do

!$omp end do

      else

        if(mpopt.eq.0.or.mpopt.eq.10) then

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            wf(i,j,3)=wf(i,j,3)                                         &
     &        +.25e0*aa(i,j,3)*(mf(i,j)*(j31u2(i,j)+j31u2(i+1,j))       &
     &        +(j32v2(i,j)+j32v2(i,j+1)))
          end do
          end do

!$omp end do

        else if(mpopt.eq.5) then

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            wf(i,j,3)=wf(i,j,3)                                         &
     &        +.25e0*aa(i,j,3)*((j31u2(i,j)+j31u2(i+1,j))               &
     &        +mf(i,j)*(j32v2(i,j)+j32v2(i,j+1)))
          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            wf(i,j,3)=wf(i,j,3)+.25e0*mf(i,j)*aa(i,j,3)                 &
     &        *((j31u2(i,j)+j31u2(i+1,j))+(j32v2(i,j)+j32v2(i,j+1)))
          end do
          end do

!$omp end do

        end if

      end if

!$omp end parallel

! -----

      end subroutine s_bbcw

!-----7--------------------------------------------------------------7--

      end module m_bbcw
