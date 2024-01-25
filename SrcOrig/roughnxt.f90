!***********************************************************************
      module m_roughnxt
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/10/15
!     Modification: 2001/10/18, 2002/04/02, 2003/02/05, 2003/04/30,
!                   2003/05/19, 2003/07/15, 2003/12/12, 2004/05/07,
!                   2004/08/20, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2009/06/16, 2009/11/13, 2011/06/01,
!                   2011/11/10, 2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     reset the roughness parameter on the sea surface to the next time
!     step.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: roughnxt, s_roughnxt

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface roughnxt

        module procedure s_roughnxt

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic max

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_roughnxt(ni,nj,land,va,cm,z0m,z0h)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

      real, intent(in) :: va(0:ni+1,0:nj+1)
                       ! Magnitude of velocity at lowest plane

      real, intent(in) :: cm(0:ni+1,0:nj+1)
                       ! Bulk coefficient for velocity

! Input and output variables

      real, intent(inout) :: z0m(0:ni+1,0:nj+1)
                       ! Roughness length for velocity

      real, intent(inout) :: z0h(0:ni+1,0:nj+1)
                       ! Roughness length for scalar

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      real ust         ! Friction velocity

!-----7--------------------------------------------------------------7--

! Calculate the roughness parameter on the sea surface to the next time
! step.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j,ust)

      do j=1,nj-1
      do i=1,ni-1

        if(land(i,j).lt.3) then

          ust=cm(i,j)*va(i,j)

          if(ust.lt.1.08e0) then

            z0m(i,j)=max(-34.7e-6+8.28e-4*ust,z0min)

          else

            z0m(i,j)=max(-.277e-2+3.39e-3*ust,z0min)

          end if

          z0h(i,j)=z0m(i,j)

        end if

      end do
      end do

!$omp end do

!$omp end parallel

! -----

      end subroutine s_roughnxt

!-----7--------------------------------------------------------------7--

      end module m_roughnxt
