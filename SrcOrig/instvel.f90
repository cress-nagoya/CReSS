!***********************************************************************
      module m_instvel
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/07/30
!     Modification: 2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2009/11/13, 2011/09/22, 2011/11/10, 2013/01/28,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the maximum instantaneous wind velocity.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getcname
      use m_inichar

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: instvel, s_instvel

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface instvel

        module procedure s_instvel

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic max
      intrinsic sqrt

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_instvel(fpdmpvar,ni,nj,nk,u,v,w,tke,maxvl)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpdmpvar
                       ! Formal parameter of unique index of dmpvar

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

      real, intent(in) :: tke(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy

! Input and output variable

      real, intent(inout) :: maxvl(0:ni+1,0:nj+1,1:nk)
                       ! Maximum instantaneous wind velocity

! Internal shared variable

      character(len=108) dmpvar
                       ! Control flag of dumped variables

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real u8s2        ! 2.0 x x components of velocity at scalar points
      real v8s2        ! 2.0 x y components of velocity at scalar points
      real w8s2        ! 2.0 x z components of velocity at scalar points

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(dmpvar)

! -----

! Get the required namelist variable.

      call getcname(fpdmpvar,dmpvar)

! -----

! Calculate the maximum instantaneous wind velocity.

      if(dmpvar(12:12).eq.'+'.or.dmpvar(12:12).eq.'-') then

!$omp parallel default(shared) private(k)

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j,u8s2,v8s2,w8s2)

          do j=1,nj-1
          do i=1,ni-1

            u8s2=u(i,j,k)+u(i+1,j,k)
            v8s2=v(i,j,k)+v(i,j+1,k)
            w8s2=w(i,j,k)+w(i,j,k+1)

            maxvl(i,j,k)=max(maxvl(i,j,k),sqrt(2.e0*tke(i,j,k))         &
     &        +.5e0*sqrt(u8s2*u8s2+v8s2*v8s2+w8s2*w8s2))

          end do
          end do

!$omp end do

        end do

!$omp end parallel

      end if

! -----

      end subroutine s_instvel

!-----7--------------------------------------------------------------7--

      end module m_instvel
