!***********************************************************************
      module m_advbspt
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/07/05,
!                   1999/08/03, 1999/09/30, 1999/10/12, 1999/11/01,
!                   2000/01/17, 2001/04/15, 2001/06/29, 2002/04/02,
!                   2003/04/30, 2003/05/19, 2003/11/28, 2003/12/12,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2010/02/01, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the base state potential temperature advection.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: advbspt, s_advbspt

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface advbspt

        module procedure s_advbspt

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
      subroutine s_advbspt(fpgwmopt,fpdziv,ni,nj,nk,ptbr,rbr,w,ptadv,   &
     &                     pta8w)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpgwmopt
                       ! Formal parameter of unique index of gwmopt

      integer, intent(in) :: fpdziv
                       ! Formal parameter of unique index of dziv

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

! Input and output variable

      real, intent(inout) :: ptadv(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature advection

! Internal shared variables

      integer gwmopt   ! Option for gravity wave mode integration

      real dziv        ! Inverse of dz

      real dziv25      ! 0.25 x dziv

      real, intent(inout) :: pta8w(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature advection
                       ! at w points

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpgwmopt,gwmopt)
      call getrname(fpdziv,dziv)

! -----

! Set the common used variable.

      dziv25=.25e0*dziv

! -----

! Calculate the base state potential temperature advection.

!$omp parallel default(shared) private(k)

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          pta8w(i,j,k)=(rbr(i,j,k-1)+rbr(i,j,k))                        &
     &      *w(i,j,k)*(ptbr(i,j,k-1)-ptbr(i,j,k))*dziv25
        end do
        end do

!$omp end do

      end do

      if(gwmopt.eq.0) then

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            ptadv(i,j,k)=ptadv(i,j,k)+(pta8w(i,j,k)+pta8w(i,j,k+1))
          end do
          end do

!$omp end do

        end do

      else

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            ptadv(i,j,k)=pta8w(i,j,k)+pta8w(i,j,k+1)
          end do
          end do

!$omp end do

        end do

      end if

!$omp end parallel

! -----

      end subroutine s_advbspt

!-----7--------------------------------------------------------------7--

      end module m_advbspt
