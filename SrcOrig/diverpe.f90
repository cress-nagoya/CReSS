!***********************************************************************
      module m_diverpe
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/07/05,
!                   1999/08/03, 1999/08/23, 1999/10/12, 1999/11/01,
!                   1999/12/17, 2000/01/17, 2000/12/18, 2001/02/24,
!                   2001/06/29, 2002/04/02, 2003/01/04, 2003/03/21,
!                   2003/04/30, 2003/05/19, 2004/06/10, 2006/11/06,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the divergence in the pressure equation with the
!     horizontally explicit and vertically explicit method.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_diver3d

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: diverpe, s_diverpe

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface diverpe

        module procedure s_diverpe

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
      subroutine s_diverpe(ni,nj,nk,jcb8u,jcb8v,jcb8w,mf,rmf,           &
     &                     rmf8u,rmf8v,rcsq,u,v,wc,pdiv,tmp1,tmp2,tmp3)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: jcb8u(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at u points

      real, intent(in) :: jcb8v(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at v points

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

      real, intent(in) :: rmf8u(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at u points

      real, intent(in) :: rmf8v(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at v points

      real, intent(in) :: rcsq(0:ni+1,0:nj+1,1:nk)
                       ! rbr x sound wave speed squared

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(in) :: wc(0:ni+1,0:nj+1,1:nk)
                       ! zeta components of contravariant velocity

! Output variable

      real, intent(out) :: pdiv(0:ni+1,0:nj+1,1:nk)
                       ! Divergence value in pressure quation

! Internal shared variables

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

! Remark

!     pdiv: This variable is also temporary.

!-----7--------------------------------------------------------------7--

! Calculate the 3 dimensional divergence.

      call diver3d(idmpopt,idmfcopt,iddxiv,iddyiv,iddziv,ni,nj,nk,      &
     &             mf,rmf,rmf8u,rmf8v,jcb8u,jcb8v,jcb8w,u,v,wc,         &
     &             pdiv,tmp1,tmp2,tmp3)

! -----

! Finally get the 3 dimensional divergence in the pressure equation.

!$omp parallel default(shared) private(k)

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          pdiv(i,j,k)=rcsq(i,j,k)*pdiv(i,j,k)
        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_diverpe

!-----7--------------------------------------------------------------7--

      end module m_diverpe
