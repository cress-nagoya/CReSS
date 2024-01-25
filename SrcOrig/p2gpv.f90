!***********************************************************************
      module m_p2gpv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/03/13
!     Modification: 2001/04/15, 2001/05/29, 2001/06/29, 2002/04/02,
!                   2002/08/15, 2002/09/09, 2002/10/31, 2003/04/30,
!                   2003/05/19, 2003/11/28, 2003/12/12, 2006/01/10,
!                   2006/05/12, 2006/09/21, 2007/05/07, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2009/03/23,
!                   2011/09/22, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform the analysis nudging to GPV data of pressure.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: p2gpv, s_p2gpv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface p2gpv

        module procedure s_p2gpv

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
      subroutine s_p2gpv(nggdmp,gtinc,ni,nj,nk,jcb,ppp,ppgpv,pptd,pfrc)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: nggdmp
                       ! Analysis nudging damping coefficient

      real, intent(in) :: gtinc
                       ! Lapse of forecast time from GPV data reading

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jabobian

      real, intent(in) :: ppp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at past

      real, intent(in) :: ppgpv(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation of GPV data
                       ! at marked time

      real, intent(in) :: pptd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! pressure perturbation of GPV data

! Input and output variable

      real, intent(inout) :: pfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in pressure equation

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Calculate the analysis nudging terms for pressure.

!$omp parallel default(shared) private(k)

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          pfrc(i,j,k)=pfrc(i,j,k)+nggdmp                                &
     &      *jcb(i,j,k)*((ppgpv(i,j,k)+pptd(i,j,k)*gtinc)-ppp(i,j,k))
        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_p2gpv

!-----7--------------------------------------------------------------7--

      end module m_p2gpv
