!***********************************************************************
      module m_getqvs
!***********************************************************************

!     Author      : Mizutani Fumihiko, Sakakibara Atsushi
!     Date        : 2004/03/16
!     Modification: 2007/05/07, 2007/10/19, 2008/05/02, 2008/07/01,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the saturation mixing ratio.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getqvs, s_getqvs

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getqvs

        module procedure s_getqvs

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic log

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_getqvs(ni,nj,nk,pbr,ptbr,pp,ptp,qvs)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation

      real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

! Output variable

      real, intent(out) :: qvs(0:ni+1,0:nj+1,1:nk)
                       ! Saturation mixing ratio

! Internal shared variables

      real rddvcp      ! rd / cp

      real p0iv        ! 1.0 / p0

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real t           ! Air temperature
      real p           ! Pressure

      real esw         ! Saturation vapor pressure for water

!-----7--------------------------------------------------------------7--

! Set the common used variables.

      rddvcp=rd/cp

      p0iv=1.e0/p0

! -----

! Calculate the saturation mixing ratio.

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j,p,t,esw)

        do j=1,nj-1
        do i=1,ni-1
          p=pbr(i,j,k)+pp(i,j,k)

          t=(ptbr(i,j,k)+ptp(i,j,k))*exp(rddvcp*log(p0iv*p))

          esw=es0*exp(17.269e0*(t-t0)/(t-35.86e0))

          qvs(i,j,k)=epsva*esw/(p-esw)

        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_getqvs

!-----7--------------------------------------------------------------7--

      end module m_getqvs
