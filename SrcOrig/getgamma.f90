!***********************************************************************
      module m_getgamma
!***********************************************************************

!     Author      : Naito Daisuke, Sakakibara Atsushi
!     Date        : 2011/03/24
!     Modification: 2011/03/29

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the value of Gamma function at specified value of
!     independent variable.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
      use m_commath

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getgamma, s_getgamma

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getgamma

        module procedure s_getgamma

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic sqrt
      intrinsic exp
      intrinsic log

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_getgamma(x,gamma)
!***********************************************************************

! Input variables

      real, intent(in) :: x
                       ! Value of independent variable

! Output variables

      real, intent(out) :: gamma
                       ! Answer value of Gamma function

! Internal shared variables

      real(kind=r8) c0 ! Coefficient of approximation formula,
                       !  1.000000000178

      real(kind=r8) c1 ! Coefficient of approximation formula,
                       !  76.180091729406

      real(kind=r8) c2 ! Coefficient of approximation formula,
                       ! -86.505320327112

      real(kind=r8) c3 ! Coefficient of approximation formula,
                       !  24.014098222230

      real(kind=r8) c4 ! Coefficient of approximation formula,
                       !  -1.231739516140

      real(kind=r8) c5 ! Coefficient of approximation formula,
                       !   0.001208580030

      real(kind=r8) c6 ! Coefficient of approximation formula,
                       !  -0.000005363820

      real(kind=r8) a  ! Temporary variable

!-----7--------------------------------------------------------------7--

! Set the coefficients of approximation formula.

      c0=1.000000000178e0
      c1=76.180091729406e0
      c2=-86.505320327112e0
      c3=24.014098222230e0
      c4=-1.231739516140e0
      c5=0.001208580030e0
      c6=-0.000005363820e0

! -----

! Get the value of Gamma function at specified value of independent
! variable.

      a=c0+c1/real(x,r8)+c2/(real(x,r8)+1.e0_r8)                        &
     &    +c3/(real(x,r8)+2.e0_r8)+c4/(real(x,r8)+3.e0_r8)              &
     &    +c5/(real(x,r8)+4.e0_r8)+c6/(real(x,r8)+5.e0_r8)

      gamma=sqrt(2.e0*cc)*exp(-real(x,r8)-4.5e0_r8)*a                   &
     &  *exp((real(x,r8)-.5e0_r8)*log(real(x,r8)+4.5e0_r8))

! -----

      end subroutine s_getgamma

!-----7--------------------------------------------------------------7--

      end module m_getgamma
