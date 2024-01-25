!***********************************************************************
      module m_commath
!***********************************************************************

!     Author      : Sakakibara Atsushi, Naito Daisuke
!     Date        : 2003/05/19
!     Modification: 2004/03/05, 2004/04/15, 2004/05/07, 2004/05/31,
!                   2004/08/01, 2004/08/20, 2004/09/01, 2004/09/10,
!                   2004/10/12, 2006/02/13, 2008/05/02, 2008/08/25,
!                   2008/10/10, 2009/02/27, 2011/03/18, 2011/03/29,
!                   2011/09/22

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the mathmatical constants.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      public

! Exceptional access control

!     none

!-----7--------------------------------------------------------------7--

! Module variables

      real, parameter :: cc=3.141592e0
                       ! Circular constant

      real, parameter :: d2r=3.141592e0/180.e0
                       ! Converter of degree to radian

      real, parameter :: r2d=180.e0/3.141592e0
                       ! Converter of radian to degree

      real, parameter :: oned3=1.e0/3.e0
                       ! 1.0 / 3.0

      real, parameter :: oned6=1.e0/6.e0
                       ! 1.0 / 6.0

      real, parameter :: oned9=1.e0/9.e0
                       ! 1.0 / 9.0

      real, parameter :: oned15=1.e0/15.e0
                       ! 1.0 / 15.0

      real, parameter :: oned24=1.e0/24.e0
                       ! 1.0 / 24.0

      real, parameter :: oned27=1.e0/27.e0
                       ! 1.0 / 27.0

      real, parameter :: oned60=1.e0/60.e0
                       ! 1.0 / 60.0

      real, parameter :: twod15=2.e0/15.e0
                       ! 2.0 / 15.0

      real, parameter :: fourd3=4.e0/3.e0
                       ! 4.0 / 3.0

      real, parameter :: sevnd3=7.e0/3.e0
                       ! 7.0 / 3.0

      real, parameter :: tend3=10.e0/3.e0
                       ! 10.0 / 3.0

      real, parameter :: tend6=10.e0/6.e0
                       ! 10.0 / 6.0

      real, parameter :: tend7=10.e0/7.e0
                       ! 10.0 / 7.0

      real, parameter :: i365=1.e0/365.e0
                       ! 1.0 / 365.0

      real, parameter :: i366=1.e0/366.e0
                       ! 1.0 / 366.0

      real, parameter :: i65536=1.e0/65536.e0
                       ! 1.0 / 65536.0

      real, parameter :: ln10=2.30259e0
                       ! ln(10.0)

      real, parameter :: tanh2=.96403e0
                       ! tanh(2.0)

      real, parameter :: gf4=6.e0
                       ! Gamma function at 4.0

      real, parameter :: gf7=720.e0
                       ! Gamma function at 7.0

      real, parameter :: gfbus=1.77245e0
                       ! Gamma function at bus

      real, parameter :: gf1buc=2.e0
                       ! Gamma function at 1 + buc

      real, parameter :: gf1bur=.931384e0
                       ! Gamma function at 1 + bur

      real, parameter :: gf1bui=1.e0
                       ! Gamma function at 1 + bui

      real, parameter :: gf1bus=.886227e0
                       ! Gamma function at 1 + bus

      real, parameter :: gf1bug=.898642e0
                       ! Gamma function at 1 + bug

      real, parameter :: gf1buh=.898642e0
                       ! Gamma function at 1 + buh

      real, parameter :: gf3bur=4.69417e0
                       ! Gamma function at 3 + bur

      real, parameter :: gf3bus=3.32335e0
                       ! Gamma function at 3 + bus

      real, parameter :: gf3bug=3.89076e0
                       ! Gamma function at 3 + bug

      real, parameter :: gf3buh=3.89076e0
                       ! Gamma function at 3 + buh

      real, parameter :: gf4buc=120.e0
                       ! Gamma function at 4 + buc

      real, parameter :: gf4bur=17.8379e0
                       ! Gamma function at 4 + bur

      real, parameter :: gf4bui=24.e0
                       ! Gamma function at 4 + bui

      real, parameter :: gf4bus=11.6314e0
                       ! Gamma function at 4 + bus

      real, parameter :: gf4bug=14.1624e0
                       ! Gamma function at 4 + bug

      real, parameter :: gf4buh=14.1624e0
                       ! Gamma function at 4 + buh

      real, parameter :: gf5bur=1.82736e0
                       ! Gamma function at (5 + bur) / 2

      real, parameter :: gf5bus=1.60836e0
                       ! Gamma function at (5 + bus) / 2

      real, parameter :: gf5bug=1.70506e0
                       ! Gamma function at (5 + bug) / 2

      real, parameter :: gf5buh=1.70506e0
                       ! Gamma function at (5 + buh) / 2

      real, parameter :: hfbus=1610.e0
                       ! Hypergeometric function at bus

      real, parameter :: eps=1.e-20
                       ! Very small constant 1.0 x 10^-20

      real, parameter :: epsn=-1.e-20
                       ! Very small negative constant - 1.0 x 10^-20

      real, parameter :: lim36=1.e36
                       ! Very large constant 1.0 x 10^36

      real, parameter :: lim34n=-1.e34
                       ! Very large negative constant - 1.0 x 10^34

      real, parameter :: lim35n=-1.e35
                       ! Very large negative constant - 1.0 x 10^35

      real, parameter :: lim36n=-1.e36
                       ! Very large negative constant - 1.0 x 10^36

! Module procedure

!     none

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

!     none

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

!     none

!-----7--------------------------------------------------------------7--

      end module m_commath
