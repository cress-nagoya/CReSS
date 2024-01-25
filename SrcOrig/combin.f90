!***********************************************************************
      module m_combin
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/09/30
!     Modification: 2007/01/20, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2009/03/12

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the array for warm bin cloud physics.

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

      real, allocatable, save :: rbw(:)
                       ! Standard radius
                       ! between adjacent water bins [cm]

      real, allocatable, save :: rrbw(:,:)
                       ! Related parameters of rbw

      real, allocatable, save :: brw(:)
                       ! Radius at water bin boundary [cm]

      real, allocatable, save :: rbrw(:,:)
                       ! Related parameters of brw

      real, allocatable, save :: bmw(:,:)
                       ! Mass at water bin boundary [g]

      real, allocatable, save :: rbmw(:,:)
                       ! Related parameters of bmw

      real, allocatable, save :: dbmw(:)
                       ! Differential between adjacent water bins [g]

      real, allocatable, save :: ewbw(:,:)
                       ! Radius weighted coalescence efficiency
                       ! between water bins [cm^2]

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

      end module m_combin
