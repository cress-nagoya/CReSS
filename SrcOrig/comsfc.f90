!***********************************************************************
      module m_comsfc
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/05/19
!     Modification: 2003/07/15, 2003/11/05, 2004/08/01, 2004/09/01,
!                   2005/01/14, 2007/01/20, 2008/05/02, 2008/08/25,
!                   2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the array for surface.

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

      integer, allocatable, save :: land(:,:)
                       ! Land use of surface

      integer, allocatable, save :: landat(:,:)
                       ! Land use in data

      real, allocatable, save :: ri(:,:)
                       ! Real indices in data region in x direction

      real, allocatable, save :: rj(:,:)
                       ! Real indices in data region in y direction

      real, allocatable, save :: albe(:,:)
                       ! Albedo

      real, allocatable, save :: beta(:,:)
                       ! Evapotranspiration efficiency

      real, allocatable, save :: z0m(:,:)
                       ! Roughness length for velocity

      real, allocatable, save :: z0h(:,:)
                       ! Roughness length for scalar

      real, allocatable, save :: cap(:,:)
                       ! Thermal capacity

      real, allocatable, save :: nuu(:,:)
                       ! Thermal diffusivity

      real, allocatable, save :: sst(:,:)
                       ! Sea surface temperature

      real, allocatable, save :: kai(:,:)
                       ! Sea ice distribution

      real, allocatable, save :: sstdat(:,:)
                       ! Sea surface temperature in data

      real, allocatable, save :: icedat(:,:)
                       ! Sea ice distribution in data

      real, allocatable, save :: tmp1(:,:)
                       ! Temporary array

      real, allocatable, save :: tmp2(:,:)
                       ! Temporary array

      real, allocatable, save :: tmp3(:,:)
                       ! Temporary array

      real, allocatable, save :: ltmp1(:,:)
                       ! Temporary array

      real, allocatable, save :: stmp1(:,:)
                       ! Temporary array

      real, allocatable, save :: itmp1(:,:)
                       ! Temporary array

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

      end module m_comsfc
