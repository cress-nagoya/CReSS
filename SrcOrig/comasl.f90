!***********************************************************************
      module m_comasl
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/09/22
!     Modification: 2011/11/10

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the array for asldata.

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

      real, allocatable, save :: z(:)
                       ! zeta coordinates

      real, allocatable, save :: zph(:,:,:)
                       ! z physical coordinates

      real, allocatable, save :: qasl(:,:,:,:)
                       ! Aerosol mixing ratio

      real, allocatable, save :: tmp1(:,:)
                       ! Temporary array

      real, allocatable, save :: tmp2(:,:)
                       ! Temporary array

      real, allocatable, save :: tmp3(:,:)
                       ! Temporary array

      real, allocatable, save :: tmp4(:,:)
                       ! Temporary array

      real, allocatable, save :: tmp5(:,:,:)
                       ! Temporary array

      real, allocatable, save :: tmp6(:,:,:)
                       ! Temporary array

      real, allocatable, save :: zdat(:,:,:)
                       ! z physical coordinates in data

      real, allocatable, save :: qadat(:,:,:,:)
                       ! Aerosol mixing ratio in data

      real, allocatable, save :: dtmp1(:,:,:)
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

      end module m_comasl
