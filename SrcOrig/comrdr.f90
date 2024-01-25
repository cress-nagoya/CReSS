!***********************************************************************
      module m_comrdr
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/05/19
!     Modification: 2003/11/05, 2003/12/12, 2007/01/20, 2008/05/02,
!                   2008/08/25, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the array for radata.

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

      real, allocatable, save :: u(:,:,:)
                       ! x components of velocity

      real, allocatable, save :: v(:,:,:)
                       ! y components of velocity

      real, allocatable, save :: w(:,:,:)
                       ! z components of velocity

      real, allocatable, save :: qp(:,:,:)
                       ! Precipitation mixing ratio

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

      real, allocatable, save :: londat(:,:)
                       ! Longitude in data

      real, allocatable, save :: zdat(:,:,:)
                       ! z physical coordinates in data

      real, allocatable, save :: udat(:,:,:)
                       ! x components of velocity in data

      real, allocatable, save :: vdat(:,:,:)
                       ! y components of velocity in data

      real, allocatable, save :: wdat(:,:,:)
                       ! z components of velocity in data

      real, allocatable, save :: qpdat(:,:,:)
                       ! Precipitation mixing ratio in data

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

      end module m_comrdr
