!***********************************************************************
      module m_comflt
!***********************************************************************

!     Author      : Satoki Tsujino
!     Date        : 2017/05/30
!     Modification: 2017/06/11

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the array for Asselin time filter procedure.

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

      real, allocatable, dimension(:,:,:), save :: dmpfltu
                       ! Filtering term in u equation

      real, allocatable, dimension(:,:,:), save :: dmpfltv
                       ! Filtering term in v equation

      real, allocatable, dimension(:,:,:), save :: dmpfltw
                       ! Filtering term in w equation

      real, allocatable, dimension(:,:,:), save :: dmpfltp
                       ! Filtering term in p equation

      real, allocatable, dimension(:,:,:), save :: dmpfltpt
                       ! Filtering term in pt equation

      real, allocatable, dimension(:,:,:), save :: dmpfltqv
                       ! Filtering term in qv equation

      real, allocatable, dimension(:,:,:), save :: netdmpu
                       ! net tendency in u equation

      real, allocatable, dimension(:,:,:), save :: netdmpv
                       ! net tendency in v equation

      real, allocatable, dimension(:,:,:), save :: netdmpw
                       ! net tendency in w equation

      real, allocatable, dimension(:,:,:), save :: netdmpp
                       ! net tendency in p equation

      real, allocatable, dimension(:,:,:), save :: netdmppt
                       ! net tendency in pt equation

      real, allocatable, dimension(:,:,:), save :: netdmpqv
                       ! net tendency in qv equation

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

      end module m_comflt
