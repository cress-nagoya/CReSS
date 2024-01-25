!***********************************************************************
      module m_comuni
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/05/19
!     Modification: 2003/11/05, 2004/07/01, 2007/01/20, 2007/04/11,
!                   2008/05/02, 2008/08/25, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the array for unite.

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

      integer, allocatable, save :: iodmp(:)
                       ! Unit number of dumped file

      real, allocatable, save :: tmp1(:)
                       ! Temporary array

      real, allocatable, save :: tmp2(:)
                       ! Temporary array

      real, allocatable, save :: tmp3(:)
                       ! Temporary array

      real, allocatable, save :: tmp4(:)
                       ! Temporary array

      real, allocatable, save :: var(:,:)
                       ! Optional variable in dumped file

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

      end module m_comuni
