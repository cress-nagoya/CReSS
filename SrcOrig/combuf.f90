!***********************************************************************
      module m_combuf
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/05/19
!     Modification: 2003/11/05, 2004/08/31, 2006/12/04, 2007/01/20,
!                   2008/05/02, 2008/08/25, 2008/10/10, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the communication buffer.

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

      integer, allocatable, save :: idxbuf(:,:)
                       ! Index of maximum or minimum value
                       ! in each processor element

      real, allocatable, save :: mxnbuf(:)
                       ! Maximum or minimum value
                       ! in each processor element

      real, allocatable, save :: sbuf(:)
                       ! Common used sending buffer

      real, allocatable, save :: rbuf(:)
                       ! Common used receiving buffer

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

      end module m_combuf
