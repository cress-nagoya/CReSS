!***********************************************************************
      module m_comgrp
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/12/04
!     Modification: 2007/01/05, 2007/01/20, 2008/05/02, 2008/08/25,
!                   2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the group domain arrangement table in entire domain.

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

      character(len=1), allocatable, save :: chrgrp(:,:)
                       ! User specified group domain arrangement
                       ! in entire domain

      integer, allocatable, save :: grpxy(:,:)
                       ! Adjusted group domain arrangement
                       ! in entire domain

      integer, allocatable, save :: xgrp(:)
                       ! Index in grpxy in x direction

      integer, allocatable, save :: ygrp(:)
                       ! Index in grpxy in x direction

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

      end module m_comgrp
