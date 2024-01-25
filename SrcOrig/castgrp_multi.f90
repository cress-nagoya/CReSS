!***********************************************************************
      module m_castgrp
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/12/04
!     Modification: 2007/01/05, 2007/01/20, 2007/08/24, 2008/05/02,
!                   2008/07/25, 2008/08/25, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     broadcast the group domain arrangemnet table from the root
!     processor element to the others.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comgrp
      use m_commpi
      use m_defmpi

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: castgrp, s_castgrp

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface castgrp

        module procedure s_castgrp

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

!     none

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_castgrp
!***********************************************************************

! Internal shared variables

      integer siz      ! Broadcasting buffer size

      integer ierr     ! Error descriptor

!-----7--------------------------------------------------------------7--

! Broadcast the group domain arrangement variable.

      siz=ngrp

      call mpi_bcast(chrgrp,siz,mpi_character,root,mpi_comm_cress,ierr)

! -----

      end subroutine s_castgrp

!-----7--------------------------------------------------------------7--

      end module m_castgrp
