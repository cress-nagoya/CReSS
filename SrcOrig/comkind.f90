!***********************************************************************
      module m_comkind
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/07/30
!     Modification: 2008/04/17, 2008/05/02, 2008/07/25, 2008/08/25,
!                   2008/10/10, 2009/01/30, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the kind of Fortran double precision integer variable.

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

!WHY? integer, parameter :: i8=selected_int_kind(2*range(0))
      integer, parameter :: i8=selected_int_kind(12)
                       ! Kind of double precision integer variable

!WHY? integer, parameter :: r8=selected_real_kind(2*precision(0.e0))
      integer, parameter :: r8=selected_real_kind(15)
                       ! Kind of double precision real variable

! Module procedure

!     none

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

!WHY? intrinsic range
!WHY? intrinsic precision
      intrinsic selected_int_kind
      intrinsic selected_real_kind

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

!     none

!-----7--------------------------------------------------------------7--

      end module m_comkind
