!***********************************************************************
      module m_chkkind
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/07/30
!     Modification: 2007/09/25, 2008/04/17, 2008/05/02, 2008/07/25,
!                   2008/08/25, 2009/01/05, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     check the kind of the Fortran variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_comkind
      use m_commpi
      use m_cpondpe
      use m_destroy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: chkkind, s_chkkind

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface chkkind

        module procedure s_chkkind

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
      subroutine s_chkkind
!***********************************************************************

! Internal shared variables

      integer stat     ! Runtime status

      integer i4rng    ! Exponent range of
                       ! single precision integer variable

!WHY? integer i8rng    ! Exponent range of
!WHY?                  ! double precision integer variable

      integer r4rng    ! Exponent range of
                       ! single precision real variable

!-----7--------------------------------------------------------------7--

! Initialize the runtime status.

      stat=0

! -----

! Check the kind of the Fortran single precision variables.

      i4rng=range(0)
      r4rng=range(0.e0)

      if(i4rng.lt.9) then
        stat=stat+1
      end if

      if(r4rng.lt.37) then
        stat=stat+1
      end if

! -----

! Check the kind of the Fortran double precision integer variable.

!WHY? if(i8.lt.0) then
!WHY?   stat=stat+1
!WHY? end if

!WHY? if(stat.eq.0) then

!WHY?   i8rng=range(0_i8)

!WHY?   if(i8rng.lt.2*i4rng) then
!WHY?     stat=stat+1
!WHY?   end if

!WHY? end if

! -----

! If error occured, call the procedure destroy.

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('chkkind ',7,'cont',7,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('chkkind ',7,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

      end subroutine s_chkkind

!-----7--------------------------------------------------------------7--

      end module m_chkkind
