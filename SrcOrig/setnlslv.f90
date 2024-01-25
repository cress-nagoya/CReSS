!***********************************************************************
      module m_setnlslv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2008/04/17
!     Modification: 2008/05/02, 2008/08/25, 2008/10/10, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read out and set the namelist variables for solver.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_chknlslv
      use m_chkstd
      use m_comkind
      use m_commpi
      use m_destroy
      use m_outstd04
      use m_outstd13
      use m_rdconf
      use m_setname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: setnlslv, s_setnlslv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface setnlslv

        module procedure s_setnlslv

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
      subroutine s_setnlslv
!***********************************************************************

! Internal shared variable

      integer stat     ! Runtime status

!-----7--------------------------------------------------------------7--

!! Read out and check the namelist variables.

! Read in the message to the standard i/o.

      if(mype.eq.root) then

        call outstd04(0,0_i8)

      end if

      call chkstd(root)

! -----

! Read out the namelist variables.

      stat=0

      if(mype.eq.root) then

        call rdconf(stat)

      end if

! -----

! If error occured, call the procedure destroy.

      call chkerr(stat)

      if(stat.lt.0) then

        call destroy('rdconf  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

! Read in the messages to the standard i/o.

      if(mype.eq.root) then

        call outstd13('setnlslv',8)

      end if

      call chkstd(root)

! -----

! Check the namelist variables.

      stat=0

      if(mype.eq.root) then

        call chknlslv('solver  ',6,stat)

      end if

! -----

! If error occured, call the procedure destroy.

      call chkerr(stat)

      if(stat.lt.0) then

        call destroy('setnlslv',8,'stop',1001,'              ',14,5,    &
     &               stat)

      end if

! -----

!! -----

! Set the table to archive namelist variables.

      call setname

! -----

      end subroutine s_setnlslv

!-----7--------------------------------------------------------------7--

      end module m_setnlslv
