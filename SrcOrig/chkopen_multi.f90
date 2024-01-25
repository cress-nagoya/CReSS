!***********************************************************************
      module m_chkopen
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/06/14
!     Modification: 1999/06/21, 2000/01/17, 2001/05/29, 2001/11/20,
!                   2002/06/18, 2003/05/19, 2003/11/05, 2004/08/31,
!                   2004/09/25, 2005/02/10, 2006/12/04, 2007/01/20,
!                   2008/05/02, 2008/07/25, 2008/08/25, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     check the common files open between all processor elements.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_commpi
      use m_cpondpe
      use m_defmpi
      use m_destroy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: chkopen, s_chkopen

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface chkopen

        module procedure s_chkopen

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
      subroutine s_chkopen(extpe,ionum)
!***********************************************************************

! Input variable

      integer, intent(in) :: extpe
                       ! Control flag of file name extension

      integer, intent(in) :: ionum
                       ! Opened file unit number

! Internal shared variables

      integer extcom   ! Control flag of common file name extension

      integer stat     ! Runtime status

      integer ierr     ! Error descriptor

!-----7--------------------------------------------------------------7--

!! Check the common files opened.

! Set the reference variable.

      extcom=extpe

! -----

! Broadcast the control flag.

      call mpi_bcast(extcom,1,mpi_integer,root,mpi_comm_cress,ierr)

! -----

! If error occured, call the procedure destroy.

      stat=extpe-extcom

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('chkopen ',7,'cont',1,'              ',14,ionum, &
     &                 stat)

        end if

        call cpondpe

        call destroy('chkopen ',7,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

!! -----

      end subroutine s_chkopen

!-----7--------------------------------------------------------------7--

      end module m_chkopen
