!***********************************************************************
      module m_outstd14
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2008/10/10
!     Modification: 2009/02/27, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read in message to standard i/o.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: outstd14, s_outstd14

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outstd14

        module procedure s_outstd14

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
      subroutine s_outstd14(sname,ncsn,crsdir,datdir,nccrs,ncdat)
!***********************************************************************

! Input variables

      character(len=8), intent(in) :: sname
                       ! Called procedure name

      character(len=108), intent(in) :: crsdir
                       ! User specified directory for CReSS files

      character(len=108), intent(in) :: datdir
                       ! User specified directory for external data

      integer, intent(in) :: ncsn
                       ! Number of character of sname

      integer, intent(in) :: nccrs
                       ! Number of character of crsdir

      integer, intent(in) :: ncdat
                       ! Number of character of datdir

!-----7--------------------------------------------------------------7--

! Read in messages to standard i/o.

      write(6,*)

      write(6,'(a,a,a)') '  messages: procedure, ',sname(1:ncsn),';'

      if(nccrs.le.21.and.ncdat.le.24) then

        write(6,'(a,a,a)') '    Specified the directory, "',            &
     &     crsdir(1:nccrs),'" for CReSS processing files.'

        write(6,'(a,a,a)') '    Specified the directory, "',            &
     &     datdir(1:ncdat),'" for external data files.'

      else

        write(6,'(a,a,a)')                                              &
     &          '    Specified the directory, "',crsdir(1:nccrs),'"'

        write(6,'(a)')                                                  &
     &          '    for CReSS processing files.'

        write(6,*)

        write(6,'(a,a,a)')                                              &
     &          '    Specified the directory, "',datdir(1:ncdat),'"'

        write(6,'(a)')                                                  &
     &          '    for external data files.'

      end if

! -----

      end subroutine s_outstd14

!-----7--------------------------------------------------------------7--

      end module m_outstd14
