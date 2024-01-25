!***********************************************************************
      module m_outstd09
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/01/15
!     Modification: 2001/03/13, 2001/04/15, 2001/06/29, 2001/12/10,
!                   2003/03/28, 2003/04/30, 2003/05/19, 2004/06/10,
!                   2006/09/21, 2007/01/05, 2008/04/17, 2008/05/02,
!                   2008/08/25, 2008/10/10, 2009/02/27, 2011/09/22

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read in the messages to standard i/o.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: outstd09, s_outstd09

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outstd09

        module procedure s_outstd09

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
      subroutine s_outstd09(vname,ncvn,vcap,ncvc,ndim,cntdmp)
!***********************************************************************

! Input variables

      character(len=6), intent(in) :: vname
                       ! Optional variable name

      character(len=60), intent(in) :: vcap
                       ! Caption for dumped variable

      integer, intent(in) :: ncvn
                       ! Number of character of vname

      integer, intent(in) :: ncvc
                       ! Number of character of vcap

      integer, intent(in) :: ndim
                       ! Input variable dimension

      integer, intent(in) :: cntdmp
                       ! Counter of dumped variables

!-----7--------------------------------------------------------------7--

! Read in the messages to standard i/o.

      if(cntdmp.eq.1) then

        write(6,*)

        write(6,'(a)') '  i/o: procedure, outdmp[2d|3d];'

      end if

      if(ndim.eq.2) then

        write(6,'(a,a,a,a,a)') '    Dumped the 2d variable, ',          &
     &            vcap(1:ncvc),', ',vname(1:ncvn),'.'

      else if(ndim.eq.3) then

        write(6,'(a,a,a,a,a)') '    Dumped the 3d variable, ',          &
     &            vcap(1:ncvc),', ',vname(1:ncvn),'.'

      end if

! -----

      end subroutine s_outstd09

!-----7--------------------------------------------------------------7--

      end module m_outstd09
