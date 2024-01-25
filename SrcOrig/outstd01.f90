!***********************************************************************
      module m_outstd01
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/04/06
!     Modification: 1999/06/14, 2000/01/17, 2001/05/29, 2003/04/30,
!                   2003/05/19, 2004/05/31, 2004/08/20, 2006/09/21,
!                   2007/01/20, 2008/05/02, 2008/08/25, 2008/10/10,
!                   2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read in the first message to standard i/o.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comstd

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: outstd01, s_outstd01

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outstd01

        module procedure s_outstd01

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
      subroutine s_outstd01(pname,ncpn)
!***********************************************************************

! Input variables

      character(len=8), intent(in) :: pname
                       ! Running program name

      integer, intent(in) :: ncpn
                       ! Number of character of pname

!-----7--------------------------------------------------------------7--

! Read in the first message to standard i/o.

      write(6,'(a,a,a)') 'Now the program, ',pname(1:ncpn),' start.'

! -----

! Initialize the module variable in the m_comstd.

      write(fstd(1:3),'(a3)') 'off'

! -----

      end subroutine s_outstd01

!-----7--------------------------------------------------------------7--

      end module m_outstd01
