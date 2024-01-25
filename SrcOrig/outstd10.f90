!***********************************************************************
      module m_outstd10
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/07/03
!     Modification: 2003/03/28, 2003/04/30, 2003/05/19, 2004/04/15,
!                   2004/06/10, 2006/09/21, 2007/01/05, 2007/04/11,
!                   2008/05/02, 2008/08/25, 2008/10/10, 2009/02/27,
!                   2011/09/22

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

      public :: outstd10, s_outstd10

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outstd10

        module procedure s_outstd10

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
      subroutine s_outstd10(vname,ncvn,vcap,ncvc,cntgeo)
!***********************************************************************

! Input variable

      character(len=6), intent(in) :: vname
                       ! Optional variable name

      character(len=60), intent(in) :: vcap
                       ! Caption for dumped variable

      integer, intent(in) :: ncvn
                       ! Number of character of vname

      integer, intent(in) :: ncvc
                       ! Number of character of vcap

      integer, intent(in) :: cntgeo
                       ! Dumped variables counter

!-----7--------------------------------------------------------------7--

! Read in the messages to standard i/o.

      if(cntgeo.eq.1) then

        write(6,*)

        write(6,'(a)') '  i/o: procedure, outgeo;'

      end if

      write(6,'(a,a,a,a,a)') '    Dumped the 2d variable, ',            &
     &          vcap(1:ncvc),', ',vname(1:ncvn),'.'

! -----

      end subroutine s_outstd10

!-----7--------------------------------------------------------------7--

      end module m_outstd10
