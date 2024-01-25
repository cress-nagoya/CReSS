!***********************************************************************
      module m_outstd06
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/03/25, 1999/04/06, 1999/05/10, 1999/06/14,
!                   1999/06/21, 1999/07/28, 1999/10/22, 2000/01/05,
!                   2000/01/17, 2002/04/02, 2002/07/15, 2003/05/19,
!                   2004/05/31, 2004/06/10, 2007/01/05, 2007/07/30,
!                   2008/05/02, 2008/08/25, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read in the message to standard i/o.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
      use m_outstd05

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: outstd06, s_outstd06

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outstd06

        module procedure s_outstd06

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_outstd06(nbstp0,ibstp,ctime)
!***********************************************************************

! Input variables

      integer(kind=i8), intent(in) :: nbstp0
                       ! Start index of large time steps

      integer(kind=i8), intent(in) :: ibstp
                       ! Index of large time steps

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

!-----7--------------------------------------------------------------7--

! Call another output procedure.

      if(ibstp.eq.nbstp0) then

        call outstd05(2)

      else

        call outstd05(1)

      end if

! -----

! Read in the messages to standard i/o.

      if(ibstp.eq.1_i8) then

        write(6,'(a,i8,a,i8,a,i3.3,a)') '  end of time step  = ',ibstp, &
     &                                  ', real time = ',ctime/1000_i8, &
     &                                  '.',mod(ctime,1000_i8),' [s]'

      else if(ibstp.ge.2) then

        write(6,'(a,i8,a,i8,a,i3.3,a)') '  end of time steps = ',ibstp, &
     &                                  ', real time = ',ctime/1000_i8, &
     &                                  '.',mod(ctime,1000_i8),' [s]'

      end if

! -----

      end subroutine s_outstd06

!-----7--------------------------------------------------------------7--

      end module m_outstd06
