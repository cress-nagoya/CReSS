!***********************************************************************
      module m_gettime
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/03/25, 1999/05/10, 1999/06/07, 1999/07/28,
!                   2000/01/05, 2000/01/17, 2001/03/13, 2001/06/29,
!                   2002/04/02, 2002/06/18, 2002/09/09, 2002/10/31,
!                   2002/12/02, 2003/01/20, 2003/05/19, 2004/05/07,
!                   2004/08/20, 2004/09/01, 2006/04/03, 2007/07/30,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the model forecast time in the large time steps
!     integration.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
      use m_commath
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: gettime, s_gettime

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface gettime

        module procedure s_gettime

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic int

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_gettime(fpadvopt,fpdtbig,                            &
     &                     ibstp,ctime,ptime,ftime,pmin)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fpdtbig
                       ! Formal parameter of unique index of dtbig

      integer(kind=i8), intent(in) :: ibstp
                       ! Index of large time steps

! Output variables

      integer(kind=i8), intent(out) :: ctime
                       ! Model current forecast time

      integer(kind=i8), intent(out) :: ptime
                       ! Model forecast time at 1 step past

      integer(kind=i8), intent(out) :: ftime
                       ! Model forecast time at 1 step future

      integer(kind=i8), intent(out) :: pmin
                       ! 60000 x (ptime / 60000)

! Internal shared variables

      integer advopt   ! Option for advection scheme

      real dtbig       ! Large time steps inteval

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpadvopt,advopt)
      call getrname(fpdtbig,dtbig)

! -----

! Calculate the model forecast time in the large time steps.

      if(advopt.le.3) then

        ctime=10_i8*int(1.e2*(dtbig+.001e0),i8)*(ibstp-1_i8)
        ftime=10_i8*int(1.e2*(dtbig+.001e0),i8)*ibstp

        if(ctime.eq.0_i8) then
          ptime=ctime
        else
          ptime=ctime-10_i8*int(1.e2*(dtbig+.001e0),i8)
        end if

        pmin=60000_i8*(ptime/60000_i8)

      else

        ptime=10_i8*int(1.e2*(dtbig+.001e0),i8)*(ibstp-1_i8)
        ctime=10_i8*int(1.e2*(dtbig+.001e0),i8)*(ibstp-1_i8)
        ftime=10_i8*int(1.e2*(dtbig+.001e0),i8)*ibstp

        pmin=60000_i8*(ptime/60000_i8)

      end if

! -----

      end subroutine s_gettime

!-----7--------------------------------------------------------------7--

      end module m_gettime
