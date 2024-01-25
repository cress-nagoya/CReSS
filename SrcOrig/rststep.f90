!***********************************************************************
      module m_rststep
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2005/07/01
!     Modification: 2005/08/05, 2007/01/20, 2007/07/30, 2008/05/02,
!                   2008/08/25, 2009/01/05, 2009/02/27, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the number of steps of the main do loop in rstdrv.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
      use m_destroy
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: rststep, s_rststep

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rststep

        module procedure s_rststep

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic int
      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_rststep(fpflitv_rst,fpstime,fpetime,nstp0,nstp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpflitv_rst
                       ! Formal  parameter of unique index of flitv_rst

      integer, intent(in) :: fpstime
                       ! Formal  parameter of unique index of stime

      integer, intent(in) :: fpetime
                       ! Formal  parameter of unique index of etime

! Output variables

      integer(kind=i8), intent(out) :: nstp0
                       ! Start index of main do loop

      integer(kind=i8), intent(out) :: nstp1
                       ! End index of main do loop

! Internal shared variable

      integer(kind=i8) fl103
                       ! 1000 x int(flitv_rst + 0.1)

      integer stat     ! Runtime status

      real stime       ! Forecast start time
      real etime       ! Forecast end time

      real flitv_rst   ! Time interval of processed file

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getrname(fpflitv_rst,flitv_rst)
      call getrname(fpstime,stime)
      call getrname(fpetime,etime)

! -----

! Calculate the number of steps of the main do loop.

      fl103=1000_i8*int(flitv_rst+.1e0,i8)

      if(mod(1000_i8*int(stime+.1e0,i8),fl103).eq.0_i8) then

        nstp0=1000_i8*int(stime+.1e0,i8)/fl103+1_i8

      else

        nstp0=1000_i8*int(stime+.1e0,i8)/fl103+2_i8

      end if

      nstp1=1000_i8*int(etime+.1e0,i8)/fl103+1_i8

      if(nstp0.gt.nstp1) then

        call destroy('rststep ',7,'stop',7,'              ',14,101,stat)

      end if

! -----

      end subroutine s_rststep

!-----7--------------------------------------------------------------7--

      end module m_rststep
