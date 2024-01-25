!***********************************************************************
      module m_unistep
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/09/09
!     Modification: 2003/04/30, 2003/05/19, 2004/08/01, 2004/08/20,
!                   2004/09/01, 2005/01/14, 2005/02/10, 2006/09/21,
!                   2006/01/20, 2007/07/30, 2007/09/28, 2008/04/17,
!                   2008/05/02, 2008/08/25, 2009/01/05, 2009/02/27,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the number of steps of the main do loop in unidrv.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
      use m_destroy
      use m_getcname
      use m_getrname
      use m_inichar

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: unistep, s_unistep

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface unistep

        module procedure s_unistep

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
      subroutine s_unistep(fpfltyp_uni,fpflitv_uni,fpstime,fpetime,     &
     &                     nstp0,nstp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpfltyp_uni
                       ! Formal parameter of unique index of fltyp_uni

      integer, intent(in) :: fpflitv_uni
                       ! Formal parameter of unique index of flitv_uni

      integer, intent(in) :: fpstime
                       ! Formal parameter of unique index of stime

      integer, intent(in) :: fpetime
                       ! Formal parameter of unique index of etime

! Output variables

      integer(kind=i8), intent(out) :: nstp0
                       ! Start index of main do loop

      integer(kind=i8), intent(out) :: nstp1
                       ! End index of main do loop

! Internal shared variable

      character(len=108) fltyp_uni
                       ! Control flag of processed file type

      integer(kind=i8) fl103
                       ! 1000 x int(flitv_uni + 0.1)

      integer stat     ! Runtime status

      real stime       ! Forecast start time
      real etime       ! Forecast end time

      real flitv_uni   ! Time interval of processed file

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(fltyp_uni)

! -----

! Get the required namelist variables.

      call getcname(fpfltyp_uni,fltyp_uni)
      call getrname(fpflitv_uni,flitv_uni)
      call getrname(fpstime,stime)
      call getrname(fpetime,etime)

! -----

! Calculate the number of steps of the main do loop.

      if(fltyp_uni(1:3).eq.'all'                                        &
     &  .or.fltyp_uni(1:3).eq.'dmp'.or.fltyp_uni(1:3).eq.'mon') then

        fl103=1000_i8*int(flitv_uni+.1e0,i8)

        if(mod(1000_i8*int(stime+.1e0,i8),fl103).eq.0_i8) then

          nstp0=1000_i8*int(stime+.1e0,i8)/fl103+1_i8

        else

          nstp0=1000_i8*int(stime+.1e0,i8)/fl103+2_i8

        end if

        nstp1=1000_i8*int(etime+.1e0,i8)/fl103+1_i8

        if(nstp0.gt.nstp1) then

          call destroy('unistep ',7,'stop',7,'              ',14,101,   &
     &                 stat)

        end if

      else if(fltyp_uni(1:3).eq.'geo') then

        nstp0=1_i8
        nstp1=1_i8

      end if

! -----

      end subroutine s_unistep

!-----7--------------------------------------------------------------7--

      end module m_unistep
