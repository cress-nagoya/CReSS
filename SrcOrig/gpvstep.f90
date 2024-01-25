!***********************************************************************
      module m_gpvstep
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/09/09
!     Modification: 2003/05/19, 2004/05/07, 2004/08/20, 2004/09/01,
!                   2005/02/10, 2006/11/27, 2007/01/20, 2007/03/10,
!                   2007/05/21, 2007/07/30, 2008/05/02, 2008/08/25,
!                   2008/12/11, 2009/01/05, 2009/02/27, 2011/09/22,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the number of steps of the main do loop in gpvdrv.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_comkind
      use m_commpi
      use m_cpondpe
      use m_destroy
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: gpvstep, s_gpvstep

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface gpvstep

        module procedure s_gpvstep

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
      subroutine s_gpvstep(fpnggopt,fpexbopt,fplspopt,fpvspopt,fpgpvitv,&
     &                     fpstime,fpetime,nstp0,nstp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpnggopt
                       ! Formal parameter of unique index of nggopt

      integer, intent(in) :: fpexbopt
                       ! Formal parameter of unique index of exbopt

      integer, intent(in) :: fplspopt
                       ! Formal parameter of unique index of lspopt

      integer, intent(in) :: fpvspopt
                       ! Formal parameter of unique index of vspopt

      integer, intent(in) :: fpgpvitv
                       ! Formal parameter of unique index of gpvitv

      integer, intent(in) :: fpstime
                       ! Formal parameter of unique index of stime

      integer, intent(in) :: fpetime
                       ! Formal parameter of unique index of etime

! Output variables

      integer(kind=i8), intent(out) :: nstp0
                       ! Start index of main do loop

      integer(kind=i8), intent(out) :: nstp1
                       ! End index of main do loop

! Internal shared variables

      integer nggopt   ! Option for analysis nudging to GPV
      integer exbopt   ! Option for external boundary forcing
      integer lspopt   ! Option for lateral sponge damping
      integer vspopt   ! Option for vertical sponge damping

      integer(kind=i8) gpv103
                       ! 1000 x int(gpvitv + 0.1)

      integer stat     ! Runtime status

      real stime       ! Forecast start time
      real etime       ! Forecast end time

      real gpvitv      ! Time interval of GPV data file

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpnggopt,nggopt)
      call getiname(fpexbopt,exbopt)
      call getiname(fplspopt,lspopt)
      call getiname(fpvspopt,vspopt)
      call getrname(fpgpvitv,gpvitv)
      call getrname(fpstime,stime)
      call getrname(fpetime,etime)

! -----

! Calculate the number of steps of the main do loop in gpvdrv.

      if(nggopt.eq.1.or.exbopt.ge.1                                     &
     &  .or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

        gpv103=1000_i8*int(gpvitv+.1e0,i8)

        nstp0=1000_i8*int(stime+.1e0,i8)/gpv103+1_i8

        if(mod(1000_i8*int(etime+.1e0,i8),gpv103).eq.0_i8) then

          nstp1=1000_i8*int(etime+.1e0,i8)/gpv103+1_i8

        else

          nstp1=1000_i8*int(etime+.1e0,i8)/gpv103+2_i8

        end if

        stat=0

        if(nstp0.gt.nstp1) then
          stat=1
        end if

        call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('gpvstep ',7,'cont',7,'              ',14,101, &
     &                   stat)

          end if

          call cpondpe

          call destroy('gpvstep ',7,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

      else

        nstp0=1_i8
        nstp1=1_i8

      end if

! -----

      end subroutine s_gpvstep

!-----7--------------------------------------------------------------7--

      end module m_gpvstep
