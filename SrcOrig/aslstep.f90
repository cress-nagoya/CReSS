!***********************************************************************
      module m_aslstep
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/09/22
!     Modification: 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the number of steps of the main do loop in asldrv.

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

      public :: aslstep, s_aslstep

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface aslstep

        module procedure s_aslstep

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
      subroutine s_aslstep(fpnggopt,fpexbopt,fplspopt,fpvspopt,fpaslitv,&
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

      integer, intent(in) :: fpaslitv
                       ! Formal parameter of unique index of aslitv

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

      integer(kind=i8) asl103
                       ! 1000 x int(aslitv + 0.1)

      integer stat     ! Runtime status

      real stime       ! Forecast start time
      real etime       ! Forecast end time

      real aslitv      ! Time interval of aerosol data file

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpnggopt,nggopt)
      call getiname(fpexbopt,exbopt)
      call getiname(fplspopt,lspopt)
      call getiname(fpvspopt,vspopt)
      call getrname(fpaslitv,aslitv)
      call getrname(fpstime,stime)
      call getrname(fpetime,etime)

! -----

! Calculate the number of steps of the main do loop in asldrv.

      if(nggopt.eq.1.or.exbopt.ge.1                                     &
     &  .or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

        asl103=1000_i8*int(aslitv+.1e0,i8)

        nstp0=1000_i8*int(stime+.1e0,i8)/asl103+1_i8

        if(mod(1000_i8*int(etime+.1e0,i8),asl103).eq.0_i8) then

          nstp1=1000_i8*int(etime+.1e0,i8)/asl103+1_i8

        else

          nstp1=1000_i8*int(etime+.1e0,i8)/asl103+2_i8

        end if

        stat=0

        if(nstp0.gt.nstp1) then
          stat=1
        end if

        call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('aslstep ',7,'cont',7,'              ',14,101, &
     &                   stat)

          end if

          call cpondpe

          call destroy('aslstep ',7,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

      else

        nstp0=1_i8
        nstp1=1_i8

      end if

! -----

      end subroutine s_aslstep

!-----7--------------------------------------------------------------7--

      end module m_aslstep
