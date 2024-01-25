!***********************************************************************
      module m_sststep
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/11/10
!     Modification: 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the number of steps of the main do loop in sstdrv.

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

      public :: sststep, s_sststep

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface sststep

        module procedure s_sststep

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
      subroutine s_sststep(fpsfcopt,fpsstitv,fpstime,fpetime,           &
     &                     nstp0,nstp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpsfcopt
                       ! Formal parameter of unique index of sfcopt

      integer, intent(in) :: fpsstitv
                       ! Formal parameter of unique index of sstitv

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

      integer sfcopt   ! Option for surface physics

      integer(kind=i8) sst103
                       ! 1000 x int(sstitv + 0.1)

      integer stat     ! Runtime status

      real stime       ! Forecast start time
      real etime       ! Forecast end time

      real sstitv      ! Time interval of
                       ! sea surface temperature data file

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpsfcopt,sfcopt)
      call getrname(fpsstitv,sstitv)
      call getrname(fpstime,stime)
      call getrname(fpetime,etime)

! -----

! Calculate the number of steps of the main do loop in sstdrv.

      if(sfcopt.eq.3.or.sfcopt.eq.13) then

        sst103=1000_i8*int(sstitv+.1e0,i8)

        nstp0=1000_i8*int(stime+.1e0,i8)/sst103+1_i8

        if(mod(1000_i8*int(etime+.1e0,i8),sst103).eq.0_i8) then

          nstp1=1000_i8*int(etime+.1e0,i8)/sst103+1_i8

        else

          nstp1=1000_i8*int(etime+.1e0,i8)/sst103+2_i8

        end if

        stat=0

        if(nstp0.gt.nstp1) then
          stat=1
        end if

        call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('sststep ',7,'cont',7,'              ',14,101, &
     &                   stat)

          end if

          call cpondpe

          call destroy('sststep ',7,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

      else

        nstp0=1_i8
        nstp1=1_i8

      end if

! -----

      end subroutine s_sststep

!-----7--------------------------------------------------------------7--

      end module m_sststep
