!***********************************************************************
      module m_rdrstep
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/09/09
!     Modification: 2002/12/02, 2003/05/19, 2004/09/01, 2005/02/10,
!                   2007/01/20, 2007/05/21, 2007/07/30, 2008/05/02,
!                   2008/08/25, 2009/01/05, 2009/02/27, 2009/11/13,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the number of steps of the main do loop in rdrdrv.

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

      public :: rdrstep, s_rdrstep

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rdrstep

        module procedure s_rdrstep

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic int
      intrinsic max
      intrinsic min
      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_rdrstep(fpngropt,fprdritv,fpngrstr,fpngrend,         &
     &                     fpstime,fpetime,nstp0,nstp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpngropt
                       ! Formal parameter of unique index of ngropt

      integer, intent(in) :: fprdritv
                       ! Formal parameter of unique index of rdritv

      integer, intent(in) :: fpngrstr
                       ! Formal parameter of unique index of ngrstr

      integer, intent(in) :: fpngrend
                       ! Formal parameter of unique index of ngrend

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

      integer ngropt   ! Option for analysis nudging to radar data

      integer(kind=i8) rdr103
                       ! 1000 x int(rdritv + 0.1)

      integer stat     ! Runtime status

      real stime       ! Forecast start time
      real etime       ! Forecast end time

      real rdritv      ! Time interval of radar data file

      real ngrstr      ! Analysis nudging start time to radar data
      real ngrend      ! Analysis nudging end time to radar data

      real stmstr      ! stime - ngrstr
      real etmstr      ! etime - ngrstr or ngrend - ngrstr

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpngropt,ngropt)
      call getrname(fprdritv,rdritv)
      call getrname(fpngrstr,ngrstr)
      call getrname(fpngrend,ngrend)
      call getrname(fpstime,stime)
      call getrname(fpetime,etime)

! -----

! Calculate the number of steps of the main do loop in rdrdrv.

      if(ngropt.ge.1) then

        rdr103=1000_i8*int(rdritv+.1e0,i8)

        stmstr=max(stime-ngrstr,0.e0)
        etmstr=max(min(etime-ngrstr,ngrend-ngrstr),0.e0)

        if(ngropt.eq.1) then

          nstp0=1000_i8*int(stmstr+.1e0,i8)/rdr103+1_i8

          if(mod(1000_i8*int(etmstr+.1e0,i8),rdr103).eq.0_i8) then

            nstp1=1000_i8*int(etmstr+.1e0,i8)/rdr103+1_i8

          else

            nstp1=1000_i8*int(etmstr+.1e0,i8)/rdr103+2_i8

          end if

        else if(ngropt.ge.2) then

          nstp0=1000_i8*int(stmstr+.1e0,i8)/rdr103+1_i8
          nstp1=1000_i8*int(etmstr+.1e0,i8)/rdr103+1_i8

        end if

        stat=0

        if(nstp0.gt.nstp1) then
          stat=1
        end if

        call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('rdrstep ',7,'cont',7,'              ',14,101, &
     &                   stat)

          end if

          call cpondpe

          call destroy('rdrstep ',7,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

      else

        nstp0=1_i8
        nstp1=1_i8

      end if

! -----

      end subroutine s_rdrstep

!-----7--------------------------------------------------------------7--

      end module m_rdrstep
