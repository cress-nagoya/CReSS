!***********************************************************************
      module m_slvstep
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/07/10
!     Modification: 2008/01/11, 2008/05/02, 2008/08/25, 2009/01/05,
!                   2009/01/30, 2009/02/27, 2009/11/13, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the number of time integration steps in slvdrv.

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

      public :: slvstep, s_slvstep

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface slvstep

        module procedure s_slvstep

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic int
      intrinsic min
      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_slvstep(fpiniopt,fpadvopt,fpcphopt,fpstime,fpetime,  &
     &                     fpdtbig,fpdtsml,fpdtvcul,fpdtcmph,           &
     &                     nbstp0,nbstp1,nsstp,nvstp,nclstp)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpiniopt
                       ! Formal parameter of unique index of iniopt

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fpstime
                       ! Formal parameter of unique index of stime

      integer, intent(in) :: fpetime
                       ! Formal parameter of unique index of etime

      integer, intent(in) :: fpdtbig
                       ! Formal parameter of unique index of dtbig

      integer, intent(in) :: fpdtsml
                       ! Formal parameter of unique index of dtsml

      integer, intent(in) :: fpdtvcul
                       ! Formal parameter of unique index of dtvcul

      integer, intent(in) :: fpdtcmph
                       ! Formal parameter of unique index of dtcmph

! Output variables

      integer(kind=i8), intent(out) :: nbstp0
                       ! Start index of large time steps

      integer(kind=i8), intent(out) :: nbstp1
                       ! End index of large time steps

      integer, intent(out) :: nsstp
                       ! Number of small time steps

      integer, intent(out) :: nvstp
                       ! Number of steps
                       ! of vertical Cubic Lagrange advection

      integer, intent(out) :: nclstp(0:3)
                       ! Number of steps of cloud micro physics

! Internal shared variables

      integer iniopt   ! Option for model initialization
      integer advopt   ! Option for advection scheme
      integer cphopt   ! Option for cloud micro physics

      integer ilcm     ! Parameter to get lowest common multiple

      integer stat     ! Runtime status

      real stime       ! Forecast start time
      real etime       ! Forecast end time

      real dtbig       ! Large time steps interval
      real dtsml       ! Small time steps interval

      real dtvcul      ! Time steps interval
                       ! of vertical Cubic Lagrange advection

      real dtcmph(1:3) ! Time interval of cloud micro physics

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpiniopt,iniopt)
      call getiname(fpadvopt,advopt)
      call getiname(fpcphopt,cphopt)
      call getrname(fpstime,stime)
      call getrname(fpetime,etime)
      call getrname(fpdtbig,dtbig)
      call getrname(fpdtsml,dtsml)
      call getrname(fpdtvcul,dtvcul)
      call getrname(fpdtcmph,dtcmph(1))
      call getrname(fpdtcmph+1,dtcmph(2))
      call getrname(fpdtcmph+2,dtcmph(3))

! -----

! Calculate the number of time integration steps.

      nbstp0=1000_i8*int(stime+.1e0,i8)                                 &
     &  /(10_i8*int(1.e2*(dtbig+.001e0),i8))+1_i8

      nbstp1=1000_i8*int(etime+.1e0,i8)                                 &
     &  /(10_i8*int(1.e2*(dtbig+.001e0),i8))

      if(advopt.le.3) then

        nsstp=20*int(1.e2*(dtbig+.001e0))/int(1.e3*(dtsml+.0001e0))

        nvstp=0

      else

        nsstp=10*int(1.e2*(dtbig+.001e0))/int(1.e3*(dtsml+.0001e0))

        nvstp=int(1.e2*(dtbig+.001e0))/int(1.e2*(dtvcul+.001e0))

      end if

      if(abs(cphopt).lt.10) then

        nclstp(0)=0
        nclstp(1)=0
        nclstp(2)=0
        nclstp(3)=0

      else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

        if(advopt.le.3) then

          nclstp(1)=20*int(1.e2*(dtbig+.001e0))                         &
     &      /int(1.e3*(dtcmph(1)+.0001e0))

          nclstp(2)=20*int(1.e2*(dtbig+.001e0))                         &
     &      /int(1.e3*(dtcmph(2)+.0001e0))

        else

          nclstp(1)=10*int(1.e2*(dtbig+.001e0))                         &
     &      /int(1.e3*(dtcmph(1)+.0001e0))

          nclstp(2)=10*int(1.e2*(dtbig+.001e0))                         &
     &      /int(1.e3*(dtcmph(2)+.0001e0))

        end if

        nclstp(0)=min(nclstp(1),nclstp(2))

        ilcm=0

        iterate: do

          ilcm=ilcm+1

          nclstp(0)=nclstp(0)*ilcm

          if(mod(nclstp(0),nclstp(1)).eq.0.and.                         &
     &       mod(nclstp(0),nclstp(2)).eq.0) then

            exit iterate

          else

            nclstp(0)=nclstp(0)/ilcm

          end if

        end do iterate

        nclstp(3)=0

      end if

! -----

! If error occured, call the procedure destroy.

      stat=0

      if(nbstp0.gt.nbstp1) then

        stat=stat+1

      end if

      if((mod(iniopt,10).ne.2.and.nbstp0.ne.1_i8)                       &
     &  .or.(mod(iniopt,10).eq.2.and.nbstp0.eq.1_i8)) then

        stat=stat+1

      end if

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('slvstep ',7,'cont',7,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('slvstep ',7,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

      end subroutine s_slvstep

!-----7--------------------------------------------------------------7--

      end module m_slvstep
