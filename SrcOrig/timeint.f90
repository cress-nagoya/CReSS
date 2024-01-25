!***********************************************************************
      module m_timeint
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2008/07/01
!     Modification: 2008/08/25, 2009/01/30, 2009/02/27, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the time interval.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: timeint, s_timeint

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface timeint

        module procedure s_timeint

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic int
      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_timeint(fpsfcopt,fpadvopt,fpcphopt,fpdtbig,fpdtsml,  &
     &                     fpdtvcul,fpdtgrd,fpdtcmph,ctime,             &
     &                     dtb,dts,dtsep,dtsoil,dtcl)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpsfcopt
                       ! Formal parameter of unique index of sfcopt

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fpdtbig
                       ! Formal parameter of unique index of dtbig

      integer, intent(in) :: fpdtsml
                       ! Formal parameter of unique index of dtsml

      integer, intent(in) :: fpdtvcul
                       ! Formal parameter of unique index of dtvcul

      integer, intent(in) :: fpdtgrd
                       ! Formal parameter of unique index of dtgrd

      integer, intent(in) :: fpdtcmph
                       ! Formal parameter of unique index of dtcmph

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

! Output variables

      real, intent(out) :: dtb
                       ! Large time steps interval at
                       ! current processed time

      real, intent(out) :: dts
                       ! Small time steps interval at
                       ! current processed time

      real, intent(out) :: dtsep
                       ! Separated time steps interval
                       ! of vertical Cubic Lagrange advection

      real, intent(out) :: dtsoil
                       ! Time interval of soil temperature calculation
                       ! at current processed time

      real, intent(out) :: dtcl(1:3)
                       ! Time interval of cloud micro physics
                       ! at current processed time

! Internal shared variables

      integer sfcopt   ! Option for surface physics
      integer advopt   ! Option for advection scheme
      integer cphopt   ! Option for cloud micro physics

      real dtbig       ! Large time steps interval
      real dtsml       ! Small time steps interval

      real dtvcul      ! Time steps interval
                       ! of vertical Cubic Lagrange advection

      real dtgrd       ! Time interval of soil temperature calculation

      real dtcmph(1:3) ! Time interval of cloud micro physics

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpsfcopt,sfcopt)
      call getiname(fpadvopt,advopt)
      call getiname(fpcphopt,cphopt)
      call getrname(fpdtbig,dtbig)
      call getrname(fpdtsml,dtsml)
      call getrname(fpdtvcul,dtvcul)
      call getrname(fpdtgrd,dtgrd)
      call getrname(fpdtcmph,dtcmph(1))
      call getrname(fpdtcmph+1,dtcmph(2))
      call getrname(fpdtcmph+2,dtcmph(3))

! -----

! Calculate the time interval.

      if(advopt.le.3.and.ctime.eq.0_i8) then

        dtb=.5e0*dtbig
        dts=.5e0*dtsml

        dtsep=-.1e0

        if(sfcopt.eq.0) then
          dtsoil=-.1e0
        else
          dtsoil=.5e0*dtgrd
        end if

        if(abs(cphopt).lt.10) then

          dtcl(1)=dtb
          dtcl(2)=dtb
          dtcl(3)=dtb

        else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

          dtcl(1)=.5e0*dtcmph(1)
          dtcl(2)=.5e0*dtcmph(2)
          dtcl(3)=dtb

        end if

      else

        dtb=dtbig
        dts=dtsml

        dtsep=dtvcul

        if(sfcopt.eq.0) then

          dtsoil=-.1e0

        else

          if(mod(ctime,10_i8*int(1.e2*(dtgrd+.001e0),i8)).eq.0_i8) then
            dtsoil=dtgrd
          else
            dtsoil=-.1e0
          end if

        end if

        if(abs(cphopt).lt.10) then

          dtcl(1)=dtb
          dtcl(2)=dtb
          dtcl(3)=dtb

        else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

          dtcl(1)=dtcmph(1)
          dtcl(2)=dtcmph(2)
          dtcl(3)=dtb

        end if

      end if

! -----

      end subroutine s_timeint

!-----7--------------------------------------------------------------7--

      end module m_timeint
