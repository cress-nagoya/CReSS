!***********************************************************************
      module m_ndgstep
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/09/09
!     Modification: 2002/10/31, 2002/12/02, 2003/03/28, 2003/05/19,
!                   2004/09/01, 2007/03/10, 2007/04/24, 2007/05/21,
!                   2007/07/30, 2008/05/02, 2008/07/25, 2008/08/25,
!                   2009/01/05, 2009/02/27, 2009/11/13, 2010/05/17,
!                   2011/01/19, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the control flag of analysis nudging.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
      use m_commath
      use m_getcname
      use m_getiname
      use m_getrname
      use m_inichar

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: ndgstep, s_ndgstep

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface ndgstep

        module procedure s_ndgstep

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic cos
      intrinsic int
      intrinsic min
      intrinsic mod
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_ndgstep(fpngrvar,fpnggopt,fpngropt,                  &
     &                     fpnggcoe,fpnggdlt,fpnggstr,fpnggend,fpnggc20,&
     &                     fpngrcoe,fpngrdlt,fpngrstr,fpngrend,fpngrc20,&
     &                     fpngraff,ctime,frdr,rtinc,nggdmp,ngrdmp)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpngrvar
                       ! Formal parameter of unique index of ngrvar

      integer, intent(in) :: fpnggopt
                       ! Formal parameter of unique index of nggopt

      integer, intent(in) :: fpngropt
                       ! Formal parameter of unique index of ngropt

      integer, intent(in) :: fpnggcoe
                       ! Formal parameter of unique index of nggcoe

      integer, intent(in) :: fpnggdlt
                       ! Formal parameter of unique index of nggdlt

      integer, intent(in) :: fpnggstr
                       ! Formal parameter of unique index of nggstr

      integer, intent(in) :: fpnggend
                       ! Formal parameter of unique index of nggend

      integer, intent(in) :: fpnggc20
                       ! Formal parameter of unique index of nggc20

      integer, intent(in) :: fpngrcoe
                       ! Formal parameter of unique index of ngrcoe

      integer, intent(in) :: fpngrdlt
                       ! Formal parameter of unique index of ngrdlt

      integer, intent(in) :: fpngrstr
                       ! Formal parameter of unique index of ngrstr

      integer, intent(in) :: fpngrend
                       ! Formal parameter of unique index of ngrend

      integer, intent(in) :: fpngrc20
                       ! Formal parameter of unique index of ngrc20

      integer, intent(in) :: fpngraff
                       ! Formal parameter of unique index of ngraff

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

      integer, intent(in) :: frdr(1:2)
                       ! Descriptor to put into motion
                       ! for radar data nudging

      real, intent(in) :: rtinc(1:2)
                       ! Lapse of forecast time from radar data reading

! Output variables

      real, intent(out) :: nggdmp
                       ! Analysis nudging damping coefficient for GPV

      real, intent(out) :: ngrdmp(1:2)
                       ! Analysis nudging damping coefficient for radar

! Internal shared variables

      character(len=108) ngrvar
                       ! Control flag of
                       ! analysis nudged variables to radar data

      integer nggopt   ! Option for analysis nudging to GPV
      integer ngropt   ! Option for analysis nudging to radar data

      real nggcoe      ! Analysis nudging coefficient to GPV
      real nggdlt      ! Time interval of analysis nudging to GPV
      real nggstr      ! Analysis nudging start time to GPV
      real nggend      ! Analysis nudging end time to GPV
      real nggc20      ! Start time to decrease nggcoe to 0 to GPV

      real ngrcoe      ! Analysis nudging coefficient to radar
      real ngrdlt      ! Time interval of analysis nudging to radar
      real ngrstr      ! Analysis nudging start time to radar
      real ngrend      ! Analysis nudging end time to radar
      real ngrc20      ! Start time to decrease ngrcoe to 0 to radar
      real ngraff      ! Analysis nudging affected time to radar

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(ngrvar)

! -----

! Get the required namelist variables.

      call getcname(fpngrvar,ngrvar)
      call getiname(fpnggopt,nggopt)
      call getiname(fpngropt,ngropt)
      call getrname(fpnggcoe,nggcoe)
      call getrname(fpnggdlt,nggdlt)
      call getrname(fpnggstr,nggstr)
      call getrname(fpnggend,nggend)
      call getrname(fpnggc20,nggc20)
      call getrname(fpngrcoe,ngrcoe)
      call getrname(fpngrdlt,ngrdlt)
      call getrname(fpngrstr,ngrstr)
      call getrname(fpngrend,ngrend)
      call getrname(fpngrc20,ngrc20)
      call getrname(fpngraff,ngraff)

! -----

! Get the control flag of analysis nudging to GPV.

      if(nggopt.eq.1) then

        if(ctime.ge.1000_i8*int(nggstr+.1e0,i8).and.                    &
     &     ctime.lt.1000_i8*int(nggend+.1e0,i8).and.                    &
     &     mod(ctime,10_i8*int(1.e2*(nggdlt+.001e0),i8)).eq.0_i8) then

          if(ctime.le.1000_i8*int(nggc20+.1e0,i8)) then

            nggdmp=nggcoe

          else

            nggdmp=min(1.e0,                                            &
     &        (.001e0*real(ctime)-nggc20)/(nggend-nggc20+eps))

            nggdmp=.5e0*nggcoe*(cos(cc*nggdmp)+1.e0)

          end if

        else

          nggdmp=-.1e0

        end if

      else

        nggdmp=-.1e0

      end if

! -----

! Get the control flag of analysis nudging to radar data.

      if(ngropt.ge.1) then

        if(ngropt.eq.1) then

          if(ctime.ge.1000_i8*int(ngrstr+.1e0,i8).and.                  &
     &       ctime.lt.1000_i8*int(ngrend+.1e0,i8).and.                  &
     &       mod(ctime,10_i8*int(1.e2*(ngrdlt+.001e0),i8)).eq.0_i8) then

            if(ctime.le.1000_i8*int(ngrc20+.1e0,i8)) then

              ngrdmp(1)=ngrcoe

            else

              ngrdmp(1)=min(1.e0,                                       &
     &          (.001e0*real(ctime)-ngrc20)/(ngrend-ngrc20+eps))

              ngrdmp(1)=.5e0*ngrcoe*(cos(cc*ngrdmp(1))+1.e0)

            end if

          else

            ngrdmp(1)=-.1e0

          end if

        else

          ngrdmp(1)=-.1e0

        end if

        if(ngrvar(1:3).ne.'xxx') then

          if(frdr(2).ge.0.and.rtinc(2).le.ngraff.and.                   &
     &       ctime.ge.1000_i8*int(ngrstr-ngraff+.1e0,i8).and.           &
     &       ctime.le.1000_i8*int(ngrend+ngraff+.1e0,i8).and.           &
     &       mod(ctime,10_i8*int(1.e2*(ngrdlt+.001e0),i8)).eq.0_i8) then

            ngrdmp(2)=min(1.e0,rtinc(2)/ngraff)

            ngrdmp(2)=.5e0*ngrcoe*(cos(cc*ngrdmp(2))+1.e0)

          else

            ngrdmp(2)=-.1e0

          end if

        else

          ngrdmp(2)=-.1e0

        end if

      else

        ngrdmp(1)=-.1e0
        ngrdmp(2)=-.1e0

      end if

! -----

      end subroutine s_ndgstep

!-----7--------------------------------------------------------------7--

      end module m_ndgstep
