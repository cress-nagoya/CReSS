!***********************************************************************
      module m_chkrdr
!***********************************************************************

!     Author      : Sakakibara Atsushi, Monoe Daisuke
!     Date        : 2004/04/15
!     Modification: 2004/05/31, 2004/06/10, 2004/08/20, 2004/08/31,
!                   2006/01/10, 2006/09/21, 2007/01/20, 2007/07/30,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2009/03/31,
!                   2013/02/05, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     check the radar data variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkmxn
      use m_comkind
      use m_getcname
      use m_inichar
      use m_outstd12

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: chkrdr, s_chkrdr

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface chkrdr

        module procedure s_chkrdr

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

!     none

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_chkrdr(fprdrvar,fpdatype_rdr,ctime,stat,             &
     &                    nid,njd,nkd,zdat,udat,vdat,wdat,qpdat)
!***********************************************************************

! Input variables

      integer, intent(in) :: fprdrvar
                       ! Formal parameter of unique index of rdrvar

      integer, intent(in) :: fpdatype_rdr
                       ! Formal parameter of unique index of datype_rdr

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

      integer, intent(in) :: nid
                       ! Data dimension in x direction

      integer, intent(in) :: njd
                       ! Data dimension in y direction

      integer, intent(in) :: nkd
                       ! Data dimension in z direction

      real, intent(in) :: zdat(1:nid,1:njd,1:nkd)
                       ! z physical coordinates in data

      real, intent(in) :: udat(1:nid,1:njd,1:nkd)
                       ! x components of velocity in data

      real, intent(in) :: vdat(1:nid,1:njd,1:nkd)
                       ! y components of velocity in data

      real, intent(in) :: wdat(1:nid,1:njd,1:nkd)
                       ! z components of velocity in data

      real, intent(in) :: qpdat(1:nid,1:njd,1:nkd)
                       ! Precipitation mixing ratio in data

! Output variable

      integer, intent(out) :: stat
                       ! Runtime status

! Internal shared variables

      character(len=108) rdrvar
                       ! Control flag of input radar data variables

      character(len=108) datype_rdr
                       ! Control flag of radar data type

!-----7--------------------------------------------------------------7--

! Initialize the character variables.

      call inichar(rdrvar)
      call inichar(datype_rdr)

! -----

! Get the required namelist variables.

      call getcname(fprdrvar,rdrvar)
      call getcname(fpdatype_rdr,datype_rdr)

! -----

! Initialize the runtime status.

      stat=0

! -----

! Read in the message to standard i/o.

      call outstd12(-1,'    ',4,'       ',ctime,0,0,1,0.e0,0,0,1,0.e0)

! -----

! Check the radar data variables.

      call chkmxn('z   ',1,'[m]    ','all  ',ctime,stat,                &
     &            -1000.e0,1000000.e0,0.e0,0.e0,nid,njd,nkd,zdat)

      if(stat.ne.0) then

        return

      end if

      if(rdrvar(1:1).eq.'o') then

        call chkmxn('u   ',1,'[m/s]  ','undef',ctime,stat,              &
     &              -200.e0,200.e0,-1.e34,1.e34,nid,njd,nkd,udat)

        if(stat.ne.0) then

          return

        end if

      end if

      if(rdrvar(2:2).eq.'o') then

        call chkmxn('v   ',1,'[m/s]  ','undef',ctime,stat,              &
     &              -200.e0,200.e0,-1.e34,1.e34,nid,njd,nkd,vdat)

        if(stat.ne.0) then

          return

        end if

      end if

      if(rdrvar(3:3).eq.'o') then

        call chkmxn('w   ',1,'[m/s]  ','undef',ctime,stat,              &
     &              -200.e0,200.e0,-1.e34,1.e34,nid,njd,nkd,wdat)

        if(stat.ne.0) then

          return

        end if

      end if

      if(rdrvar(4:4).eq.'o') then

        if(datype_rdr(1:1).eq.'m') then

          call chkmxn('qp  ',2,'[kg/kg]','undef',ctime,stat,            &
     &                0.e0,.1e0,-1.e34,1.e34,nid,njd,nkd,qpdat)

          if(stat.ne.0) then

            return

          end if

        else if(datype_rdr(1:1).eq.'r') then

          call chkmxn('qp  ',2,'[dBZe] ','undef',ctime,stat,            &
     &                -1.e0,100.e0,-1.e34,1.e34,nid,njd,nkd,qpdat)

          if(stat.ne.0) then

            return

          end if

        end if

      end if

! -----

      end subroutine s_chkrdr

!-----7--------------------------------------------------------------7--

      end module m_chkrdr
