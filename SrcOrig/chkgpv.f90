!***********************************************************************
      module m_chkgpv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/04/15
!     Modification: 2004/05/31, 2004/06/10, 2004/08/20, 2006/01/10,
!                   2006/08/18, 2006/09/21, 2007/01/20, 2007/07/30,
!                   2008/05/02, 2008/08/25, 2008/12/11, 2009/02/27,
!                   2009/03/31, 2010/12/28, 2011/09/22, 2013/02/13,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     check the GPV data variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkmxn
      use m_comkind
      use m_getcname
      use m_getiname
      use m_inichar
      use m_outstd12

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: chkgpv, s_chkgpv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface chkgpv

        module procedure s_chkgpv

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_chkgpv(fpgpvvar,fpdatype_gpv,fpexbopt,ctime,stat,    &
     &                    nid,njd,nkd,htdat,zdat,udat,vdat,wdat,        &
     &                    pdat,ptdat,qvdat,qcdat,qrdat,                 &
     &                    qidat,qsdat,qgdat,qhdat)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpgpvvar
                       ! Formal parameter of unique index of gpvvar

      integer, intent(in) :: fpdatype_gpv
                       ! Formal parameter of unique index of datype_gpv

      integer, intent(in) :: fpexbopt
                       ! Formal parameter of unique index of exbopt

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

      integer, intent(in) :: nid
                       ! Data dimension in x direction

      integer, intent(in) :: njd
                       ! Data dimension in y direction

      integer, intent(in) :: nkd
                       ! Data dimension in z direction

      real, intent(in) :: htdat(1:nid,1:njd)
                       ! Terrain height in data

      real, intent(in) :: zdat(1:nid,1:njd,1:nkd)
                       ! z physical coordinates in data

      real, intent(in) :: udat(1:nid,1:njd,1:nkd)
                       ! x components of velocity in data

      real, intent(in) :: vdat(1:nid,1:njd,1:nkd)
                       ! y components of velocity in data

      real, intent(in) :: wdat(1:nid,1:njd,1:nkd)
                       ! z components of velocity in data

      real, intent(in) :: pdat(1:nid,1:njd,1:nkd)
                       ! Pressure in data

      real, intent(in) :: ptdat(1:nid,1:njd,1:nkd)
                       ! Potential temperature in data

      real, intent(in) :: qvdat(1:nid,1:njd,1:nkd)
                       ! Water vapor mixing ratio in data

      real, intent(in) :: qcdat(1:nid,1:njd,1:nkd)
                       ! Cloud water mixing ratio in data

      real, intent(in) :: qrdat(1:nid,1:njd,1:nkd)
                       ! Rain water mixing ratio in data

      real, intent(in) :: qidat(1:nid,1:njd,1:nkd)
                       ! Cloud ice mixing ratio in data

      real, intent(in) :: qsdat(1:nid,1:njd,1:nkd)
                       ! Snow mixing ratio in data

      real, intent(in) :: qgdat(1:nid,1:njd,1:nkd)
                       ! Graupel mixing ratio in data

      real, intent(in) :: qhdat(1:nid,1:njd,1:nkd)
                       ! Hail mixing ratio in data

! Output variable

      integer, intent(out) :: stat
                       ! Runtime status

! Internal shared variables

      character(len=108) gpvvar
                       ! Control flag of input GPV data variables

      character(len=108) datype_gpv
                       ! Control flag of GPV data type

      integer exbopt   ! Option for external boundary forcing

!-----7--------------------------------------------------------------7--

! Initialize the character variables.

      call inichar(gpvvar)
      call inichar(datype_gpv)

! -----

! Get the required namelist variables.

      call getcname(fpgpvvar,gpvvar)
      call getcname(fpdatype_gpv,datype_gpv)
      call getiname(fpexbopt,exbopt)

! -----

! Initialize the runtime status.

      stat=0

! -----

! Read in the message to standard i/o.

      call outstd12(-1,'    ',4,'       ',ctime,0,0,1,0.e0,0,0,1,0.e0)

! -----

! Check the GPV data variables.

      if(mod(exbopt,10).eq.2) then

        call s_chkmxn('ht  ',2,'[m]    ','all  ',ctime,stat,            &
     &                -1000.e0,10000.e0,0.e0,0.e0,nid,njd,1,htdat)

        if(stat.ne.0) then

          return

        end if

      end if

      call chkmxn('z   ',1,'[m]    ','all  ',ctime,stat,                &
     &            -1000.e0,1000000.e0,0.e0,0.e0,nid,njd,nkd,zdat)

      if(stat.ne.0) then

        return

      end if

      call chkmxn('u   ',1,'[m/s]  ','all  ',ctime,stat,                &
     &            -300.e0,300.e0,0.e0,0.e0,nid,njd,nkd,udat)

      if(stat.ne.0) then

        return

      end if

      call chkmxn('v   ',1,'[m/s]  ','all  ',ctime,stat,                &
     &            -300.e0,300.e0,0.e0,0.e0,nid,njd,nkd,vdat)

      if(stat.ne.0) then

        return

      end if

      if(gpvvar(1:1).eq.'o') then

        call chkmxn('w   ',1,'[m/s]  ','all  ',ctime,stat,              &
     &              -200.e0,200.e0,0.e0,0.e0,nid,njd,nkd,wdat)

        if(stat.ne.0) then

          return

        end if

      end if

      call chkmxn('p   ',1,'[Pa]   ','all  ',ctime,stat,                &
     &            0.e0,200000.e0,0.e0,0.e0,nid,njd,nkd,pdat)

      if(stat.ne.0) then

        return

      end if

      if(datype_gpv(1:1).eq.'p') then

        call chkmxn('pt  ',2,'[K]    ','all  ',ctime,stat,              &
     &              123.16e0,2273.16e0,0.e0,0.e0,nid,njd,nkd,ptdat)

        if(stat.ne.0) then

          return

        end if

      else if(datype_gpv(1:1).eq.'t') then

        call chkmxn('t   ',1,'[K]    ','all  ',ctime,stat,              &
     &              123.16e0,373.16e0,0.e0,0.e0,nid,njd,nkd,ptdat)

        if(stat.ne.0) then

          return

        end if

      end if

      if(gpvvar(2:2).eq.'o') then

        if(datype_gpv(2:2).eq.'m') then

          call chkmxn('qv  ',2,'[kg/kg]','all  ',ctime,stat,0.e0,.1e0,  &
     &                0.e0,0.e0,nid,njd,nkd,qvdat)

          if(stat.ne.0) then

            return

          end if

        else if(datype_gpv(2:2).eq.'r') then

          call chkmxn('qv  ',2,'[%]    ','all  ',ctime,stat,            &
     &                0.e0,100.e0,0.e0,0.e0,nid,njd,nkd,qvdat)

          if(stat.ne.0) then

            return

          end if

        end if

      end if

      if(gpvvar(3:3).eq.'o') then

        call chkmxn('qc  ',2,'[kg/kg]','all  ',ctime,stat,0.e0,.1e0,    &
     &              0.e0,0.e0,nid,njd,nkd,qcdat)

        if(stat.ne.0) then

          return

        end if

      end if

      if(gpvvar(4:4).eq.'o') then

        call chkmxn('qr  ',2,'[kg/kg]','all  ',ctime,stat,0.e0,.1e0,    &
     &              0.e0,0.e0,nid,njd,nkd,qrdat)

        if(stat.ne.0) then

          return

        end if

      end if

      if(gpvvar(5:5).eq.'o') then

        call chkmxn('qi  ',2,'[kg/kg]','all  ',ctime,stat,0.e0,.1e0,    &
     &              0.e0,0.e0,nid,njd,nkd,qidat)

        if(stat.ne.0) then

          return

        end if

      end if

      if(gpvvar(6:6).eq.'o') then

        call chkmxn('qs  ',2,'[kg/kg]','all  ',ctime,stat,0.e0,.1e0,    &
     &              0.e0,0.e0,nid,njd,nkd,qsdat)

        if(stat.ne.0) then

          return

        end if

      end if

      if(gpvvar(7:7).eq.'o') then

        call chkmxn('qg  ',2,'[kg/kg]','all  ',ctime,stat,0.e0,.1e0,    &
     &              0.e0,0.e0,nid,njd,nkd,qgdat)

        if(stat.ne.0) then

          return

        end if

      end if

      if(gpvvar(8:8).eq.'o') then

        call chkmxn('qh  ',2,'[kg/kg]','all  ',ctime,stat,0.e0,.1e0,    &
     &              0.e0,0.e0,nid,njd,nkd,qhdat)

        if(stat.ne.0) then

          return

        end if

      end if

! -----

      end subroutine s_chkgpv

!-----7--------------------------------------------------------------7--

      end module m_chkgpv
