!***********************************************************************
      module m_chknlslv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/01/20
!     Modification: 2007/03/10, 2007/04/24, 2007/05/21, 2007/06/27,
!                   2007/07/30, 2007/11/26, 2008/01/11, 2008/03/12,
!                   2008/04/17, 2008/05/02, 2008/07/01, 2008/08/25,
!                   2008/10/10, 2008/12/11, 2009/01/05, 2009/01/30,
!                   2009/02/27, 2009/03/23, 2009/11/13, 2009/12/18,
!                   2010/05/17, 2011/01/19, 2011/05/16, 2011/07/15,
!                   2011/08/09, 2011/08/18, 2011/09/22, 2011/11/10,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     check the namelist variables for solver.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
      use m_commath
      use m_commpi
      use m_comphy
      use m_defname
      use m_destroy
      use m_numchar
      use m_outstd14
      use m_outstd15

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: chknlslv, s_chknlslv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface chknlslv

        module procedure s_chknlslv

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic aint
      intrinsic ichar
      intrinsic int
      intrinsic mod
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_chknlslv(pname,ncpn,stat)
!***********************************************************************

! Input variables

      character(len=8), intent(in) :: pname
                       ! Running program name

      integer, intent(in) :: ncpn
                       ! Number of character of pname

! Output variable

      integer, intent(out) :: stat
                       ! Runtime status

! Internal shared variables

      integer ierr     ! Error descriptor

      integer fmsg     ! Control flag of message type for standard i/o

      integer ncspc    ! Number of space character of exprim and prvres

      integer ichar0   ! Character code of '0'
      integer ichar9   ! Character code of '9'

      real r4scl       ! Optional real scalar variable

!-----7--------------------------------------------------------------7--

! Set the word length for direct access file.

      r4scl=0.e0

      inquire(iolength=wlngth) r4scl

! -----

!! Check the namelist variables.

! Initialize the controler.

      stat=0

      ierr=0

      fmsg=0

! -----

! For the section, sysdep.

      if(savmem.ne.0.and.savmem.ne.1) then

        call destroy('chknlslv',8,'cont',201,'savmem        ',6,101,    &
     &               stat)

      end if

! -----

! For the section, runame.

      call numchar(exprim,1,ncexp,ncspc)

      if(ncexp.gt.64.or.ncspc.ne.0) then

        call destroy('chknlslv',8,'cont',201,'exprim        ',6,101,    &
     &               stat)

      end if

      if(pname(1:ncpn).ne.'check') then

        write(exprim(ncexp+1:ncexp+1),'(a1)') '.'

        ncexp=ncexp+1

      end if

! -----

! For the section, drname.

      call numchar(crsdir,1,nccrs,ncspc)

      if(crsdir(nccrs:nccrs).ne.'/') then

        nccrs=nccrs+1

        if(nccrs.le.108) then

          write(crsdir(nccrs:nccrs),'(a1)') '/'

        end if

      end if

      if(nccrs.gt.108.or.ncspc.ne.0) then

        call destroy('chknlslv',8,'cont',201,'crsdir        ',6,101,    &
     &               stat)

      end if

      call numchar(datdir,1,ncdat,ncspc)

      if(datdir(ncdat:ncdat).ne.'/') then

        ncdat=ncdat+1

        if(ncdat.le.108) then

          write(datdir(ncdat:ncdat),'(a1)') '/'

        end if

      end if

      if(ncdat.gt.108.or.ncspc.ne.0) then

        call destroy('chknlslv',8,'cont',201,'datdir        ',6,101,    &
     &               stat)

      end if

      if(pname(1:ncpn).ne.'check') then

        call outstd14('chknlslv',8,crsdir,datdir,nccrs,ncdat)

      end if

! -----

! For the section, dimset.

      if(xdim.lt.4) then

        call destroy('chknlslv',8,'cont',201,'xdim          ',4,101,    &
     &               stat)

      end if

      if(mpopt.eq.5) then

        if(ydim.ne.4) then

          call destroy('chknlslv',8,'cont',201,'ydim          ',4,101,  &
     &                 stat)

        end if

      else

        if(ydim.lt.4) then

          call destroy('chknlslv',8,'cont',201,'ydim          ',4,101,  &
     &                 stat)

        end if

      end if

      if(zdim.lt.7) then

        call destroy('chknlslv',8,'cont',201,'zdim          ',4,101,    &
     &               stat)

      end if

! -----

! For the section, project.

      if(mpopt.ne.0.and.mpopt.ne.1.and.mpopt.ne.2.and.mpopt.ne.3.and.   &
     &   mpopt.ne.4.and.mpopt.ne.5.and.mpopt.ne.10.and.mpopt.ne.13) then

        call destroy('chknlslv',8,'cont',201,'mpopt         ',5,101,    &
     &               stat)

      end if

      if(mpopt.eq.1.or.mpopt.eq.2.or.mpopt.eq.3.or.mpopt.eq.13) then

        if(nspol.ne.1.and.nspol.ne.-1) then

          call destroy('chknlslv',8,'cont',201,'nspol         ',5,101,  &
     &                 stat)

        end if

        if(tlat1.lt.-90.e0.or.tlat1.gt.90.e0) then

          call destroy('chknlslv',8,'cont',201,'tlat1         ',5,101,  &
     &                 stat)

        end if

      end if

      if(mpopt.eq.2) then

        if(tlat2.lt.-90.e0.or.tlat2.gt.90.e0) then

          call destroy('chknlslv',8,'cont',201,'tlat2         ',5,101,  &
     &                 stat)

        end if

      end if

      if(mpopt.eq.1.or.mpopt.eq.2.or.mpopt.eq.4) then

        if(tlon.lt.-180.e0.or.tlon.gt.180.e0) then

          call destroy('chknlslv',8,'cont',201,'tlon          ',4,101,  &
     &                 stat)

        end if

      end if

      if(ulat.lt.-90.e0.or.ulat.gt.90.e0) then

        call destroy('chknlslv',8,'cont',201,'ulat          ',4,101,    &
     &               stat)

      end if

      if(ulon.lt.-180.e0.or.ulon.gt.180.e0) then

        call destroy('chknlslv',8,'cont',201,'ulon          ',4,101,    &
     &               stat)

      end if

      if(mpopt.eq.5) then

        if(disr.lt.eps) then

          call destroy('chknlslv',8,'cont',201,'disr          ',4,101,  &
     &                 stat)

        end if

      end if

! -----

! For the section, gridset.

      if(mpopt.lt.10) then

        if(dx.lt.eps) then

          call destroy('chknlslv',8,'cont',201,'dx            ',2,101,  &
     &                 stat)

        else

          if(pname(1:ncpn).ne.'check') then

            if(mpopt.eq.0) then
              dx=rearth*d2r*dx
            end if

            dxiv=1.e0/dx

          end if

        end if

      else

        if(pname(1:ncpn).ne.'check') then

          if(xdim.ge.4) then

            if(mpopt.eq.10) then

              fmsg=fmsg+1

              dx=360.e0/real(xdim)

              call outstd15('dx    ',2,fmsg,2,dx)

              dx=rearth*d2r*dx

            end if

            if(mpopt.eq.13) then

              fmsg=fmsg+1

              dx=2.e0*cc*rearth/real(xdim-3)

              call outstd15('dx    ',2,fmsg,1,dx)

            end if

            dxiv=1.e0/dx

          end if

        end if

      end if

      if(mpopt.eq.5) then

        if(2.e0*dx.lt.disr) then

          call destroy('chknlslv',8,'cont',201,'dx            ',2,101,  &
     &                 stat)

        end if

      end if

      if(dy.lt.eps) then

        call destroy('chknlslv',8,'cont',201,'dy            ',2,101,    &
     &               stat)

      else

        if(pname(1:ncpn).ne.'check') then

          if(mpopt.eq.0.or.mpopt.eq.10) then
            dy=rearth*d2r*dy
          end if

          dyiv=1.e0/dy

        end if

      end if

      if(dz.lt.eps) then

        call destroy('chknlslv',8,'cont',201,'dz            ',2,101,    &
     &               stat)

      else

        if(pname(1:ncpn).ne.'check') then
          dziv=1.e0/dz
        end if

      end if

! -----

! For the section, gridsth.

      if(zsfc.lt.-1000.e0) then

        call destroy('chknlslv',8,'cont',201,'zsfc          ',4,101,    &
     &               stat)

      end if

      if(zflat.lt.zsfc) then

        call destroy('chknlslv',8,'cont',201,'zflat         ',5,101,    &
     &               stat)

      end if

      if(sthopt.ne.0.and.sthopt.ne.1.and.sthopt.ne.2) then

        call destroy('chknlslv',8,'cont',201,'sthopt        ',6,101,    &
     &               stat)

      end if

      if(sthopt.eq.0) then

        if(pname(1:ncpn).ne.'check') then
          dzmin=dz
        end if

      else

        if(dzmin.lt.eps.or.dzmin.gt.dz) then

          call destroy('chknlslv',8,'cont',201,'dzmin         ',5,101,  &
     &                 stat)

        end if

        if(layer1.lt.zsfc) then

          call destroy('chknlslv',8,'cont',201,'layer1        ',6,101,  &
     &                 stat)

        end if

        if(layer2.lt.layer1) then

          call destroy('chknlslv',8,'cont',201,'layer2        ',6,101,  &
     &                 stat)

        end if

      end if

! -----

! For the section, terrain.

      if(trnopt.ne.0.and.trnopt.ne.1.and.trnopt.ne.2) then

        call destroy('chknlslv',8,'cont',201,'trnopt        ',6,101,    &
     &               stat)

      end if

      if(trnopt.eq.0.or.trnopt.eq.1) then

        if(mnthgh(1).lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'mnthgh(1)     ',9,101,  &
     &                 stat)

        end if

      end if

      if(trnopt.eq.1) then

        if(mntwx.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'mntwx         ',5,101,  &
     &                 stat)

        end if

        if(mntwy.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'mntwy         ',5,101,  &
     &                 stat)

        end if

      end if

! -----

! For the section, flength.

      if(pname(1:ncpn).ne.'check') then
        idate(1:4)=sfcast(1:4)
        idate(5:6)=sfcast(6:7)
        idate(7:8)=sfcast(9:10)
        idate(9:10)=sfcast(12:13)
        idate(11:12)=sfcast(15:16)
      end if

      if(stime.lt.0.e0) then

        call destroy('chknlslv',8,'cont',201,'stime         ',5,101,    &
     &               stat)

      end if

      if(etime.lt.0.e0) then

        call destroy('chknlslv',8,'cont',201,'etime         ',5,101,    &
     &               stat)

      end if

      if(etime.lt.stime) then

        call destroy('chknlslv',8,'cont',201,'etime         ',5,101,    &
     &               stat)

      end if

! -----

! For the section, boundry.

      if(wbc.ne.-1.and.wbc.ne.1.and.wbc.ne.2.and.wbc.ne.3.and.          &
     &   wbc.ne.4.and.wbc.ne.5.and.wbc.ne.6.and.wbc.ne.7.and.           &
     &   wbc.ne.14.and.wbc.ne.15.and.wbc.ne.16) then

        call destroy('chknlslv',8,'cont',201,'wbc           ',3,101,    &
     &               stat)

      end if

      if(ebc.ne.-1.and.ebc.ne.1.and.ebc.ne.2.and.ebc.ne.3.and.          &
     &   ebc.ne.4.and.ebc.ne.5.and.ebc.ne.6.and.ebc.ne.7.and.           &
     &   ebc.ne.14.and.ebc.ne.15.and.ebc.ne.16) then

        call destroy('chknlslv',8,'cont',201,'ebc           ',3,101,    &
     &               stat)

      end if

      if(sbc.ne.-1.and.sbc.ne.1.and.sbc.ne.2.and.sbc.ne.3.and.          &
     &   sbc.ne.4.and.sbc.ne.5.and.sbc.ne.6.and.sbc.ne.7.and.           &
     &   sbc.ne.14.and.sbc.ne.15.and.sbc.ne.16) then

        call destroy('chknlslv',8,'cont',201,'sbc           ',3,101,    &
     &               stat)

      end if

      if(nbc.ne.-1.and.nbc.ne.1.and.nbc.ne.2.and.nbc.ne.3.and.          &
     &   nbc.ne.4.and.nbc.ne.5.and.nbc.ne.6.and.nbc.ne.7.and.           &
     &   nbc.ne.14.and.nbc.ne.15.and.nbc.ne.16) then

        call destroy('chknlslv',8,'cont',201,'nbc           ',3,101,    &
     &               stat)

      end if

      if(bbc.ne.2.and.bbc.ne.3) then

        call destroy('chknlslv',8,'cont',201,'bbc           ',3,101,    &
     &               stat)

      end if

      if(tbc.ne.2.and.tbc.ne.3) then

        call destroy('chknlslv',8,'cont',201,'tbc           ',3,101,    &
     &               stat)

      end if

      if(mpopt.eq.1.or.mpopt.eq.2.or.mpopt.eq.3) then

        if(abs(wbc).eq.1.or.abs(ebc).eq.1) then

          call destroy('chknlslv',8,'cont',201,'wbc           ',3,101,  &
     &                 stat)

          call destroy('chknlslv',8,'cont',201,'ebc           ',3,101,  &
     &                 stat)

        end if

      end if

      if(mpopt.eq.1.or.mpopt.eq.2.or.mpopt.eq.3                         &
     &  .or.mpopt.eq.10.or.mpopt.eq.13) then

        if(abs(sbc).eq.1.or.abs(nbc).eq.1) then

          call destroy('chknlslv',8,'cont',201,'sbc           ',3,101,  &
     &                 stat)

          call destroy('chknlslv',8,'cont',201,'nbc           ',3,101,  &
     &                 stat)

        end if

      end if

      if(mpopt.eq.5) then

        if(abs(wbc).eq.1.or.abs(ebc).eq.1) then

          call destroy('chknlslv',8,'cont',201,'wbc           ',3,101,  &
     &                 stat)

          call destroy('chknlslv',8,'cont',201,'ebc           ',3,101,  &
     &                 stat)

        end if

        if(abs(sbc).ne.1.or.abs(nbc).ne.1) then

          call destroy('chknlslv',8,'cont',201,'sbc           ',3,101,  &
     &                 stat)

          call destroy('chknlslv',8,'cont',201,'nbc           ',3,101,  &
     &                 stat)

        end if

      end if

      if(abs(wbc).eq.1.or.abs(ebc).eq.1) then

        if(wbc.ne.ebc) then

          call destroy('chknlslv',8,'cont',201,'wbc           ',3,101,  &
     &                 stat)

          call destroy('chknlslv',8,'cont',201,'ebc           ',3,101,  &
     &                 stat)

        end if

      end if

      if(xdim.eq.4) then

        if(wbc.ne.ebc) then

          call destroy('chknlslv',8,'cont',201,'wbc           ',3,101,  &
     &                 stat)

          call destroy('chknlslv',8,'cont',201,'ebc           ',3,101,  &
     &                 stat)

        end if

      end if

      if(abs(sbc).eq.1.or.abs(nbc).eq.1) then

        if(sbc.ne.nbc) then

          call destroy('chknlslv',8,'cont',201,'sbc           ',3,101,  &
     &                 stat)

          call destroy('chknlslv',8,'cont',201,'nbc           ',3,101,  &
     &                 stat)

        end if

      end if

      if(ydim.eq.4) then

        if(sbc.ne.nbc) then

          call destroy('chknlslv',8,'cont',201,'sbc           ',3,101,  &
     &                 stat)

          call destroy('chknlslv',8,'cont',201,'nbc           ',3,101,  &
     &                 stat)

        end if

      end if

      if(numpe.ne.xsub*ysub*xgroup*ygroup) then

        if(wbc.eq.1.or.ebc.eq.1) then

          call destroy('chknlslv',8,'cont',201,'wbc           ',3,101,  &
     &                 stat)

          call destroy('chknlslv',8,'cont',201,'ebc           ',3,101,  &
     &                 stat)

        end if

        if(sbc.eq.1.or.nbc.eq.1) then

          call destroy('chknlslv',8,'cont',201,'sbc           ',3,101,  &
     &                 stat)

          call destroy('chknlslv',8,'cont',201,'nbc           ',3,101,  &
     &                 stat)

        end if

      end if

      if(numpe.eq.1) then

        if(pname(1:ncpn).ne.'check') then

          if(wbc.eq.1) then
            wbc=-1
          end if

          if(ebc.eq.1) then
            ebc=-1
          end if

          if(sbc.eq.1) then
            sbc=-1
          end if

          if(nbc.eq.1) then
            nbc=-1
          end if

        end if

      end if

      if(xgroup.eq.1) then

        if(pname(1:ncpn).ne.'check') then

          if(wbc.eq.1) then
            wbc=-1
          end if

          if(ebc.eq.1) then
            ebc=-1
          end if

        end if

      end if

      if(ygroup.eq.1) then

        if(pname(1:ncpn).ne.'check') then

          if(sbc.eq.1) then
            sbc=-1
          end if

          if(nbc.eq.1) then
            nbc=-1
          end if

        end if

      end if

! -----

! For the section, gpvpram.

      ierr=0

      if(gpvvar(1:1).ne.'o'.and.gpvvar(1:1).ne.'x') then
        ierr=ierr+1
      end if

      if(gpvvar(2:2).ne.'o'.and.gpvvar(2:2).ne.'x') then
        ierr=ierr+1
      end if

      if(gpvvar(3:3).ne.'o'.and.gpvvar(3:3).ne.'x') then
        ierr=ierr+1
      end if

      if(gpvvar(4:4).ne.'o'.and.gpvvar(4:4).ne.'x') then
        ierr=ierr+1
      end if

      if(gpvvar(5:5).ne.'o'.and.gpvvar(5:5).ne.'x') then
        ierr=ierr+1
      end if

      if(gpvvar(6:6).ne.'o'.and.gpvvar(6:6).ne.'x') then
        ierr=ierr+1
      end if

      if(gpvvar(7:7).ne.'o'.and.gpvvar(7:7).ne.'x') then
        ierr=ierr+1
      end if

      if(gpvvar(8:8).ne.'o'.and.gpvvar(8:8).ne.'x') then
        ierr=ierr+1
      end if

      if(ierr.ne.0) then

        call destroy('chknlslv',8,'cont',201,'gpvvar        ',6,101,    &
     &               stat)

      end if

      if(gsmopt.ne.0.and.gsmopt.ne.1) then

        call destroy('chknlslv',8,'cont',201,'gsmopt        ',6,101,    &
     &               stat)

      end if

      if(gsmopt.eq.1) then

        if(gsmcnt.lt.0) then

          call destroy('chknlslv',8,'cont',201,'gsmcnt        ',6,101,  &
     &                 stat)

        end if

        if(gsmcoe.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'gsmcoe        ',6,101,  &
     &                 stat)

        end if

      end if

      if(nggopt.ne.0.and.nggopt.ne.1) then

        call destroy('chknlslv',8,'cont',201,'nggopt        ',6,101,    &
     &               stat)

      end if

      if(exbopt.ne.0.and.exbopt.ne.1.and.exbopt.ne.2                    &
     &  .and.exbopt.ne.11.and.exbopt.ne.12) then

        call destroy('chknlslv',8,'cont',201,'exbopt        ',6,101,    &
     &               stat)

      end if

      if(mod(exbopt,10).eq.2) then

        if(exbwid.lt.2) then

          call destroy('chknlslv',8,'cont',201,'exbwid        ',6,101,  &
     &                 stat)

        end if

        if((xdim.ne.4.and.2*exbwid.gt.xdim)                             &
     &    .or.(ydim.ne.4.and.2*exbwid.gt.ydim)) then

          call destroy('chknlslv',8,'cont',201,'exbwid        ',6,101,  &
     &                 stat)

        end if

      end if

      if(lspopt.ne.0.and.lspopt.ne.1.and.lspopt.ne.2                    &
     &  .and.lspopt.ne.11.and.lspopt.ne.12) then

        call destroy('chknlslv',8,'cont',201,'lspopt        ',6,101,    &
     &               stat)

      end if

      if(vspopt.ne.0.and.vspopt.ne.1.and.vspopt.ne.2) then

        call destroy('chknlslv',8,'cont',201,'vspopt        ',6,101,    &
     &               stat)

      end if

      if(nggopt.eq.1.or.exbopt.ge.1                                     &
     &  .or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

        if(int(gpvitv+.1e0).lt.60) then

          call destroy('chknlslv',8,'cont',201,'gpvitv        ',6,101,  &
     &                 stat)

        end if

        if(aslopt.ge.1) then

          if(int(aslitv+.1e0).lt.60) then

            call destroy('chknlslv',8,'cont',201,'aslitv        ',6,101,&
     &                   stat)

          end if

        end if

      end if

      if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

        ierr=0

        write(lbcvar(8:10),'(a3)') 'xxx'

        if(lbcvar(1:1).ne.'o'.and.lbcvar(1:1).ne.'x') then
          ierr=ierr+1
        end if

        if(lbcvar(2:2).ne.'o'.and.lbcvar(2:2).ne.'x') then
          ierr=ierr+1
        end if

        if(lbcvar(3:3).ne.'o'.and.lbcvar(3:3).ne.'x') then
          ierr=ierr+1
        end if

        if(lbcvar(4:4).ne.'o'.and.lbcvar(4:4).ne.'x') then
          ierr=ierr+1
        end if

        if(lbcvar(5:5).ne.'o'.and.lbcvar(5:5).ne.'x') then
          ierr=ierr+1
        end if

        if(lbcvar(6:6).ne.'o'.and.lbcvar(6:6).ne.'x') then
          ierr=ierr+1
        end if

        if(lbcvar(7:7).ne.'o'.and.lbcvar(7:7).ne.'x') then
          ierr=ierr+1
        end if

        if(lbcvar(8:8).ne.'o'.and.lbcvar(8:8).ne.'x') then
          ierr=ierr+1
        end if

        if(lbcvar(9:9).ne.'o'.and.lbcvar(9:9).ne.'x') then
          ierr=ierr+1
        end if

        if(lbcvar(10:10).ne.'o'.and.lbcvar(10:10).ne.'x') then
          ierr=ierr+1
        end if

        if(ierr.ne.0) then

          call destroy('chknlslv',8,'cont',201,'lbcvar        ',6,101,  &
     &                 stat)

        end if

        if(lbnews.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'lbnews        ',6,101,  &
     &                 stat)

        end if

        if(lbnorm.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'lbnorm        ',6,101,  &
     &                 stat)

        end if

      end if

      if(wbc.ge.6.or.ebc.ge.6.or.sbc.ge.6.or.nbc.ge.6) then

        if(gwave.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'gwave         ',5,101,  &
     &                 stat)

        end if

      end if

      if(nggopt.eq.1) then

        if(wbc.eq.2) then

          call destroy('chknlslv',8,'cont',201,'wbc           ',3,101,  &
     &                 stat)

        end if

        if(ebc.eq.2) then

          call destroy('chknlslv',8,'cont',201,'ebc           ',3,101,  &
     &                 stat)

        end if

        if(sbc.eq.2) then

          call destroy('chknlslv',8,'cont',201,'sbc           ',3,101,  &
     &                 stat)

        end if

        if(nbc.eq.2) then

          call destroy('chknlslv',8,'cont',201,'nbc           ',3,101,  &
     &                 stat)

        end if

        ierr=0

        if(nggvar(1:1).ne.'o'.and.nggvar(1:1).ne.'x') then
          ierr=ierr+1
        end if

        if(nggvar(2:2).ne.'o'.and.nggvar(2:2).ne.'x') then
          ierr=ierr+1
        end if

        if(nggvar(3:3).ne.'o'.and.nggvar(3:3).ne.'x') then
          ierr=ierr+1
        end if

        if(nggvar(4:4).ne.'o'.and.nggvar(4:4).ne.'x') then
          ierr=ierr+1
        end if

        if(nggvar(5:5).ne.'o'.and.nggvar(5:5).ne.'x') then
          ierr=ierr+1
        end if

        if(nggvar(6:6).ne.'o'.and.nggvar(6:6).ne.'x') then
          ierr=ierr+1
        end if

        if(nggvar(7:7).ne.'o'.and.nggvar(7:7).ne.'x') then
          ierr=ierr+1
        end if

        if(nggvar(8:8).ne.'o'.and.nggvar(8:8).ne.'x') then
          ierr=ierr+1
        end if

        if(gpvvar(1:1).eq.'x'.and.nggvar(3:3).eq.'o') then
          ierr=ierr+1
        end if

        if(gpvvar(2:2).eq.'x'.and.nggvar(6:6).eq.'o') then
          ierr=ierr+1
        end if

        if(gpvvar(3:8).eq.'xxxxxx'.and.nggvar(7:7).eq.'o') then
          ierr=ierr+1
        end if

        if(ierr.ne.0) then

          call destroy('chknlslv',8,'cont',201,'nggvar        ',6,101,  &
     &                 stat)

        end if

        if(nggcoe.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'nggcoe        ',6,101,  &
     &                 stat)

        end if

        if(nggdlt.lt..01e0) then

          call destroy('chknlslv',8,'cont',201,'nggdlt        ',6,101,  &
     &                 stat)

        end if

        if(dtbig.ge..01e0) then

          if(mod(10*int(1.e2*(nggdlt+.001e0)),                          &
     &      10*int(1.e2*(dtbig+.001e0))).ne.0) then

            call destroy('chknlslv',8,'cont',201,'nggdlt        ',6,101,&
     &                   stat)

          end if

        end if

        if(nggstr.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'nggstr        ',6,101,  &
     &                 stat)

        end if

        if(nggend.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'nggend        ',6,101,  &
     &                 stat)

        end if

        if(nggend.lt.nggstr) then

          call destroy('chknlslv',8,'cont',201,'nggend        ',6,101,  &
     &                 stat)

        end if

        if(nggc20.lt.nggstr) then

          call destroy('chknlslv',8,'cont',201,'nggc20        ',6,101,  &
     &                 stat)

        end if

        if(nggc20.gt.nggend) then

          call destroy('chknlslv',8,'cont',201,'nggc20        ',6,101,  &
     &                 stat)

        end if

        if(pname(1:ncpn).ne.'check') then
          nggstr=aint(nggstr+.1e0)
          nggend=aint(nggend+.1e0)
          nggc20=aint(nggc20+.1e0)
        end if

      end if

      if(exbopt.ge.1) then

        if(wbc.eq.2.or.wbc.eq.3) then

          call destroy('chknlslv',8,'cont',201,'wbc           ',3,101,  &
     &                 stat)

        end if

        if(ebc.eq.2.or.ebc.eq.3) then

          call destroy('chknlslv',8,'cont',201,'ebc           ',3,101,  &
     &                 stat)

        end if

        if(sbc.le.3) then

          call destroy('chknlslv',8,'cont',201,'sbc           ',3,101,  &
     &                 stat)

        end if

        if(nbc.le.3) then

          call destroy('chknlslv',8,'cont',201,'nbc           ',3,101,  &
     &                 stat)

        end if

        ierr=0

        if(exbvar(1:1).ne.'o'.and.exbvar(1:1).ne.'+'                    &
     &    .and.exbvar(1:1).ne.'-'.and.exbvar(1:1).ne.'x') then

          ierr=ierr+1

        end if

        if(exbvar(2:2).ne.'o'.and.exbvar(2:2).ne.'+'                    &
     &    .and.exbvar(2:2).ne.'-'.and.exbvar(2:2).ne.'x') then

          ierr=ierr+1

        end if

        if(exbvar(3:3).ne.'o'                                           &
     &    .and.exbvar(3:3).ne.'-'.and.exbvar(3:3).ne.'x') then

          ierr=ierr+1

        end if

        if(exbvar(4:4).ne.'o'                                           &
     &    .and.exbvar(4:4).ne.'-'.and.exbvar(4:4).ne.'x') then

          ierr=ierr+1

        end if

        if(exbvar(5:5).ne.'o'                                           &
     &    .and.exbvar(5:5).ne.'-'.and.exbvar(5:5).ne.'x') then

          ierr=ierr+1

        end if

        if(exbvar(6:6).ne.'o'                                           &
     &    .and.exbvar(6:6).ne.'-'.and.exbvar(6:6).ne.'x') then

          ierr=ierr+1

        end if

        if(exbvar(7:7).ne.'o'                                           &
     &    .and.exbvar(7:7).ne.'-'.and.exbvar(7:7).ne.'x') then

          ierr=ierr+1

        end if

        if(exbvar(8:8).ne.'o'                                           &
     &    .and.exbvar(8:8).ne.'-'.and.exbvar(8:8).ne.'x') then

          ierr=ierr+1

        end if

        if(gpvvar(1:1).eq.'x'.and.exbvar(3:3).ne.'x') then
          ierr=ierr+1
        end if

        if(gpvvar(2:2).eq.'x'.and.exbvar(6:6).ne.'x') then
          ierr=ierr+1
        end if

        if(gpvvar(3:8).eq.'xxxxxx'.and.exbvar(7:7).ne.'x') then
          ierr=ierr+1
        end if

        if(ierr.ne.0) then

          call destroy('chknlslv',8,'cont',201,'exbvar        ',6,101,  &
     &                 stat)

        end if

        if(exnews.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'exnews        ',6,101,  &
     &                 stat)

        end if

        if(exnorm.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'exnorm        ',6,101,  &
     &                 stat)

        end if

      end if

      if(mod(lspopt,10).eq.2) then

        if(vspopt.eq.1) then

          call destroy('chknlslv',8,'cont',201,'vspopt        ',6,101,  &
     &                 stat)

        end if

      end if

      if(nggopt.eq.1.or.exbopt.ge.1) then

        if(mod(lspopt,10).eq.2) then

          call destroy('chknlslv',8,'cont',201,'lspopt        ',6,101,  &
     &                 stat)

        end if

      end if

      if(lspopt.ge.1) then

        ierr=0

        if(lspvar(1:1).ne.'o'.and.lspvar(1:1).ne.'+'                    &
     &    .and.lspvar(1:1).ne.'-'.and.lspvar(1:1).ne.'x') then

          ierr=ierr+1

        end if

        if(lspvar(2:2).ne.'o'.and.lspvar(2:2).ne.'+'                    &
     &    .and.lspvar(2:2).ne.'-'.and.lspvar(2:2).ne.'x') then

          ierr=ierr+1

        end if

        if(lspvar(3:3).ne.'o'.and.lspvar(3:3).ne.'x') then
          ierr=ierr+1
        end if

        if(lspvar(4:4).ne.'o'.and.lspvar(4:4).ne.'x') then
          ierr=ierr+1
        end if

        if(lspvar(5:5).ne.'o'.and.lspvar(5:5).ne.'x') then
          ierr=ierr+1
        end if

        if(lspvar(6:6).ne.'o'.and.lspvar(6:6).ne.'x') then
          ierr=ierr+1
        end if

        if(lspvar(7:7).ne.'o'.and.lspvar(7:7).ne.'x') then
          ierr=ierr+1
        end if

        if(lspvar(8:8).ne.'o'.and.lspvar(8:8).ne.'x') then
          ierr=ierr+1
        end if

        if(lspvar(9:9).ne.'o'.and.lspvar(9:9).ne.'x') then
          ierr=ierr+1
        end if

        if(lspvar(10:10).ne.'o'.and.lspvar(10:10).ne.'x') then
          ierr=ierr+1
        end if

        if(ierr.ne.0) then

          call destroy('chknlslv',8,'cont',201,'lspvar        ',6,101,  &
     &                 stat)

        end if

        if(lspsmt.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'lspsmt        ',6,101,  &
     &                 stat)

        end if

        if(lsnews.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'lsnews        ',6,101,  &
     &                 stat)

        end if

        if(lsnorm.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'lsnorm        ',6,101,  &
     &                 stat)

        end if

        if(wdnews.lt.0) then

          call destroy('chknlslv',8,'cont',201,'wdnews        ',6,101,  &
     &                 stat)

        end if

        if(wdnorm.lt.0) then

          call destroy('chknlslv',8,'cont',201,'wdnorm        ',6,101,  &
     &                 stat)

        end if

        if((xdim.ne.4.and.2*(wdnews+2).gt.xdim+1)                       &
     &    .or.(ydim.ne.4.and.2*(wdnews+2).gt.ydim+1)) then

          call destroy('chknlslv',8,'cont',201,'wdnews        ',6,101,  &
     &                 stat)

        end if

        if((xdim.ne.4.and.2*(wdnorm+2).gt.xdim+1)                       &
     &    .or.(ydim.ne.4.and.2*(wdnorm+2).gt.ydim+1)) then

          call destroy('chknlslv',8,'cont',201,'wdnorm        ',6,101,  &
     &                 stat)

        end if

      end if

      if(vspopt.ge.1) then

        if(tbc.ne.2) then

          call destroy('chknlslv',8,'cont',201,'tbc           ',3,101,  &
     &                 stat)

        end if

        ierr=0

        if(vspvar(1:1).ne.'o'.and.vspvar(1:1).ne.'x') then
          ierr=ierr+1
        end if

        if(vspvar(2:2).ne.'o'.and.vspvar(2:2).ne.'x') then
          ierr=ierr+1
        end if

        if(vspvar(3:3).ne.'o'.and.vspvar(3:3).ne.'x') then
          ierr=ierr+1
        end if

        if(vspvar(4:4).ne.'o'.and.vspvar(4:4).ne.'x') then
          ierr=ierr+1
        end if

        if(vspvar(5:5).ne.'o'.and.vspvar(5:5).ne.'x') then
          ierr=ierr+1
        end if

        if(vspvar(6:6).ne.'o'.and.vspvar(6:6).ne.'x') then
          ierr=ierr+1
        end if

        if(vspvar(7:7).ne.'o'.and.vspvar(7:7).ne.'x') then
          ierr=ierr+1
        end if

        if(vspvar(8:8).ne.'o'.and.vspvar(8:8).ne.'x') then
          ierr=ierr+1
        end if

        if(vspvar(9:9).ne.'o'.and.vspvar(9:9).ne.'x') then
          ierr=ierr+1
        end if

        if(vspvar(10:10).ne.'o'.and.vspvar(10:10).ne.'x') then
          ierr=ierr+1
        end if

        if(ierr.ne.0) then

          call destroy('chknlslv',8,'cont',201,'vspvar        ',6,101,  &
     &                 stat)

        end if

        if(vspopt.eq.1) then

          if(vspgpv.lt.0.e0) then

            call destroy('chknlslv',8,'cont',201,'vspgpv        ',6,101,&
     &                   stat)

          end if

          if(botgpv.lt.zsfc) then

            call destroy('chknlslv',8,'cont',201,'botgpv        ',6,101,&
     &                   stat)

          end if

        end if

        if(vspbar.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'vspbar        ',6,101,  &
     &                 stat)

        end if

        if(botbar.lt.zsfc) then

          call destroy('chknlslv',8,'cont',201,'botbar        ',6,101,  &
     &                 stat)

        end if

      end if

      if(pname(1:ncpn).ne.'check') then
        gpvitv=aint(gpvitv+.1e0)
        aslitv=aint(aslitv+.1e0)
      end if

! -----

! For the section, rdrpram.

      ierr=0

      if(rdrvar(1:1).ne.'o'.and.rdrvar(1:1).ne.'x') then
        ierr=ierr+1
      end if

      if(rdrvar(2:2).ne.'o'.and.rdrvar(2:2).ne.'x') then
        ierr=ierr+1
      end if

      if(rdrvar(3:3).ne.'o'.and.rdrvar(3:3).ne.'x') then
        ierr=ierr+1
      end if

      if(rdrvar(4:4).ne.'o'.and.rdrvar(4:4).ne.'x') then
        ierr=ierr+1
      end if

      if(ierr.ne.0) then

        call destroy('chknlslv',8,'cont',201,'rdrvar        ',6,101,    &
     &               stat)

      end if

      if(ngropt.ne.0.and.ngropt.ne.1                                    &
     &  .and.ngropt.ne.2.and.ngropt.ne.12) then

        call destroy('chknlslv',8,'cont',201,'ngropt        ',6,101,    &
     &               stat)

      end if

      if(ngropt.ge.1) then

        if(rdrvar(1:4).eq.'xxxx') then

          call destroy('chknlslv',8,'cont',201,'rdrvar        ',6,101,  &
     &                 stat)

        end if

        if(int(rdritv+.1e0).lt.1.or.rdritv.lt.4.e0*dtbig) then

          call destroy('chknlslv',8,'cont',201,'rdritv        ',6,101,  &
     &                 stat)

        end if

        ierr=0

        if(ngrvar(1:1).ne.'o'.and.ngrvar(1:1).ne.'x') then
          ierr=ierr+1
        end if

        if(ngrvar(2:2).ne.'o'.and.ngrvar(2:2).ne.'x') then
          ierr=ierr+1
        end if

        if(ngrvar(3:3).ne.'o'.and.ngrvar(3:3).ne.'x') then
          ierr=ierr+1
        end if

        if(ngrvar(4:4).ne.'o'.and.ngrvar(4:4).ne.'x') then
          ierr=ierr+1
        end if

        if(ngrvar(5:5).ne.'o'.and.ngrvar(5:5).ne.'x') then
          ierr=ierr+1
        end if

        if(rdrvar(1:1).eq.'x'.and.ngrvar(1:1).eq.'o') then
          ierr=ierr+1
        end if

        if(rdrvar(2:2).eq.'x'.and.ngrvar(2:2).eq.'o') then
          ierr=ierr+1
        end if

        if(rdrvar(3:3).eq.'x'.and.ngrvar(3:3).eq.'o') then
          ierr=ierr+1
        end if

        if(rdrvar(4:4).eq.'x'.and.ngrvar(4:4).eq.'o') then
          ierr=ierr+1
        end if

        if(rdrvar(4:4).eq.'x'.and.ngrvar(5:5).eq.'o') then
          ierr=ierr+1
        end if

        if(ierr.ne.0) then

          call destroy('chknlslv',8,'cont',201,'ngrvar        ',6,101,  &
     &                 stat)

        end if

        if(ngrcoe.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'ngrcoe        ',6,101,  &
     &                 stat)

        end if

        if(ngrdlt.lt..01e0) then

          call destroy('chknlslv',8,'cont',201,'ngrdlt        ',6,101,  &
     &                 stat)

        end if

        if(dtbig.ge..01e0) then

          if(mod(10*int(1.e2*(ngrdlt+.001e0)),                          &
     &      10*int(1.e2*(dtbig+.001e0))).ne.0) then

            call destroy('chknlslv',8,'cont',201,'ngrdlt        ',6,101,&
     &                   stat)

          end if

        end if

        if(ngrstr.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'ngrstr        ',6,101,  &
     &                 stat)

        end if

        if(ngrend.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'ngrend        ',6,101,  &
     &                 stat)

        end if

        if(ngrend.lt.ngrstr) then

          call destroy('chknlslv',8,'cont',201,'ngrend        ',6,101,  &
     &                 stat)

        end if

        if(pname(1:ncpn).ne.'check') then
          ngrstr=aint(ngrstr+.1e0)
          ngrend=aint(ngrend+.1e0)
        end if

        if(ngropt.eq.1) then

          if(ngrc20.lt.ngrstr) then

            call destroy('chknlslv',8,'cont',201,'ngrc20        ',6,101,&
     &                   stat)

          end if

          if(ngrc20.gt.ngrend) then

            call destroy('chknlslv',8,'cont',201,'ngrc20        ',6,101,&
     &                   stat)

          end if

          if(pname(1:ncpn).ne.'check') then
            ngrc20=aint(ngrc20+.1e0)
          end if

        end if

        if(ngropt.ge.2.or.ngrvar(1:3).ne.'xxx') then

          if(ngraff.lt.dtbig.or.2.e0*ngraff.gt.rdritv) then

            call destroy('chknlslv',8,'cont',201,'ngraff        ',6,101,&
     &                   stat)

          end if

        end if

      end if

      if(pname(1:ncpn).ne.'check') then
        rdritv=aint(rdritv+.1e0)
      end if

! -----

! For the section, sfcphys.

      ierr=0

      if(sfcdat(1:1).ne.'o'.and.sfcdat(1:1).ne.'x') then
        ierr=ierr+1
      end if

      if(sfcdat(2:2).ne.'o'.and.sfcdat(2:2).ne.'x') then
        ierr=ierr+1
      end if

      if(sfcdat(3:3).ne.'o'.and.sfcdat(3:3).ne.'x') then
        ierr=ierr+1
      end if

      if(ierr.ne.0) then

        call destroy('chknlslv',8,'cont',201,'sfcdat        ',6,101,    &
     &               stat)

      end if

      if(sfcopt.ne.0.and.sfcopt.ne.1.and.sfcopt.ne.2.and.sfcopt.ne.3    &
     &  .and.sfcopt.ne.11.and.sfcopt.ne.12.and.sfcopt.ne.13) then

        call destroy('chknlslv',8,'cont',201,'sfcopt        ',6,101,    &
     &               stat)

      end if

      if(sfcopt.ge.1) then

        if(dzmin.gt.zamax) then

          call destroy('chknlslv',8,'cont',201,'dz            ',2,101,  &
     &                 stat)

          call destroy('chknlslv',8,'cont',201,'dzmin         ',5,101,  &
     &                 stat)

        end if

        if(tubopt.eq.0) then

          if(levpbl.lt.1.or.levpbl.gt.zdim-2) then

            call destroy('chknlslv',8,'cont',201,'levpbl        ',6,101,&
     &                   stat)

          end if

        else

          if(pname(1:ncpn).ne.'check') then
            levpbl=1
          end if

        end if

        if(levund.lt.4.or.levund.gt.zdim) then

          call destroy('chknlslv',8,'cont',201,'levund        ',6,101,  &
     &                 stat)

        end if

        if(dtgrd.lt..01e0) then

          call destroy('chknlslv',8,'cont',201,'dtgrd         ',5,101,  &
     &                 stat)

        end if

        if(dtbig.ge..01e0) then

          if(mod(10*int(1.e2*(dtgrd+.001e0)),                           &
     &      10*int(1.e2*(dtbig+.001e0))).ne.0) then

            call destroy('chknlslv',8,'cont',201,'dtgrd         ',5,101,&
     &                   stat)

          end if

        end if

        if(dzgrd.lt.eps) then

          call destroy('chknlslv',8,'cont',201,'dzgrd         ',5,101,  &
     &                 stat)

        end if

        if(dzsea.lt.eps) then

          call destroy('chknlslv',8,'cont',201,'dzsea         ',5,101,  &
     &                 stat)

        end if

        if(tgdeep.lt.173.16e0.or.tgdeep.gt.373.16e0) then

          call destroy('chknlslv',8,'cont',201,'tgdeep        ',6,101,  &
     &                 stat)

        end if

        if(sfcopt.gt.10) then

          call numchar(prvres,1,ncprv,ncspc)

          if(ncspc.ne.0) then

            call destroy('chknlslv',8,'cont',201,'prvres        ',6,101,&
     &                   stat)

          end if

          ierr=0

          if(prvres(ncprv-11:ncprv-8).ne.'.res') then
            ierr=1
          end if

          if(ierr.ne.0) then

            call destroy('chknlslv',8,'cont',201,'prvres        ',6,101,&
     &                   stat)

          end if

          ichar0=ichar('0')
          ichar9=ichar('9')

          if(ichar0.lt.ichar9) then

            if(ichar(prvres(ncprv:ncprv)).lt.ichar0                     &
     &        .or.ichar(prvres(ncprv:ncprv)).gt.ichar9) then

              ierr=ierr+1

            end if

            if(ichar(prvres(ncprv-1:ncprv-1)).lt.ichar0                 &
     &        .or.ichar(prvres(ncprv-1:ncprv-1)).gt.ichar9) then

              ierr=ierr+1

            end if

            if(ichar(prvres(ncprv-2:ncprv-2)).lt.ichar0                 &
     &        .or.ichar(prvres(ncprv-2:ncprv-2)).gt.ichar9) then

              ierr=ierr+1

            end if

            if(ichar(prvres(ncprv-3:ncprv-3)).lt.ichar0                 &
     &        .or.ichar(prvres(ncprv-3:ncprv-3)).gt.ichar9) then

              ierr=ierr+1

            end if

            if(ichar(prvres(ncprv-4:ncprv-4)).lt.ichar0                 &
     &        .or.ichar(prvres(ncprv-4:ncprv-4)).gt.ichar9) then

              ierr=ierr+1

            end if

            if(ichar(prvres(ncprv-5:ncprv-5)).lt.ichar0                 &
     &        .or.ichar(prvres(ncprv-5:ncprv-5)).gt.ichar9) then

              ierr=ierr+1

            end if

            if(ichar(prvres(ncprv-6:ncprv-6)).lt.ichar0                 &
     &        .or.ichar(prvres(ncprv-6:ncprv-6)).gt.ichar9) then

              ierr=ierr+1

            end if

            if(ichar(prvres(ncprv-7:ncprv-7)).lt.ichar0                 &
     &        .or.ichar(prvres(ncprv-7:ncprv-7)).gt.ichar9) then

              ierr=ierr+1

            end if

          else

            if(ichar(prvres(ncprv:ncprv)).lt.ichar9                     &
     &        .or.ichar(prvres(ncprv:ncprv)).gt.ichar0) then

              ierr=ierr+1

            end if

            if(ichar(prvres(ncprv-1:ncprv-1)).lt.ichar9                 &
     &        .or.ichar(prvres(ncprv-1:ncprv-1)).gt.ichar0) then

              ierr=ierr+1

            end if

            if(ichar(prvres(ncprv-2:ncprv-2)).lt.ichar9                 &
     &        .or.ichar(prvres(ncprv-2:ncprv-2)).gt.ichar0) then

              ierr=ierr+1

            end if

            if(ichar(prvres(ncprv-3:ncprv-3)).lt.ichar9                 &
     &        .or.ichar(prvres(ncprv-3:ncprv-3)).gt.ichar0) then

              ierr=ierr+1

            end if

            if(ichar(prvres(ncprv-4:ncprv-4)).lt.ichar9                 &
     &        .or.ichar(prvres(ncprv-4:ncprv-4)).gt.ichar0) then

              ierr=ierr+1

            end if

            if(ichar(prvres(ncprv-5:ncprv-5)).lt.ichar9                 &
     &        .or.ichar(prvres(ncprv-5:ncprv-5)).gt.ichar0) then

              ierr=ierr+1

            end if

            if(ichar(prvres(ncprv-6:ncprv-6)).lt.ichar9                 &
     &        .or.ichar(prvres(ncprv-6:ncprv-6)).gt.ichar0) then

              ierr=ierr+1

            end if

            if(ichar(prvres(ncprv-7:ncprv-7)).lt.ichar9                 &
     &        .or.ichar(prvres(ncprv-7:ncprv-7)).gt.ichar0) then

              ierr=ierr+1

            end if

          end if

          if(ierr.ne.0) then

            call destroy('chknlslv',8,'cont',201,'prvres        ',6,101,&
     &                   stat)

          end if

        end if

        if(sfcdat(1:1).eq.'x') then

          if(lnduse.ne.1.and.lnduse.ne.2.and.lnduse.ne.3) then

            call destroy('chknlslv',8,'cont',201,'lnduse        ',6,101,&
     &                   stat)

          end if

          if(gralbe.lt.0.e0.or.gralbe.gt.1.e0) then

            call destroy('chknlslv',8,'cont',201,'gralbe        ',6,101,&
     &                   stat)

          end if

          if(grbeta.lt.0.e0.or.grbeta.gt.1.e0) then

            call destroy('chknlslv',8,'cont',201,'grbeta        ',6,101,&
     &                   stat)

          end if

          if(grz0m.lt.0.e0) then

            call destroy('chknlslv',8,'cont',201,'grz0m         ',5,101,&
     &                   stat)

          end if

          if(grz0h.lt.0.e0) then

            call destroy('chknlslv',8,'cont',201,'grz0h         ',5,101,&
     &                   stat)

          end if

          if(grcap.lt.0.e0) then

            call destroy('chknlslv',8,'cont',201,'grcap         ',5,101,&
     &                   stat)

          end if

          if(grnuu.lt.0.e0) then

            call destroy('chknlslv',8,'cont',201,'grnuu         ',5,101,&
     &                   stat)

          end if

        end if

        if(sfcdat(2:2).eq.'x') then

          if(sfcopt.eq.3.or.sfcopt.eq.13) then

            call destroy('chknlslv',8,'cont',201,'sfcopt        ',6,101,&
     &                   stat)

          end if

          if(sstcst.lt.268.16e0.or.sstcst.gt.323.16e0) then

            call destroy('chknlslv',8,'cont',201,'sstcst        ',6,101,&
     &                   stat)

          end if

        else if(sfcdat(2:2).eq.'o') then

          if(sfcopt.eq.3.or.sfcopt.eq.13) then

            if(int(sstitv+.1e0).lt.60) then

              call destroy('chknlslv',8,'cont',201,'sstitv        ',6,  &
     &                     101,stat)

            end if

            if(pname(1:ncpn).ne.'check') then
              sstitv=aint(sstitv+.1e0)
            end if

          end if

        end if

        if(sfcdat(3:3).eq.'o') then

          if(dstopt.ne.1.and.dstopt.ne.2) then

            call destroy('chknlslv',8,'cont',201,'dstopt        ',6,101,&
     &                   stat)

          end if

        end if

      end if

! -----

! For the section, initype.

      if(iniopt.ne.1.and.iniopt.ne.2                                    &
     &  .and.iniopt.ne.3.and.iniopt.ne.12) then

        call destroy('chknlslv',8,'cont',201,'iniopt        ',6,101,    &
     &               stat)

      end if

      if(iniopt.eq.1) then

        if(nggopt.eq.1) then

          call destroy('chknlslv',8,'cont',201,'nggopt        ',6,101,  &
     &                 stat)

        end if

        if(mod(exbopt,10).eq.1.or.mod(exbopt,10).eq.2) then

          call destroy('chknlslv',8,'cont',201,'exbopt        ',6,101,  &
     &                 stat)

        end if

        if(mod(lspopt,10).eq.1) then

          call destroy('chknlslv',8,'cont',201,'lspopt        ',6,101,  &
     &                 stat)

        end if

        if(vspopt.eq.1) then

          call destroy('chknlslv',8,'cont',201,'vspopt        ',6,101,  &
     &                 stat)

        end if

        if(snddim.lt.1) then

          call destroy('chknlslv',8,'cont',201,'snddim        ',6,101,  &
     &                 stat)

        end if

        if(sndtyp(1:1).ne.'z'.and.sndtyp(1:1).ne.'p') then

          call destroy('chknlslv',8,'cont',201,'sndtyp        ',6,101,  &
     &                 stat)

        end if

        if(sndtyp(2:2).ne.'t'.and.sndtyp(2:2).ne.'p') then

          call destroy('chknlslv',8,'cont',201,'sndtyp        ',6,101,  &
     &                 stat)

        end if

        if(sndtyp(3:3).ne.'m'.and.sndtyp(3:3).ne.'r') then

          call destroy('chknlslv',8,'cont',201,'sndtyp        ',6,101,  &
     &                 stat)

        end if

        if(zsnd0.lt.-1000.e0) then

          call destroy('chknlslv',8,'cont',201,'zsnd0         ',5,101,  &
     &                 stat)

        end if

        if(psnd0.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'psnd0         ',5,101,  &
     &                 stat)

        end if

      end if

      if(masopt.ne.0.and.masopt.ne.1) then

        call destroy('chknlslv',8,'cont',201,'masopt        ',6,101,    &
     &               stat)

      end if

      if(masopt.eq.1) then

        if(maseps.lt.eps) then

          call destroy('chknlslv',8,'cont',201,'maseps        ',6,101,  &
     &                 stat)

        end if

        if(alpha1.lt.eps) then

          call destroy('chknlslv',8,'cont',201,'alpha1        ',6,101,  &
     &                 stat)

        end if

        if(alpha2.lt.eps) then

          call destroy('chknlslv',8,'cont',201,'alpha2        ',6,101,  &
     &                 stat)

        end if

      end if

! -----

! For the section, gridmove.

      if(movopt.ne.0.and.movopt.ne.1) then

        call destroy('chknlslv',8,'cont',201,'movopt        ',6,101,    &
     &               stat)

      end if

      if(movopt.eq.1) then

        if(iniopt.eq.3) then

          call destroy('chknlslv',8,'cont',201,'iniopt        ',6,101,  &
     &                 stat)

        end if

        if(umove.lt.-100.e0.or.umove.gt.100.e0) then

          call destroy('chknlslv',8,'cont',201,'umove         ',5,101,  &
     &                 stat)

        end if

        if(vmove.lt.-100.e0.or.vmove.gt.100.e0) then

          call destroy('chknlslv',8,'cont',201,'vmove         ',5,101,  &
     &                 stat)

        end if

      end if

      if(sfcopt.ge.1) then

        if(movopt.eq.1) then

          call destroy('chknlslv',8,'cont',201,'movopt        ',6,101,  &
     &                 stat)

        end if

      end if

! -----

! For the section, ptinicon.

      if(pt0opt.ne.0.and.pt0opt.ne.1.and.pt0opt.ne.2                    &
     &  .and.pt0opt.ne.3.and.pt0opt.ne.4.and.pt0opt.ne.5) then

        call destroy('chknlslv',8,'cont',201,'pt0opt        ',6,101,    &
     &               stat)

      end if

      if(pt0opt.eq.1.or.pt0opt.eq.2) then

        if(pt0num.lt.0.or.pt0num.gt.256) then

          call destroy('chknlslv',8,'cont',201,'pt0num        ',6,101,  &
     &                 stat)

        end if

      end if

! -----

! For the section, integrat.

      if(dtbig.lt..01e0) then

        call destroy('chknlslv',8,'cont',201,'dtbig         ',5,101,    &
     &               stat)

      else

        if(mod(1000_i8*int(stime+.1e0,i8),                              &
     &    10_i8*int(1.e2*(dtbig+.001e0),i8)).ne.0_i8) then

          call destroy('chknlslv',8,'cont',201,'stime         ',5,101,  &
     &                 stat)

        end if

        if(mod(1000_i8*int(etime+.1e0,i8),                              &
     &    10_i8*int(1.e2*(dtbig+.001e0),i8)).ne.0_i8) then

          call destroy('chknlslv',8,'cont',201,'etime         ',5,101,  &
     &                 stat)

        end if

      end if

      if(dtsml.lt..001e0) then

        call destroy('chknlslv',8,'cont',201,'dtsml         ',5,101,    &
     &               stat)

      end if

      if(gwmopt.ne.0.and.gwmopt.ne.1) then

        call destroy('chknlslv',8,'cont',201,'gwmopt        ',6,101,    &
     &               stat)

      end if

      if(impopt.ne.0.and.impopt.ne.1                                    &
     &  .and.impopt.ne.2.and.impopt.ne.3) then

        call destroy('chknlslv',8,'cont',201,'impopt        ',6,101,    &
     &               stat)

      end if

      if(impopt.ge.1) then

        if(bbc.ne.2) then

          call destroy('chknlslv',8,'cont',201,'bbc           ',3,101,  &
     &                 stat)

        end if

        if(tbc.ne.2) then

          call destroy('chknlslv',8,'cont',201,'tbc           ',3,101,  &
     &                 stat)

        end if

        if(weicoe.lt.0.e0.or.weicoe.gt.1.e0) then

          call destroy('chknlslv',8,'cont',201,'weicoe        ',6,101,  &
     &                 stat)

        end if

      end if

      if(impopt.eq.3) then

        if(gsdeps.lt.eps) then

          call destroy('chknlslv',8,'cont',201,'gsdeps        ',6,101,  &
     &                 stat)

        end if

      end if

      if(advopt.ne.1.and.advopt.ne.2                                    &
     &  .and.advopt.ne.3.and.advopt.ne.4.and.advopt.ne.5) then

        call destroy('chknlslv',8,'cont',201,'advopt        ',6,101,    &
     &               stat)

      end if

      if(advopt.le.3) then

        if(filcoe.lt.0.e0.or.filcoe.gt.1.e0) then

          call destroy('chknlslv',8,'cont',201,'filcoe        ',6,101,  &
     &                 stat)

        end if

        if(dtsml.ge..001e0) then

          if(mod(20*int(1.e2*(dtbig+.001e0)),                           &
     &      int(1.e3*(dtsml+.0001e0))).ne.0) then

            call destroy('chknlslv',8,'cont',201,'dtbig         ',5,101,&
     &                   stat)

            call destroy('chknlslv',8,'cont',201,'dtsml         ',5,101,&
     &                   stat)

          end if

        end if

      else

        if(pname(1:ncpn).ne.'check') then

          if(wbc.eq.4.or.wbc.eq.5) then
            wbc=6
          else if(wbc.eq.14.or.wbc.eq.15) then
            wbc=16
          end if

          if(ebc.eq.4.or.ebc.eq.5) then
            ebc=6
          else if(ebc.eq.14.or.ebc.eq.15) then
            ebc=16
          end if

          if(sbc.eq.4.or.sbc.eq.5) then
            sbc=6
          else if(sbc.eq.14.or.sbc.eq.15) then
            sbc=16
          end if

          if(nbc.eq.4.or.nbc.eq.5) then
            nbc=6
          else if(nbc.eq.14.or.nbc.eq.15) then
            nbc=16
          end if

        end if

        if(dtsml.ge..001e0) then

          if(mod(10*int(1.e2*(dtbig+.001e0)),                           &
     &      int(1.e3*(dtsml+.0001e0))).ne.0) then

            call destroy('chknlslv',8,'cont',201,'dtbig         ',5,101,&
     &                   stat)

            call destroy('chknlslv',8,'cont',201,'dtsml         ',5,101,&
     &                   stat)

          end if

        end if

        if(dtvcul.lt..01e0) then

          call destroy('chknlslv',8,'cont',201,'dtvcul        ',6,101,  &
     &                 stat)

        else

          if(mod(int(1.e2*(dtbig+.001e0)),                              &
     &      int(1.e2*(dtvcul+.001e0))).ne.0) then

            call destroy('chknlslv',8,'cont',201,'dtbig         ',5,101,&
     &                   stat)

            call destroy('chknlslv',8,'cont',201,'dtvcul        ',6,101,&
     &                   stat)

          end if

        end if

      end if

! -----

! For the section, smoother.

      if(smtopt.ne.0.and.smtopt.ne.1.and.smtopt.ne.2.and.smtopt.ne.3    &
     &  .and.smtopt.ne.11.and.smtopt.ne.12.and.smtopt.ne.13) then

        call destroy('chknlslv',8,'cont',201,'smtopt        ',6,101,    &
     &               stat)

      end if

      if(smtopt.gt.0) then

        if(smhcoe.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'smhcoe        ',6,101,  &
     &                 stat)

        end if

        if(smvcoe.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'smvcoe        ',6,101,  &
     &                 stat)

        end if

      end if

      if(smtopt.gt.10) then

        if(nlhcoe.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'nlhcoe        ',6,101,  &
     &                 stat)

        end if

        if(nlvcoe.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'nlvcoe        ',6,101,  &
     &                 stat)

        end if

      end if

! -----

! For the sections, mapfcter, coriolis, earthcrv and buoyancy.

      if(mfcopt.ne.0.and.mfcopt.ne.1) then

        call destroy('chknlslv',8,'cont',201,'mfcopt        ',6,101,    &
     &               stat)

      end if

      if(coropt.ne.0.and.coropt.ne.1.and.coropt.ne.2) then

        call destroy('chknlslv',8,'cont',201,'coropt        ',6,101,    &
     &               stat)

      end if

      if(crvopt.ne.0.and.crvopt.ne.1) then

        call destroy('chknlslv',8,'cont',201,'crvopt        ',6,101,    &
     &               stat)

      end if

      if(buyopt.ne.0.and.buyopt.ne.1) then

        call destroy('chknlslv',8,'cont',201,'buyopt        ',6,101,    &
     &               stat)

      end if

! -----

! For the section, diabatic.

      if(diaopt.ne.0.and.diaopt.ne.1) then

        call destroy('chknlslv',8,'cont',201,'diaopt        ',6,101,    &
     &               stat)

      end if

      if(diaopt.eq.1) then

        if(advopt.ge.4) then

          call destroy('chknlslv',8,'cont',201,'advopt        ',6,101,  &
     &                 stat)

        end if

      end if

! -----

! For the section, ddamping.

      if(divopt.ne.0.and.divopt.ne.1.and.divopt.ne.2) then

        call destroy('chknlslv',8,'cont',201,'divopt        ',6,101,    &
     &               stat)

      end if

! -----

! For the section, cloudphy.

      if(cphopt.ne.-4.and.cphopt.ne.-3.and.cphopt.ne.-2                 &
     &  .and.cphopt.ne.0.and.cphopt.ne.1.and.cphopt.ne.2                &
     &  .and.cphopt.ne.3.and.cphopt.ne.4.and.cphopt.ne.11) then

        call destroy('chknlslv',8,'cont',201,'cphopt        ',6,101,    &
     &               stat)

      end if

      if(abs(cphopt).eq.0) then

        if(ngropt.ge.1) then

          call destroy('chknlslv',8,'cont',201,'ngropt        ',6,101,  &
     &                 stat)

        end if

      end if

      if(abs(cphopt).ge.1) then

        if(thresq.lt.eps) then

          call destroy('chknlslv',8,'cont',201,'thresq        ',6,101,  &
     &                 stat)

        end if

      end if

      if(abs(cphopt).lt.10) then

        if(abs(cphopt).ge.2) then

          if(savmem.eq.0) then

            if(zdim.lt.10) then

              call destroy('chknlslv',8,'cont',201,'zdim          ',4,  &
     &                     101,stat)

            end if

          end if

          if(haiopt.ne.0.and.haiopt.ne.1) then

            call destroy('chknlslv',8,'cont',201,'haiopt        ',6,101,&
     &                   stat)

          end if

        end if

      end if

      if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

        if(dtcmph(1).lt..001e0) then

          call destroy('chknlslv',8,'cont',201,'dtcmph(1)     ',9,101,  &
     &                 stat)

        end if

        if(dtcmph(2).lt..001e0) then

          call destroy('chknlslv',8,'cont',201,'dtcmph(2)     ',9,101,  &
     &                 stat)

        end if

        if(advopt.le.3) then

          if(dtcmph(1).ge..001e0) then

            if(mod(20*int(1.e2*(dtbig+.001e0)),                         &
     &        int(1.e3*(dtcmph(1)+.0001e0))).ne.0) then

              call destroy('chknlslv',8,'cont',201,'dtbig         ',5,  &
     &                     101,stat)

              call destroy('chknlslv',8,'cont',201,'dtcmph(1)     ',9,  &
     &                     101,stat)

            end if

          end if

          if(dtcmph(2).ge..001e0) then

            if(mod(20*int(1.e2*(dtbig+.001e0)),                         &
     &        int(1.e3*(dtcmph(2)+.0001e0))).ne.0) then

              call destroy('chknlslv',8,'cont',201,'dtbig         ',5,  &
     &                     101,stat)

              call destroy('chknlslv',8,'cont',201,'dtcmph(2)     ',9,  &
     &                     101,stat)

            end if

          end if

        else

          if(dtcmph(1).ge..001e0) then

            if(mod(10*int(1.e2*(dtbig+.001e0)),                         &
     &        int(1.e3*(dtcmph(1)+.0001e0))).ne.0) then

              call destroy('chknlslv',8,'cont',201,'dtbig         ',5,  &
     &                     101,stat)

              call destroy('chknlslv',8,'cont',201,'dtcmph(1)     ',9,  &
     &                     101,stat)

            end if

          end if

          if(dtcmph(2).ge..001e0) then

            if(mod(10*int(1.e2*(dtbig+.001e0)),                         &
     &        int(1.e3*(dtcmph(2)+.0001e0))).ne.0) then

              call destroy('chknlslv',8,'cont',201,'dtbig         ',5,  &
     &                     101,stat)

              call destroy('chknlslv',8,'cont',201,'dtcmph(2)     ',9,  &
     &                     101,stat)

            end if

          end if

        end if

        if(ncbinw.lt.2) then

          call destroy('chknlslv',8,'cont',201,'ncbinw        ',6,101,  &
     &                 stat)

        end if

        if(bbinw.le.1.e0) then

          call destroy('chknlslv',8,'cont',201,'bbinw         ',5,101,  &
     &                 stat)

        end if

        if(sbinw.le.0.e0) then

          call destroy('chknlslv',8,'cont',201,'sbinw         ',5,101,  &
     &                 stat)

        end if

      end if

      if(abs(cphopt).gt.10) then

        if(nggopt.eq.1) then

          call destroy('chknlslv',8,'cont',201,'nggopt        ',6,101,  &
     &                 stat)

        end if

        if(exbopt.ge.1) then

          call destroy('chknlslv',8,'cont',201,'exbopt        ',6,101,  &
     &                 stat)

        end if

        if(mod(lspopt,10).eq.1) then

          call destroy('chknlslv',8,'cont',201,'lspopt        ',6,101,  &
     &                 stat)

        end if

        if(vspopt.eq.1) then

          call destroy('chknlslv',8,'cont',201,'vspopt        ',6,101,  &
     &                 stat)

        end if

        if(ngropt.ge.1) then

          call destroy('chknlslv',8,'cont',201,'ngropt        ',6,101,  &
     &                 stat)

        end if

      end if

      if(cphopt.lt.0) then

        if(qcgopt.ne.1.and.qcgopt.ne.2) then

          call destroy('chknlslv',8,'cont',201,'qcgopt        ',6,101,  &
     &                 stat)

        end if

        if(eledlt.lt..01e0) then

          call destroy('chknlslv',8,'cont',201,'eledlt        ',6,101,  &
     &                 stat)

        end if

        if(dtbig.ge..01e0) then

          if(mod(10*int(1.e2*(eledlt+.001e0)),                          &
     &      10*int(1.e2*(dtbig+.001e0))).ne.0) then

            call destroy('chknlslv',8,'cont',201,'eledlt        ',6,101,&
     &                   stat)

          end if

        end if

      end if

! -----

! For the section, asolproc.

      if(aslopt.ne.0.and.aslopt.ne.1) then

        call destroy('chknlslv',8,'cont',201,'aslopt        ',6,101,    &
     &               stat)

      end if

      if(aslopt.ge.1) then

        if(abs(cphopt).lt.2.or.abs(cphopt).ge.10) then

          call destroy('chknlslv',8,'cont',201,'cphopt        ',6,101,  &
     &                 stat)

        end if

      end if

! -----

! For the section, mixtrace.

      if(trkopt.ne.0.and.trkopt.ne.1.and.trkopt.ne.2) then

        call destroy('chknlslv',8,'cont',201,'trkopt        ',6,101,    &
     &               stat)

      end if

      if(trkopt.eq.1) then

        if(qt0opt.ne.1.and.qt0opt.ne.2                                  &
     &    .and.qt0opt.ne.3.and.qt0opt.ne.4.and.qt0opt.ne.5) then

          call destroy('chknlslv',8,'cont',201,'qt0opt        ',6,101,  &
     &                 stat)

        end if

        if(qt0opt.eq.1.or.qt0opt.eq.2) then

          if(qt0num.lt.0.or.qt0num.gt.256) then

            call destroy('chknlslv',8,'cont',201,'qt0num        ',6,101,&
     &                   stat)

          end if

        end if

      end if

      if(trkopt.eq.2) then

        if(qt0opt.ne.1.and.qt0opt.ne.2) then

          call destroy('chknlslv',8,'cont',201,'qt0opt        ',6,101,  &
     &                 stat)

        end if

        if(qt0num.lt.0.or.qt0num.gt.256) then

          call destroy('chknlslv',8,'cont',201,'qt0num        ',6,101,  &
     &                 stat)

        end if

        if(qtdt.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'qtdt          ',4,101,  &
     &                 stat)

        end if

        if(qtstr.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'qtstr         ',5,101,  &
     &                 stat)

        end if

        if(qtend.lt.0.e0) then

          call destroy('chknlslv',8,'cont',201,'qtend         ',5,101,  &
     &                 stat)

        end if

        if(qtend.lt.qtstr) then

          call destroy('chknlslv',8,'cont',201,'qtend         ',5,101,  &
     &                 stat)

        end if

        if(pname(1:ncpn).ne.'check') then
          qtstr=aint(qtstr+.1e0)
          qtend=aint(qtend+.1e0)
        end if

      end if

! -----

! For the section, turbulen.

      if(tubopt.ne.0.and.tubopt.ne.1                                    &
     &  .and.tubopt.ne.2.and.tubopt.ne.3) then

        call destroy('chknlslv',8,'cont',201,'tubopt        ',6,101,    &
     &               stat)

      end if

      if(tubopt.ge.1) then

        if(isoopt.ne.1.and.isoopt.ne.2) then

          call destroy('chknlslv',8,'cont',201,'isoopt        ',6,101,  &
     &                 stat)

        end if

      end if

! -----

! For the section, outfomat.

      if(dmpfmt.ne.1.and.dmpfmt.ne.2) then

        call destroy('chknlslv',8,'cont',201,'dmpfmt        ',6,101,    &
     &               stat)

      end if

      if(dmplev.ne.1.and.dmplev.ne.2.and.dmplev.ne.3                    &
     &  .and.dmplev.ne.11.and.dmplev.ne.12.and.dmplev.ne.13) then

        call destroy('chknlslv',8,'cont',201,'dmplev        ',6,101,    &
     &               stat)

      end if

      if(mod(dmplev,10).eq.2.or.mod(dmplev,10).eq.3) then

        if(zdim.lt.10) then

          call destroy('chknlslv',8,'cont',201,'zdim          ',4,101,  &
     &                 stat)

        end if

      end if

      if(dmpmon.ne.0.and.dmpmon.ne.1) then

        call destroy('chknlslv',8,'cont',201,'dmpmon        ',6,101,    &
     &               stat)

      end if

      ierr=0

      if(dmpvar(1:1).ne.'o'                                             &
     &  .and.dmpvar(1:1).ne.'-'.and.dmpvar(1:1).ne.'x') then

        ierr=ierr+1

      end if

      if(dmpvar(2:2).ne.'o'                                             &
     &  .and.dmpvar(2:2).ne.'-'.and.dmpvar(2:2).ne.'x') then

        ierr=ierr+1

      end if

      if(dmpvar(3:3).ne.'o'.and.dmpvar(3:3).ne.'x') then
        ierr=ierr+1
      end if

      if(dmpvar(4:4).ne.'o'                                             &
     &  .and.dmpvar(4:4).ne.'-'.and.dmpvar(4:4).ne.'x') then

        ierr=ierr+1

      end if

      if(dmpvar(5:5).ne.'o'                                             &
     &  .and.dmpvar(5:5).ne.'-'.and.dmpvar(5:5).ne.'x') then

        ierr=ierr+1

      end if

      if(dmpvar(6:6).ne.'o'                                             &
     &  .and.dmpvar(6:6).ne.'-'.and.dmpvar(6:6).ne.'x') then

        ierr=ierr+1

      end if

      if(dmpvar(7:7).ne.'o'.and.dmpvar(7:7).ne.'x') then
        ierr=ierr+1
      end if

      if(dmpvar(8:8).ne.'o'.and.dmpvar(8:8).ne.'x') then
        ierr=ierr+1
      end if

      if(dmpvar(9:9).ne.'o'.and.dmpvar(9:9).ne.'x') then
        ierr=ierr+1
      end if

      if(dmpvar(10:10).ne.'o'.and.dmpvar(10:10).ne.'x') then
        ierr=ierr+1
      end if

      if(dmpvar(11:11).ne.'o'.and.dmpvar(11:11).ne.'x') then
        ierr=ierr+1
      end if

      if(dmpvar(12:12).ne.'+'.and.dmpvar(12:12).ne.'o'                  &
     &  .and.dmpvar(12:12).ne.'-'.and.dmpvar(12:12).ne.'x') then

        ierr=ierr+1

      end if

      if(dmpvar(13:13).ne.'o'.and.dmpvar(13:13).ne.'x') then
        ierr=ierr+1
      end if

      if(dmpvar(14:14).ne.'o'.and.dmpvar(14:14).ne.'x') then
        ierr=ierr+1
      end if

      if(dmpvar(15:15).ne.'o'.and.dmpvar(15:15).ne.'x') then
        ierr=ierr+1
      end if

      if(ierr.ne.0) then

        call destroy('chknlslv',8,'cont',201,'dmpvar        ',6,101,    &
     &               stat)

      end if

      if(dmpitv.lt.1.e0) then

        call destroy('chknlslv',8,'cont',201,'dmpitv        ',6,101,    &
     &               stat)

      end if

      if(dmpmon.eq.1) then

        if(monitv.lt.1.e0) then

          call destroy('chknlslv',8,'cont',201,'monitv        ',6,101,  &
     &                 stat)

        end if

      end if

      if(resopt.ne.0.and.resopt.ne.1) then

        call destroy('chknlslv',8,'cont',201,'resopt        ',6,101,    &
     &               stat)

      end if

      if(resopt.eq.1) then

        if(resitv.lt.1.e0) then

          call destroy('chknlslv',8,'cont',201,'resitv        ',6,101,  &
     &                 stat)

        end if

      end if

      if(mxnopt.ne.0.and.mxnopt.ne.1) then

        call destroy('chknlslv',8,'cont',201,'mxnopt        ',6,101,    &
     &               stat)

      end if

      if(mxnopt.eq.1) then

        ierr=0

        if(mxnvar(1:1).ne.'o'.and.mxnvar(1:1).ne.'x') then
          ierr=ierr+1
        end if

        if(mxnvar(2:2).ne.'o'.and.mxnvar(2:2).ne.'x') then
          ierr=ierr+1
        end if

        if(mxnvar(3:3).ne.'o'.and.mxnvar(3:3).ne.'x') then
          ierr=ierr+1
        end if

        if(mxnvar(4:4).ne.'o'.and.mxnvar(4:4).ne.'x') then
          ierr=ierr+1
        end if

        if(mxnvar(5:5).ne.'o'.and.mxnvar(5:5).ne.'x') then
          ierr=ierr+1
        end if

        if(mxnvar(6:6).ne.'o'.and.mxnvar(6:6).ne.'x') then
          ierr=ierr+1
        end if

        if(mxnvar(7:7).ne.'o'.and.mxnvar(7:7).ne.'x') then
          ierr=ierr+1
        end if

        if(mxnvar(8:8).ne.'o'.and.mxnvar(8:8).ne.'x') then
          ierr=ierr+1
        end if

        if(mxnvar(9:9).ne.'o'.and.mxnvar(9:9).ne.'x') then
          ierr=ierr+1
        end if

        if(mxnvar(10:10).ne.'o'.and.mxnvar(10:10).ne.'x') then
          ierr=ierr+1
        end if

        if(mxnvar(11:11).ne.'o'.and.mxnvar(11:11).ne.'x') then
          ierr=ierr+1
        end if

        if(mxnvar(12:12).ne.'o'.and.mxnvar(12:12).ne.'x') then
          ierr=ierr+1
        end if

        if(ierr.ne.0) then

          call destroy('chknlslv',8,'cont',201,'mxnvar        ',6,101,  &
     &                 stat)

        end if

        if(mxnitv.lt.1.e0) then

          call destroy('chknlslv',8,'cont',201,'mxnitv        ',6,101,  &
     &                 stat)

        end if

      end if

      if(pname(1:ncpn).ne.'check') then
        dmpitv=aint(dmpitv+.1e0)
        monitv=aint(monitv+.1e0)
        resitv=aint(resitv+.1e0)
        mxnitv=aint(mxnitv+.1e0)
      end if

! -----

! For the sections, project_rdr, gridset_rdr and datconf_rdr.

      if(ngropt.eq.12) then

        if(mpopt_rdr.ne.0.and.                                          &
     &     mpopt_rdr.ne.1.and.mpopt_rdr.ne.2.and.                       &
     &     mpopt_rdr.ne.3.and.mpopt_rdr.ne.4.and.                       &
     &     mpopt_rdr.ne.10.and.mpopt_rdr.ne.13) then

          call destroy('chknlslv',8,'cont',201,'mpopt_rdr     ',9,101,  &
     &                 stat)

        end if

        if(mpopt.ge.10) then

          if(mpopt_rdr.lt.10) then

            call destroy('chknlslv',8,'cont',201,'mpopt_rdr     ',9,101,&
     &                   stat)

          end if

        end if

        if(mpopt_rdr.eq.1.or.mpopt_rdr.eq.2                             &
     &    .or.mpopt_rdr.eq.3.or.mpopt_rdr.eq.13) then

          if(nspol_rdr.ne.1.and.nspol_rdr.ne.-1) then

            call destroy('chknlslv',8,'cont',201,'nspol_rdr     ',9,101,&
     &                   stat)

          end if

          if(tlat1_rdr.lt.-90.e0.or.tlat1_rdr.gt.90.e0) then

            call destroy('chknlslv',8,'cont',201,'tlat1_rdr     ',9,101,&
     &                   stat)

          end if

        end if

        if(mpopt_rdr.eq.2) then

          if(tlat2_rdr.lt.-90.e0.or.tlat2_rdr.gt.90.e0) then

            call destroy('chknlslv',8,'cont',201,'tlat2_rdr     ',9,101,&
     &                   stat)

          end if

        end if

        if(mpopt_rdr.eq.1.or.mpopt_rdr.eq.2.or.mpopt_rdr.eq.4) then

          if(tlon_rdr.lt.-180.e0.or.tlon_rdr.gt.180.e0) then

            call destroy('chknlslv',8,'cont',201,'tlon_rdr      ',8,101,&
     &                   stat)

          end if

        end if

        if(ulat_rdr.lt.-90.e0.or.ulat_rdr.gt.90.e0) then

          call destroy('chknlslv',8,'cont',201,'ulat_rdr      ',8,101,  &
     &                 stat)

        end if

        if(ulon_rdr.lt.-180.e0.or.ulon_rdr.gt.180.e0) then

          call destroy('chknlslv',8,'cont',201,'ulon_rdr      ',8,101,  &
     &                 stat)

        end if

        if(xdim_rdr.lt.2) then

          call destroy('chknlslv',8,'cont',201,'xdim_rdr      ',8,101,  &
     &                 stat)

        end if

        if(ydim_rdr.lt.2) then

          call destroy('chknlslv',8,'cont',201,'ydim_rdr      ',8,101,  &
     &                 stat)

        end if

        if(zdim_rdr.lt.2) then

          call destroy('chknlslv',8,'cont',201,'zdim_rdr      ',8,101,  &
     &                 stat)

        end if

        if(mpopt_rdr.lt.10) then

          if(dx_rdr.lt.eps) then

            call destroy('chknlslv',8,'cont',201,'dx_rdr        ',6,101,&
     &                   stat)

          else

            if(pname(1:ncpn).ne.'check') then
              dxiv_rdr=1.e0/dx_rdr
            end if

          end if

        else

          if(pname(1:ncpn).ne.'check') then

            if(xdim_rdr.ge.4) then

              if(mpopt_rdr.eq.10) then

                fmsg=fmsg+1

                dx_rdr=360.e0/real(xdim_rdr)

                call outstd15('dx_rdr',6,fmsg,2,dx_rdr)

              end if

              if(mpopt_rdr.eq.13) then

                fmsg=fmsg+1

                dx_rdr=2.e0*cc*rearth/real(xdim_rdr)

                call outstd15('dx_rdr',6,fmsg,1,dx_rdr)

              end if

              dxiv_rdr=1.e0/dx_rdr

            end if

          end if

        end if

        if(dy_rdr.lt.eps) then

          call destroy('chknlslv',8,'cont',201,'dy_rdr        ',6,101,  &
     &                 stat)

        else

          if(pname(1:ncpn).ne.'check') then
            dyiv_rdr=1.e0/dy_rdr
          end if

        end if

        if(rotopt_rdr.ne.0.and.rotopt_rdr.ne.1) then

          call destroy('chknlslv',8,'cont',201,'rotopt_rdr    ',10,101, &
     &                 stat)

        end if

        if(datype_rdr(1:1).ne.'r'.and.datype_rdr(1:1).ne.'m') then

          call destroy('chknlslv',8,'cont',201,'datype_rdr    ',10,101, &
     &                 stat)

        end if

      end if

! -----

!! -----

! Set the constant value.

      if(pname(1:ncpn).ne.'check') then

        write(gpvvar(9:9),'(a1)') 'o'
        write(exbvar(9:9),'(a1)') 'x'

      end if

! -----

! Check the parameters of parallelizing.

      if(npe.lt.0) then
        npe=numpe
      end if

      if(npe.ne.numpe) then

        stat=stat+1

        call destroy('chknlslv',8,'cont',11,'              ',14,101,    &
     &               stat)

      end if

      if(xgroup.eq.0.or.xsub.eq.0.or.ygroup.eq.0.or.ysub.eq.0) then

        stat=stat+1

        call destroy('chknlslv',8,'cont',11,'              ',14,101,    &
     &               stat)

      else

        if(mod((xdim-3),xgroup*xsub).ne.0                               &
     &    .or.mod((ydim-3),ygroup*ysub).ne.0) then

          stat=stat+1

          call destroy('chknlslv',8,'cont',11,'              ',14,101,  &
     &                 stat)

        end if

      end if

      if(numpe.lt.xsub*ysub.or.numpe.gt.xsub*ysub*xgroup*ygroup) then

        stat=stat+1

        call destroy('chknlslv',8,'cont',11,'              ',14,101,    &
     &               stat)

      end if

! -----

      end subroutine s_chknlslv

!-----7--------------------------------------------------------------7--

      end module m_chknlslv
