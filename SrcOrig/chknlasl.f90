!***********************************************************************
      module m_chknlasl
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/09/22
!     Modification: 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     check the namelist variables for asldata.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
      use m_commath
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

      public :: chknlasl, s_chknlasl

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface chknlasl

        module procedure s_chknlasl

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic aint
      intrinsic int
      intrinsic mod
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_chknlasl(pname,ncpn,stat)
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

      integer fmsg     ! Control flag of message type for standard i/o

      integer ncspc    ! Number of space character of exprim and prvres

      real r4scl       ! Optional real scalar variable

!-----7--------------------------------------------------------------7--

! Set the word length for direct access file.

      r4scl=0.e0

      inquire(iolength=wlngth) r4scl

! -----

!! Check the namelist variables.

! Initialize the controler.

      stat=0

      fmsg=0

! -----

! For the section, runame.

      call numchar(exprim,1,ncexp,ncspc)

      if(ncexp.gt.64.or.ncspc.ne.0) then

        call destroy('chknlasl',8,'cont',201,'exprim        ',6,101,    &
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

        call destroy('chknlasl',8,'cont',201,'crsdir        ',6,101,    &
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

        call destroy('chknlasl',8,'cont',201,'datdir        ',6,101,    &
     &               stat)

      end if

      if(pname(1:ncpn).ne.'check') then

        call outstd14('chknlasl',8,crsdir,datdir,nccrs,ncdat)

      end if

! -----

! For the section, dimset.

      if(xdim.le.4) then

        call destroy('chknlasl',8,'cont',201,'xdim          ',4,101,    &
     &               stat)

      end if

      if(mpopt.eq.5) then

        if(ydim.ne.4) then

          call destroy('chknlasl',8,'cont',201,'ydim          ',4,101,  &
     &                 stat)

        end if

      else

        if(ydim.le.4) then

          call destroy('chknlasl',8,'cont',201,'ydim          ',4,101,  &
     &                 stat)

        end if

      end if

      if(zdim.lt.7) then

        call destroy('chknlasl',8,'cont',201,'zdim          ',4,101,    &
     &               stat)

      end if

! -----

! For the section, project.

      if(mpopt.ne.0.and.mpopt.ne.1.and.mpopt.ne.2.and.mpopt.ne.3.and.   &
     &   mpopt.ne.4.and.mpopt.ne.5.and.mpopt.ne.10.and.mpopt.ne.13) then

        call destroy('chknlasl',8,'cont',201,'mpopt         ',5,101,    &
     &               stat)

      end if

      if(mpopt.eq.1.or.mpopt.eq.2.or.mpopt.eq.3.or.mpopt.eq.13) then

        if(nspol.ne.1.and.nspol.ne.-1) then

          call destroy('chknlasl',8,'cont',201,'nspol         ',5,101,  &
     &                 stat)

        end if

        if(tlat1.lt.-90.e0.or.tlat1.gt.90.e0) then

          call destroy('chknlasl',8,'cont',201,'tlat1         ',5,101,  &
     &                 stat)

        end if

      end if

      if(mpopt.eq.2) then

        if(tlat2.lt.-90.e0.or.tlat2.gt.90.e0) then

          call destroy('chknlasl',8,'cont',201,'tlat2         ',5,101,  &
     &                 stat)

        end if

      end if

      if(mpopt.eq.1.or.mpopt.eq.2.or.mpopt.eq.4) then

        if(tlon.lt.-180.e0.or.tlon.gt.180.e0) then

          call destroy('chknlasl',8,'cont',201,'tlon          ',4,101,  &
     &                 stat)

        end if

      end if

      if(ulat.lt.-90.e0.or.ulat.gt.90.e0) then

        call destroy('chknlasl',8,'cont',201,'ulat          ',4,101,    &
     &               stat)

      end if

      if(ulon.lt.-180.e0.or.ulon.gt.180.e0) then

        call destroy('chknlasl',8,'cont',201,'ulon          ',4,101,    &
     &               stat)

      end if

      if(mpopt.eq.5) then

        if(disr.lt.eps) then

          call destroy('chknlasl',8,'cont',201,'disr          ',4,101,  &
     &                 stat)

        end if

      end if

! -----

! For the section, gridset.

      if(mpopt.lt.10) then

        if(dx.lt.eps) then

          call destroy('chknlasl',8,'cont',201,'dx            ',2,101,  &
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

          call destroy('chknlasl',8,'cont',201,'dx            ',2,101,  &
     &                 stat)

        end if

      end if

      if(dy.lt.eps) then

        call destroy('chknlasl',8,'cont',201,'dy            ',2,101,    &
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

        call destroy('chknlasl',8,'cont',201,'dz            ',2,101,    &
     &               stat)

      else

        if(pname(1:ncpn).ne.'check') then
          dziv=1.e0/dz
        end if

      end if

! -----

! For the section, gridsth.

      if(zsfc.lt.-1000.e0) then

        call destroy('chknlasl',8,'cont',201,'zsfc          ',4,101,    &
     &               stat)

      end if

      if(zflat.lt.zsfc) then

        call destroy('chknlasl',8,'cont',201,'zflat         ',5,101,    &
     &               stat)

      end if

      if(sthopt.ne.0.and.sthopt.ne.1.and.sthopt.ne.2) then

        call destroy('chknlasl',8,'cont',201,'sthopt        ',6,101,    &
     &               stat)

      end if

      if(sthopt.eq.0) then

        if(pname(1:ncpn).ne.'check') then
          dzmin=dz
        end if

      else

        if(dzmin.lt.eps.or.dzmin.gt.dz) then

          call destroy('chknlasl',8,'cont',201,'dzmin         ',5,101,  &
     &                 stat)

        end if

        if(layer1.lt.zsfc) then

          call destroy('chknlasl',8,'cont',201,'layer1        ',6,101,  &
     &                 stat)

        end if

        if(layer2.lt.layer1) then

          call destroy('chknlasl',8,'cont',201,'layer2        ',6,101,  &
     &                 stat)

        end if

      end if

! -----

! For the section, terrain.

      if(trnopt.ne.0.and.trnopt.ne.1.and.trnopt.ne.2) then

        call destroy('chknlasl',8,'cont',201,'trnopt        ',6,101,    &
     &               stat)

      end if

      if(trnopt.eq.0.or.trnopt.eq.1) then

        if(mnthgh(1).lt.0.e0) then

          call destroy('chknlasl',8,'cont',201,'mnthgh(1)     ',9,101,  &
     &                 stat)

        end if

      end if

      if(trnopt.eq.1) then

        if(mntwx.lt.0.e0) then

          call destroy('chknlasl',8,'cont',201,'mntwx         ',5,101,  &
     &                 stat)

        end if

        if(mntwy.lt.0.e0) then

          call destroy('chknlasl',8,'cont',201,'mntwy         ',5,101,  &
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

        call destroy('chknlasl',8,'cont',201,'stime         ',5,101,    &
     &               stat)

      end if

      if(etime.lt.0.e0) then

        call destroy('chknlasl',8,'cont',201,'etime         ',5,101,    &
     &               stat)

      end if

      if(etime.lt.stime) then

        call destroy('chknlasl',8,'cont',201,'etime         ',5,101,    &
     &               stat)

      end if

      if(1000_i8*int(stime+.1e0,i8).gt.1000_i8*int(etime+.1e0,i8)) then

        call destroy('chknlasl',8,'cont',201,'stime         ',5,101,    &
     &               stat)

      end if

! -----

! For the section, boundry.

      if(wbc.ne.-1.and.wbc.ne.1.and.wbc.ne.2.and.wbc.ne.3.and.          &
     &   wbc.ne.4.and.wbc.ne.5.and.wbc.ne.6.and.wbc.ne.7.and.           &
     &   wbc.ne.14.and.wbc.ne.15.and.wbc.ne.16) then

        call destroy('chknlasl',8,'cont',201,'wbc           ',3,101,    &
     &               stat)

      end if

      if(ebc.ne.-1.and.ebc.ne.1.and.ebc.ne.2.and.ebc.ne.3.and.          &
     &   ebc.ne.4.and.ebc.ne.5.and.ebc.ne.6.and.ebc.ne.7.and.           &
     &   ebc.ne.14.and.ebc.ne.15.and.ebc.ne.16) then

        call destroy('chknlasl',8,'cont',201,'ebc           ',3,101,    &
     &               stat)

      end if

      if(sbc.ne.-1.and.sbc.ne.1.and.sbc.ne.2.and.sbc.ne.3.and.          &
     &   sbc.ne.4.and.sbc.ne.5.and.sbc.ne.6.and.sbc.ne.7.and.           &
     &   sbc.ne.14.and.sbc.ne.15.and.sbc.ne.16) then

        call destroy('chknlasl',8,'cont',201,'sbc           ',3,101,    &
     &               stat)

      end if

      if(nbc.ne.-1.and.nbc.ne.1.and.nbc.ne.2.and.nbc.ne.3.and.          &
     &   nbc.ne.4.and.nbc.ne.5.and.nbc.ne.6.and.nbc.ne.7.and.           &
     &   nbc.ne.14.and.nbc.ne.15.and.nbc.ne.16) then

        call destroy('chknlasl',8,'cont',201,'nbc           ',3,101,    &
     &               stat)

      end if

      if(bbc.ne.2.and.bbc.ne.3) then

        call destroy('chknlasl',8,'cont',201,'bbc           ',3,101,    &
     &               stat)

      end if

      if(tbc.ne.2.and.tbc.ne.3) then

        call destroy('chknlasl',8,'cont',201,'tbc           ',3,101,    &
     &               stat)

      end if

      if(mpopt.eq.1.or.mpopt.eq.2.or.mpopt.eq.3) then

        if(abs(wbc).eq.1.or.abs(ebc).eq.1) then

          call destroy('chknlasl',8,'cont',201,'wbc           ',3,101,  &
     &                 stat)

          call destroy('chknlasl',8,'cont',201,'ebc           ',3,101,  &
     &                 stat)

        end if

      end if

      if(mpopt.eq.1.or.mpopt.eq.2.or.mpopt.eq.3                         &
     &  .or.mpopt.eq.10.or.mpopt.eq.13) then

        if(abs(sbc).eq.1.or.abs(nbc).eq.1) then

          call destroy('chknlasl',8,'cont',201,'sbc           ',3,101,  &
     &                 stat)

          call destroy('chknlasl',8,'cont',201,'nbc           ',3,101,  &
     &                 stat)

        end if

      end if

      if(mpopt.eq.5) then

        if(abs(wbc).eq.1.or.abs(ebc).eq.1) then

          call destroy('chknlasl',8,'cont',201,'wbc           ',3,101,  &
     &                 stat)

          call destroy('chknlasl',8,'cont',201,'ebc           ',3,101,  &
     &                 stat)

        end if

        if(abs(sbc).ne.1.or.abs(nbc).ne.1) then

          call destroy('chknlasl',8,'cont',201,'sbc           ',3,101,  &
     &                 stat)

          call destroy('chknlasl',8,'cont',201,'nbc           ',3,101,  &
     &                 stat)

        end if

      end if

      if(abs(wbc).eq.1.or.abs(ebc).eq.1) then

        if(wbc.ne.ebc) then

          call destroy('chknlasl',8,'cont',201,'wbc           ',3,101,  &
     &                 stat)

          call destroy('chknlasl',8,'cont',201,'ebc           ',3,101,  &
     &                 stat)

        end if

      end if

      if(xdim.eq.4) then

        if(wbc.ne.ebc) then

          call destroy('chknlasl',8,'cont',201,'wbc           ',3,101,  &
     &                 stat)

          call destroy('chknlasl',8,'cont',201,'ebc           ',3,101,  &
     &                 stat)

        end if

      end if

      if(abs(sbc).eq.1.or.abs(nbc).eq.1) then

        if(sbc.ne.nbc) then

          call destroy('chknlasl',8,'cont',201,'sbc           ',3,101,  &
     &                 stat)

          call destroy('chknlasl',8,'cont',201,'nbc           ',3,101,  &
     &                 stat)

        end if

      end if

      if(ydim.eq.4) then

        if(sbc.ne.nbc) then

          call destroy('chknlasl',8,'cont',201,'sbc           ',3,101,  &
     &                 stat)

          call destroy('chknlasl',8,'cont',201,'nbc           ',3,101,  &
     &                 stat)

        end if

      end if

      if(numpe.ne.xsub*ysub*xgroup*ygroup) then

        if(wbc.eq.1.or.ebc.eq.1) then

          call destroy('chknlasl',8,'cont',201,'wbc           ',3,101,  &
     &                 stat)

          call destroy('chknlasl',8,'cont',201,'ebc           ',3,101,  &
     &                 stat)

        end if

        if(sbc.eq.1.or.nbc.eq.1) then

          call destroy('chknlasl',8,'cont',201,'sbc           ',3,101,  &
     &                 stat)

          call destroy('chknlasl',8,'cont',201,'nbc           ',3,101,  &
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

      if(gsmopt.ne.0.and.gsmopt.ne.1) then

        call destroy('chknlasl',8,'cont',201,'gsmopt        ',6,101,    &
     &               stat)

      end if

      if(nggopt.ne.0.and.nggopt.ne.1) then

        call destroy('chknlasl',8,'cont',201,'nggopt        ',6,101,    &
     &               stat)

      end if

      if(exbopt.ne.0.and.exbopt.ne.1.and.exbopt.ne.2                    &
     &  .and.exbopt.ne.11.and.exbopt.ne.12) then

        call destroy('chknlasl',8,'cont',201,'exbopt        ',6,101,    &
     &               stat)

      end if

      if(mod(exbopt,10).eq.2) then

        if(exbwid.lt.2) then

          call destroy('chknlasl',8,'cont',201,'exbwid        ',6,101,  &
     &                 stat)

        end if

        if((xdim.ne.4.and.2*exbwid.gt.xdim)                             &
     &    .or.(ydim.ne.4.and.2*exbwid.gt.ydim)) then

          call destroy('chknlasl',8,'cont',201,'exbwid        ',6,101,  &
     &                 stat)

        end if

      end if

      if(lspopt.ne.0.and.lspopt.ne.1.and.lspopt.ne.2                    &
     &  .and.lspopt.ne.11.and.lspopt.ne.12) then

        call destroy('chknlasl',8,'cont',201,'lspopt        ',6,101,    &
     &               stat)

      end if

      if(vspopt.ne.0.and.vspopt.ne.1.and.vspopt.ne.2) then

        call destroy('chknlasl',8,'cont',201,'vspopt        ',6,101,    &
     &               stat)

      end if

      if(nggopt.eq.1.or.exbopt.ge.1                                     &
     &  .or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

        if(int(aslitv+.1e0).lt.60) then

          call destroy('chknlasl',8,'cont',201,'aslitv        ',6,101,  &
     &                 stat)

        end if

      end if

      if(pname(1:ncpn).ne.'check') then
        aslitv=aint(aslitv+.1e0)
      end if

      if(nggopt.eq.0.and.exbopt.eq.0                                    &
     &  .and.lspopt.eq.0.and.vspopt.eq.0) then

        call destroy('chknlasl',8,'cont',201,'nggopt        ',6,101,    &
     &               stat)

        call destroy('chknlasl',8,'cont',201,'exbopt        ',6,101,    &
     &               stat)

        call destroy('chknlasl',8,'cont',201,'lspopt        ',6,101,    &
     &               stat)

        call destroy('chknlasl',8,'cont',201,'vspopt        ',6,101,    &
     &               stat)

      end if

! -----

! For the section, asolproc.

      if(aslopt.eq.0) then

        call destroy('chknlslv',8,'cont',201,'aslopt        ',6,101,    &
     &               stat)

      end if

! -----

! For the sections, project_asl, gridset_asl and datconf_asl.

      if(mpopt_asl.ne.0.and.                                            &
     &   mpopt_asl.ne.1.and.mpopt_asl.ne.2.and.                         &
     &   mpopt_asl.ne.3.and.mpopt_asl.ne.4.and.                         &
     &   mpopt_asl.ne.10.and.mpopt_asl.ne.13) then

        call destroy('chknlasl',8,'cont',201,'mpopt_asl     ',9,101,    &
     &               stat)

      end if

      if(mpopt.ge.10) then

        if(mpopt_asl.lt.10) then

          call destroy('chknlasl',8,'cont',201,'mpopt_asl     ',9,101,  &
     &                 stat)

        end if

      end if

      if(mpopt_asl.eq.1.or.mpopt_asl.eq.2                               &
     &  .or.mpopt_asl.eq.3.or.mpopt_asl.eq.13) then

        if(nspol_asl.ne.1.and.nspol_asl.ne.-1) then

          call destroy('chknlasl',8,'cont',201,'nspol_asl     ',9,101,  &
     &                 stat)

        end if

        if(tlat1_asl.lt.-90.e0.or.tlat1_asl.gt.90.e0) then

          call destroy('chknlasl',8,'cont',201,'tlat1_asl     ',9,101,  &
     &                 stat)

        end if

      end if

      if(mpopt_asl.eq.2) then

        if(tlat2_asl.lt.-90.e0.or.tlat2_asl.gt.90.e0) then

          call destroy('chknlasl',8,'cont',201,'tlat2_asl     ',9,101,  &
     &                 stat)

        end if

      end if

      if(mpopt_asl.eq.1.or.mpopt_asl.eq.2.or.mpopt_asl.eq.4) then

        if(tlon_asl.lt.-180.e0.or.tlon_asl.gt.180.e0) then

          call destroy('chknlasl',8,'cont',201,'tlon_asl      ',8,101,  &
     &                 stat)

        end if

      end if

      if(ulat_asl.lt.-90.e0.or.ulat_asl.gt.90.e0) then

        call destroy('chknlasl',8,'cont',201,'ulat_asl      ',8,101,    &
     &               stat)

      end if

      if(ulon_asl.lt.-180.e0.or.ulon_asl.gt.180.e0) then

        call destroy('chknlasl',8,'cont',201,'ulon_asl      ',8,101,    &
     &               stat)

      end if

      if(xdim_asl.lt.2) then

        call destroy('chknlasl',8,'cont',201,'xdim_asl      ',8,101,    &
     &               stat)

      end if

      if(ydim_asl.lt.2) then

        call destroy('chknlasl',8,'cont',201,'ydim_asl      ',8,101,    &
     &               stat)

      end if

      if(zdim_asl.lt.2) then

        call destroy('chknlasl',8,'cont',201,'zdim_asl      ',8,101,    &
     &               stat)

      end if

      if(mpopt_asl.lt.10) then

        if(dx_asl.lt.eps) then

          call destroy('chknlasl',8,'cont',201,'dx_asl        ',6,101,  &
     &                 stat)

        else

          if(pname(1:ncpn).ne.'check') then
            dxiv_asl=1.e0/dx_asl
          end if

        end if

      else

        if(pname(1:ncpn).ne.'check') then

          if(xdim_asl.ge.4) then

            if(mpopt_asl.eq.10) then

              fmsg=fmsg+1

              dx_asl=360.e0/real(xdim_asl)

              call outstd15('dx_asl',6,fmsg,2,dx_asl)

            end if

            if(mpopt_asl.eq.13) then

              fmsg=fmsg+1

              dx_asl=2.e0*cc*rearth/real(xdim_asl)

              call outstd15('dx_asl',6,fmsg,1,dx_asl)

            end if

            dxiv_asl=1.e0/dx_asl

          end if

        end if

      end if

      if(dy_asl.lt.eps) then

        call destroy('chknlasl',8,'cont',201,'dy_asl        ',6,101,    &
     &               stat)

      else

        if(pname(1:ncpn).ne.'check') then
          dyiv_asl=1.e0/dy_asl
        end if

      end if

      if(intopt_asl.ne.1.and.intopt_asl.ne.2) then

        call destroy('chknlasl',8,'cont',201,'intopt_asl    ',10,101,   &
     &               stat)

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

      if(xgroup.eq.0.or.xsub.eq.0.or.ygroup.eq.0.or.ysub.eq.0) then

        stat=stat+1

        call destroy('chknlasl',8,'cont',11,'              ',14,101,    &
     &               stat)

      else

        if(mod((xdim-3),xgroup*xsub).ne.0                               &
     &    .or.mod((ydim-3),ygroup*ysub).ne.0) then

          stat=stat+1

          call destroy('chknlasl',8,'cont',11,'              ',14,101,  &
     &                 stat)

        end if

      end if

      if(numpe.lt.xsub*ysub.or.numpe.gt.xsub*ysub*xgroup*ygroup) then

        stat=stat+1

        call destroy('chknlasl',8,'cont',11,'              ',14,101,    &
     &               stat)

      end if

! -----

      end subroutine s_chknlasl

!-----7--------------------------------------------------------------7--

      end module m_chknlasl
