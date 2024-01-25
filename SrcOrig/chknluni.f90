!***********************************************************************
      module m_chknluni
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/01/20
!     Modification: 2007/04/11, 2007/07/30, 2007/09/28, 2007/11/26,
!                   2008/01/11, 2008/03/12, 2008/04/17, 2008/05/02,
!                   2008/08/25, 2008/10/10, 2009/01/05, 2009/01/30,
!                   2009/02/27, 2009/03/23, 2011/08/09, 2011/08/18,
!                   2011/09/22, 2011/11/10, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     check the namelist variables for unite.

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

      public :: chknluni, s_chknluni

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface chknluni

        module procedure s_chknluni

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
      subroutine s_chknluni(pname,ncpn,stat)
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

! For the section, runame.

      call numchar(exprim,1,ncexp,ncspc)

      if(ncexp.gt.64.or.ncspc.ne.0) then

        call destroy('chknluni',8,'cont',201,'exprim        ',6,101,    &
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

        call destroy('chknluni',8,'cont',201,'crsdir        ',6,101,    &
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

        call destroy('chknluni',8,'cont',201,'datdir        ',6,101,    &
     &               stat)

      end if

      if(pname(1:ncpn).ne.'check') then

        call outstd14('chknluni',8,crsdir,datdir,nccrs,ncdat)

      end if

! -----

! For the section, dimset.

      if(xdim.lt.4) then

        call destroy('chknluni',8,'cont',201,'xdim          ',4,101,    &
     &               stat)

      end if

      if(mpopt.eq.5) then

        if(ydim.ne.4) then

          call destroy('chknluni',8,'cont',201,'ydim          ',4,101,  &
     &                 stat)

        end if

      else

        if(ydim.lt.4) then

          call destroy('chknluni',8,'cont',201,'ydim          ',4,101,  &
     &                 stat)

        end if

      end if

      if(zdim.lt.7) then

        call destroy('chknluni',8,'cont',201,'zdim          ',4,101,    &
     &               stat)

      end if

! -----

! For the section, project.

      if(mpopt.ne.0.and.mpopt.ne.1.and.mpopt.ne.2.and.mpopt.ne.3.and.   &
     &   mpopt.ne.4.and.mpopt.ne.5.and.mpopt.ne.10.and.mpopt.ne.13) then

        call destroy('chknluni',8,'cont',201,'mpopt         ',5,101,    &
     &               stat)

      end if

      if(mpopt.eq.1.or.mpopt.eq.2.or.mpopt.eq.3.or.mpopt.eq.13) then

        if(nspol.ne.1.and.nspol.ne.-1) then

          call destroy('chknluni',8,'cont',201,'nspol         ',5,101,  &
     &                 stat)

        end if

        if(tlat1.lt.-90.e0.or.tlat1.gt.90.e0) then

          call destroy('chknluni',8,'cont',201,'tlat1         ',5,101,  &
     &                 stat)

        end if

      end if

      if(mpopt.eq.2) then

        if(tlat2.lt.-90.e0.or.tlat2.gt.90.e0) then

          call destroy('chknluni',8,'cont',201,'tlat2         ',5,101,  &
     &                 stat)

        end if

      end if

      if(mpopt.eq.1.or.mpopt.eq.2.or.mpopt.eq.4) then

        if(tlon.lt.-180.e0.or.tlon.gt.180.e0) then

          call destroy('chknluni',8,'cont',201,'tlon          ',4,101,  &
     &                 stat)

        end if

      end if

      if(ulat.lt.-90.e0.or.ulat.gt.90.e0) then

        call destroy('chknluni',8,'cont',201,'ulat          ',4,101,    &
     &               stat)

      end if

      if(ulon.lt.-180.e0.or.ulon.gt.180.e0) then

        call destroy('chknluni',8,'cont',201,'ulon          ',4,101,    &
     &               stat)

      end if

      if(mpopt.eq.5) then

        if(disr.lt.eps) then

          call destroy('chknluni',8,'cont',201,'disr          ',4,101,  &
     &                 stat)

        end if

      end if

! -----

! For the section, gridset.

      if(mpopt.lt.10) then

        if(dx.lt.eps) then

          call destroy('chknluni',8,'cont',201,'dx            ',2,101,  &
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

          call destroy('chknluni',8,'cont',201,'dx            ',2,101,  &
     &                 stat)

        end if

      end if

      if(dy.lt.eps) then

        call destroy('chknluni',8,'cont',201,'dy            ',2,101,    &
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

        call destroy('chknluni',8,'cont',201,'dz            ',2,101,    &
     &               stat)

      else

        if(pname(1:ncpn).ne.'check') then
          dziv=1.e0/dz
        end if

      end if

! -----

! For the section, gridsth.

      if(zsfc.lt.-1000.e0) then

        call destroy('chknluni',8,'cont',201,'zsfc          ',4,101,    &
     &               stat)

      end if

      if(zflat.lt.zsfc) then

        call destroy('chknluni',8,'cont',201,'zflat         ',5,101,    &
     &               stat)

      end if

      if(sthopt.ne.0.and.sthopt.ne.1.and.sthopt.ne.2) then

        call destroy('chknluni',8,'cont',201,'sthopt        ',6,101,    &
     &               stat)

      end if

      if(sthopt.eq.0) then

        if(pname(1:ncpn).ne.'check') then
          dzmin=dz
        end if

      else

        if(dzmin.lt.eps.or.dzmin.gt.dz) then

          call destroy('chknluni',8,'cont',201,'dzmin         ',5,101,  &
     &                 stat)

        end if

        if(layer1.lt.zsfc) then

          call destroy('chknluni',8,'cont',201,'layer1        ',6,101,  &
     &                 stat)

        end if

        if(layer2.lt.layer1) then

          call destroy('chknluni',8,'cont',201,'layer2        ',6,101,  &
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

        call destroy('chknluni',8,'cont',201,'stime         ',5,101,    &
     &               stat)

      end if

      if(etime.lt.0.e0) then

        call destroy('chknluni',8,'cont',201,'etime         ',5,101,    &
     &               stat)

      end if

      if(etime.lt.stime) then

        call destroy('chknluni',8,'cont',201,'etime         ',5,101,    &
     &               stat)

      end if

      if(1000_i8*int(stime+.1e0,i8).gt.1000_i8*int(etime+.1e0,i8)) then

        call destroy('chknluni',8,'cont',201,'stime         ',5,101,    &
     &               stat)

      end if

! -----

! For the sections, mapfcter, coriolis, sfcphys, turbulen, cloudphy,
! asolproc and mixtrace.

      if(mfcopt.ne.0.and.mfcopt.ne.1) then

        call destroy('chknluni',8,'cont',201,'mfcopt        ',6,101,    &
     &               stat)

      end if

      if(coropt.ne.0.and.coropt.ne.1.and.coropt.ne.2) then

        call destroy('chknluni',8,'cont',201,'coropt        ',6,101,    &
     &               stat)

      end if

      if(sfcopt.ne.0.and.sfcopt.ne.1.and.sfcopt.ne.2.and.sfcopt.ne.3    &
     &  .and.sfcopt.ne.11.and.sfcopt.ne.12.and.sfcopt.ne.13) then

        call destroy('chknluni',8,'cont',201,'sfcopt        ',6,101,    &
     &               stat)

      end if

      if(cphopt.ne.-4.and.cphopt.ne.-3.and.cphopt.ne.-2                 &
     &  .and.cphopt.ne.0.and.cphopt.ne.1.and.cphopt.ne.2                &
     &  .and.cphopt.ne.3.and.cphopt.ne.4.and.cphopt.ne.11) then

        call destroy('chknluni',8,'cont',201,'cphopt        ',6,101,    &
     &               stat)

      end if

      if(abs(cphopt).lt.10) then

        if(abs(cphopt).ge.2) then

          if(haiopt.ne.0.and.haiopt.ne.1) then

            call destroy('chknluni',8,'cont',201,'haiopt        ',6,101,&
     &                   stat)

          end if

        end if

      end if

      if(cphopt.lt.0) then

        if(qcgopt.ne.1.and.qcgopt.ne.2) then

          call destroy('chknluni',8,'cont',201,'qcgopt        ',6,101,  &
     &                 stat)

        end if

      end if

      if(aslopt.ne.0.and.aslopt.ne.1) then

        call destroy('chknluni',8,'cont',201,'aslopt        ',6,101,    &
     &               stat)

      end if

      if(trkopt.ne.0.and.trkopt.ne.1.and.trkopt.ne.2) then

        call destroy('chknluni',8,'cont',201,'trkopt        ',6,101,    &
     &               stat)

      end if

      if(tubopt.ne.0.and.tubopt.ne.1                                    &
     &  .and.tubopt.ne.2.and.tubopt.ne.3) then

        call destroy('chknluni',8,'cont',201,'tubopt        ',6,101,    &
     &               stat)

      end if

! -----

! For the section, outfomat.

      if(dmpfmt.ne.2) then

        call destroy('chknluni',8,'cont',201,'dmpfmt        ',6,101,    &
     &               stat)

      end if

      if(dmplev.ne.1.and.dmplev.ne.2.and.dmplev.ne.3                    &
     &  .and.dmplev.ne.11.and.dmplev.ne.12.and.dmplev.ne.13) then

        call destroy('chknluni',8,'cont',201,'dmplev        ',6,101,    &
     &               stat)

      end if

      if(mod(dmplev,10).eq.2.or.mod(dmplev,10).eq.3) then

        if(zdim.lt.10) then

          call destroy('chknluni',8,'cont',201,'zdim          ',4,101,  &
     &                 stat)

        end if

      end if

      if(dmpmon.ne.0.and.dmpmon.ne.1) then

        call destroy('chknluni',8,'cont',201,'dmpmon        ',6,101,    &
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

        call destroy('chknluni',8,'cont',201,'dmpvar        ',6,101,    &
     &               stat)

      end if

! -----

! For the section, uniconf_uni.

      if(fltyp_uni(1:3).ne.'all'.and.fltyp_uni(1:3).ne.'dmp'            &
     &  .and.fltyp_uni(1:3).ne.'mon'.and.fltyp_uni(1:3).ne.'geo') then

        call destroy('chknluni',8,'cont',201,'fltyp_uni     ',9,101,    &
     &               stat)

      end if

      if(dmpmon.eq.0) then

        if(fltyp_uni(1:3).eq.'mon') then

          call destroy('chknluni',8,'cont',201,'fltyp_uni     ',9,101,  &
     &                 stat)

        end if

      end if

      if(rmopt_uni.ne.0.and.rmopt_uni.ne.1) then

        call destroy('chknluni',8,'cont',201,'rmopt_uni     ',9,101,    &
     &               stat)

      end if

      ierr=0

      if(uniopt_uni.ge.0                                                &
     &  .and.uniopt_uni.ne.1.and.uniopt_uni.ne.2                        &
     &  .and.uniopt_uni.ne.3.and.uniopt_uni.ne.4                        &
     &  .and.uniopt_uni.ne.5.and.uniopt_uni.ne.6) then

        ierr=ierr+1

      end if

      if(uniopt_uni.ge.-10.and.uniopt_uni.lt.0                          &
     &  .and.uniopt_uni.ne.-1.and.uniopt_uni.ne.-2                      &
     &  .and.uniopt_uni.ne.-3.and.uniopt_uni.ne.-4                      &
     &  .and.uniopt_uni.ne.-5.and.uniopt_uni.ne.-6) then

        ierr=ierr+1

      end if

      if(uniopt_uni.ge.-20.and.uniopt_uni.lt.-10                        &
     &  .and.uniopt_uni.ne.-11.and.uniopt_uni.ne.-12                    &
     &  .and.uniopt_uni.ne.-13.and.uniopt_uni.ne.-14                    &
     &  .and.uniopt_uni.ne.-15.and.uniopt_uni.ne.-16) then

        ierr=ierr+1

      end if

      if(uniopt_uni.lt.-20                                              &
     &  .and.uniopt_uni.ne.-21.and.uniopt_uni.ne.-22                    &
     &  .and.uniopt_uni.ne.-23.and.uniopt_uni.ne.-24                    &
     &  .and.uniopt_uni.ne.-25.and.uniopt_uni.ne.-26) then

        ierr=ierr+1

      end if

      if(ierr.ne.0) then

        call destroy('chknluni',8,'cont',201,'uniopt_uni    ',10,101,   &
     &               stat)

      end if

      if(abs(mod(uniopt_uni,10)).eq.6) then

        if(ugroup_uni.lt.0.or.ugroup_uni.ge.xgroup*ygroup) then

          call destroy('chknluni',8,'cont',201,'ugroup_uni    ',10,101, &
     &                 stat)

        end if

      end if

      if(xgroup*ygroup.eq.1) then

        if(abs(mod(uniopt_uni,10)).le.4) then

          call destroy('chknluni',8,'cont',201,'uniopt_uni    ',10,101, &
     &                 stat)

        end if

      end if

      if(xsub*ysub.eq.1) then

        if(abs(mod(uniopt_uni,10)).ge.5) then

          call destroy('chknluni',8,'cont',201,'uniopt_uni    ',10,101, &
     &                 stat)

        end if

      end if

      if(flitv_uni.lt.0.e0) then

        call destroy('chknluni',8,'cont',201,'flitv_uni     ',9,101,    &
     &               stat)

      end if

      if(pname(1:ncpn).ne.'check') then
        flitv_uni=aint(flitv_uni+.1e0)
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

        call destroy('chknluni',8,'cont',11,'              ',14,101,    &
     &               stat)

      else

        if(mod((xdim-3),xgroup*xsub).ne.0                               &
     &    .or.mod((ydim-3),ygroup*ysub).ne.0) then

          stat=stat+1

          call destroy('chknluni',8,'cont',11,'              ',14,101,  &
     &                 stat)

        end if

      end if

      if(numpe.lt.xsub*ysub.or.numpe.gt.xsub*ysub*xgroup*ygroup) then

        stat=stat+1

        call destroy('chknluni',8,'cont',11,'              ',14,101,    &
     &               stat)

      end if

! -----

      end subroutine s_chknluni

!-----7--------------------------------------------------------------7--

      end module m_chknluni
