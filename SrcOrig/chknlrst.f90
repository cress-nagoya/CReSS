!***********************************************************************
      module m_chknlrst
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/01/20
!     Modification: 2007/03/10, 2007/03/23, 2007/07/30, 2007/09/04,
!                   2007/09/28, 2007/11/26, 2008/01/11, 2008/03/12,
!                   2008/04/17, 2008/05/02, 2008/08/25, 2008/10/10,
!                   2009/01/30, 2009/02/27, 2011/01/19, 2011/07/15,
!                   2011/08/09, 2011/08/18, 2011/09/22, 2011/11/10,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     check the namelist variables for rstruct.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
      use m_defname
      use m_destroy
      use m_numchar
      use m_outstd14

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: chknlrst, s_chknlrst

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface chknlrst

        module procedure s_chknlrst

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic aint
      intrinsic int
      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_chknlrst(pname,ncpn,stat)
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

! -----

! For the section, sysdep.

      if(savmem.ne.0.and.savmem.ne.1) then

        call destroy('chknlrst',8,'cont',201,'savmem        ',6,101,    &
     &               stat)

      end if

! -----

! For the section, runame.

      call numchar(exprim,1,ncexp,ncspc)

      if(ncexp.gt.64.or.ncspc.ne.0) then

        call destroy('chknlrst',8,'cont',201,'exprim        ',6,101,    &
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

        call destroy('chknlrst',8,'cont',201,'crsdir        ',6,101,    &
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

        call destroy('chknlrst',8,'cont',201,'datdir        ',6,101,    &
     &               stat)

      end if

      if(pname(1:ncpn).ne.'check') then

        call outstd14('chknlrst',8,crsdir,datdir,nccrs,ncdat)

      end if

! -----

! For the section, dimset.

      if(xdim.lt.4) then

        call destroy('chknlrst',8,'cont',201,'xdim          ',4,101,    &
     &               stat)

      end if

      if(mpopt.eq.5) then

        if(ydim.ne.4) then

          call destroy('chknlrst',8,'cont',201,'ydim          ',4,101,  &
     &                 stat)

        end if

      else

        if(ydim.lt.4) then

          call destroy('chknlrst',8,'cont',201,'ydim          ',4,101,  &
     &                 stat)

        end if

      end if

      if(zdim.lt.7) then

        call destroy('chknlrst',8,'cont',201,'zdim          ',4,101,    &
     &               stat)

      end if

! -----

! For the section, project.

      if(mpopt.ne.0.and.mpopt.ne.1.and.mpopt.ne.2.and.mpopt.ne.3.and.   &
     &   mpopt.ne.4.and.mpopt.ne.5.and.mpopt.ne.10.and.mpopt.ne.13) then

        call destroy('chknlrst',8,'cont',201,'mpopt         ',5,101,    &
     &               stat)

      end if

! -----

! For the section, terrain.

      if(trnopt.ne.0.and.trnopt.ne.1.and.trnopt.ne.2) then

        call destroy('chknlrst',8,'cont',201,'trnopt        ',6,101,    &
     &               stat)

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

        call destroy('chknlrst',8,'cont',201,'stime         ',5,101,    &
     &               stat)

      end if

      if(etime.lt.0.e0) then

        call destroy('chknlrst',8,'cont',201,'etime         ',5,101,    &
     &               stat)

      end if

      if(etime.lt.stime) then

        call destroy('chknlrst',8,'cont',201,'etime         ',5,101,    &
     &               stat)

      end if

      if(1000_i8*int(stime+.1e0,i8).gt.1000_i8*int(etime+.1e0,i8)) then

        call destroy('chknlrst',8,'cont',201,'stime         ',5,101,    &
     &               stat)

      end if

! -----

! For the section, boundry.

      if(wbc.ne.-1.and.wbc.ne.1.and.wbc.ne.2.and.wbc.ne.3.and.          &
     &   wbc.ne.4.and.wbc.ne.5.and.wbc.ne.6.and.wbc.ne.7.and.           &
     &   wbc.ne.14.and.wbc.ne.15.and.wbc.ne.16) then

        call destroy('chknlrst',8,'cont',201,'wbc           ',3,101,    &
     &               stat)

      end if

      if(ebc.ne.-1.and.ebc.ne.1.and.ebc.ne.2.and.ebc.ne.3.and.          &
     &   ebc.ne.4.and.ebc.ne.5.and.ebc.ne.6.and.ebc.ne.7.and.           &
     &   ebc.ne.14.and.ebc.ne.15.and.ebc.ne.16) then

        call destroy('chknlrst',8,'cont',201,'ebc           ',3,101,    &
     &               stat)

      end if

      if(sbc.ne.-1.and.sbc.ne.1.and.sbc.ne.2.and.sbc.ne.3.and.          &
     &   sbc.ne.4.and.sbc.ne.5.and.sbc.ne.6.and.sbc.ne.7.and.           &
     &   sbc.ne.14.and.sbc.ne.15.and.sbc.ne.16) then

        call destroy('chknlrst',8,'cont',201,'sbc           ',3,101,    &
     &               stat)

      end if

      if(nbc.ne.-1.and.nbc.ne.1.and.nbc.ne.2.and.nbc.ne.3.and.          &
     &   nbc.ne.4.and.nbc.ne.5.and.nbc.ne.6.and.nbc.ne.7.and.           &
     &   nbc.ne.14.and.nbc.ne.15.and.nbc.ne.16) then

        call destroy('chknlrst',8,'cont',201,'nbc           ',3,101,    &
     &               stat)

      end if

      if(bbc.ne.2.and.bbc.ne.3) then

        call destroy('chknlrst',8,'cont',201,'bbc           ',3,101,    &
     &               stat)

      end if

      if(tbc.ne.2.and.tbc.ne.3) then

        call destroy('chknlrst',8,'cont',201,'tbc           ',3,101,    &
     &               stat)

      end if

      if(mpopt.eq.1.or.mpopt.eq.2.or.mpopt.eq.3) then

        if(abs(wbc).eq.1.or.abs(ebc).eq.1) then

          call destroy('chknlrst',8,'cont',201,'wbc           ',3,101,  &
     &                 stat)

          call destroy('chknlrst',8,'cont',201,'ebc           ',3,101,  &
     &                 stat)

        end if

      end if

      if(mpopt.eq.1.or.mpopt.eq.2.or.mpopt.eq.3                         &
     &  .or.mpopt.eq.10.or.mpopt.eq.13) then

        if(abs(sbc).eq.1.or.abs(nbc).eq.1) then

          call destroy('chknlrst',8,'cont',201,'sbc           ',3,101,  &
     &                 stat)

          call destroy('chknlrst',8,'cont',201,'nbc           ',3,101,  &
     &                 stat)

        end if

      end if

      if(mpopt.eq.5) then

        if(abs(wbc).eq.1.or.abs(ebc).eq.1) then

          call destroy('chknlrst',8,'cont',201,'wbc           ',3,101,  &
     &                 stat)

          call destroy('chknlrst',8,'cont',201,'ebc           ',3,101,  &
     &                 stat)

        end if

        if(abs(sbc).ne.1.or.abs(nbc).ne.1) then

          call destroy('chknlrst',8,'cont',201,'sbc           ',3,101,  &
     &                 stat)

          call destroy('chknlrst',8,'cont',201,'nbc           ',3,101,  &
     &                 stat)

        end if

      end if

      if(abs(wbc).eq.1.or.abs(ebc).eq.1) then

        if(wbc.ne.ebc) then

          call destroy('chknlrst',8,'cont',201,'wbc           ',3,101,  &
     &                 stat)

          call destroy('chknlrst',8,'cont',201,'ebc           ',3,101,  &
     &                 stat)

        end if

      end if

      if(xdim.eq.4) then

        if(wbc.ne.ebc) then

          call destroy('chknlrst',8,'cont',201,'wbc           ',3,101,  &
     &                 stat)

          call destroy('chknlrst',8,'cont',201,'ebc           ',3,101,  &
     &                 stat)

        end if

      end if

      if(abs(sbc).eq.1.or.abs(nbc).eq.1) then

        if(sbc.ne.nbc) then

          call destroy('chknlrst',8,'cont',201,'sbc           ',3,101,  &
     &                 stat)

          call destroy('chknlrst',8,'cont',201,'nbc           ',3,101,  &
     &                 stat)

        end if

      end if

      if(ydim.eq.4) then

        if(sbc.ne.nbc) then

          call destroy('chknlrst',8,'cont',201,'sbc           ',3,101,  &
     &                 stat)

          call destroy('chknlrst',8,'cont',201,'nbc           ',3,101,  &
     &                 stat)

        end if

      end if

      if(numpe.ne.xsub*ysub*xgroup*ygroup) then

        if(wbc.eq.1.or.ebc.eq.1) then

          call destroy('chknlrst',8,'cont',201,'wbc           ',3,101,  &
     &                 stat)

          call destroy('chknlrst',8,'cont',201,'ebc           ',3,101,  &
     &                 stat)

        end if

        if(sbc.eq.1.or.nbc.eq.1) then

          call destroy('chknlrst',8,'cont',201,'sbc           ',3,101,  &
     &                 stat)

          call destroy('chknlrst',8,'cont',201,'nbc           ',3,101,  &
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

! For the section, sfcphys.

      if(sfcopt.ne.0.and.sfcopt.ne.1.and.sfcopt.ne.2.and.sfcopt.ne.3    &
     &  .and.sfcopt.ne.11.and.sfcopt.ne.12.and.sfcopt.ne.13) then

        call destroy('chknlrst',8,'cont',201,'sfcopt        ',6,101,    &
     &               stat)

      end if

      if(sfcopt.ge.1) then

        if(levund.lt.4.or.levund.gt.zdim) then

          call destroy('chknlrst',8,'cont',201,'levund        ',6,101,  &
     &                 stat)

        end if

      end if

! -----

! For the section, integrat.

      if(advopt.ne.1.and.advopt.ne.2                                    &
     &  .and.advopt.ne.3.and.advopt.ne.4.and.advopt.ne.5) then

        call destroy('chknlrst',8,'cont',201,'advopt        ',6,101,    &
     &               stat)

      end if

! -----

! For the section, diabatic.

      if(diaopt.ne.0.and.diaopt.ne.1) then

        call destroy('chknlrst',8,'cont',201,'diaopt        ',6,101,    &
     &               stat)

      end if

      if(diaopt.eq.1) then

        if(advopt.ge.4) then

          call destroy('chknlrst',8,'cont',201,'advopt        ',6,101,  &
     &                 stat)

        end if

      end if

! -----

! For the section, cloudphy.

      if(cphopt.ne.-4.and.cphopt.ne.-3.and.cphopt.ne.-2                 &
     &  .and.cphopt.ne.0.and.cphopt.ne.1.and.cphopt.ne.2                &
     &  .and.cphopt.ne.3.and.cphopt.ne.4.and.cphopt.ne.11) then

        call destroy('chknlrst',8,'cont',201,'cphopt        ',6,101,    &
     &               stat)

      end if

      if(abs(cphopt).lt.10) then

        if(abs(cphopt).ge.2) then

          if(zdim.lt.10) then

            call destroy('chknlrst',8,'cont',201,'zdim          ',4,101,&
     &                   stat)

          end if

          if(haiopt.ne.0.and.haiopt.ne.1) then

            call destroy('chknlrst',8,'cont',201,'haiopt        ',6,101,&
     &                   stat)

          end if

        end if

      end if

      if(cphopt.lt.0) then

        if(qcgopt.ne.1.and.qcgopt.ne.2) then

          call destroy('chknlrst',8,'cont',201,'qcgopt        ',6,101,  &
     &                 stat)

        end if

      end if

! -----

! For the section, asolproc.

      if(aslopt.ne.0.and.aslopt.ne.1) then

        call destroy('chknlrst',8,'cont',201,'aslopt        ',6,101,    &
     &               stat)

      end if

! -----

! For the section, mixtrace.

      if(trkopt.ne.0.and.trkopt.ne.1.and.trkopt.ne.2) then

        call destroy('chknlrst',8,'cont',201,'trkopt        ',6,101,    &
     &               stat)

      end if

! -----

! For the section, turbulen.

      if(tubopt.ne.0.and.tubopt.ne.1                                    &
     &  .and.tubopt.ne.2.and.tubopt.ne.3) then

        call destroy('chknlrst',8,'cont',201,'tubopt        ',6,101,    &
     &               stat)

      end if

! -----

! For the section, outfomat.

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

        call destroy('chknlrst',8,'cont',201,'dmpvar        ',6,101,    &
     &               stat)

      end if

! -----

! For the section, rstconf_rst.

      if(flitv_rst.lt.0.e0) then

        call destroy('chknlrst',8,'cont',201,'flitv_rst     ',9,101,    &
     &               stat)

      end if

      if(pname(1:ncpn).ne.'check') then
        flitv_rst=aint(flitv_rst+.1e0)
      end if

      if(rmopt_rst.ne.0.and.rmopt_rst.ne.1) then

        call destroy('chknlrst',8,'cont',201,'rmopt_rst     ',9,101,    &
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

        call destroy('chknlrst',8,'cont',11,'              ',14,101,    &
     &               stat)

      else

        if(mod((xdim-3),xgroup*xsub).ne.0                               &
     &    .or.mod((ydim-3),ygroup*ysub).ne.0) then

          stat=stat+1

          call destroy('chknlrst',8,'cont',11,'              ',14,101,  &
     &                 stat)

        end if

      end if

      if(numpe.lt.xsub*ysub.or.numpe.gt.xsub*ysub*xgroup*ygroup) then

        stat=stat+1

        call destroy('chknlrst',8,'cont',11,'              ',14,101,    &
     &               stat)

      end if

      if(xsub_rst.eq.0.or.ysub_rst.eq.0) then

        stat=stat+1

        call destroy('chknlrst',8,'cont',11,'              ',14,101,    &
     &               stat)

      else

        if(mod((xdim-3),xgroup*xsub_rst).ne.0                           &
     &    .or.mod((ydim-3),ygroup*ysub_rst).ne.0) then

          stat=stat+1

          call destroy('chknlrst',8,'cont',11,'              ',14,101,  &
     &                 stat)

        end if

      end if

! -----

      end subroutine s_chknlrst

!-----7--------------------------------------------------------------7--

      end module m_chknlrst
