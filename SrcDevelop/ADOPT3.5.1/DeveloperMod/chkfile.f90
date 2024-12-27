!***********************************************************************
      module m_chkfile
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/01/25, 1999/03/25, 1999/04/06, 1999/05/10,
!                   1999/05/20, 1999/06/07, 1999/07/05, 1999/08/03,
!                   1999/08/18, 1999/09/16, 1999/09/30, 1999/10/12,
!                   1999/11/01, 1999/11/24, 1999/12/17, 2000/01/17,
!                   2000/04/18, 2000/06/01, 2000/12/18, 2001/01/15,
!                   2001/03/13, 2001/05/29, 2001/07/13, 2001/08/07,
!                   2001/10/18, 2001/11/20, 2002/02/05, 2002/06/18,
!                   2002/07/03, 2002/08/15, 2002/09/09, 2002/10/31,
!                   2002/11/11, 2002/12/02, 2003/02/13, 2003/03/13,
!                   2003/03/28, 2003/04/30, 2003/05/19, 2003/07/15,
!                   2003/08/08, 2003/10/10, 2003/11/05, 2004/01/09,
!                   2004/03/05, 2004/04/15, 2004/05/07, 2004/05/31,
!                   2004/08/01, 2004/08/20, 2004/09/01, 2004/09/10,
!                   2005/01/14, 2005/02/10, 2005/04/04, 2005/08/05,
!                   2005/12/13, 2006/01/10, 2006/04/03, 2006/09/21,
!                   2006/09/30, 2006/11/06, 2006/11/27, 2006/12/04,
!                   2007/01/05, 2007/01/20, 2007/03/23, 2007/04/11,
!                   2007/03/23, 2007/04/24, 2007/05/21, 2007/06/27,
!                   2007/11/26, 2008/01/11, 2008/03/12, 2008/04/17,
!                   2008/05/02, 2008/06/09, 2008/07/01, 2008/08/19,
!                   2008/08/25, 2008/12/11, 2009/01/05, 2009/01/30,
!                   2009/02/27, 2010/05/17, 2011/08/18, 2011/09/22,
!                   2011/11/10, 2013/03/27
!     Modification: 2024/12/25 (Satoki Tsujino)

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     check the namelist variables which are read out from the input
!     files.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_commath
      use m_destroy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: chkfile, s_chkfile

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface chkfile

        module procedure s_chkfile

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic mod
      intrinsic sign

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_chkfile(fname,stat,ncn,nin,nrn,                      &
     &                     cname,iname,rname,rcname,riname,rrname)
!***********************************************************************

! Input variables

      character(len=3), intent(in) :: fname
                       ! Character variable of checked file name

      integer, intent(in) :: ncn
                       ! Dimension of character namelist table

      integer, intent(in) :: nin
                       ! Dimension of integer namelist table

      integer, intent(in) :: nrn
                       ! Dimension of real namelist table

      character(len=108), intent(in) :: cname(1:ncn)
                       ! Character namelist table

      character(len=108), intent(in) :: rcname(1:ncn)
                       ! Read out character namelist table

      integer, intent(in) :: iname(-1:nin)
                       ! Integer namelist table

      integer, intent(in) :: riname(1:nin)
                       ! Read out integer namelist table

      real, intent(in) :: rname(1:nrn)
                       ! Real namelist table

      real, intent(in) :: rrname(1:nrn)
                       ! Read out real namelist table

! Output variable

      integer, intent(out) :: stat
                       ! Runtime status

! Internal shared variables

      integer ierr_lnd ! Error descriptor for lnduse_lnd

      integer ierr_albe
                       ! Error descriptor for albe_lnd

      integer ierr_beta
                       ! Error descriptor for beta_lnd

      integer ierr_z0m ! Error descriptor for z0m_lnd
      integer ierr_z0h ! Error descriptor for z0h_lnd
      integer ierr_cap ! Error descriptor for cap_lnd
      integer ierr_nuu ! Error descriptor for nuu_lnd

      real chkeps      ! Very small constant
                       ! for real variables comparison

      real crn         ! Temporary variable
      real crrn        ! Temporary variable

! Internal private variables

      integer iid      ! Index of do loops

      real crn_sub     ! Substitute for crn
      real crrn_sub    ! Substitute for crrn

!-----7--------------------------------------------------------------7--

! Do nothing.

      if(fname(1:3).eq.'ctl') then

        return

      end if

! -----

! Set the common used variables.

      stat=0

      chkeps=1.e-5

! -----

! Check the read out namelist variables from the all input files.

      if(fname(1:3).ne.'und') then

        if(rcname(idexprim)(1:108).ne.cname(idexprim)(1:108)) then

          call destroy('chkfile ',7,'cont',202,'exprim        ',6,101,  &
     &                 stat)

        end if

      end if

      if(riname(idxdim).ne.iname(idxdim)) then

        call destroy('chkfile ',7,'cont',202,'xdim          ',4,101,    &
     &               stat)

      end if

      if(riname(idydim).ne.iname(idydim)) then

        call destroy('chkfile ',7,'cont',202,'ydim          ',4,101,    &
     &               stat)

      end if

      if(fname(1:3).ne.'und'                                            &
     &  .and.fname(1:3).ne.'geo'.and.fname(1:3).ne.'mon') then

        if(riname(idzdim).ne.iname(idzdim)) then

          call destroy('chkfile ',7,'cont',202,'zdim          ',4,101,  &
     &                 stat)

        end if

      end if

      if(riname(idnumpe).ne.iname(idnumpe)) then

        call destroy('chkfile ',7,'cont',202,'numpe         ',5,101,    &
     &               stat)

      end if

      if(riname(idxgroup).ne.iname(idxgroup)) then

        call destroy('chkfile ',7,'cont',202,'xgroup        ',6,101,    &
     &               stat)

      end if

      if(riname(idygroup).ne.iname(idygroup)) then

        call destroy('chkfile ',7,'cont',202,'ygroup        ',6,101,    &
     &               stat)

      end if

      if(riname(idxsub).ne.iname(idxsub)) then

        call destroy('chkfile ',7,'cont',202,'xsub          ',4,101,    &
     &               stat)

      end if

      if(riname(idysub).ne.iname(idysub)) then

        call destroy('chkfile ',7,'cont',202,'ysub          ',4,101,    &
     &               stat)

      end if

      if(riname(idmpopt).ne.iname(idmpopt)) then

        call destroy('chkfile ',7,'cont',202,'mpopt         ',5,101,    &
     &               stat)

      end if

      if(fname(1:3).ne.'rst') then

        if(riname(idmpopt).eq.1.or.riname(idmpopt).eq.2                 &
     &    .or.riname(idmpopt).eq.3.or.riname(idmpopt).eq.13) then

          if(riname(idnspol).ne.iname(idnspol)) then

            call destroy('chkfile ',7,'cont',202,'nspol         ',5,101,&
     &                   stat)

          end if

          crn=rname(idtlat1)+sign(eps,rname(idtlat1))
          crrn=rrname(idtlat1)+sign(eps,rrname(idtlat1))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'tlat1         ',5,101,&
     &                   stat)

          end if

        end if

        if(riname(idmpopt).eq.2) then

          crn=rname(idtlat2)+sign(eps,rname(idtlat2))
          crrn=rrname(idtlat2)+sign(eps,rrname(idtlat2))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'tlat2         ',5,101,&
     &                   stat)

          end if

        end if

        if(riname(idmpopt).eq.1                                         &
     &    .or.riname(idmpopt).eq.2.or.riname(idmpopt).eq.4) then

          crn=rname(idtlon)+sign(eps,rname(idtlon))
          crrn=rrname(idtlon)+sign(eps,rrname(idtlon))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'tlon          ',4,101,&
     &                   stat)

          end if

        end if

        crn=rname(idulat)+sign(eps,rname(idulat))
        crrn=rrname(idulat)+sign(eps,rrname(idulat))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'ulat          ',4,101,  &
     &                 stat)

        end if

        crn=rname(idulon)+sign(eps,rname(idulon))
        crrn=rrname(idulon)+sign(eps,rrname(idulon))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'ulon          ',4,101,  &
     &                 stat)

        end if

        crn=rname(idriu)+sign(eps,rname(idriu))
        crrn=rrname(idriu)+sign(eps,rrname(idriu))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'riu           ',3,101,  &
     &                 stat)

        end if

        crn=rname(idrju)+sign(eps,rname(idrju))
        crrn=rrname(idrju)+sign(eps,rrname(idrju))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'rju           ',3,101,  &
     &                 stat)

        end if

        crn=rname(iddx)+sign(eps,rname(iddx))
        crrn=rrname(iddx)+sign(eps,rrname(iddx))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'dx            ',2,101,  &
     &                 stat)

        end if

        crn=rname(iddy)+sign(eps,rname(iddy))
        crrn=rrname(iddy)+sign(eps,rrname(iddy))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'dy            ',2,101,  &
     &                 stat)

        end if

        if(fname(1:3).ne.'und') then

          crn=rname(iddz)+sign(eps,rname(iddz))
          crrn=rrname(iddz)+sign(eps,rrname(iddz))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'dz            ',2,101,&
     &                   stat)

          end if

        end if

      end if

! -----

! Check the read out namelist variables from the restart, the GPV data,
! the aerosol data, the radar data files or the dumped file.

      if(fname(1:3).eq.'res'.or.fname(1:3).eq.'gpv'.or.                 &
     &   fname(1:3).eq.'asl'.or.fname(1:3).eq.'rdr'.or.                 &
     &   fname(1:3).eq.'geo'.or.fname(1:3).eq.'dmp'.or.                 &
     &   fname(1:3).eq.'mon') then

        crn=rname(idzflat)+sign(eps,rname(idzflat))
        crrn=rrname(idzflat)+sign(eps,rrname(idzflat))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'zflat         ',5,101,  &
     &                 stat)

        end if

        crn=rname(idzsfc)+sign(eps,rname(idzsfc))
        crrn=rrname(idzsfc)+sign(eps,rrname(idzsfc))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'zsfc          ',4,101,  &
     &                 stat)

        end if

        if(riname(idsthopt).ne.iname(idsthopt)) then

          call destroy('chkfile ',7,'cont',202,'sthopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idsthopt).eq.1.or.riname(idsthopt).eq.2) then

          crn=rname(iddzmin)+sign(eps,rname(iddzmin))
          crrn=rrname(iddzmin)+sign(eps,rrname(iddzmin))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'dzmin         ',5,101,&
     &                   stat)

          end if

          crn=rname(idlayer1)+sign(eps,rname(idlayer1))
          crrn=rrname(idlayer1)+sign(eps,rrname(idlayer1))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'layer1        ',6,101,&
     &                   stat)

          end if

          crn=rname(idlayer2)+sign(eps,rname(idlayer2))
          crrn=rrname(idlayer2)+sign(eps,rrname(idlayer2))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'layer2        ',6,101,&
     &                   stat)

          end if

        end if

      end if

! -----

! Check the read out namelist variables from the restart, the GPV data,
! the aerosol data or the radar data files.

      if(fname(1:3).eq.'res'.or.fname(1:3).eq.'gpv'                     &
     &  .or.fname(1:3).eq.'asl'.or.fname(1:3).eq.'rdr') then

        if(riname(idtrnopt).ne.iname(idtrnopt)) then

          call destroy('chkfile ',7,'cont',202,'trnopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idtrnopt).eq.0.or.riname(idtrnopt).eq.1) then

          crn=rname(idmnthgh)+sign(eps,rname(idmnthgh))
          crrn=rrname(idmnthgh)+sign(eps,rrname(idmnthgh))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'mnthgh(1)     ',9,101,&
     &                   stat)

          end if

          crn=rname(idmnthgh+1)+sign(eps,rname(idmnthgh+1))
          crrn=rrname(idmnthgh+1)+sign(eps,rrname(idmnthgh+1))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'mnthgh(2)     ',9,101,&
     &                   stat)

          end if

        end if

        if(riname(idtrnopt).eq.1) then

          crn=rname(idmntwx)+sign(eps,rname(idmntwx))
          crrn=rrname(idmntwx)+sign(eps,rrname(idmntwx))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'mntwx         ',5,101,&
     &                   stat)

          end if

          crn=rname(idmntwy)+sign(eps,rname(idmntwy))
          crrn=rrname(idmntwy)+sign(eps,rrname(idmntwy))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'mntwy         ',5,101,&
     &                   stat)

          end if

          crn=rname(idmntcx)+sign(eps,rname(idmntcx))
          crrn=rrname(idmntcx)+sign(eps,rrname(idmntcx))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'mntcx         ',5,101,&
     &                   stat)

          end if

          crn=rname(idmntcy)+sign(eps,rname(idmntcy))
          crrn=rrname(idmntcy)+sign(eps,rrname(idmntcy))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'mntcy         ',5,101,&
     &                   stat)

          end if

        end if

        if(rcname(ididate).ne.cname(ididate)) then

          call destroy('chkfile ',7,'cont',202,'idate         ',5,101,  &
     &                 stat)

        end if

      end if

! -----

! Check the read out namelist variables from the restart or the GPV data
! files.

      if(fname(1:3).eq.'res'.or.fname(1:3).eq.'gpv') then

        if(rcname(idgpvvar)(1:8).ne.cname(idgpvvar)(1:8)) then

          call destroy('chkfile ',7,'cont',202,'gpvvar        ',6,101,  &
     &                 stat)

        end if

        if(riname(idnggopt).eq.1.or.riname(idexbopt).ge.1.or.           &
     &     mod(riname(idlspopt),10).eq.1.or.riname(idvspopt).eq.1) then

          crn=rname(idgpvitv)+sign(eps,rname(idgpvitv))
          crrn=rrname(idgpvitv)+sign(eps,rrname(idgpvitv))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'gpvitv        ',6,101,&
     &                   stat)

          end if

        end if

        if(riname(idmpopt_gpv).ne.iname(idmpopt_gpv)) then

          call destroy('chkfile ',7,'cont',202,'mpopt_gpv     ',9,101,  &
     &                 stat)

        end if

        if(riname(idmpopt_gpv).eq.1.or.riname(idmpopt_gpv).eq.2.or.     &
     &    riname(idmpopt_gpv).eq.3.or.riname(idmpopt_gpv).eq.13) then

          if(riname(idnspol_gpv).ne.iname(idnspol_gpv)) then

            call destroy('chkfile ',7,'cont',202,'nspol_gpv     ',9,101,&
     &                   stat)

          end if

          crn=rname(idtlat1_gpv)+sign(eps,rname(idtlat1_gpv))
          crrn=rrname(idtlat1_gpv)+sign(eps,rrname(idtlat1_gpv))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'tlat1_gpv     ',9,101,&
     &                   stat)

          end if

        end if

        if(riname(idmpopt_gpv).eq.2) then

          crn=rname(idtlat2_gpv)+sign(eps,rname(idtlat2_gpv))
          crrn=rrname(idtlat2_gpv)+sign(eps,rrname(idtlat2_gpv))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'tlat2_gpv     ',9,101,&
     &                   stat)

          end if

        end if

        if(riname(idmpopt_gpv).eq.1.or.                                 &
     &     riname(idmpopt_gpv).eq.2.or.riname(idmpopt_gpv).eq.4) then

          crn=rname(idtlon_gpv)+sign(eps,rname(idtlon_gpv))
          crrn=rrname(idtlon_gpv)+sign(eps,rrname(idtlon_gpv))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'tlon_gpv      ',8,101,&
     &                   stat)

          end if

        end if

        if(riname(idxdim_gpv).ne.iname(idxdim_gpv)) then

          call destroy('chkfile ',7,'cont',202,'xdim_gpv      ',8,101,  &
     &                 stat)

        end if

        if(riname(idydim_gpv).ne.iname(idydim_gpv)) then

          call destroy('chkfile ',7,'cont',202,'ydim_gpv      ',8,101,  &
     &                 stat)

        end if

        if(riname(idzdim_gpv).ne.iname(idzdim_gpv)) then

          call destroy('chkfile ',7,'cont',202,'zdim_gpv      ',8,101,  &
     &                 stat)

        end if

        crn=rname(iddx_gpv)+sign(eps,rname(iddx_gpv))
        crrn=rrname(iddx_gpv)+sign(eps,rrname(iddx_gpv))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'dx_gpv        ',6,101,  &
     &                 stat)

        end if

        crn=rname(iddy_gpv)+sign(eps,rname(iddy_gpv))
        crrn=rrname(iddy_gpv)+sign(eps,rrname(iddy_gpv))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'dy_gpv        ',6,101,  &
     &                 stat)

        end if

        crn=rname(idulat_gpv)+sign(eps,rname(idulat_gpv))
        crrn=rrname(idulat_gpv)+sign(eps,rrname(idulat_gpv))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'ulat_gpv      ',8,101,  &
     &                 stat)

        end if

        crn=rname(idulon_gpv)+sign(eps,rname(idulon_gpv))
        crrn=rrname(idulon_gpv)+sign(eps,rrname(idulon_gpv))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'ulon_gpv      ',8,101,  &
     &                 stat)

        end if

        crn=rname(idriu_gpv)+sign(eps,rname(idriu_gpv))
        crrn=rrname(idriu_gpv)+sign(eps,rrname(idriu_gpv))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'riu_gpv       ',7,101,  &
     &                 stat)

        end if

        crn=rname(idrju_gpv)+sign(eps,rname(idrju_gpv))
        crrn=rrname(idrju_gpv)+sign(eps,rrname(idrju_gpv))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'rju_gpv       ',7,101,  &
     &                 stat)

        end if

        if(riname(idintopt_gpv).ne.iname(idintopt_gpv)) then

          call destroy('chkfile ',7,'cont',202,'intopt_gpv    ',10,101, &
     &                 stat)

        end if

        if(riname(idrotopt_gpv).ne.iname(idrotopt_gpv)) then

          call destroy('chkfile ',7,'cont',202,'rotopt_gpv    ',10,101, &
     &                 stat)

        end if

        if(rcname(iddatype_gpv)(1:2).ne.cname(iddatype_gpv)(1:2)) then

          call destroy('chkfile ',7,'cont',202,'datype_gpv    ',10,101, &
     &                 stat)

        end if

        if(riname(idrefsfc_gpv).ne.iname(idrefsfc_gpv)) then

          call destroy('chkfile ',7,'cont',202,'refsfc_gpv    ',10,101, &
     &                 stat)

        end if

        if(rcname(idetrvar_gpv)(1:7).ne.cname(idetrvar_gpv)(1:7)) then

          call destroy('chkfile ',7,'cont',202,'etrvar_gpv    ',10,101, &
     &                 stat)

        end if

      end if

! -----

! Check the read out namelist variables from the restart or the aerosol
! data files.

      if(fname(1:3).eq.'res'.or.fname(1:3).eq.'asl') then

        if(riname(idnggopt).eq.1.or.riname(idexbopt).ge.1.or.           &
     &     mod(riname(idlspopt),10).eq.1.or.riname(idvspopt).eq.1) then

          crn=rname(idaslitv)+sign(eps,rname(idaslitv))
          crrn=rrname(idaslitv)+sign(eps,rrname(idaslitv))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'aslitv        ',6,101,&
     &                   stat)

          end if

        end if

        if(riname(idmpopt_asl).ne.iname(idmpopt_asl)) then

          call destroy('chkfile ',7,'cont',202,'mpopt_asl     ',9,101,  &
     &                 stat)

        end if

        if(riname(idmpopt_asl).eq.1.or.riname(idmpopt_asl).eq.2.or.     &
     &    riname(idmpopt_asl).eq.3.or.riname(idmpopt_asl).eq.13) then

          if(riname(idnspol_asl).ne.iname(idnspol_asl)) then

            call destroy('chkfile ',7,'cont',202,'nspol_asl     ',9,101,&
     &                   stat)

          end if

          crn=rname(idtlat1_asl)+sign(eps,rname(idtlat1_asl))
          crrn=rrname(idtlat1_asl)+sign(eps,rrname(idtlat1_asl))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'tlat1_asl     ',9,101,&
     &                   stat)

          end if

        end if

        if(riname(idmpopt_asl).eq.2) then

          crn=rname(idtlat2_asl)+sign(eps,rname(idtlat2_asl))
          crrn=rrname(idtlat2_asl)+sign(eps,rrname(idtlat2_asl))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'tlat2_asl     ',9,101,&
     &                   stat)

          end if

        end if

        if(riname(idmpopt_asl).eq.1.or.                                 &
     &     riname(idmpopt_asl).eq.2.or.riname(idmpopt_asl).eq.4) then

          crn=rname(idtlon_asl)+sign(eps,rname(idtlon_asl))
          crrn=rrname(idtlon_asl)+sign(eps,rrname(idtlon_asl))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'tlon_asl      ',8,101,&
     &                   stat)

          end if

        end if

        if(riname(idxdim_asl).ne.iname(idxdim_asl)) then

          call destroy('chkfile ',7,'cont',202,'xdim_asl      ',8,101,  &
     &                 stat)

        end if

        if(riname(idydim_asl).ne.iname(idydim_asl)) then

          call destroy('chkfile ',7,'cont',202,'ydim_asl      ',8,101,  &
     &                 stat)

        end if

        if(riname(idzdim_asl).ne.iname(idzdim_asl)) then

          call destroy('chkfile ',7,'cont',202,'zdim_asl      ',8,101,  &
     &                 stat)

        end if

        crn=rname(iddx_asl)+sign(eps,rname(iddx_asl))
        crrn=rrname(iddx_asl)+sign(eps,rrname(iddx_asl))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'dx_asl        ',6,101,  &
     &                 stat)

        end if

        crn=rname(iddy_asl)+sign(eps,rname(iddy_asl))
        crrn=rrname(iddy_asl)+sign(eps,rrname(iddy_asl))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'dy_asl        ',6,101,  &
     &                 stat)

        end if

        crn=rname(idulat_asl)+sign(eps,rname(idulat_asl))
        crrn=rrname(idulat_asl)+sign(eps,rrname(idulat_asl))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'ulat_asl      ',8,101,  &
     &                 stat)

        end if

        crn=rname(idulon_asl)+sign(eps,rname(idulon_asl))
        crrn=rrname(idulon_asl)+sign(eps,rrname(idulon_asl))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'ulon_asl      ',8,101,  &
     &                 stat)

        end if

        crn=rname(idriu_asl)+sign(eps,rname(idriu_asl))
        crrn=rrname(idriu_asl)+sign(eps,rrname(idriu_asl))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'riu_asl       ',7,101,  &
     &                 stat)

        end if

        crn=rname(idrju_asl)+sign(eps,rname(idrju_asl))
        crrn=rrname(idrju_asl)+sign(eps,rrname(idrju_asl))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'rju_asl       ',7,101,  &
     &                 stat)

        end if

        if(riname(idintopt_asl).ne.iname(idintopt_asl)) then

          call destroy('chkfile ',7,'cont',202,'intopt_asl    ',10,101, &
     &                 stat)

        end if

      end if

! -----

! Check the read out namelist variables from the restart or the radar
! data files.

      if(fname(1:3).eq.'res'.or.fname(1:3).eq.'rdr') then

        if(rcname(idrdrvar)(1:4).ne.cname(idrdrvar)(1:4)) then

          call destroy('chkfile ',7,'cont',202,'rdrvar        ',6,101,  &
     &                 stat)

        end if

        if(riname(idngropt).ge.1) then

          crn=rname(idrdritv)+sign(eps,rname(idrdritv))
          crrn=rrname(idrdritv)+sign(eps,rrname(idrdritv))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'rdritv        ',6,101,&
     &                   stat)

          end if

        end if

        if(riname(idmpopt_rdr).ne.iname(idmpopt_rdr)) then

          call destroy('chkfile ',7,'cont',202,'mpopt_rdr     ',9,101,  &
     &                 stat)

        end if

        if(riname(idmpopt_rdr).eq.1.or.riname(idmpopt_rdr).eq.2.or.     &
     &    riname(idmpopt_rdr).eq.3.or.riname(idmpopt_rdr).eq.13) then

          if(riname(idnspol_rdr).ne.iname(idnspol_rdr)) then

            call destroy('chkfile ',7,'cont',202,'nspol_rdr     ',9,101,&
     &                   stat)

          end if

          crn=rname(idtlat1_rdr)+sign(eps,rname(idtlat1_rdr))
          crrn=rrname(idtlat1_rdr)+sign(eps,rrname(idtlat1_rdr))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'tlat1_rdr     ',9,101,&
     &                   stat)

          end if

        end if

        if(riname(idmpopt_rdr).eq.2) then

          crn=rname(idtlat2_rdr)+sign(eps,rname(idtlat2_rdr))
          crrn=rrname(idtlat2_rdr)+sign(eps,rrname(idtlat2_rdr))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'tlat2_rdr     ',9,101,&
     &                   stat)

          end if

        end if

        if(riname(idmpopt_rdr).eq.1.or.                                 &
     &     riname(idmpopt_rdr).eq.2.or.riname(idmpopt_rdr).eq.4) then

          crn=rname(idtlon_rdr)+sign(eps,rname(idtlon_rdr))
          crrn=rrname(idtlon_rdr)+sign(eps,rrname(idtlon_rdr))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'tlon_rdr      ',8,101,&
     &                   stat)

          end if

        end if

        if(riname(idxdim_rdr).ne.iname(idxdim_rdr)) then

          call destroy('chkfile ',7,'cont',202,'xdim_rdr      ',8,101,  &
     &                 stat)

        end if

        if(riname(idydim_rdr).ne.iname(idydim_rdr)) then

          call destroy('chkfile ',7,'cont',202,'ydim_rdr      ',8,101,  &
     &                 stat)

        end if

        if(riname(idzdim_rdr).ne.iname(idzdim_rdr)) then

          call destroy('chkfile ',7,'cont',202,'zdim_rdr      ',8,101,  &
     &                 stat)

        end if

        crn=rname(iddx_rdr)+sign(eps,rname(iddx_rdr))
        crrn=rrname(iddx_rdr)+sign(eps,rrname(iddx_rdr))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'dx_rdr        ',6,101,  &
     &                 stat)

        end if

        crn=rname(iddy_rdr)+sign(eps,rname(iddy_rdr))
        crrn=rrname(iddy_rdr)+sign(eps,rrname(iddy_rdr))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'dy_rdr        ',6,101,  &
     &                 stat)

        end if

        crn=rname(idulat_rdr)+sign(eps,rname(idulat_rdr))
        crrn=rrname(idulat_rdr)+sign(eps,rrname(idulat_rdr))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'ulat_rdr      ',8,101,  &
     &                 stat)

        end if

        crn=rname(idulon_rdr)+sign(eps,rname(idulon_rdr))
        crrn=rrname(idulon_rdr)+sign(eps,rrname(idulon_rdr))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'ulon_rdr      ',8,101,  &
     &                 stat)

        end if

        crn=rname(idriu_rdr)+sign(eps,rname(idriu_rdr))
        crrn=rrname(idriu_rdr)+sign(eps,rrname(idriu_rdr))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'riu_rdr       ',7,101,  &
     &                 stat)

        end if

        crn=rname(idrju_rdr)+sign(eps,rname(idrju_rdr))
        crrn=rrname(idrju_rdr)+sign(eps,rrname(idrju_rdr))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'rju_rdr       ',7,101,  &
     &                 stat)

        end if

        if(riname(idrotopt_rdr).ne.iname(idrotopt_rdr)) then

          call destroy('chkfile ',7,'cont',202,'rotopt_rdr    ',10,101, &
     &                 stat)

        end if

        if(rcname(iddatype_rdr)(1:1).ne.cname(iddatype_rdr)(1:1)) then

          call destroy('chkfile ',7,'cont',202,'datype_rdr    ',10,101, &
     &                 stat)

        end if

        if(rcname(iddatype_rdr)(1:1).eq.'r') then

          crn=rname(idrdrcoe_rdr)+sign(eps,rname(idrdrcoe_rdr))
          crrn=rrname(idrdrcoe_rdr)+sign(eps,rrname(idrdrcoe_rdr))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

           call destroy('chkfile ',7,'cont',202,'rdrcoe_rdr    ',10,101,&
     &                  stat)

          end if

          crn=rname(idrdrexp_rdr)+sign(eps,rname(idrdrexp_rdr))
          crrn=rrname(idrdrexp_rdr)+sign(eps,rrname(idrdrexp_rdr))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

           call destroy('chkfile ',7,'cont',202,'rdrexp_rdr    ',10,101,&
     &                  stat)

          end if

        end if

      end if

! -----

! Check the read out namelist variables from the restart or the terrain
! data files.

      if(fname(1:3).eq.'res'.or.fname(1:3).eq.'trn') then

        if(riname(idmpopt_trn).ne.iname(idmpopt_trn)) then

          call destroy('chkfile ',7,'cont',202,'mpopt_trn     ',9,101,  &
     &                 stat)

        end if

        if(riname(idmpopt_trn).eq.1.or.riname(idmpopt_trn).eq.2.or.     &
     &    riname(idmpopt_trn).eq.3.or.riname(idmpopt_trn).eq.13) then

          if(riname(idnspol_trn).ne.iname(idnspol_trn)) then

            call destroy('chkfile ',7,'cont',202,'nspol_trn     ',9,101,&
     &                   stat)

          end if

          crn=rname(idtlat1_trn)+sign(eps,rname(idtlat1_trn))
          crrn=rrname(idtlat1_trn)+sign(eps,rrname(idtlat1_trn))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'tlat1_trn     ',9,101,&
     &                   stat)

          end if

        end if

        if(riname(idmpopt_trn).eq.2) then

          crn=rname(idtlat2_trn)+sign(eps,rname(idtlat2_trn))
          crrn=rrname(idtlat2_trn)+sign(eps,rrname(idtlat2_trn))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'tlat2_trn     ',9,101,&
     &                   stat)

          end if

        end if

        if(riname(idmpopt_trn).eq.1.or.                                 &
     &     riname(idmpopt_trn).eq.2.or.riname(idmpopt_trn).eq.4) then

          crn=rname(idtlon_trn)+sign(eps,rname(idtlon_trn))
          crrn=rrname(idtlon_trn)+sign(eps,rrname(idtlon_trn))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'tlon_trn      ',8,101,&
     &                   stat)

          end if

        end if

        if(riname(idxdim_trn).ne.iname(idxdim_trn)) then

          call destroy('chkfile ',7,'cont',202,'xdim_trn      ',8,101,  &
     &                 stat)

        end if

        if(riname(idydim_trn).ne.iname(idydim_trn)) then

          call destroy('chkfile ',7,'cont',202,'ydim_trn      ',8,101,  &
     &                 stat)

        end if

        crn=rname(iddx_trn)+sign(eps,rname(iddx_trn))
        crrn=rrname(iddx_trn)+sign(eps,rrname(iddx_trn))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'dx_trn        ',6,101,  &
     &                 stat)

        end if

        crn=rname(iddy_trn)+sign(eps,rname(iddy_trn))
        crrn=rrname(iddy_trn)+sign(eps,rrname(iddy_trn))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'dy_trn        ',6,101,  &
     &                 stat)

        end if

        crn=rname(idulat_trn)+sign(eps,rname(idulat_trn))
        crrn=rrname(idulat_trn)+sign(eps,rrname(idulat_trn))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'ulat_trn      ',8,101,  &
     &                 stat)

        end if

        crn=rname(idulon_trn)+sign(eps,rname(idulon_trn))
        crrn=rrname(idulon_trn)+sign(eps,rrname(idulon_trn))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'ulon_trn      ',8,101,  &
     &                 stat)

        end if

        crn=rname(idriu_trn)+sign(eps,rname(idriu_trn))
        crrn=rrname(idriu_trn)+sign(eps,rrname(idriu_trn))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'riu_trn       ',7,101,  &
     &                 stat)

        end if

        crn=rname(idrju_trn)+sign(eps,rname(idrju_trn))
        crrn=rrname(idrju_trn)+sign(eps,rrname(idrju_trn))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'rju_trn       ',7,101,  &
     &                 stat)

        end if

        if(riname(idintopt_trn).ne.iname(idintopt_trn)) then

          call destroy('chkfile ',7,'cont',202,'intopt_trn    ',10,101, &
     &                 stat)

        end if

      end if

! -----

! Check the read out namelist variables from the restart or the surface
! data files.

      if(fname(1:3).eq.'res'.or.fname(1:3).eq.'sfc') then

        if(rcname(idsfcdat)(1:1).ne.cname(idsfcdat)(1:1)) then

          call destroy('chkfile ',7,'cont',202,'sfcdat        ',6,101,  &
     &                 stat)

        end if

        if(rcname(idsfcdat)(3:3).ne.cname(idsfcdat)(3:3)) then

          call destroy('chkfile ',7,'cont',202,'sfcdat        ',6,101,  &
     &                 stat)

        end if

        if(rcname(idsfcdat)(1:1).eq.'o') then

          if(riname(idmpopt_lnd).ne.iname(idmpopt_lnd)) then

            call destroy('chkfile ',7,'cont',202,'mpopt_lnd     ',9,101,&
     &                   stat)

          end if

          if(riname(idmpopt_lnd).eq.1.or.riname(idmpopt_lnd).eq.2.or.   &
     &      riname(idmpopt_lnd).eq.3.or.riname(idmpopt_lnd).eq.13) then

            if(riname(idnspol_lnd).ne.iname(idnspol_lnd)) then

              call destroy('chkfile ',7,'cont',202,'nspol_lnd     ',9,  &
     &                     101,stat)

            end if

            crn=rname(idtlat1_lnd)+sign(eps,rname(idtlat1_lnd))
            crrn=rrname(idtlat1_lnd)+sign(eps,rrname(idtlat1_lnd))

            if(abs(crrn/crn-1.e0).gt.chkeps) then

              call destroy('chkfile ',7,'cont',202,'tlat1_lnd     ',9,  &
     &                     101,stat)

            end if

          end if

          if(riname(idmpopt_lnd).eq.2) then

            crn=rname(idtlat2_lnd)+sign(eps,rname(idtlat2_lnd))
            crrn=rrname(idtlat2_lnd)+sign(eps,rrname(idtlat2_lnd))

            if(abs(crrn/crn-1.e0).gt.chkeps) then

              call destroy('chkfile ',7,'cont',202,'tlat2_lnd     ',9,  &
     &                     101,stat)

            end if

          end if

          if(riname(idmpopt_lnd).eq.1.or.                               &
     &       riname(idmpopt_lnd).eq.2.or.riname(idmpopt_lnd).eq.4) then

            crn=rname(idtlon_lnd)+sign(eps,rname(idtlon_lnd))
            crrn=rrname(idtlon_lnd)+sign(eps,rrname(idtlon_lnd))

            if(abs(crrn/crn-1.e0).gt.chkeps) then

              call destroy('chkfile ',7,'cont',202,'tlon_lnd      ',8,  &
     &                     101,stat)

            end if

          end if

          if(riname(idxdim_lnd).ne.iname(idxdim_lnd)) then

            call destroy('chkfile ',7,'cont',202,'xdim_lnd      ',8,101,&
     &                   stat)

          end if

          if(riname(idydim_lnd).ne.iname(idydim_lnd)) then

            call destroy('chkfile ',7,'cont',202,'ydim_lnd      ',8,101,&
     &                   stat)

          end if

          crn=rname(iddx_lnd)+sign(eps,rname(iddx_lnd))
          crrn=rrname(iddx_lnd)+sign(eps,rrname(iddx_lnd))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'dx_lnd        ',6,101,&
     &                   stat)

          end if

          crn=rname(iddy_lnd)+sign(eps,rname(iddy_lnd))
          crrn=rrname(iddy_lnd)+sign(eps,rrname(iddy_lnd))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'dy_lnd        ',6,101,&
     &                   stat)

          end if

          crn=rname(idulat_lnd)+sign(eps,rname(idulat_lnd))
          crrn=rrname(idulat_lnd)+sign(eps,rrname(idulat_lnd))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'ulat_lnd      ',8,101,&
     &                   stat)

          end if

          crn=rname(idulon_lnd)+sign(eps,rname(idulon_lnd))
          crrn=rrname(idulon_lnd)+sign(eps,rrname(idulon_lnd))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'ulon_lnd      ',8,101,&
     &                   stat)

          end if

          crn=rname(idriu_lnd)+sign(eps,rname(idriu_lnd))
          crrn=rrname(idriu_lnd)+sign(eps,rrname(idriu_lnd))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'riu_lnd       ',7,101,&
     &                   stat)

          end if

          crn=rname(idrju_lnd)+sign(eps,rname(idrju_lnd))
          crrn=rrname(idrju_lnd)+sign(eps,rrname(idrju_lnd))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'rju_lnd       ',7,101,&
     &                   stat)

          end if

          if(riname(idintopt_lnd).ne.iname(idintopt_lnd)) then

            call destroy('chkfile ',7,'cont',202,'intopt_lnd    ',10,   &
     &                   101,stat)

          end if

          if(riname(idnumctg_lnd).ne.iname(idnumctg_lnd)) then

            call destroy('chkfile ',7,'cont',202,'numctg_lnd    ',10,   &
     &                   101,stat)

          end if

        end if

        if(rcname(idsfcdat)(3:3).eq.'o') then

          if(riname(idmpopt_ice).ne.iname(idmpopt_ice)) then

            call destroy('chkfile ',7,'cont',202,'mpopt_ice     ',9,101,&
     &                   stat)

          end if

          if(riname(idmpopt_ice).eq.1.or.riname(idmpopt_ice).eq.2.or.   &
     &      riname(idmpopt_ice).eq.3.or.riname(idmpopt_ice).eq.13) then

            if(riname(idnspol_ice).ne.iname(idnspol_ice)) then

              call destroy('chkfile ',7,'cont',202,'nspol_ice     ',9,  &
     &                     101,stat)

            end if

            crn=rname(idtlat1_ice)+sign(eps,rname(idtlat1_ice))
            crrn=rrname(idtlat1_ice)+sign(eps,rrname(idtlat1_ice))

            if(abs(crrn/crn-1.e0).gt.chkeps) then

              call destroy('chkfile ',7,'cont',202,'tlat1_ice     ',9,  &
     &                     101,stat)

            end if

          end if

          if(riname(idmpopt_ice).eq.2) then

            crn=rname(idtlat2_ice)+sign(eps,rname(idtlat2_ice))
            crrn=rrname(idtlat2_ice)+sign(eps,rrname(idtlat2_ice))

            if(abs(crrn/crn-1.e0).gt.chkeps) then

              call destroy('chkfile ',7,'cont',202,'tlat2_ice     ',9,  &
     &                     101,stat)

            end if

          end if

          if(riname(idmpopt_ice).eq.1.or.                               &
     &       riname(idmpopt_ice).eq.2.or.riname(idmpopt_ice).eq.4) then

            crn=rname(idtlon_ice)+sign(eps,rname(idtlon_ice))
            crrn=rrname(idtlon_ice)+sign(eps,rrname(idtlon_ice))

            if(abs(crrn/crn-1.e0).gt.chkeps) then

              call destroy('chkfile ',7,'cont',202,'tlon_ice      ',8,  &
     &                     101,stat)

            end if

          end if

          if(riname(idxdim_ice).ne.iname(idxdim_ice)) then

            call destroy('chkfile ',7,'cont',202,'xdim_ice      ',8,101,&
     &                   stat)

          end if

          if(riname(idydim_ice).ne.iname(idydim_ice)) then

            call destroy('chkfile ',7,'cont',202,'ydim_ice      ',8,101,&
     &                   stat)

          end if

          crn=rname(iddx_ice)+sign(eps,rname(iddx_ice))
          crrn=rrname(iddx_ice)+sign(eps,rrname(iddx_ice))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'dx_ice        ',6,101,&
     &                   stat)

          end if

          crn=rname(iddy_ice)+sign(eps,rname(iddy_ice))
          crrn=rrname(iddy_ice)+sign(eps,rrname(iddy_ice))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'dy_ice        ',6,101,&
     &                   stat)

          end if

          crn=rname(idulat_ice)+sign(eps,rname(idulat_ice))
          crrn=rrname(idulat_ice)+sign(eps,rrname(idulat_ice))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'ulat_ice      ',8,101,&
     &                   stat)

          end if

          crn=rname(idulon_ice)+sign(eps,rname(idulon_ice))
          crrn=rrname(idulon_ice)+sign(eps,rrname(idulon_ice))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'ulon_ice      ',8,101,&
     &                   stat)

          end if

          crn=rname(idriu_ice)+sign(eps,rname(idriu_ice))
          crrn=rrname(idriu_ice)+sign(eps,rrname(idriu_ice))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'riu_ice       ',7,101,&
     &                   stat)

          end if

          crn=rname(idrju_ice)+sign(eps,rname(idrju_ice))
          crrn=rrname(idrju_ice)+sign(eps,rrname(idrju_ice))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'rju_ice       ',7,101,&
     &                   stat)

          end if

        end if

      end if

! -----

! Check the read out namelist variables from the restart or the sea
! surface temperature data files.

      if(fname(1:3).eq.'res'.or.fname(1:3).eq.'sst') then

        if(rcname(idsfcdat)(2:2).ne.cname(idsfcdat)(2:2)) then

          call destroy('chkfile ',7,'cont',202,'sfcdat        ',6,101,  &
     &                 stat)

        end if

        if(rcname(idsfcdat)(2:2).eq.'o') then

          if(riname(idsfcopt).eq.3.or.riname(idsfcopt).eq.13) then

            crn=rname(idsstitv)+sign(eps,rname(idsstitv))
            crrn=rrname(idsstitv)+sign(eps,rrname(idsstitv))

            if(abs(crrn/crn-1.e0).gt.chkeps) then

              call destroy('chkfile ',7,'cont',202,'sstitv        ',6,  &
     &                     101,stat)

            end if

          end if

          if(riname(idmpopt_sst).ne.iname(idmpopt_sst)) then

            call destroy('chkfile ',7,'cont',202,'mpopt_sst     ',9,101,&
     &                   stat)

          end if

          if(riname(idmpopt_sst).eq.1.or.riname(idmpopt_sst).eq.2.or.   &
     &      riname(idmpopt_sst).eq.3.or.riname(idmpopt_sst).eq.13) then

            if(riname(idnspol_sst).ne.iname(idnspol_sst)) then

              call destroy('chkfile ',7,'cont',202,'nspol_sst     ',9,  &
     &                     101,stat)

            end if

            crn=rname(idtlat1_sst)+sign(eps,rname(idtlat1_sst))
            crrn=rrname(idtlat1_sst)+sign(eps,rrname(idtlat1_sst))

            if(abs(crrn/crn-1.e0).gt.chkeps) then

              call destroy('chkfile ',7,'cont',202,'tlat1_sst     ',9,  &
     &                     101,stat)

            end if

          end if

          if(riname(idmpopt_sst).eq.2) then

            crn=rname(idtlat2_sst)+sign(eps,rname(idtlat2_sst))
            crrn=rrname(idtlat2_sst)+sign(eps,rrname(idtlat2_sst))

            if(abs(crrn/crn-1.e0).gt.chkeps) then

              call destroy('chkfile ',7,'cont',202,'tlat2_sst     ',9,  &
     &                     101,stat)

            end if

          end if

          if(riname(idmpopt_sst).eq.1.or.                               &
     &       riname(idmpopt_sst).eq.2.or.riname(idmpopt_sst).eq.4) then

            crn=rname(idtlon_sst)+sign(eps,rname(idtlon_sst))
            crrn=rrname(idtlon_sst)+sign(eps,rrname(idtlon_sst))

            if(abs(crrn/crn-1.e0).gt.chkeps) then

              call destroy('chkfile ',7,'cont',202,'tlon_sst      ',8,  &
     &                     101,stat)

            end if

          end if

          if(riname(idxdim_sst).ne.iname(idxdim_sst)) then

            call destroy('chkfile ',7,'cont',202,'xdim_sst      ',8,101,&
     &                   stat)

          end if

          if(riname(idydim_sst).ne.iname(idydim_sst)) then

            call destroy('chkfile ',7,'cont',202,'ydim_sst      ',8,101,&
     &                   stat)

          end if

          crn=rname(iddx_sst)+sign(eps,rname(iddx_sst))
          crrn=rrname(iddx_sst)+sign(eps,rrname(iddx_sst))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'dx_sst        ',6,101,&
     &                   stat)

          end if

          crn=rname(iddy_sst)+sign(eps,rname(iddy_sst))
          crrn=rrname(iddy_sst)+sign(eps,rrname(iddy_sst))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'dy_sst        ',6,101,&
     &                   stat)

          end if

          crn=rname(idulat_sst)+sign(eps,rname(idulat_sst))
          crrn=rrname(idulat_sst)+sign(eps,rrname(idulat_sst))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'ulat_sst      ',8,101,&
     &                   stat)

          end if

          crn=rname(idulon_sst)+sign(eps,rname(idulon_sst))
          crrn=rrname(idulon_sst)+sign(eps,rrname(idulon_sst))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'ulon_sst      ',8,101,&
     &                   stat)

          end if

          crn=rname(idriu_sst)+sign(eps,rname(idriu_sst))
          crrn=rrname(idriu_sst)+sign(eps,rrname(idriu_sst))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'riu_sst       ',7,101,&
     &                   stat)

          end if

          crn=rname(idrju_sst)+sign(eps,rname(idrju_sst))
          crrn=rrname(idrju_sst)+sign(eps,rrname(idrju_sst))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'rju_sst       ',7,101,&
     &                   stat)

          end if

        end if

      end if

! -----

! Check the read out namelist variables from the restart or the dumped
! files.

      if(fname(1:3).eq.'res'.or.fname(1:3).eq.'geo'                     &
     &  .or.fname(1:3).eq.'dmp'.or.fname(1:3).eq.'mon') then

        if(riname(iddmpfmt).ne.iname(iddmpfmt)) then

          call destroy('chkfile ',7,'cont',202,'dmpfmt        ',6,101,  &
     &                 stat)

        end if

        if(riname(iddmplev).ne.iname(iddmplev)) then

          call destroy('chkfile ',7,'cont',202,'dmplev        ',6,101,  &
     &                 stat)

        end if

        if(riname(iddmpmon).ne.iname(iddmpmon)) then

          call destroy('chkfile ',7,'cont',202,'dmpmon        ',6,101,  &
     &                 stat)

        end if

        if(rcname(iddmpvar)(1:18).ne.cname(iddmpvar)(1:18)) then

          call destroy('chkfile ',7,'cont',202,'dmpvar        ',6,101,  &
     &                 stat)

        end if

      end if

! -----

! Check the read out namelist variables from the restart files.

      if(fname(1:3).eq.'res') then

        if(riname(idwbc).ne.iname(idwbc)) then

          call destroy('chkfile ',7,'cont',202,'wbc           ',3,101,  &
     &                 stat)

        end if

        if(riname(idebc).ne.iname(idebc)) then

          call destroy('chkfile ',7,'cont',202,'ebc           ',3,101,  &
     &                 stat)

        end if

        if(riname(idsbc).ne.iname(idsbc)) then

          call destroy('chkfile ',7,'cont',202,'sbc           ',3,101,  &
     &                 stat)

        end if

        if(riname(idnbc).ne.iname(idnbc)) then

          call destroy('chkfile ',7,'cont',202,'nbc           ',3,101,  &
     &                 stat)

        end if

        if(riname(idbbc).ne.iname(idbbc)) then

          call destroy('chkfile ',7,'cont',202,'bbc           ',3,101,  &
     &                 stat)

        end if

        if(riname(idtbc).ne.iname(idtbc)) then

          call destroy('chkfile ',7,'cont',202,'tbc           ',3,101,  &
     &                 stat)

        end if

        if(riname(idwbc).ge.4.or.riname(idebc).ge.4                     &
     &    .or.riname(idsbc).ge.4.or.riname(idnbc).ge.4) then

          if(rcname(idlbcvar)(1:10).ne.cname(idlbcvar)(1:10)) then

            call destroy('chkfile ',7,'cont',202,'lbcvar        ',6,101,&
     &                   stat)

          end if

          crn=rname(idlbnews)+sign(eps,rname(idlbnews))
          crrn=rrname(idlbnews)+sign(eps,rrname(idlbnews))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'lbnews        ',6,101,&
     &                   stat)

          end if

          crn=rname(idlbnorm)+sign(eps,rname(idlbnorm))
          crrn=rrname(idlbnorm)+sign(eps,rrname(idlbnorm))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'lbnorm        ',6,101,&
     &                   stat)

          end if

        end if

        if(mod(riname(idwbc),10).ge.6.or.mod(riname(idebc),10).ge.6.or. &
     &    mod(riname(idsbc),10).ge.6.or.mod(riname(idnbc),10).ge.6) then

          crn=rname(idgwave)+sign(eps,rname(idgwave))
          crrn=rrname(idgwave)+sign(eps,rrname(idgwave))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'gwave         ',5,101,&
     &                   stat)

          end if

        end if

        if(riname(idgsmopt).ne.iname(idgsmopt)) then

          call destroy('chkfile ',7,'cont',202,'gsmopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idgsmopt).eq.1) then

          if(riname(idgsmcnt).ne.iname(idgsmcnt)) then

            call destroy('chkfile ',7,'cont',202,'gsmcnt        ',6,101,&
     &                   stat)

          end if

          crn=rname(idgsmcoe)+sign(eps,rname(idgsmcoe))
          crrn=rrname(idgsmcoe)+sign(eps,rrname(idgsmcoe))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'gsmcoe        ',6,101,&
     &                   stat)

          end if

        end if

        if(riname(idnggopt).ne.iname(idnggopt)) then

          call destroy('chkfile ',7,'cont',202,'nggopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idnggopt).eq.1) then

          if(rcname(idnggvar)(1:8).ne.cname(idnggvar)(1:8)) then

            call destroy('chkfile ',7,'cont',202,'nggvar        ',6,101,&
     &                   stat)

          end if

          crn=rname(idnggcoe)+sign(eps,rname(idnggcoe))
          crrn=rrname(idnggcoe)+sign(eps,rrname(idnggcoe))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'nggcoe        ',6,101,&
     &                   stat)

          end if

          crn=rname(idnggdlt)+sign(eps,rname(idnggdlt))
          crrn=rrname(idnggdlt)+sign(eps,rrname(idnggdlt))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'nggdlt        ',6,101,&
     &                   stat)

          end if

          crn=rname(idnggstr)+sign(eps,rname(idnggstr))
          crrn=rrname(idnggstr)+sign(eps,rrname(idnggstr))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'nggstr        ',6,101,&
     &                   stat)

          end if

          crn=rname(idnggend)+sign(eps,rname(idnggend))
          crrn=rrname(idnggend)+sign(eps,rrname(idnggend))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'nggend        ',6,101,&
     &                   stat)

          end if

          crn=rname(idnggc20)+sign(eps,rname(idnggc20))
          crrn=rrname(idnggc20)+sign(eps,rrname(idnggc20))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'nggc20        ',6,101,&
     &                   stat)

          end if

        end if

        if(riname(idexbopt).ne.iname(idexbopt)) then

          call destroy('chkfile ',7,'cont',202,'exbopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idexbopt).ge.1) then

          if(rcname(idexbvar)(1:8).ne.cname(idexbvar)(1:8)) then

            call destroy('chkfile ',7,'cont',202,'exbvar        ',6,101,&
     &                   stat)

          end if

          crn=rname(idexnews)+sign(eps,rname(idexnews))
          crrn=rrname(idexnews)+sign(eps,rrname(idexnews))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'exnews        ',6,101,&
     &                   stat)

          end if

          crn=rname(idexnorm)+sign(eps,rname(idexnorm))
          crrn=rrname(idexnorm)+sign(eps,rrname(idexnorm))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'exnorm        ',6,101,&
     &                   stat)

          end if

        end if

        if(riname(idlspopt).ne.iname(idlspopt)) then

          call destroy('chkfile ',7,'cont',202,'lspopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idlspopt).ge.1) then

          if(rcname(idlspvar)(1:10).ne.cname(idlspvar)(1:10)) then

            call destroy('chkfile ',7,'cont',202,'lspvar        ',6,101,&
     &                   stat)

          end if

          crn=rname(idlspsmt)+sign(eps,rname(idlspsmt))
          crrn=rrname(idlspsmt)+sign(eps,rrname(idlspsmt))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'lspsmt        ',6,101,&
     &                   stat)

          end if

          crn=rname(idlsnews)+sign(eps,rname(idlsnews))
          crrn=rrname(idlsnews)+sign(eps,rrname(idlsnews))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'lsnews        ',6,101,&
     &                   stat)

          end if

          crn=rname(idlsnorm)+sign(eps,rname(idlsnorm))
          crrn=rrname(idlsnorm)+sign(eps,rrname(idlsnorm))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'lsnorm        ',6,101,&
     &                   stat)

          end if

          if(riname(idwdnews).ne.iname(idwdnews)) then

            call destroy('chkfile ',7,'cont',202,'wdnews        ',6,101,&
     &                   stat)

          end if

          if(riname(idwdnorm).ne.iname(idwdnorm)) then

            call destroy('chkfile ',7,'cont',202,'wdnorm        ',6,101,&
     &                   stat)

          end if

        end if

        if(riname(idvspopt).ne.iname(idvspopt)) then

          call destroy('chkfile ',7,'cont',202,'vspopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idvspopt).ge.1) then

          if(rcname(idvspvar)(1:10).ne.cname(idvspvar)(1:10)) then

            call destroy('chkfile ',7,'cont',202,'vspvar        ',6,101,&
     &                   stat)

          end if

          if(riname(idvspopt).eq.1) then

            crn=rname(idvspgpv)+sign(eps,rname(idvspgpv))
            crrn=rrname(idvspgpv)+sign(eps,rrname(idvspgpv))

            if(abs(crrn/crn-1.e0).gt.chkeps) then

              call destroy('chkfile ',7,'cont',202,'vspgpv        ',6,  &
     &                     101,stat)

            end if

            crn=rname(idbotgpv)+sign(eps,rname(idbotgpv))
            crrn=rrname(idbotgpv)+sign(eps,rrname(idbotgpv))

            if(abs(crrn/crn-1.e0).gt.chkeps) then

              call destroy('chkfile ',7,'cont',202,'botgpv        ',6,  &
     &                     101,stat)

            end if

          end if

          crn=rname(idvspbar)+sign(eps,rname(idvspbar))
          crrn=rrname(idvspbar)+sign(eps,rrname(idvspbar))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'vspbar        ',6,101,&
     &                   stat)

          end if

          crn=rname(idbotbar)+sign(eps,rname(idbotbar))
          crrn=rrname(idbotbar)+sign(eps,rrname(idbotbar))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'botbar        ',6,101,&
     &                   stat)

          end if

        end if

        if(riname(idngropt).ne.iname(idngropt)) then

          call destroy('chkfile ',7,'cont',202,'ngropt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idngropt).ge.1) then

          if(rcname(idngrvar)(1:5).ne.cname(idngrvar)(1:5)) then

            call destroy('chkfile ',7,'cont',202,'ngrvar        ',6,101,&
     &                   stat)

          end if

          crn=rname(idngrcoe)+sign(eps,rname(idngrcoe))
          crrn=rrname(idngrcoe)+sign(eps,rrname(idngrcoe))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'ngrcoe        ',6,101,&
     &                   stat)

          end if

          crn=rname(idngrdlt)+sign(eps,rname(idngrdlt))
          crrn=rrname(idngrdlt)+sign(eps,rrname(idngrdlt))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'ngrdlt        ',6,101,&
     &                   stat)

          end if

          crn=rname(idngrstr)+sign(eps,rname(idngrstr))
          crrn=rrname(idngrstr)+sign(eps,rrname(idngrstr))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'ngrstr        ',6,101,&
     &                   stat)

          end if

          crn=rname(idngrend)+sign(eps,rname(idngrend))
          crrn=rrname(idngrend)+sign(eps,rrname(idngrend))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'ngrend        ',6,101,&
     &                   stat)

          end if

          if(riname(idngropt).eq.1) then

            crn=rname(idngrc20)+sign(eps,rname(idngrc20))
            crrn=rrname(idngrc20)+sign(eps,rrname(idngrc20))

            if(abs(crrn/crn-1.e0).gt.chkeps) then

              call destroy('chkfile ',7,'cont',202,'ngrc20        ',6,  &
     &                     101,stat)

            end if

          end if

          if(riname(idngropt).ge.2                                      &
     &      .or.rcname(idngrvar)(1:3).ne.'xxx') then

            crn=rname(idngraff)+sign(eps,rname(idngraff))
            crrn=rrname(idngraff)+sign(eps,rrname(idngraff))

            if(abs(crrn/crn-1.e0).gt.chkeps) then

              call destroy('chkfile ',7,'cont',202,'ngraff        ',6,  &
     &                     101,stat)

            end if

          end if

        end if

        if(riname(idsfcopt).ne.iname(idsfcopt)) then

          call destroy('chkfile ',7,'cont',202,'sfcopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idsfcopt).ge.1) then

          if(riname(idtubopt).eq.0) then

            if(riname(idlevpbl).ne.iname(idlevpbl)) then

              call destroy('chkfile ',7,'cont',202,'levpbl        ',6,  &
     &                     101,stat)

            end if

          end if

          if(riname(idlevund).ne.iname(idlevund)) then

            call destroy('chkfile ',7,'cont',202,'levund        ',6,101,&
     &                   stat)

          end if

          crn=rname(iddtgrd)+sign(eps,rname(iddtgrd))
          crrn=rrname(iddtgrd)+sign(eps,rrname(iddtgrd))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'dtgrd         ',5,101,&
     &                   stat)

          end if

          crn=rname(iddzgrd)+sign(eps,rname(iddzgrd))
          crrn=rrname(iddzgrd)+sign(eps,rrname(iddzgrd))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'dzgrd         ',5,101,&
     &                   stat)

          end if

          crn=rname(iddzsea)+sign(eps,rname(iddzsea))
          crrn=rrname(iddzsea)+sign(eps,rrname(iddzsea))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'dzsea         ',5,101,&
     &                   stat)

          end if

          crn=rname(idtgdeep)+sign(eps,rname(idtgdeep))
          crrn=rrname(idtgdeep)+sign(eps,rrname(idtgdeep))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'tgdeep        ',6,101,&
     &                   stat)

          end if

          if(rcname(idsfcdat)(1:1).eq.'x') then

            if(riname(idlnduse).ne.iname(idlnduse)) then

              call destroy('chkfile ',7,'cont',202,'lnduse        ',6,  &
     &                     101,stat)

            end if

            if(riname(idlnduse).eq.1) then

              crn=rname(idgralbe)+sign(eps,rname(idgralbe))
              crrn=rrname(idgralbe)+sign(eps,rrname(idgralbe))

              if(abs(crrn/crn-1.e0).gt.chkeps) then

                call destroy('chkfile ',7,'cont',202,'gralbe        ',6,&
     &                       101,stat)

              end if

              crn=rname(idgrbeta)+sign(eps,rname(idgrbeta))
              crrn=rrname(idgrbeta)+sign(eps,rrname(idgrbeta))

              if(abs(crrn/crn-1.e0).gt.chkeps) then

                call destroy('chkfile ',7,'cont',202,'grbeta        ',6,&
     &                       101,stat)

              end if

              crn=rname(idgrz0m)+sign(eps,rname(idgrz0m))
              crrn=rrname(idgrz0m)+sign(eps,rrname(idgrz0m))

              if(abs(crrn/crn-1.e0).gt.chkeps) then

                call destroy('chkfile ',7,'cont',202,'grz0m         ',5,&
     &                       101,stat)

              end if

              crn=rname(idgrz0h)+sign(eps,rname(idgrz0h))
              crrn=rrname(idgrz0h)+sign(eps,rrname(idgrz0h))

              if(abs(crrn/crn-1.e0).gt.chkeps) then

                call destroy('chkfile ',7,'cont',202,'grz0h         ',5,&
     &                       101,stat)

              end if

              crn=rname(idgrcap)+sign(eps,rname(idgrcap))
              crrn=rrname(idgrcap)+sign(eps,rrname(idgrcap))

              if(abs(crrn/crn-1.e0).gt.chkeps) then

                call destroy('chkfile ',7,'cont',202,'grcap         ',5,&
     &                       101,stat)

              end if

              crn=rname(idgrnuu)+sign(eps,rname(idgrnuu))
              crrn=rrname(idgrnuu)+sign(eps,rrname(idgrnuu))

              if(abs(crrn/crn-1.e0).gt.chkeps) then

                call destroy('chkfile ',7,'cont',202,'grnuu         ',5,&
     &                       101,stat)

              end if

            end if

          end if

          if(rcname(idsfcdat)(2:2).eq.'x') then

            crn=rname(idsstcst)+sign(eps,rname(idsstcst))
            crrn=rrname(idsstcst)+sign(eps,rrname(idsstcst))

            if(abs(crrn/crn-1.e0).gt.chkeps) then

              call destroy('chkfile ',7,'cont',202,'sstcst        ',6,  &
     &                     101,stat)

            end if

          end if

          if(rcname(idsfcdat)(3:3).eq.'o') then

            if(riname(iddstopt).ne.iname(iddstopt)) then

              call destroy('chkfile ',7,'cont',202,'dstopt        ',6,  &
     &                     101,stat)

            end if

          end if

        end if

        if(riname(idmovopt).ne.iname(idmovopt)) then

          call destroy('chkfile ',7,'cont',202,'movopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idmovopt).eq.1) then

          crn=rname(idumove)+sign(eps,rname(idumove))
          crrn=rrname(idumove)+sign(eps,rrname(idumove))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'umove         ',5,101,&
     &                   stat)

          end if

          crn=rname(idvmove)+sign(eps,rname(idvmove))
          crrn=rrname(idvmove)+sign(eps,rrname(idvmove))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'vmove         ',5,101,&
     &                   stat)

          end if

        end if

        if(riname(idgwmopt).ne.iname(idgwmopt)) then

          call destroy('chkfile ',7,'cont',202,'gwmopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idimpopt).ne.iname(idimpopt)) then

          call destroy('chkfile ',7,'cont',202,'impopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idimpopt).eq.1.or.riname(idimpopt).eq.2) then

          crn=rname(idweicoe)+sign(eps,rname(idweicoe))
          crrn=rrname(idweicoe)+sign(eps,rrname(idweicoe))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'weicoe        ',6,101,&
     &                   stat)

          end if

        end if

        if(riname(idimpopt).eq.2) then

          crn=rname(idgsdeps)+sign(eps,rname(idgsdeps))
          crrn=rrname(idgsdeps)+sign(eps,rrname(idgsdeps))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'gsdeps        ',6,101,&
     &                   stat)

          end if

        end if

        crn=rname(idfilcoe)+sign(eps,rname(idfilcoe))
        crrn=rrname(idfilcoe)+sign(eps,rrname(idfilcoe))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'filcoe        ',6,101,  &
     &                 stat)

        end if

        if(riname(idadvopt).ne.iname(idadvopt)) then

          call destroy('chkfile ',7,'cont',202,'advopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idsmtopt).ne.iname(idsmtopt)) then

          call destroy('chkfile ',7,'cont',202,'smtopt        ',6,101,  &
     &                 stat)

        end if

        if(mod(riname(idsmtopt),10).eq.1                                &
     &    .or.mod(riname(idsmtopt),10).eq.2                             &
     &    .or.mod(riname(idsmtopt),10).eq.3) then

          crn=rname(idsmhcoe)+sign(eps,rname(idsmhcoe))
          crrn=rrname(idsmhcoe)+sign(eps,rrname(idsmhcoe))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'smhcoe        ',6,101,&
     &                   stat)

          end if

          crn=rname(idsmvcoe)+sign(eps,rname(idsmvcoe))
          crrn=rrname(idsmvcoe)+sign(eps,rrname(idsmvcoe))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'smvcoe        ',6,101,&
     &                   stat)

          end if

        end if

        if(riname(idsmtopt).ge.11) then

          crn=rname(idnlhcoe)+sign(eps,rname(idnlhcoe))
          crrn=rrname(idnlhcoe)+sign(eps,rrname(idnlhcoe))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'nlhcoe        ',6,101,&
     &                   stat)

          end if

          crn=rname(idnlvcoe)+sign(eps,rname(idnlvcoe))
          crrn=rrname(idnlvcoe)+sign(eps,rrname(idnlvcoe))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'nlvcoe        ',6,101,&
     &                   stat)

          end if

        end if

        if(riname(idmfcopt).ne.iname(idmfcopt)) then

          call destroy('chkfile ',7,'cont',202,'mfcopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idcoropt).ne.iname(idcoropt)) then

          call destroy('chkfile ',7,'cont',202,'coropt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idcrvopt).ne.iname(idcrvopt)) then

          call destroy('chkfile ',7,'cont',202,'crvopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idbuyopt).ne.iname(idbuyopt)) then

          call destroy('chkfile ',7,'cont',202,'buyopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(iddiaopt).ne.iname(iddiaopt)) then

          call destroy('chkfile ',7,'cont',202,'diaopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(iddivopt).ne.iname(iddivopt)) then

          call destroy('chkfile ',7,'cont',202,'divopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idcphopt).ne.iname(idcphopt)) then

          call destroy('chkfile ',7,'cont',202,'cphopt        ',6,101,  &
     &                 stat)

        end if

        if(abs(riname(idcphopt)).lt.10) then

          if(abs(riname(idcphopt)).ge.2) then

            if(riname(idhaiopt).ne.iname(idhaiopt)) then

              call destroy('chkfile ',7,'cont',202,'haiopt        ',6,  &
     &                     101,stat)

            end if

          end if

        end if

        if(abs(riname(idcphopt)).ge.1) then

          crn=rname(idthresq)+sign(eps,rname(idthresq))
          crrn=rrname(idthresq)+sign(eps,rrname(idthresq))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'thresq        ',6,101,&
     &                   stat)

          end if

        end if

        if(abs(riname(idcphopt)).gt.10                                  &
     &    .and.abs(riname(idcphopt)).lt.20) then

          crn=rname(iddtcmph)+sign(eps,rname(iddtcmph))
          crrn=rrname(iddtcmph)+sign(eps,rrname(iddtcmph))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'dtcmph(1)     ',9,101,&
     &                   stat)

          end if

          crn=rname(iddtcmph+1)+sign(eps,rname(iddtcmph+1))
          crrn=rrname(iddtcmph+1)+sign(eps,rrname(iddtcmph+1))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'dtcmph(2)     ',9,101,&
     &                   stat)

          end if

          if(riname(idncbinw).ne.iname(idncbinw)) then

            call destroy('chkfile ',7,'cont',202,'ncbinw        ',6,101,&
     &                   stat)

          end if

          crn=rname(idbbinw)+sign(eps,rname(idbbinw))
          crrn=rrname(idbbinw)+sign(eps,rrname(idbbinw))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'bbinw         ',5,101,&
     &                   stat)

          end if

          crn=rname(idsbinw)+sign(eps,rname(idsbinw))
          crrn=rrname(idsbinw)+sign(eps,rrname(idsbinw))

          if(abs(crrn/crn-1.e0).gt.chkeps) then

            call destroy('chkfile ',7,'cont',202,'sbinw         ',5,101,&
     &                   stat)

          end if

        end if

        if(riname(idcphopt).lt.0) then

          if(riname(idqcgopt).ne.iname(idqcgopt)) then

            call destroy('chkfile ',7,'cont',202,'qcgopt        ',6,101,&
     &                   stat)

          end if

          if(abs(riname(idqcgopt)).ge.1) then

            crn=rname(ideledlt)+sign(eps,rname(ideledlt))
            crrn=rrname(ideledlt)+sign(eps,rrname(ideledlt))

            if(abs(crrn/crn-1.e0).gt.chkeps) then

              call destroy('chkfile ',7,'cont',202,'eledlt        ',6,  &
     &                     101,stat)

            end if

          end if

        end if

        if(riname(idaslopt).ne.iname(idaslopt)) then

          call destroy('chkfile ',7,'cont',202,'aslopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idaslopt).ge.1) then

          if(abs(riname(idcphopt)).lt.2                                 &
     &      .or.abs(riname(idcphopt)).ge.10) then

            call destroy('chkfile ',7,'cont',202,'cphopt        ',6,101,&
     &                   stat)

          end if

        end if

        if(riname(idtrkopt).ne.iname(idtrkopt)) then

          call destroy('chkfile ',7,'cont',202,'trkopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idtubopt).ne.iname(idtubopt)) then

          call destroy('chkfile ',7,'cont',202,'tubopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idtubopt).ge.1) then

          if(riname(idisoopt).ne.iname(idisoopt)) then

            call destroy('chkfile ',7,'cont',202,'isoopt        ',6,101,&
     &                   stat)

          end if

        end if

      end if

! -----

! Check the read out namelist variables from the GPV data files.

      if(fname(1:3).eq.'gpv') then

        if(riname(idnggopt).eq.0.and.riname(idexbopt).eq.0              &
     &    .and.riname(idlspopt).eq.0.and.riname(idvspopt).eq.0) then

          if(iname(idnggopt).ne.0) then

            call destroy('chkfile ',7,'cont',202,'nggopt        ',6,101,&
     &                   stat)

          end if

          if(iname(idexbopt).ne.0) then

            call destroy('chkfile ',7,'cont',202,'exbopt        ',6,101,&
     &                   stat)

          end if

          if(iname(idlspopt).ne.0) then

            call destroy('chkfile ',7,'cont',202,'lspopt        ',6,101,&
     &                   stat)

          end if

          if(iname(idvspopt).ne.0) then

            call destroy('chkfile ',7,'cont',202,'vspopt        ',6,101,&
     &                   stat)

          end if

        end if

        if(riname(idexbopt).eq.0) then

          if(riname(idwbc).ne.iname(idwbc)) then

            call destroy('chkfile ',7,'cont',202,'wbc           ',3,101,&
     &                   stat)

          end if

          if(riname(idebc).ne.iname(idebc)) then

            call destroy('chkfile ',7,'cont',202,'ebc           ',3,101,&
     &                   stat)

          end if

          if(riname(idsbc).ne.iname(idsbc)) then

            call destroy('chkfile ',7,'cont',202,'sbc           ',3,101,&
     &                   stat)

          end if

          if(riname(idnbc).ne.iname(idnbc)) then

            call destroy('chkfile ',7,'cont',202,'nbc           ',3,101,&
     &                   stat)

          end if

        end if

        if(riname(idexbopt).eq.2.or.riname(idexbopt).eq.12) then

          if(rcname(idsfcdat)(1:1).ne.cname(idsfcdat)(1:1)) then

            call destroy('chkfile ',7,'cont',202,'sfcdat        ',6,101,&
     &                   stat)

          end if

          if(riname(idsfcopt).ne.iname(idsfcopt)) then

            call destroy('chkfile ',7,'cont',202,'sfcopt        ',6,101,&
     &                   stat)

          end if

        end if

      end if

! -----

! Check the read out namelist variables from the aerosol data files.

      if(fname(1:3).eq.'asl') then

        if(riname(idnggopt).eq.0.and.riname(idexbopt).eq.0              &
     &    .and.riname(idlspopt).eq.0.and.riname(idvspopt).eq.0) then

          if(iname(idnggopt).ne.0) then

            call destroy('chkfile ',7,'cont',202,'nggopt        ',6,101,&
     &                   stat)

          end if

          if(iname(idexbopt).ne.0) then

            call destroy('chkfile ',7,'cont',202,'exbopt        ',6,101,&
     &                   stat)

          end if

          if(iname(idlspopt).ne.0) then

            call destroy('chkfile ',7,'cont',202,'lspopt        ',6,101,&
     &                   stat)

          end if

          if(iname(idvspopt).ne.0) then

            call destroy('chkfile ',7,'cont',202,'vspopt        ',6,101,&
     &                   stat)

          end if

        end if

        if(riname(idexbopt).eq.0) then

          if(riname(idwbc).ne.iname(idwbc)) then

            call destroy('chkfile ',7,'cont',202,'wbc           ',3,101,&
     &                   stat)

          end if

          if(riname(idebc).ne.iname(idebc)) then

            call destroy('chkfile ',7,'cont',202,'ebc           ',3,101,&
     &                   stat)

          end if

          if(riname(idsbc).ne.iname(idsbc)) then

            call destroy('chkfile ',7,'cont',202,'sbc           ',3,101,&
     &                   stat)

          end if

          if(riname(idnbc).ne.iname(idnbc)) then

            call destroy('chkfile ',7,'cont',202,'nbc           ',3,101,&
     &                   stat)

          end if

        end if

      end if

! -----

! Check the read out namelist variables from the surface data files.

      if(fname(1:3).eq.'sfc') then

        if(rcname(idsfcdat)(1:1).eq.'o'                                 &
     &    .and.cname(idsfcdat)(1:1).eq.'x') then

          call destroy('chkfile ',7,'cont',202,'sfcdat        ',6,101,  &
     &                 stat)

        end if

        if(rcname(idsfcdat)(3:3).eq.'o'                                 &
     &    .and.cname(idsfcdat)(3:3).eq.'x') then

          call destroy('chkfile ',7,'cont',202,'sfcdat        ',6,101,  &
     &                 stat)

        end if

      end if

! -----

! Check the read out namelist variables from the sea surface temperature
! data files.

      if(fname(1:3).eq.'sst') then

        if(rcname(idsfcdat)(2:2).eq.'o'                                 &
     &    .and.cname(idsfcdat)(2:2).eq.'x') then

          call destroy('chkfile ',7,'cont',202,'sfcdat        ',6,101,  &
     &                 stat)

        end if

      end if

! -----

! Check the read out namelist variables from the restart files in the
! case only soil and sea temperature are read.

      if(fname(1:3).eq.'und') then

        if(riname(idsfcopt).lt.1) then

          call destroy('chkfile ',7,'cont',202,'sfcopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idlevund).ne.iname(idlevund)) then

          call destroy('chkfile ',7,'cont',202,'levund        ',6,101,  &
     &                 stat)

        end if

        crn=rname(iddtgrd)+sign(eps,rname(iddtgrd))
        crrn=rrname(iddtgrd)+sign(eps,rrname(iddtgrd))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'dtgrd         ',5,101,  &
     &                 stat)

        end if

        crn=rname(iddzgrd)+sign(eps,rname(iddzgrd))
        crrn=rrname(iddzgrd)+sign(eps,rrname(iddzgrd))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'dzgrd         ',5,101,  &
     &                 stat)

        end if

        crn=rname(iddzsea)+sign(eps,rname(iddzsea))
        crrn=rrname(iddzsea)+sign(eps,rrname(iddzsea))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'dzsea         ',5,101,  &
     &                 stat)

        end if

        crn=rname(idtgdeep)+sign(eps,rname(idtgdeep))
        crrn=rrname(idtgdeep)+sign(eps,rrname(idtgdeep))

        if(abs(crrn/crn-1.e0).gt.chkeps) then

          call destroy('chkfile ',7,'cont',202,'tgdeep        ',6,101,  &
     &                 stat)

        end if

      end if

! -----

! Check the read out namelist variables from the dumped files.

      if(fname(1:3).eq.'geo'                                            &
     &  .or.fname(1:3).eq.'dmp'.or.fname(1:3).eq.'mon') then

        if(riname(idmfcopt).ne.iname(idmfcopt)) then

          call destroy('chkfile ',7,'cont',202,'mfcopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idcoropt).ne.iname(idcoropt)) then

          call destroy('chkfile ',7,'cont',202,'coropt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idsfcopt).ne.iname(idsfcopt)) then

          call destroy('chkfile ',7,'cont',202,'sfcopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idcphopt).ne.iname(idcphopt)) then

          call destroy('chkfile ',7,'cont',202,'cphopt        ',6,101,  &
     &                 stat)

        end if

        if(abs(riname(idcphopt)).lt.10) then

          if(abs(riname(idcphopt)).ge.2) then

            if(riname(idhaiopt).ne.iname(idhaiopt)) then

              call destroy('chkfile ',7,'cont',202,'haiopt        ',6,  &
     &                     101,stat)

            end if

          end if

        end if

        if(riname(idcphopt).lt.0) then

          if(riname(idqcgopt).ne.iname(idqcgopt)) then

            call destroy('chkfile ',7,'cont',202,'qcgopt        ',6,101,&
     &                   stat)

          end if

        end if

        if(riname(idaslopt).ne.iname(idaslopt)) then

          call destroy('chkfile ',7,'cont',202,'aslopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idtrkopt).ne.iname(idtrkopt)) then

          call destroy('chkfile ',7,'cont',202,'trkopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idtubopt).ne.iname(idtubopt)) then

          call destroy('chkfile ',7,'cont',202,'tubopt        ',6,101,  &
     &                 stat)

        end if

      end if

! -----

! Check the read out namelist variables from the restructed restart
! files.

      if(fname(1:3).eq.'rst') then

        if(riname(idtrnopt).ne.iname(idtrnopt)) then

          call destroy('chkfile ',7,'cont',202,'trnopt        ',6,101,  &
     &                 stat)

        end if

        if(rcname(ididate).ne.cname(ididate)) then

          call destroy('chkfile ',7,'cont',202,'idate         ',5,101,  &
     &                 stat)

        end if

        if(riname(idwbc).ne.iname(idwbc)) then

          call destroy('chkfile ',7,'cont',202,'wbc           ',3,101,  &
     &                 stat)

        end if

        if(riname(idebc).ne.iname(idebc)) then

          call destroy('chkfile ',7,'cont',202,'ebc           ',3,101,  &
     &                 stat)

        end if

        if(riname(idsbc).ne.iname(idsbc)) then

          call destroy('chkfile ',7,'cont',202,'sbc           ',3,101,  &
     &                 stat)

        end if

        if(riname(idnbc).ne.iname(idnbc)) then

          call destroy('chkfile ',7,'cont',202,'nbc           ',3,101,  &
     &                 stat)

        end if

        if(riname(idbbc).ne.iname(idbbc)) then

          call destroy('chkfile ',7,'cont',202,'bbc           ',3,101,  &
     &                 stat)

        end if

        if(riname(idtbc).ne.iname(idtbc)) then

          call destroy('chkfile ',7,'cont',202,'tbc           ',3,101,  &
     &                 stat)

        end if

        if(riname(idsfcopt).ne.iname(idsfcopt)) then

          call destroy('chkfile ',7,'cont',202,'sfcopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idsfcopt).ge.1) then

          if(riname(idlevund).ne.iname(idlevund)) then

            call destroy('chkfile ',7,'cont',202,'levund        ',6,101,&
     &                   stat)

          end if

        end if

        if(riname(idadvopt).ne.iname(idadvopt)) then

          call destroy('chkfile ',7,'cont',202,'advopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(iddiaopt).ne.iname(iddiaopt)) then

          call destroy('chkfile ',7,'cont',202,'diaopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idcphopt).ne.iname(idcphopt)) then

          call destroy('chkfile ',7,'cont',202,'cphopt        ',6,101,  &
     &                 stat)

        end if

        if(abs(riname(idcphopt)).lt.10) then

          if(abs(riname(idcphopt)).ge.2) then

            if(riname(idhaiopt).ne.iname(idhaiopt)) then

              call destroy('chkfile ',7,'cont',202,'haiopt        ',6,  &
     &                     101,stat)

            end if

          end if

        end if

        if(riname(idcphopt).lt.0) then

          if(riname(idqcgopt).ne.iname(idqcgopt)) then

            call destroy('chkfile ',7,'cont',202,'qcgopt        ',6,101,&
     &                   stat)

          end if

        end if

        if(riname(idaslopt).ne.iname(idaslopt)) then

          call destroy('chkfile ',7,'cont',202,'aslopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idtrkopt).ne.iname(idtrkopt)) then

          call destroy('chkfile ',7,'cont',202,'trkopt        ',6,101,  &
     &                 stat)

        end if

        if(riname(idtubopt).ne.iname(idtubopt)) then

          call destroy('chkfile ',7,'cont',202,'tubopt        ',6,101,  &
     &                 stat)

        end if

      end if

! -----

!! Check the read out namelist variables from the restart or the surface
!! data files in parallel processing.

      if(fname(1:3).eq.'res'.or.fname(1:3).eq.'sfc') then

        if(rcname(idsfcdat)(1:1).eq.'o') then

! Initialize the error descriptors.

          ierr_lnd=0
          ierr_albe=0
          ierr_beta=0
          ierr_z0m=0
          ierr_z0h=0
          ierr_cap=0
          ierr_nuu=0

! -----

! Perform checking.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(iid) reduction(+: ierr_lnd)

          do iid=0,riname(idnumctg_lnd)-1

            if(riname(idlnduse_lnd+iid).ne.iname(idlnduse_lnd+iid)) then
              ierr_lnd=ierr_lnd+1
            end if

          end do

!$omp end do

!$omp do schedule(runtime)                                              &
!$omp&   private(iid,crn_sub,crrn_sub) reduction(+: ierr_albe)

          do iid=0,riname(idnumctg_lnd)-1

            crn_sub=rname(idalbe_lnd+iid)                               &
     &        +sign(eps,rname(idalbe_lnd+iid))

            crrn_sub=rrname(idalbe_lnd+iid)                             &
     &        +sign(eps,rrname(idalbe_lnd+iid))

            if(abs(crrn_sub/crn_sub-1.e0).gt.chkeps) then
              ierr_albe=ierr_albe+1
            end if

          end do

!$omp end do

!$omp do schedule(runtime)                                              &
!$omp&   private(iid,crn_sub,crrn_sub) reduction(+: ierr_beta)

          do iid=0,riname(idnumctg_lnd)-1

            crn_sub=rname(idbeta_lnd+iid)                               &
     &        +sign(eps,rname(idbeta_lnd+iid))

            crrn_sub=rrname(idbeta_lnd+iid)                             &
     &        +sign(eps,rrname(idbeta_lnd+iid))

            if(abs(crrn_sub/crn_sub-1.e0).gt.chkeps) then
              ierr_beta=ierr_beta+1
            end if

          end do

!$omp end do

!$omp do schedule(runtime)                                              &
!$omp&   private(iid,crn_sub,crrn_sub) reduction(+: ierr_z0m)

          do iid=0,riname(idnumctg_lnd)-1

            crn_sub=rname(idz0m_lnd+iid)                                &
     &        +sign(eps,rname(idz0m_lnd+iid))

            crrn_sub=rrname(idz0m_lnd+iid)                              &
     &        +sign(eps,rrname(idz0m_lnd+iid))

            if(abs(crrn_sub/crn_sub-1.e0).gt.chkeps) then
              ierr_z0m=ierr_z0m+1
            end if

          end do

!$omp end do

!$omp do schedule(runtime)                                              &
!$omp&   private(iid,crn_sub,crrn_sub) reduction(+: ierr_z0h)

          do iid=0,riname(idnumctg_lnd)-1

            crn_sub=rname(idz0h_lnd+iid)                                &
     &        +sign(eps,rname(idz0h_lnd+iid))

            crrn_sub=rrname(idz0h_lnd+iid)                              &
     &        +sign(eps,rrname(idz0h_lnd+iid))

            if(abs(crrn_sub/crn_sub-1.e0).gt.chkeps) then
              ierr_z0h=ierr_z0h+1
            end if

          end do

!$omp end do

!$omp do schedule(runtime)                                              &
!$omp&   private(iid,crn_sub,crrn_sub) reduction(+: ierr_cap)

          do iid=0,riname(idnumctg_lnd)-1

            crn_sub=rname(idcap_lnd+iid)                                &
     &        +sign(eps,rname(idcap_lnd+iid))

            crrn_sub=rrname(idcap_lnd+iid)                              &
     &        +sign(eps,rrname(idcap_lnd+iid))

            if(abs(crrn_sub/crn_sub-1.e0).gt.chkeps) then
              ierr_cap=ierr_cap+1
            end if

          end do

!$omp end do

!$omp do schedule(runtime)                                              &
!$omp&   private(iid,crn_sub,crrn_sub) reduction(+: ierr_nuu)

          do iid=0,riname(idnumctg_lnd)-1

            crn_sub=rname(idnuu_lnd+iid)                                &
     &        +sign(eps,rname(idnuu_lnd+iid))

            crrn_sub=rrname(idnuu_lnd+iid)                              &
     &        +sign(eps,rrname(idnuu_lnd+iid))

            if(abs(crrn_sub/crn_sub-1.e0).gt.chkeps) then
              ierr_nuu=ierr_nuu+1
            end if

          end do

!$omp end do

!$omp end parallel

! -----

! If error occured, call the procedure destroy.

          if(ierr_lnd.ne.0) then

            call destroy('chkfile ',7,'cont',202,'lnduse_lnd    ',10,   &
     &                   101,stat)

          end if

          if(ierr_albe.ne.0) then

            call destroy('chkfile ',7,'cont',202,'albe_lnd      ',8,101,&
     &                   stat)

          end if

          if(ierr_beta.ne.0) then

            call destroy('chkfile ',7,'cont',202,'beta_lnd      ',8,101,&
     &                   stat)

          end if

          if(ierr_z0m.ne.0) then

            call destroy('chkfile ',7,'cont',202,'z0m_lnd       ',7,101,&
     &                   stat)

          end if

          if(ierr_z0h.ne.0) then

            call destroy('chkfile ',7,'cont',202,'z0h_lnd       ',7,101,&
     &                   stat)

          end if

          if(ierr_cap.ne.0) then

            call destroy('chkfile ',7,'cont',202,'cap_lnd       ',7,101,&
     &                   stat)

          end if

          if(ierr_nuu.ne.0) then

            call destroy('chkfile ',7,'cont',202,'nuu_lnd       ',7,101,&
     &                   stat)

          end if

! -----

        end if

      end if

!! -----

      end subroutine s_chkfile

!-----7--------------------------------------------------------------7--

      end module m_chkfile
