!***********************************************************************
      module m_ini0name
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/01/20
!     Modification: 2007/01/31, 2007/03/23, 2007/04/24, 2007/05/21,
!                   2007/07/30, 2007/10/19, 2008/01/11, 2008/04/17,
!                   2008/05/02, 2008/07/01, 2008/08/25, 2008/10/10,
!                   2008/12/11, 2009/01/05, 2009/01/30, 2009/02/27,
!                   2011/05/16, 2011/08/09, 2011/08/18, 2011/09/22,
!                   2011/11/10, 2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     initialize the namelist variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_defname
      use m_inichar

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: ini0name, s_ini0name

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface ini0name

        module procedure s_ini0name

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
      subroutine s_ini0name
!***********************************************************************

! Internal private variable

      integer iid      ! Index of do loops

!-----7--------------------------------------------------------------7--

!!! Initialize the namelist variables.

!! For user specified variables.

! For the character variables.

      call inichar(exprim)
      call inichar(crsdir)
      call inichar(datdir)
      call inichar(lbcvar)
      call inichar(gpvvar)
      call inichar(nggvar)
      call inichar(exbvar)
      call inichar(lspvar)
      call inichar(vspvar)
      call inichar(rdrvar)
      call inichar(ngrvar)
      call inichar(sfcdat)
      call inichar(prvres)
      call inichar(sndtyp)
      call inichar(dmpvar)
      call inichar(mxnvar)
      call inichar(datype_gpv)
      call inichar(etrvar_gpv)
      call inichar(datype_rdr)
      call inichar(fltyp_uni)

      write(sfcast(1:22),'(a22)') '                      '

! -----

! For the integer variables.

      savmem=0

      numpe=0

      xgroup=0
      ygroup=0

      xsub=0
      ysub=0

      xdim=0
      ydim=0
      zdim=0

      mpopt=0
      nspol=0

      sthopt=0

      trnopt=0

      wbc=0
      ebc=0
      sbc=0
      nbc=0
      bbc=0
      tbc=0

      gsmopt=0
      gsmcnt=0
      nggopt=0
      exbopt=0
      exbwid=0
      lspopt=0
      wdnews=0
      wdnorm=0
      vspopt=0

      ngropt=0

      sfcopt=0
      levpbl=0
      levund=0
      lnduse=0
      dstopt=0

      iniopt=0
      snddim=0
      masopt=0

      movopt=0

      pt0opt=0
      pt0num=0

      gwmopt=0
      impopt=0
      advopt=0

      smtopt=0

      mfcopt=0

      coropt=0

      crvopt=0

      buyopt=0

      diaopt=0

      divopt=0

      cphopt=0
      haiopt=0
      ncbinw=0
      qcgopt=0

      aslopt=0

      trkopt=0
      qt0opt=0
      qt0num=0

      tubopt=0
      isoopt=0

      dmpfmt=0
      dmplev=0
      dmpmon=0
      resopt=0
      mxnopt=0

      mpopt_gpv=0
      nspol_gpv=0

      xdim_gpv=0
      ydim_gpv=0
      zdim_gpv=0

      intopt_gpv=0

      rotopt_gpv=0

      refsfc_gpv=0

      mpopt_asl=0
      nspol_asl=0

      xdim_asl=0
      ydim_asl=0
      zdim_asl=0

      intopt_asl=0

      mpopt_rdr=0
      nspol_rdr=0

      xdim_rdr=0
      ydim_rdr=0
      zdim_rdr=0

      rotopt_rdr=0

      mpopt_trn=0
      nspol_trn=0

      xdim_trn=0
      ydim_trn=0

      intopt_trn=0

      mpopt_lnd=0
      nspol_lnd=0

      xdim_lnd=0
      ydim_lnd=0

      intopt_lnd=0
      numctg_lnd=0

      mpopt_sst=0
      nspol_sst=0

      xdim_sst=0
      ydim_sst=0

      mpopt_ice=0
      nspol_ice=0

      xdim_ice=0
      ydim_ice=0

      rmopt_uni=0
      uniopt_uni=0
      ugroup_uni=0

      xsub_rst=0
      ysub_rst=0
      rmopt_rst=0

! -----

! For the real variables.

      tlat1=0.e0
      tlat2=0.e0
      tlon=0.e0
      disr=0.e0

      dx=0.e0
      dy=0.e0
      dz=0.e0
      ulat=0.e0
      ulon=0.e0
      riu=0.e0
      rju=0.e0

      zsfc=0.e0
      zflat=0.e0
      dzmin=0.e0
      layer1=0.e0
      layer2=0.e0

      mnthgh(1)=0.e0
      mnthgh(2)=0.e0

      mntwx=0.e0
      mntwy=0.e0
      mntcx=0.e0
      mntcy=0.e0

      stime=0.e0
      etime=0.e0

      lbnews=0.e0
      lbnorm=0.e0

      gwave=0.e0

      gpvitv=0.e0
      aslitv=0.e0
      gsmcoe=0.e0
      nggcoe=0.e0
      nggdlt=0.e0
      nggstr=0.e0
      nggend=0.e0
      nggc20=0.e0

      exnews=0.e0
      exnorm=0.e0

      lspsmt=0.e0
      lsnews=0.e0
      lsnorm=0.e0

      vspgpv=0.e0
      vspbar=0.e0
      botgpv=0.e0
      botbar=0.e0

      rdritv=0.e0
      ngrcoe=0.e0
      ngrdlt=0.e0
      ngrstr=0.e0
      ngrend=0.e0
      ngrc20=0.e0
      ngraff=0.e0

      dtgrd=0.e0
      dzgrd=0.e0
      dzsea=0.e0
      tgdeep=0.e0
      gralbe=0.e0
      grbeta=0.e0
      grz0m=0.e0
      grz0h=0.e0
      grcap=0.e0
      grnuu=0.e0
      sstcst=0.e0
      sstitv=0.e0

      zsnd0=0.e0
      psnd0=0.e0
      maseps=0.e0
      alpha1=0.e0
      alpha2=0.e0

      umove=0.e0
      vmove=0.e0

      ptp0=0.e0
      pt0rx=0.e0
      pt0ry=0.e0
      pt0rz=0.e0
      pt0cx=0.e0
      pt0cy=0.e0
      pt0cz=0.e0
      pt0ds=0.e0

      dtbig=0.e0
      dtsml=0.e0
      gsdeps=0.e0
      weicoe=0.e0
      filcoe=0.e0
      dtvcul=0.e0

      smhcoe=0.e0
      smvcoe=0.e0
      nlhcoe=0.e0
      nlvcoe=0.e0

      thresq=0.e0
      dtcmph(1)=0.e0
      dtcmph(2)=0.e0
      dtcmph(3)=0.e0
      dtcmph(4)=0.e0
      dtcmph(5)=0.e0
      dtcmph(6)=0.e0
      bbinw=0.e0
      sbinw=0.e0
      eledlt=0.e0

      qt0=0.e0
      qt0rx=0.e0
      qt0ry=0.e0
      qt0rz=0.e0
      qt0cx=0.e0
      qt0cy=0.e0
      qt0cz=0.e0
      qt0ds=0.e0
      qtdt=0.e0
      qtstr=0.e0
      qtend=0.e0

      dmpitv=0.e0
      monitv=0.e0
      resitv=0.e0
      mxnitv=0.e0

      tlat1_gpv=0.e0
      tlat2_gpv=0.e0
      tlon_gpv=0.e0

      dx_gpv=0.e0
      dy_gpv=0.e0

      ulat_gpv=0.e0
      ulon_gpv=0.e0

      riu_gpv=0.e0
      rju_gpv=0.e0

      tlat1_asl=0.e0
      tlat2_asl=0.e0
      tlon_asl=0.e0

      dx_asl=0.e0
      dy_asl=0.e0

      ulat_asl=0.e0
      ulon_asl=0.e0

      riu_asl=0.e0
      rju_asl=0.e0

      tlat1_rdr=0.e0
      tlat2_rdr=0.e0
      tlon_rdr=0.e0

      dx_rdr=0.e0
      dy_rdr=0.e0

      ulat_rdr=0.e0
      ulon_rdr=0.e0

      riu_rdr=0.e0
      rju_rdr=0.e0

      rdrcoe_rdr=0.e0
      rdrexp_rdr=0.e0

      tlat1_trn=0.e0
      tlat2_trn=0.e0
      tlon_trn=0.e0

      dx_trn=0.e0
      dy_trn=0.e0

      ulat_trn=0.e0
      ulon_trn=0.e0

      riu_trn=0.e0
      rju_trn=0.e0

      tlat1_lnd=0.e0
      tlat2_lnd=0.e0
      tlon_lnd=0.e0

      dx_lnd=0.e0
      dy_lnd=0.e0

      ulat_lnd=0.e0
      ulon_lnd=0.e0

      riu_lnd=0.e0
      rju_lnd=0.e0

      tlat1_sst=0.e0
      tlat2_sst=0.e0
      tlon_sst=0.e0

      dx_sst=0.e0
      dy_sst=0.e0

      ulat_sst=0.e0
      ulon_sst=0.e0

      riu_sst=0.e0
      rju_sst=0.e0

      tlat1_ice=0.e0
      tlat2_ice=0.e0
      tlon_ice=0.e0

      dx_ice=0.e0
      dy_ice=0.e0

      ulat_ice=0.e0
      ulon_ice=0.e0

      riu_ice=0.e0
      rju_ice=0.e0

      flitv_uni=0.e0

      flitv_rst=0.e0

! -----

! For the table.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(iid)

      do iid=1,100
        lnduse_lnd(iid)=0

        albe_lnd(iid)=0.e0
        beta_lnd(iid)=0.e0
        z0m_lnd(iid)=0.e0
        z0h_lnd(iid)=0.e0
        cap_lnd(iid)=0.e0
        nuu_lnd(iid)=0.e0

      end do

!$omp end do

!$omp end parallel

! -----

!! -----

! For subordinate variables.

      write(idate(1:12),'(a12)') '            '

      wlngth=0

      ncexp=0
      nccrs=0
      ncdat=0
      ncprv=0

      iwest=0
      ieast=0
      jsouth=0
      jnorth=0

      dxiv=0.e0
      dyiv=0.e0
      dziv=0.e0

      dxiv_gpv=0.e0
      dyiv_gpv=0.e0

      dxiv_asl=0.e0
      dyiv_asl=0.e0

      dxiv_rdr=0.e0
      dyiv_rdr=0.e0

      dxiv_trn=0.e0
      dyiv_trn=0.e0

      dxiv_lnd=0.e0
      dyiv_lnd=0.e0

      dxiv_sst=0.e0
      dyiv_sst=0.e0

      dxiv_ice=0.e0
      dyiv_ice=0.e0

! -----

!!! -----

      end subroutine s_ini0name

!-----7--------------------------------------------------------------7--

      end module m_ini0name
