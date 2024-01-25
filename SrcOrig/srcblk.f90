!***********************************************************************
      module m_srcblk
!***********************************************************************

!     Author      : Sakakibara Atsushi, Naito Daisuke
!     Date        : 2003/05/19
!     Modification: 2003/12/12, 2004/03/22, 2004/04/01, 2004/05/31,
!                   2004/07/10, 2004/09/01, 2004/09/25, 2004/10/12,
!                   2004/12/17, 2005/01/07, 2005/04/04, 2005/08/05,
!                   2005/09/30, 2005/10/05, 2005/11/22, 2006/01/10,
!                   2006/02/13, 2006/04/03, 2006/07/21, 2006/09/30,
!                   2007/11/26, 2008/01/11, 2008/05/02, 2008/07/01,
!                   2008/08/25, 2008/10/10, 2009/02/27, 2011/03/18,
!                   2011/09/22, 2013/01/28, 2013/02/13, 2013/10/08

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the inferior procedures for source term of the bulk cold
!     rain cloud phycics.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_aggregat
      use m_charging
      use m_collect
      use m_convers
      use m_copy3d
      use m_depsit
      use m_distrpg
      use m_freezing
      use m_melting
      use m_more0q
      use m_newblk
      use m_newblk_noevap
      use m_nuc1stc
      use m_nuc1stv
      use m_nuc2nd
      use m_prodctwg
      use m_setblk
      use m_setcst3d
      use m_shedding
      use m_temparam

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: srcblk, s_srcblk

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface srcblk

        module procedure s_srcblk

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
      subroutine s_srcblk(cphopt,haiopt,qcgopt,dtb,thresq,              &
     &                    ni,nj,nk,ptbr,rbr,rbv,pi,p,ptpp,qvp,qcp,qrp,  &
     &                    qip,qsp,qgp,qhp,nccp,ncrp,ncip,ncsp,ncgp,nchp,&
     &                    qallp,urq,usq,ugq,uhq,urn,usn,ugn,uhn,ptpf,   &
     &                    qvf,qcf,qrf,qif,qsf,qgf,qhf,nccf,ncrf,ncif,   &
     &                    ncsf,ncgf,nchf,qccf,qrcf,qicf,qscf,qgcf,qhcf, &
     &                    t,tcel,qvsst0,qvsw,qvsi,lv,ls,lf,kp,mu,dv,mi, &
     &                    diaqc,diaqr,diaqi,diaqs,diaqg,vntr,vnts,vntg, &
     &                    nuvi,nuci,clcr,clcs,clcg,clri,clrs,clrg,      &
     &                    clir,clis,clig,clsr,clsg,clrsg,clrin,         &
     &                    clrsn,clsrn,clsgn,agcn,agrn,agin,agsn,        &
     &                    vdvr,vdvi,vdvs,vdvg,cncr,cnis,cnsg,cnsgn,     &
     &                    spsi,spgi,mlic,mlsr,mlgr,frrg,frrgn,          &
     &                    shsr,shgr,pgwet,ecs)
!***********************************************************************

! Input variables

      integer, intent(in) :: cphopt
                       ! Option for cloud micro physics

      integer, intent(in) :: haiopt
                       ! Option for additional hail processes

      integer, intent(in) :: qcgopt
                       ! Option for charging distribution

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: thresq
                       ! Minimum threshold value of mixing ratio

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rbv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of base state density

      real, intent(in) :: pi(0:ni+1,0:nj+1,1:nk)
                       ! Exnar function

      real, intent(in) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure

      real, intent(in) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(in) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(in) :: qcp(0:ni+1,0:nj+1,1:nk)
                       ! Cloud water mixing ratio at past

      real, intent(in) :: qrp(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixing ratio at past

      real, intent(in) :: qip(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixing ratio at past

      real, intent(in) :: qsp(0:ni+1,0:nj+1,1:nk)
                       ! Snow mixing ratio at past

      real, intent(in) :: qgp(0:ni+1,0:nj+1,1:nk)
                       ! Graupel mixing ratio at past

      real, intent(in) :: qhp(0:ni+1,0:nj+1,1:nk)
                       ! Hail mixing ratio at past

      real, intent(in) :: nccp(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud water at past

      real, intent(in) :: ncrp(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of rain water at past

      real, intent(in) :: ncip(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud ice at past

      real, intent(in) :: ncsp(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of snow at past

      real, intent(in) :: ncgp(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of graupel at past

      real, intent(in) :: nchp(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of hail at past

      real, intent(in) :: qallp(0:ni+1,0:nj+1,1:nk)
                       ! Total water and ice mixing ratio at past

      real, intent(in) :: urq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of rain water

      real, intent(in) :: usq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of snow

      real, intent(in) :: ugq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of graupel

      real, intent(in) :: uhq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of hail

      real, intent(in) :: urn(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of rain water concentrations

      real, intent(in) :: usn(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of snow concentrations

      real, intent(in) :: ugn(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of graupel concentrations

      real, intent(in) :: uhn(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of hail concentrations

! Input and output variables

      real, intent(inout) :: ptpf(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at future

      real, intent(inout) :: qvf(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at future

      real, intent(inout) :: qcf(0:ni+1,0:nj+1,1:nk)
                       ! Cloud water mixing ratio at future

      real, intent(inout) :: qrf(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixing ratio at future

      real, intent(inout) :: qif(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixing ratio at future

      real, intent(inout) :: qsf(0:ni+1,0:nj+1,1:nk)
                       ! Snow mixing ratio at future

      real, intent(inout) :: qgf(0:ni+1,0:nj+1,1:nk)
                       ! Graupel mixing ratio at future

      real, intent(inout) :: qhf(0:ni+1,0:nj+1,1:nk)
                       ! Hail mixing ratio at future

      real, intent(inout) :: nccf(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud water at future

      real, intent(inout) :: ncrf(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of rain water at future

      real, intent(inout) :: ncif(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud ice at future

      real, intent(inout) :: ncsf(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of snow at future

      real, intent(inout) :: ncgf(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of graupel at future

      real, intent(inout) :: nchf(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of hail at future

      real, intent(inout) :: qccf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of cloud water at future

      real, intent(inout) :: qrcf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of rain water at future

      real, intent(inout) :: qicf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of cloud ice at future

      real, intent(inout) :: qscf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of snow at future

      real, intent(inout) :: qgcf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of graupel at future

      real, intent(inout) :: qhcf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of hail at future

! Internal shared variables

      real, intent(inout) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

      real, intent(inout) :: tcel(0:ni+1,0:nj+1,1:nk)
                       ! Ambient air temperature

      real, intent(inout) :: qvsst0(0:ni+1,0:nj+1,1:nk)
                       ! Super saturation mixing ratio at melting point

      real, intent(inout) :: qvsw(0:ni+1,0:nj+1,1:nk)
                       ! Saturation mixing ratio for water

      real, intent(inout) :: qvsi(0:ni+1,0:nj+1,1:nk)
                       ! Saturation mixing ratio for ice

      real, intent(inout) :: lv(0:ni+1,0:nj+1,1:nk)
                       ! Latent heat of evapolation

      real, intent(inout) :: ls(0:ni+1,0:nj+1,1:nk)
                       ! Latent heat of sublimation

      real, intent(inout) :: lf(0:ni+1,0:nj+1,1:nk)
                       ! Latent heat of fusion

      real, intent(inout) :: kp(0:ni+1,0:nj+1,1:nk)
                       ! Thermal conductivity of air

      real, intent(inout) :: mu(0:ni+1,0:nj+1,1:nk)
                       ! Viscosity of air

      real, intent(inout) :: dv(0:ni+1,0:nj+1,1:nk)
                       ! Molecular diffusivity of water

      real, intent(inout) :: mi(0:ni+1,0:nj+1,1:nk)
                       ! Mean mass of cloud ice

      real, intent(inout) :: diaqc(0:ni+1,0:nj+1,1:nk)
                       ! Mean diameter of cloud water

      real, intent(inout) :: diaqr(0:ni+1,0:nj+1,1:nk)
                       ! Mean diameter of rain water

      real, intent(inout) :: diaqi(0:ni+1,0:nj+1,1:nk)
                       ! Mean diameter of cloud ice

      real, intent(inout) :: diaqs(0:ni+1,0:nj+1,1:nk)
                       ! Mean diameter of snow

      real, intent(inout) :: diaqg(0:ni+1,0:nj+1,1:nk)
                       ! Mean diameter of graupel

      real, intent(inout) :: vntr(0:ni+1,0:nj+1,1:nk)
                       ! Ventilation factor for rain water

      real, intent(inout) :: vnts(0:ni+1,0:nj+1,1:nk)
                       ! Ventilation factor for snow

      real, intent(inout) :: vntg(0:ni+1,0:nj+1,1:nk)
                       ! Ventilation factor for graupel

      real, intent(inout) :: nuvi(0:ni+1,0:nj+1,1:nk)
                       ! Nucleation rate of deposition or sorption

      real, intent(inout) :: nuci(0:ni+1,0:nj+1,1:nk)
                       ! Nucleation rate
                       ! of condensation, contact and homogeneous

      real, intent(inout) :: clcr(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between cloud water and rain water

      real, intent(inout) :: clcs(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between cloud water and snow

      real, intent(inout) :: clcg(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between cloud water and graupel

      real, intent(inout) :: clri(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between rain water and cloud ice

      real, intent(inout) :: clrs(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! from rain water to snow

      real, intent(inout) :: clrg(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between rain water and graupel

      real, intent(inout) :: clir(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between rain water and cloud ice

      real, intent(inout) :: clis(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between cloud ice and snow

      real, intent(inout) :: clig(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between cloud ice and graupel

      real, intent(inout) :: clsr(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! from snow to rain water

      real, intent(inout) :: clsg(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between snow and graupel

      real, intent(inout) :: clrsg(0:ni+1,0:nj+1,1:nk)
                       ! Production rate of graupel
                       ! from collection rate form rain to snow

      real, intent(inout) :: clrin(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate for concentrations
                       ! between rain water and cloud ice

      real, intent(inout) :: clrsn(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate for concentrations
                       ! between rain water and snow

      real, intent(inout) :: clsrn(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate for concentrations
                       ! between rain water and snow

      real, intent(inout) :: clsgn(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate for concentrations
                       ! between snow and graupel

      real, intent(inout) :: agcn(0:ni+1,0:nj+1,1:nk)
                       ! Aggregation rate for cloud water

      real, intent(inout) :: agrn(0:ni+1,0:nj+1,1:nk)
                       ! Aggregation rate for rain water

      real, intent(inout) :: agin(0:ni+1,0:nj+1,1:nk)
                       ! Aggregation rate for cloud ice

      real, intent(inout) :: agsn(0:ni+1,0:nj+1,1:nk)
                       ! Aggregation rate for snow

      real, intent(inout) :: vdvr(0:ni+1,0:nj+1,1:nk)
                       ! Evaporation rate from rain water to water vapor

      real, intent(inout) :: vdvi(0:ni+1,0:nj+1,1:nk)
                       ! Deposiotn rate from water vapor to cloud ice

      real, intent(inout) :: vdvs(0:ni+1,0:nj+1,1:nk)
                       ! Deposiotn rate from water vapor to snow

      real, intent(inout) :: vdvg(0:ni+1,0:nj+1,1:nk)
                       ! Deposiotn rate from water vapor to graupel

      real, intent(inout) :: cncr(0:ni+1,0:nj+1,1:nk)
                       ! Conversion rate from cloud water to rain water

      real, intent(inout) :: cnis(0:ni+1,0:nj+1,1:nk)
                       ! Conversion rate from cloud ice to snow

      real, intent(inout) :: cnsg(0:ni+1,0:nj+1,1:nk)
                       ! Conversion rate from snow to graupel

      real, intent(inout) :: cnsgn(0:ni+1,0:nj+1,1:nk)
                       ! Conversion rate for concentrations
                       ! from snow to graupel

      real, intent(inout) :: spsi(0:ni+1,0:nj+1,1:nk)
                       ! Secondary nucleation rate from snow

      real, intent(inout) :: spgi(0:ni+1,0:nj+1,1:nk)
                       ! Secondary nucleation rate from graupel

      real, intent(inout) :: mlic(0:ni+1,0:nj+1,1:nk)
                       ! Melting rate from cloud ice to cloud water

      real, intent(inout) :: mlsr(0:ni+1,0:nj+1,1:nk)
                       ! Melting rate from snow to rain water

      real, intent(inout) :: mlgr(0:ni+1,0:nj+1,1:nk)
                       ! Melting rate from graupel to rain water

      real, intent(inout) :: frrg(0:ni+1,0:nj+1,1:nk)
                       ! Freezing rate from rain water to graupel

      real, intent(inout) :: frrgn(0:ni+1,0:nj+1,1:nk)
                       ! Freezing rate for concentrations
                       ! from rain water to graupel

      real, intent(inout) :: shsr(0:ni+1,0:nj+1,1:nk)
                       ! Shedding rate of liquid water from snow

      real, intent(inout) :: shgr(0:ni+1,0:nj+1,1:nk)
                       ! Shedding rate of liquid water from graupel

      real, intent(inout) :: pgwet(0:ni+1,0:nj+1,1:nk)
                       ! Graupel produntion rate for moist process

      real, intent(inout) :: ecs(0:ni+1,0:nj+1,1:nk)
                       ! Collection efficiency of snow for cloud water

!-----7--------------------------------------------------------------7--

! Calculate the air temperature and suchlike.

      call s_setblk(cphopt,thresq,ni,nj,nk,ptbr,rbr,rbv,pi,p,ptpp,qvp,  &
     &              qcp,qrp,qip,qsp,qgp,nccp,ncrp,ncip,ncsp,ncgp,       &
     &              t,tcel,qvsst0,qvsw,qvsi,lv,ls,lf,kp,mu,dv,mi,       &
     &              diaqc,diaqr,diaqi,diaqs,diaqg,vntr,vnts,vntg,       &
     &              ecs(0,0,1))

! -----

! Calculate the nucleation rate of the deposition or sorption.

      call nuc1stv(thresq,ni,nj,nk,rbv,qvp,qip,t,qvsi,nuvi)

! -----

! Calculate the nucleation rate of the condensation, contact and
! homogeneous.

      call nuc1stc(dtb,thresq,ni,nj,nk,p,qcp,nccp,t,tcel,lv,kp,mu,diaqc,&
     &             nuci)

! -----

! Calculate the collection rate.

      call collect(cphopt,dtb,thresq,ni,nj,nk,rbr,qcp,qrp,qip,qsp,qgp,  &
     &            ncrp,ncsp,ncgp,urq,usq,ugq,urn,usn,ugn,t,mu,mi,       &
     &            diaqc,diaqr,diaqs,diaqg,clcr,clcs,clcg,clri,clrs,clrg,&
     &            clir,clis,clig,clsr,clsg,clrin,clrsn,clsrn,clsgn,ecs)

! -----

! Calculate the graupel production rate.

      call prodctwg(cphopt,dtb,thresq,ni,nj,nk,rbr,                     &
     &             qip,qsp,qgp,ncsp,tcel,qvsst0,lv,lf,kp,dv,vntg,       &
     &             clcg,clrg,clir,clis,clig,clsr,clsg,clsrn,clsgn,pgwet)

! -----

! Calculate the distribution ratio at which the collisions between rain
! water and snow and reset the collection rate.

      call distrpg(cphopt,thresq,ni,nj,nk,qsp,t,diaqr,diaqs,clrs,clsr,  &
     &             clrsn,clsrn,clrsg)

! -----

! Calculate the aggregation rate for the cloud ice and snow.

      call aggregat(cphopt,dtb,thresq,ni,nj,nk,rbr,rbv,qcp,qrp,qip,qsp, &
     &              nccp,ncrp,ncip,ncsp,diaqc,diaqr,agcn,agrn,agin,agsn)

! -----

! Calculate the melting rate.

      call melting(dtb,thresq,ni,nj,nk,rbr,rbv,                         &
     &             qip,qsp,qgp,tcel,qvsst0,lv,lf,kp,dv,                 &
     &             vnts,vntg,clcs,clcg,clrs,clrg,mlic,mlsr,mlgr)

! -----

! Calculate the evaporation and deposition rate.

      call depsit(dtb,thresq,ni,nj,nk,rbr,rbv,qvp,qrp,qip,qsp,qgp,ncip, &
     &            t,tcel,qvsst0,qvsw,qvsi,lv,ls,lf,kp,dv,mi,            &
     &            vntr,vnts,vntg,clcs,clcg,mlsr,mlgr,                   &
     &            vdvr,vdvi,vdvs,vdvg)

! -----

! Calculate the conversion rate.

      call convers(cphopt,dtb,thresq,ni,nj,nk,rbr,rbv,qcp,qip,qsp,      &
     &             nccp,ncsp,t,mu,mi,diaqi,diaqs,clcs,vdvi,vdvs,        &
     &             ecs,cncr,cnis,cnsg,cnsgn)

! -----

! Calculate the secondary nucleation rate.

      if(nuc2nd_opt.eq.1) then

        call nuc2nd(ni,nj,nk,rbv,t,clcs,clcg,pgwet,spsi,spgi)

      else

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,spsi)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,spgi)

      end if

! -----

! Calculate the freezing rate.

      call freezing(cphopt,dtb,thresq,ni,nj,nk,qrp,ncrp,tcel,diaqr,     &
     &              frrg,frrgn)

! -----

! Calculate the shedding rate.

      call shedding(thresq,ni,nj,nk,qsp,qgp,t,clcs,clcg,clrs,clrg,clig, &
     &              clsg,pgwet,shsr,shgr)

! -----

! Force the mixing ratio to not being less than 0.

      call more0q(thresq,ni,nj,nk,qcp,qrp,qip,qsp,qgp,                  &
     &            qcf,qrf,qif,qsf,qgf,nuvi,nuci,clcr,clcs,clcg,         &
     &            clri,clrs,clrg,clir,clis,clig,clsr,clsg,clrsg,        &
     &            vdvr,vdvi,vdvs,vdvg,cncr,cnis,cnsg,spsi,spgi,         &
     &            mlic,mlsr,mlgr,frrg,shsr,shgr)

! -----

! Solve the new potential temperature perturbation, the mixing ratio and
! concentrations.

      if(evapor_opt.ge.0) then

        call newblk(cphopt,thresq,ni,nj,nk,pi,qcp,qrp,qip,qsp,qgp,      &
     &            nccp,ncsp,ncgp,lv,ls,lf,mi,nuvi,nuci,clcr,clcs,clcg,  &
     &            clri,clrs,clrg,clir,clis,clig,clsr,clsg,clrsg,        &
     &            clrin,clrsn,clsrn,clsgn,agcn,agrn,agin,agsn,          &
     &            vdvr,vdvi,vdvs,vdvg,cncr,cnis,cnsg,cnsgn,             &
     &            spsi,spgi,mlic,mlsr,mlgr,frrg,frrgn,shsr,shgr,        &
     &            ptpf,qvf,qcf,qrf,qif,qsf,qgf,nccf,ncrf,ncif,ncsf,ncgf)

      else

        call newblk_noevap(cphopt,thresq,                               &
     &            ni,nj,nk,pi,qcp,qrp,qip,qsp,qgp,                      &
     &            nccp,ncsp,ncgp,lv,ls,lf,mi,nuvi,nuci,clcr,clcs,clcg,  &
     &            clri,clrs,clrg,clir,clis,clig,clsr,clsg,clrsg,        &
     &            clrin,clrsn,clsrn,clsgn,agcn,agrn,agin,agsn,          &
     &            vdvr,vdvi,vdvs,vdvg,cncr,cnis,cnsg,cnsgn,             &
     &            spsi,spgi,mlic,mlsr,mlgr,frrg,frrgn,shsr,shgr,        &
     &            ptpf,qvf,qcf,qrf,qif,qsf,qgf,nccf,ncrf,ncif,ncsf,ncgf)

      end if

! -----

! Solve the new charging distributions.

      if(cphopt.lt.0) then

        call charging(haiopt,qcgopt,dtb,thresq,                         &
     &                ni,nj,nk,qcp,qrp,qip,qsp,qgp,qhp,                 &
     &                qallp,qccf,qrcf,qicf,qscf,qgcf,qhcf)

      end if

! -----

! Now underconstruction, set array with zero or use copy procedure
! to avoid warning.

      if(haiopt.eq.1) then

        call copy3d(0,ni+1,0,nj+1,1,nk,qhp,qhf)
        call copy3d(0,ni+1,0,nj+1,1,nk,nchp,nchf)

        call copy3d(0,ni+1,0,nj+1,1,nk,uhq,qhf)
        call copy3d(0,ni+1,0,nj+1,1,nk,uhn,nchf)

        if(cphopt.lt.0) then

          call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qhcf)

        end if

      end if

! -----

      end subroutine s_srcblk

!-----7--------------------------------------------------------------7--

      end module m_srcblk
