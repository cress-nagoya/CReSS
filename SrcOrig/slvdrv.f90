!***********************************************************************
      module m_slvdrv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/01/25, 1999/03/25, 1999/04/06, 1999/05/10,
!                   1999/05/20, 1999/06/07, 1999/06/14, 1999/06/21,
!                   1999/07/05, 1999/07/28, 1999/08/03, 1999/08/23,
!                   1999/09/06, 1999/09/16, 1999/09/30, 1999/10/12,
!                   1999/10/22, 1999/11/01, 1999/11/19, 1999/11/24,
!                   1999/12/17, 2000/01/05, 2000/01/17, 2000/03/08,
!                   2000/03/23, 2000/04/18, 2000/06/01, 2000/07/05,
!                   2000/12/18, 2001/01/15, 2001/02/13, 2001/03/13,
!                   2001/04/15, 2001/05/29, 2001/06/29, 2001/07/13,
!                   2001/08/07, 2001/09/13, 2001/10/18, 2001/11/20,
!                   2001/12/11, 2002/01/15, 2002/04/02, 2002/07/03,
!                   2002/07/15, 2002/07/23, 2002/08/15, 2002/09/09,
!                   2002/10/15, 2002/10/31, 2002/12/02, 2003/01/04,
!                   2003/01/20, 2003/02/13, 2003/03/13, 2003/03/21,
!                   2003/04/30, 2003/05/19, 2003/07/15, 2003/09/01,
!                   2003/10/06, 2003/11/05, 2003/12/12, 2004/01/09,
!                   2004/02/01, 2004/03/05, 2004/04/01, 2004/04/15,
!                   2004/05/07, 2004/05/31, 2004/06/10, 2004/07/01,
!                   2004/08/01, 2004/08/20, 2004/09/25, 2004/12/17,
!                   2005/01/14, 2005/04/04, 2005/08/05, 2006/01/10,
!                   2006/02/13, 2006/04/03, 2006/06/21, 2006/07/21,
!                   2006/09/21, 2006/09/30, 2006/11/06, 2006/11/27,
!                   2007/01/20, 2007/03/10, 2007/04/11, 2007/04/24,
!                   2007/05/07, 2007/05/14, 2007/05/21, 2007/07/30,
!                   2008/01/11, 2008/04/17, 2008/05/02, 2008/06/09,
!                   2008/07/01, 2008/07/25, 2008/08/25, 2008/10/10,
!                   2008/12/11, 2009/01/05, 2009/01/30, 2009/02/27,
!                   2009/03/23, 2010/02/01, 2010/05/17, 2011/08/18,
!                   2011/09/22, 2011/11/10, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the inferior procedures for the time steps integration.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_abortslv
      use m_chkstd
      use m_closedmp
      use m_comindx
      use m_comkind
      use m_commpi
      use m_getcname
      use m_getdate
      use m_getiname
      use m_gettime
      use m_inichar
      use m_instvel
      use m_intgdrv
      use m_masscon
      use m_mxndrv
      use m_ndgstep
      use m_opendmp
      use m_outcheck
      use m_outdmp
      use m_outres
      use m_outstd05
      use m_outstd06
      use m_rdaslnxt
      use m_rdgpvnxt
      use m_rdrdrmrk
      use m_rdrdrnxt
      use m_rdsstnxt
      use m_setrdr
      use m_sfcphy
      use m_slvstep
      use m_timeint

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: slvdrv, s_slvdrv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface slvdrv

        module procedure s_slvdrv

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
      subroutine s_slvdrv(fpidate,fpnggopt,fpexbopt,fplspopt,fpvspopt,  &
     &                    fpngropt,fpsfcopt,fpmasopt,fpadvopt,fpaslopt, &
     &                    fptubopt,fpresopt,fpmxnopt,fmois,ksp0,area,   &
     &                    ni,nj,nk,nqw,nnw,nqi,nni,km,nqa,nund,zph,     &
     &                    zsth,lat,lon,j31,j32,jcb,jcb8u,jcb8v,jcb8w,   &
     &                    mf,mf8u,mf8v,rmf,rmf8u,rmf8v,fc,ubr,vbr,      &
     &                    pbr,ptbr,qvbr,rbr,rst,rst8u,rst8v,rst8w,rcsq, &
     &                    rbcx,rbcy,rbcxy,rbct,land,albe,beta,cap,nuu,  &
     &                    kai,u,up,v,vp,w,wp,pp,ppp,ptp,ptpp,qv,qvp,    &
     &                    qwtr,qwtrp,nwtr,nwtrp,qice,qicep,nice,nicep,  &
     &                    qcwtr,qcwtrp,qcice,qcicep,qasl,qaslp,qt,qtp,  &
     &                    tke,tkep,ucpx,ucpy,vcpx,vcpy,wcpx,wcpy,       &
     &                    pcpx,pcpy,ptcpx,ptcpy,qvcpx,qvcpy,qwcpx,qwcpy,&
     &                    nwcpx,nwcpy,qicpx,qicpy,nicpx,nicpy,          &
     &                    qcwcpx,qcwcpy,qcicpx,qcicpy,qacpx,qacpy,      &
     &                    qtcpx,qtcpy,tkecpx,tkecpy,ugpv,utd,vgpv,vtd,  &
     &                    wgpv,wtd,ppgpv,pptd,ptpgpv,ptptd,qvgpv,qvtd,  &
     &                    qwgpv,qwtd,qigpv,qitd,qagpv,qatd,urdr,vrdr,   &
     &                    wrdr,qwrdr,qwrtd,qirdr,qirtd,sst,sstd,        &
     &                    maxvl,prwtr,price,qall,qallp,pdia,z0m,z0h,    &
     &                    tund,tundp,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,     &
     &                    tmp7,tmp8,tmp9,tmp10,tmp11,tmp12,tmp13,       &
     &                    tmp14,tmp15,tmp16,tmp17,tmp18,tmp19,          &
     &                    qwtmp,nwtmp,qitmp,nitmp,qcwtmp,qcitmp,qatmp,  &
     &                    qttmp,tketmp,tutmp,nid_rdr,njd_rdr,nkd_rdr,   &
     &                    km_rdr,lon_rdr,z_rdr,u_rdr,v_rdr,w_rdr,qp_rdr,&
     &                    tmp1_rdr)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: fpidate
                       ! Formal parameter of unique index of idate

      integer, intent(in) :: fpnggopt
                       ! Formal parameter of unique index of nggopt

      integer, intent(in) :: fpexbopt
                       ! Formal parameter of unique index of exbopt

      integer, intent(in) :: fplspopt
                       ! Formal parameter of unique index of lspopt

      integer, intent(in) :: fpvspopt
                       ! Formal parameter of unique index of vspopt

      integer, intent(in) :: fpngropt
                       ! Formal parameter of unique index of ngropt

      integer, intent(in) :: fpsfcopt
                       ! Formal parameter of unique index of sfcopt

      integer, intent(in) :: fpmasopt
                       ! Formal parameter of unique index of masopt

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fpaslopt
                       ! Formal parameter of unique index of aslopt

      integer, intent(in) :: fptubopt
                       ! Formal parameter of unique index of tubopt

      integer, intent(in) :: fpresopt
                       ! Formal parameter of unique index of resopt

      integer, intent(in) :: fpmxnopt
                       ! Formal parameter of unique index of mxnopt

      integer, intent(in) :: ksp0(1:2)
                       ! Index of lowest vertical sponge level

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

      integer, intent(in) :: nnw
                       ! Number of categories of water concentrations

      integer, intent(in) :: nqi
                       ! Number of categories of ice hydrometeor

      integer, intent(in) :: nni
                       ! Number of categories of ice concentrations

      integer, intent(in) :: km
                       ! Dimension of max(nk, nqw, nqi)

      integer, intent(in) :: nqa(0:4)
                       ! Number of types of aerosol

      integer, intent(in) :: nund
                       ! Number of soil and sea layers

      integer, intent(in) :: nid_rdr
                       ! Radar data dimension in x direction

      integer, intent(in) :: njd_rdr
                       ! Radar data dimension in y direction

      integer, intent(in) :: nkd_rdr
                       ! Radar data dimension in z direction

      integer, intent(in) :: km_rdr
                       ! Dimension of max(nk, nkd_rdr)

      integer, intent(in) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

      real, intent(in) :: area(0:4)
                       ! Area of each boundary plane

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(in) :: zsth(1:nk)
                       ! 1 dimensional stretched z coordinates

      real, intent(in) :: lat(0:ni+1,0:nj+1)
                       ! Latitude

      real, intent(in) :: lon(0:ni+1,0:nj+1)
                       ! Longitude

      real, intent(in) :: j31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of Jacobian

      real, intent(in) :: j32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of Jacobian

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: jcb8u(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at u points

      real, intent(in) :: jcb8v(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at v points

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: mf8u(0:ni+1,0:nj+1)
                       ! Map scale factors at u points

      real, intent(in) :: mf8v(0:ni+1,0:nj+1)
                       ! Map scale factors at v points

      real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

      real, intent(in) :: rmf8u(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at u points

      real, intent(in) :: rmf8v(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at v points

      real, intent(in) :: fc(0:ni+1,0:nj+1,1:2)
                       ! 0.25 x Coriolis parameters

      real, intent(in) :: ubr(0:ni+1,0:nj+1,1:nk)
                       ! Base state x components of velocity

      real, intent(in) :: vbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state y components of velocity

      real, intent(in) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: qvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state water vapor mixing ratio

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: rst8u(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at u points

      real, intent(in) :: rst8v(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at v points

      real, intent(in) :: rst8w(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at w points

      real, intent(in) :: rcsq(0:ni+1,0:nj+1,1:nk)
                       ! rbr x sound wave speed squared

      real, intent(in) :: rbcx(1:ni)
                       ! Relaxed lateral sponge damping coefficients
                       ! in x direction

      real, intent(in) :: rbcy(1:nj)
                       ! Relaxed lateral sponge damping coefficients
                       ! in y direction

      real, intent(in) :: rbcxy(1:ni,1:nj)
                       ! Relaxed lateral sponge damping coefficients

      real, intent(in) :: rbct(1:ni,1:nj,1:nk,1:2)
                       ! Relaxed top sponge damping coefficients

      real, intent(in) :: albe(0:ni+1,0:nj+1)
                       ! Albedo

      real, intent(in) :: beta(0:ni+1,0:nj+1)
                       ! Evapotranspiration efficiency

      real, intent(in) :: cap(0:ni+1,0:nj+1)
                       ! Thermal capacity

      real, intent(in) :: nuu(0:ni+1,0:nj+1)
                       ! Thermal diffusivity

      real, intent(in) :: kai(0:ni+1,0:nj+1)
                       ! Sea ice distribution

! Input and output variables

      real, intent(inout) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at present

      real, intent(inout) :: up(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at past

      real, intent(inout) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at present

      real, intent(inout) :: vp(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at past

      real, intent(inout) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at present

      real, intent(inout) :: wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

      real, intent(inout) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at present

      real, intent(inout) :: ppp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at past

      real, intent(inout) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at present

      real, intent(inout) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(inout) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at present

      real, intent(inout) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(inout) :: qwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at present

      real, intent(inout) :: qwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at past

      real, intent(inout) :: nwtr(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at present

      real, intent(inout) :: nwtrp(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at past

      real, intent(inout) :: qice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at present

      real, intent(inout) :: qicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at past

      real, intent(inout) :: nice(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at present

      real, intent(inout) :: nicep(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at past

      real, intent(inout) :: qcwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at present

      real, intent(inout) :: qcwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at past

      real, intent(inout) :: qcice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at present

      real, intent(inout) :: qcicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at past

      real, intent(inout) :: qasl(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at present

      real, intent(inout) :: qaslp(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at past

      real, intent(inout) :: qt(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at present

      real, intent(inout) :: qtp(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at past

      real, intent(inout) :: tke(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at present

      real, intent(inout) :: tkep(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at past

      real, intent(inout) :: ucpx(1:nj,1:nk,1:2)
                       ! Phase speed of x components of velocity
                       ! on west and east boundary

      real, intent(inout) :: ucpy(1:ni,1:nk,1:2)
                       ! Phase speed of x components of velocity
                       ! on south and north boundary

      real, intent(inout) :: vcpx(1:nj,1:nk,1:2)
                       ! Phase speed of y components of velocity
                       ! on west and east boundary

      real, intent(inout) :: vcpy(1:ni,1:nk,1:2)
                       ! Phase speed of y components of velocity
                       ! on south and north boundary

      real, intent(inout) :: wcpx(1:nj,1:nk,1:2)
                       ! Phase speed of z components of velocity
                       ! on west and east boundary

      real, intent(inout) :: wcpy(1:ni,1:nk,1:2)
                       ! Phase speed of z components of velocity
                       ! on south and north boundary

      real, intent(inout) :: pcpx(1:nj,1:nk,1:2)
                       ! Phase speed of pressure
                       ! on west and east boundary

      real, intent(inout) :: pcpy(1:ni,1:nk,1:2)
                       ! Phase speed of pressure
                       ! on south and north boundary

      real, intent(inout) :: ptcpx(1:nj,1:nk,1:2)
                       ! Phase speed of potential temperature
                       ! on west and east boundary

      real, intent(inout) :: ptcpy(1:ni,1:nk,1:2)
                       ! Phase speed of potential temperature
                       ! on south and north boundary

      real, intent(inout) :: qvcpx(1:nj,1:nk,1:2)
                       ! Phase speed of water vapor mixing ratio
                       ! on west and east boundary

      real, intent(inout) :: qvcpy(1:ni,1:nk,1:2)
                       ! Phase speed of water vapor mixing ratio
                       ! on south and north boundary

      real, intent(inout) :: qwcpx(1:nj,1:nk,1:2,1:nqw)
                       ! Phase speed of water hydrometeor
                       ! on west and east boundary

      real, intent(inout) :: qwcpy(1:ni,1:nk,1:2,1:nqw)
                       ! Phase speed of water hydrometeor
                       ! on south and north boundary

      real, intent(inout) :: nwcpx(1:nj,1:nk,1:2,1:nnw)
                       ! Phase speed of water concentrations
                       ! on west and east boundary

      real, intent(inout) :: nwcpy(1:ni,1:nk,1:2,1:nnw)
                       ! Phase speed of water concentrations
                       ! on south and north boundary

      real, intent(inout) :: qicpx(1:nj,1:nk,1:2,1:nqi)
                       ! Phase speed of ice hydrometeor
                       ! on west and east boundary

      real, intent(inout) :: qicpy(1:ni,1:nk,1:2,1:nqi)
                       ! Phase speed of ice hydrometeor
                       ! on south and north boundary

      real, intent(inout) :: nicpx(1:nj,1:nk,1:2,1:nni)
                       ! Phase speed of ice concentrations
                       ! on west and east boundary

      real, intent(inout) :: nicpy(1:ni,1:nk,1:2,1:nni)
                       ! Phase speed of ice concentrations
                       ! on south and north boundary

      real, intent(inout) :: qcwcpx(1:nj,1:nk,1:2,1:nqw)
                       ! Phase speed of charging distribution for water
                       ! on west and east boundary

      real, intent(inout) :: qcwcpy(1:ni,1:nk,1:2,1:nqw)
                       ! Phase speed of charging distribution for water
                       ! on south and north boundary

      real, intent(inout) :: qcicpx(1:nj,1:nk,1:2,1:nqi)
                       ! Phase speed of charging distribution for ice
                       ! on west and east boundary

      real, intent(inout) :: qcicpy(1:ni,1:nk,1:2,1:nqi)
                       ! Phase speed of charging distribution for ice
                       ! on south and north boundary

      real, intent(inout) :: qacpx(1:nj,1:nk,1:2,1:nqa(0))
                       ! Phase speed of aerosol mixing ratio
                       ! on west and east boundary

      real, intent(inout) :: qacpy(1:ni,1:nk,1:2,1:nqa(0))
                       ! Phase speed of aerosol mixing ratio
                       ! on south and north boundary

      real, intent(inout) :: qtcpx(1:nj,1:nk,1:2)
                       ! Phase speed of tracer mixing ratio
                       ! on west and east boundary

      real, intent(inout) :: qtcpy(1:ni,1:nk,1:2)
                       ! Phase speed of tracer mixing ratio
                       ! on south and north boundary

      real, intent(inout) :: tkecpx(1:nj,1:nk,1:2)
                       ! Phase speed of turbulent kinetic energy
                       ! on west and east boundary

      real, intent(inout) :: tkecpy(1:ni,1:nk,1:2)
                       ! Phase speed of turbulent kinetic energy
                       ! on south and north boundary

      real, intent(inout) :: ugpv(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity of GPV data
                       ! at marked time

      real, intent(inout) :: utd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! x components of velocity of GPV data

      real, intent(inout) :: vgpv(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity of GPV data
                       ! at marked time

      real, intent(inout) :: vtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! y components of velocity of GPV data

      real, intent(inout) :: wgpv(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity of GPV data
                       ! at marked time

      real, intent(inout) :: wtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! z components of velocity of GPV data

      real, intent(inout) :: ppgpv(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation of GPV data
                       ! at marked time

      real, intent(inout) :: pptd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! pressure perturbation of GPV data

      real, intent(inout) :: ptpgpv(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation of GPV data
                       ! at marked time

      real, intent(inout) :: ptptd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! potential temperature perturbation of GPV data

      real, intent(inout) :: qvgpv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio of GPV data
                       ! at marked time

      real, intent(inout) :: qvtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! water vapor mixing ratio of GPV data

      real, intent(inout) :: qwgpv(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor of GPV data at marked time

      real, intent(inout) :: qwtd(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Time tendency of water hydrometeor of GPV data

      real, intent(inout) :: qigpv(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor of GPV data at marked time

      real, intent(inout) :: qitd(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Time tendency of ice hydrometeor of GPV data

      real, intent(inout) :: qagpv(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Mixing ratio of aerosol data at marked time

      real, intent(inout) :: qatd(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Time tendency of mixing ratio of aerosol data

      real, intent(inout) :: urdr(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity of radar data
                       ! at marked time

      real, intent(inout) :: vrdr(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity of radar data
                       ! at marked time

      real, intent(inout) :: wrdr(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity of radar data
                       ! at marked time

      real, intent(inout) :: qwrdr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor of radar data
                       ! at marked time

      real, intent(inout) :: qwrtd(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Time tendency of
                       ! water hydrometeor of radar data

      real, intent(inout) :: qirdr(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor of radar data
                       ! at marked time

      real, intent(inout) :: qirtd(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Time tendency of
                       ! ice hydrometeor of radar data

      real, intent(inout) :: sst(0:ni+1,0:nj+1)
                       ! Sea surface temperature of external data
                       ! at marked time

      real, intent(inout) :: sstd(0:ni+1,0:nj+1)
                       ! Time tendency of
                       ! sea surface temperature of external data

      real, intent(inout) :: maxvl(0:ni+1,0:nj+1,1:nk)
                       ! Maximum instantaneous wind velocity

      real, intent(inout) :: prwtr(0:ni+1,0:nj+1,1:2,1:nqw)
                       ! Precipitation and accumulation for water

      real, intent(inout) :: price(0:ni+1,0:nj+1,1:2,1:nqi)
                       ! Precipitation and accumulation for ice

      real, intent(inout) :: qall(0:ni+1,0:nj+1,1:nk)
                       ! Total water and ice mixing ratio at present

      real, intent(inout) :: qallp(0:ni+1,0:nj+1,1:nk)
                       ! Total water and ice mixing ratio at past

      real, intent(inout) :: pdia(0:ni+1,0:nj+1,1:nk)
                       ! Diabatic value

      real, intent(inout) :: z0m(0:ni+1,0:nj+1)
                       ! Roughness length for velocity

      real, intent(inout) :: z0h(0:ni+1,0:nj+1)
                       ! Roughness length for scalar

      real, intent(inout) :: tund(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at present

      real, intent(inout) :: tundp(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at past

! Internal shared variables

      character(len=108) idate
                       ! Forecast start date
                       ! with Gregorian calendar, yyyymmddhhmm

      character(len=12) pdate
                       ! Forecast date at 1 step past
                       ! with Gregorian calendar, yyyymmddhhmm

      integer nggopt   ! Option for analysis nudging to GPV
      integer exbopt   ! Option for external boundary forcing
      integer lspopt   ! Option for lateral sponge damping
      integer vspopt   ! Option for vertical sponge damping
      integer ngropt   ! Option for analysis nudging to radar
      integer sfcopt   ! Option for surface physics
      integer masopt   ! Option for masscon model
      integer advopt   ! Option for advection scheme
      integer aslopt   ! Option for aerosol processes
      integer tubopt   ! Option for turbulent mixing
      integer resopt   ! Option for restart output
      integer mxnopt   ! Option for maxmum and minimum output

      integer(kind=i8) ibstp
                       ! Index of large time steps integration

      integer(kind=i8) nbstp0
                       ! Start index of large time steps

      integer(kind=i8) nbstp1
                       ! End index of large time steps

      integer nsstp    ! Number of small time steps

      integer nvstp    ! Number of steps
                       ! of vertical Cubic Lagrange advection

      integer nclstp(0:3)
                       ! Number of steps of cloud micro physics

      integer(kind=i8) ctime
                       ! Model current forecast time

      integer(kind=i8) ptime
                       ! Model forecast time at 1 step past

      integer(kind=i8) ftime
                       ! Model forecast time at 1 step future

      integer(kind=i8) pmin
                       ! 60000 x (ptime / 60000)

      integer fgpv     ! Descriptor to put into motion
                       ! for GPV nudging

      integer fasl     ! Descriptor to put into motion
                       ! for aerosol nudging

      integer frdr(1:2)
                       ! Descriptor to put into motion
                       ! for radar nudging

      integer fsst     ! Descriptor to put into motion
                       ! for sea surface temperature data reading

      real dtb         ! Large time steps interval
      real dts         ! Small time steps interval

      real dtsep       ! Time steps interval
                       ! of vertical Cubic Lagrange advection

      real dtsoil      ! Time interval of soil temperature calculation

      real dtcl(1:3)   ! Time interval of cloud micro physics

      real nggdmp      ! Analysis nudging damping coefficient for GPV

      real ngrdmp(1:2) ! Analysis nudging damping coefficient for radar

      real gtinc       ! Lapse of forecast time
                       ! from GPV data reading

      real atinc       ! Lapse of forecast time
                       ! from aerosol data reading

      real rtinc(1:2)  ! Lapse of forecast time
                       ! from radar data reading

      real stinc       ! Lapse of forecast time
                       ! from sea surface temperature data reading

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp4(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp5(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp6(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp7(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp8(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp9(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp10(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp11(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp12(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp13(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp14(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp15(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp16(0:ni+1,0:nj+1,1:km)
                       ! Temporary array

      real, intent(inout) :: tmp17(0:ni+1,0:nj+1,1:km)
                       ! Temporary array

      real, intent(inout) :: tmp18(0:ni+1,0:nj+1,1:km)
                       ! Temporary array

      real, intent(inout) :: tmp19(0:ni+1,0:nj+1,1:km)
                       ! Temporary array

      real, intent(inout) :: qwtmp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Temporary array

      real, intent(inout) :: nwtmp(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Temporary array

      real, intent(inout) :: qitmp(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Temporary array

      real, intent(inout) :: nitmp(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Temporary array

      real, intent(inout) :: qcwtmp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Temporary array

      real, intent(inout) :: qcitmp(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Temporary array

      real, intent(inout) :: qatmp(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Temporary array

      real, intent(inout) :: qttmp(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tketmp(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tutmp(0:ni+1,0:nj+1,1:nund)
                       ! Temporary array

      real, intent(inout) :: lon_rdr(1:nid_rdr,1:njd_rdr)
                       ! Longitude in radar data

      real, intent(inout) :: z_rdr(1:nid_rdr,1:njd_rdr,1:nkd_rdr)
                       ! z physical coordinates in radar data

      real, intent(inout) :: u_rdr(1:nid_rdr,1:njd_rdr,1:km_rdr)
                       ! x components of velocity in radar data

      real, intent(inout) :: v_rdr(1:nid_rdr,1:njd_rdr,1:km_rdr)
                       ! y components of velocity in radar data

      real, intent(inout) :: w_rdr(1:nid_rdr,1:njd_rdr,1:km_rdr)
                       ! z components of velocity in radar data

      real, intent(inout) :: qp_rdr(1:nid_rdr,1:njd_rdr,1:km_rdr)
                       ! Precipitation mixing ratio in radar data

      real, intent(inout) :: tmp1_rdr(1:nid_rdr,1:njd_rdr,1:nk)
                       ! Temporary array

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(idate)

! -----

! Get the required namelist variables.

      call getcname(fpidate,idate)
      call getiname(fpnggopt,nggopt)
      call getiname(fpexbopt,exbopt)
      call getiname(fplspopt,lspopt)
      call getiname(fpvspopt,vspopt)
      call getiname(fpngropt,ngropt)
      call getiname(fpsfcopt,sfcopt)
      call getiname(fpmasopt,masopt)
      call getiname(fpadvopt,advopt)
      call getiname(fpaslopt,aslopt)
      call getiname(fptubopt,tubopt)
      call getiname(fpresopt,resopt)
      call getiname(fpmxnopt,mxnopt)

! -----

! Read in the message to standard i/o.

      if(mype.eq.root) then

        call outstd05(0)

      end if

      call chkstd(root)

! -----

! Initialize the runtime variables.

      fgpv=0

      fasl=0

      frdr(1)=0
      frdr(2)=0

      fsst=0

      gtinc=0.e0

      atinc=0.e0

      rtinc(1)=0.e0
      rtinc(2)=0.e0

      stinc=0.e0

! -----

! Calculate the number of time integration steps.

      call slvstep(idiniopt,idadvopt,idcphopt,idstime,idetime,          &
     &             iddtbig,iddtsml,iddtvcul,iddtcmph,                   &
     &             nbstp0,nbstp1,nsstp,nvstp,nclstp)

! -----

!!!! The loop for the large time steps integration.

      do ibstp=nbstp0,nbstp1

! Get the current integration time and forecast date.

        call gettime(idadvopt,iddtbig,ibstp,ctime,ptime,ftime,pmin)

        call getdate(idate,pmin,pdate)

! -----

! Calculate the time interval.

        call timeint(idsfcopt,idadvopt,idcphopt,iddtbig,iddtsml,        &
     &               iddtvcul,iddtgrd,iddtcmph,ctime,                   &
     &               dtb,dts,dtsep,dtsoil,dtcl)

! -----

!!! Set the nudged variables.

! Read out the data from the interpolated GPV data file.

        if(nggopt.eq.1.or.exbopt.ge.1                                   &
     &    .or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

          if(ibstp-nbstp0.eq.0_i8) then
            fgpv=0
          else
            fgpv=1
          end if

          if(ctime.eq.0_i8) then

            call rdgpvnxt(idexprim,idcrsdir,idgpvvar,idncexp,idnccrs,   &
     &                    idwlngth,idgsmopt,idcphopt,idhaiopt,idetime,  &
     &                    idgpvitv,fgpv,ptime,ctime,gtinc,ni,nj,nk,     &
     &                    nqw,nqi,pbr,ptbr,ugpv,utd,vgpv,vtd,wgpv,wtd,  &
     &                    ppgpv,pptd,ptpgpv,ptptd,qvgpv,qvtd,           &
     &                    qwgpv,qwtd,qigpv,qitd,tmp1)

          else

            call rdgpvnxt(idexprim,idcrsdir,idgpvvar,idncexp,idnccrs,   &
     &                    idwlngth,idgsmopt,idcphopt,idhaiopt,idetime,  &
     &                    idgpvitv,fgpv,ctime,ftime,gtinc,ni,nj,nk,     &
     &                    nqw,nqi,pbr,ptbr,ugpv,utd,vgpv,vtd,wgpv,wtd,  &
     &                    ppgpv,pptd,ptpgpv,ptptd,qvgpv,qvtd,           &
     &                    qwgpv,qwtd,qigpv,qitd,tmp1)

          end if

        end if

! -----

! Read out the data from the interpolated aerosol data file.

        if(aslopt.ge.1) then

          if(nggopt.eq.1.or.exbopt.ge.1                                 &
     &      .or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

            if(ibstp-nbstp0.eq.0_i8) then
              fasl=0
            else
              fasl=1
            end if

            if(ctime.eq.0_i8) then

              call rdaslnxt(idexprim,idcrsdir,idncexp,idnccrs,idwlngth, &
     &                      idgsmopt,idetime,idaslitv,fasl,ptime,ctime, &
     &                      atinc,ni,nj,nk,nqa,qagpv,qatd,tmp1)

            else

              call rdaslnxt(idexprim,idcrsdir,idncexp,idnccrs,idwlngth, &
     &                      idgsmopt,idetime,idaslitv,fasl,ctime,ftime, &
     &                      atinc,ni,nj,nk,nqa,qagpv,qatd,tmp1)

            end if

          end if

        end if

! -----

!! Initialize the radar data nudging.

! Read out the data from the interpolated radar data file.

        if(ngropt.eq.1.or.ngropt.eq.2) then

          if(ngropt.eq.1) then

            if(ibstp-nbstp0.eq.0_i8) then
              frdr(1)=0
            else
              frdr(1)=1
            end if

            if(ctime.eq.0_i8) then

              call rdrdrnxt(idexprim,idcrsdir,idrdrvar,idncexp,idnccrs, &
     &                     idwlngth,idcphopt,idrdritv,idngrstr,idngrend,&
     &                     frdr,ptime,ctime,rtinc,ni,nj,nk,nqw,nqi,rbr, &
     &                     qwtrp,qicep,qwrdr,qwrtd,qirdr,qirtd,tmp1)

            else

              call rdrdrnxt(idexprim,idcrsdir,idrdrvar,idncexp,idnccrs, &
     &                     idwlngth,idcphopt,idrdritv,idngrstr,idngrend,&
     &                     frdr,ctime,ftime,rtinc,ni,nj,nk,nqw,nqi,rbr, &
     &                     qwtrp,qicep,qwrdr,qwrtd,qirdr,qirtd,tmp1)

            end if

          end if

          if(ibstp-nbstp0.eq.0_i8) then
            frdr(2)=0
          else
            frdr(2)=1
          end if

          if(ctime.eq.0_i8) then

            call rdrdrmrk(idexprim,idcrsdir,idrdrvar,idncexp,idnccrs,   &
     &                    idwlngth,idngropt,idcphopt,idrdritv,          &
     &                    idngrstr,idngrend,idngraff,frdr,ptime,ctime,  &
     &                    rtinc,ni,nj,nk,nqw,nqi,rbr,qwtrp,qicep,       &
     &                    urdr,vrdr,wrdr,qwrdr,qirdr,tmp1)

          else

            call rdrdrmrk(idexprim,idcrsdir,idrdrvar,idncexp,idnccrs,   &
     &                    idwlngth,idngropt,idcphopt,idrdritv,          &
     &                    idngrstr,idngrend,idngraff,frdr,ctime,ftime,  &
     &                    rtinc,ni,nj,nk,nqw,nqi,rbr,qwtrp,qicep,       &
     &                    urdr,vrdr,wrdr,qwrdr,qirdr,tmp1)

          end if

! -----

! Interpolate radar data.

        else if(ngropt.eq.12) then

          if(frdr(2).ge.0) then

            if(ibstp-nbstp0.eq.0_i8) then
              frdr(2)=0
            else
              frdr(2)=1
            end if

          end if

          if(ctime.eq.0_i8) then

            call s_setrdr(idrdrvar,iddatype_rdr,ididate,idcphopt,       &
     &                    idrotopt_rdr,idrdritv,idngrstr,idngrend,      &
     &                    idngraff,ptime,ctime,frdr,rtinc,ni,nj,nk,     &
     &                    nqw,nqi,zph,rbr,qwtrp,qicep,urdr,vrdr,wrdr,   &
     &                    qwrdr,qirdr,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,    &
     &                    tmp7,tmp8,nid_rdr,njd_rdr,nkd_rdr,km_rdr,     &
     &                    lon_rdr,z_rdr,u_rdr,v_rdr,w_rdr,qp_rdr,       &
     &                    tmp1_rdr)

          else

            call s_setrdr(idrdrvar,iddatype_rdr,ididate,idcphopt,       &
     &                    idrotopt_rdr,idrdritv,idngrstr,idngrend,      &
     &                    idngraff,ctime,ftime,frdr,rtinc,ni,nj,nk,     &
     &                    nqw,nqi,zph,rbr,qwtrp,qicep,urdr,vrdr,wrdr,   &
     &                    qwrdr,qirdr,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,    &
     &                    tmp7,tmp8,nid_rdr,njd_rdr,nkd_rdr,km_rdr,     &
     &                    lon_rdr,z_rdr,u_rdr,v_rdr,w_rdr,qp_rdr,       &
     &                    tmp1_rdr)

          end if

        end if

! -----

!! -----

! Get the control flag of analysis nudging.

        call ndgstep(idngrvar,idnggopt,idngropt,                        &
     &            idnggcoe,idnggdlt,idnggstr,idnggend,idnggc20,         &
     &            idngrcoe,idngrdlt,idngrstr,idngrend,idngrc20,idngraff,&
     &            ctime,frdr,rtinc,nggdmp,ngrdmp)

! -----

!!! -----

! Set and read out the external sea surface temperature.

        if(sfcopt.eq.3.or.sfcopt.eq.13) then

          if(ibstp-nbstp0.eq.0_i8) then
            fsst=0
          else
            fsst=1
          end if

          if(ctime.eq.0_i8) then

            call rdsstnxt(idexprim,idcrsdir,idncexp,idnccrs,idwlngth,   &
     &                    idetime,idsstitv,fsst,ptime,ctime,stinc,      &
     &                    ni,nj,sst,sstd)

          else

            call rdsstnxt(idexprim,idcrsdir,idncexp,idnccrs,idwlngth,   &
     &                    idetime,idsstitv,fsst,ctime,ftime,stinc,      &
     &                    ni,nj,sst,sstd)

          end if

        end if

! -----

!! Fit the x, y and z components of velocity to the mass consistent
!! equation and read in data to the dumped file at forecast start time.

        if(ctime.eq.0_i8) then

! Fit the x, y and z components of velocity to the mass consistent
! equation.

          if(masopt.eq.1) then

           call masscon(idadvopt,idsmtopt,idtubopt,iddxiv,iddyiv,iddziv,&
     &                  idmaseps,idalpha1,idalpha2,ctime,ni,nj,nk,      &
     &                  j31,j32,jcb8u,jcb8v,jcb8w,mf,u,up,v,vp,w,wp,    &
     &                  tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)

          end if

! -----

! Read in the data to the dumped file at forecast start time.

          call s_opendmp(idexprim,idcrsdir,idncexp,idnccrs,             &
     &                  idwlngth,iddmpfmt,iddmplev,iddmpmon,idetime,    &
     &                  iddmpitv,idmonitv,iddz,ctime,ni,nj,nk,zsth,tmp1)

          call outcheck(iddmpmon,ctime)

          if(sfcopt.ge.1) then

            call sfcphy(idadvopt,fmois,pdate,ctime,dtb,dtsoil,stinc,    &
     &                  ni,nj,nk,nqw,nqi,nund,zph,lat,lon,j31,j32,jcb8w,&
     &                  pbr,ptbr,rbr,rst,rst8u,rst8v,up,vp,wp,ppp,ptpp, &
     &                  qvp,prwtr,price,qallp,land,albe,beta,cap,nuu,   &
     &                  kai,sst,sstd,tmp1,tmp2,tmp3,tmp4,z0m,z0h,       &
     &                  tundp,tutmp,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,     &
     &                  tmp11,tmp12,tmp13,tmp14,tmp15,tmp16,tmp17,tmp18)

          end if

          call s_outdmp(iddmpvar,idcphopt,idhaiopt,idqcgopt,            &
     &                  idaslopt,idtrkopt,idtubopt,iddmplev,iddz,       &
     &                  fmois,ni,nj,nk,nqw,nnw,nqi,nni,nqa,zph,zsth,    &
     &                  lon,ubr,vbr,pbr,ptbr,qvbr,up,vp,wp,ppp,ptpp,    &
     &                  qvp,qwtrp,nwtrp,qicep,nicep,qcwtrp,qcicep,      &
     &                  qaslp,qtp,tkep,maxvl,prwtr,price,tmp1,tmp2,     &
     &                  tmp3,tmp4,tmp5,tmp6,tmp7)

          call closedmp(iddmpmon,ctime)

! -----

! Read in the maximum and minimum value of optional prognostic variable
! to the standard i/o at forecast start time.

          if(mxnopt.eq.1) then

            call mxndrv(idmxnvar,idcphopt,idhaiopt,idqcgopt,            &
     &                  idaslopt,idtrkopt,idtubopt,idetime,idmxnitv,    &
     &                  fmois,ctime,ni,nj,nk,nqw,nnw,nqi,nni,nqa,       &
     &                  up,vp,wp,ppp,ptpp,qvp,qwtrp,nwtrp,qicep,nicep,  &
     &                  qcwtrp,qcicep,qaslp,qtp,tkep,tmp1)

          end if

! -----

! Check the file to abort the solver.

          call abortslv(idcrsdir,idnccrs,iddmpmon,idresopt,idmxnopt,    &
     &                  iddmpitv,idmonitv,idresitv,idmxnitv,ctime)

! -----

! Read in the message to standard i/o.

          if(mype.eq.root) then

            call outstd05(0)

          end if

          call chkstd(root)

! -----

        end if

!! -----

! Open the dumped file when the current forecast time reaches marked
! time.

        call s_opendmp(idexprim,idcrsdir,idncexp,idnccrs,               &
     &                 idwlngth,iddmpfmt,iddmplev,iddmpmon,idetime,     &
     &                 iddmpitv,idmonitv,iddz,ftime,ni,nj,nk,zsth,tmp1)

        call outcheck(iddmpmon,ctime)

! -----

! Start the time integration driver routine.

        call intgdrv(idexbopt,idsfcopt,idadvopt,                        &
     &               idcphopt,idaslopt,idtubopt,iddiaopt,               &
     &               fmois,pdate,ksp0,nsstp,nvstp,nclstp,               &
     &               ctime,ftime,dtb,dts,dtsep,dtsoil,dtcl,             &
     &               nggdmp,ngrdmp,gtinc,atinc,rtinc,stinc,area,        &
     &               ni,nj,nk,nqw,nnw,nqi,nni,km,nqa,nund,zph,lat,lon,  &
     &               j31,j32,jcb,jcb8u,jcb8v,jcb8w,mf,mf8u,mf8v,        &
     &               rmf,rmf8u,rmf8v,fc,ubr,vbr,pbr,ptbr,qvbr,rbr,      &
     &               rst,rst8u,rst8v,rst8w,rcsq,rbcx,rbcy,rbcxy,        &
     &               rbct,ugpv,utd,vgpv,vtd,wgpv,wtd,ppgpv,pptd,        &
     &               ptpgpv,ptptd,qvgpv,qvtd,qwgpv,qwtd,qigpv,qitd,     &
     &               qagpv,qatd,urdr,vrdr,wrdr,qwrdr,qwrtd,qirdr,qirtd, &
     &               land,albe,beta,cap,nuu,kai,sst,sstd,u,up,v,vp,     &
     &               w,wp,pp,ppp,ptp,ptpp,qv,qvp,qwtr,qwtrp,nwtr,nwtrp, &
     &               qice,qicep,nice,nicep,qcwtr,qcwtrp,qcice,qcicep,   &
     &               qasl,qaslp,qt,qtp,tke,tkep,ucpx,ucpy,vcpx,vcpy,    &
     &               wcpx,wcpy,pcpx,pcpy,ptcpx,ptcpy,qvcpx,qvcpy,       &
     &               qwcpx,qwcpy,nwcpx,nwcpy,qicpx,qicpy,nicpx,nicpy,   &
     &               qcwcpx,qcwcpy,qcicpx,qcicpy,qacpx,qacpy,           &
     &               qtcpx,qtcpy,tkecpx,tkecpy,prwtr,price,qall,qallp,  &
     &               pdia,z0m,z0h,tund,tundp,tmp1,tmp2,tmp3,tmp4,tmp5,  &
     &               tmp6,tmp7,qwtmp,nwtmp,qitmp,nitmp,qcwtmp,qcitmp,   &
     &               qatmp,qttmp,tketmp,tutmp,tmp8,tmp9,tmp10,tmp11,    &
     &               tmp12,tmp13,tmp14,tmp15,tmp16,tmp17,tmp18,tmp19)

! -----

! Calculate the maximum instantaneous wind velocity.

        if(tubopt.ge.2) then

          if(advopt.le.3) then

            call instvel(iddmpvar,ni,nj,nk,u,v,w,tke,maxvl)

          else

            call instvel(iddmpvar,ni,nj,nk,up,vp,wp,tkep,maxvl)

          end if

        end if

! -----

! Read in the data to the dumped file when the current forecast time
! reaches marked time.

        if(advopt.le.3) then

          call s_outdmp(iddmpvar,idcphopt,idhaiopt,idqcgopt,            &
     &                  idaslopt,idtrkopt,idtubopt,iddmplev,iddz,       &
     &                  fmois,ni,nj,nk,nqw,nnw,nqi,nni,nqa,zph,zsth,    &
     &                  lon,ubr,vbr,pbr,ptbr,qvbr,u,v,w,pp,ptp,         &
     &                  qv,qwtr,nwtr,qice,nice,qcwtr,qcice,             &
     &                  qasl,qt,tke,maxvl,prwtr,price,tmp1,tmp2,        &
     &                  tmp3,tmp4,tmp5,tmp6,tmp7)

        else

          call s_outdmp(iddmpvar,idcphopt,idhaiopt,idqcgopt,            &
     &                  idaslopt,idtrkopt,idtubopt,iddmplev,iddz,       &
     &                  fmois,ni,nj,nk,nqw,nnw,nqi,nni,nqa,zph,zsth,    &
     &                  lon,ubr,vbr,pbr,ptbr,qvbr,up,vp,wp,ppp,ptpp,    &
     &                  qvp,qwtrp,nwtrp,qicep,nicep,qcwtrp,qcicep,      &
     &                  qaslp,qtp,tkep,maxvl,prwtr,price,tmp1,tmp2,     &
     &                  tmp3,tmp4,tmp5,tmp6,tmp7)

        end if

        call closedmp(iddmpmon,ftime)

! -----

! Read in the data to the restart file when the current forecast time
! reaches marked time.

        if(resopt.eq.1) then

          call outres(idexprim,idcrsdir,iddmpvar,idncexp,idnccrs,       &
     &                idwbc,idebc,idsbc,idnbc,idsfcopt,idadvopt,        &
     &                idcphopt,idqcgopt,idaslopt,idtrkopt,idtubopt,     &
     &                iddiaopt,idetime,idresitv,fmois,ftime,ni,nj,nk,   &
     &                nqw,nnw,nqi,nni,nqa,nund,ubr,vbr,pbr,ptbr,qvbr,   &
     &                u,up,v,vp,w,wp,pp,ppp,ptp,ptpp,qv,qvp,            &
     &                qwtr,qwtrp,nwtr,nwtrp,qice,qicep,nice,nicep,      &
     &                qcwtr,qcwtrp,qcice,qcicep,qasl,qaslp,qt,qtp,      &
     &                tke,tkep,ucpx,ucpy,vcpx,vcpy,wcpx,wcpy,pcpx,pcpy, &
     &                ptcpx,ptcpy,qvcpx,qvcpy,qwcpx,qwcpy,nwcpx,nwcpy,  &
     &                qicpx,qicpy,nicpx,nicpy,qcwcpx,qcwcpy,            &
     &                qcicpx,qcicpy,qacpx,qacpy,qtcpx,qtcpy,            &
     &                tkecpx,tkecpy,maxvl,prwtr,price,pdia,             &
     &                z0m,z0h,tund,tundp)

        end if

! -----

! Read in the maximum and minimum value of optional prognostic variable
! to the standard i/o when the current forecast time reaches marked
! time.

        if(mxnopt.eq.1) then

          if(advopt.le.3) then

            call mxndrv(idmxnvar,idcphopt,idhaiopt,idqcgopt,            &
     &                  idaslopt,idtrkopt,idtubopt,idetime,idmxnitv,    &
     &                  fmois,ftime,ni,nj,nk,nqw,nnw,nqi,nni,nqa,       &
     &                  u,v,w,pp,ptp,qv,qwtr,nwtr,qice,nice,            &
     &                  qcwtr,qcice,qasl,qt,tke,tmp1)

          else

            call mxndrv(idmxnvar,idcphopt,idhaiopt,idqcgopt,            &
     &                  idaslopt,idtrkopt,idtubopt,idetime,idmxnitv,    &
     &                  fmois,ftime,ni,nj,nk,nqw,nnw,nqi,nni,nqa,       &
     &                  up,vp,wp,ppp,ptpp,qvp,qwtrp,nwtrp,qicep,nicep,  &
     &                  qcwtrp,qcicep,qaslp,qtp,tkep,tmp1)

          end if

        end if

! -----

! Check the file to abort the solver.

        call abortslv(idcrsdir,idnccrs,iddmpmon,idresopt,idmxnopt,      &
     &                iddmpitv,idmonitv,idresitv,idmxnitv,ftime)

! -----

! Read in the message to standard i/o.

        if(mype.eq.root) then

          call outstd06(nbstp0,ibstp,ftime)

        end if

        call chkstd(root)

! -----

      end do

!!!! -----

      end subroutine s_slvdrv

!-----7--------------------------------------------------------------7--

      end module m_slvdrv
