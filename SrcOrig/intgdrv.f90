!***********************************************************************
      module m_intgdrv
!***********************************************************************

!     Author      : Sakakibara Atsushi, Naito Daisuke
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/05/20,
!                   1999/06/07, 1999/06/21, 1999/07/05, 1999/08/03,
!                   1999/08/18, 1999/08/23, 1999/09/16, 1999/09/30,
!                   1999/10/12, 1999/11/01, 1999/11/24, 1999/12/06,
!                   1999/12/20, 2000/01/05, 2000/01/17, 2000/02/02,
!                   2000/02/07, 2000/03/08, 2000/03/17, 2000/04/18,
!                   2000/06/01, 2000/07/05, 2000/12/19, 2001/01/15,
!                   2001/03/13, 2001/04/15, 2001/05/29, 2001/06/06,
!                   2001/07/13, 2001/08/07, 2001/09/13, 2001/10/18,
!                   2001/11/20, 2001/12/07, 2002/01/07, 2002/01/15,
!                   2002/04/02, 2002/06/06, 2002/07/03, 2002/07/23,
!                   2002/08/15, 2002/09/09, 2002/10/31, 2002/11/11,
!                   2002/12/02, 2003/01/04, 2003/01/20, 2003/02/13,
!                   2003/03/13, 2003/03/21, 2003/04/30, 2003/05/19,
!                   2003/07/15, 2003/08/01, 2003/10/31, 2003/11/05,
!                   2003/11/28, 2003/12/12, 2003/12/26, 2004/02/01,
!                   2004/03/05, 2004/03/22, 2004/04/01, 2004/04/15,
!                   2004/05/31, 2004/06/10, 2004/08/01, 2004/08/20,
!                   2004/09/25, 2004/10/12, 2004/10/12, 2004/12/17,
!                   2005/01/07, 2005/01/14, 2005/01/31, 2005/02/10,
!                   2005/04/04, 2005/08/05, 2005/10/05, 2006/01/10,
!                   2006/02/13, 2006/04/03, 2006/05/12, 2006/06/21,
!                   2006/07/21, 2006/09/21, 2006/09/30, 2006/11/06,
!                   2007/01/05, 2007/01/20, 2007/01/31, 2007/05/07,
!                   2007/05/21, 2007/07/30, 2007/11/26, 2008/01/11,
!                   2008/05/02, 2008/07/01, 2008/08/25, 2008/10/10,
!                   2008/12/11, 2009/01/30, 2009/02/27, 2009/03/23,
!                   2010/02/01, 2011/01/19, 2011/03/29, 2011/07/15,
!                   2011/08/18, 2011/09/22, 2011/11/10, 2012/06/19,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the inferior procedures for the time steps integration.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_adjstuv
      use m_aerophy
      use m_cloudphy
      use m_comindx
      use m_comkind
      use m_culintg
      use m_diabat
      use m_fdmdrv
      use m_getiname
      use m_phasev
      use m_phasevbc
      use m_resetag
      use m_satadjst
      use m_sfcphy
      use m_shift2nd
      use m_swp2nxt
      use m_timeflt
      use m_totalqwi

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: intgdrv, s_intgdrv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface intgdrv

        module procedure s_intgdrv

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_intgdrv(fpexbopt,fpsfcopt,fpadvopt,                  &
     &                fpcphopt,fpaslopt,fptubopt,fpdiaopt,              &
     &                fmois,pdate,ksp0,nsstp,nvstp,nclstp,              &
     &                ctime,ftime,dtb,dts,dtsep,dtsoil,dtcl,            &
     &                nggdmp,ngrdmp,gtinc,atinc,rtinc,stinc,area,       &
     &                ni,nj,nk,nqw,nnw,nqi,nni,km,nqa,nund,zph,lat,lon, &
     &                j31,j32,jcb,jcb8u,jcb8v,jcb8w,mf,mf8u,mf8v,       &
     &                rmf,rmf8u,rmf8v,fc,ubr,vbr,pbr,ptbr,qvbr,rbr,     &
     &                rst,rst8u,rst8v,rst8w,rcsq,rbcx,rbcy,rbcxy,       &
     &                rbct,ugpv,utd,vgpv,vtd,wgpv,wtd,ppgpv,pptd,       &
     &                ptpgpv,ptptd,qvgpv,qvtd,qwgpv,qwtd,qigpv,qitd,    &
     &                qagpv,qatd,urdr,vrdr,wrdr,qwrdr,qwrtd,qirdr,qirtd,&
     &                land,albe,beta,cap,nuu,kai,sst,sstd,u,up,v,vp,    &
     &                w,wp,pp,ppp,ptp,ptpp,qv,qvp,qwtr,qwtrp,nwtr,nwtrp,&
     &                qice,qicep,nice,nicep,qcwtr,qcwtrp,qcice,qcicep,  &
     &                qasl,qaslp,qt,qtp,tke,tkep,ucpx,ucpy,vcpx,vcpy,   &
     &                wcpx,wcpy,pcpx,pcpy,ptcpx,ptcpy,qvcpx,qvcpy,      &
     &                qwcpx,qwcpy,nwcpx,nwcpy,qicpx,qicpy,nicpx,nicpy,  &
     &                qcwcpx,qcwcpy,qcicpx,qcicpy,qacpx,qacpy,          &
     &                qtcpx,qtcpy,tkecpx,tkecpy,prwtr,price,qall,qallp, &
     &                pfrc,z0m,z0h,tund,tundp,uf,vf,wf,wc,ppf,ptpf,     &
     &                qvf,qwtrf,nwtrf,qicef,nicef,qcwtrf,qcicef,qaslf,  &
     &                qtf,tkef,tundf,ufrc,vfrc,wfrc,ptfrc,qvfrc,        &
     &                tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      character(len=12), intent(in) :: pdate
                       ! Forecast date at 1 step past
                       ! with Gregorian calendar, yyyymmddhhmm

      integer, intent(in) :: fpexbopt
                       ! Formal parameter of unique index of exbopt

      integer, intent(in) :: fpsfcopt
                       ! Formal parameter of unique index of sfcopt

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fpaslopt
                       ! Formal parameter of unique index of aslopt

      integer, intent(in) :: fptubopt
                       ! Formal parameter of unique index of tubopt

      integer, intent(in) :: fpdiaopt
                       ! Formal parameter of unique index of diaopt

      integer, intent(in) :: ksp0(1:2)
                       ! Index of lowest vertical sponge level

      integer, intent(in) :: nsstp
                       ! Number of small time steps

      integer, intent(in) :: nvstp
                       ! Number of steps
                       ! of vertical Cubic Lagrange advection

      integer, intent(in) :: nclstp(0:3)
                       ! Number of steps of cloud micro physics

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

      integer(kind=i8), intent(in) :: ftime
                       ! Model forecast time at 1 step future

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

      integer, intent(in) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: dts
                       ! Small time steps interval

      real, intent(in) :: dtsep
                       ! Time steps interval
                       ! of vertical Cubic Lagrange advection

      real, intent(in) :: dtsoil
                       ! Time interval of soil temperature calculation

      real, intent(in) :: dtcl(1:3)
                       ! Time interval of cloud micro physics

      real, intent(in) :: nggdmp
                       ! Analysis nudging damping coefficient for GPV

      real, intent(in) :: ngrdmp(1:2)
                       ! Analysis nudging damping coefficient for radar

      real, intent(in) :: gtinc
                       ! Lapse of forecast time
                       ! from GPV data reading

      real, intent(in) :: atinc
                       ! Lapse of forecast time
                       ! from aerosol data reading

      real, intent(in) :: rtinc(1:2)
                       ! Lapse of forecast time
                       ! from radar data reading

      real, intent(in) :: stinc
                       ! Lapse of forecast time
                       ! from sea surface temperature data reading

      real, intent(in) :: area(0:4)
                       ! Area of each boundary plane

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

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

      real, intent(in) :: ugpv(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: utd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! x components of velocity of GPV data

      real, intent(in) :: vgpv(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: vtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! y components of velocity of GPV data

      real, intent(in) :: wgpv(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: wtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! z components of velocity of GPV data

      real, intent(in) :: ppgpv(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation of GPV data
                       ! at marked time

      real, intent(in) :: pptd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! pressure perturbation of GPV data

      real, intent(in) :: ptpgpv(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation of GPV data
                       ! at marked time

      real, intent(in) :: ptptd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! potential temperature perturbation of GPV data

      real, intent(in) :: qvgpv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio of GPV data
                       ! at marked time

      real, intent(in) :: qvtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! water vapor mixing ratio of GPV data

      real, intent(in) :: qwgpv(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor of GPV data at marked time

      real, intent(in) :: qwtd(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Time tendency of water hydrometeor of GPV data

      real, intent(in) :: qigpv(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor of GPV data at marked time

      real, intent(in) :: qitd(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Time tendency of ice hydrometeor of GPV data

      real, intent(in) :: qagpv(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Mixing ratio of aerosol data at marked time

      real, intent(in) :: qatd(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Time tendency of mixing ratio of aerosol data

      real, intent(in) :: urdr(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity of radar data
                       ! at marked time

      real, intent(in) :: vrdr(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity of radar data
                       ! at marked time

      real, intent(in) :: wrdr(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity of radar data
                       ! at marked time

      real, intent(in) :: qwrdr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor of radar data
                       ! at marked time

      real, intent(in) :: qwrtd(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Time tendency of
                       ! water hydrometeor of radar data

      real, intent(in) :: qirdr(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor of radar data
                       ! at marked time

      real, intent(in) :: qirtd(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Time tendency of
                       ! ice hydrometeor of radar data

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

      real, intent(in) :: sst(0:ni+1,0:nj+1)
                       ! Sea surface temperature of external data
                       ! at marked time

      real, intent(in) :: sstd(0:ni+1,0:nj+1)
                       ! Time tendency of
                       ! sea surface temperature of external data

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

      real, intent(inout) :: prwtr(0:ni+1,0:nj+1,1:2,1:nqw)
                       ! Precipitation and accumulation for water

      real, intent(inout) :: price(0:ni+1,0:nj+1,1:2,1:nqi)
                       ! Precipitation and accumulation for ice

      real, intent(inout) :: qall(0:ni+1,0:nj+1,1:nk)
                       ! Total water and ice mixing ratio at present

      real, intent(inout) :: qallp(0:ni+1,0:nj+1,1:nk)
                       ! Total water and ice mixing ratio at past

      real, intent(inout) :: pfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in pressure equation

      real, intent(inout) :: z0m(0:ni+1,0:nj+1)
                       ! Roughness length for velocity

      real, intent(inout) :: z0h(0:ni+1,0:nj+1)
                       ! Roughness length for scalar

      real, intent(inout) :: tund(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at present

      real, intent(inout) :: tundp(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at past

! Internal shared variables

      integer exbopt   ! Option for external boundary forcing
      integer sfcopt   ! Option for surface physics
      integer advopt   ! Option for advection scheme
      integer cphopt   ! Option for cloud micro physics
      integer aslopt   ! Option for aerosol processes
      integer tubopt   ! Option for turbulent mixing
      integer diaopt   ! Option for diabatic calculation

      real, intent(inout) :: uf(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at future

      real, intent(inout) :: vf(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at future

      real, intent(inout) :: wf(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at future

      real, intent(inout) :: wc(0:ni+1,0:nj+1,1:nk)
                       ! zeta components of contravariant velocity

      real, intent(inout) :: ppf(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at future

      real, intent(inout) :: ptpf(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at future

      real, intent(inout) :: qvf(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at future

      real, intent(inout) :: qwtrf(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at future

      real, intent(inout) :: nwtrf(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at future

      real, intent(inout) :: qicef(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at future

      real, intent(inout) :: nicef(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at future

      real, intent(inout) :: qcwtrf(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at future

      real, intent(inout) :: qcicef(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at future

      real, intent(inout) :: qaslf(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at future

      real, intent(inout) :: qtf(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at future

      real, intent(inout) :: tkef(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at future

      real, intent(inout) :: tundf(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at future

      real, intent(inout) :: ufrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in u equation

      real, intent(inout) :: vfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in v equation

      real, intent(inout) :: wfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in w equation

      real, intent(inout) :: ptfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in potential temperature equation

      real, intent(inout) :: qvfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in water vapor mixing ratio equation

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp4(0:ni+1,0:nj+1,1:km)
                       ! Temporary array

      real, intent(inout) :: tmp5(0:ni+1,0:nj+1,1:km)
                       ! Temporary array

      real, intent(inout) :: tmp6(0:ni+1,0:nj+1,1:km)
                       ! Temporary array

      real, intent(inout) :: tmp7(0:ni+1,0:nj+1,1:km)
                       ! Temporary array

! Remark

!     u,v,w,pp,ptp,qv,pfrc: These variables are also temporary.

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpexbopt,exbopt)
      call getiname(fpsfcopt,sfcopt)
      call getiname(fpadvopt,advopt)
      call getiname(fpcphopt,cphopt)
      call getiname(fpaslopt,aslopt)
      call getiname(fptubopt,tubopt)
      call getiname(fpdiaopt,diaopt)

! -----

! Get the total water and ice mixing ratio at past.

      if(fmois(1:5).eq.'moist') then

        if(advopt.le.3) then

          call totalqwi(idcphopt,idhaiopt,ni,nj,nk,nqw,nqi,qwtr,qice,   &
     &                  qall)

        end if

        call totalqwi(idcphopt,idhaiopt,ni,nj,nk,nqw,nqi,qwtrp,qicep,   &
     &                qallp)

      end if

! -----

! Perform the surface physics.

      if(sfcopt.ge.1.and.tubopt.ge.1) then

        call sfcphy(idadvopt,fmois,pdate,ftime,dtb,dtsoil,stinc,        &
     &              ni,nj,nk,nqw,nqi,nund,zph,lat,lon,j31,j32,jcb8w,    &
     &              pbr,ptbr,rbr,rst,rst8u,rst8v,up,vp,wp,ppp,ptpp,qvp, &
     &              prwtr,price,qallp,land,albe,beta,cap,nuu,kai,       &
     &              sst,sstd,uf,vf,ptpf,qvf,z0m,z0h,tundp,tundf,        &
     &              ufrc,vfrc,ptfrc,qvfrc,tmp1,tmp2,tmp3,tmp4,          &
     &              tmp5,tmp6,tmp7,wf,wc,wfrc)

      end if

! -----

! Perform the time integration.

      call fdmdrv(fmois,ksp0,nsstp,ctime,                               &
     &            dtb,dts,nggdmp,ngrdmp,gtinc,atinc,rtinc,              &
     &            ni,nj,nk,nqw,nnw,nqi,nni,nqa,zph,                     &
     &            j31,j32,jcb,jcb8u,jcb8v,jcb8w,mf,mf8u,mf8v,           &
     &            rmf,rmf8u,rmf8v,fc,ubr,vbr,pbr,ptbr,qvbr,rbr,         &
     &            rst,rst8u,rst8v,rst8w,rcsq,u,up,v,vp,w,wp,            &
     &            pp,ppp,ptp,ptpp,qv,qvp,qwtr,qwtrp,nwtr,nwtrp,         &
     &            qice,qicep,nice,nicep,qcwtr,qcwtrp,qcice,qcicep,      &
     &            qasl,qaslp,qt,qtp,ucpx,ucpy,vcpx,vcpy,                &
     &            wcpx,wcpy,pcpx,pcpy,ptcpx,ptcpy,qvcpx,qvcpy,          &
     &            qwcpx,qwcpy,nwcpx,nwcpy,qicpx,qicpy,nicpx,nicpy,      &
     &            qcwcpx,qcwcpy,qcicpx,qcicpy,qacpx,qacpy,              &
     &            qtcpx,qtcpy,tkecpx,tkecpy,rbcx,rbcy,rbcxy,rbct,       &
     &            ugpv,utd,vgpv,vtd,wgpv,wtd,ppgpv,pptd,                &
     &            ptpgpv,ptptd,qvgpv,qvtd,qwgpv,qwtd,qigpv,qitd,        &
     &            qagpv,qatd,urdr,vrdr,wrdr,qwrdr,qwrtd,qirdr,qirtd,    &
     &            qall,qallp,tke,tkep,ufrc,vfrc,pfrc,ptfrc,qvfrc,       &
     &            uf,vf,wf,wc,ppf,ptpf,qvf,qwtrf,nwtrf,qicef,nicef,     &
     &            qcwtrf,qcicef,qaslf,qtf,tkef,wfrc,tmp1,tmp2,tmp3,     &
     &            tmp4,tmp5,tmp6,tmp7)

! -----

! Perform the Cubic Lagrange scheme.

      if(advopt.ge.4) then

       call shift2nd(idwbc,idebc,idsbc,idnbc,idadvopt,idsmtopt,idcphopt,&
     &               idhaiopt,idqcgopt,idaslopt,idtrkopt,idtubopt,fmois,&
     &               'oooooooooo',ni,nj,nk,nqw,nnw,nqi,nni,nqa,uf,vf,wf,&
     &               ppf,ptpf,qvf,qwtrf,nwtrf,qicef,nicef,qcwtrf,qcicef,&
     &               qaslf,qtf,tkef)

       call culintg(fmois,nvstp,dtb,dtsep,ni,nj,nk,nqw,nnw,nqi,nni,nqa, &
     &              j31,j32,jcb8w,mf,mf8u,mf8v,uf,vf,wf,ppf,ptpf,qvf,   &
     &              qwtrf,nwtrf,qicef,nicef,qcwtrf,qcicef,qaslf,qtf,    &
     &              tkef,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,            &
     &              ufrc,vfrc,wfrc)

      end if

! -----

! Adjust the x and y components of velocity.

      if(exbopt.ge.11) then

        call adjstuv(idwbc,idebc,idadvopt,idmpopt,idmfcopt,             &
     &               iddx,iddy,iddz,dtb,gtinc,area,ni,nj,nk,            &
     &               rmf,rmf8u,rmf8v,rst8u,rst8v,ppp,ppf,               &
     &               ugpv,utd,vgpv,vtd,uf,vf)

      end if

! -----

! Perform the surface physics.

      if(sfcopt.ge.1.and.tubopt.eq.0) then

        call sfcphy(idadvopt,fmois,pdate,ftime,dtb,dtsoil,stinc,        &
     &              ni,nj,nk,nqw,nqi,nund,zph,lat,lon,j31,j32,jcb8w,    &
     &              pbr,ptbr,rbr,rst,rst8u,rst8v,up,vp,wp,ppp,ptpp,qvp, &
     &              prwtr,price,qallp,land,albe,beta,cap,nuu,kai,       &
     &              sst,sstd,uf,vf,ptpf,qvf,z0m,z0h,tundp,tundf,        &
     &              ufrc,vfrc,ptfrc,qvfrc,tmp1,tmp2,tmp3,tmp4,          &
     &              tmp5,tmp6,tmp7,wc,wfrc,pfrc)

      end if

! -----

! Perform the cloud physics.

      if(fmois(1:5).eq.'moist') then

        call cloudphy(idadvopt,idcphopt,idhaiopt,nclstp,dtb,dtcl,       &
     &               ni,nj,nk,nqw,nnw,nqi,nni,km,jcb,jcb8w,pbr,ptbr,rbr,&
     &               rst,wf,ppp,ptpp,qvp,qwtrp,qicep,qallp,ptpf,qvf,    &
     &               qwtrf,nwtrf,nwtrp,qicef,nicef,nicep,qcwtrf,qcicef, &
     &               prwtr,price,ufrc,vfrc,wfrc,pfrc,ptfrc,qvfrc,wc,    &
     &               tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7)

      end if

! -----

! Perform the aerosol processes.

      if(aslopt.ge.1) then

        call aerophy(idadvopt,idcphopt,dtb,ni,nj,nk,                    &
     &               nqw,nnw,nqi,nni,nqa,jcb,pbr,ptbr,rbr,rst,          &
     &               ppp,ptpp,qvp,qwtrp,nwtrp,qicep,nicep,qaslp,        &
     &               qaslf,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7)

      end if

! -----

! Perform the Asselin time filter.

      if(advopt.le.3.and.abs(cphopt).lt.20.and.ctime.ne.0_i8) then

        call timeflt(idsfcopt,idcphopt,idhaiopt,                        &
     &               idqcgopt,idaslopt,idtrkopt,idtubopt,idfilcoe,      &
     &               fmois,dtsoil,ni,nj,nk,nqw,nnw,nqi,nni,nqa,nund,    &
     &               land,up,uf,vp,vf,wp,wf,ppp,ppf,ptpp,ptpf,qvp,qvf,  &
     &               qwtrp,qwtrf,nwtrp,nwtrf,qicep,qicef,nicep,nicef,   &
     &               qcwtrp,qcwtrf,qcicep,qcicef,qaslp,qaslf,           &
     &               qtp,qtf,tkep,tkef,tundp,tundf,u,v,w,pp,ptp,qv,     &
     &               qwtr,nwtr,qice,nice,qcwtr,qcice,qasl,qt,tke,tund)

      end if

! -----

! Perform the saturation adjustment.

      if(fmois(1:5).eq.'moist') then

        if(advopt.le.3) then

          call satadjst(idcphopt,ni,nj,nk,nqw,nnw,nqi,nni,pbr,ptbr,     &
     &                  w,pp,ptp,qv,qwtr,nwtr,qice,nice,tmp1,tmp2)

        end if

        call satadjst(idcphopt,ni,nj,nk,nqw,nnw,nqi,nni,pbr,ptbr,       &
     &                wf,ppf,ptpf,qvf,qwtrf,nwtrf,qicef,nicef,tmp1,tmp2)

      end if

! -----

! Exchange the value in the case the 4th order calculation is performed.

      if(advopt.le.3) then

       call shift2nd(idwbc,idebc,idsbc,idnbc,idadvopt,idsmtopt,idcphopt,&
     &               idhaiopt,idqcgopt,idaslopt,idtrkopt,idtubopt,fmois,&
     &               'oooooooooo',ni,nj,nk,nqw,nnw,nqi,nni,nqa,u,v,w,   &
     &               pp,ptp,qv,qwtr,nwtr,qice,nice,qcwtr,qcice,         &
     &               qasl,qt,tke)

      end if

      call shift2nd(idwbc,idebc,idsbc,idnbc,idadvopt,idsmtopt,idcphopt, &
     &              idhaiopt,idqcgopt,idaslopt,idtrkopt,idtubopt,fmois, &
     &              'oooooooooo',ni,nj,nk,nqw,nnw,nqi,nni,nqa,uf,vf,wf, &
     &              ppf,ptpf,qvf,qwtrf,nwtrf,qicef,nicef,qcwtrf,qcicef, &
     &              qaslf,qtf,tkef)

! -----

! Calculate the phase speed for the open boundary conditions.

      if(advopt.le.3) then

        call phasev(idgwmopt,idcphopt,idhaiopt,idqcgopt,idaslopt,       &
     &           idtrkopt,idtubopt,fmois,dtb,dts,dtsep,ni,nj,nk,        &
     &           nqw,nnw,nqi,nni,nqa,rmf,rmf8u,rmf8v,u,up,uf,v,vp,vf,   &
     &           w,wp,wf,pp,ppp,ppf,ptp,ptpp,ptpf,qv,qvp,qvf,           &
     &           qwtr,qwtrp,qwtrf,nwtr,nwtrp,nwtrf,qice,qicep,qicef,    &
     &           nice,nicep,nicef,qcwtr,qcwtrp,qcwtrf,                  &
     &           qcice,qcicep,qcicef,qasl,qaslp,qaslf,qt,qtp,qtf,       &
     &           tke,tkep,tkef,ucpx,ucpy,vcpx,vcpy,wcpx,wcpy,pcpx,pcpy, &
     &           ptcpx,ptcpy,qvcpx,qvcpy,qwcpx,qwcpy,nwcpx,nwcpy,       &
     &           qicpx,qicpy,nicpx,nicpy,qcwcpx,qcwcpy,qcicpx,qcicpy,   &
     &           qacpx,qacpy,qtcpx,qtcpy,tkecpx,tkecpy,                 &
     &           tmp1,tmp2,tmp3,tmp4)

      else

        call phasev(idgwmopt,idcphopt,idhaiopt,idqcgopt,idaslopt,       &
     &           idtrkopt,idtubopt,fmois,dtb,dts,dtsep,ni,nj,nk,        &
     &           nqw,nnw,nqi,nni,nqa,rmf,rmf8u,rmf8v,up,up,uf,vp,vp,vf, &
     &           wp,wp,wf,ppp,ppp,ppf,ptpp,ptpp,ptpf,qvp,qvp,qvf,       &
     &           qwtrp,qwtrp,qwtrf,nwtrp,nwtrp,nwtrf,qicep,qicep,qicef, &
     &           nicep,nicep,nicef,qcwtrp,qcwtrp,qcwtrf,                &
     &           qcicep,qcicep,qcicef,qaslp,qaslp,qaslf,qtp,qtp,qtf,    &
     &           tkep,tkep,tkef,ucpx,ucpy,vcpx,vcpy,wcpx,wcpy,pcpx,pcpy,&
     &           ptcpx,ptcpy,qvcpx,qvcpy,qwcpx,qwcpy,nwcpx,nwcpy,       &
     &           qicpx,qicpy,nicpx,nicpy,qcwcpx,qcwcpy,qcicpx,qcicpy,   &
     &           qacpx,qacpy,qtcpx,qtcpy,tkecpx,tkecpy,                 &
     &           tmp1,tmp2,tmp3,tmp4)

      end if

! -----

! Calculate the diffrential phase speed term between the external
! boundary and model grid variables.

      if(exbopt.ge.1) then

        if(advopt.le.3) then

          call phasevbc(idexbvar,idgwmopt,idcphopt,idhaiopt,idaslopt,   &
     &                  fmois,dtb,dts,dtsep,gtinc,ni,nj,nk,nqw,nqi,nqa, &
     &                  rmf,rmf8u,rmf8v,u,up,uf,v,vp,vf,w,wp,wf,        &
     &                  pp,ppp,ppf,ptp,ptpp,ptpf,qv,qvp,qvf,            &
     &                  qwtr,qwtrp,qwtrf,qice,qicep,qicef,              &
     &                  qasl,qaslp,qaslf,ugpv,utd,vgpv,vtd,wgpv,wtd,    &
     &                  ppgpv,pptd,ptpgpv,ptptd,qvgpv,qvtd,qwgpv,qwtd,  &
     &                  qigpv,qitd,qagpv,qatd,ucpx,ucpy,vcpx,vcpy,      &
     &                  wcpx,wcpy,pcpx,pcpy,ptcpx,ptcpy,qvcpx,qvcpy,    &
     &                  qwcpx,qwcpy,qicpx,qicpy,qacpx,qacpy,            &
     &                  tmp1,tmp2,tmp3,tmp4)

        else

          call phasevbc(idexbvar,idgwmopt,idcphopt,idhaiopt,idaslopt,   &
     &                  fmois,dtb,dts,dtsep,gtinc,ni,nj,nk,nqw,nqi,nqa, &
     &                  rmf,rmf8u,rmf8v,up,up,uf,vp,vp,vf,wp,wp,wf,     &
     &                  ppp,ppp,ppf,ptpp,ptpp,ptpf,qvp,qvp,qvf,         &
     &                  qwtrp,qwtrp,qwtrf,qicep,qicep,qicef,            &
     &                  qaslp,qaslp,qaslf,ugpv,utd,vgpv,vtd,wgpv,wtd,   &
     &                  ppgpv,pptd,ptpgpv,ptptd,qvgpv,qvtd,qwgpv,qwtd,  &
     &                  qigpv,qitd,qagpv,qatd,ucpx,ucpy,vcpx,vcpy,      &
     &                  wcpx,wcpy,pcpx,pcpy,ptcpx,ptcpy,qvcpx,qvcpy,    &
     &                  qwcpx,qwcpy,qicpx,qicpy,qacpx,qacpy,            &
     &                  tmp1,tmp2,tmp3,tmp4)

        end if

      end if

! -----

! Calculate the diabatic to the next large time step.

      if(advopt.le.3.and.diaopt.eq.1) then

        if(fmois(1:5).eq.'moist') then

          call totalqwi(idcphopt,idhaiopt,ni,nj,nk,nqw,nqi,qwtrf,qicef, &
     &                  tmp1)

        end if

        call diabat(idsmtopt,idcphopt,idiwest,idieast,idjsouth,idjnorth,&
     &              fmois,dtb,ni,nj,nk,j31,j32,jcb,jcb8u,jcb8v,jcb8w,mf,&
     &              mf8u,mf8v,ptbr,rcsq,u,v,w,ptp,ptpp,ptpf,qv,qvp,qvf, &
     &              qall,qallp,tmp1,pfrc,wc,tmp2,tmp3,tmp4,             &
     &              ufrc,vfrc,wfrc,ptfrc,qvfrc)

      end if

! -----

! Swap the prognostic variables to the next time step.

      call swp2nxt(idsfcopt,idadvopt,idcphopt,                          &
     &             idhaiopt,idqcgopt,idaslopt,idtrkopt,idtubopt,        &
     &             idiwest,idieast,idjsouth,idjnorth,fmois,dtsoil,      &
     &             ni,nj,nk,nqw,nnw,nqi,nni,nqa,nund,uf,vf,wf,ppf,ptpf, &
     &             qvf,qwtrf,nwtrf,qicef,nicef,qcwtrf,qcicef,qaslf,qtf, &
     &             tkef,tundf,u,v,w,pp,ptp,qv,qwtr,nwtr,qice,nice,      &
     &             qcwtr,qcice,qasl,qt,tke,tund,up,vp,wp,ppp,ptpp,      &
     &             qvp,qwtrp,nwtrp,qicep,nicep,qcwtrp,qcicep,           &
     &             qaslp,qtp,tkep,tundp)

! -----

! Reset the message tag.

      call resetag

! -----

      end subroutine s_intgdrv

!-----7--------------------------------------------------------------7--

      end module m_intgdrv
