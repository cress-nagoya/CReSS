!***********************************************************************
      module m_inivar
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/01/25, 1999/03/25, 1999/04/06, 1999/05/10,
!                   1999/05/20, 1999/06/07, 1999/06/21, 1999/06/28,
!                   1999/07/05, 1999/07/23, 1999/08/03, 1999/08/09,
!                   1999/08/23, 1999/09/06, 1999/09/30, 1999/10/12,
!                   1999/11/01, 1999/11/19, 1999/12/20, 2000/01/05,
!                   2000/01/17, 2000/02/02, 2000/02/07, 2000/03/08,
!                   2000/04/18, 2000/06/01, 2000/12/18, 2001/01/15,
!                   2001/03/13, 2001/06/06, 2001/06/29, 2001/07/13,
!                   2001/08/07, 2001/09/13, 2001/10/18, 2001/11/14,
!                   2002/01/15, 2002/02/05, 2002/04/02, 2002/06/18,
!                   2002/07/03, 2002/07/15, 2002/07/23, 2002/08/15,
!                   2002/08/27, 2002/10/31, 2003/01/04, 2003/03/21,
!                   2003/05/19, 2003/07/15, 2003/08/08, 2003/10/10,
!                   2003/11/05, 2003/12/12, 2004/01/09, 2004/02/01,
!                   2004/03/05, 2004/04/01, 2004/04/15, 2004/05/31,
!                   2004/06/10, 2004/07/01, 2004/07/10, 2004/08/01,
!                   2004/08/31, 2004/09/01, 2004/09/10, 2004/09/25,
!                   2005/01/14, 2005/02/10, 2005/08/05, 2006/01/10,
!                   2006/02/13, 2006/04/03, 2006/07/21, 2006/09/30,
!                   2006/11/06, 2006/11/27, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/04/11, 2007/05/07, 2007/05/14,
!                   2007/05/21, 2007/07/30, 2007/08/24, 2008/01/11,
!                   2008/05/02, 2008/07/01, 2008/08/25, 2008/10/10,
!                   2008/12/11, 2009/01/30, 2009/02/27, 2011/08/18,
!                   2011/09/22, 2011/11/10, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the initial conditions.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_baserho
      use m_castvar
      use m_chkerr
      use m_chkstd
      use m_comindx
      use m_commpi
      use m_destroy
      use m_getiname
      use m_inisfc
      use m_initund
      use m_intrpsnd
      use m_lspdmp
      use m_move2d
      use m_outgeo
      use m_rdaslini
      use m_rdgpvini
      use m_rdres
      use m_rdsnd
      use m_set1d
      use m_setbase
      use m_setbin
      use m_setvar1d
      use m_setvar3d
      use m_sndwave
      use m_vintsnd
      use m_vspdmp

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: inivar, s_inivar

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface inivar

        module procedure s_inivar

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_inivar(fpiniopt,fpmovopt,fplspopt,fpvspopt,          &
     &                    fpsfcopt,fpcphopt,fpaslopt,fmois,ksp0,        &
     &                    ni,nj,nk,nqw,nnw,nqi,nni,nqa,nund,zph,        &
     &                    lat,lon,j31,j32,jcb,jcb8w,mf,fc,ubr,vbr,      &
     &                    pbr,ptbr,qvbr,rbr,rst,rst8u,rst8v,rst8w,      &
     &                    rcsq,u,up,v,vp,w,wp,pp,ppp,ptp,ptpp,qv,qvp,   &
     &                    qwtr,qwtrp,nwtr,nwtrp,qice,qicep,nice,nicep,  &
     &                    qcwtr,qcwtrp,qcice,qcicep,qasl,qaslp,         &
     &                    qt,qtp,tke,tkep,rbcx,rbcy,rbcxy,rbct,         &
     &                    ucpx,ucpy,vcpx,vcpy,wcpx,wcpy,pcpx,pcpy,      &
     &                    ptcpx,ptcpy,qvcpx,qvcpy,qwcpx,qwcpy,          &
     &                    nwcpx,nwcpy,qicpx,qicpy,nicpx,nicpy,          &
     &                    qcwcpx,qcwcpy,qcicpx,qcicpy,qacpx,qacpy,      &
     &                    qtcpx,qtcpy,tkecpx,tkecpy,maxvl,prwtr,price,  &
     &                    pdia,land,albe,beta,z0m,z0h,cap,nuu,kai,      &
     &                    tund,tundp,tmp1,tmp2,tmp3,nlev,z1d,u1d,v1d,   &
     &                    p1d,pt1d,qv1d,ltmp1,ltmp2,ltmp3)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpiniopt
                       ! Formal parameter of unique index of iniopt

      integer, intent(in) :: fpmovopt
                       ! Formal parameter of unique index of movopt

      integer, intent(in) :: fplspopt
                       ! Formal parameter of unique index of lspopt

      integer, intent(in) :: fpvspopt
                       ! Formal parameter of unique index of vspopt

      integer, intent(in) :: fpsfcopt
                       ! Formal parameter of unique index of sfcopt

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fpaslopt
                       ! Formal parameter of unique index of aslopt

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

      integer, intent(in) :: nqa(0:4)
                       ! Number of types of aerosol

      integer, intent(in) :: nund
                       ! Number of soil and sea layers

      integer, intent(in) :: nlev
                       ! Horizontally averaged vertical dimension

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

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: fc(0:ni+1,0:nj+1,1:2)
                       ! 0.25 x Coriolis parameters

! Output variables

      character(len=5), intent(out) :: fmois
                       ! Control flag of air moisture

      integer, intent(out) :: ksp0(1:2)
                       ! Index of lowest vertical sponge level

      integer, intent(out) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

      real, intent(out) :: ubr(0:ni+1,0:nj+1,1:nk)
                       ! Base state x components of velocity

      real, intent(out) :: vbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state y components of velocity

      real, intent(out) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(out) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(out) :: qvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state water vapor mixing ratio

      real, intent(out) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(out) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(out) :: rst8u(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at u points

      real, intent(out) :: rst8v(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at v points

      real, intent(out) :: rst8w(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at w points

      real, intent(out) :: rcsq(0:ni+1,0:nj+1,1:nk)
                       ! rbr x sound wave speed squared

      real, intent(out) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at present

      real, intent(out) :: up(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at past

      real, intent(out) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at present

      real, intent(out) :: vp(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at past

      real, intent(out) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at present

      real, intent(out) :: wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

      real, intent(out) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at present

      real, intent(out) :: ppp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at past

      real, intent(out) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at present

      real, intent(out) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(out) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at present

      real, intent(out) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(out) :: qwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at present

      real, intent(out) :: qwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at past

      real, intent(out) :: nwtr(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at present

      real, intent(out) :: nwtrp(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at past

      real, intent(out) :: qice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at present

      real, intent(out) :: qicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at past

      real, intent(out) :: nice(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at present

      real, intent(out) :: nicep(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at past

      real, intent(out) :: qcwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at present

      real, intent(out) :: qcwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at past

      real, intent(out) :: qcice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at present

      real, intent(out) :: qcicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at past

      real, intent(out) :: qasl(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at present

      real, intent(out) :: qaslp(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at past

      real, intent(out) :: qt(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at present

      real, intent(out) :: qtp(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at past

      real, intent(out) :: tke(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at present

      real, intent(out) :: tkep(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at past

      real, intent(out) :: rbcx(1:ni)
                       ! Relaxed lateral sponge damping coefficients
                       ! in x direction

      real, intent(out) :: rbcy(1:nj)
                       ! Relaxed lateral sponge damping coefficients
                       ! in y direction

      real, intent(out) :: rbcxy(1:ni,1:nj)
                       ! Relaxed lateral sponge damping coefficients

      real, intent(out) :: rbct(1:ni,1:nj,1:nk,1:2)
                       ! Relaxed top sponge damping coefficients

      real, intent(out) :: ucpx(1:nj,1:nk,1:2)
                       ! Phase speed of x components of velocity
                       ! on west and east boundary

      real, intent(out) :: ucpy(1:ni,1:nk,1:2)
                       ! Phase speed of x components of velocity
                       ! on south and north boundary

      real, intent(out) :: vcpx(1:nj,1:nk,1:2)
                       ! Phase speed of y components of velocity
                       ! on west and east boundary

      real, intent(out) :: vcpy(1:ni,1:nk,1:2)
                       ! Phase speed of y components of velocity
                       ! on south and north boundary

      real, intent(out) :: wcpx(1:nj,1:nk,1:2)
                       ! Phase speed of z components of velocity
                       ! on west and east boundary

      real, intent(out) :: wcpy(1:ni,1:nk,1:2)
                       ! Phase speed of z components of velocity
                       ! on south and north boundary

      real, intent(out) :: pcpx(1:nj,1:nk,1:2)
                       ! Phase speed of pressure
                       ! on west and east boundary

      real, intent(out) :: pcpy(1:ni,1:nk,1:2)
                       ! Phase speed of pressure
                       ! on south and north boundary

      real, intent(out) :: ptcpx(1:nj,1:nk,1:2)
                       ! Phase speed of potential temperature
                       ! on west and east boundary

      real, intent(out) :: ptcpy(1:ni,1:nk,1:2)
                       ! Phase speed of potential temperature
                       ! on south and north boundary

      real, intent(out) :: qvcpx(1:nj,1:nk,1:2)
                       ! Phase speed of water vapor mixing ratio
                       ! on west and east boundary

      real, intent(out) :: qvcpy(1:ni,1:nk,1:2)
                       ! Phase speed of water vapor mixing ratio
                       ! on south and north boundary

      real, intent(out) :: qwcpx(1:nj,1:nk,1:2,1:nqw)
                       ! Phase speed of water hydrometeor
                       ! on west and east boundary

      real, intent(out) :: qwcpy(1:ni,1:nk,1:2,1:nqw)
                       ! Phase speed of water hydrometeor
                       ! on south and north boundary

      real, intent(out) :: nwcpx(1:nj,1:nk,1:2,1:nnw)
                       ! Phase speed of water concentrations
                       ! on west and east boundary

      real, intent(out) :: nwcpy(1:ni,1:nk,1:2,1:nnw)
                       ! Phase speed of water concentrations
                       ! on south and north boundary

      real, intent(out) :: qicpx(1:nj,1:nk,1:2,1:nqi)
                       ! Phase speed of ice hydrometeor
                       ! on west and east boundary

      real, intent(out) :: qicpy(1:ni,1:nk,1:2,1:nqi)
                       ! Phase speed of ice hydrometeor
                       ! on south and north boundary

      real, intent(out) :: nicpx(1:nj,1:nk,1:2,1:nni)
                       ! Phase speed of ice concentrations
                       ! on west and east boundary

      real, intent(out) :: nicpy(1:ni,1:nk,1:2,1:nni)
                       ! Phase speed of ice concentrations
                       ! on south and north boundary

      real, intent(out) :: qcwcpx(1:nj,1:nk,1:2,1:nqw)
                       ! Phase speed of charging distribution for water
                       ! on west and east boundary

      real, intent(out) :: qcwcpy(1:ni,1:nk,1:2,1:nqw)
                       ! Phase speed of charging distribution for water
                       ! on south and north boundary

      real, intent(out) :: qcicpx(1:nj,1:nk,1:2,1:nqi)
                       ! Phase speed of charging distribution for ice
                       ! on west and east boundary

      real, intent(out) :: qcicpy(1:ni,1:nk,1:2,1:nqi)
                       ! Phase speed of charging distribution for ice
                       ! on south and north boundary

      real, intent(out) :: qacpx(1:nj,1:nk,1:2,1:nqa(0))
                       ! Phase speed of aerosol mixing ratio
                       ! on west and east boundary

      real, intent(out) :: qacpy(1:ni,1:nk,1:2,1:nqa(0))
                       ! Phase speed of aerosol mixing ratio
                       ! on south and north boundary

      real, intent(out) :: qtcpx(1:nj,1:nk,1:2)
                       ! Phase speed of tracer mixing ratio
                       ! on west and east boundary

      real, intent(out) :: qtcpy(1:ni,1:nk,1:2)
                       ! Phase speed of tracer mixing ratio
                       ! on south and north boundary

      real, intent(out) :: tkecpx(1:nj,1:nk,1:2)
                       ! Phase speed of turbulent kinetic energy
                       ! on west and east boundary

      real, intent(out) :: tkecpy(1:ni,1:nk,1:2)
                       ! Phase speed of turbulent kinetic energy
                       ! on south and north boundary

      real, intent(out) :: maxvl(0:ni+1,0:nj+1,1:nk)
                       ! Maximum instantaneous wind velocity

      real, intent(out) :: prwtr(0:ni+1,0:nj+1,1:2,1:nqw)
                       ! Precipitation and accumulation for water

      real, intent(out) :: price(0:ni+1,0:nj+1,1:2,1:nqi)
                       ! Precipitation and accumulation for ice

      real, intent(out) :: pdia(0:ni+1,0:nj+1,1:nk)
                       ! Diabatic value

      real, intent(out) :: albe(0:ni+1,0:nj+1)
                       ! Albedo

      real, intent(out) :: beta(0:ni+1,0:nj+1)
                       ! Evapotranspiration efficiency

      real, intent(out) :: z0m(0:ni+1,0:nj+1)
                       ! Roughness length for velocity

      real, intent(out) :: z0h(0:ni+1,0:nj+1)
                       ! Roughness length for scalar

      real, intent(out) :: cap(0:ni+1,0:nj+1)
                       ! Thermal capacity

      real, intent(out) :: nuu(0:ni+1,0:nj+1)
                       ! Thermal diffusivity

      real, intent(out) :: kai(0:ni+1,0:nj+1)
                       ! Sea ice distribution

      real, intent(out) :: tund(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at present

      real, intent(out) :: tundp(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at past

! Internal shared variables

      integer iniopt   ! Option for model initialization
      integer movopt   ! Option for grid moving
      integer lspopt   ! Option for lateral sponge damping
      integer vspopt   ! Option for vertical sponge damping
      integer sfcopt   ! Option for surface physics
      integer cphopt   ! Option for cloud micro physics
      integer aslopt   ! Option for aerosol processes

      integer stat     ! Runtime status

      real, intent(inout) :: z1d(0:nlev)
                       ! Horizontally averaged z physical coordinates

      real, intent(inout) :: u1d(0:nlev)
                       ! Horizontally averaged x components of velocity

      real, intent(inout) :: v1d(0:nlev)
                       ! Horizontally averaged y components of velocity

      real, intent(inout) :: p1d(0:nlev)
                       ! Horizontally averaged pressure

      real, intent(inout) :: pt1d(0:nlev)
                       ! Horizontally averaged potential temrerature

      real, intent(inout) :: qv1d(0:nlev)
                       ! Horizontally averaged water vapor mixing ratio

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: ltmp1(1:nlev)
                       ! Temporary array

      real, intent(inout) :: ltmp2(1:nlev)
                       ! Temporary array

      real, intent(inout) :: ltmp3(1:nlev)
                       ! Temporary array

! Remark

!     rcsq: This variable is also temporary.

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpiniopt,iniopt)
      call getiname(fpmovopt,movopt)
      call getiname(fplspopt,lspopt)
      call getiname(fpvspopt,vspopt)
      call getiname(fpsfcopt,sfcopt)
      call getiname(fpcphopt,cphopt)
      call getiname(fpaslopt,aslopt)

! -----

! Initialize the surface physical parameters.

      if(sfcopt.ge.1) then

        call s_inisfc(idsfcdat,idwbc,idebc,idexbopt,                    &
     &                idlnduse,iddstopt,idzsfc,idgralbe,idgrbeta,       &
     &                idgrz0m,idgrz0h,idgrcap,idgrnuu,ni,nj,nk,zph,     &
     &                land,albe,beta,z0m,z0h,cap,nuu,kai,tmp1)

      end if

! -----

!! Initialize the variables with reading out the data from the sounding
!! file.

      if(iniopt.eq.1) then

! Read out the data from the sounding file.

        stat=0

        if(mype.eq.root) then

          call s_rdsnd(idexprim,idcrsdir,idncexp,idnccrs,               &
     &                 stat,nlev/4,ltmp1,u1d(1),v1d(1),pt1d(1),qv1d(1))

        end if

        call chkerr(stat)

        if(stat.lt.0) then

          call destroy('rdsnd   ',5,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        call chkstd(root)

! -----

! Broadcast the horizontally averaged variables to the other processor
! elements.

        call s_castvar(1,nlev/4,1,1,1,1,ltmp1)
        call s_castvar(1,nlev/4,1,1,1,1,u1d(1))
        call s_castvar(1,nlev/4,1,1,1,1,v1d(1))
        call s_castvar(1,nlev/4,1,1,1,1,pt1d(1))
        call s_castvar(1,nlev/4,1,1,1,1,qv1d(1))

! -----

! Interpolate the sounding data to fine interval levels vartically.

        call s_vintsnd(idsndtyp,nlev/4,ltmp1,ltmp2,                     &
     &                 nlev,u1d(1),v1d(1),pt1d(1),qv1d(1),z1d(1))

! -----

! Be averaged the sounding variables to the horizontally averaged
! variables.

        call set1d(idsndtyp,idzsnd0,idpsnd0,idthresq,fmois,nlev,z1d,    &
     &             u1d,v1d,pt1d,qv1d,p1d,ltmp1,ltmp2,ltmp3)

! -----

! Set the grid moving velocity.

        if(movopt.eq.1) then

          call move2d(idumove,idvmove,nlev,u1d,v1d)

        end if

! -----

! Extract the base state variables from interpolated sounding data.

        call intrpsnd(idsndtyp,ni,nj,nk,zph,ubr,vbr,pbr,ptbr,qvbr,      &
     &                tmp1,tmp2,tmp3,nlev,z1d,u1d,v1d,p1d,pt1d,qv1d)

! -----

! Set the base state variables.

        call setbase(ni,nj,nk,zph,ubr,vbr,pbr,ptbr,qvbr,rbr,            &
     &               tmp1,tmp2,tmp3)

! -----

! Calculate the base state density multiplyed by Jacobian.

        call baserho(idadvopt,idsmtopt,ni,nj,nk,jcb,rbr,rst,            &
     &               rst8u,rst8v,rst8w)

! -----

! Set the initial velocity and the scalar perturbations.

        call setvar1d(iddmpvar,idadvopt,idcphopt,                       &
     &                idqcgopt,idaslopt,idtrkopt,idtubopt,fmois,        &
     &                ni,nj,nk,nqw,nnw,nqi,nni,nqa,zph,ubr,vbr,qvbr,    &
     &                u,up,v,vp,w,wp,pp,ppp,ptp,ptpp,qv,qvp,            &
     &                qwtr,qwtrp,nwtr,nwtrp,qice,qicep,nice,nicep,      &
     &                qcwtr,qcwtrp,qcice,qcicep,qasl,qaslp,qt,qtp,      &
     &                tke,tkep,maxvl,prwtr,price)

! -----

!! -----

!! Initialize the variables with reading out the data from the restart
!! file.

      else if(iniopt.eq.2.or.iniopt.eq.12) then

! Read out the data from the restart file.

        if(iniopt.eq.2) then

          call rdres(idexprim,idcrsdir,iddmpvar,idncexp,idnccrs,        &
     &               idwbc,idebc,idsbc,idnbc,idsfcopt,idadvopt,idcphopt,&
     &               idqcgopt,idaslopt,idtrkopt,idtubopt,iddiaopt,      &
     &               idstime,'check',fmois,ni,nj,nk,nqw,nnw,nqi,nni,nqa,&
     &               nund,ubr,vbr,pbr,ptbr,qvbr,u,up,v,vp,w,wp,         &
     &               pp,ppp,ptp,ptpp,qv,qvp,qwtr,qwtrp,nwtr,nwtrp,      &
     &               qice,qicep,nice,nicep,qcwtr,qcwtrp,qcice,qcicep,   &
     &               qasl,qaslp,qt,qtp,tke,tkep,ucpx,ucpy,vcpx,vcpy,    &
     &               wcpx,wcpy,pcpx,pcpy,ptcpx,ptcpy,qvcpx,qvcpy,       &
     &               qwcpx,qwcpy,nwcpx,nwcpy,qicpx,qicpy,nicpx,nicpy,   &
     &               qcwcpx,qcwcpy,qcicpx,qcicpy,qacpx,qacpy,           &
     &               qtcpx,qtcpy,tkecpx,tkecpy,maxvl,prwtr,price,pdia,  &
     &               z0m,z0h,tund,tundp)

        else

          call rdres(idexprim,idcrsdir,iddmpvar,idncexp,idnccrs,        &
     &               idwbc,idebc,idsbc,idnbc,idsfcopt,idadvopt,idcphopt,&
     &               idqcgopt,idaslopt,idtrkopt,idtubopt,iddiaopt,      &
     &               idstime,'skip ',fmois,ni,nj,nk,nqw,nnw,nqi,nni,nqa,&
     &               nund,ubr,vbr,pbr,ptbr,qvbr,u,up,v,vp,w,wp,         &
     &               pp,ppp,ptp,ptpp,qv,qvp,qwtr,qwtrp,nwtr,nwtrp,      &
     &               qice,qicep,nice,nicep,qcwtr,qcwtrp,qcice,qcicep,   &
     &               qasl,qaslp,qt,qtp,tke,tkep,ucpx,ucpy,vcpx,vcpy,    &
     &               wcpx,wcpy,pcpx,pcpy,ptcpx,ptcpy,qvcpx,qvcpy,       &
     &               qwcpx,qwcpy,nwcpx,nwcpy,qicpx,qicpy,nicpx,nicpy,   &
     &               qcwcpx,qcwcpy,qcicpx,qcicpy,qacpx,qacpy,           &
     &               qtcpx,qtcpy,tkecpx,tkecpy,maxvl,prwtr,price,pdia,  &
     &               z0m,z0h,tund,tundp)

        end if

! -----

! Set the base state variables.

        call setbase(ni,nj,nk,zph,ubr,vbr,pbr,ptbr,qvbr,rbr,            &
     &               tmp1,tmp2,tmp3)

! -----

! Calculate the base state density multiplyed by Jacobian.

        call baserho(idadvopt,idsmtopt,ni,nj,nk,jcb,rbr,rst,            &
     &               rst8u,rst8v,rst8w)

! -----

!! -----

!! Initialize the variables with reading out the data from the
!! 3 dimensional initial file.

      else if(iniopt.eq.3) then

! Read out the data from the interpolated GPV data file.

        call rdgpvini(idexprim,idcrsdir,idgpvvar,idncexp,idnccrs,       &
     &                idwlngth,idgsmopt,idcphopt,idhaiopt,idstime,      &
     &                fmois,ni,nj,nk,nqw,nqi,ubr,vbr,pbr,ptbr,qvbr,     &
     &                up,vp,wp,ppp,ptpp,qvp,qwtrp,qicep,tmp1)

! -----

! Read out the data from the interpolated aerosol data file.

        if(aslopt.ge.1) then

          call rdaslini(idexprim,idcrsdir,idncexp,idnccrs,idwlngth,     &
     &                  idgsmopt,idstime,ni,nj,nk,nqa,qaslp,tmp1)

        end if

! -----

! Set the base state variables.

        call setbase(ni,nj,nk,zph,ubr,vbr,pbr,ptbr,qvbr,rbr,            &
     &               tmp1,tmp2,tmp3)

! -----

! Calculate the base state density multiplyed by Jacobian.

        call baserho(idadvopt,idsmtopt,ni,nj,nk,jcb,rbr,rst,            &
     &               rst8u,rst8v,rst8w)

! -----

! Set the initial velocity and the scalar perturbations.

        call setvar3d(idgpvvar,iddmpvar,idadvopt,idcphopt,idhaiopt,     &
     &                idqcgopt,idaslopt,idtrkopt,idtubopt,fmois,        &
     &                ni,nj,nk,nqw,nnw,nqi,nni,nqa,zph,j31,j32,jcb8w,   &
     &                mf,rbr,up,vp,wp,ppp,ptpp,qvp,qwtrp,qicep,u,v,w,   &
     &                pp,ptp,qv,qwtr,nwtr,nwtrp,qice,nice,nicep,        &
     &                qcwtr,qcwtrp,qcice,qcicep,qasl,qaslp,qt,qtp,      &
     &                tke,tkep,maxvl,prwtr,price,tmp1,tmp2,tmp3,rcsq)

! -----

      end if

!! -----

! Calculate the lateral sponge damping coefficients.

      if(lspopt.ge.1) then

        call lspdmp(idwbc,idebc,idsbc,idnbc,idwdnews,idwdnorm,ni,nj,    &
     &              rbcx,rbcy,rbcxy)

      end if

! -----

! Calculate the vertical sponge damping coefficients.

      if(vspopt.ge.1) then

        call s_vspdmp(idvspopt,idvspgpv,idvspbar,idbotgpv,idbotbar,ksp0,&
     &                ni,nj,nk,zph,rbct,tmp1)

      end if

! -----

! Initialize the soil and sea temperature.

      if(sfcopt.ge.1.and.mod(iniopt,10).ne.2) then

        call s_initund(idsfcdat,idsfcopt,idadvopt,                      &
     &                 iddzgrd,idtgdeep,idsstcst,ni,nj,nk,nund,         &
     &                 pbr,ptbr,ppp,ptpp,land,tund,tundp,tmp1,tmp2)

      end if

! -----

! Calculate the sound wave speed squared.

      call sndwave(ni,nj,nk,pbr,rcsq)

! -----

! Set the bin parameters.

      if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

        call setbin(idbbinw,idsbinw,nqw)

      end if

! -----

! Read in 2 dimensional grid data to the geography file.

      if(mod(iniopt,10).ne.2) then

        call s_outgeo(idexprim,idcrsdir,idncexp,idnccrs,idwlngth,       &
     &                idmfcopt,idcoropt,idsfcopt,iddmpfmt,ni,nj,        &
     &                zph(0,0,2),lat,lon,mf,fc,land)

      end if

! -----

      end subroutine s_inivar

!-----7--------------------------------------------------------------7--

      end module m_inivar
