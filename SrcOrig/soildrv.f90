!***********************************************************************
      module m_soildrv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/10/16
!     Modification: 2001/10/18, 2001/11/14, 2001/12/10, 2002/01/15,
!                   2002/04/02, 2002/07/03, 2002/08/15, 2003/01/20,
!                   2003/03/13, 2003/04/30, 2003/05/19, 2003/07/15,
!                   2003/09/01, 2003/10/31, 2003/11/28, 2003/12/12,
!                   2004/02/01, 2004/03/05, 2004/04/01, 2004/05/07,
!                   2004/05/31, 2004/06/10, 2004/07/01, 2004/08/01,
!                   2004/08/20, 2004/09/01, 2004/09/10, 2004/12/17,
!                   2005/01/14, 2005/01/31, 2005/04/04, 2005/06/10,
!                   2005/08/05, 2006/01/10, 2006/02/13, 2006/04/03,
!                   2006/08/08, 2006/09/21, 2007/01/20, 2007/05/14,
!                   2007/05/21, 2007/06/27, 2007/07/30, 2008/03/12,
!                   2008/04/17, 2008/05/02, 2008/07/01, 2008/08/25,
!                   2009/02/27, 2009/08/20, 2010/02/01, 2011/09/22,
!                   2011/11/10, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the inferior procedures for the surface physics.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkrain
      use m_cloudcov
      use m_comindx
      use m_comkind
      use m_forcesfc
      use m_getiname
      use m_getrich
      use m_getzlow
      use m_heatsfc
      use m_outpbl
      use m_pbldrv
      use m_radiat
      use m_roughitr
      use m_roughnxt
      use m_setsfc
      use m_sfcflx
      use m_steptund

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: soildrv, s_soildrv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface soildrv

        module procedure s_soildrv

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
      subroutine s_soildrv(fptubopt,fmois,pdate,ctime,dtb,dtsoil,stinc, &
     &                     ni,nj,nk,nqw,nqi,nund,zph,lat,lon,j31,j32,   &
     &                     jcb8w,pbr,ptbr,rbr,rst,rst8u,rst8v,up,vp,wp, &
     &                     ppp,ptpp,qvp,prwtr,price,qallp,land,albe,    &
     &                     beta,cap,nuu,kai,sst,sstd,uf,vf,ptpf,qvf,    &
     &                     z0m,z0h,tundp,tundf,ufrc,vfrc,ptfrc,qvfrc,   &
     &                     p,t,ptv,qvsfc,tice,va,rch,cm,ch,ce,ct,cq,    &
     &                     hs,le,rgd,rsd,rld,rlu,cdl,cdm,cdh,fall,za,   &
     &                     tmp1,tmp2,tmp3,tmp4)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      character(len=12), intent(in) :: pdate
                       ! Forecast date at 1 step past
                       ! with Gregorian calendar, yyyymmddhhmm

      integer, intent(in) :: fptubopt
                       ! Formal parameter of unique index of tubopt

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

      integer, intent(in) :: nqi
                       ! Number of categories of ice hydrometeor

      integer, intent(in) :: nund
                       ! Number of soil and sea layers

      integer, intent(in) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: dtsoil
                       ! Time interval of soil temperature calculation

      real, intent(in) :: stinc
                       ! Lapse of forecast time
                       ! from sea surface temperature data reading

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

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: rst8u(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at u points

      real, intent(in) :: rst8v(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at v points

      real, intent(in) :: up(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at past

      real, intent(in) :: vp(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at past

      real, intent(in) :: wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

      real, intent(in) :: ppp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at past

      real, intent(in) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbarion at past

      real, intent(in) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(in) :: prwtr(0:ni+1,0:nj+1,1:2,1:nqw)
                       ! Precipitation and accumulation for water

      real, intent(in) :: price(0:ni+1,0:nj+1,1:2,1:nqi)
                       ! Precipitation and accumulation for ice

      real, intent(in) :: qallp(0:ni+1,0:nj+1,1:nk)
                       ! Total water and ice mixing ratio at past

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

      real, intent(inout) :: uf(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at future

      real, intent(inout) :: vf(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at future

      real, intent(inout) :: ptpf(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbarion at future

      real, intent(inout) :: qvf(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at future

      real, intent(inout) :: z0m(0:ni+1,0:nj+1)
                       ! Roughness length for velocity

      real, intent(inout) :: z0h(0:ni+1,0:nj+1)
                       ! Roughness length for scalar

      real, intent(inout) :: tundp(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at past

      real, intent(inout) :: tundf(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at future

! Output variables

      real, intent(out) :: ufrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in u equation

      real, intent(out) :: vfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in v equation

      real, intent(out) :: ptfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in potential temperature equation

      real, intent(out) :: qvfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in water vapor mixing ratio equation

! Internal shared variables

      integer tubopt   ! Option for turbulent mixing

      real, intent(inout) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure

      real, intent(inout) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

      real, intent(inout) :: ptv(0:ni+1,0:nj+1,1:nk)
                       ! Virtual potential temperature

      real, intent(inout) :: qvsfc(0:ni+1,0:nj+1)
                       ! Water vapor mixing ratio on surface

      real, intent(inout) :: tice(0:ni+1,0:nj+1)
                       ! Mixed ice surface temperature

      real, intent(inout) :: va(0:ni+1,0:nj+1)
                       ! Magnitude of velocity at lowest plane

      real, intent(inout) :: rch(0:ni+1,0:nj+1)
                       ! Bulk Richardson number

      real, intent(inout) :: cm(0:ni+1,0:nj+1)
                       ! Bulk coefficient for velocity

      real, intent(inout) :: ch(0:ni+1,0:nj+1)
                       ! Bulk coefficient for scalar

      real, intent(inout) :: ce(0:ni+1,0:nj+1)
                       ! Exchange coefficient of surface momentum flux

      real, intent(inout) :: ct(0:ni+1,0:nj+1)
                       ! Exchange coefficient of surface heat flux

      real, intent(inout) :: cq(0:ni+1,0:nj+1)
                       ! Exchange coefficient of surface moisture flux

      real, intent(inout) :: hs(0:ni+1,0:nj+1)
                       ! Sensible heat

      real, intent(inout) :: le(0:ni+1,0:nj+1)
                       ! Latent heat

      real, intent(inout) :: rgd(0:ni+1,0:nj+1)
                       ! Global solar radiation

      real, intent(inout) :: rsd(0:ni+1,0:nj+1)
                       ! Net downward short wave radiation

      real, intent(inout) :: rld(0:ni+1,0:nj+1)
                       ! Downward long wave radiation

      real, intent(inout) :: rlu(0:ni+1,0:nj+1)
                       ! Upward long wave radiation

      real, intent(inout) :: cdl(0:ni+1,0:nj+1)
                       ! Cloud cover in lower layer

      real, intent(inout) :: cdm(0:ni+1,0:nj+1)
                       ! Cloud cover in middle layer

      real, intent(inout) :: cdh(0:ni+1,0:nj+1)
                       ! Cloud cover in upper layer

      real, intent(inout) :: fall(0:ni+1,0:nj+1)
                       ! Precipitation flag

      real, intent(inout) :: za(0:ni+1,0:nj+1)
                       ! z physical coordinates at lowest plane

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp4(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Remark

!     ufrc,vfrc: These variables are also temporary, because they are
!                not used again.

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fptubopt,tubopt)

! -----

! Calculate the z physical coordinates at lowest plane.

      call getzlow(ni,nj,nk,zph,za)

! -----

! Check the precipitation on the surface.

      call chkrain(idcphopt,idhaiopt,fmois,ni,nj,nqw,nqi,prwtr,price,   &
     &             fall)

! -----

! Calculate the magnitude of velocity, virtual potential temperature and
! water vapor mixing ratio on the surface.

      call s_setsfc(idlevpbl,idtubopt,idcphopt,fmois,                   &
     &              ni,nj,nk,nund,pbr,ptbr,up,vp,wp,ppp,ptpp,qvp,       &
     &              land,beta,kai,tundp,fall,p,t,ptv,qvsfc,tice,va,     &
     &              tmp1,tmp2(0,0,1),tmp2(0,0,2),tmp2(0,0,3))

! -----

! Calculate the cloud cover.

      call s_cloudcov(idcphopt,iddz,fmois,'rh ',ni,nj,nk,               &
     &                zph,rst,p,t,qvp,qallp,cdl,cdm,cdh,tmp1,           &
     &                tmp2(0,0,1),tmp2(0,0,2),tmp2(0,0,3),tmp2(0,0,4),  &
     &                tmp3,tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3))

! -----

! Calculate the short and long wave radiation.

      call s_radiat(idcphopt,fmois,pdate,ni,nj,nk,nund,zph,lat,lon,p,t, &
     &              qvp,land,albe,kai,tundp,tice,cdl,cdm,cdh,fall,      &
     &              rgd,rsd,rld,rlu,tmp1,tmp2(0,0,1),tmp2(0,0,2))

! -----

! Calculate the roughness length for velocity on the sea surface by
! iteration at forecast start time.

      if(ctime.eq.0_i8) then

        call s_roughitr(ni,nj,nk,za,land,kai,ptv,va,z0m,z0h,rch,cm,ch,  &
     &                  tmp1(0,0,1))

      end if

! -----

! Calculate the bulk Richardson number on the surface.

      call getrich(ni,nj,nk,za,land,kai,z0m,z0h,ptv,va,rch)

! -----

! Calculate the friction velocity and the common coefficients of surface
! flux calculation.

      call sfcflx(ni,nj,nk,za,rbr,land,kai,z0m,z0h,va,rch,cm,ch,        &
     &            ce,ct,cq)

! -----

! Calculate the sensible and latent heat on the surface.

      call heatsfc(fmois,ni,nj,nk,nund,t,qvp,qvsfc,ct,cq,land,kai,      &
     &             tundp,tice,hs,le)

! -----

! Reset the roughness parameter on the sea surface to the next time
! step.

      call roughnxt(ni,nj,land,va,cm,z0m,z0h)

! -----

! Get the surface flux to bottom boundary.

      call forcesfc(fmois,ni,nj,nk,j31,j32,ptbr,up,vp,wp,ptpp,qvp,      &
     &              ptv,qvsfc,ce,ct,cq,ufrc,vfrc,ptfrc,qvfrc)

! -----

! Read in the surface variables to the dumped file.

      call s_outpbl(iddmpvar,iddmplev,fmois,                            &
     &              ni,nj,nk,nund,za,lon,p,up,vp,qvp,                   &
     &              ufrc,vfrc,ptfrc,qvfrc,land,kai,z0m,z0h,ptv,         &
     &              qvsfc,rch,cm,ch,tundp,tice,hs,le,rgd,rsd,rld,rlu,   &
     &              cdl,cdm,cdh,tmp1(0,0,1),tmp1(0,0,2),tmp1(0,0,3),    &
     &              tmp2(0,0,1),tmp2(0,0,2),tmp2(0,0,3),tmp2(0,0,4),    &
     &              tmp3(0,0,1),tmp3(0,0,2),tmp3(0,0,3),tmp3(0,0,4),    &
     &              tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4))

! -----

!! Solve the surface processes.

      if(ctime.ne.0_i8) then

! Solve the soil and sea temperature to the next time step.

        if(dtsoil.gt.0.e0) then

          call steptund(idsfcopt,iddzgrd,iddzsea,dtsoil,stinc,          &
     &                  ni,nj,nk,nund,t,land,cap,nuu,sst,sstd,          &
     &                  hs,le,rsd,rld,rlu,tundp,tundf,                  &
     &                  tmp1,tmp2,tmp3,tmp4)

        end if

! -----

! Calculate the diffusion term in planetary boundary layer.

        if(tubopt.eq.0) then

          call pbldrv(idlevpbl,fmois,dtb,ni,nj,nk,zph,jcb8w,            &
     &                ptbr,rbr,rst,rst8u,rst8v,ce,ct,cq,qvsfc,ptv,      &
     &                uf,vf,ptpf,qvf,tmp1,tmp2,tmp3,tmp4,ufrc,vfrc)

        end if

! -----

      end if

!! -----

      end subroutine s_soildrv

!-----7--------------------------------------------------------------7--

      end module m_soildrv
