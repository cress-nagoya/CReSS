!***********************************************************************
      module m_sfcphy
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/10/16
!     Modification: 2001/10/18, 2001/11/14, 2001/12/10, 2002/01/15,
!                   2002/04/02, 2002/07/03, 2002/08/15, 2003/01/20,
!                   2003/03/13, 2003/04/30, 2003/05/19, 2003/07/15,
!                   2003/09/01, 2003/10/31, 2003/11/28, 2003/12/12,
!                   2004/02/01, 2004/03/05, 2004/04/01, 2004/08/01,
!                   2004/08/20, 2004/09/10, 2004/12/17, 2005/01/14,
!                   2005/08/05, 2006/01/10, 2006/02/13, 2006/04/03,
!                   2006/09/21, 2007/01/20, 2007/07/30, 2008/05/02,
!                   2008/07/01, 2008/08/25, 2009/02/27, 2009/08/20,
!                   2010/02/01, 2011/09/22, 2011/11/10, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     call the surface process driver routine to bridge the complicated
!     procedure.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_comkind
      use m_getiname
      use m_soildrv

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: sfcphy, s_sfcphy

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface sfcphy

        module procedure s_sfcphy

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
      subroutine s_sfcphy(fpadvopt,fmois,pdate,ctime,dtb,dtsoil,stinc,  &
     &                    ni,nj,nk,nqw,nqi,nund,zph,lat,lon,j31,j32,    &
     &                    jcb8w,pbr,ptbr,rbr,rst,rst8u,rst8v,up,vp,wp,  &
     &                    ppp,ptpp,qvp,prwtr,price,qallp,land,albe,beta,&
     &                    cap,nuu,kai,sst,sstd,uf,vf,ptpf,qvf,z0m,z0h,  &
     &                    tundp,tundf,ufrc,vfrc,ptfrc,qvfrc,            &
     &                    tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,                &
     &                    tmp7,tmp8,tmp9,tmp10)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      character(len=12), intent(in) :: pdate
                       ! Forecast date at 1 step past
                       ! with Gregorian calendar, yyyymmddhhmm

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

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

      integer advopt   ! Option for advection scheme

      real dtb_sub     ! Substitute for dtb
      real dtsoil_sub  ! Substitute for dtsoil

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

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpadvopt,advopt)

! -----

! Reset the large time steps interval.

      if(advopt.le.3) then

        dtb_sub=2.e0*dtb
        dtsoil_sub=2.e0*dtsoil

      else

        dtb_sub=dtb
        dtsoil_sub=dtsoil

      end if

! -----

! Call the surface process driver routine to bridge the complicated
! procedure.

      call s_soildrv(idtubopt,fmois,pdate,ctime,                        &
     &               dtb_sub,dtsoil_sub,stinc,ni,nj,nk,nqw,nqi,nund,zph,&
     &               lat,lon,j31,j32,jcb8w,pbr,ptbr,rbr,rst,rst8u,rst8v,&
     &               up,vp,wp,ppp,ptpp,qvp,prwtr,price,qallp,land,albe, &
     &               beta,cap,nuu,kai,sst,sstd,uf,vf,ptpf,qvf,z0m,z0h,  &
     &               tundp,tundf,ufrc,vfrc,ptfrc,qvfrc,tmp1,tmp2,tmp3,  &
     &               tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),   &
     &               tmp4(0,0,5),tmp4(0,0,6),tmp4(0,0,7),tmp5(0,0,1),   &
     &               tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4),tmp5(0,0,5),   &
     &               tmp5(0,0,6),tmp5(0,0,7),tmp6(0,0,1),tmp6(0,0,2),   &
     &               tmp6(0,0,3),tmp6(0,0,4),tmp6(0,0,5),tmp6(0,0,6),   &
     &               tmp7,tmp8,tmp9,tmp10)

! -----

      end subroutine s_sfcphy

!-----7--------------------------------------------------------------7--

      end module m_sfcphy
