!***********************************************************************
      module m_comslv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/05/19
!     Modification: 2003/07/15, 2003/11/05, 2003/12/12, 2004/02/01,
!                   2004/04/15, 2004/05/31, 2004/06/10, 2004/08/01,
!                   2004/08/20, 2005/01/14, 2006/01/10, 2006/02/13,
!                   2006/06/21, 2006/07/21, 2006/11/06, 2007/01/20,
!                   2007/07/30, 2008/05/02, 2008/07/01, 2008/08/25,
!                   2008/12/11, 2009/02/27, 2009/03/23, 2011/08/18,
!                   2011/09/22, 2011/11/10

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the array and variable for solver.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      public

! Exceptional access control

!     none

!-----7--------------------------------------------------------------7--

! Module variables

      character(len=5) fmois
                       ! Control flag of air moisture

      integer ksp0(1:2)
                       ! Index of lowest vertical sponge level

      real area(0:4)   ! Area of each boundary plane

      integer, allocatable, save :: land(:,:)
                       ! Land use of surface

      real, allocatable, save :: zph(:,:,:)
                       ! z physical coordinates

      real, allocatable, save :: zsth(:)
                       ! 1 dimensional stretched z coordinates

      real, allocatable, save :: lat(:,:)
                       ! Latitude

      real, allocatable, save :: lon(:,:)
                       ! Longitude

      real, allocatable, save :: j31(:,:,:)
                       ! z-x components of Jacobian

      real, allocatable, save :: j32(:,:,:)
                       ! z-y components of Jacobian

      real, allocatable, save :: jcb(:,:,:)
                       ! Jacobian

      real, allocatable, save :: jcb8u(:,:,:)
                       ! Jacobian at u points

      real, allocatable, save :: jcb8v(:,:,:)
                       ! Jacobian at v points

      real, allocatable, save :: jcb8w(:,:,:)
                       ! Jacobian at w points

      real, allocatable, save :: mf(:,:)
                       ! Map scale factors

      real, allocatable, save :: mf8u(:,:)
                       ! Map scale factors at u points

      real, allocatable, save :: mf8v(:,:)
                       ! Map scale factors at v points

      real, allocatable, save :: rmf(:,:,:)
                       ! Related parameters of map scale factors

      real, allocatable, save :: rmf8u(:,:,:)
                       ! Related parameters of map scale factors
                       ! at u points

      real, allocatable, save :: rmf8v(:,:,:)
                       ! Related parameters of map scale factors
                       ! at v points

      real, allocatable, save :: fc(:,:,:)
                       ! 0.25 x Coriolis parameters

      real, allocatable, save :: ubr(:,:,:)
                       ! Base state x components of velocity

      real, allocatable, save :: vbr(:,:,:)
                       ! Base state y components of velocity

      real, allocatable, save :: pbr(:,:,:)
                       ! Base state pressure

      real, allocatable, save :: ptbr(:,:,:)
                       ! Base state potential temperature

      real, allocatable, save :: qvbr(:,:,:)
                       ! Base state water vapor mixing ratio

      real, allocatable, save :: rbr(:,:,:)
                       ! Base state density

      real, allocatable, save :: rst(:,:,:)
                       ! Base state density x Jacobian

      real, allocatable, save :: rst8u(:,:,:)
                       ! Base state density x Jacobian at u points

      real, allocatable, save :: rst8v(:,:,:)
                       ! Base state density x Jacobian at v points

      real, allocatable, save :: rst8w(:,:,:)
                       ! Base state density x Jacobian at w points

      real, allocatable, save :: rcsq(:,:,:)
                       ! rbr x sound wave speed squared

      real, allocatable, save :: u(:,:,:)
                       ! x components of velocity at present

      real, allocatable, save :: up(:,:,:)
                       ! x components of velocity at past

      real, allocatable, save :: v(:,:,:)
                       ! y components of velocity at present

      real, allocatable, save :: vp(:,:,:)
                       ! y components of velocity at past

      real, allocatable, save :: w(:,:,:)
                       ! z components of velocity at present

      real, allocatable, save :: wp(:,:,:)
                       ! z components of velocity at past

      real, allocatable, save :: pp(:,:,:)
                       ! Pressure perturbation at present

      real, allocatable, save :: ppp(:,:,:)
                       ! Pressure perturbation at past

      real, allocatable, save :: ptp(:,:,:)
                       ! Potential temperature perturbation at present

      real, allocatable, save :: ptpp(:,:,:)
                       ! Potential temperature perturbation at past

      real, allocatable, save :: qv(:,:,:)
                       ! Water vapor mixing ratio at present

      real, allocatable, save :: qvp(:,:,:)
                       ! Water vapor mixing ratio at past

      real, allocatable, save :: qwtr(:,:,:,:)
                       ! Water hydrometeor at present

      real, allocatable, save :: qwtrp(:,:,:,:)
                       ! Water hydrometeor at past

      real, allocatable, save :: nwtr(:,:,:,:)
                       ! Water concentrations at present

      real, allocatable, save :: nwtrp(:,:,:,:)
                       ! Water concentrations at past

      real, allocatable, save :: qice(:,:,:,:)
                       ! Ice hydrometeor at present

      real, allocatable, save :: qicep(:,:,:,:)
                       ! Ice hydrometeor at past

      real, allocatable, save :: nice(:,:,:,:)
                       ! Ice concentrations at present

      real, allocatable, save :: nicep(:,:,:,:)
                       ! Ice concentrations at past

      real, allocatable, save :: qcwtr(:,:,:,:)
                       ! Charging distribution for water at present

      real, allocatable, save :: qcwtrp(:,:,:,:)
                       ! Charging distribution for water at past

      real, allocatable, save :: qcice(:,:,:,:)
                       ! Charging distribution for ice at present

      real, allocatable, save :: qcicep(:,:,:,:)
                       ! Charging distribution for ice at past

      real, allocatable, save :: qasl(:,:,:,:)
                       ! Aerosol mixing ratio at present

      real, allocatable, save :: qaslp(:,:,:,:)
                       ! Aerosol mixing ratio at past

      real, allocatable, save :: qt(:,:,:)
                       ! Tracer mixing ratio at present

      real, allocatable, save :: qtp(:,:,:)
                       ! Tracer mixing ratio at past

      real, allocatable, save :: tke(:,:,:)
                       ! Turbulent kinetic energy at present

      real, allocatable, save :: tkep(:,:,:)
                       ! Turbulent kinetic energy at past

      real, allocatable, save :: ucpx(:,:,:)
                       ! Phase speed of x components of velocity
                       ! on west and east boundary

      real, allocatable, save :: ucpy(:,:,:)
                       ! Phase speed of x components of velocity
                       ! on south and north boundary

      real, allocatable, save :: vcpx(:,:,:)
                       ! Phase speed of y components of velocity
                       ! on west and east boundary

      real, allocatable, save :: vcpy(:,:,:)
                       ! Phase speed of y components of velocity
                       ! on south and north boundary

      real, allocatable, save :: wcpx(:,:,:)
                       ! Phase speed of z components of velocity
                       ! on west and east boundary

      real, allocatable, save :: wcpy(:,:,:)
                       ! Phase speed of z components of velocity
                       ! on south and north boundary

      real, allocatable, save :: pcpx(:,:,:)
                       ! Phase speed of pressure
                       ! on west and east boundary

      real, allocatable, save :: pcpy(:,:,:)
                       ! Phase speed of pressure
                       ! on south and north boundary

      real, allocatable, save :: ptcpx(:,:,:)
                       ! Phase speed of potential temperature
                       ! on west and east boundary

      real, allocatable, save :: ptcpy(:,:,:)
                       ! Phase speed of potential temperature
                       ! on south and north boundary

      real, allocatable, save :: qvcpx(:,:,:)
                       ! Phase speed of water vapor mixing ratio
                       ! on west and east boundary

      real, allocatable, save :: qvcpy(:,:,:)
                       ! Phase speed of water vapor mixing ratio
                       ! on south and north boundary

      real, allocatable, save :: qwcpx(:,:,:,:)
                       ! Phase speed of water hydrometeor
                       ! on west and east boundary

      real, allocatable, save :: qwcpy(:,:,:,:)
                       ! Phase speed of water hydrometeor
                       ! on south and north boundary

      real, allocatable, save :: nwcpx(:,:,:,:)
                       ! Phase speed of water concentrations
                       ! on west and east boundary

      real, allocatable, save :: nwcpy(:,:,:,:)
                       ! Phase speed of water concentrations
                       ! on south and north boundary

      real, allocatable, save :: qicpx(:,:,:,:)
                       ! Phase speed of ice hydrometeor
                       ! on west and east boundary

      real, allocatable, save :: qicpy(:,:,:,:)
                       ! Phase speed of ice hydrometeor
                       ! on south and north boundary

      real, allocatable, save :: nicpx(:,:,:,:)
                       ! Phase speed of ice concentrations
                       ! on west and east boundary

      real, allocatable, save :: nicpy(:,:,:,:)
                       ! Phase speed of ice concentrations
                       ! on south and north boundary

      real, allocatable, save :: qcwcpx(:,:,:,:)
                       ! Phase speed of charging distribution for water
                       ! on west and east boundary

      real, allocatable, save :: qcwcpy(:,:,:,:)
                       ! Phase speed of charging distribution for water
                       ! on south and north boundary

      real, allocatable, save :: qcicpx(:,:,:,:)
                       ! Phase speed of charging distribution for ice
                       ! on west and east boundary

      real, allocatable, save :: qcicpy(:,:,:,:)
                       ! Phase speed of charging distribution for ice
                       ! on south and north boundary

      real, allocatable, save :: qacpx(:,:,:,:)
                       ! Phase speed of aerosol mixing ratio
                       ! on west and east boundary

      real, allocatable, save :: qacpy(:,:,:,:)
                       ! Phase speed of aerosol mixing ratio
                       ! on south and north boundary

      real, allocatable, save :: qtcpx(:,:,:)
                       ! Phase speed of tracer mixing ratio
                       ! on west and east boundary

      real, allocatable, save :: qtcpy(:,:,:)
                       ! Phase speed of tracer mixing ratio
                       ! on south and north boundary

      real, allocatable, save :: tkecpx(:,:,:)
                       ! Phase speed of turbulent kinetic energy
                       ! on west and east boundary

      real, allocatable, save :: tkecpy(:,:,:)
                       ! Phase speed of turbulent kinetic energy
                       ! on south and north boundary

      real, allocatable, save :: rbcx(:)
                       ! Relaxed lateral sponge damping coefficients
                       ! in x direction

      real, allocatable, save :: rbcy(:)
                       ! Relaxed lateral sponge damping coefficients
                       ! in y direction

      real, allocatable, save :: rbcxy(:,:)
                       ! Relaxed lateral sponge damping coefficients

      real, allocatable, save :: rbct(:,:,:,:)
                       ! Relaxed top sponge damping coefficients

      real, allocatable, save :: ugpv(:,:,:)
                       ! x components of velocity of GPV data
                       ! at marked time

      real, allocatable, save :: utd(:,:,:)
                       ! Time tendency of
                       ! x components of velocity of GPV data

      real, allocatable, save :: vgpv(:,:,:)
                       ! y components of velocity of GPV data
                       ! at marked time

      real, allocatable, save :: vtd(:,:,:)
                       ! Time tendency of
                       ! y components of velocity of GPV data

      real, allocatable, save :: wgpv(:,:,:)
                       ! z components of velocity of GPV data
                       ! at marked time

      real, allocatable, save :: wtd(:,:,:)
                       ! Time tendency of
                       ! z components of velocity of GPV data

      real, allocatable, save :: ppgpv(:,:,:)
                       ! Pressure perturbation of GPV data
                       ! at marked time

      real, allocatable, save :: pptd(:,:,:)
                       ! Time tendency of
                       ! pressure perturbation of GPV data

      real, allocatable, save :: ptpgpv(:,:,:)
                       ! Potential temperature perturbation of GPV data
                       ! at marked time

      real, allocatable, save :: ptptd(:,:,:)
                       ! Time tendency of
                       ! potential temperature perturbation of GPV data

      real, allocatable, save :: qvgpv(:,:,:)
                       ! Water vapor mixing ratio of GPV data
                       ! at marked time

      real, allocatable, save :: qvtd(:,:,:)
                       ! Time tendency of
                       ! water vapor mixing ratio of GPV data

      real, allocatable, save :: qwgpv(:,:,:,:)
                       ! Water hydrometeor of GPV data at marked time

      real, allocatable, save :: qwtd(:,:,:,:)
                       ! Time tendency of water hydrometeor of GPV data

      real, allocatable, save :: qigpv(:,:,:,:)
                       ! Ice hydrometeor of GPV data at marked time

      real, allocatable, save :: qitd(:,:,:,:)
                       ! Time tendency of ice hydrometeor of GPV data

      real, allocatable, save :: qagpv(:,:,:,:)
                       ! Mixing ratio of aerosol data at marked time

      real, allocatable, save :: qatd(:,:,:,:)
                       ! Time tendency of mixing ratio of aerosol data

      real, allocatable, save :: urdr(:,:,:)
                       ! x components of velocity of radar data
                       ! at marked time

      real, allocatable, save :: vrdr(:,:,:)
                       ! y components of velocity of radar data
                       ! at marked time

      real, allocatable, save :: wrdr(:,:,:)
                       ! z components of velocity of radar data
                       ! at marked time

      real, allocatable, save :: qwrdr(:,:,:,:)
                       ! Water hydrometeor of radar data
                       ! at marked time

      real, allocatable, save :: qwrtd(:,:,:,:)
                       ! Time tendency of
                       ! water hydrometeor of radar data

      real, allocatable, save :: qirdr(:,:,:,:)
                       ! Ice hydrometeor of radar data
                       ! at marked time

      real, allocatable, save :: qirtd(:,:,:,:)
                       ! Time tendency of
                       ! ice hydrometeor of radar data

      real, allocatable, save :: maxvl(:,:,:)
                       ! Maximum instantaneous wind velocity

      real, allocatable, save :: prwtr(:,:,:,:)
                       ! Precipitation and accumulation for water

      real, allocatable, save :: price(:,:,:,:)
                       ! Precipitation and accumulation for ice

      real, allocatable, save :: qall(:,:,:)
                       ! Total water and ice mixing ratio at present

      real, allocatable, save :: qallp(:,:,:)
                       ! Total water and ice mixing ratio at past

      real, allocatable, save :: pdia(:,:,:)
                       ! Diabatic value

      real, allocatable, save :: albe(:,:)
                       ! Albedo

      real, allocatable, save :: beta(:,:)
                       ! Evapotranspiration efficiency

      real, allocatable, save :: z0m(:,:)
                       ! Roughness length for velocity

      real, allocatable, save :: z0h(:,:)
                       ! Roughness length for scalar

      real, allocatable, save :: cap(:,:)
                       ! Thermal capacity

      real, allocatable, save :: nuu(:,:)
                       ! Thermal diffusivity

      real, allocatable, save :: kai(:,:)
                       ! Sea ice distribution

      real, allocatable, save :: sst(:,:)
                       ! Sea surface temperature of external data
                       ! at marked time

      real, allocatable, save :: sstd(:,:)
                       ! Time tendency of
                       ! sea surface temperature of external data

      real, allocatable, save :: tund(:,:,:)
                       ! Soil and sea temperature at present

      real, allocatable, save :: tundp(:,:,:)
                       ! Soil and sea temperature at past

      real, allocatable, save :: z1d(:)
                       ! Horizontally averaged z physical coordinates

      real, allocatable, save :: u1d(:)
                       ! Horizontally averaged x components of velocity

      real, allocatable, save :: v1d(:)
                       ! Horizontally averaged y components of velocity

      real, allocatable, save :: p1d(:)
                       ! Horizontally averaged pressure

      real, allocatable, save :: pt1d(:)
                       ! Horizontally averaged potential temrerature

      real, allocatable, save :: qv1d(:)
                       ! Horizontally averaged water vapor mixing ratio

      real, allocatable, save :: lon_rdr(:,:)
                       ! Longitude in radar data

      real, allocatable, save :: z_rdr(:,:,:)
                       ! z physical coordinates in radar data

      real, allocatable, save :: u_rdr(:,:,:)
                       ! x components of velocity in radar data

      real, allocatable, save :: v_rdr(:,:,:)
                       ! y components of velocity in radar data

      real, allocatable, save :: w_rdr(:,:,:)
                       ! z components of velocity in radar data

      real, allocatable, save :: qp_rdr(:,:,:)
                       ! Precipitation mixing ratio in radar data

      real, allocatable, save :: tmp1(:,:,:)
                       ! Temporary array

      real, allocatable, save :: tmp2(:,:,:)
                       ! Temporary array

      real, allocatable, save :: tmp3(:,:,:)
                       ! Temporary array

      real, allocatable, save :: tmp4(:,:,:)
                       ! Temporary array

      real, allocatable, save :: tmp5(:,:,:)
                       ! Temporary array

      real, allocatable, save :: tmp6(:,:,:)
                       ! Temporary array

      real, allocatable, save :: tmp7(:,:,:)
                       ! Temporary array

      real, allocatable, save :: tmp8(:,:,:)
                       ! Temporary array

      real, allocatable, save :: tmp9(:,:,:)
                       ! Temporary array

      real, allocatable, save :: tmp10(:,:,:)
                       ! Temporary array

      real, allocatable, save :: tmp11(:,:,:)
                       ! Temporary array

      real, allocatable, save :: tmp12(:,:,:)
                       ! Temporary array

      real, allocatable, save :: tmp13(:,:,:)
                       ! Temporary array

      real, allocatable, save :: tmp14(:,:,:)
                       ! Temporary array

      real, allocatable, save :: tmp15(:,:,:)
                       ! Temporary array

      real, allocatable, save :: tmp16(:,:,:)
                       ! Temporary array

      real, allocatable, save :: tmp17(:,:,:)
                       ! Temporary array

      real, allocatable, save :: tmp18(:,:,:)
                       ! Temporary array

      real, allocatable, save :: tmp19(:,:,:)
                       ! Temporary array

      real, allocatable, save :: qwtmp(:,:,:,:)
                       ! Temporary array

      real, allocatable, save :: nwtmp(:,:,:,:)
                       ! Temporary array

      real, allocatable, save :: qitmp(:,:,:,:)
                       ! Temporary array

      real, allocatable, save :: nitmp(:,:,:,:)
                       ! Temporary array

      real, allocatable, save :: qcwtmp(:,:,:,:)
                       ! Temporary array

      real, allocatable, save :: qcitmp(:,:,:,:)
                       ! Temporary array

      real, allocatable, save :: qatmp(:,:,:,:)
                       ! Temporary array

      real, allocatable, save :: qttmp(:,:,:)
                       ! Temporary array

      real, allocatable, save :: tketmp(:,:,:)
                       ! Temporary array

      real, allocatable, save :: tutmp(:,:,:)
                       ! Temporary array

      real, allocatable, save :: ltmp1(:)
                       ! Temporary array

      real, allocatable, save :: ltmp2(:)
                       ! Temporary array

      real, allocatable, save :: ltmp3(:)
                       ! Temporary array

      real, allocatable, save :: tmp1_rdr(:,:,:)
                       ! Temporary array

! Module procedure

!     none

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

!     none

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

!     none

!-----7--------------------------------------------------------------7--

      end module m_comslv
