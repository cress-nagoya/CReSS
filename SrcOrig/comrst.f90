!***********************************************************************
      module m_comrst
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2005/07/01
!     Modification: 2005/08/05, 2006/01/10, 2006/02/13, 2006/07/21,
!                   2007/01/20, 2007/05/07, 2007/07/30, 2008/05/02,
!                   2008/07/01, 2008/08/25, 2009/02/27, 2011/08/18,
!                   2011/09/22, 2011/11/10

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the array for rstruct.

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

      real, allocatable, save :: maxvl(:,:,:)
                       ! Maximum instantaneous wind velocity

      real, allocatable, save :: prwtr(:,:,:,:)
                       ! Precipitation and accumulation for water

      real, allocatable, save :: price(:,:,:,:)
                       ! Precipitation and accumulation for ice

      real, allocatable, save :: pdia(:,:,:)
                       ! Diabatic value

      real, allocatable, save :: z0m(:,:)
                       ! Roughness length for velocity

      real, allocatable, save :: z0h(:,:)
                       ! Roughness length for scalar

      real, allocatable, save :: tund(:,:,:)
                       ! Soil and sea temperature at present

      real, allocatable, save :: tundp(:,:,:)
                       ! Soil and sea temperature at past

      real, allocatable, save :: ubr_rst(:,:,:)
                       ! ubr in restructed domain

      real, allocatable, save :: vbr_rst(:,:,:)
                       ! vbr in restructed domain

      real, allocatable, save :: pbr_rst(:,:,:)
                       ! pbr in restructed domain

      real, allocatable, save :: ptbr_rst(:,:,:)
                       ! ptbr in restructed domain

      real, allocatable, save :: qvbr_rst(:,:,:)
                       ! qvbr in restructed domain

      real, allocatable, save :: u_rst(:,:,:)
                       ! u in restructed domain

      real, allocatable, save :: up_rst(:,:,:)
                       ! up in restructed domain

      real, allocatable, save :: v_rst(:,:,:)
                       ! v in restructed domain

      real, allocatable, save :: vp_rst(:,:,:)
                       ! vp in restructed domain

      real, allocatable, save :: w_rst(:,:,:)
                       ! w in restructed domain

      real, allocatable, save :: wp_rst(:,:,:)
                       ! wp in restructed domain

      real, allocatable, save :: pp_rst(:,:,:)
                       ! pp in restructed domain

      real, allocatable, save :: ppp_rst(:,:,:)
                       ! ppp in restructed domain

      real, allocatable, save :: ptp_rst(:,:,:)
                       ! ptp in restructed domain

      real, allocatable, save :: ptpp_rst(:,:,:)
                       ! ptpp in restructed domain

      real, allocatable, save :: qv_rst(:,:,:)
                       ! qv in restructed domain

      real, allocatable, save :: qvp_rst(:,:,:)
                       ! qvp in restructed domain

      real, allocatable, save :: qwtr_rst(:,:,:,:)
                       ! qwtr in restructed domain

      real, allocatable, save :: qwtrp_rst(:,:,:,:)
                       ! qwtrp in restructed domain

      real, allocatable, save :: nwtr_rst(:,:,:,:)
                       ! nwtr in restructed domain

      real, allocatable, save :: nwtrp_rst(:,:,:,:)
                       ! nwtrp in restructed domain

      real, allocatable, save :: qice_rst(:,:,:,:)
                       ! qice in restructed domain

      real, allocatable, save :: qicep_rst(:,:,:,:)
                       ! qicep in restructed domain

      real, allocatable, save :: nice_rst(:,:,:,:)
                       ! nice in restructed domain

      real, allocatable, save :: nicep_rst(:,:,:,:)
                       ! nicep in restructed domain

      real, allocatable, save :: qcwtr_rst(:,:,:,:)
                       ! qcwtr in restructed domain

      real, allocatable, save :: qcwtrp_rst(:,:,:,:)
                       ! qcwtrp in restructed domain

      real, allocatable, save :: qcice_rst(:,:,:,:)
                       ! qcice in restructed domain

      real, allocatable, save :: qcicep_rst(:,:,:,:)
                       ! qcicep in restructed domain

      real, allocatable, save :: qasl_rst(:,:,:,:)
                       ! qasl in restructed domain

      real, allocatable, save :: qaslp_rst(:,:,:,:)
                       ! qaslp in restructed domain

      real, allocatable, save :: qt_rst(:,:,:)
                       ! qt in restructed domain

      real, allocatable, save :: qtp_rst(:,:,:)
                       ! qtp in restructed domain

      real, allocatable, save :: tke_rst(:,:,:)
                       ! tke in restructed domain

      real, allocatable, save :: tkep_rst(:,:,:)
                       ! tkep in restructed domain

      real, allocatable, save :: ucpx_rst(:,:,:)
                       ! ucpx in restructed domain

      real, allocatable, save :: ucpy_rst(:,:,:)
                       ! ucpy in restructed domain

      real, allocatable, save :: vcpx_rst(:,:,:)
                       ! vcpx in restructed domain

      real, allocatable, save :: vcpy_rst(:,:,:)
                       ! vcpy in restructed domain

      real, allocatable, save :: wcpx_rst(:,:,:)
                       ! wcpx in restructed domain

      real, allocatable, save :: wcpy_rst(:,:,:)
                       ! wcpy in restructed domain

      real, allocatable, save :: pcpx_rst(:,:,:)
                       ! pcpx in restructed domain

      real, allocatable, save :: pcpy_rst(:,:,:)
                       ! pcpy in restructed domain

      real, allocatable, save :: ptcpx_rst(:,:,:)
                       ! ptcpx in restructed domain

      real, allocatable, save :: ptcpy_rst(:,:,:)
                       ! ptcpy in restructed domain

      real, allocatable, save :: qvcpx_rst(:,:,:)
                       ! qvcpx in restructed domain

      real, allocatable, save :: qvcpy_rst(:,:,:)
                       ! qvcpy in restructed domain

      real, allocatable, save :: qwcpx_rst(:,:,:,:)
                       ! qwcpx in restructed domain

      real, allocatable, save :: qwcpy_rst(:,:,:,:)
                       ! qwcpy in restructed domain

      real, allocatable, save :: nwcpx_rst(:,:,:,:)
                       ! nwcpx in restructed domain

      real, allocatable, save :: nwcpy_rst(:,:,:,:)
                       ! nwcpy in restructed domain

      real, allocatable, save :: qicpx_rst(:,:,:,:)
                       ! qicpx in restructed domain

      real, allocatable, save :: qicpy_rst(:,:,:,:)
                       ! qicpy in restructed domain

      real, allocatable, save :: nicpx_rst(:,:,:,:)
                       ! nicpx in restructed domain

      real, allocatable, save :: nicpy_rst(:,:,:,:)
                       ! nicpy in restructed domain

      real, allocatable, save :: qcwcpx_rst(:,:,:,:)
                       ! qcwcpx in restructed domain

      real, allocatable, save :: qcwcpy_rst(:,:,:,:)
                       ! qcwcpy in restructed domain

      real, allocatable, save :: qcicpx_rst(:,:,:,:)
                       ! qcicpx in restructed domain

      real, allocatable, save :: qcicpy_rst(:,:,:,:)
                       ! qcicpy in restructed domain

      real, allocatable, save :: qacpx_rst(:,:,:,:)
                       ! qacpx in restructed domain

      real, allocatable, save :: qacpy_rst(:,:,:,:)
                       ! qacpy in restructed domain

      real, allocatable, save :: qtcpx_rst(:,:,:)
                       ! qtcpx in restructed domain

      real, allocatable, save :: qtcpy_rst(:,:,:)
                       ! qtcpy in restructed domain

      real, allocatable, save :: tkecpx_rst(:,:,:)
                       ! tkecpx in restructed domain

      real, allocatable, save :: tkecpy_rst(:,:,:)
                       ! tkecpy in restructed domain

      real, allocatable, save :: maxvl_rst(:,:,:)
                       ! maxvl in restructed domain

      real, allocatable, save :: prwtr_rst(:,:,:,:)
                       ! prwtr in restructed domain

      real, allocatable, save :: price_rst(:,:,:,:)
                       ! price in restructed domain

      real, allocatable, save :: pdia_rst(:,:,:)
                       ! pdia in restructed domain

      real, allocatable, save :: z0m_rst(:,:)
                       ! z0m in restructed domain

      real, allocatable, save :: z0h_rst(:,:)
                       ! z0h in restructed domain

      real, allocatable, save :: tund_rst(:,:,:)
                       ! tund in restructed domain

      real, allocatable, save :: tundp_rst(:,:,:)
                       ! tundp in restructed domain

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

      end module m_comrst
