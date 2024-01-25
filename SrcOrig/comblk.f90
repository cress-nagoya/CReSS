!***********************************************************************
      module m_comblk
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/05/19
!     Modification: 2003/11/05, 2004/09/25, 2004/10/12, 2004/12/17,
!                   2005/04/04, 2006/01/10, 2007/01/20, 2007/11/26,
!                   2008/01/11, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2011/01/14, 2011/09/22

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the array for bulk cold rain model.

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

      character(len=6) fmem
                       ! Control flag of memory spacing

      real, allocatable, save :: ucq(:,:,:)
                       ! Terminal velocity of cloud water

      real, allocatable, save :: urq(:,:,:)
                       ! Terminal velocity of rain water

      real, allocatable, save :: uiq(:,:,:)
                       ! Terminal velocity of cloud ice

      real, allocatable, save :: usq(:,:,:)
                       ! Terminal velocity of snow

      real, allocatable, save :: ugq(:,:,:)
                       ! Terminal velocity of graupel

      real, allocatable, save :: uhq(:,:,:)
                       ! Terminal velocity of hail

      real, allocatable, save :: ucn(:,:,:)
                       ! Terminal velocity of cloud water concentrations

      real, allocatable, save :: urn(:,:,:)
                       ! Terminal velocity of rain water concentrations

      real, allocatable, save :: uin(:,:,:)
                       ! Terminal velocity of cloud ice concentrations

      real, allocatable, save :: usn(:,:,:)
                       ! Terminal velocity of snow concentrations

      real, allocatable, save :: ugn(:,:,:)
                       ! Terminal velocity of graupel concentrations

      real, allocatable, save :: uhn(:,:,:)
                       ! Terminal velocity of hail concentrations

      real, allocatable, save :: t(:,:,:)
                       ! Air temperature

      real, allocatable, save :: tcel(:,:,:)
                       ! Ambient air temperature

      real, allocatable, save :: qvsst0(:,:,:)
                       ! Super saturation mixing ratio at melting point

      real, allocatable, save :: qvsw(:,:,:)
                       ! Saturation mixing ratio for water

      real, allocatable, save :: qvsi(:,:,:)
                       ! Saturation mixing ratio for ice

      real, allocatable, save :: lv(:,:,:)
                       ! Latent heat of evapolation

      real, allocatable, save :: ls(:,:,:)
                       ! Latent heat of sublimation

      real, allocatable, save :: lf(:,:,:)
                       ! Latent heat of fusion

      real, allocatable, save :: kp(:,:,:)
                       ! Thermal conductivity of air

      real, allocatable, save :: mu(:,:,:)
                       ! Viscosity of air

      real, allocatable, save :: dv(:,:,:)
                       ! Molecular diffusivity of water

      real, allocatable, save :: mi(:,:,:)
                       ! Mean mass of cloud ice

      real, allocatable, save :: diaqc(:,:,:)
                       ! Mean diameter of cloud water

      real, allocatable, save :: diaqr(:,:,:)
                       ! Mean diameter of rain water

      real, allocatable, save :: diaqi(:,:,:)
                       ! Mean diameter of cloud ice

      real, allocatable, save :: diaqs(:,:,:)
                       ! Mean diameter of snow

      real, allocatable, save :: diaqg(:,:,:)
                       ! Mean diameter of graupel

      real, allocatable, save :: vntr(:,:,:)
                       ! Ventilation factor for rain water

      real, allocatable, save :: vnts(:,:,:)
                       ! Ventilation factor for snow

      real, allocatable, save :: vntg(:,:,:)
                       ! Ventilation factor for graupel

      real, allocatable, save :: nuvi(:,:,:)
                       ! Nucleation rate
                       ! of deposition or sorption

      real, allocatable, save :: nuci(:,:,:)
                       ! Nucleation rate
                       ! of condensation, contact and homogeneous

      real, allocatable, save :: clcr(:,:,:)
                       ! Collection rate
                       ! between cloud water and rain water

      real, allocatable, save :: clcs(:,:,:)
                       ! Collection rate
                       ! between cloud water and snow

      real, allocatable, save :: clcg(:,:,:)
                       ! Collection rate
                       ! between cloud water and graupel

      real, allocatable, save :: clri(:,:,:)
                       ! Collection rate
                       ! between rain water and cloud ice

      real, allocatable, save :: clrs(:,:,:)
                       ! Collection rate
                       ! from rain water to snow

      real, allocatable, save :: clrg(:,:,:)
                       ! Collection rate
                       ! between rain water and graupel

      real, allocatable, save :: clir(:,:,:)
                       ! Collection rate
                       ! between rain water and cloud ice

      real, allocatable, save :: clis(:,:,:)
                       ! Collection rate
                       ! between cloud ice and snow

      real, allocatable, save :: clig(:,:,:)
                       ! Collection rate
                       ! between cloud ice and graupel

      real, allocatable, save :: clsr(:,:,:)
                       ! Collection rate
                       ! from snow to rain water

      real, allocatable, save :: clsg(:,:,:)
                       ! Collection rate
                       ! between snow and graupel

      real, allocatable, save :: clrsg(:,:,:)
                       ! Production rate of graupel
                       ! from collection rate form rain to snow

      real, allocatable, save :: clrin(:,:,:)
                       ! Collection rate for concentrations
                       ! between rain water and cloud ice

      real, allocatable, save :: clrsn(:,:,:)
                       ! Collection rate for concentrations
                       ! between rain water and snow

      real, allocatable, save :: clsrn(:,:,:)
                       ! Collection rate for concentrations
                       ! between rain water and snow

      real, allocatable, save :: clsgn(:,:,:)
                       ! Collection rate for concentrations
                       ! between snow and graupel

      real, allocatable, save :: agcn(:,:,:)
                       ! Aggregation rate for cloud water

      real, allocatable, save :: agrn(:,:,:)
                       ! Aggregation rate for rain water

      real, allocatable, save :: agin(:,:,:)
                       ! Aggregation rate for cloud ice

      real, allocatable, save :: agsn(:,:,:)
                       ! Aggregation rate for snow

      real, allocatable, save :: vdvr(:,:,:)
                       ! Evaporation rate from rain water to water vapor

      real, allocatable, save :: vdvi(:,:,:)
                       ! Deposiotn rate from water vapor to cloud ice

      real, allocatable, save :: vdvs(:,:,:)
                       ! Deposiotn rate from water vapor to snow

      real, allocatable, save :: vdvg(:,:,:)
                       ! Deposiotn rate from water vapor to graupel

      real, allocatable, save :: cncr(:,:,:)
                       ! Conversion rate from cloud water to rain water

      real, allocatable, save :: cnis(:,:,:)
                       ! Conversion rate from cloud ice to snow

      real, allocatable, save :: cnsg(:,:,:)
                       ! Conversion rate from snow to graupel

      real, allocatable, save :: cnsgn(:,:,:)
                       ! Conversion rate for concentrations
                       ! from snow to graupel

      real, allocatable, save :: spsi(:,:,:)
                       ! Secondary nucleation rate from snow

      real, allocatable, save :: spgi(:,:,:)
                       ! Secondary nucleation rate from graupel

      real, allocatable, save :: mlic(:,:,:)
                       ! Melting rate from cloud ice to cloud water

      real, allocatable, save :: mlsr(:,:,:)
                       ! Melting rate from snow to rain water

      real, allocatable, save :: mlgr(:,:,:)
                       ! Melting rate from graupel to rain water

      real, allocatable, save :: frrg(:,:,:)
                       ! Freezing rate from rain water to graupel

      real, allocatable, save :: frrgn(:,:,:)
                       ! Freezing rate for concentrations
                       ! from rain water to graupel

      real, allocatable, save :: shsr(:,:,:)
                       ! Shedding rate of liquid water from snow

      real, allocatable, save :: shgr(:,:,:)
                       ! Shedding rate of liquid water from graupel

      real, allocatable, save :: pgwet(:,:,:)
                       ! Graupel produntion rate for moist process

      real, allocatable, save :: ecs(:,:,:)
                       ! Collection efficiency of snow for cloud water

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

      end module m_comblk
