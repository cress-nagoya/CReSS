!***********************************************************************
      module m_comphy
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/05/19
!     Modification: 2003/07/15, 2003/10/31, 2003/11/05, 2004/02/01,
!                   2004/03/22, 2004/04/01, 2004/05/07, 2004/05/31,
!                   2004/07/01, 2004/08/01, 2004/08/20, 2004/09/01,
!                   2004/09/10, 2004/09/25, 2004/10/12, 2004/12/17,
!                   2005/01/07, 2005/01/14, 2005/04/04, 2005/08/05,
!                   2005/09/30, 2005/11/22, 2006/02/13, 2006/03/06,
!                   2006/08/08, 2006/09/30, 2006/12/04, 2007/04/11,
!                   2007/05/07, 2007/09/04, 2007/11/26, 2008/01/11,
!                   2008/05/02, 2008/07/01, 2008/08/25, 2008/10/10,
!                   2009/01/30, 2009/03/12, 2009/06/16, 2011/03/18,
!                   2011/04/08, 2011/06/01, 2011/09/22, 2011/11/10,
!                   2011/12/17, 2013/01/27, 2013/02/05

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the physical constants.

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

      real, parameter :: g=9.8e0
                       ! Gravity acceleration

      real, parameter :: omega=7.292e-5
                       ! Angular velocity of earth rotation

      real, parameter :: rearth=6371000.e0
                       ! Earth radius

      real, parameter :: tclat=23.44e0
                       ! Tropic of Cancer

      real, parameter :: cp=1004.e0
                       ! Specific heat of dry air at constant pressure

      real, parameter :: cv=717.e0
                       ! Specific heat of dry air at constant volume

      real, parameter :: cw=4.218e3
                       ! Specific heat of water

      real, parameter :: ci=2.106e3
                       ! Specific heat of ice

      real, parameter :: rd=287.e0
                       ! Gas constant for dry air

      real, parameter :: rv=461.e0
                       ! Gas constant for water vapor

      real, parameter :: epsva=.622e0
                       ! Morecular weight ratio of vapor / air

      real, parameter :: epsav=1.6077e0
                       ! Morecular weight ratio of air / vapor

      real, parameter :: es0=610.78e0
                       ! Saturation vapor pressure at 0 [deg C]

      real, parameter :: r0=1.225e0
                       ! Reference air density

      real, parameter :: p0=1.e5
                       ! Reference pressure

      real, parameter :: t0=273.16e0
                       ! Temperature at melting point

      real, parameter :: ttp=273.16e0
                       ! Temperature at triple point

      real, parameter :: tlow=233.16e0
                       ! Temperature at lowest super cooled point

      real, parameter :: t0cel=0.e0
                       ! Ambient air temperature at melting point

      real, parameter :: t20=293.16e0
                       ! Temperature at 20.0 [deg C]

      real, parameter :: t23=296.16e0
                       ! Temperature at 23.0 [deg C]

      real, parameter :: p20=101325.e0
                       ! Pressure for standard atmosphere at t20

      real, parameter :: lmb20=6.6e-8
                       ! Mean free path at t20 and p20

      real, parameter :: diasl=6.e-7
                       ! Diameter of aerosol particle

      real, parameter :: boltz=1.38066e-23
                       ! Boltzmann constant

      real, parameter :: sc=.6e0
                       ! Schmidt number

      real, parameter :: lv0=2.50078e6
                       ! Latent heat of evaporation

      real, parameter :: lf0=3.34e5
                       ! Latent heat of fusion

      real, parameter :: kp0=10.79e0
                       ! Standard value of thermal conductivity of air

      real, parameter :: kpa=2.4e-1
                       ! Thermal conductivity of aerosol

      real, parameter :: nu0=1.328e-5
                       ! Kinetic viscosity of air at 0 [deg C]

      real, parameter :: da0=1.87e-5
                       ! Molecular diffusivity of air at 0 [deg C]

      real, parameter :: dv0=2.23e-5
                       ! Molecular diffusivity of water at 0 [deg C]

      real, parameter :: xl=.25e0
                       ! Dispersion of fall velocity spectrum
                       ! of ice crystals

      real, parameter :: nc0=3.e12
                       ! Parameter of cloud water size distribution

      real, parameter :: nr0=8.e6
                       ! Parameter of rain water size distribution

      real, parameter :: ni0=3.e11
                       ! Parameter of cloud ice size distribution

      real, parameter :: ns0=1.8e6
                       ! Parameter of snow size distribution

      real, parameter :: ng0=1.1e6
                       ! Parameter of graupel size distribution

      real, parameter :: nh0=1.1e6
                       ! Parameter of hail size distribution

      real, parameter :: auc=2.98e7
                       ! Constant in empirical formula for cloud water

      real, parameter :: aur=842.e0
                       ! Constant in empirical formula for rain water

      real, parameter :: aui=770.e0
                       ! Constant in empirical formula for cloud ice

      real, parameter :: aus=27.e0
                       ! Constant in empirical formula for snow

      real, parameter :: aug=132.e0
                       ! Constant in empirical formula for graupel

      real, parameter :: auh=132.e0
                       ! Constant in empirical formula for hail

      real, parameter :: buc=2.e0
                       ! Constant in empirical formula for cloud water

      real, parameter :: bur=.8e0
                       ! Constant in empirical formula for rain water

      real, parameter :: bui=1.e0
                       ! Constant in empirical formula for cloud ice

      real, parameter :: bus=.5e0
                       ! Constant in empirical formula for snow

      real, parameter :: bug=.64e0
                       ! Constant in empirical formula for graupel

      real, parameter :: buh=.64e0
                       ! Constant in empirical formula for hail

      real, parameter :: guc=1.e0
                       ! Constant in empirical formula for cloud water

      real, parameter :: gur=.5e0
                       ! Constant in empirical formula for rain water

      real, parameter :: gui=.33e0
                       ! Constant in empirical formula for cloud ice

      real, parameter :: gus=.5e0
                       ! Constant in empirical formula for snow

      real, parameter :: gug=.5e0
                       ! Constant in empirical formula for graupel

      real, parameter :: guh=.5e0
                       ! Constant in empirical formula for hail

      real, parameter :: rhow=1.e3
                       ! Density of water

      real, parameter :: rhoi=2.e2
                       ! Density of cloud ice

      real, parameter :: rhos=1.e2
                       ! Density of snow

      real, parameter :: rhog=4.e2
                       ! Density of graupel

      real, parameter :: rhoh=4.e2
                       ! Density of hail

      real, parameter :: ers=1.e0
                       ! Collection efficiency of snow
                       ! for rain water

      real, parameter :: erg=1.e0
                       ! Collection efficiency of graupel
                       ! for rain water

      real, parameter :: eir=1.e0
                       ! Collection efficiency of rain water
                       ! for cloud ice

      real, parameter :: eis=.1e0
                       ! Collection efficiency of snow
                       ! for cloud ice

      real, parameter :: eig=.1e0
                       ! Collection efficiency of graupel
                       ! for cloud ice

      real, parameter :: esr=1.e0
                       ! Collection efficiency of rain water
                       ! for snow

      real, parameter :: esg=.001e0
                       ! Collection efficiency of graupel
                       ! for snow

      real, parameter :: ecc=.55e0
                       ! Collection efficiency among cloud water

      real, parameter :: eii=.1e0
                       ! Collection efficiency among cloud ice

      real, parameter :: ess=.1e0
                       ! Collection efficiency among snow

      real, parameter :: mc0=.5e-12
                       ! Mass of minimum cloud water

      real, parameter :: mr0=2.7e-10
                       ! Mass of minimum rain water

      real, parameter :: mi0=1.e-12
                       ! Mass of minimum cloud ice

      real, parameter :: ms0=1.76e-10
                       ! Mass of minimum snow

      real, parameter :: mg0=7.e-10
                       ! Mass of minimum graupel

      real, parameter :: mh0=7.e-10
                       ! Mass of minimum hail

      real, parameter :: mcmax=2.7e-10
                       ! Mass of maximum cloud water

      real, parameter :: mrmax=2.68e-4
                       ! Mass of maximum rain water

      real, parameter :: mimax=1.4e-10
                       ! Mass of maximum cloud ice

      real, parameter :: msmax=6.e-4
                       ! Mass of maximum snow

      real, parameter :: mgmax=8.e-3
                       ! Mass of maximum graupel

      real, parameter :: mhmax=8.e-3
                       ! Mass of maximum hail

      real, parameter :: nclcst=1.e8
                       ! Constant concentrations of cloud

      real, parameter :: diaqcm=2.e-5
                       ! Critical mean diameter of cloud water

      real, parameter :: qccrit=1.7e-5
                       ! Critical mean total water mixing ratio

      real, parameter :: dc0=.036e0
                       ! Deposition coefficient of water bin

      real, parameter :: eta0=1.818e-4
                       ! Constant of terminal velocity of water bin

      real, parameter :: wsten=75.64e0
                       ! Surface tension of water [dyn/cm^2]

      real, parameter :: rbwmin=1.e-4
                       ! Radius of minimum water bin [cm]

      real, parameter :: rbwmax=.35e0
                       ! Radius of maximum water bin for
                       ! terminal fall velocity [cm]

      real, parameter :: dbmw0=1.e-10
                       ! Minimum differential
                       ! between adjacent water bins [g]

      real, parameter :: ckm=.1e0
                       ! Constant of eddy viscosity

      real, parameter :: ckmin=1.e-6
                       ! Coefficient of lower limit of eddy viscosity

      real, parameter :: csnum=0.21e0
                       ! Smagorinsky constant

      real, parameter :: prnum=0.33e0
                       ! Constant turbulent Prandtl number

      real, parameter :: prnumg=.74e0
                       ! Turbulent Prandtl number on soil

      real, parameter :: prnumw=1.e0
                       ! Turbulent Prandtl number on sea

      real, parameter :: kappa=.4e0
                       ! Karman constant

      real, parameter :: wkappa=.35e0
                       ! Karman constant on sea surface

      real, parameter :: sigma=5.67e-8
                       ! Stefan Boltzman constant

      real, parameter :: sun0=1367.e0
                       ! Sun light constant

      real, parameter :: el=.95e0
                       ! Emissivity on surface

      real, parameter :: sealbe=.05e0
                       ! Minimum albedo on sea surface

      real, parameter :: icalbe=.5e0
                       ! Albedo on ice surface

      real, parameter :: snalbe=.6e0
                       ! Albedo on snow surface

      real, parameter :: sebeta=1.e0
                       ! Evapotranspiration efficiency on sea surface

      real, parameter :: icbeta=1.e0
                       ! Evapotranspiration efficiency on ice surface

      real, parameter :: snbeta=1.e0
                       ! Evapotranspiration efficiency on snow surface

!ORIG real, parameter :: sez0m=5.e-5
      real, parameter :: sez0m=2.e-4
                       ! Constant roughness length for velocity
                       ! on sea surface

      real, parameter :: icz0m=5.e-4
                       ! Constant roughness length for velocity
                       ! on ice surface

      real, parameter :: snz0m=1.e-3
                       ! Constant roughness length for velocity
                       ! on snow surface

!ORIG real, parameter :: sez0h=1.e-5
      real, parameter :: sez0h=2.e-4
                       ! Constant roughness length for scalar
                       ! on sea surface

      real, parameter :: icz0h=1.e-4
                       ! Constant roughness length for scalar
                       ! on ice surface

      real, parameter :: snz0h=2.e-4
                       ! Constant roughness length for scalar
                       ! on snow surface

      real, parameter :: secap=4.18e6
                       ! Thermal capacity of sea water

!ORIG real, parameter :: senuu=2.e-3
      real, parameter :: senuu=1.e-3
                       ! Thermal diffusivity of sea water

      real, parameter :: rmg=15.e0
                       ! Constant of bulk coefficients for soil surface

      real, parameter :: rhg=9.e0
                       ! Constant of bulk coefficients for soil surface

      real, parameter :: rms=16.e0
                       ! Constant of bulk coefficients for sea surface

      real, parameter :: rhs=16.e0
                       ! Constant of bulk coefficients for sea surface

      real, parameter :: zarad=100.e0
                       ! Reference height for downward radiation

      real, parameter :: zamax=600.e0
                       ! Maximum z physical coordinates
                       ! of planetary boundary layer

      real, parameter :: vamin=.4e0
                       ! Minimum magnitude of velocity to avoid
                       ! to be divided by 0.0

!ORIG real, parameter :: z0min=1.e-5
      real, parameter :: z0min=1.5e-5
                       ! Minimum roughness length on sea surface

!ORIG real, parameter :: rchmin=-100.e0
      real, parameter :: rchmin=-10.e0
                       ! Minimum bulk Richardson number

      real, parameter :: prmin=1.389e-7
                       ! Lower limit of precipitation

      real, parameter :: divndc=0.05e0
                       ! No dimensional divergence damping coefficient

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

      end module m_comphy
