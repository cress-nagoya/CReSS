!***********************************************************************
      module m_defname
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/05/10,
!                   1999/05/20, 1999/07/05, 1999/07/23, 1999/08/03,
!                   1999/08/18, 1999/09/16, 1999/09/30, 1999/11/01,
!                   1999/11/19, 1999/11/24, 1999/12/17, 2000/01/17,
!                   2000/04/18, 2000/06/01, 2000/07/05, 2000/08/10,
!                   2000/12/18, 2001/01/15, 2001/03/13, 2001/04/15,
!                   2001/06/06, 2001/07/13, 2001/08/07, 2001/10/18,
!                   2001/11/20, 2002/02/05, 2002/06/18, 2002/07/03,
!                   2002/08/15, 2002/08/27, 2002/09/02, 2002/09/09,
!                   2002/10/15, 2002/10/31, 2002/11/11, 2003/03/28,
!                   2003/04/30, 2003/05/19, 2003/07/15, 2003/08/08,
!                   2003/09/01, 2003/10/10, 2003/11/05, 2003/12/12,
!                   2004/01/09, 2004/02/01, 2004/03/05, 2004/04/15,
!                   2004/05/07, 2004/05/31, 2004/08/01, 2004/08/20,
!                   2004/09/01, 2004/09/10, 2004/10/12, 2004/12/17,
!                   2005/01/14, 2005/02/10, 2005/08/05, 2005/12/13,
!                   2006/04/03, 2006/07/21, 2006/09/21, 2006/09/30,
!                   2006/11/27, 2006/12/04, 2007/01/05, 2007/01/20,
!                   2007/03/23, 2007/04/24, 2007/05/21, 2007/06/21,
!                   2008/01/11, 2008/04/17, 2008/05/02, 2008/07/01,
!                   2008/08/25, 2008/10/10, 2008/12/11, 2009/01/05,
!                   2009/01/30, 2009/02/27, 2011/05/16, 2011/08/09,
!                   2011/08/18, 2011/09/22, 2011/11/10, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the namelist variables and blocks.

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

! Module variables (namelist variables)

      character(len=108) exprim
                       ! Optional run name

      character(len=108) crsdir
                       ! User specified directory for CReSS files

      character(len=108) datdir
                       ! User specified directory for external data

      character(len=22) sfcast
                       ! Forecast start date
                       ! with Gregorian calendar, "yyyy/mm/dd hh:mm"

      character(len=108) lbcvar
                       ! Control flag of
                       ! lateral boundary forced variables

      character(len=108) gpvvar
                       ! Control flag of
                       ! input GPV data variables

      character(len=108) nggvar
                       ! Control flag of
                       ! analysis nudged variables to GPV

      character(len=108) exbvar
                       ! Control flag of
                       ! external boundary forced variables

      character(len=108) lspvar
                       ! Control flag of
                       ! lateral sponge damped variables

      character(len=108) vspvar
                       ! Control flag of
                       ! vertical sponge damped variables

      character(len=108) rdrvar
                       ! Control flag of
                       ! input radar data variables

      character(len=108) ngrvar
                       ! Control flag of
                       ! analysis nudged variables to radar

      character(len=108) sfcdat
                       ! Control flag of input surface data type

      character(len=108) prvres
                       ! Restart file name without extension
                       ! of previous running

      character(len=108) sndtyp
                       ! Control flag of sounding data type

      character(len=108) dmpvar
                       ! Control flag of dumped variables

      character(len=108) mxnvar
                       ! Control flag of maximum and mininum output

      character(len=108) datype_gpv
                       ! Control flag of GPV data type

      character(len=108) etrvar_gpv
                       ! Control flag of extrapolating method
                       ! for GPV data

      character(len=108) datype_rdr
                       ! Control flag of radar data type

      character(len=108) fltyp_uni
                       ! Control flag of processed file type

      integer savmem   ! Option for memory saving

      integer numpe    ! Total number of processor elements

      integer xgroup   ! Number of group domain
                       ! in entire domain in x direction

      integer ygroup   ! Number of group domain
                       ! in entire domain in y direction

      integer xsub     ! Number of sub domain
                       ! in group domain in x direction

      integer ysub     ! Number of sub domain
                       ! in group domain in y direction

      integer xdim     ! Model dimension in x direction
      integer ydim     ! Model dimension in y direction
      integer zdim     ! Model dimension in z direction

      integer mpopt    ! Option for map projection
      integer nspol    ! Option for projected region

      integer sthopt   ! Option for vertical grid stretching

      integer trnopt   ! Option for terrain height setting

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions
      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

      integer bbc      ! Option for bottom boundary condition
      integer tbc      ! Option for top boundary condition

      integer gsmopt   ! Option for interpolated GPV data smoothing
      integer gsmcnt   ! Iteration count for GPV data smoothing

      integer nggopt   ! Option for analysis nudging to GPV

      integer exbopt   ! Option for external boundary forcing

      integer exbwid   ! Corresponded thickness
                       ! between model and data height

      integer lspopt   ! Option for lateral sponge damping

      integer wdnews   ! Lateral sponge damping thickness

      integer wdnorm   ! Lateral sponge damping thickness
                       ! for u and v in normal

      integer vspopt   ! Option for vertical sponge damping

      integer ngropt   ! Option for analysis nudging to radar

      integer sfcopt   ! Option for surface physics

      integer levpbl   ! Number of planetary boundary layer
      integer levund   ! Number of soil and sea layers

      integer lnduse   ! User specified land use category

      integer dstopt   ! Option for sea ice distribution

      integer iniopt   ! Option for model initialization

      integer snddim   ! Sounding data dimension

      integer masopt   ! Option for masscon model

      integer movopt   ! Option for grid moving

      integer pt0opt   ! Option for
                       ! initial potential temperature perturbation

      integer pt0num   ! Number of buble shaped
                       ! initial potential temperature perturbation

      integer gwmopt   ! Option for gravity wave mode integration

      integer impopt   ! Option for vertical implicit time integration

      integer advopt   ! Option for advection scheme

      integer smtopt   ! Option for numerical smoothing

      integer mfcopt   ! Option for map scale factor

      integer coropt   ! Option for Coriolis force

      integer crvopt   ! Option for earth curvature

      integer buyopt   ! Option for buoyancy calculation

      integer diaopt   ! Option for diabatic calculation

      integer divopt   ! Option for divergence damping

      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes

      integer ncbinw   ! Number of categories for warm bin

      integer qcgopt   ! Option for charging distribution

      integer aslopt   ! Option for aerosol processes

      integer trkopt   ! Option for mixing ratio tracking

      integer qt0opt   ! Option for initial tracer location
      integer qt0num   ! Number of buble shaped initial tracer

      integer tubopt   ! Option for turbulent mixing
      integer isoopt   ! Option for grid shape

      integer dmpfmt   ! Option for dumped file format
      integer dmplev   ! Option for z coordinates of dumped variable
      integer dmpmon   ! Option for dumped monitor variables output

      integer resopt   ! Option for restart output

      integer mxnopt   ! Option for maxmum and minimum output

      integer xdim_gpv ! GPV data dimension in x direction
      integer ydim_gpv ! GPV data dimension in y direction
      integer zdim_gpv ! GPV data dimension in z direction

      integer mpopt_gpv
                       ! Option for map projection for GPV data

      integer nspol_gpv
                       ! Option for projected region for GPV data

      integer intopt_gpv
                       ! Option for GPV data interpolating

      integer rotopt_gpv
                       ! Option for rotation of wind direction
                       ! in GPV data

      integer refsfc_gpv
                       ! Option for surface data reference
                       ! in interpolating for GPV

      integer xdim_asl ! Aerosol data dimension in x direction
      integer ydim_asl ! Aerosol data dimension in y direction
      integer zdim_asl ! Aerosol data dimension in z direction

      integer mpopt_asl
                       ! Option for map projection for aerosol data

      integer nspol_asl
                       ! Option for projected region for aerosol data

      integer intopt_asl
                       ! Option for aerosol data interpolating

      integer xdim_rdr ! Radar data dimension in x direction
      integer ydim_rdr ! Radar data dimension in y direction
      integer zdim_rdr ! Radar data dimension in z direction

      integer mpopt_rdr
                       ! Option for map projection for radar data

      integer nspol_rdr
                       ! Option for projected region for radar data

      integer rotopt_rdr
                       ! Option for rotation of wind direction
                       ! in radar data

      integer xdim_trn ! Terrain data dimension in x direction
      integer ydim_trn ! Terrain data dimension in y direction

      integer mpopt_trn
                       ! Option for map projection for terrain data

      integer nspol_trn
                       ! Option for projected region for terrain data

      integer intopt_trn
                       ! Option for terrain data interpolating

      integer xdim_lnd ! Land use data dimension in x direction
      integer ydim_lnd ! Land use data dimension in y direction

      integer mpopt_lnd
                       ! Option for map projection for land use data

      integer nspol_lnd
                       ! Option for projected region for land use data

      integer intopt_lnd
                       ! Option for linear interpolating
                       ! for land use data

      integer numctg_lnd
                       ! Number of land use data categories

      integer lnduse_lnd(1:100)
                       ! Land use data

      integer xdim_sst ! Sea surface temperature data dimension
                       ! in x direction

      integer ydim_sst ! Sea surface temperature data dimension
                       ! in y direction

      integer mpopt_sst
                       ! Option for map projection
                       ! for sea surface temperature data

      integer nspol_sst
                       ! Option for projected region
                       ! for sea surface temperature data

      integer xdim_ice ! Sea ice distribution data dimension
                       ! in x direction

      integer ydim_ice ! Sea ice distribution data dimension
                       ! in y direction

      integer mpopt_ice
                       ! Option for map projection
                       ! for sea ice distribution data

      integer nspol_ice
                       ! Option for projected region
                       ! for sea ice distribution data

      integer rmopt_uni
                       ! Option for original dumped files removing

      integer uniopt_uni
                       ! Option for uniting process

      integer ugroup_uni
                       ! User specified and united group number
                       ! in entire domain

      integer xsub_rst ! Number of sub domain in group domain
                       ! in x direction for restructed files

      integer ysub_rst ! Number of sub domain in group domain
                       ! in y direction for restructed files

      integer rmopt_rst
                       ! Option for original restart files removing

      real tlat1       ! True latitude
      real tlat2       ! True latitude
      real tlon        ! True longitude

      real disr        ! Distance from center of circular cylinder
                       ! to origin of calculation domain

      real dx          ! Grid distance in x direction
      real dy          ! Grid distance in y direction
      real dz          ! Grid distance in z direction

      real ulat        ! User specified latitude
      real ulon        ! User specified longitude

      real riu         ! User specified real index in x direction
      real rju         ! User specified real index in y direction

      real zsfc        ! Sea surface terrain height

      real zflat       ! Lowest flat level

      real dzmin       ! Minimum dz in lowest layer

      real layer1      ! Lowest stretching level
      real layer2      ! Highest stretching level

      real mnthgh(1:2) ! Flat or bell shaped mountain height
                       ! and base level height

      real mntwx       ! Bell shaped mountain width in x direction
      real mntwy       ! Bell shaped mountain width in y direction

      real mntcx       ! Center in x coordinates of
                       ! bell shaped mountain

      real mntcy       ! Center in y coordinates of
                       ! bell shaped mountain

      real stime       ! Forecast start time
      real etime       ! Forecast stop time

      real lbnews      ! Lateral damping coefficient
                       ! of open boundary

      real lbnorm      ! Lateral damping coefficient
                       ! of open boundary for u and v in normal

      real gwave       ! Fastest gravity wave speed

      real gpvitv      ! Time interval of GPV data file
      real aslitv      ! Time interval of aerosol data file

      real gsmcoe      ! Interpolated GPV data smoothing coefficient

      real nggcoe      ! Analysis nudging coefficient to GPV
      real nggdlt      ! Time interval of analysis nudging to GPV
      real nggstr      ! Analysis nudging start time to GPV
      real nggend      ! Analysis nudging end time to GPV
      real nggc20      ! Start time to decrease nggcoe to 0 to GPV

      real exnews      ! External boundary damping coefficient

      real exnorm      ! External boundary damping coefficient
                       ! for u and v in normal

      real lspsmt      ! Lateral sponge smoothing coefficient

      real lsnews      ! Lateral sponge damping coefficient

      real lsnorm      ! Lateral sponge damping coefficient
                       ! for u and v in normal

      real vspgpv      ! Vertical sponge damping coefficient
                       ! for external data

      real vspbar      ! Vertical sponge damping coefficient
                       ! for base state or 0

      real botgpv      ! Lowest height of vertical sponge damping
                       ! for external data

      real botbar      ! Lowest height of vertical sponge damping
                       ! for base state or 0

      real rdritv      ! Time interval of radar data file

      real ngrcoe      ! Analysis nudging coefficient to radar
      real ngrdlt      ! Time interval of analysis nudging to radar
      real ngrstr      ! Analysis nudging start time to radar
      real ngrend      ! Analysis nudging end time to radar
      real ngrc20      ! Start time to decrease ngrcoe to 0 to radar
      real ngraff      ! Analysis nudging affected time to radar

      real dtgrd       ! Time interval of soil temperature calculation

      real dzgrd       ! Grid distance in soil layers in z direction
      real dzsea       ! Grid distance in sea layers in z direction

      real tgdeep      ! Constant soil temperature in deepest layer

      real gralbe      ! Albedo on soil surface
      real grbeta      ! Evapotranspiration efficiency on soil surface
      real grz0m       ! Roughness length for velocity on soil surface
      real grz0h       ! Roughness length for scalar on soil surface
      real grcap       ! Thermal capacity of soil
      real grnuu       ! Thermal diffusivity of soil

      real sstcst      ! Constant sea surface temperature

      real sstitv      ! Time interval of
                       ! sea surface temperature data file

      real zsnd0       ! Reference height of sounding data
      real psnd0       ! Reference pressure of sounding data

      real maseps      ! Value of convergence of iteration for masscon

      real alpha1      ! Weighting coefficinent
                       ! for fitting to mass consistent equation

      real alpha2      ! Weighting coefficinent
                       ! for fitting to mass consistent equation

      real umove       ! x components of grid moving speed
      real vmove       ! y components of grid moving speed

      real ptp0        ! Magnitude or amplitude
                       ! of potential temperature perturbation

      real pt0rx       ! Potential temperature perturbation radius
                       ! or length in x direction

      real pt0ry       ! Potential temperature perturbation radius
                       ! or length in y direction

      real pt0rz       ! Potential temperature perturbation radius
                       ! or length in z direction

      real pt0cx       ! Center or origin in x coordinates
                       ! of potential temperature perturbation

      real pt0cy       ! Center or origin in y coordinates
                       ! of potential temperature perturbation

      real pt0cz       ! Center or origin in z coordinates
                       ! of potential temperature perturbation

      real pt0ds       ! Distance between each perturbation buble

      real dtbig       ! Large time steps interval
      real dtsml       ! Small time steps interval

      real gsdeps      ! Value of convergence of iteration
                       ! for Gauss-Seidel method

      real weicoe      ! Weighting coefficient for impilicit method

      real filcoe      ! Coefficient of Asselin time filter

      real dtvcul      ! Time steps interval
                       ! of vertical Cubic Lagrange advection

      real smhcoe      ! Horizontal smoothig coefficient
      real smvcoe      ! Vertical smoothig coefficient

      real nlhcoe      ! Horizontal non linear smoothig coefficient
      real nlvcoe      ! Vertical non linear smoothig coefficient

      real thresq      ! Minimum threshold value of mixing ratio

      real dtcmph(1:6) ! Time interval of cloud micro physics

      real bbinw       ! Base of exponential function
                       ! in mass components

      real sbinw       ! Exponent of exponential function
                       ! in mass components

      real eledlt      ! Time interval of electric field calculation

      real qt0         ! Magnitude or amplitude of tracer

      real qt0rx       ! Tracer radius or length in x direction
      real qt0ry       ! Tracer radius or length in y direction
      real qt0rz       ! Tracer radius or length in z direction

      real qt0cx       ! Center or origin in x coordinates of tracer
      real qt0cy       ! Center or origin in y coordinates of tracer
      real qt0cz       ! Center or origin in z coordinates of tracer

      real qt0ds       ! Distance between each tracer buble

      real qtdt        ! Emitted tracer intensity

      real qtstr       ! User specified emitting start time
      real qtend       ! User specified emitting end time

      real dmpitv      ! Time interval of dumped files

      real monitv      ! Time interval of dumped files
                       ! for monitor variables

      real resitv      ! Time interval of restart files

      real mxnitv      ! Time interval of maximum and minimum output
                       ! to standard i/o

      real tlat1_gpv   ! True latitude of GPV data
      real tlat2_gpv   ! True latitude of GPV data
      real tlon_gpv    ! True longitude of GPV data

      real dx_gpv      ! Grid distance in x direction of GPV data
      real dy_gpv      ! Grid distance in y direction of GPV data

      real ulat_gpv    ! User specified latitude on GPV data
      real ulon_gpv    ! User specified longitude on GPV data

      real riu_gpv     ! User specified real index on GPV data
                       ! in x direction

      real rju_gpv     ! User specified real index on GPV data
                       ! in y direction

      real tlat1_asl   ! True latitude of aerosol data
      real tlat2_asl   ! True latitude of aerosol data
      real tlon_asl    ! True longitude of aerosol data

      real dx_asl      ! Grid distance in x direction of aerosol data
      real dy_asl      ! Grid distance in y direction of aerosol data

      real ulat_asl    ! User specified latitude on aerosol data
      real ulon_asl    ! User specified longitude on aerosol data

      real riu_asl     ! User specified real index on aerosol data
                       ! in x direction

      real rju_asl     ! User specified real index on aerosol data
                       ! in y direction

      real tlat1_rdr   ! True latitude of radar data
      real tlat2_rdr   ! True latitude of radar data
      real tlon_rdr    ! True longitude of radar data

      real dx_rdr      ! Grid distance in x direction of radar data
      real dy_rdr      ! Grid distance in y direction of radar data

      real ulat_rdr    ! User specified latitude on radar data
      real ulon_rdr    ! User specified longitude on radar data

      real riu_rdr     ! User specified real index on radar data
                       ! in x direction

      real rju_rdr     ! User specified real index on radar data
                       ! in y direction

      real rdrcoe_rdr  ! Coeffient of converter of [dBZe] to [kg/m^3]
      real rdrexp_rdr  ! Exponent of converter of [dBZe] to [kg/m^3]

      real tlat1_trn   ! True latitude of terrain data
      real tlat2_trn   ! True latitude of terrain data
      real tlon_trn    ! True longitude of terrain data

      real dx_trn      ! Grid distance in x direction of terrain data
      real dy_trn      ! Grid distance in y direction of terrain data

      real ulat_trn    ! User specified latitude on terrain data
      real ulon_trn    ! User specified longitude on terrain data

      real riu_trn     ! User specified real index on terrain data
                       ! in x direction

      real rju_trn     ! User specified real index on terrain data
                       ! in y direction

      real tlat1_lnd   ! True latitude of land use data
      real tlat2_lnd   ! True latitude of land use data
      real tlon_lnd    ! True longitude of land use data

      real dx_lnd      ! Grid distance in x direction of land use data
      real dy_lnd      ! Grid distance in y direction of land use data

      real ulat_lnd    ! User specified latitude on land use data
      real ulon_lnd    ! User specified longitude on land use data

      real riu_lnd     ! User specified real index on land use data
                       ! in x direction

      real rju_lnd     ! User specified real index on land use data
                       ! in y direction

      real albe_lnd(1:100)
                       ! Albedo data

      real beta_lnd(1:100)
                       ! Evapotranspiration efficiency data

      real z0m_lnd(1:100)
                       ! Roughness length data for velocity

      real z0h_lnd(1:100)
                       ! Roughness length data for scalar

      real cap_lnd(1:100)
                       ! Thermal capacity data

      real nuu_lnd(1:100)
                       ! Thermal diffusivity data

      real tlat1_sst   ! True latitude of sea surface temperature data
      real tlat2_sst   ! True latitude of sea surface temperature data
      real tlon_sst    ! True longitude of sea surface temperature data

      real dx_sst      ! Grid distance in x direction of
                       ! sea surface temperature data

      real dy_sst      ! Grid distance in y direction of
                       ! sea surface temperature data

      real ulat_sst    ! User specified latitude on
                       ! sea surface temperature data

      real ulon_sst    ! User specified longitude on
                       ! sea surface temperature data

      real riu_sst     ! User specified real index on
                       ! sea surface temperature data in x direction

      real rju_sst     ! User specified real index on
                       ! sea surface temperature data in y direction

      real tlat1_ice   ! True latitude of sea ice distribution data
      real tlat2_ice   ! True latitude of sea ice distribution data
      real tlon_ice    ! True longitude of sea ice distribution data

      real dx_ice      ! Grid distance in x direction of
                       ! sea ice distribution data

      real dy_ice      ! Grid distance in y direction of
                       ! sea ice distribution data

      real ulat_ice    ! User specified latitude on
                       ! sea ice distribution data

      real ulon_ice    ! User specified longitude on
                       ! sea ice distribution data

      real riu_ice     ! User specified real index on
                       ! sea ice distribution data in x direction

      real rju_ice     ! User specified real index on
                       ! sea ice distribution data in y direction

      real flitv_uni   ! Time interval of original dumped file
      real flitv_rst   ! Time interval of original restart file

! Module variables (subordinate variables)

      character(len=12) idate
                       ! Forecast start date
                       ! with Gregorian calendar, yyyymmddhhmm

      integer wlngth   ! Word length of direct access file

      integer ncexp    ! Number of character of exprim
      integer nccrs    ! Number of character of crsdir
      integer ncdat    ! Number of character of datdir
      integer ncprv    ! Number of character of prvres

      integer iwest    ! Added index on west boundary
      integer ieast    ! Subtracted index on east boundary
      integer jsouth   ! Added index on south boundary
      integer jnorth   ! Subtracted index on north boundary

      real dxiv        ! Inverse of dx
      real dyiv        ! Inverse of dy
      real dziv        ! Inverse of dz

      real dxiv_gpv    ! Inverse of dx_gpv
      real dyiv_gpv    ! Inverse of dy_gpv

      real dxiv_asl    ! Inverse of dx_asl
      real dyiv_asl    ! Inverse of dy_asl

      real dxiv_rdr    ! Inverse of dx_rdr
      real dyiv_rdr    ! Inverse of dy_rdr

      real dxiv_trn    ! Inverse of dx_trn
      real dyiv_trn    ! Inverse of dy_trn

      real dxiv_lnd    ! Inverse of dx_lnd
      real dyiv_lnd    ! Inverse of dy_lnd

      real dxiv_sst    ! Inverse of dx_sst
      real dyiv_sst    ! Inverse of dy_sst

      real dxiv_ice    ! Inverse of dx_ice
      real dyiv_ice    ! Inverse of dy_ice

! Module procedure

!     none

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

!     none

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Namelist statements

      namelist /sysdep/ savmem

      namelist /runame/ exprim

      namelist /drname/ crsdir,datdir

      namelist /dimset/ numpe,xgroup,ygroup,xsub,ysub,xdim,ydim,zdim

      namelist /project/ mpopt,nspol,tlat1,tlat2,tlon,disr

      namelist /gridset/ dx,dy,dz,ulat,ulon,riu,rju

      namelist /gridsth/ zsfc,zflat,sthopt,dzmin,layer1,layer2

      namelist /terrain/ trnopt,mnthgh,mntwx,mntwy,mntcx,mntcy

      namelist /flength/ sfcast,stime,etime

      namelist /boundry/ wbc,ebc,sbc,nbc,bbc,tbc,lbcvar,lbnews,lbnorm,  &
     &                   gwave

      namelist /gpvpram/ gpvvar,gpvitv,aslitv,                          &
     &                   gsmopt,gsmcnt,gsmcoe,nggopt,nggvar,nggcoe,     &
     &                   nggdlt,nggstr,nggend,nggc20,exbopt,exbvar,     &
     &                   exnews,exnorm,exbwid,lspopt,lspvar,lspsmt,     &
     &                   lsnews,lsnorm,wdnews,wdnorm,vspopt,vspvar,     &
     &                   vspgpv,vspbar,botgpv,botbar

      namelist /rdrpram/ rdrvar,rdritv,ngropt,ngrvar,ngrcoe,ngrdlt,     &
     &                   ngrstr,ngrend,ngrc20,ngraff

      namelist /sfcphys/ sfcdat,sfcopt,levpbl,levund,dtgrd,dzgrd,dzsea, &
     &                   tgdeep,prvres,lnduse,gralbe,grbeta,grz0m,grz0h,&
     &                   grcap,grnuu,sstcst,sstitv,dstopt

      namelist /initype/ iniopt,snddim,sndtyp,zsnd0,psnd0,              &
     &                   masopt,maseps,alpha1,alpha2

      namelist /gridmove/ movopt,umove,vmove

      namelist /ptinicon/ pt0opt,pt0num,ptp0,pt0rx,pt0ry,pt0rz,pt0cx,   &
     &                    pt0cy,pt0cz,pt0ds

      namelist /integrat/ dtbig,dtsml,gwmopt,impopt,advopt,             &
     &                    gsdeps,weicoe,filcoe,dtvcul

      namelist /smoother/ smtopt,smhcoe,smvcoe,nlhcoe,nlvcoe

      namelist /mapfcter/ mfcopt

      namelist /coriolis/ coropt

      namelist /earthcrv/ crvopt

      namelist /buoyancy/ buyopt

      namelist /diabatic/ diaopt

      namelist /ddamping/ divopt

      namelist /cloudphy/ cphopt,haiopt,thresq,                         &
     &                    dtcmph,ncbinw,bbinw,sbinw,qcgopt,eledlt

      namelist /asolproc/ aslopt

      namelist /mixtrace/ trkopt,qt0opt,qt0num,qt0,qt0rx,qt0ry,qt0rz,   &
     &                    qt0cx,qt0cy,qt0cz,qt0ds,qtdt,qtstr,qtend

      namelist /turbulen/ tubopt,isoopt

      namelist /outfomat/ dmpfmt,dmplev,dmpmon,dmpvar,dmpitv,monitv,    &
     &                    resopt,resitv,mxnopt,mxnvar,mxnitv

      namelist /project_gpv/ mpopt_gpv,nspol_gpv,                       &
     &                       tlat1_gpv,tlat2_gpv,tlon_gpv

      namelist /gridset_gpv/ xdim_gpv,ydim_gpv,zdim_gpv,dx_gpv,dy_gpv,  &
     &                       ulat_gpv,ulon_gpv,riu_gpv,rju_gpv

      namelist /datconf_gpv/ intopt_gpv,rotopt_gpv,datype_gpv,          &
     &                       refsfc_gpv,etrvar_gpv

      namelist /project_asl/ mpopt_asl,nspol_asl,                       &
     &                       tlat1_asl,tlat2_asl,tlon_asl

      namelist /gridset_asl/ xdim_asl,ydim_asl,zdim_asl,dx_asl,dy_asl,  &
     &                       ulat_asl,ulon_asl,riu_asl,rju_asl

      namelist /datconf_asl/ intopt_asl

      namelist /project_rdr/ mpopt_rdr,nspol_rdr,                       &
     &                       tlat1_rdr,tlat2_rdr,tlon_rdr

      namelist /gridset_rdr/ xdim_rdr,ydim_rdr,zdim_rdr,dx_rdr,dy_rdr,  &
     &                       ulat_rdr,ulon_rdr,riu_rdr,rju_rdr

      namelist /datconf_rdr/ rotopt_rdr,datype_rdr,rdrcoe_rdr,rdrexp_rdr

      namelist /project_trn/ mpopt_trn,nspol_trn,                       &
     &                       tlat1_trn,tlat2_trn,tlon_trn

      namelist /gridset_trn/ xdim_trn,ydim_trn,dx_trn,dy_trn,           &
     &                       ulat_trn,ulon_trn,riu_trn,rju_trn

      namelist /datconf_trn/ intopt_trn

      namelist /project_lnd/ mpopt_lnd,nspol_lnd,                       &
     &                       tlat1_lnd,tlat2_lnd,tlon_lnd

      namelist /gridset_lnd/ xdim_lnd,ydim_lnd,dx_lnd,dy_lnd,           &
     &                       ulat_lnd,ulon_lnd,riu_lnd,rju_lnd

      namelist /datconf_lnd/ intopt_lnd,numctg_lnd,lnduse_lnd,albe_lnd, &
     &                       beta_lnd,z0m_lnd,z0h_lnd,cap_lnd,nuu_lnd

      namelist /project_sst/ mpopt_sst,nspol_sst,                       &
     &                       tlat1_sst,tlat2_sst,tlon_sst

      namelist /gridset_sst/ xdim_sst,ydim_sst,dx_sst,dy_sst,           &
     &                       ulat_sst,ulon_sst,riu_sst,rju_sst

      namelist /project_ice/ mpopt_ice,nspol_ice,                       &
     &                       tlat1_ice,tlat2_ice,tlon_ice

      namelist /gridset_ice/ xdim_ice,ydim_ice,dx_ice,dy_ice,           &
     &                       ulat_ice,ulon_ice,riu_ice,rju_ice

      namelist /uniconf_uni/ fltyp_uni,flitv_uni,                       &
     &                       rmopt_uni,uniopt_uni,ugroup_uni

      namelist /rstconf_rst/ xsub_rst,ysub_rst,flitv_rst,rmopt_rst

!-----7--------------------------------------------------------------7--

! Internal module procedure

!     none

!-----7--------------------------------------------------------------7--

      end module m_defname
