!***********************************************************************
      module m_defdim
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/05/19
!     Modification: 2003/07/15, 2003/12/12, 2004/07/01, 2004/08/01,
!                   2005/08/05, 2006/01/10, 2006/09/30, 2007/01/20,
!                   2008/05/02, 2008/07/01, 2008/08/25, 2008/12/11,
!                   2009/02/27, 2011/08/18, 2011/09/22

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the model dimension variables.

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

      integer ni       ! Model dimension in x direction
      integer nj       ! Model dimension in y direction
      integer nk       ! Model dimension in z direction

      integer nqw      ! Number of categories of water hydrometeor
      integer nnw      ! Number of categories of water concentrations
      integer nqi      ! Number of categories of ice hydrometeor
      integer nni      ! Number of categories of ice concentrations

      integer km       ! Dimension of max(nk, nqw, nqi)

      integer nqa(0:4) ! Number of types of aerosol

      integer nund     ! Number of soil and sea layers

      integer nlev     ! Horizontally averaged vertical dimension

      integer nid_gpv  ! GPV data dimension in x direction
      integer njd_gpv  ! GPV data dimension in y direction
      integer nkd_gpv  ! GPV data dimension in z direction

      integer km_gpv   ! Dimension of max(nk, nkd_gpv)

      integer nid_asl  ! Aerosol data dimension in x direction
      integer njd_asl  ! Aerosol data dimension in y direction
      integer nkd_asl  ! Aerosol data dimension in z direction

      integer km_asl   ! Dimension of max(nk, nkd_asl)

      integer nid_rdr  ! Radar data dimension in x direction
      integer njd_rdr  ! Radar data dimension in y direction
      integer nkd_rdr  ! Radar data dimension in z direction

      integer km_rdr   ! Dimension of max(nk, nkd_rdr)

      integer nid_trn  ! Terrain data dimension in x direction
      integer njd_trn  ! Terrain data dimension in y direction

      integer nid_lnd  ! Land use data dimension in x direction
      integer njd_lnd  ! Land use data dimension in y direction

      integer nid_sst  ! Sea surface temperature data dimension
                       ! in x direction

      integer njd_sst  ! Sea surface temperature data dimension
                       ! in y direction

      integer nid_ice  ! Sea ice distribution data dimension
                       ! in x direction

      integer njd_ice  ! Sea ice distribution data dimension
                       ! in y direction

      integer ni_uni   ! Model dimension of unite in x direction
      integer nj_uni   ! Model dimension of unite in y direction

      integer nio_uni  ! Maximum unit number of unite

      integer ni_rst   ! Restructed files dimension in x direction
      integer nj_rst   ! Restructed files dimension in y direction

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

      end module m_defdim
