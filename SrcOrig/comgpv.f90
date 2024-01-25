!***********************************************************************
      module m_comgpv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/05/19
!     Modification: 2003/11/05, 2003/12/12, 2004/01/09, 2004/07/01,
!                   2006/02/03, 2007/01/20, 2008/05/02, 2008/08/19,
!                   2008/08/25, 2008/12/11, 2009/02/27, 2011/09/22

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the array for gridata.

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

      integer, allocatable, save :: land(:,:)
                       ! Land use of surface

      real, allocatable, save :: z(:)
                       ! zeta coordinates

      real, allocatable, save :: zph(:,:,:)
                       ! z physical coordinates

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
                       ! x components of velocity

      real, allocatable, save :: v(:,:,:)
                       ! y components of velocity

      real, allocatable, save :: w(:,:,:)
                       ! z components of velocity

      real, allocatable, save :: pp(:,:,:)
                       ! Pressure perturbation

      real, allocatable, save :: ptp(:,:,:)
                       ! Potential temperature perturbation

      real, allocatable, save :: qv(:,:,:)
                       ! Water vapor mixing ratio

      real, allocatable, save :: qc(:,:,:)
                       ! Cloud water mixing ratio

      real, allocatable, save :: qr(:,:,:)
                       ! Rain water mixing ratio

      real, allocatable, save :: qi(:,:,:)
                       ! Cloud ice mixing ratio

      real, allocatable, save :: qs(:,:,:)
                       ! Snow mixing ratio

      real, allocatable, save :: qg(:,:,:)
                       ! Graupel mixing ratio

      real, allocatable, save :: qh(:,:,:)
                       ! Hail mixing ratio

      real, allocatable, save :: z1d(:)
                       ! Horizontally averaged z coordinates

      real, allocatable, save :: u1d(:)
                       ! Horizontally averaged x components of velocity

      real, allocatable, save :: v1d(:)
                       ! Horizontally averaged y components of velocity

      real, allocatable, save :: p1d(:)
                       ! Horizontally averaged pressure

      real, allocatable, save :: pt1d(:)
                       ! Horizontally averaged potential temperature

      real, allocatable, save :: qv1d(:)
                       ! Horizontally averaged water vapor mixing ratio

      real, allocatable, save :: tmp1(:)
                       ! Temporary array

      real, allocatable, save :: tmp2(:,:)
                       ! Temporary array

      real, allocatable, save :: tmp3(:,:)
                       ! Temporary array

      real, allocatable, save :: tmp4(:,:)
                       ! Temporary array

      real, allocatable, save :: tmp5(:,:)
                       ! Temporary array

      real, allocatable, save :: tmp6(:,:,:)
                       ! Temporary array

      real, allocatable, save :: tmp7(:,:,:)
                       ! Temporary array

      real, allocatable, save :: londat(:,:)
                       ! Longitude in data

      real, allocatable, save :: htdat(:,:)
                       ! Terrain height in data

      real, allocatable, save :: zdat(:,:,:)
                       ! z physical coordinates in data

      real, allocatable, save :: udat(:,:,:)
                       ! x components of velocity in data

      real, allocatable, save :: vdat(:,:,:)
                       ! y components of velocity in data

      real, allocatable, save :: wdat(:,:,:)
                       ! z components of velocity in data

      real, allocatable, save :: pdat(:,:,:)
                       ! Pressure in data

      real, allocatable, save :: ptdat(:,:,:)
                       ! Potential temperature in data

      real, allocatable, save :: qvdat(:,:,:)
                       ! Water vapor mixing ratio in data

      real, allocatable, save :: qcdat(:,:,:)
                       ! Cloud water mixing ratio in data

      real, allocatable, save :: qrdat(:,:,:)
                       ! Rain water mixing ratio in data

      real, allocatable, save :: qidat(:,:,:)
                       ! Cloud ice mixing ratio in data

      real, allocatable, save :: qsdat(:,:,:)
                       ! Snow mixing ratio in data

      real, allocatable, save :: qgdat(:,:,:)
                       ! Graupel mixing ratio in data

      real, allocatable, save :: qhdat(:,:,:)
                       ! Hail mixing ratio in data

      real, allocatable, save :: pbdat(:,:,:)
                       ! Base state pressure in data

      real, allocatable, save :: ptbdat(:,:,:)
                       ! Base state potential temperature in data

      real, allocatable, save :: dtmp1(:,:)
                       ! Temporary array

      real, allocatable, save :: dtmp2(:,:,:)
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

      end module m_comgpv
