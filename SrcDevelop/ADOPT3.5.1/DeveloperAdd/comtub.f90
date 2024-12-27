!***********************************************************************
      module m_comtub
!***********************************************************************

!     Author      : Satoki Tsujino
!     Date        : 2014/01/28
!     Modification: 2017/06/02   Adding turbulence terms for pt and qv
!                   2017/06/18   Adding numerical diffusion terms
!                   2018/01/05   Adding eddy vis and dif coefficients

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the array for tubulence.

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

      real, allocatable, save :: turbu(:,:,:)
                       ! x component of frictional force of tubulence process

      real, allocatable, save :: turbv(:,:,:)
                       ! y component of frictional force of tubulence process

      real, allocatable, save :: turbpt(:,:,:)
                       ! heating of tubulence process

      real, allocatable, save :: turbqv(:,:,:)
                       ! moistening of tubulence process

      real, allocatable, save :: numdu(:,:,:)
                       ! x component of numerical diffusion process

      real, allocatable, save :: numdv(:,:,:)
                       ! y component of numerical diffusion process

      real, allocatable, save :: numdw(:,:,:)
                       ! z component of numerical diffusion process

      real, allocatable, save :: numdpt(:,:,:)
                       ! heating of numerical diffusion process

      real, allocatable, save :: numdqv(:,:,:)
                       ! moistening of numerical diffusion process

      real, allocatable, save :: difh(:,:,:)
                       ! eddy diffusion coefficient in horizontal

      real, allocatable, save :: difv(:,:,:)
                       ! eddy diffusion coefficient in vertical

      real, allocatable, save :: vish(:,:,:)
                       ! eddy viscousity coefficient in horizontal

      real, allocatable, save :: visv(:,:,:)
                       ! eddy viscousity coefficient in vertical

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

      end module m_comtub
