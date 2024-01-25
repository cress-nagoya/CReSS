!***********************************************************************
      module m_comaero
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/09/22
!     Modification:

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the aerosol physical constants.

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

      real, parameter :: rhodu=2.5e3
                       ! Density of soil dust

      real, parameter :: rhoca(1:2)=(/1.5e3,1.25e3/)
                       ! Density of carbonaceoues

      real, parameter :: rhosu=1.769e3
                       ! Density of sulfate

      real, parameter :: rhosa=2.25e3
                       ! Density of sea salt

      real, parameter :: diadu(1:6)=(/2.60e-7,6.60e-7,1.64e-6,2.54e-6,  &
     &                                6.40e-6,1.60e-5/)
                       ! Mean diameter of soil dust

      real, parameter :: diaca(1:4)=(/3.60e-7,7.54e-7,2.00e-7,6.24e-7/)
                       ! Mean diameter of cobonaceous

      real, parameter :: diasu(1:4)=(/2.78e-7,9.24e-7,1.39e-7,4.62e-7/)
                       ! Mean diameter of sulfate

      real, parameter :: diasa(1:4)=(/3.60e-7,1.12e-6,3.56e-6,1.12e-5/)
                       ! Mean diameter of sea salt

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

      end module m_comaero
