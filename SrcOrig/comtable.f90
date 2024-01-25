!***********************************************************************
      module m_comtable
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/05/19
!     Modification: 2003/11/05, 2004/10/12, 2006/03/06, 2006/08/08,
!                   2006/09/30, 2007/01/20, 2007/07/30, 2008/03/12,
!                   2008/05/02, 2008/08/25, 2008/10/10, 2009/02/27,
!                   2009/03/12, 2009/08/20, 2009/11/05

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the referenced table.

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

      real rcdl(0:101) ! Standard referenced cloud cover in low layer
      real rcdm(0:101) ! Standard referenced cloud cover in middle layer
      real rcdh(0:101) ! Standard referenced cloud cover in high layer

      real cdrat(1:3)  ! Distribution ratio of each cloud cover

      real ckoe(0:40)  ! Koenig temperature dependent parameter
      real pkoe(0:40)  ! Koenig temperature dependent parameter

      real refch(0:11,0:17)
                       ! Standard referenced charging rate

      real nccn(1:5)   ! Cloud condensation nuclei [1/cm^3]

      real rrdbw(1:15) ! Standard referenced water bin radius [cm]

      real rrcbw(1:21) ! Ratio of standard referenced water bin radius

      real rewbw(1:15,1:21)
                       ! Standard referenced coalescence efficiency
                       ! between water bins

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

      end module m_comtable
