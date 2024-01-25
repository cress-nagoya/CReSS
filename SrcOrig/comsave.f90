!***********************************************************************
      module m_comsave
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2005/02/10
!     Modification: 2006/12/04, 2008/05/02, 2008/08/25, 2008/12/11,
!                   2009/02/27, 2011/09/22, 2011/11/10

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the external referenced variables to avoid save attribute
!     statement.

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

      integer nxtgpv   ! Next time to read out GPV data file

      integer nxtasl   ! Next time to read out aerosol data file

      integer nxtrdr(1:2)
                       ! Next time to read out radar data file

      integer nxtsst   ! Next time to read out
                       ! sea surface temperature data file

      integer extcom   ! Control flag of common file name extension

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

      end module m_comsave
