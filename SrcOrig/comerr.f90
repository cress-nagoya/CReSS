!***********************************************************************
      module m_comerr
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/04/06, 1999/05/10, 1999/05/20, 1999/06/07,
!                   1999/07/05, 1999/07/23, 1999/07/28, 1999/09/16,
!                   1999/09/30, 1999/10/22, 1999/11/01, 1999/11/19,
!                   1999/11/24, 1999/12/17, 2000/01/17, 2000/04/18,
!                   2000/06/01, 2000/07/05, 2000/12/18, 2001/01/15,
!                   2001/03/13, 2001/05/29, 2001/06/06, 2001/07/13,
!                   2001/08/07, 2001/10/18, 2001/11/20, 2002/06/18,
!                   2002/07/03, 2002/08/15, 2002/08/27, 2002/09/02,
!                   2002/09/09, 2002/10/15, 2002/10/31, 2002/11/11,
!                   2003/04/30, 2003/05/19, 2003/07/15, 2003/09/01,
!                   2003/10/10, 2003/11/05, 2003/12/12, 2004/01/09,
!                   2006/09/21, 2007/01/20, 2008/05/02, 2008/08/25,
!                   2008/10/10, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the error list table.

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

      integer, parameter :: nerr=200
                       ! Dimension of error list table

      character(len=14) errlst(1:nerr)
                       ! Table to check namelist errors

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

      end module m_comerr
