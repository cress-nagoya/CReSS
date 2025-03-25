!***********************************************************************
      module m_comname
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
!                   2006/09/21, 2007/01/20, 2007/04/24, 2007/05/21,
!                   2007/07/30, 2007/09/04, 2008/01/11, 2008/04/17,
!                   2008/05/02, 2008/07/01, 2008/08/25, 2008/10/10,
!                   2009/01/05, 2009/01/30, 2009/02/27, 2011/05/16,
!                   2011/08/09, 2011/08/18, 2011/09/22, 2011/11/10,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the table to archive the namelist variables.

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

      integer, parameter :: ncn=21
                       ! Dimension of character namelist table

      integer, parameter :: nin=219
                       ! Dimension of integer namelist table

      integer, parameter :: nrn=798
                       ! Dimension of real namelist table

      character(len=108) cname(1:ncn)
                       ! Character namelist table

      character(len=108) rcname(1:ncn)
                       ! Read out character namelist table

      integer iname(-1:nin)
                       ! Integer namelist table

      integer riname(1:nin)
                       ! Read out integer namelist table

      real rname(1:nrn)
                       ! Real namelist table

      real rrname(1:nrn)
                       ! Read out real namelist table

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

      end module m_comname
