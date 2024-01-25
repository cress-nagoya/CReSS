!***********************************************************************
      module m_chkname
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/01/25, 1999/03/16, 1999/03/25,
!                   1999/04/06, 1999/05/10, 1999/05/20, 1999/06/07,
!                   1999/06/21, 1999/07/05, 1999/07/23, 1999/07/28,
!                   1999/08/03, 1999/08/09, 1999/08/18, 1999/08/23,
!                   1999/09/06, 1999/09/16, 1999/09/30, 1999/10/12,
!                   1999/10/22, 1999/10/27, 1999/11/01, 1999/11/19,
!                   1999/11/24, 1999/12/17, 2000/01/05, 2000/01/17,
!                   2000/02/02, 2000/03/17, 2000/04/18, 2000/07/05,
!                   2000/08/10, 2000/12/18, 2001/01/15, 2001/02/13,
!                   2001/03/13, 2001/04/15, 2001/05/29, 2001/06/06,
!                   2001/07/13, 2001/08/07, 2001/09/13, 2001/10/18,
!                   2001/11/20, 2002/01/07, 2002/02/05, 2002/04/02,
!                   2002/06/18, 2002/07/03, 2002/07/15, 2002/07/23,
!                   2002/08/15, 2002/08/27, 2002/09/02, 2002/09/09,
!                   2002/10/15, 2002/10/31, 2002/11/11, 2002/12/02,
!                   2002/12/17, 2003/02/05, 2003/03/13, 2003/03/28,
!                   2003/04/30, 2003/05/19, 2003/06/27, 2003/07/15,
!                   2003/08/08, 2003/09/01, 2003/10/10, 2003/10/31,
!                   2003/11/05, 2003/12/12, 2004/01/09, 2004/03/05,
!                   2004/04/01, 2004/04/15, 2004/05/07, 2004/05/31,
!                   2004/06/10, 2004/07/01, 2004/07/10, 2004/08/01,
!                   2004/08/20, 2004/09/01, 2004/09/10, 2004/09/25,
!                   2004/12/17, 2005/01/14, 2005/02/10, 2005/04/04,
!                   2005/08/05, 2005/10/05, 2005/11/22, 2005/12/13,
!                   2006/01/10, 2006/02/13, 2006/04/03, 2006/05/12,
!                   2006/06/21, 2006/07/21, 2006/09/21, 2006/09/30,
!                   2006/11/06, 2006/11/27, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/08/24, 2008/03/12, 2008/05/02,
!                   2008/08/25, 2008/10/10, 2008/12/11, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read out and set the namelist variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_chknlgpv
      use m_chknlrdr
      use m_chknlrst
      use m_chknlsfc
      use m_chknlslv
      use m_chknltrn
      use m_chknluni
      use m_chkstd
      use m_comkind
      use m_commpi
      use m_destroy
      use m_outstd04
      use m_outstd13
      use m_outstd16
      use m_rdconf

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: chkname, s_chkname

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface chkname

        module procedure s_chkname

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

!     none

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_chkname
!***********************************************************************

! Internal shared variable

      integer stat     ! Runtime status

!-----7--------------------------------------------------------------7--

!! Read out and check the namelist variables.

! Read in the message to the standard i/o.

      if(mype.eq.root) then

        call outstd04(0,0_i8)

      end if

      call chkstd(root)

! -----

! Read out the namelist variables.

      stat=0

      if(mype.eq.root) then

        call rdconf(stat)

      end if

! -----

! If error occured, call the procedure destroy.

      call chkerr(stat)

      if(stat.lt.0) then

        call destroy('rdconf  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

! Read in the messages to the standard i/o.

      if(mype.eq.root) then

        call outstd13('chkname ',7)

      end if

      call chkstd(root)

! -----

! Check the namelist variables.

      stat=0

      if(mype.eq.root) then

        call chknlslv('check   ',5,stat)

        call outstd16('solver  ',6,stat)

        call chknlgpv('check   ',5,stat)

        call outstd16('gridata ',7,stat)

        call chknlrdr('check   ',5,stat)

        call outstd16('radata  ',6,stat)

        call chknltrn('check   ',5,stat)

        call outstd16('terrain ',7,stat)

        call chknlsfc('check   ',5,stat)

        call outstd16('surface ',7,stat)

        call chknluni('check   ',5,stat)

        call outstd16('unite   ',5,stat)

        call chknlrst('check   ',5,stat)

        call outstd16('rstruct ',7,stat)

      end if

! -----

!! -----

      end subroutine s_chkname

!-----7--------------------------------------------------------------7--

      end module m_chkname
