!***********************************************************************
      module m_currpe
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/06/07
!     Modification: 2000/01/05, 2000/01/17, 2000/03/23, 2000/04/18,
!                   2001/02/13, 2002/06/18, 2002/08/15, 2003/02/05,
!                   2003/05/19, 2004/08/20, 2005/08/05, 2006/12/04,
!                   2007/01/05, 2007/01/20, 2007/04/11, 2008/05/02,
!                   2008/08/25, 2008/10/10, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the parameters of parallelizing from user arrangement.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comgrp
      use m_commpi

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: currpe, s_currpe

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface currpe

        module procedure s_currpe

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_currpe(pname,ncpn,fproc)
!***********************************************************************

! Input variable

      character(len=8), intent(in) :: pname
                       ! Running program name

      character(len=5), intent(in) :: fproc
                       ! Control flag of processing type

      integer, intent(in) :: ncpn
                       ! Number of character of pname

!-----7--------------------------------------------------------------7--

!! Set the parameters of parallelizing from user arrangement.

! For the all programs.

      if(pname(1:ncpn).eq.'all') then

        mysub=mod(mype,nsub)

        isub=mod(mysub,nisub)
        jsub=mysub/nisub

        mysrl=mype/nsub

        igrp=xgrp(mysrl)-1
        jgrp=ygrp(mysrl)-1

        mygrp=nigrp*jgrp+igrp

        if(grpxy(igrp,jgrp+1).lt.0) then
          ebw=1
        else
          ebw=0
        end if

        if(grpxy(igrp+2,jgrp+1).lt.0) then
          ebe=1
        else
          ebe=0
        end if

        if(grpxy(igrp+1,jgrp).lt.0) then
          ebs=1
        else
          ebs=0
        end if

        if(grpxy(igrp+1,jgrp+2).lt.0) then
          ebn=1
        else
          ebn=0
        end if

        if(grpxy(igrp,jgrp).lt.0.and.ebw.eq.0.and.ebs.eq.0) then
          ebsw=1
        else
          ebsw=0
        end if

        if(grpxy(igrp+2,jgrp).lt.0.and.ebe.eq.0.and.ebs.eq.0) then
          ebse=1
        else
          ebse=0
        end if

        if(grpxy(igrp,jgrp+2).lt.0.and.ebw.eq.0.and.ebn.eq.0) then
          ebnw=1
        else
          ebnw=0
        end if

        if(grpxy(igrp+2,jgrp+2).lt.0.and.ebe.eq.0.and.ebn.eq.0) then
          ebne=1
        else
          ebne=0
        end if

! -----

! For the program, unite.

      else if(pname(1:ncpn).eq.'unite') then

        if(fproc(1:5).eq.'mygrp') then

          igrp=xgrp(mysrl)-1
          jgrp=ygrp(mysrl)-1

          mygrp=nigrp*jgrp+igrp

        else if(fproc(1:5).eq.'ijgrp') then

          igrp=mod(mygrp,nigrp)
          jgrp=mygrp/nigrp

        else if(fproc(1:5).eq.'mysub') then

          mysub=nisub*jsub+isub

        else if(fproc(1:5).eq.'ijsub') then

          isub=mod(mysub,nisub)
          jsub=mysub/nisub

        else if(fproc(1:5).eq.'ijred') then

          ired=mod(myred,nired)
          jred=myred/nired

          igrp=ired+iwred
          jgrp=jred+jsred

          mygrp=nigrp*jgrp+igrp

        end if

! -----

! For the program, rstruct.

      else if(pname(1:ncpn).eq.'rstruct') then

        if(fproc(1:5).eq.'ijgrp') then

          igrp=xgrp(mysrl)-1
          jgrp=ygrp(mysrl)-1

          mygrp=nigrp*jgrp+igrp

          if(grpxy(igrp,jgrp+1).lt.0) then
            ebw=1
          else
            ebw=0
          end if

          if(grpxy(igrp+2,jgrp+1).lt.0) then
            ebe=1
          else
            ebe=0
          end if

          if(grpxy(igrp+1,jgrp).lt.0) then
            ebs=1
          else
            ebs=0
          end if

          if(grpxy(igrp+1,jgrp+2).lt.0) then
            ebn=1
          else
            ebn=0
          end if

        else if(fproc(1:5).eq.'ijsub') then

          isub=mod(mysub,nisub)
          jsub=mysub/nisub

          mype=mysrl*nsrl+(nisub*jsub+isub)

        else if(fproc(1:5).eq.'ijrst') then

          isub_rst=mod(mysub_rst,nisub_rst)
          jsub_rst=mysub_rst/nisub_rst

          mype=mysrl*nsrl+(nisub_rst*jsub_rst+isub_rst)

        end if

      end if

! -----

!! -----

      end subroutine s_currpe

!-----7--------------------------------------------------------------7--

      end module m_currpe
