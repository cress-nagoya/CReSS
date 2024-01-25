!***********************************************************************
      module m_outstd12
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/01/15
!     Modification: 2001/04/15, 2001/06/29, 2002/06/18, 2002/07/15,
!                   2003/03/28, 2003/04/30, 2003/05/19, 2003/09/01,
!                   2004/04/15, 2004/05/31, 2004/06/10, 2006/01/10,
!                   2006/09/21, 2007/01/05, 2007/07/30, 2008/05/02,
!                   2008/08/25, 2008/10/10, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read in the messages to standard i/o.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
      use m_outstd04

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: outstd12, s_outstd12

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outstd12

        module procedure s_outstd12

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_outstd12(fmsg,vname,ncvn,cunit,ctime,                &
     &                      maxi,maxj,maxk,maxvl,mini,minj,mink,minvl)
!***********************************************************************

! Input variables

      character(len=4), intent(in) :: vname
                       ! Optional variable name

      character(len=7), intent(in) :: cunit
                       ! Optional dumped variable unit

      integer, intent(in) :: fmsg
                       ! Control flag of message type

      integer, intent(in) :: ncvn
                       ! Number of character of vname

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

      integer, intent(in) :: maxi
                       ! Index of maximum value point in x direction

      integer, intent(in) :: maxj
                       ! Index of maxinum value point in y direction

      integer, intent(in) :: maxk
                       ! Index of maxinum value point in z direction

      integer, intent(in) :: mini
                       ! Index of mininum value point in x direction

      integer, intent(in) :: minj
                       ! Index of mininum value point in y direction

      integer, intent(in) :: mink
                       ! Index of mininum value point in z direction

      real, intent(in) :: maxvl
                       ! Maximum value

      real, intent(in) :: minvl
                       ! Minimum value

! Internal shared variable

      integer absmsg   ! abs(fmsg)

!-----7--------------------------------------------------------------7--

! Read in the messages to standard i/o.

      absmsg=abs(fmsg)

      if(absmsg.le.1) then

        call outstd04(absmsg,ctime)

        write(6,*)

        if(fmsg.eq.1) then

          write(6,'(a)') '  maximum and minimum: procedure, outmxn;'

        else

          write(6,'(a)') '  maximum and minimum: procedure, chkmxn;'

        end if

      else if(absmsg.eq.2) then

        if(ncvn.eq.1) then

          write(6,'(a,a,a,e15.8e2,a,a,a,i6,a,i6,a,i6,a)') '    ',       &
     &          vname(1:ncvn),'max    = ',maxvl,' ',cunit(1:7),         &
     &          ' at (',maxi,',',maxj,',',maxk,')'

          write(6,'(a,a,a,e15.8e2,a,a,a,i6,a,i6,a,i6,a)') '    ',       &
     &          vname(1:ncvn),'min    = ',minvl,' ',cunit(1:7),         &
     &          ' at (',mini,',',minj,',',mink,')'

        else if(ncvn.eq.2) then

          write(6,'(a,a,a,e15.8e2,a,a,a,i6,a,i6,a,i6,a)') '    ',       &
     &          vname(1:ncvn),'max   = ',maxvl,' ',cunit(1:7),          &
     &          ' at (',maxi,',',maxj,',',maxk,')'

          write(6,'(a,a,a,e15.8e2,a,a,a,i6,a,i6,a,i6,a)') '    ',       &
     &          vname(1:ncvn),'min   = ',minvl,' ',cunit(1:7),          &
     &          ' at (',mini,',',minj,',',mink,')'

        else if(ncvn.eq.3) then

          write(6,'(a,a,a,e15.8e2,a,a,a,i6,a,i6,a,i6,a)') '    ',       &
     &          vname(1:ncvn),'max  = ',maxvl,' ',cunit(1:7),           &
     &          ' at (',maxi,',',maxj,',',maxk,')'

          write(6,'(a,a,a,e15.8e2,a,a,a,i6,a,i6,a,i6,a)') '    ',       &
     &          vname(1:ncvn),'min  = ',minvl,' ',cunit(1:7),           &
     &          ' at (',mini,',',minj,',',mink,')'

        else if(ncvn.eq.4) then

          write(6,'(a,a,a,e15.8e2,a,a,a,i6,a,i6,a,i6,a)') '    ',       &
     &          vname(1:ncvn),'max = ',maxvl,' ',cunit(1:7),            &
     &          ' at (',maxi,',',maxj,',',maxk,')'

          write(6,'(a,a,a,e15.8e2,a,a,a,i6,a,i6,a,i6,a)') '    ',       &
     &          vname(1:ncvn),'min = ',minvl,' ',cunit(1:7),            &
     &          ' at (',mini,',',minj,',',mink,')'

        end if

      else

        write(6,'(a,a)') '    Because of wrong option of mxnvar,',      &
     &                   ' no maximum and minimum value is output.'

      end if

! -----

      end subroutine s_outstd12

!-----7--------------------------------------------------------------7--

      end module m_outstd12
