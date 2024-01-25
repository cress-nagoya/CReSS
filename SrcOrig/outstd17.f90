!***********************************************************************
      module m_outstd17
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/04/11
!     Modification: 2008/05/02, 2008/08/25, 2008/10/10, 2009/02/27,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read in message to standard i/o.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: outstd17, s_outstd17

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outstd17

        module procedure s_outstd17

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
      subroutine s_outstd17(latsw,lonsw)
!***********************************************************************

! Input variables

      real, intent(in) :: latsw
                       ! Latitude at south-west corner

      real, intent(in) :: lonsw
                       ! Longitude at south-west corner

!-----7--------------------------------------------------------------7--

! Read in messages to standard i/o.

      write(6,*)

      write(6,'(a)') '  messages: procedure, outllsw;'

      if(latsw.le.-10.e0) then

        write(6,'(a,f10.5,a)')                                          &
     &    '    latitude at south-west corner, lat =',latsw,' [degree].'

      else if(latsw.gt.-10.e0.and.latsw.lt.0.e0) then

        write(6,'(a,f9.5,a)')                                           &
     &    '    latitude at south-west corner, lat =',latsw,' [degree].'

      else if(latsw.ge.0.e0.and.latsw.lt.10.e0) then

        write(6,'(a,f8.5,a)')                                           &
     &    '    latitude at south-west corner, lat =',latsw,' [degree].'

      else

        write(6,'(a,f9.5,a)')                                           &
     &    '    latitude at south-west corner, lat =',latsw,' [degree].'

      end if

      if(lonsw.le.-100.e0) then

        write(6,'(a,f11.5,a)')                                          &
     &    '    longitude at south-west corner, lon =',lonsw,' [degree].'

      else if(lonsw.gt.-100.e0.and.lonsw.le.-10.e0) then

        write(6,'(a,f10.5,a)')                                          &
     &    '    longitude at south-west corner, lon =',lonsw,' [degree].'

      else if(lonsw.gt.-10.e0.and.lonsw.lt.0.e0) then

        write(6,'(a,f9.5,a)')                                           &
     &    '    longitude at south-west corner, lon =',lonsw,' [degree].'

      else if(lonsw.ge.0.e0.and.lonsw.lt.10.e0) then

        write(6,'(a,f8.5,a)')                                           &
     &    '    longitude at south-west corner, lon =',lonsw,' [degree].'

      else if(lonsw.ge.10.e0.and.lonsw.lt.100.e0) then

        write(6,'(a,f9.5,a)')                                           &
     &    '    longitude at south-west corner, lon =',lonsw,' [degree].'

      else

        write(6,'(a,f10.5,a)')                                          &
     &    '    longitude at south-west corner, lon =',lonsw,' [degree].'

      end if

! -----

      end subroutine s_outstd17

!-----7--------------------------------------------------------------7--

      end module m_outstd17
