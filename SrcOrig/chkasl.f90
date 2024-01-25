!***********************************************************************
      module m_chkasl
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/09/22
!     Modification: 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     check the aerosol data variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkmxn
      use m_comkind
      use m_outstd12

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: chkasl, s_chkasl

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface chkasl

        module procedure s_chkasl

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
      subroutine s_chkasl(ctime,stat,nid,njd,nkd,nqa,zdat,qadat)
!***********************************************************************

! Input variables

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

      integer, intent(in) :: nid
                       ! Data dimension in x direction

      integer, intent(in) :: njd
                       ! Data dimension in y direction

      integer, intent(in) :: nkd
                       ! Data dimension in z direction

      integer, intent(in) :: nqa(0:4)
                       ! Number of types of aerosol

      real, intent(in) :: zdat(1:nid,1:njd,1:nkd)
                       ! z physical coordinates in data

      real, intent(in) :: qadat(1:nid,1:njd,1:nkd,1:nqa(0))
                       ! Aerosol mixing ratio in data

! Output variable

      integer, intent(out) :: stat
                       ! Runtime status

! Internal shared variables

      integer istr     ! Start index of types of aerosol array
      integer iend     ! End index of types of aerosol array

      integer n        ! Array index in 4th direction

!-----7--------------------------------------------------------------7--

! Initialize the runtime status.

      stat=0

! -----

! Read in the message to standard i/o.

      call outstd12(-1,'    ',4,'       ',ctime,0,0,1,0.e0,0,0,1,0.e0)

! -----

! Check the aerosol data variables.

      call chkmxn('z   ',1,'[m]    ','all  ',ctime,stat,                &
     &            -1000.e0,1000000.e0,0.e0,0.e0,nid,njd,nkd,zdat)

      if(stat.ne.0) then

        return

      end if

      istr=1
      iend=nqa(1)

      do n=istr,iend

        call s_chkmxn('dust',4,'[kg/kg]','all  ',ctime,stat,            &
     &                0.e0,.1e0,0.e0,0.e0,nid,njd,nkd,qadat(1,1,1,n))

        if(stat.ne.0) then

          return

        end if

      end do

      istr=iend+1
      iend=iend+nqa(2)

      do n=istr,iend

        call s_chkmxn('cbon',4,'[kg/kg]','all  ',ctime,stat,            &
     &                0.e0,.1e0,0.e0,0.e0,nid,njd,nkd,qadat(1,1,1,n))

        if(stat.ne.0) then

          return

        end if

      end do

      istr=iend+1
      iend=iend+nqa(3)

      do n=istr,iend

        call s_chkmxn('sulf',4,'[kg/kg]','all  ',ctime,stat,            &
     &                0.e0,.1e0,0.e0,0.e0,nid,njd,nkd,qadat(1,1,1,n))

        if(stat.ne.0) then

          return

        end if

      end do

      istr=iend+1
      iend=iend+nqa(4)

      do n=istr,iend

        call s_chkmxn('salt',4,'[kg/kg]','all  ',ctime,stat,            &
     &                0.e0,.1e0,0.e0,0.e0,nid,njd,nkd,qadat(1,1,1,n))

        if(stat.ne.0) then

          return

        end if

      end do

! -----

      end subroutine s_chkasl

!-----7--------------------------------------------------------------7--

      end module m_chkasl
