!***********************************************************************
      module m_vextund
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/06/27
!     Modification: 2008/05/02, 2008/08/25, 2008/12/11, 2009/01/05,
!                   2009/02/27, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     extrapolate the base state variables in undefined layer.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_commath
      use m_commpi
      use m_cpondpe
      use m_destroy
      use m_getcname
      use m_inichar

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: vextund, s_vextund

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vextund

        module procedure s_vextund

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
      subroutine s_vextund(fpgpvvar,fpetrvar_gpv,nk,z,u1d,v1d,pt1d,qv1d)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpgpvvar
                       ! Formal parameter of unique index of gpvvar

      integer, intent(in) :: fpetrvar_gpv
                       ! Formal parameter of unique index of etrvar_gpv

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: z(1:nk)
                       ! zeta coordinates

! Input and output variables

      real, intent(inout) :: u1d(1:nk)
                       ! Horizontally averaged x components of velocity

      real, intent(inout) :: v1d(1:nk)
                       ! Horizontally averaged y components of velocity

      real, intent(inout) :: pt1d(1:nk)
                       ! Horizontally averaged potential temperature

      real, intent(inout) :: qv1d(1:nk)
                       ! Horizontally averaged water vapor mixing ratio

! Internal shared variables

      character(len=108) gpvvar
                       ! Control flag of input GPV data variables

      character(len=108) etrvar_gpv
                       ! Control flag of extrapolating method

      integer k        ! Array index in z direction

      integer stat     ! Runtime status

!-----7--------------------------------------------------------------7--

! Initialize the character variables.

      call inichar(gpvvar)
      call inichar(etrvar_gpv)

! -----

! Get the required namelist variables.

      call getcname(fpgpvvar,gpvvar)
      call getcname(fpetrvar_gpv,etrvar_gpv)

! -----

! Check errors.

      stat=0

      if(etrvar_gpv(1:1).eq.'o') then

        if(u1d(nk-1).lt.lim34n.or.u1d(nk).lt.lim34n) then
          stat=stat+1
        end if

      else

        if(u1d(nk).lt.lim34n) then
          stat=stat+1
        end if

      end if

      if(etrvar_gpv(2:2).eq.'o') then

        if(v1d(nk-1).lt.lim34n.or.v1d(nk).lt.lim34n) then
          stat=stat+1
        end if

      else

        if(v1d(nk).lt.lim34n) then
          stat=stat+1
        end if

      end if

      if(etrvar_gpv(5:5).eq.'o') then

        if(pt1d(nk-1).lt.lim34n.or.pt1d(nk).lt.lim34n) then
          stat=stat+1
        end if

      else

        if(pt1d(nk).lt.lim34n) then
          stat=stat+1
        end if

      end if

      if(gpvvar(2:2).eq.'o') then

        if(etrvar_gpv(6:6).eq.'o') then

          if(qv1d(nk-1).lt.lim34n.or.qv1d(nk).lt.lim34n) then
            stat=stat+1
          end if

        else

          if(qv1d(nk).lt.lim34n) then
            stat=stat+1
          end if

        end if

      end if

! -----

! If error occured, call the procedure destroy.

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('vextund ',7,'cont',7,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('vextund ',7,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

! Perform extrapolation in the undefined layer.

      if(etrvar_gpv(1:1).eq.'o') then

        do k=nk-2,1,-1

          if(u1d(k).lt.lim34n) then

            u1d(k)=((z(k+2)-z(k))*u1d(k+1)-(z(k+1)-z(k))*u1d(k+2))      &
     &        /(z(k+2)-z(k+1))

          end if

        end do

      else

        do k=nk-1,1,-1

          if(u1d(k).lt.lim34n) then

            u1d(k)=u1d(k+1)

          end if

        end do

      end if

      if(etrvar_gpv(2:2).eq.'o') then

        do k=nk-2,1,-1

          if(v1d(k).lt.lim34n) then

            v1d(k)=((z(k+2)-z(k))*v1d(k+1)-(z(k+1)-z(k))*v1d(k+2))      &
     &        /(z(k+2)-z(k+1))

          end if

        end do

      else

        do k=nk-1,1,-1

          if(v1d(k).lt.lim34n) then

            v1d(k)=v1d(k+1)

          end if

        end do

      end if

      if(etrvar_gpv(5:5).eq.'o') then

        do k=nk-2,1,-1

          if(pt1d(k).lt.lim34n) then

            pt1d(k)=((z(k+2)-z(k))*pt1d(k+1)-(z(k+1)-z(k))*pt1d(k+2))   &
     &        /(z(k+2)-z(k+1))

          end if

        end do

      else

        do k=nk-1,1,-1

          if(pt1d(k).lt.lim34n) then

            pt1d(k)=pt1d(k+1)

          end if

        end do

      end if

      if(gpvvar(2:2).eq.'o') then

        if(etrvar_gpv(6:6).eq.'o') then

          do k=nk-2,1,-1

            if(qv1d(k).lt.lim34n) then

              qv1d(k)=((z(k+2)-z(k))*qv1d(k+1)-(z(k+1)-z(k))*qv1d(k+2)) &
     &          /(z(k+2)-z(k+1))

            end if

          end do

        else

          do k=nk-1,1,-1

            if(qv1d(k).lt.lim34n) then

              qv1d(k)=qv1d(k+1)

            end if

          end do

        end if

      end if

! -----

      end subroutine s_vextund

!-----7--------------------------------------------------------------7--

      end module m_vextund
