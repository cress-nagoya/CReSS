!***********************************************************************
      module m_setproj
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/03/17, 1999/03/25, 1999/06/07, 1999/09/30,
!                   1999/10/12, 1999/11/01, 2000/01/17, 2000/07/05,
!                   2001/01/09, 2001/02/13, 2001/11/20, 2002/04/02,
!                   2003/05/19, 2003/09/01, 2004/05/31, 2004/09/01,
!                   2005/02/10, 2006/11/06, 2008/05/02, 2008/08/25,
!                   2008/10/10, 2009/02/27, 2011/08/09, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the map projection parameters.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comphy
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: setproj, s_setproj, s_setproj_2

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface setproj

        module procedure s_setproj, s_setproj_2

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic cos
      intrinsic sin
      intrinsic tan
      intrinsic exp
      intrinsic log
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_setproj(fpmpopt,fpnspol,fpdx,fpdy,fpulat,fpulon,     &
     &                     fpriu,fprju,fptlat1,fptlat2,fptlon,          &
     &                     pname,ncpn,x0,y0,cpj)
!***********************************************************************

! Input variables

      character(len=8), intent(in) :: pname
                       ! Running program name

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpnspol
                       ! Formal parameter of unique index of nspol

      integer, intent(in) :: fpdx
                       ! Formal parameter of unique index of dx

      integer, intent(in) :: fpdy
                       ! Formal parameter of unique index of dy

      integer, intent(in) :: fpulat
                       ! Formal parameter of unique index of ulat

      integer, intent(in) :: fpulon
                       ! Formal parameter of unique index of ulon

      integer, intent(in) :: fpriu
                       ! Formal parameter of unique index of riu

      integer, intent(in) :: fprju
                       ! Formal parameter of unique index of rju

      integer, intent(in) :: fptlat1
                       ! Formal parameter of unique index of tlat1

      integer, intent(in) :: fptlat2
                       ! Formal parameter of unique index of tlat2

      integer, intent(in) :: fptlon
                       ! Formal parameter of unique index of tlon

      integer, intent(in) :: ncpn
                       ! Number of character of pname

! Output variables

      real, intent(out) :: x0
                       ! x origin

      real, intent(out) :: y0
                       ! y origin

      real, intent(out) :: cpj(1:7)
                       ! Map projection parameters

! Internal shared variables

      integer mpopt    ! Option for map projection
      integer nspol    ! Option for projected region

      real dx          ! Grid distance in x direction
      real dy          ! Grid distance in y direction

      real ulat        ! User specified latitude
      real ulon        ! User specified longitude

      real riu         ! User specified real index in x direction
      real rju         ! User specified real index in y direction

      real tlat1       ! True latitude 1
      real tlat2       ! True latitude 2

      real tlon        ! True longitude

      real ru          ! User specified radius
                       ! on map coordinates system

      real xu          ! User specified x coordinate
                       ! on map coordinates system

      real yu          ! User specified y coordinate
                       ! on map coordinates system

      real dlon        ! Distance of longitude

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpmpopt,mpopt)
      call getrname(fpdx,dx)
      call getrname(fpdy,dy)
      call getrname(fpulat,ulat)
      call getrname(fpulon,ulon)
      call getrname(fpriu,riu)
      call getrname(fprju,rju)

! -----

!! Calculate the parameters for latitude and longitude coordinates.

      if(mpopt.eq.0.or.mpopt.eq.10) then

! For the program, solver.

        if(pname(1:ncpn).eq.'solver') then

          cpj(1)=1.e0
          cpj(2)=rearth
          cpj(3)=1.e0/cpj(2)
          cpj(4)=0.e0
          cpj(5)=0.e0
          cpj(6)=0.e0
          cpj(7)=0.e0

          xu=cpj(2)*ulon*d2r
          yu=cpj(2)*ulat*d2r

          x0=xu-(riu-2.e0)*dx
          y0=yu-(rju-2.e0)*dy

! -----

! For the pre and post processors.

        else

          cpj(1)=0.e0
          cpj(2)=0.e0
          cpj(3)=0.e0
          cpj(4)=0.e0
          cpj(5)=0.e0
          cpj(6)=0.e0
          cpj(7)=0.e0

          x0=ulon-(riu-2.e0)*dx
          y0=ulat-(rju-2.e0)*dy

        end if

! -----

!! -----

!! Calculate the map projection parameters for the Polar Stereographic
!! projection method.

      else if(mpopt.eq.1) then

! Get the required namelist variables.

        call getiname(fpnspol,nspol)
        call getrname(fptlat1,tlat1)
        call getrname(fptlon,tlon)

! -----

! Get the map projection parameters.

        cpj(1)=1.e0+sin(real(nspol)*tlat1*d2r)
        cpj(2)=rearth*cpj(1)
        cpj(3)=1.e0/cpj(2)
        cpj(4)=0.e0
        cpj(5)=0.e0
        cpj(6)=0.e0
        cpj(7)=0.e0

        dlon=ulon-tlon

        ru=cpj(2)*cos(ulat*d2r)/((1.e0+sin(real(nspol)*ulat*d2r))+eps)

        xu=ru*sin(dlon*d2r)
        yu=-ru*cos(dlon*d2r)

        x0=xu-(riu-2.e0)*dx
        y0=yu-real(nspol)*(rju-2.e0)*dy

! -----

!! -----

!! Calculate the map projection parameters for the Lambert Conformal
!! Conic projection method.

      else if(mpopt.eq.2) then

! Get the required namelist variables.

        call getiname(fpnspol,nspol)
        call getrname(fptlat1,tlat1)
        call getrname(fptlat2,tlat2)
        call getrname(fptlon,tlon)

! -----

! Get the map projection parameters.

        cpj(1)=cos(tlat1*d2r)

        cpj(2)=tan((45.e0-.5e0*real(nspol)*tlat1)*d2r)

        cpj(3)=1.e0/cpj(2)

        cpj(4)=log(cpj(1)/cos(tlat2*d2r))                               &
     &    /log(cpj(2)/tan((45.e0-.5e0*real(nspol)*tlat2)*d2r))

        cpj(5)=1.e0/cpj(4)

        cpj(6)=rearth*cpj(5)*cpj(1)

        cpj(7)=1.e0/cpj(6)

        dlon=ulon-tlon

        ru=cpj(6)*exp(cpj(4)                                            &
     &    *log(tan((45.e0-.5e0*real(nspol)*ulat)*d2r)*cpj(3)+eps))

        xu=ru*sin(cpj(4)*dlon*d2r)
        yu=-ru*cos(cpj(4)*dlon*d2r)

        x0=xu-(riu-2.e0)*dx
        y0=yu-real(nspol)*(rju-2.e0)*dy

! -----

!! -----

!! Calculate the map projection parameters for the Mercator projection
!! method.

      else if(mpopt.eq.3.or.mpopt.eq.13) then

! Get the required namelist variables.

        call getiname(fpnspol,nspol)
        call getrname(fptlat1,tlat1)

! -----

! Get the map projection parameters.

        cpj(1)=cos(tlat1*d2r)
        cpj(2)=rearth*cpj(1)
        cpj(3)=1.e0/cpj(2)
        cpj(4)=0.e0
        cpj(5)=0.e0
        cpj(6)=0.e0
        cpj(7)=0.e0

        xu=cpj(2)*ulon*d2r

        yu=-real(nspol)                                                 &
     &    *cpj(2)*log(tan((45.e0-.5e0*real(nspol)*ulat)*d2r)+eps)

        x0=xu-(riu-2.e0)*dx
        y0=yu-(rju-2.e0)*dy

! -----

!! -----

!! Calculate the parameters without any projection method.

      else if(mpopt.eq.4) then

! Get the required namelist variable.

        call getrname(fptlon,tlon)

! -----

! Get the map projection parameters.

        cpj(1)=d2r*rearth
        cpj(2)=r2d/rearth
        cpj(3)=0.e0
        cpj(4)=0.e0
        cpj(5)=0.e0
        cpj(6)=0.e0
        cpj(7)=0.e0

        dlon=ulon-tlon

        if(dlon.lt.-180.e0) then
          dlon=dlon+360.e0
        end if

        if(dlon.gt.180.e0) then
          dlon=dlon-360.e0
        end if

        xu=cpj(1)*dlon*cos(ulat*d2r)
        yu=cpj(1)*ulat

        x0=xu-(riu-2.e0)*dx
        y0=yu-(rju-2.e0)*dy

! -----

!! -----

! Calculate the parameters with the circular cylinder coordinates
! system.

      else if(mpopt.eq.5) then

        cpj(1)=d2r*rearth
        cpj(2)=r2d/rearth
        cpj(3)=0.e0
        cpj(4)=0.e0
        cpj(5)=0.e0
        cpj(6)=0.e0
        cpj(7)=0.e0

        xu=cpj(1)*ulon*cos(ulat*d2r)
        yu=cpj(1)*ulat

        x0=xu-(riu-2.e0)*dx
        y0=yu-(rju-2.e0)*dy

      end if

! -----

      end subroutine s_setproj

!***********************************************************************
      subroutine s_setproj_2(fpmpopt,fpnspol,fptlat1,fptlat2,cpj)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpnspol
                       ! Formal parameter of unique index of nspol

      integer, intent(in) :: fptlat1
                       ! Formal parameter of unique index of tlat1

      integer, intent(in) :: fptlat2
                       ! Formal parameter of unique index of tlat2

! Output variable

      real, intent(out) :: cpj(1:7)
                       ! Map projection parameters

! Internal shared variables

      integer mpopt    ! Option for map projection
      integer nspol    ! Option for projected region

      real tlat1       ! True latitude 1
      real tlat2       ! True latitude 2

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpmpopt,mpopt)
      call getiname(fpnspol,nspol)
      call getrname(fptlat1,tlat1)
      call getrname(fptlat2,tlat2)

! -----

! Calculate the map projection parameters.

      if(mpopt.eq.0.or.mpopt.eq.10) then

        cpj(1)=1.e0
        cpj(2)=rearth
        cpj(3)=1.e0/cpj(2)
        cpj(4)=0.e0
        cpj(5)=0.e0
        cpj(6)=0.e0
        cpj(7)=0.e0

      else if(mpopt.eq.1) then

        cpj(1)=1.e0+sin(real(nspol)*tlat1*d2r)
        cpj(2)=rearth*cpj(1)
        cpj(3)=1.e0/cpj(2)
        cpj(4)=0.e0
        cpj(5)=0.e0
        cpj(6)=0.e0
        cpj(7)=0.e0

      else if(mpopt.eq.2) then

        cpj(1)=cos(tlat1*d2r)

        cpj(2)=tan((45.e0-.5e0*real(nspol)*tlat1)*d2r)

        cpj(3)=1.e0/cpj(2)

        cpj(4)=log(cpj(1)/cos(tlat2*d2r))                               &
     &    /log(cpj(2)/tan((45.e0-.5e0*real(nspol)*tlat2)*d2r))

        cpj(5)=1.e0/cpj(4)

        cpj(6)=rearth*cpj(5)*cpj(1)

        cpj(7)=1.e0/cpj(6)

      else if(mpopt.eq.3.or.mpopt.eq.13) then

        cpj(1)=cos(tlat1*d2r)
        cpj(2)=rearth*cpj(1)
        cpj(3)=1.e0/cpj(2)
        cpj(4)=0.e0
        cpj(5)=0.e0
        cpj(6)=0.e0
        cpj(7)=0.e0

      else if(mpopt.eq.4) then

        cpj(1)=d2r*rearth
        cpj(2)=r2d/rearth
        cpj(3)=0.e0
        cpj(4)=0.e0
        cpj(5)=0.e0
        cpj(6)=0.e0
        cpj(7)=0.e0

      else if(mpopt.eq.5) then

        cpj(1)=d2r*rearth
        cpj(2)=r2d/rearth
        cpj(3)=0.e0
        cpj(4)=0.e0
        cpj(5)=0.e0
        cpj(6)=0.e0
        cpj(7)=0.e0

      end if

! -----

      end subroutine s_setproj_2

!-----7--------------------------------------------------------------7--

      end module m_setproj
