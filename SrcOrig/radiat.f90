!***********************************************************************
      module m_radiat
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/08/31
!     Modification: 2001/10/15, 2001/12/11, 2002/04/02, 2002/07/03,
!                   2002/07/15, 2002/12/02, 2003/04/30, 2003/05/19,
!                   2003/07/15, 2003/10/31, 2003/11/05, 2003/12/12,
!                   2004/04/15, 2004/05/07, 2004/08/01, 2004/09/01,
!                   2004/09/10, 2005/01/07, 2005/01/31, 2005/08/05,
!                   2006/08/08, 2006/09/21, 2007/01/05, 2007/01/20,
!                   2007/06/27, 2007/07/30, 2007/09/04, 2007/10/19,
!                   2008/03/12, 2008/05/02, 2008/07/01, 2008/08/25,
!                   2008/10/10, 2009/02/27, 2009/08/20, 2009/11/13,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the short and long wave radiation.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comdays
      use m_commath
      use m_comphy
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: radiat, s_radiat

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface radiat

        module procedure s_radiat

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic cos
      intrinsic sin
      intrinsic exp
      intrinsic log10
      intrinsic max
      intrinsic min
      intrinsic mod
      intrinsic real
      intrinsic sqrt

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_radiat(fpcphopt,fmois,cdate,ni,nj,nk,nund,zph,       &
     &                    lat,lon,p,t,qv,land,albe,kai,tund,tice,       &
     &                    cdl,cdm,cdh,fall,rgd,rsd,rld,rlu,             &
     &                    zph8s,zref,coseta)
!***********************************************************************

! Input variable

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      character(len=12), intent(in) :: cdate
                       ! Current forecast date
                       ! with Gregorian calendar, yyyymmddhhmm

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nund
                       ! Number of soil and sea layers

      integer, intent(in) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(in) :: lat(0:ni+1,0:nj+1)
                       ! Latitude

      real, intent(in) :: lon(0:ni+1,0:nj+1)
                       ! Longitude

      real, intent(in) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure

      real, intent(in) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

      real, intent(in) :: albe(0:ni+1,0:nj+1)
                       ! Albedo

      real, intent(in) :: kai(0:ni+1,0:nj+1)
                       ! Sea ice distribution

      real, intent(in) :: tund(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at present

      real, intent(in) :: tice(0:ni+1,0:nj+1)
                       ! Mixed ice surface temperature

      real, intent(in) :: cdl(0:ni+1,0:nj+1)
                       ! Cloud cover in lower layer

      real, intent(in) :: cdm(0:ni+1,0:nj+1)
                       ! Cloud cover in middle layer

      real, intent(in) :: cdh(0:ni+1,0:nj+1)
                       ! Cloud cover in upper layer

      real, intent(in) :: fall(0:ni+1,0:nj+1)
                       ! Precipitation flag

! Output variables

      real, intent(out) :: rgd(0:ni+1,0:nj+1)
                       ! Global solar radiation

      real, intent(out) :: rsd(0:ni+1,0:nj+1)
                       ! Net downward short wave radiation

      real, intent(out) :: rld(0:ni+1,0:nj+1)
                       ! Downward long wave radiation

      real, intent(out) :: rlu(0:ni+1,0:nj+1)
                       ! Upward long wave radiation

! Internal shared variables

      integer cphopt   ! Option for cloud micro physics

      integer nkm1     ! nk - 1

      integer cyr      ! Year of current forecast date
      integer cmo      ! Month of current forecast date
      integer cdy      ! Day of current forecast date
      integer chr      ! Hour of current forecast date
      integer cmn      ! Minite of current forecast date

      real ln1013      ! - 0.13 x ln10

      real esgm        ! el x sigma
      real esgm51      ! 0.51 x el x sigma

      real jday        ! Number of elapse of days from start of year

      real eqt         ! Equation of local time

      real phs         ! Solar angle

      real rchr        ! real(chr)
      real rcmn        ! real(cmn) / 60.0

      real sinphs      ! sin(phs)
      real cosphs      ! cos(phs)

      real, intent(inout) :: zph8s(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates at scalar points

      real, intent(inout) :: zref(0:ni+1,0:nj+1)
                       ! Reference z physical coordinates
                       ! for downward radiation

      real, intent(inout) :: coseta(0:ni+1,0:nj+1)
                       ! cos(Zenith angle)

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real tlc         ! Local time

      real pa          ! Pressure at specified plane
      real ta          ! Temperature at specified plane
      real ea          ! Pertial vapor pressure at specified plane

      real cdall       ! Total cloud cover

      real absrp       ! Absorption rate
                       ! of downward short wave radiation

      real dk          ! Distance in z direction for interpolating

      real a           ! Temporary variable
      real b           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpcphopt,cphopt)

! -----

! Read out the integer variables from the input current forecast date
! with Gregorian calendar, yyyymmddhhmm.

      read(cdate(1:12),'(i4.4,4i2.2)') cyr,cmo,cdy,chr,cmn

! -----

! Calculate the solar angle.

      if(mod(cyr,400).eq.0                                              &
     &  .or.(mod(cyr,4).eq.0.and.mod(cyr,100).ne.0)) then

        jday=2.e0*cc*i366*real(elaitc(cmo-1)+cdy-1)

      else

        jday=2.e0*cc*i365*real(ela(cmo-1)+cdy-1)

      end if

      eqt=.000075e0+.001868e0*cos(jday)-.032077e0*sin(jday)             &
     &  -.014615e0*cos(2.e0*jday)-.040849e0*sin(2.e0*jday)

      phs=.006918e0-.399912e0*cos(jday)+.070257e0*sin(jday)             &
     &  -.006758e0*cos(2.e0*jday)+.000907e0*sin(2.e0*jday)              &
     &  -.002697e0*cos(3.e0*jday)+.001480e0*sin(3.e0*jday)

! -----

! Set the common used variables.

      nkm1=nk-1

      ln1013=-.13e0*ln10

      esgm=el*sigma
      esgm51=.51e0*el*sigma

      rchr=real(chr)
      rcmn=oned60*real(cmn)

      sinphs=sin(phs)
      cosphs=cos(phs)

! -----

!!!!! Calculte the zenith angle, the global solar radiation, the net
!!!!! downward short wave radiation and the upward and downward long
!!!!! wave radiation.

!$omp parallel default(shared) private(k)

! Calculte the zenith angle.

!$omp do schedule(runtime) private(i,j,tlc)

      do j=1,nj-1
      do i=1,ni-1

        tlc=rchr+rcmn+oned15*lon(i,j)

        coseta(i,j)=sinphs*sin(lat(i,j)*d2r)                            &
     &    +cosphs*cos(lat(i,j)*d2r)*cos(eqt+15.e0*(tlc-12.e0)*d2r)

      end do
      end do

!$omp end do

! -----

! Get the z physical coordinates at scalar points and the reference
! z physical coordinates for downward radiation.

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          zph8s(i,j,k)=.5e0*(zph(i,j,k)+zph(i,j,k+1))
        end do
        end do

!$omp end do

      end do

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=1,ni-1
        zref(i,j)=min(zarad+zph(i,j,2),zph8s(i,j,nkm1))
      end do
      end do

!$omp end do

! -----

!!! In the case of dry air.

      if(fmois(1:3).eq.'dry') then

        do k=1,nk-2

!$omp do schedule(runtime) private(i,j,dk,pa,ta,absrp,a)

          do j=1,nj-1
          do i=1,ni-1

!! Perform calculation for specified plane.

            if(zref(i,j).gt.zph8s(i,j,k)                                &
     &        .and.zref(i,j).le.zph8s(i,j,k+1)) then

! Get the temperature and pressure at the specified plane.

              dk=(zref(i,j)-zph8s(i,j,k))/(zph8s(i,j,k+1)-zph8s(i,j,k))

              pa=(1.e0-dk)*p(i,j,k)+dk*p(i,j,k+1)
              ta=(1.e0-dk)*t(i,j,k)+dk*t(i,j,k+1)

! -----

! Calculate the global solar radiation and the net downward short wave
! radiation.

              if(coseta(i,j).gt.0.e0) then

                if(land(i,j).lt.0) then

                  absrp=1.e0                                            &
     &              -(9.e0*(1.e0-coseta(i,j))*albe(i,j)+albe(i,j))

                else if(land(i,j).eq.1) then

                  absrp=1.e0-(kai(i,j)*icalbe+(1.e0-kai(i,j))           &
     &              *(9.e0*(1.e0-coseta(i,j))*albe(i,j)+albe(i,j)))

                else

                  absrp=max(1.e0                                        &
     &              -(.5e0*(1.e0-coseta(i,j))*albe(i,j)+albe(i,j)),0.e0)

                end if

                rgd(i,j)=sun0*(.554e0                                   &
     &            +.43e0*exp(ln1013/(coseta(i,j)+eps)))*coseta(i,j)

                rsd(i,j)=absrp*rgd(i,j)

              else

                rgd(i,j)=0.e0
                rsd(i,j)=0.e0

              end if

! -----

! Calculate the upward and downward long wave radiation.

              a=ta*ta

              rld(i,j)=esgm51*a*a

              if(land(i,j).eq.1) then

                a=kai(i,j)*tice(i,j)+(1.e0-kai(i,j))*tund(i,j,1)

                a=a*a

              else

                a=tund(i,j,1)*tund(i,j,1)

              end if

              rlu(i,j)=esgm*a*a

! -----

            end if

!! -----

          end do
          end do

!$omp end do

        end do

!!! -----

!!!! In the case of moist air.

      else if(fmois(1:5).eq.'moist') then

!!! In the case of no cloud physics.

        if(abs(cphopt).eq.0) then

          do k=1,nk-2

!$omp do schedule(runtime) private(i,j,dk,pa,ta,ea,cdall,absrp,a,b)

            do j=1,nj-1
            do i=1,ni-1

!! Perform calculation for specified plane.

              if(zref(i,j).gt.zph8s(i,j,k)                              &
     &          .and.zref(i,j).le.zph8s(i,j,k+1)) then

! Set the common used variables.

                dk=(zref(i,j)-zph8s(i,j,k))                             &
     &            /(zph8s(i,j,k+1)-zph8s(i,j,k))

                pa=(1.e0-dk)*p(i,j,k)+dk*p(i,j,k+1)
                ta=(1.e0-dk)*t(i,j,k)+dk*t(i,j,k+1)
                ea=(1.e0-dk)*qv(i,j,k)+dk*qv(i,j,k+1)

                ea=pa*ea/(epsva+ea)

                cdall=cdl(i,j)+cdm(i,j)+cdh(i,j)

! -----

! Calculate the global solar radiation and the net downward short wave
! radiation.

                if(coseta(i,j).gt.0.e0) then

                  if(land(i,j).lt.0) then

                    absrp=1.e0-((9.e0-3.e0*cdall)                       &
     &                *(1.e0-coseta(i,j))*albe(i,j)+albe(i,j))

                  else if(land(i,j).eq.1) then

                    absrp=1.e0-(kai(i,j)*icalbe                         &
     &                +(1.e0-kai(i,j))*((9.e0-3.e0*cdall)               &
     &                *(1.e0-coseta(i,j))*albe(i,j)+albe(i,j)))

                  else

                    absrp=max(1.e0-((.5e0-oned6*cdall)                  &
     &                *(1.e0-coseta(i,j))*albe(i,j)+albe(i,j)),0.e0)

                  end if

                  b=.43e0+.00016e0*ea

                  if(ea.gt.3000.e0) then

                    a=0.e0

                  else if(ea.gt.100.e0.and.ea.le.3000.e0) then

                    a=1.12e0-b-.06e0*log10(ea)

                  else

                    a=.554e0

                  end if

                  rgd(i,j)=sun0*(a+b*exp(ln1013/(coseta(i,j)+eps)))     &
     &              *(1.e0-.7e0*cdl(i,j))*(1.e0-.6e0*cdm(i,j))          &
     &              *(1.e0-.3e0*cdh(i,j))*coseta(i,j)

                  rsd(i,j)=absrp*rgd(i,j)

                else

                  rgd(i,j)=0.e0
                  rsd(i,j)=0.e0

                end if

! -----

! Calculate the upward and downward long wave radiation.

                a=cdl(i,j)+.85e0*cdm(i,j)+.5e0*cdh(i,j)

                b=ta*ta

                rld(i,j)=esgm*b*b*(1.e0+(.66e-2*sqrt(ea)-.49e0)         &
     &            *(1.e0-(.75e0-.5e-4*ea)*a))

                if(land(i,j).eq.1) then

                  b=kai(i,j)*tice(i,j)+(1.e0-kai(i,j))*tund(i,j,1)

                  b=b*b

                else

                  b=tund(i,j,1)*tund(i,j,1)

                end if

                rlu(i,j)=esgm*b*b

! -----

              end if

!! -----

            end do
            end do

!$omp end do

          end do

!!! -----

!!! In the case of performing cloud physics.

        else

          do k=1,nk-2

!$omp do schedule(runtime) private(i,j,dk,pa,ta,ea,cdall,absrp,a,b)

            do j=1,nj-1
            do i=1,ni-1

!! Perform calculation for specified plane.

              if(zref(i,j).gt.zph8s(i,j,k)                              &
     &          .and.zref(i,j).le.zph8s(i,j,k+1)) then

! Set the common used variables.

                dk=(zref(i,j)-zph8s(i,j,k))                             &
     &            /(zph8s(i,j,k+1)-zph8s(i,j,k))

                pa=(1.e0-dk)*p(i,j,k)+dk*p(i,j,k+1)
                ta=(1.e0-dk)*t(i,j,k)+dk*t(i,j,k+1)
                ea=(1.e0-dk)*qv(i,j,k)+dk*qv(i,j,k+1)

                ea=pa*ea/(epsva+ea)

                cdall=cdl(i,j)+cdm(i,j)+cdh(i,j)

! -----

! Calculate the global solar radiation and the net downward short wave
! radiation.

                if(coseta(i,j).gt.0.e0) then

                  if(land(i,j).lt.0) then

                    absrp=1.e0-((9.e0-3.e0*cdall)                       &
     &                *(1.e0-coseta(i,j))*albe(i,j)+albe(i,j))

                  else if(land(i,j).eq.1) then

                    absrp=1.e0-(kai(i,j)*icalbe                         &
     &                +(1.e0-kai(i,j))*((9.e0-3.e0*cdall)               &
     &                *(1.e0-coseta(i,j))*albe(i,j)+albe(i,j)))

                  else

                    absrp=max(1.e0-((.5e0-oned6*cdall)                  &
     &                *(1.e0-coseta(i,j))*albe(i,j)+albe(i,j)),0.e0)

                  end if

                  b=.43e0+.00016e0*ea

                  if(ea.gt.3000.e0) then

                    a=0.e0

                  else if(ea.gt.100.e0.and.ea.le.3000.e0) then

                    a=1.12e0-b-.06e0*log10(ea)

                  else

                    a=.554e0

                  end if

                  rgd(i,j)=sun0*(a+b*exp(ln1013/(coseta(i,j)+eps)))     &
     &              *(1.e0-.7e0*cdl(i,j))*(1.e0-.6e0*cdm(i,j))          &
     &              *(1.e0-.3e0*cdh(i,j))*coseta(i,j)

                  rsd(i,j)=absrp*rgd(i,j)

                else

                  rgd(i,j)=0.e0
                  rsd(i,j)=0.e0

                end if

! -----

! Calculate the upward and downward long wave radiation.

                if(fall(i,j).gt.0.e0) then

                  a=cdl(i,j)+.85e0*cdm(i,j)+.5e0*cdh(i,j)+.1e0*cdall

                else

                  a=cdl(i,j)+.85e0*cdm(i,j)+.5e0*cdh(i,j)

                end if

                b=ta*ta

                rld(i,j)=esgm*b*b*(1.e0+(.66e-2*sqrt(ea)-.49e0)         &
     &            *(1.e0-(.75e0-.5e-4*ea)*a))

                if(land(i,j).eq.1) then

                  b=kai(i,j)*tice(i,j)+(1.e0-kai(i,j))*tund(i,j,1)

                  b=b*b

                else

                  b=tund(i,j,1)*tund(i,j,1)

                end if

                rlu(i,j)=esgm*b*b

! -----

              end if

!! -----

            end do
            end do

!$omp end do

          end do

        end if

!!! -----

      end if

!!!! -----

!$omp end parallel

!!!!! -----

      end subroutine s_radiat

!-----7--------------------------------------------------------------7--

      end module m_radiat
