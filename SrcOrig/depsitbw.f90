!***********************************************************************
      module m_depsitbw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/08/08
!     Modification: 2006/09/30, 2007/04/11, 2007/10/19, 2008/01/11,
!                   2008/05/02, 2008/07/01, 2008/08/25, 2009/02/27,
!                   2009/11/13, 2011/03/18, 2011/08/18, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform deposition for water bin.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comphy
      use m_remapbw

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: depsitbw, s_depsitbw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface depsitbw

        module procedure s_depsitbw

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic log
      intrinsic max
      intrinsic sqrt

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_depsitbw(k,dtcond,ni,nj,nk,nqw,nnw,rbv,pi,p,t,       &
     &                      brw,rbrw,bmw,rbmw,dbmw,ptp,qv,mwbin,nwbin,  &
     &                      bmws,mws,nws,ssw,lv,kp,dv,dm,               &
     &                      cn,cd1,cd2,c1,c2,c3)
!***********************************************************************

! Input variables

      integer, intent(in) :: k
                       ! Array index in z direction

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

      integer, intent(in) :: nnw
                       ! Number of categories of water concentrations

      real, intent(in) :: dtcond
                       ! Time interval of condensation processes

      real, intent(in) :: rbv(0:ni+1,0:nj+1)
                       ! Inverse of base state density [cm^3/g]

      real, intent(in) :: pi(0:ni+1,0:nj+1)
                       ! Exnar function

      real, intent(in) :: p(0:ni+1,0:nj+1)
                       ! Pressure [dyn/cm^2]

      real, intent(in) :: t(0:ni+1,0:nj+1)
                       ! Air temperature

      real, intent(in) :: brw(1:nqw+1)
                       ! Radius of water bin boundary [cm]

      real, intent(in) :: rbrw(1:nqw+1,1:8)
                       ! Related parameters of brw

      real, intent(in) :: bmw(1:nqw+1,1:3)
                       ! Mass at water bin boundary [g]

      real, intent(in) :: rbmw(1:nqw,1:2)
                       ! Related parameters of bmw

      real, intent(in) :: dbmw(1:nqw)
                       ! Differential between adjacent water bins [g]

! Input and output variables

      real, intent(inout) :: ptp(0:ni+1,0:nj+1)
                       ! Potential temperature perturbation

      real, intent(inout) :: qv(0:ni+1,0:nj+1)
                       ! Water vapor mixing ratio [g/g]

      real, intent(inout) :: mwbin(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Total water mass [g/cm^3]

      real, intent(inout) :: nwbin(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations [1/cm^3]

! Internal shared variables

      integer nqw_sub  ! Substitute for nqw
      integer nnw_sub  ! Substitute for nnw

      real cc43rw      ! 3.0 x 10^3 / (4.0 x cc x rhow)

      real cc2rv4      ! sqrt(0.5 x 10^-4 x cc / rv)
      real cc4rv4      ! 4.0 x 10^4 x cc x rv

      real ccrdcp      ! 4.0 x 10^-5
                       ! x sqrt(0.125 x cc / rd) / (7.0 x cp)

      real cp4         ! 1.0 x 10^4 x cp

      real rv28        ! 1.0 x 10^8 x rv x rv

      real p205        ! 1.0 x 10^5 x p20

      real lv04        ! 1.0 x 10^4 x lv0

      real kp05        ! 1.0 x 10^5 x kp0

      real es01        ! 10.0 x es0

      real t23iv       ! 1.0 / t23

      real rwrvws      ! 0.2 x wsten / (rhow x rv)

      real dc0iv2      ! 2.0 / dc0

      real cdv         ! dv0^0.55249 / t0

      real, intent(inout) :: bmws(0:ni+1,0:nj+1,1:nqw+1)
                       ! Shifted mass at water bin boundary

      real, intent(inout) :: mws(0:ni+1,0:nj+1,1:nqw)
                       ! Shifted water mass

      real, intent(inout) :: nws(0:ni+1,0:nj+1,1:nnw)
                       ! Shifted water concentrations

      real, intent(inout) :: ssw(0:ni+1,0:nj+1)
                       ! Super saturation indices

      real, intent(inout) :: lv(0:ni+1,0:nj+1)
                       ! Latent heat of evapolation

      real, intent(inout) :: kp(0:ni+1,0:nj+1)
                       ! Thermal conductivity of water

      real, intent(inout) :: dv(0:ni+1,0:nj+1)
                       ! Molecular diffusivity of water

      real, intent(inout) :: dm(0:ni+1,0:nj+1)
                       ! Total alternations of water mass

      real, intent(inout) :: cn(0:ni+1,0:nj+1)
                       ! Common used temporary array

      real, intent(inout) :: cd1(0:ni+1,0:nj+1)
                       ! Common used temporary array

      real, intent(inout) :: cd2(0:ni+1,0:nj+1)
                       ! Common used temporary array

      real, intent(inout) :: c1(0:ni+1,0:nj+1)
                       ! Common used temporary array

      real, intent(inout) :: c2(0:ni+1,0:nj+1)
                       ! Common used temporary array

      real, intent(inout) :: c3(0:ni+1,0:nj+1)
                       ! Common used temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      integer n        ! Array index in water bin categories

      real esw         ! Saturation vapor pressure for water
      real qvsw        ! Saturation mixing ratio for water

      real mwbr        ! Mean water mass
      real rwbr        ! Mean water radius

      real rwbriv      ! Inverse of rwbr

      real sswr        ! Temporary variable
      real dvm         ! Temporary variable
      real kpm         ! Temporary variable
      real ttmp        ! Temporary variable

      real a           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Set the substituted variables.

      nqw_sub=nqw
      nnw_sub=nnw

! -----

! Set the common used variables.

      cc43rw=3.e3/(4.e0*cc*rhow)

      cc2rv4=sqrt(.5e-4*cc/rv)
      cc4rv4=4.e4*cc*rv

      ccrdcp=4.e-5*sqrt(.125e0*cc/rd)/(7.e0*cp)

      cp4=1.e4*cp

      rv28=1.e8*rv*rv

      p205=1.e5*p20

      lv04=1.e4*lv0

      kp05=1.e5*kp0

      es01=10.e0*es0

      t23iv=1.e0/t23

      rwrvws=.2e0*wsten/(rhow*rv)

      dc0iv2=2.e0/dc0

      cdv=exp(.55249e0*log(dv0))/t0

! -----

!! Perform deposition for water bin.

!$omp parallel default(shared) private(n)

! Set the common used variables.

!$omp do schedule(runtime) private(i,j,esw,qvsw,ttmp)

      do j=1,nj-1
      do i=1,ni-1

        esw=es01*exp(17.269e0*(t(i,j)-t0)/(t(i,j)-35.86e0))

        qvsw=epsva*esw/(p(i,j)-esw)

        ssw(i,j)=qv(i,j)/qvsw

        lv(i,j)=lv04*exp((.167e0+3.67e-4*t(i,j))*log(t0/t(i,j)))

        kp(i,j)=(kp05/(120.e0+t(i,j)))*exp(1.5e0*log(t23iv*t(i,j)))

        dv(i,j)=p205*exp(1.81e0*log(cdv*t(i,j)))/p(i,j)

        ttmp=t(i,j)*t(i,j)

        cn(i,j)=cc4rv4*esw*ttmp

        cd1(i,j)=rv28*ttmp*t(i,j)
        cd2(i,j)=esw*lv(i,j)*lv(i,j)

        ttmp=1.e0/t(i,j)

        c1(i,j)=rwrvws*ttmp

        ttmp=sqrt(ttmp)

        c2(i,j)=cc2rv4*ttmp*dv(i,j)
        c3(i,j)=ccrdcp*ttmp*kp(i,j)*rbv(i,j)

        dm(i,j)=0.e0

      end do
      end do

!$omp end do

! -----

! Shift the bin boundaries.

      do n=1,nqw+1

!$omp do schedule(runtime) private(i,j,sswr,dvm,kpm)

        do j=1,nj-1
        do i=1,ni-1

          sswr=exp(c1(i,j)*rbrw(n,6))

          dvm=dv(i,j)/(c2(i,j)*rbrw(n,7)+brw(n)/(brw(n)+c2(i,j)))
          kpm=kp(i,j)/(rbrw(n,8)+c3(i,j)*rbrw(n,6))

          bmws(i,j,n)=bmw(n,1)+dvm*kpm*brw(n)*cn(i,j)*(ssw(i,j)-sswr)   &
     &      /(kpm*cd1(i,j)+dvm*cd2(i,j))*dtcond

        end do
        end do

!$omp end do

      end do

! -----

! Shift the mean mass.

      do n=1,nqw

!$omp do schedule(runtime) private(i,j,mwbr,rwbr,rwbriv,sswr,dvm,kpm)

        do j=1,nj-1
        do i=1,ni-1

          if(mwbin(i,j,k,n).gt.0.e0) then

            mwbr=mwbin(i,j,k,n)/nwbin(i,j,k,n)
            rwbr=exp(oned3*log(cc43rw*mwbr))

            rwbriv=1.e0/rwbr

            sswr=exp(c1(i,j)*rwbriv)

            dvm=dv(i,j)/(dc0iv2*rwbriv*c2(i,j)+rwbr/(rwbr+c2(i,j)))
            kpm=kp(i,j)/(rwbr/(rwbr+2.16e-5)+c3(i,j)*rwbriv)

            mws(i,j,n)=(mwbr+dvm*kpm*rwbr*cn(i,j)*(ssw(i,j)-sswr)       &
     &        /(kpm*cd1(i,j)+dvm*cd2(i,j))*dtcond)*nwbin(i,j,k,n)

            if(mws(i,j,n).gt.0.e0) then

              if(bmws(i,j,n+1).gt.bmw(1,1)) then

                nws(i,j,n)=nwbin(i,j,k,n)

                dm(i,j)=dm(i,j)+(mws(i,j,n)-mwbin(i,j,k,n))

              else

                mws(i,j,n)=0.e0

                dm(i,j)=dm(i,j)-mwbin(i,j,k,n)

              end if

            else

              mws(i,j,n)=0.e0

              dm(i,j)=dm(i,j)-mwbin(i,j,k,n)

            end if

            mwbin(i,j,k,n)=0.e0
            nwbin(i,j,k,n)=0.e0

          else

            mws(i,j,n)=0.e0

          end if

        end do
        end do

!$omp end do

      end do

! -----

! Solve the new potential temperature and water vapor mixing ratio.

!$omp do schedule(runtime) private(i,j,a)

      do j=1,nj-1
      do i=1,ni-1

        a=rbv(i,j)*dm(i,j)

        ptp(i,j)=ptp(i,j)+a*lv(i,j)/(cp4*pi(i,j))

        qv(i,j)=max(qv(i,j)-a,0.e0)

      end do
      end do

!$omp end do

! -----

!$omp end parallel

!! -----

! Perform remapping.

      call remapbw(k,2,ni,nj,nk,nqw,nnw,nqw_sub,nnw_sub,bmw,rbmw,dbmw,  &
     &             bmws,mws,nws,mwbin,nwbin,cn,c1,c2,c3)

! -----

      end subroutine s_depsitbw

!-----7--------------------------------------------------------------7--

      end module m_depsitbw
