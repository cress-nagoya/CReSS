!***********************************************************************
      module m_setblk
!***********************************************************************

!     Author      : Sakakibara Atsushi, Naito Daisuke
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2000/10/18, 2000/10/27, 2000/11/17,
!                   2001/05/29, 2001/06/29, 2001/10/18, 2001/11/20,
!                   2001/12/11, 2002/01/15, 2002/04/02, 2003/03/28,
!                   2003/04/30, 2003/05/19, 2003/10/31, 2003/11/05,
!                   2003/12/12, 2004/03/22, 2004/04/01, 2004/04/15,
!                   2004/05/31, 2004/06/10, 2004/08/01, 2004/09/01,
!                   2004/09/10, 2004/09/25, 2004/10/12, 2004/12/17,
!                   2005/01/07, 2006/09/30, 2007/10/19, 2007/11/26,
!                   2008/01/11, 2008/05/02, 2008/07/01, 2008/08/25,
!                   2009/02/27, 2011/03/18, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the air temperature, saturation mixing ratio, latent
!     heat, thermal conductivity of air, viscosity of air, molecular
!     diffusivity of water, mean mass of cloud ice, mean diameters and
!     ventilation factors.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comphy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: setblk, s_setblk

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface setblk

        module procedure s_setblk

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic log
      intrinsic sqrt

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_setblk(cphopt,thresq,ni,nj,nk,ptbr,rbr,rbv,pi,p,ptp, &
     &                    qv,qc,qr,qi,qs,qg,ncc,ncr,nci,ncs,ncg,t,tcel, &
     &                    qvsst0,qvsw,qvsi,lv,ls,lf,kp,mu,dv,mi,        &
     &                    diaqc,diaqr,diaqi,diaqs,diaqg,                &
     &                    vntr,vnts,vntg,nu)
!***********************************************************************

! Input variables

      integer, intent(in) :: cphopt
                       ! Option for cloud micro physics

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: thresq
                       ! Minimum threshold value of mixing ratio

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rbv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of base state density

      real, intent(in) :: pi(0:ni+1,0:nj+1,1:nk)
                       ! Exner function

      real, intent(in) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure

      real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

      real, intent(in) :: qc(0:ni+1,0:nj+1,1:nk)
                       ! Cloud water mixing ratio

      real, intent(in) :: qr(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixing ratio

      real, intent(in) :: qi(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixing rati

      real, intent(in) :: qs(0:ni+1,0:nj+1,1:nk)
                       ! Snow mixing ratio

      real, intent(in) :: qg(0:ni+1,0:nj+1,1:nk)
                       ! Graupel mixing ratio

      real, intent(in) :: ncc(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud water

      real, intent(in) :: ncr(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of rain water

      real, intent(in) :: nci(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud ice

      real, intent(in) :: ncs(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of snow

      real, intent(in) :: ncg(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of graupel

! Output variables

      real, intent(out) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

      real, intent(out) :: tcel(0:ni+1,0:nj+1,1:nk)
                       ! Ambient air temperature

      real, intent(out) :: qvsst0(0:ni+1,0:nj+1,1:nk)
                       ! Super saturation mixing ratio at melting point

      real, intent(out) :: qvsw(0:ni+1,0:nj+1,1:nk)
                       ! Saturation mixing ratio for water

      real, intent(out) :: qvsi(0:ni+1,0:nj+1,1:nk)
                       ! Saturation mixing ratio for ice

      real, intent(out) :: lv(0:ni+1,0:nj+1,1:nk)
                       ! Latent heat of evapolation

      real, intent(out) :: ls(0:ni+1,0:nj+1,1:nk)
                       ! Latent heat of sublimation

      real, intent(out) :: lf(0:ni+1,0:nj+1,1:nk)
                       ! Latent heat of fusion

      real, intent(out) :: kp(0:ni+1,0:nj+1,1:nk)
                       ! Thermal conductivity of air

      real, intent(out) :: mu(0:ni+1,0:nj+1,1:nk)
                       ! Viscosity of air

      real, intent(out) :: dv(0:ni+1,0:nj+1,1:nk)
                       ! Molecular diffusivity of water

      real, intent(out) :: mi(0:ni+1,0:nj+1,1:nk)
                       ! Mean mass of cloud ice

      real, intent(out) :: diaqc(0:ni+1,0:nj+1,1:nk)
                       ! Mean diameter of cloud water

      real, intent(out) :: diaqr(0:ni+1,0:nj+1,1:nk)
                       ! Mean diameter of rain water

      real, intent(out) :: diaqi(0:ni+1,0:nj+1,1:nk)
                       ! Mean diameter of cloud ice

      real, intent(out) :: diaqs(0:ni+1,0:nj+1,1:nk)
                       ! Mean diameter of snow

      real, intent(out) :: diaqg(0:ni+1,0:nj+1,1:nk)
                       ! Mean diameter of graupel

      real, intent(out) :: vntr(0:ni+1,0:nj+1,1:nk)
                       ! Ventilation factor for rain water

      real, intent(out) :: vnts(0:ni+1,0:nj+1,1:nk)
                       ! Ventilation factor for snow

      real, intent(out) :: vntg(0:ni+1,0:nj+1,1:nk)
                       ! Ventilation factor for graupel

! Internal shared variables

      real cwmci       ! cw - ci

      real epses0      ! epsva x es0

      real t23iv       ! 1.0 / t23

      real cnu         ! nu0^0.57013 / t0
      real cdv         ! dv0^0.55249 / t0

      real pdiaqr      ! (5.0 + bur) / 2.0
      real pdiaqs      ! (5.0 + bus) / 2.0 - 1.0
      real pdiaqg      ! (5.0 + bug) / 2.0 - 1.0

      real ccrw6       ! 6.0 / (cc x rhow)
      real ccri6       ! 6.0 / (cc x rhoi)

      real cdiaqc      ! Coefficient of mean diameter of cloud water
      real cdiaqr      ! Coefficient of mean diameter of rain water
      real cdiaqs      ! Coefficient of mean diameter of snow
      real cdiaqg      ! Coefficient of mean diameter of graupel

      real cvntr       ! Coefficient of ventilation facter
                       ! for rain water

      real cvnts       ! Coefficient of ventilation facter
                       ! for snow

      real cvntg       ! Coefficient of ventilation facter
                       ! for graupel

      real, intent(inout) :: nu(0:ni+1,0:nj+1)
                       ! Kinetic viscosity of air

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real esw         ! Saturation vapor pressure for water
      real esi         ! Saturation vapor pressure for ice

      real p20dvp      ! p20 / p

      real cvnt        ! Coefficient of ventilation facter

!-----7--------------------------------------------------------------7--

! Set the common used variables.

      cwmci=cw-ci

      epses0=epsva*es0

      t23iv=1.e0/t23

      cdv=1.e0/t0

      cnu=exp(.57013e0*log(nu0))*cdv
      cdv=exp(.55249e0*log(dv0))*cdv

      pdiaqr=.5e0*(5.e0+bur)-1.e0
      pdiaqs=.5e0*(5.e0+bus)-1.e0
      pdiaqg=.5e0*(5.e0+bug)-1.e0

      ccrw6=6.e0/(cc*rhow)
      ccri6=6.e0/(cc*rhoi)

      cdiaqc=1.e0/(cc*rhow)
      cdiaqr=1.e0/(cc*rhow)
      cdiaqs=1.e0/(cc*rhos)
      cdiaqg=1.e0/(cc*rhog)

      cvntr=.31e0*sqrt(aur)*exp(oned3*log(sc))*gf5bur
      cvnts=.31e0*sqrt(aus)*exp(oned3*log(sc))*gf5bus
      cvntg=.31e0*sqrt(aug)*exp(oned3*log(sc))*gf5bug

! -----

!!!!! Calculate the air temperature, saturation mixing ratio, latent
!!!!! heat, thermal conductivity of air, viscosity of air, molecular
!!!!! diffusivity of water, mean mass of cloud ice, mean diameters
!!!!! and ventilation factors.

!$omp parallel default(shared) private(k)

!!!! In the case nk = 1.

      if(nk.eq.1) then

!!! Perform calculating in the case the option abs(cphopt) is less
!!! than 3.

        if(abs(cphopt).le.3) then

!! Calculate the air temperature, saturation mixing ratio, latent heat,
!! thermal conductivity of air, viscosity of air and molecular
!! diffusivity of water.

!$omp do schedule(runtime) private(i,j,esw,esi,p20dvp)

          do j=1,nj-1
          do i=1,ni-1

! Calculate the air temperature.

            t(i,j,1)=(ptbr(i,j,1)+ptp(i,j,1))*pi(i,j,1)

            tcel(i,j,1)=t(i,j,1)-t0

! -----

! Calculate the saturation mixing ratio.

            esw=es0*exp(17.269e0*tcel(i,j,1)/(t(i,j,1)-35.86e0))
            esi=es0*exp(21.875e0*tcel(i,j,1)/(t(i,j,1)-7.66e0))

            qvsst0(i,j,1)=qv(i,j,1)-epses0/(p(i,j,1)-es0)

            qvsw(i,j,1)=epsva*esw/(p(i,j,1)-esw)
            qvsi(i,j,1)=epsva*esi/(p(i,j,1)-esi)

! -----

! Calculate the latent heat of evapolation, sublimation and fusion.

            lv(i,j,1)=lv0                                               &
     &        *exp((.167e0+3.67e-4*t(i,j,1))*log(t0/t(i,j,1)))

            lf(i,j,1)=lf0+cwmci*tcel(i,j,1)

            ls(i,j,1)=lv(i,j,1)+lf(i,j,1)

! -----

! Calculate the viscosity and kinetic viscosity of air and molecular
! diffusivity of water vapor.

            p20dvp=p20/p(i,j,1)

            nu(i,j)=p20dvp*exp(1.754e0*log(cnu*t(i,j,1)))

            kp(i,j,1)                                                   &
     &        =(kp0/(120.e0+t(i,j,1)))*exp(1.5e0*log(t23iv*t(i,j,1)))

            mu(i,j,1)=rbr(i,j,1)*nu(i,j)

            dv(i,j,1)=p20dvp*exp(1.81e0*log(cdv*t(i,j,1)))

! -----

          end do
          end do

!$omp end do

!! -----

! Calculate the mean mass of cloud ice, mean diameters and ventilation
! factors.

!$omp do schedule(runtime) private(i,j,cvnt)

          do j=1,nj-1
          do i=1,ni-1
            cvnt=sqrt(sqrt(r0*rbv(i,j,1))/nu(i,j))

            if(qc(i,j,1).gt.thresq) then

              diaqc(i,j,1)=exp(oned3*log(ccrw6*qc(i,j,1)/ncc(i,j,1)))

            else

              diaqc(i,j,1)=0.e0

            end if

            if(qr(i,j,1).gt.thresq) then

              diaqr(i,j,1)=exp(oned3*log(cdiaqr*qr(i,j,1)/ncr(i,j,1)))

              vntr(i,j,1)=rbr(i,j,1)*ncr(i,j,1)*(.78e0*diaqr(i,j,1)     &
     &          +cvnt*cvntr*exp(pdiaqr*log(diaqr(i,j,1))))

            else

              diaqr(i,j,1)=0.e0

              vntr(i,j,1)=0.e0

            end if

            if(qi(i,j,1).gt.thresq) then

              mi(i,j,1)=qi(i,j,1)/nci(i,j,1)

              diaqi(i,j,1)=exp(oned3*log(ccri6*mi(i,j,1)))

            else

              mi(i,j,1)=0.e0

              diaqi(i,j,1)=0.e0

            end if

            if(qs(i,j,1).gt.thresq) then

              diaqs(i,j,1)=exp(oned3*log(cdiaqs*qs(i,j,1)/ncs(i,j,1)))

              vnts(i,j,1)=rbr(i,j,1)*ncs(i,j,1)*(.78e0*diaqs(i,j,1)     &
     &          +cvnt*cvnts*exp(pdiaqs*log(diaqs(i,j,1))))

            else

              diaqs(i,j,1)=0.e0

              vnts(i,j,1)=0.e0

            end if

            if(qg(i,j,1).gt.thresq) then

              diaqg(i,j,1)=exp(oned3*log(cdiaqg*qg(i,j,1)/ncg(i,j,1)))

              vntg(i,j,1)=rbr(i,j,1)*ncg(i,j,1)*(.78e0*diaqg(i,j,1)     &
     &          +cvnt*cvntg*exp(pdiaqg*log(diaqg(i,j,1))))

            else

              diaqg(i,j,1)=0.e0

              vntg(i,j,1)=0.e0

            end if

          end do
          end do

!$omp end do

! -----

!!! -----

!!! Perform calculating in the case the option abs(cphopt) is equal
!!! to 4.

        else

!! Calculate the air temperature, saturation mixing ratio, latent heat,
!! thermal conductivity of air, viscosity of air and molecular
!! diffusivity of water.

!$omp do schedule(runtime) private(i,j,esw,esi,p20dvp)

          do j=1,nj-1
          do i=1,ni-1

! Calculate the air temperature.

            t(i,j,1)=(ptbr(i,j,1)+ptp(i,j,1))*pi(i,j,1)

            tcel(i,j,1)=t(i,j,1)-t0

! -----

! Calculate the saturation mixing ratio.

            esw=es0*exp(17.269e0*tcel(i,j,1)/(t(i,j,1)-35.86e0))
            esi=es0*exp(21.875e0*tcel(i,j,1)/(t(i,j,1)-7.66e0))

            qvsst0(i,j,1)=qv(i,j,1)-epses0/(p(i,j,1)-es0)

            qvsw(i,j,1)=epsva*esw/(p(i,j,1)-esw)
            qvsi(i,j,1)=epsva*esi/(p(i,j,1)-esi)

! -----

! Calculate the latent heat of evapolation, sublimation and fusion.

            lv(i,j,1)=lv0                                               &
     &        *exp((.167e0+3.67e-4*t(i,j,1))*log(t0/t(i,j,1)))

            lf(i,j,1)=lf0+cwmci*tcel(i,j,1)

            ls(i,j,1)=lv(i,j,1)+lf(i,j,1)

! -----

! Calculate the viscosity and kinetic viscosity of air and molecular
! diffusivity of water vapor.

            p20dvp=p20/p(i,j,1)

            nu(i,j)=p20dvp*exp(1.754e0*log(cnu*t(i,j,1)))

            kp(i,j,1)                                                   &
     &        =(kp0/(120.e0+t(i,j,1)))*exp(1.5e0*log(t23iv*t(i,j,1)))

            mu(i,j,1)=rbr(i,j,1)*nu(i,j)

            dv(i,j,1)=p20dvp*exp(1.81e0*log(cdv*t(i,j,1)))

! -----

          end do
          end do

!$omp end do

!! -----

! Calculate the mean mass of cloud ice, mean diameters and ventilation
! factors.

!$omp do schedule(runtime) private(i,j,cvnt)

          do j=1,nj-1
          do i=1,ni-1
            cvnt=sqrt(sqrt(r0*rbv(i,j,1))/nu(i,j))

            if(qc(i,j,1).gt.thresq) then

              diaqc(i,j,1)=exp(oned3*log(cdiaqc*qc(i,j,1)/ncc(i,j,1)))

            else

              diaqc(i,j,1)=0.e0

            end if

            if(qr(i,j,1).gt.thresq) then

              diaqr(i,j,1)=exp(oned3*log(cdiaqr*qr(i,j,1)/ncr(i,j,1)))

              vntr(i,j,1)=rbr(i,j,1)*ncr(i,j,1)*(.78e0*diaqr(i,j,1)     &
     &          +cvnt*cvntr*exp(pdiaqr*log(diaqr(i,j,1))))

            else

              diaqr(i,j,1)=0.e0

              vntr(i,j,1)=0.e0

            end if

            if(qi(i,j,1).gt.thresq) then

              mi(i,j,1)=qi(i,j,1)/nci(i,j,1)

              diaqi(i,j,1)=exp(oned3*log(ccri6*mi(i,j,1)))

            else

              mi(i,j,1)=0.e0

              diaqi(i,j,1)=0.e0

            end if

            if(qs(i,j,1).gt.thresq) then

              diaqs(i,j,1)=exp(oned3*log(cdiaqs*qs(i,j,1)/ncs(i,j,1)))

              vnts(i,j,1)=rbr(i,j,1)*ncs(i,j,1)*(.78e0*diaqs(i,j,1)     &
     &          +cvnt*cvnts*exp(pdiaqs*log(diaqs(i,j,1))))

            else

              diaqs(i,j,1)=0.e0

              vnts(i,j,1)=0.e0

            end if

            if(qg(i,j,1).gt.thresq) then

              diaqg(i,j,1)=exp(oned3*log(cdiaqg*qg(i,j,1)/ncg(i,j,1)))

              vntg(i,j,1)=rbr(i,j,1)*ncg(i,j,1)*(.78e0*diaqg(i,j,1)     &
     &          +cvnt*cvntg*exp(pdiaqg*log(diaqg(i,j,1))))

            else

              diaqg(i,j,1)=0.e0

              vntg(i,j,1)=0.e0

            end if

          end do
          end do

!$omp end do

! -----

        end if

!!! -----

!!!! -----

!!!! In the case nk > 1.

      else

!!! Perform calculating in the case the option abs(cphopt) is less
!!! than 3.

        if(abs(cphopt).le.3) then

          do k=1,nk-1

!! Calculate the air temperature, saturation mixing ratio, latent heat,
!! thermal conductivity of air, viscosity of air and molecular
!! diffusivity of water.

!$omp do schedule(runtime) private(i,j,esw,esi,p20dvp)

            do j=1,nj-1
            do i=1,ni-1

! Calculate the air temperature.

              t(i,j,k)=(ptbr(i,j,k)+ptp(i,j,k))*pi(i,j,k)

              tcel(i,j,k)=t(i,j,k)-t0

! -----

! Calculate the saturation mixing ratio.

              esw=es0*exp(17.269e0*tcel(i,j,k)/(t(i,j,k)-35.86e0))
              esi=es0*exp(21.875e0*tcel(i,j,k)/(t(i,j,k)-7.66e0))

              qvsst0(i,j,k)=qv(i,j,k)-epses0/(p(i,j,k)-es0)

              qvsw(i,j,k)=epsva*esw/(p(i,j,k)-esw)
              qvsi(i,j,k)=epsva*esi/(p(i,j,k)-esi)

! -----

! Calculate the latent heat of evapolation, sublimation and fusion.

              lv(i,j,k)=lv0                                             &
     &          *exp((.167e0+3.67e-4*t(i,j,k))*log(t0/t(i,j,k)))

              lf(i,j,k)=lf0+cwmci*tcel(i,j,k)

              ls(i,j,k)=lv(i,j,k)+lf(i,j,k)

! -----

! Calculate the viscosity and kinetic viscosity of air and molecular
! diffusivity of water vapor.

              p20dvp=p20/p(i,j,k)

              nu(i,j)=p20dvp*exp(1.754e0*log(cnu*t(i,j,k)))

              kp(i,j,k)                                                 &
     &          =(kp0/(120.e0+t(i,j,k)))*exp(1.5e0*log(t23iv*t(i,j,k)))

              mu(i,j,k)=rbr(i,j,k)*nu(i,j)

              dv(i,j,k)=p20dvp*exp(1.81e0*log(cdv*t(i,j,k)))

! -----

            end do
            end do

!$omp end do

!! -----

! Calculate the mean mass of cloud ice, mean diameters and ventilation
! factors.

!$omp do schedule(runtime) private(i,j,cvnt)

            do j=1,nj-1
            do i=1,ni-1
              cvnt=sqrt(sqrt(r0*rbv(i,j,k))/nu(i,j))

              if(qc(i,j,k).gt.thresq) then

                diaqc(i,j,k)=exp(oned3*log(ccrw6*qc(i,j,k)/ncc(i,j,k)))

              else

                diaqc(i,j,k)=0.e0

              end if

              if(qr(i,j,k).gt.thresq) then

                diaqr(i,j,k)=exp(oned3*log(cdiaqr*qr(i,j,k)/ncr(i,j,k)))

                vntr(i,j,k)=rbr(i,j,k)*ncr(i,j,k)*(.78e0*diaqr(i,j,k)   &
     &            +cvnt*cvntr*exp(pdiaqr*log(diaqr(i,j,k))))

              else

                diaqr(i,j,k)=0.e0

                vntr(i,j,k)=0.e0

              end if

              if(qi(i,j,k).gt.thresq) then

                mi(i,j,k)=qi(i,j,k)/nci(i,j,k)

                diaqi(i,j,k)=exp(oned3*log(ccri6*mi(i,j,k)))

              else

                mi(i,j,k)=0.e0

                diaqi(i,j,k)=0.e0

              end if

              if(qs(i,j,k).gt.thresq) then

                diaqs(i,j,k)=exp(oned3*log(cdiaqs*qs(i,j,k)/ncs(i,j,k)))

                vnts(i,j,k)=rbr(i,j,k)*ncs(i,j,k)*(.78e0*diaqs(i,j,k)   &
     &            +cvnt*cvnts*exp(pdiaqs*log(diaqs(i,j,k))))

              else

                diaqs(i,j,k)=0.e0

                vnts(i,j,k)=0.e0

              end if

              if(qg(i,j,k).gt.thresq) then

                diaqg(i,j,k)=exp(oned3*log(cdiaqg*qg(i,j,k)/ncg(i,j,k)))

                vntg(i,j,k)=rbr(i,j,k)*ncg(i,j,k)*(.78e0*diaqg(i,j,k)   &
     &            +cvnt*cvntg*exp(pdiaqg*log(diaqg(i,j,k))))

              else

                diaqg(i,j,k)=0.e0

                vntg(i,j,k)=0.e0

              end if

            end do
            end do

!$omp end do

! -----

          end do

!!! -----

!!! Perform calculating in the case the option abs(cphopt) is equal
!!! to 4.

        else

          do k=1,nk-1

!! Calculate the air temperature, saturation mixing ratio, latent heat,
!! thermal conductivity of air, viscosity of air and molecular
!! diffusivity of water.

!$omp do schedule(runtime) private(i,j,esw,esi,p20dvp)

            do j=1,nj-1
            do i=1,ni-1

! Calculate the air temperature.

              t(i,j,k)=(ptbr(i,j,k)+ptp(i,j,k))*pi(i,j,k)

              tcel(i,j,k)=t(i,j,k)-t0

! -----

! Calculate the saturation mixing ratio.

              esw=es0*exp(17.269e0*tcel(i,j,k)/(t(i,j,k)-35.86e0))
              esi=es0*exp(21.875e0*tcel(i,j,k)/(t(i,j,k)-7.66e0))

              qvsst0(i,j,k)=qv(i,j,k)-epses0/(p(i,j,k)-es0)

              qvsw(i,j,k)=epsva*esw/(p(i,j,k)-esw)
              qvsi(i,j,k)=epsva*esi/(p(i,j,k)-esi)

! -----

! Calculate the latent heat of evapolation, sublimation and fusion.

              lv(i,j,k)=lv0                                             &
     &          *exp((.167e0+3.67e-4*t(i,j,k))*log(t0/t(i,j,k)))

              lf(i,j,k)=lf0+cwmci*tcel(i,j,k)

              ls(i,j,k)=lv(i,j,k)+lf(i,j,k)

! -----

! Calculate the viscosity and kinetic viscosity of air and molecular
! diffusivity of water vapor.

              p20dvp=p20/p(i,j,k)

              nu(i,j)=p20dvp*exp(1.754e0*log(cnu*t(i,j,k)))

              kp(i,j,k)                                                 &
     &          =(kp0/(120.e0+t(i,j,k)))*exp(1.5e0*log(t23iv*t(i,j,k)))

              mu(i,j,k)=rbr(i,j,k)*nu(i,j)

              dv(i,j,k)=p20dvp*exp(1.81e0*log(cdv*t(i,j,k)))

! -----

            end do
            end do

!$omp end do

!! -----

! Calculate the mean mass of cloud ice, mean diameters and ventilation
! factors.

!$omp do schedule(runtime) private(i,j,cvnt)

            do j=1,nj-1
            do i=1,ni-1
              cvnt=sqrt(sqrt(r0*rbv(i,j,k))/nu(i,j))

              if(qc(i,j,k).gt.thresq) then

                diaqc(i,j,k)=exp(oned3*log(cdiaqc*qc(i,j,k)/ncc(i,j,k)))

              else

                diaqc(i,j,k)=0.e0

              end if

              if(qr(i,j,k).gt.thresq) then

                diaqr(i,j,k)=exp(oned3*log(cdiaqr*qr(i,j,k)/ncr(i,j,k)))

                vntr(i,j,k)=rbr(i,j,k)*ncr(i,j,k)*(.78e0*diaqr(i,j,k)   &
     &            +cvnt*cvntr*exp(pdiaqr*log(diaqr(i,j,k))))

              else

                diaqr(i,j,k)=0.e0

                vntr(i,j,k)=0.e0

              end if

              if(qi(i,j,k).gt.thresq) then

                mi(i,j,k)=qi(i,j,k)/nci(i,j,k)

                diaqi(i,j,k)=exp(oned3*log(ccri6*mi(i,j,k)))

              else

                mi(i,j,k)=0.e0

                diaqi(i,j,k)=0.e0

              end if

              if(qs(i,j,k).gt.thresq) then

                diaqs(i,j,k)=exp(oned3*log(cdiaqs*qs(i,j,k)/ncs(i,j,k)))

                vnts(i,j,k)=rbr(i,j,k)*ncs(i,j,k)*(.78e0*diaqs(i,j,k)   &
     &            +cvnt*cvnts*exp(pdiaqs*log(diaqs(i,j,k))))

              else

                diaqs(i,j,k)=0.e0

                vnts(i,j,k)=0.e0

              end if

              if(qg(i,j,k).gt.thresq) then

                diaqg(i,j,k)=exp(oned3*log(cdiaqg*qg(i,j,k)/ncg(i,j,k)))

                vntg(i,j,k)=rbr(i,j,k)*ncg(i,j,k)*(.78e0*diaqg(i,j,k)   &
     &            +cvnt*cvntg*exp(pdiaqg*log(diaqg(i,j,k))))

              else

                diaqg(i,j,k)=0.e0

                vntg(i,j,k)=0.e0

              end if

            end do
            end do

!$omp end do

! -----

          end do

        end if

!!! -----

      end if

!!!! -----

!$omp end parallel

!!!!! -----

      end subroutine s_setblk

!-----7--------------------------------------------------------------7--

      end module m_setblk
