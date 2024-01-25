!***********************************************************************
      module m_nuc1stc
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2000/10/18, 2001/06/29, 2001/10/18,
!                   2001/11/14, 2001/11/20, 2001/12/11, 2002/01/15,
!                   2002/04/02, 2002/07/03, 2002/12/06, 2003/01/04,
!                   2003/03/21, 2003/03/28, 2003/04/30, 2003/05/19,
!                   2003/10/31, 2003/12/12, 2004/03/22, 2004/04/01,
!                   2004/04/15, 2004/08/01, 2004/09/01, 2004/09/10,
!                   2004/10/12, 2004/12/17, 2005/09/30, 2006/04/03,
!                   2007/10/19, 2007/11/26, 2008/01/11, 2008/05/02,
!                   2008/07/01, 2008/08/25, 2009/02/27, 2009/11/13,
!                   2011/03/18, 2012/06/19, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the nucleation rate of the condensation, contact and
!     homogeneous.

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

      public :: nuc1stc, s_nuc1stc

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface nuc1stc

        module procedure s_nuc1stc

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic log
      intrinsic min

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_nuc1stc(dtb,thresq,ni,nj,nk,p,qc,ncc,t,tcel,lv,kp,mu,&
     &                     diaqc,nuci)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: thresq
                       ! Minimum threshold value of mixing ratio

      real, intent(in) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure

      real, intent(in) :: qc(0:ni+1,0:nj+1,1:nk)
                       ! Cloud water mixing ratio

      real, intent(in) :: ncc(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud water

      real, intent(in) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

      real, intent(in) :: tcel(0:ni+1,0:nj+1,1:nk)
                       ! Ambient air temperature

      real, intent(in) :: lv(0:ni+1,0:nj+1,1:nk)
                       ! Latent heat of evapolation

      real, intent(in) :: kp(0:ni+1,0:nj+1,1:nk)
                       ! Thermal conductivity of air

      real, intent(in) :: mu(0:ni+1,0:nj+1,1:nk)
                       ! Viscosity of air

      real, intent(in) :: diaqc(0:ni+1,0:nj+1,1:nk)
                       ! Mean diameter of cloud water

! Output variable

      real, intent(out) :: nuci(0:ni+1,0:nj+1,1:nk)
                       ! Nucleation rate
                       ! of condensation, contact and homogeneous

! Internal shared variables

      real rwdt2       ! 100.0 / rhow x dtb

      real cc45        ! 200000.0 x 2.0 x cc

      real cknd        ! 2.0 x lmb20 x p20 / (t20 x diasl)

      real cdar        ! boltz / (3.0 x cc x diasl)

      real kpa25       ! 2.5 x kpa
      real kpa50       ! 5.0 x kpa

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real tc          ! Cloud water temperature

      real piv         ! Inverse of pressure

      real knd         ! Knudsen number

      real dar         ! Molecular diffusivity of aerosol

      real f1          ! cc45 x (270.16 - tc)^1.3 x diaqc
      real f2          ! kpa x piv x (t - tc)

      real ft          ! Function of Knudsen number

      real nufci       ! Nucleation rate of condensation
      real nucci       ! Nucleation rate of contact
      real nuhci       ! Nucleation rate of homogeneous

!-----7--------------------------------------------------------------7--

! Set the common used variables.

      rwdt2=(100.e0/rhow)*dtb

      cc45=4.e5*cc

      cknd=2.e0*lmb20*p20/(t20*diasl)

      cdar=boltz/(3.e0*cc*diasl)

      kpa25=2.5e0*kpa
      kpa50=5.e0*kpa

! -----

!!!! Calculate the nucleation rate of the condensation, contact and
!!!! homogeneous.

!$omp parallel default(shared) private(k)

!!! In the case nk = 1.

      if(nk.eq.1) then

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,tc,piv,knd,dar,f1,f2,ft,nufci,nucci,nuhci)

        do j=1,nj-1
        do i=1,ni-1

!! Calculate the nucleation rate in the case the cloud water is grater
!! than threshold value.

          if(qc(i,j,1).gt.thresq) then

! Calculate the nucleation rate of the condensation.

            if(t(i,j,1).gt.tlow.and.t(i,j,1).lt.t0) then

              nufci=min(rwdt2*(exp(-.66e0*tcel(i,j,1))-1.e0)            &
     &          *qc(i,j,1)*qc(i,j,1)/ncc(i,j,1),qc(i,j,1))

            else

              nufci=0.e0

            end if

! -----

! Calculate the nucleation rate of the contact.

            tc=t(i,j,1)

            if(tc.gt.tlow.and.tc.lt.270.16e0) then

              piv=1.e0/p(i,j,1)

              knd=cknd*piv*t(i,j,1)

              dar=cdar*tc*(1.e0+knd)/mu(i,j,1)

              f1=cc45*exp(1.3e0*log(270.16e0-tc))*diaqc(i,j,1)

              if(tc.lt.t(i,j,1)) then

                f2=kpa*piv*(t(i,j,1)-tc)

                ft=(kp(i,j,1)+kpa25*knd)                                &
     &            *(.4e0+.58e0*knd+.16e0*exp(-1.e0/knd))                &
     &            /((1.e0+3.e0*knd)*(2.e0*kp(i,j,1)+kpa50*knd+kpa))

                nucci=f1*f2*(rv*t(i,j,1)/lv(i,j,1)+ft)

              else

                nucci=0.e0

              end if

              nucci=(nucci+f1*dar)*qc(i,j,1)*dtb

            else

              nucci=0.e0

            end if

! -----

! Calculate the nucleation rate of the homogeneous.

            if(t(i,j,1).le.tlow) then

              nuhci=qc(i,j,1)

            else

              nuhci=0.e0

            end if

! -----

! Finally get the total nucleation rate.

            nuci(i,j,1)=nufci+nucci+nuhci

! -----

!! -----

! Fill in the array nuci with 0 in the case the cloud water is lower
! than threshold value.

          else

            nuci(i,j,1)=0.e0

          end if

! -----

        end do
        end do

!$omp end do

!!! -----

!!! In the case nk > 1.

      else

        do k=1,nk-1

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,tc,piv,knd,dar,f1,f2,ft,nufci,nucci,nuhci)

          do j=1,nj-1
          do i=1,ni-1

!! Calculate the nucleation rate in the case the cloud water is grater
!! than threshold value.

            if(qc(i,j,k).gt.thresq) then

! Calculate the nucleation rate of the condensation.

              if(t(i,j,k).gt.tlow.and.t(i,j,k).lt.t0) then

                nufci=min(rwdt2*(exp(-.66e0*tcel(i,j,k))-1.e0)          &
     &            *qc(i,j,k)*qc(i,j,k)/ncc(i,j,k),qc(i,j,k))

              else

                nufci=0.e0

              end if

! -----

! Calculate the nucleation rate of the contact.

              tc=t(i,j,k)

              if(tc.gt.tlow.and.tc.lt.270.16e0) then

                piv=1.e0/p(i,j,k)

                knd=cknd*piv*t(i,j,k)

                dar=cdar*tc*(1.e0+knd)/mu(i,j,k)

                f1=cc45*exp(1.3e0*log(270.16e0-tc))*diaqc(i,j,k)

                if(tc.lt.t(i,j,k)) then

                  f2=kpa*piv*(t(i,j,k)-tc)

                  ft=(kp(i,j,k)+kpa25*knd)                              &
     &              *(.4e0+.58e0*knd+.16e0*exp(-1.e0/knd))              &
     &              /((1.e0+3.e0*knd)*(2.e0*kp(i,j,k)+kpa50*knd+kpa))

                  nucci=f1*f2*(rv*t(i,j,k)/lv(i,j,k)+ft)

                else

                  nucci=0.e0

                end if

                nucci=(nucci+f1*dar)*qc(i,j,k)*dtb

              else

                nucci=0.e0

              end if

! -----

! Calculate the nucleation rate of the homogeneous.

              if(t(i,j,k).le.tlow) then

                nuhci=qc(i,j,k)

              else

                nuhci=0.e0

              end if

! -----

! Finally get the total nucleation rate.

              nuci(i,j,k)=nufci+nucci+nuhci

! -----

!! -----

! Fill in the array nuci with 0 in the case the cloud water is lower
! than threshold value.

            else

              nuci(i,j,k)=0.e0

            end if

! -----

          end do
          end do

!$omp end do

        end do

      end if

!!! -----

!$omp end parallel

!!!! -----

      end subroutine s_nuc1stc

!-----7--------------------------------------------------------------7--

      end module m_nuc1stc
