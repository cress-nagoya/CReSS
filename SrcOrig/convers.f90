!***********************************************************************
      module m_convers
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2000/10/18, 2000/11/17, 2001/05/29,
!                   2001/06/29, 2001/10/18, 2001/11/20, 2001/12/11,
!                   2002/01/15, 2002/04/02, 2002/07/03, 2002/12/06,
!                   2003/03/13, 2003/03/28, 2003/04/30, 2003/05/19,
!                   2003/10/31, 2003/12/12, 2004/03/22, 2004/04/01,
!                   2004/04/15, 2004/05/31, 2004/06/10, 2004/07/01,
!                   2004/08/01, 2004/09/01, 2004/09/10, 2004/09/25,
!                   2004/10/12, 2004/12/17, 2005/01/07, 2005/04/04,
!                   2005/08/05, 2005/09/30, 2006/04/03, 2007/07/30,
!                   2007/10/19, 2007/11/26, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2009/11/13, 2011/03/18, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the conversion rate from the cloud water to the rain
!     water, from the cloud ice to the snow and from the snow to the
!     graupel.

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

      public :: convers, s_convers

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface convers

        module procedure s_convers

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic exp
      intrinsic log
      intrinsic min
      intrinsic sqrt

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_convers(cphopt,dtb,thresq,ni,nj,nk,rbr,rbv,          &
     &                     qc,qi,qs,ncc,ncs,t,mu,mi,diaqi,diaqs,        &
     &                     clcs,vdvi,vdvs,ecs,cncr,cnis,cnsg,cnsgn)
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

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: thresq
                       ! Minimum threshold value of mixing ratio

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rbv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of base state density

      real, intent(in) :: qc(0:ni+1,0:nj+1,1:nk)
                       ! Cloud water mixing ratio

      real, intent(in) :: qi(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixing ratio

      real, intent(in) :: qs(0:ni+1,0:nj+1,1:nk)
                       ! Snow mixing ratio

      real, intent(in) :: ncc(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud water

      real, intent(in) :: ncs(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of snow

      real, intent(in) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

      real, intent(in) :: mu(0:ni+1,0:nj+1,1:nk)
                       ! Viscosity of air

      real, intent(in) :: mi(0:ni+1,0:nj+1,1:nk)
                       ! Mean mass of cloud ice

      real, intent(in) :: diaqi(0:ni+1,0:nj+1,1:nk)
                       ! Mean diameter of cloud ice

      real, intent(in) :: diaqs(0:ni+1,0:nj+1,1:nk)
                       ! Mean diameter of snow

      real, intent(in) :: clcs(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate from cloud water to snow

      real, intent(in) :: vdvi(0:ni+1,0:nj+1,1:nk)
                       ! Deposiotn rate from water vapor to cloud ice

      real, intent(in) :: vdvs(0:ni+1,0:nj+1,1:nk)
                       ! Deposiotn rate from water vapor to snow

      real, intent(in) :: ecs(0:ni+1,0:nj+1,1:nk)
                       ! Collection efficiency of snow for cloud water

! Output variables

      real, intent(out) :: cncr(0:ni+1,0:nj+1,1:nk)
                       ! Conversion rate from cloud water to rain water

      real, intent(out) :: cnis(0:ni+1,0:nj+1,1:nk)
                       ! Conversion rate from cloud ice to snow

      real, intent(out) :: cnsg(0:ni+1,0:nj+1,1:nk)
                       ! Conversion rate from snow to graupel

      real, intent(out) :: cnsgn(0:ni+1,0:nj+1,1:nk)
                       ! Conversion rate for concentrations
                       ! from snow to graupel

! Internal shared variables

      real busm1       ! bus - 1.0

      real ms05        ! 0.5 x ms0

      real qccm        ! Critical concentrations of cloud water

      real diaqs0      ! Diameter of minimum snow

      real cagin       ! Coefficient of aggregation rate for cloud ice

      real ccncr       ! Coefficient of conversion rate
      real ccnsg       ! Coefficient of conversion rate

      real ccnsgn      ! Coefficient of conversion rate
                       ! for concentrations

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real qcm         ! Critical mixing ratio of cloud water

!-----7--------------------------------------------------------------7--

! Set the common used variables.

      busm1=bus-1.e0

      ms05=.5e0*ms0

      qccm=oned6*rhow*cc*diaqcm*diaqcm*diaqcm

      diaqs0=exp(oned3*log(6.e0*ms0/(cc*rhos)))

      cagin=r0*exp(3.e0*log(.5e0*oned3*aui*eii*xl/rhoi*dtb))

      ccncr=.104e0*g*ecc*exp(-oned3*log(rhow))*dtb

      ccnsgn=1.e0/(rhog-rhos)

      ccnsg=rhog*ccnsgn
      ccnsgn=1.5e0*aus*gfbus*sqrt(r0)*ccnsgn*dtb

! -----

!!!!! Calculate the conversion rate.

!$omp parallel default(shared) private(k)

!!!! In the case nk = 1.

      if(nk.eq.1) then

!!! Perform calculating in the case the option abs(cphopt) is equal
!!! to 2.

        if(abs(cphopt).eq.2) then

!$omp do schedule(runtime) private(i,j,qcm)

          do j=1,nj-1
          do i=1,ni-1

! Calculate the conversion rate from the cloud water to the rain water.

            if(qc(i,j,1).gt.thresq) then

              qcm=qccm*ncc(i,j,1)

              if(t(i,j,1).gt.tlow.and.qc(i,j,1).ge.qcm) then

                cncr(i,j,1)=ccncr*rbr(i,j,1)*exp(-oned3*log(ncc(i,j,1)))&
     &            *exp(sevnd3*log(qc(i,j,1)))/mu(i,j,1)

              else

                cncr(i,j,1)=0.e0

              end if

            else

              cncr(i,j,1)=0.e0

            end if

! -----

!! Calculate the conversion rate from the cloud ice to the snow and the
!! snow to the graupel.

            if(t(i,j,1).lt.t0) then

! Calculate the conversion rate from the cloud ice to the snow.

              if(qi(i,j,1).gt.thresq) then

                cnis(i,j,1)=rbr(i,j,1)*exp(oned3*log(cagin*rbv(i,j,1))) &
     &            *qi(i,j,1)*qi(i,j,1)/log(diaqs0/diaqi(i,j,1))

                if(vdvi(i,j,1).ge.0.e0) then

                  if(mi(i,j,1).lt.ms05) then

                    cnis(i,j,1)=cnis(i,j,1)                             &
     &                +mi(i,j,1)/(ms0-mi(i,j,1))*vdvi(i,j,1)

                  else

                    cnis(i,j,1)=cnis(i,j,1)                             &
     &                +(vdvi(i,j,1)+(1.e0-ms05/mi(i,j,1))*qi(i,j,1))

                  end if

                end if

                cnis(i,j,1)=min(cnis(i,j,1),qi(i,j,1))

              else

                cnis(i,j,1)=0.e0

              end if

! -----

! Calculate the conversion rate from the snow to the graupel.

              if(qs(i,j,1).gt.thresq) then

                if(vdvs(i,j,1).gt.0.e0                                  &
     &            .and.vdvs(i,j,1).lt.clcs(i,j,1)) then

                  cnsg(i,j,1)=min(ccnsg*clcs(i,j,1),qs(i,j,1))

                else

                  cnsg(i,j,1)=0.e0

                end if

              else

                cnsg(i,j,1)=0.e0

              end if

! -----

!! -----

! Fill in the array with 0 in the case the air temperature is greater
! than melting points.

            else

              cnis(i,j,1)=0.e0
              cnsg(i,j,1)=0.e0

            end if

! -----

          end do
          end do

!$omp end do

!!! -----

!!! Perform calculating in the case the option abs(cphopt) is greater
!!! than 2.

        else if(abs(cphopt).ge.3) then

!$omp do schedule(runtime) private(i,j,qcm)

          do j=1,nj-1
          do i=1,ni-1

! Calculate the conversion rate from the cloud water to the rain water.

            if(qc(i,j,1).gt.thresq) then

              qcm=qccm*ncc(i,j,1)

              if(t(i,j,1).gt.tlow.and.qc(i,j,1).ge.qcm) then

                cncr(i,j,1)=ccncr*rbr(i,j,1)*exp(-oned3*log(ncc(i,j,1)))&
     &            *exp(sevnd3*log(qc(i,j,1)))/mu(i,j,1)

              else

                cncr(i,j,1)=0.e0

              end if

            else

              cncr(i,j,1)=0.e0

            end if

! -----

!! Calculate the conversion rate from the cloud ice to the snow and the
!! snow to the graupel.

            if(t(i,j,1).lt.t0) then

! Calculate the conversion rate from the cloud ice to the snow.

              if(qi(i,j,1).gt.thresq) then

                cnis(i,j,1)=rbr(i,j,1)*exp(oned3*log(cagin*rbv(i,j,1))) &
     &            *qi(i,j,1)*qi(i,j,1)/log(diaqs0/diaqi(i,j,1))

                if(vdvi(i,j,1).ge.0.e0) then

                  if(mi(i,j,1).lt.ms05) then

                    cnis(i,j,1)=cnis(i,j,1)                             &
     &                +mi(i,j,1)/(ms0-mi(i,j,1))*vdvi(i,j,1)

                  else

                    cnis(i,j,1)=cnis(i,j,1)                             &
     &                +(vdvi(i,j,1)+(1.e0-ms05/mi(i,j,1))*qi(i,j,1))

                  end if

                end if

                cnis(i,j,1)=min(cnis(i,j,1),qi(i,j,1))

              else

                cnis(i,j,1)=0.e0

              end if

! -----

! Calculate the conversion rate from the snow to the graupel.

              if(qs(i,j,1).gt.thresq) then

                if(vdvs(i,j,1).gt.0.e0                                  &
     &            .and.vdvs(i,j,1).lt.clcs(i,j,1)) then

                  cnsg(i,j,1)=min(ccnsg*clcs(i,j,1),qs(i,j,1))

                  cnsgn(i,j,1)=min(ccnsgn*ecs(i,j,1)                    &
     &              *qc(i,j,1)*ncs(i,j,1)*sqrt(rbr(i,j,1))              &
     &              *exp(busm1*log(diaqs(i,j,1))),ncs(i,j,1))

                else

                  cnsg(i,j,1)=0.e0
                  cnsgn(i,j,1)=0.e0

                end if

              else

                cnsg(i,j,1)=0.e0
                cnsgn(i,j,1)=0.e0

              end if

! -----

!! -----

! Fill in the array with 0 in the case the air temperature is greater
! than melting points.

            else

              cnis(i,j,1)=0.e0

              cnsg(i,j,1)=0.e0
              cnsgn(i,j,1)=0.e0

            end if

! -----

          end do
          end do

!$omp end do

        end if

!!! -----

!!!! -----

!!!! In the case nk > 1.

      else

!!! Perform calculating in the case the option abs(cphopt) is equal
!!! to 2.

        if(abs(cphopt).eq.2) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j,qcm)

            do j=1,nj-1
            do i=1,ni-1

! Calculate the conversion rate from the cloud water to the rain water.

              if(qc(i,j,k).gt.thresq) then

                qcm=qccm*ncc(i,j,k)

                if(t(i,j,k).gt.tlow.and.qc(i,j,k).ge.qcm) then

                  cncr(i,j,k)=ccncr*rbr(i,j,k)                          &
     &              *exp(-oned3*log(ncc(i,j,k)))                        &
     &              *exp(sevnd3*log(qc(i,j,k)))/mu(i,j,k)

                else

                  cncr(i,j,k)=0.e0

                end if

              else

                cncr(i,j,k)=0.e0

              end if

! -----

!! Calculate the conversion rate from the cloud ice to the snow and the
!! snow to the graupel.

              if(t(i,j,k).lt.t0) then

! Calculate the conversion rate from the cloud ice to the snow.

                if(qi(i,j,k).gt.thresq) then

                  cnis(i,j,k)=rbr(i,j,k)                                &
     &              *exp(oned3*log(cagin*rbv(i,j,k)))                   &
     &              *qi(i,j,k)*qi(i,j,k)/log(diaqs0/diaqi(i,j,k))

                  if(vdvi(i,j,k).ge.0.e0) then

                    if(mi(i,j,k).lt.ms05) then

                      cnis(i,j,k)=cnis(i,j,k)                           &
     &                  +mi(i,j,k)/(ms0-mi(i,j,k))*vdvi(i,j,k)

                    else

                      cnis(i,j,k)=cnis(i,j,k)                           &
     &                  +(vdvi(i,j,k)+(1.e0-ms05/mi(i,j,k))*qi(i,j,k))

                    end if

                  end if

                  cnis(i,j,k)=min(cnis(i,j,k),qi(i,j,k))

                else

                  cnis(i,j,k)=0.e0

                end if

! -----

! Calculate the conversion rate from the snow to the graupel.

                if(qs(i,j,k).gt.thresq) then

                  if(vdvs(i,j,k).gt.0.e0                                &
     &              .and.vdvs(i,j,k).lt.clcs(i,j,k)) then

                    cnsg(i,j,k)=min(ccnsg*clcs(i,j,k),qs(i,j,k))

                  else

                    cnsg(i,j,k)=0.e0

                  end if

                else

                  cnsg(i,j,k)=0.e0

                end if

! -----

!! -----

! Fill in the array with 0 in the case the air temperature is greater
! than melting points.

              else

                cnis(i,j,k)=0.e0
                cnsg(i,j,k)=0.e0

              end if

! -----

            end do
            end do

!$omp end do

          end do

!!! -----

!!! Perform calculating in the case the option abs(cphopt) is greater
!!! than 2.

        else if(abs(cphopt).ge.3) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j,qcm)

            do j=1,nj-1
            do i=1,ni-1

! Calculate the conversion rate from the cloud water to the rain water.

              if(qc(i,j,k).gt.thresq) then

                qcm=qccm*ncc(i,j,k)

                if(t(i,j,k).gt.tlow.and.qc(i,j,k).ge.qcm) then

                  cncr(i,j,k)=ccncr*rbr(i,j,k)                          &
     &              *exp(-oned3*log(ncc(i,j,k)))                        &
     &              *exp(sevnd3*log(qc(i,j,k)))/mu(i,j,k)

                else

                  cncr(i,j,k)=0.e0

                end if

              else

                cncr(i,j,k)=0.e0

              end if

! -----

!! Calculate the conversion rate from the cloud ice to the snow and the
!! snow to the graupel.

              if(t(i,j,k).lt.t0) then

! Calculate the conversion rate from the cloud ice to the snow.

                if(qi(i,j,k).gt.thresq) then

                  cnis(i,j,k)=rbr(i,j,k)                                &
     &              *exp(oned3*log(cagin*rbv(i,j,k)))                   &
     &              *qi(i,j,k)*qi(i,j,k)/log(diaqs0/diaqi(i,j,k))

                  if(vdvi(i,j,k).ge.0.e0) then

                    if(mi(i,j,k).lt.ms05) then

                      cnis(i,j,k)=cnis(i,j,k)                           &
     &                  +mi(i,j,k)/(ms0-mi(i,j,k))*vdvi(i,j,k)

                    else

                      cnis(i,j,k)=cnis(i,j,k)                           &
     &                  +(vdvi(i,j,k)+(1.e0-ms05/mi(i,j,k))*qi(i,j,k))

                    end if

                  end if

                  cnis(i,j,k)=min(cnis(i,j,k),qi(i,j,k))

                else

                  cnis(i,j,k)=0.e0

                end if

! -----

! Calculate the conversion rate from the snow to the graupel.

                if(qs(i,j,k).gt.thresq) then

                  if(vdvs(i,j,k).gt.0.e0                                &
     &              .and.vdvs(i,j,k).lt.clcs(i,j,k)) then

                    cnsg(i,j,k)=min(ccnsg*clcs(i,j,k),qs(i,j,k))

                    cnsgn(i,j,k)=min(ccnsgn*ecs(i,j,k)                  &
     &                *qc(i,j,k)*ncs(i,j,k)*sqrt(rbr(i,j,k))            &
     &                *exp(busm1*log(diaqs(i,j,k))),ncs(i,j,k))

                  else

                    cnsg(i,j,k)=0.e0
                    cnsgn(i,j,k)=0.e0

                  end if

                else

                  cnsg(i,j,k)=0.e0
                  cnsgn(i,j,k)=0.e0

                end if

! -----

!! -----

! Fill in the array with 0 in the case the air temperature is greater
! than melting points.

              else

                cnis(i,j,k)=0.e0

                cnsg(i,j,k)=0.e0
                cnsgn(i,j,k)=0.e0

              end if

! -----

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

      end subroutine s_convers

!-----7--------------------------------------------------------------7--

      end module m_convers
