!***********************************************************************
      module m_depsit
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2000/10/18, 2000/11/17, 2001/06/29,
!                   2001/10/18, 2001/11/20, 2001/12/11, 2002/01/15,
!                   2002/04/01, 2002/04/02, 2003/03/28, 2003/04/30,
!                   2003/05/19, 2003/09/01, 2003/10/31, 2003/11/05,
!                   2003/12/12, 2004/04/01, 2004/04/15, 2004/06/10,
!                   2004/08/01, 2004/09/01, 2004/09/25, 2004/10/12,
!                   2004/12/17, 2005/01/31, 2005/04/04, 2005/10/05,
!                   2006/02/13, 2006/03/06, 2006/04/03, 2006/09/30,
!                   2007/10/19, 2008/01/11, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2009/11/13, 2011/04/08, 2011/06/01,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the evaporation rate from the rain water to the water
!     vapor and the deposition rate from the water vapor to the ice
!     hydrometeor.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comphy
      use m_comtable

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: depsit, s_depsit

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface depsit

        module procedure s_depsit

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic int
      intrinsic log
      intrinsic max

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_depsit(dtb,thresq,ni,nj,nk,rbr,rbv,qv,qr,qi,qs,qg,   &
     &                    nci,t,tcel,qvsst0,qvsw,qvsi,lv,ls,lf,kp,dv,   &
     &                    mi,vntr,vnts,vntg,clcs,clcg,mlsr,mlgr,        &
     &                    vdvr,vdvi,vdvs,vdvg)
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

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rbv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of base state density

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

      real, intent(in) :: qr(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixing ratio

      real, intent(in) :: qi(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixing ratio

      real, intent(in) :: qs(0:ni+1,0:nj+1,1:nk)
                       ! Snow mixing ratio

      real, intent(in) :: qg(0:ni+1,0:nj+1,1:nk)
                       ! Graupel mixing ratio

      real, intent(in) :: nci(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud ice

      real, intent(in) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

      real, intent(in) :: tcel(0:ni+1,0:nj+1,1:nk)
                       ! Ambient air temperature

      real, intent(in) :: qvsst0(0:ni+1,0:nj+1,1:nk)
                       ! Super saturation mixing ratio at melting point

      real, intent(in) :: qvsw(0:ni+1,0:nj+1,1:nk)
                       ! Saturation mixing ratio for water

      real, intent(in) :: qvsi(0:ni+1,0:nj+1,1:nk)
                       ! Saturation mixing ratio for ice

      real, intent(in) :: lv(0:ni+1,0:nj+1,1:nk)
                       ! Latent heat of evapolation

      real, intent(in) :: ls(0:ni+1,0:nj+1,1:nk)
                       ! Latent heat of sublimation

      real, intent(in) :: lf(0:ni+1,0:nj+1,1:nk)
                       ! Latent heat of fusion

      real, intent(in) :: kp(0:ni+1,0:nj+1,1:nk)
                       ! Thermal conductivity of air

      real, intent(in) :: dv(0:ni+1,0:nj+1,1:nk)
                       ! Molecular diffusivity of water

      real, intent(in) :: mi(0:ni+1,0:nj+1,1:nk)
                       ! Mean mass of cloud ice

      real, intent(in) :: vntr(0:ni+1,0:nj+1,1:nk)
                       ! Ventilation factor for rain water

      real, intent(in) :: vnts(0:ni+1,0:nj+1,1:nk)
                       ! Ventilation factor for snow

      real, intent(in) :: vntg(0:ni+1,0:nj+1,1:nk)
                       ! Ventilation factor for graupel

      real, intent(in) :: clcs(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate between cloud water and snow

      real, intent(in) :: clcg(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate between cloud water and graupel

      real, intent(in) :: mlsr(0:ni+1,0:nj+1,1:nk)
                       ! Melting rate from snow to rain water

      real, intent(in) :: mlgr(0:ni+1,0:nj+1,1:nk)
                       ! Melting rate from graupel to rain water

! Output variables

      real, intent(out) :: vdvr(0:ni+1,0:nj+1,1:nk)
                       ! Evaporation rate from rain water to water vapor

      real, intent(out) :: vdvi(0:ni+1,0:nj+1,1:nk)
                       ! Deposiotn rate from water vapor to cloud ice

      real, intent(out) :: vdvs(0:ni+1,0:nj+1,1:nk)
                       ! Deposiotn rate from water vapor to snow

      real, intent(out) :: vdvg(0:ni+1,0:nj+1,1:nk)
                       ! Deposiotn rate from water vapor to graupel

! Internal shared variables

      real ccdtb4      ! 2.0 x cc x dtb

      real t27311      ! t0 - 0.05

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer itsc     ! int(super cooled air temperature)

      real qvssw       ! Super saturation mixing ratio for water
      real qvssi       ! Super saturation mixing ratio for ice

      real gi          ! Thermodynamical function Gi(T,P)

      real cvdvx1      ! Coefficient of deposition rate
      real cvdvx2      ! Coefficient of deposition rate
      real cvdvx3      ! Coefficient of deposition rate
      real cvdvx4      ! Coefficient of deposition rate

      real sink        ! Sink amount in deposition

      real a           ! Temporary variable
      real b           ! Temporary variable
      real c           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Set the common used variables.

      ccdtb4=2.e0*cc*dtb

      t27311=t0-.01e0

! -----

!!!! Calculate the evaporation and deposition rate.

!$omp parallel default(shared) private(k)

!!! In the case nk = 1.

      if(nk.eq.1) then

!$omp do schedule(runtime) private(i,j,itsc)                            &
!$omp&   private(qvssw,qvssi,gi,cvdvx1,cvdvx2,cvdvx3,cvdvx4,sink,a,b,c)

        do j=1,nj-1
        do i=1,ni-1

!! Perform calculating in the case the air temperature is lower than
!! lowest super cooled point.

          if(t(i,j,1).gt.tlow) then

! Set the common used variables.

            qvssw=qv(i,j,1)-qvsw(i,j,1)
            qvssi=qv(i,j,1)-qvsi(i,j,1)

            a=1.e0/(rv*kp(i,j,1)*t(i,j,1)*t(i,j,1))
            b=dv(i,j,1)*rbr(i,j,1)
            c=ccdtb4*rbv(i,j,1)

            gi=1.e0/(ls(i,j,1)*ls(i,j,1)*a+1.e0/(b*qvsi(i,j,1)))

            cvdvx1=c*(qv(i,j,1)/qvsw(i,j,1)-1.e0)                       &
     &        /(lv(i,j,1)*lv(i,j,1)*a+1.e0/(b*qvsw(i,j,1)))

            cvdvx2=ccdtb4*dv(i,j,1)*qvsst0(i,j,1)
            cvdvx3=c*(qv(i,j,1)/qvsi(i,j,1)-1.e0)
            cvdvx4=a*ls(i,j,1)*lf(i,j,1)

! -----

! Calculate the evaporation rate from the rain water to the water vapor.

            if(qr(i,j,1).gt.thresq) then

              if(qv(i,j,1).lt.qvsw(i,j,1)) then

                vdvr(i,j,1)=max(cvdvx1*vntr(i,j,1),qvssw)

              else

                vdvr(i,j,1)=0.e0

              end if

            else

              vdvr(i,j,1)=0.e0

            end if

! -----

! Calculate the deposition rate from the water vapor to the cloud ice.

            if(qi(i,j,1).gt.thresq) then

              if(t(i,j,1).lt.t27311) then

                itsc=int(-tcel(i,j,1))

                vdvi(i,j,1)=ckoe(itsc)*(qv(i,j,1)-qvsi(i,j,1))          &
     &            /(qvsw(i,j,1)-qvsi(i,j,1))*nci(i,j,1)                 &
     &            *exp(pkoe(itsc)*log(mi(i,j,1)))*dtb

              else

                vdvi(i,j,1)=0.e0

              end if

            else

              vdvi(i,j,1)=0.e0

            end if

! -----

! Calculate the deposition rate from the water vapor to the snow.

            if(qs(i,j,1).gt.thresq) then

              if(t(i,j,1).gt.t0) then

                if(mlsr(i,j,1).gt.0.e0) then

                  vdvs(i,j,1)=cvdvx2*vnts(i,j,1)

                else

                  vdvs(i,j,1)=cvdvx1*vnts(i,j,1)

                end if

              else

                vdvs(i,j,1)=gi*(cvdvx3*vnts(i,j,1)-cvdvx4*clcs(i,j,1))

              end if

            else

              vdvs(i,j,1)=0.e0

            end if

! -----

! Calculate the deposition rate from the water vapor to the graupel.

            if(qg(i,j,1).gt.thresq) then

              if(t(i,j,1).gt.t0) then

                if(mlgr(i,j,1).gt.0.e0) then

                  vdvg(i,j,1)=cvdvx2*vntg(i,j,1)

                else

                  vdvg(i,j,1)=cvdvx1*vntg(i,j,1)

                end if

              else

                vdvg(i,j,1)=gi*(cvdvx3*vntg(i,j,1)-cvdvx4*clcg(i,j,1))

              end if

            else

              vdvg(i,j,1)=0.e0

            end if

! -----

! Adjust the deposition rate from the water vapor to the ice
! hydrometeor.

            sink=vdvi(i,j,1)+vdvs(i,j,1)+vdvg(i,j,1)

            if((qvssi.lt.sink.and.qvssi.gt.0.e0)                        &
     &        .or.(qvssi.gt.sink.and.qvssi.lt.0.e0)) then

              a=qvssi/sink

              vdvi(i,j,1)=vdvi(i,j,1)*a
              vdvs(i,j,1)=vdvs(i,j,1)*a
              vdvg(i,j,1)=vdvg(i,j,1)*a

            end if

! -----

!! -----

! Fill in the array with 0 in the case the air temperature is lower than
! lowest super cooled point.

          else

            vdvr(i,j,1)=0.e0
            vdvi(i,j,1)=0.e0
            vdvs(i,j,1)=0.e0
            vdvg(i,j,1)=0.e0

          end if

! -----

        end do
        end do

!$omp end do

!!! -----

!!! In the case nk > 1.

      else

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j,itsc)                            &
!$omp&   private(qvssw,qvssi,gi,cvdvx1,cvdvx2,cvdvx3,cvdvx4,sink,a,b,c)

          do j=1,nj-1
          do i=1,ni-1

!! Perform calculating in the case the air temperature is lower than
!! lowest super cooled point.

            if(t(i,j,k).gt.tlow) then

! Set the common used variables.

              qvssw=qv(i,j,k)-qvsw(i,j,k)
              qvssi=qv(i,j,k)-qvsi(i,j,k)

              a=1.e0/(rv*kp(i,j,k)*t(i,j,k)*t(i,j,k))
              b=dv(i,j,k)*rbr(i,j,k)
              c=ccdtb4*rbv(i,j,k)

              gi=1.e0/(ls(i,j,k)*ls(i,j,k)*a+1.e0/(b*qvsi(i,j,k)))

              cvdvx1=c*(qv(i,j,k)/qvsw(i,j,k)-1.e0)                     &
     &          /(lv(i,j,k)*lv(i,j,k)*a+1.e0/(b*qvsw(i,j,k)))

              cvdvx2=ccdtb4*dv(i,j,k)*qvsst0(i,j,k)
              cvdvx3=c*(qv(i,j,k)/qvsi(i,j,k)-1.e0)
              cvdvx4=a*ls(i,j,k)*lf(i,j,k)

! -----

! Calculate the evaporation rate from the rain water to the water vapor.

              if(qr(i,j,k).gt.thresq) then

                if(qv(i,j,k).lt.qvsw(i,j,k)) then

                  vdvr(i,j,k)=max(cvdvx1*vntr(i,j,k),qvssw)

                else

                  vdvr(i,j,k)=0.e0

                end if

              else

                vdvr(i,j,k)=0.e0

              end if

! -----

! Calculate the deposition rate from the water vapor to the cloud ice.

              if(qi(i,j,k).gt.thresq) then

                if(t(i,j,k).lt.t27311) then

                  itsc=int(-tcel(i,j,k))

                  vdvi(i,j,k)=ckoe(itsc)*(qv(i,j,k)-qvsi(i,j,k))        &
     &              /(qvsw(i,j,k)-qvsi(i,j,k))*nci(i,j,k)               &
     &              *exp(pkoe(itsc)*log(mi(i,j,k)))*dtb

                else

                  vdvi(i,j,k)=0.e0

                end if

              else

                vdvi(i,j,k)=0.e0

              end if

! -----

! Calculate the deposition rate from the water vapor to the snow.

              if(qs(i,j,k).gt.thresq) then

                if(t(i,j,k).gt.t0) then

                  if(mlsr(i,j,k).gt.0.e0) then

                    vdvs(i,j,k)=cvdvx2*vnts(i,j,k)

                  else

                    vdvs(i,j,k)=cvdvx1*vnts(i,j,k)

                  end if

                else

                  vdvs(i,j,k)=gi*(cvdvx3*vnts(i,j,k)-cvdvx4*clcs(i,j,k))

                end if

              else

                vdvs(i,j,k)=0.e0

              end if

! -----

! Calculate the deposition rate from the water vapor to the graupel.

              if(qg(i,j,k).gt.thresq) then

                if(t(i,j,k).gt.t0) then

                  if(mlgr(i,j,k).gt.0.e0) then

                    vdvg(i,j,k)=cvdvx2*vntg(i,j,k)

                  else

                    vdvg(i,j,k)=cvdvx1*vntg(i,j,k)

                  end if

                else

                  vdvg(i,j,k)=gi*(cvdvx3*vntg(i,j,k)-cvdvx4*clcg(i,j,k))

                end if

              else

                vdvg(i,j,k)=0.e0

              end if

! -----

! Adjust the deposition rate from the water vapor to the ice
! hydrometeor.

              sink=vdvi(i,j,k)+vdvs(i,j,k)+vdvg(i,j,k)

              if((qvssi.lt.sink.and.qvssi.gt.0.e0)                      &
     &          .or.(qvssi.gt.sink.and.qvssi.lt.0.e0)) then

                a=qvssi/sink

                vdvi(i,j,k)=vdvi(i,j,k)*a
                vdvs(i,j,k)=vdvs(i,j,k)*a
                vdvg(i,j,k)=vdvg(i,j,k)*a

              end if

! -----

!! -----

! Fill in the array with 0 in the case the air temperature is lower than
! lowest super cooled point.

            else

              vdvr(i,j,k)=0.e0
              vdvi(i,j,k)=0.e0
              vdvs(i,j,k)=0.e0
              vdvg(i,j,k)=0.e0

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

      end subroutine s_depsit

!-----7--------------------------------------------------------------7--

      end module m_depsit
