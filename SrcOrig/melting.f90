!***********************************************************************
      module m_melting
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2000/10/18, 2001/06/29, 2001/10/18,
!                   2001/11/20, 2001/12/11, 2002/01/15, 2002/04/02,
!                   2003/03/28, 2003/04/30, 2003/05/19, 2003/12/12,
!                   2004/03/22, 2004/04/01, 2004/06/10, 2004/09/01,
!                   2004/09/25, 2004/10/12, 2004/12/17, 2006/04/03,
!                   2007/10/19, 2008/01/11, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2009/11/13, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the melting rate from the cloud ice to the cloud water
!     and from the snow and graupel to the rain water.

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

      public :: melting, s_melting

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface melting

        module procedure s_melting

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic max

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_melting(dtb,thresq,ni,nj,nk,rbr,rbv,qi,qs,qg,tcel,   &
     &                     qvsst0,lv,lf,kp,dv,vnts,vntg,clcs,clcg,      &
     &                     clrs,clrg,mlic,mlsr,mlgr)
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

      real, intent(in) :: qi(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixing ratio

      real, intent(in) :: qs(0:ni+1,0:nj+1,1:nk)
                       ! Snow mixing ratio

      real, intent(in) :: qg(0:ni+1,0:nj+1,1:nk)
                       ! Graupel mixing ratio

      real, intent(in) :: tcel(0:ni+1,0:nj+1,1:nk)
                       ! Ambient air temperature

      real, intent(in) :: qvsst0(0:ni+1,0:nj+1,1:nk)
                       ! Super saturation mixing ratio at melting point

      real, intent(in) :: lv(0:ni+1,0:nj+1,1:nk)
                       ! Latent heat of evapolation

      real, intent(in) :: lf(0:ni+1,0:nj+1,1:nk)
                       ! Latent heat of fusion

      real, intent(in) :: kp(0:ni+1,0:nj+1,1:nk)
                       ! Thermal conductivity of air

      real, intent(in) :: dv(0:ni+1,0:nj+1,1:nk)
                       ! Molecular diffusivity of water

      real, intent(in) :: vnts(0:ni+1,0:nj+1,1:nk)
                       ! Ventilation factor for snow

      real, intent(in) :: vntg(0:ni+1,0:nj+1,1:nk)
                       ! Ventilation factor for graupel

      real, intent(in) :: clcs(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate between cloud water and snow

      real, intent(in) :: clcg(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate between cloud water and graupel

      real, intent(in) :: clrs(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate from rain water to snow

      real, intent(in) :: clrg(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate between rain water and graupel

! Output variables

      real, intent(out) :: mlic(0:ni+1,0:nj+1,1:nk)
                       ! Melting rate from cloud ice to cloud water

      real, intent(out) :: mlsr(0:ni+1,0:nj+1,1:nk)
                       ! Melting rate from snow to rain water

      real, intent(out) :: mlgr(0:ni+1,0:nj+1,1:nk)
                       ! Melting rate from graupel to rain water

! Internal shared variable

      real cc2dt       ! 2.0 x cc x dtb

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real cmlxr1      ! Coefficient of melting rate
      real cmlxr2      ! Coefficient of melting rate

!-----7--------------------------------------------------------------7--

! Set the common used variable.

      cc2dt=2.e0*cc*dtb

! -----

!!! Calculate the melting rate from the cloud ice to the cloud water and
!!! from the snow and graupel to the rain water.

!$omp parallel default(shared) private(k)

!! In the case nk = 1.

      if(nk.eq.1) then

!$omp do schedule(runtime) private(i,j,cmlxr1,cmlxr2)

        do j=1,nj-1
        do i=1,ni-1

! Set the common used variables.

          cmlxr1=cc2dt*rbv(i,j,1)*(kp(i,j,1)*tcel(i,j,1)                &
     &      +lv(i,j,1)*dv(i,j,1)*rbr(i,j,1)*qvsst0(i,j,1))

          cmlxr2=cw*tcel(i,j,1)

! -----

! Calculate the melting rate from the cloud ice to the cloud water.

          if(qi(i,j,1).gt.thresq) then

            if(tcel(i,j,1).ge.t0cel) then

              mlic(i,j,1)=qi(i,j,1)

            else

              mlic(i,j,1)=0.e0

            end if

          else

            mlic(i,j,1)=0.e0

          end if

! -----

! Calculate the melting rate from the snow to the rain water.

          if(qs(i,j,1).gt.thresq) then

            if(tcel(i,j,1).ge.t0cel) then

              mlsr(i,j,1)=max((cmlxr1*vnts(i,j,1)                       &
     &          +cmlxr2*(clcs(i,j,1)+clrs(i,j,1)))/lf(i,j,1),0.e0)

            else

              mlsr(i,j,1)=0.e0

            end if

          else

            mlsr(i,j,1)=0.e0

          end if

! -----

! Calculate the melting rate from the graupel to the rain water.

          if(qg(i,j,1).gt.thresq) then

            if(tcel(i,j,1).ge.t0cel) then

              mlgr(i,j,1)=max((cmlxr1*vntg(i,j,1)                       &
     &          +cmlxr2*(clcg(i,j,1)+clrg(i,j,1)))/lf(i,j,1),0.e0)

            else

              mlgr(i,j,1)=0.e0

            end if

          else

            mlgr(i,j,1)=0.e0

          end if

! -----

        end do
        end do

!$omp end do

!! -----

!! In the case nk > 1.

      else

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j,cmlxr1,cmlxr2)

          do j=1,nj-1
          do i=1,ni-1

! Set the common used variables.

            cmlxr1=cc2dt*rbv(i,j,k)*(kp(i,j,k)*tcel(i,j,k)              &
     &        +lv(i,j,k)*dv(i,j,k)*rbr(i,j,k)*qvsst0(i,j,k))

            cmlxr2=cw*tcel(i,j,k)

! -----

! Calculate the melting rate from the cloud ice to the cloud water.

            if(qi(i,j,k).gt.thresq) then

              if(tcel(i,j,k).ge.t0cel) then

                mlic(i,j,k)=qi(i,j,k)

              else

                mlic(i,j,k)=0.e0

              end if

            else

              mlic(i,j,k)=0.e0

            end if

! -----

! Calculate the melting rate from the snow to the rain water.

            if(qs(i,j,k).gt.thresq) then

              if(tcel(i,j,k).ge.t0cel) then

                mlsr(i,j,k)=max((cmlxr1*vnts(i,j,k)                     &
     &            +cmlxr2*(clcs(i,j,k)+clrs(i,j,k)))/lf(i,j,k),0.e0)

              else

                mlsr(i,j,k)=0.e0

              end if

            else

              mlsr(i,j,k)=0.e0

            end if

! -----

! Calculate the melting rate from the graupel to the rain water.

            if(qg(i,j,k).gt.thresq) then

              if(tcel(i,j,k).ge.t0cel) then

                mlgr(i,j,k)=max((cmlxr1*vntg(i,j,k)                     &
     &            +cmlxr2*(clcg(i,j,k)+clrg(i,j,k)))/lf(i,j,k),0.e0)

              else

                mlgr(i,j,k)=0.e0

              end if

            else

              mlgr(i,j,k)=0.e0

            end if

! -----

          end do
          end do

!$omp end do

        end do

      end if

!! -----

!$omp end parallel

!!! -----

      end subroutine s_melting

!-----7--------------------------------------------------------------7--

      end module m_melting
