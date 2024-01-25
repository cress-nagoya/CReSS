!***********************************************************************
      module m_setsfc
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/10/15
!     Modification: 2001/10/18, 2001/11/14, 2002/01/15, 2002/04/02,
!                   2002/07/03, 2002/12/02, 2002/12/27, 2003/03/28,
!                   2003/04/30, 2003/05/19, 2003/07/15, 2003/10/31,
!                   2003/12/12, 2004/03/05, 2004/04/15, 2004/05/07,
!                   2004/08/01, 2004/09/01, 2004/09/10, 2005/01/31,
!                   2007/01/20, 2007/07/30, 2007/10/19, 2008/05/02,
!                   2008/07/01, 2008/08/25, 2009/02/27, 2009/11/13,
!                   2011/01/05, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the magnitude of velocity, virtual potential temperature
!     and water vapor mixing ratio on the surface.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comphy
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: setsfc, s_setsfc

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface setsfc

        module procedure s_setsfc

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic exp
      intrinsic log
      intrinsic max
      intrinsic min
      intrinsic sqrt

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_setsfc(fplevpbl,fptubopt,fpcphopt,fmois,             &
     &                    ni,nj,nk,nund,pbr,ptbr,u,v,w,pp,ptp,qv,       &
     &                    land,beta,kai,tund,fall,p,t,ptv,qvsfc,tice,va,&
     &                    pt,ps,qvsts,qvsice)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: fplevpbl
                       ! Formal parameter of unique index of levpbl

      integer, intent(in) :: fptubopt
                       ! Formal parameter of unique index of tubopt

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

      real, intent(in) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! w components of velocity

      real, intent(in) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation

      real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

      real, intent(in) :: beta(0:ni+1,0:nj+1)
                       ! Evapotranspiration efficiency

      real, intent(in) :: kai(0:ni+1,0:nj+1)
                       ! Sea ice distribution

      real, intent(in) :: tund(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature

      real, intent(in) :: fall(0:ni+1,0:nj+1)
                       ! Precipitation flag

! Output variables

      real, intent(out) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure

      real, intent(out) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

      real, intent(out) :: ptv(0:ni+1,0:nj+1,1:nk)
                       ! Virtual potential temperature

      real, intent(out) :: qvsfc(0:ni+1,0:nj+1)
                       ! Water vapor mixing ratio on surface

      real, intent(out) :: tice(0:ni+1,0:nj+1)
                       ! Mixed ice surface temperature

      real, intent(out) :: va(0:ni+1,0:nj+1)
                       ! Magnitude of velocity at lowest plane

! Internal shared variables

      integer levpbl   ! Number of planetary boundary layer

      integer tubopt   ! Option for turbulent mixing
      integer cphopt   ! Option for cloud micro physics

      real rddvcp      ! rd / cp

      real p0iv        ! 1.0 / p0

      real, intent(inout) :: pt(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature

      real, intent(inout) :: ps(0:ni+1,0:nj+1)
                       ! Surface pressure

      real, intent(inout) :: qvsts(0:ni+1,0:nj+1)
                       ! Saturation mixing ratio
                       ! at surface temperature

      real, intent(inout) :: qvsice(0:ni+1,0:nj+1)
                       ! Saturation mixing ratio
                       ! at ice surface temperature

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real ests        ! Saturation vapor pressure
                       ! at surface temperature

      real esice       ! Saturation vapor pressure
                       ! at ice surface temperature

      real u8s         ! x components of velocity at scalar points
      real v8s         ! y components of velocity at scalar points
      real w8s         ! z components of velocity at scalar points

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fplevpbl,levpbl)
      call getiname(fptubopt,tubopt)
      call getiname(fpcphopt,cphopt)

! -----

! Set the common used variables.

      rddvcp=rd/cp

      p0iv=1.e0/p0

! -----

!! Calculate the surface parameters.

!$omp parallel default(shared) private(k)

! Get the pressure, potential temperature and air temperature.

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          p(i,j,k)=pbr(i,j,k)+pp(i,j,k)
          pt(i,j,k)=ptbr(i,j,k)+ptp(i,j,k)

          t(i,j,k)=pt(i,j,k)*exp(rddvcp*log(p0iv*p(i,j,k)))

        end do
        end do

!$omp end do

      end do

! -----

! Get the surface pressure and ice surface temperature.

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=1,ni-1

        ps(i,j)=.5e0*(p(i,j,1)+p(i,j,2))

        if(land(i,j).eq.1) then

          tice(i,j)=min(.5e0*(t(i,j,1)+t(i,j,2)),t0)

        else

          tice(i,j)=lim35n

        end if

      end do
      end do

!$omp end do

! -----

! Calculate the water vapor mixing ratio on the surface.

      if(fmois(1:3).eq.'dry') then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          qvsfc(i,j)=0.e0
        end do
        end do

!$omp end do

      else if(fmois(1:5).eq.'moist') then

!$omp do schedule(runtime) private(i,j,ests,esice)

        do j=1,nj-1
        do i=1,ni-1

          if(land(i,j).lt.0) then

            ests=es0                                                    &
     &        *exp(17.269e0*(tund(i,j,1)-t0)/(tund(i,j,1)-35.86e0))

            qvsts(i,j)=epsva*ests/(ps(i,j)-ests)

          else if(land(i,j).eq.1) then

            if(tice(i,j).gt.tlow) then

              ests=es0                                                  &
     &          *exp(17.269e0*(tund(i,j,1)-t0)/(tund(i,j,1)-35.86e0))

              esice=es0                                                 &
     &          *exp(17.269e0*(tice(i,j)-t0)/(tice(i,j)-35.86e0))

              qvsts(i,j)=epsva*ests/(ps(i,j)-ests)

              qvsice(i,j)=epsva*esice/(ps(i,j)-esice)

            else

              ests=es0                                                  &
     &          *exp(17.269e0*(tund(i,j,1)-t0)/(tund(i,j,1)-35.86e0))

              esice=es0                                                 &
     &          *exp(21.875e0*(tice(i,j)-t0)/(tice(i,j)-7.66e0))

              qvsts(i,j)=epsva*ests/(ps(i,j)-ests)

              qvsice(i,j)=epsva*esice/(ps(i,j)-esice)

            end if

          else

            if(tund(i,j,1).gt.tlow) then

              ests=es0                                                  &
     &          *exp(17.269e0*(tund(i,j,1)-t0)/(tund(i,j,1)-35.86e0))

              qvsts(i,j)=epsva*ests/(ps(i,j)-ests)

            else

              ests=es0                                                  &
     &          *exp(21.875e0*(tund(i,j,1)-t0)/(tund(i,j,1)-7.66e0))

              qvsts(i,j)=epsva*ests/(ps(i,j)-ests)

            end if

          end if

        end do
        end do

!$omp end do

        if(abs(cphopt).eq.0) then

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1

            if(land(i,j).eq.1) then

              qvsfc(i,j)=qv(i,j,2)                                      &
     &          +icbeta*kai(i,j)*(qvsice(i,j)-qv(i,j,2))                &
     &          +beta(i,j)*(1.e0-kai(i,j))*(qvsts(i,j)-qv(i,j,2))

            else

              qvsfc(i,j)=beta(i,j)*(qvsts(i,j)-qv(i,j,2))+qv(i,j,2)

            end if

          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1

            if(land(i,j).eq.1) then

              if(fall(i,j).gt.0.e0) then

                qvsfc(i,j)                                              &
     &            =kai(i,j)*qvsice(i,j)+(1.e0-kai(i,j))*qvsts(i,j)

              else

                qvsfc(i,j)=qv(i,j,2)                                    &
     &            +icbeta*kai(i,j)*(qvsice(i,j)-qv(i,j,2))              &
     &            +beta(i,j)*(1.e0-kai(i,j))*(qvsts(i,j)-qv(i,j,2))

              end if

            else

              if(fall(i,j).gt.0.e0) then

!ORIG           qvsfc(i,j)=qvsts(i,j)
                qvsfc(i,j)=beta(i,j)*(qvsts(i,j)-qv(i,j,2))+qv(i,j,2)

              else

                qvsfc(i,j)=beta(i,j)*(qvsts(i,j)-qv(i,j,2))+qv(i,j,2)

              end if

            end if

          end do
          end do

!$omp end do

        end if

      end if

! -----

! Calculate the virtual potential temperature.

      if(fmois(1:3).eq.'dry') then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1

          if(land(i,j).eq.1) then

            ptv(i,j,1)=exp(rddvcp*log(p0/ps(i,j)))                      &
     &        *(kai(i,j)*tice(i,j)+(1.e0-kai(i,j))*tund(i,j,1))

          else

            ptv(i,j,1)=exp(rddvcp*log(p0/ps(i,j)))*tund(i,j,1)

          end if

        end do
        end do

!$omp end do

        do k=2,levpbl+1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            ptv(i,j,k)=pt(i,j,k)
          end do
          end do

!$omp end do

        end do

      else if(fmois(1:5).eq.'moist') then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1

          if(land(i,j).eq.1) then

            ptv(i,j,1)=exp(rddvcp*log(p0/ps(i,j)))                      &
     &        *(1.e0+epsav*qvsfc(i,j))/(1.e0+qvsfc(i,j))                &
     &        *(kai(i,j)*tice(i,j)+(1.e0-kai(i,j))*tund(i,j,1))

          else

            ptv(i,j,1)=exp(rddvcp*log(p0/ps(i,j)))                      &
     &        *tund(i,j,1)*(1.e0+epsav*qvsfc(i,j))/(1.e0+qvsfc(i,j))

          end if

        end do
        end do

!$omp end do

        do k=2,levpbl+1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            ptv(i,j,k)=pt(i,j,k)*(1.e0+epsav*qv(i,j,k))/(1.e0+qv(i,j,k))
          end do
          end do

!$omp end do

        end do

      end if

! -----

! Calculate the magnitude of velocity.

      if(tubopt.eq.0) then

!$omp do schedule(runtime) private(i,j,u8s,v8s)

        do j=1,nj-1
        do i=1,ni-1
          u8s=u(i,j,2)+u(i+1,j,2)
          v8s=v(i,j,2)+v(i,j+1,2)

          va(i,j)=max(.5e0*sqrt(u8s*u8s+v8s*v8s),vamin)

        end do
        end do

!$omp end do

      else

!$omp do schedule(runtime) private(i,j,u8s,v8s,w8s)

        do j=1,nj-1
        do i=1,ni-1
          u8s=u(i,j,2)+u(i+1,j,2)
          v8s=v(i,j,2)+v(i,j+1,2)
          w8s=w(i,j,2)+w(i,j,3)

          va(i,j)=max(.5e0*sqrt(u8s*u8s+v8s*v8s+w8s*w8s),vamin)

        end do
        end do

!$omp end do

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_setsfc

!-----7--------------------------------------------------------------7--

      end module m_setsfc
