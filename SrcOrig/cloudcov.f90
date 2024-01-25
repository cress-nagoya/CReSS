!***********************************************************************
      module m_cloudcov
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/08/31
!     Modification: 2001/10/18, 2001/12/07, 2002/04/02, 2002/12/02,
!                   2003/04/30, 2003/05/19, 2003/10/31, 2003/11/28,
!                   2003/12/12, 2004/05/31, 2004/08/01, 2006/01/10,
!                   2007/01/20, 2007/01/31, 2007/06/27, 2007/07/30,
!                   2007/10/19, 2008/03/12, 2008/05/02, 2008/07/01,
!                   2008/08/25, 2009/02/27, 2009/08/20, 2009/11/05,
!                   2011/04/08, 2011/11/10, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the cloud cover.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comphy
      use m_comtable
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: cloudcov, s_cloudcov

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface cloudcov

        module procedure s_cloudcov

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic aint
      intrinsic exp
      intrinsic int
      intrinsic max
      intrinsic min

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_cloudcov(fpcphopt,fpdz,fmois,fproc,ni,nj,nk,         &
     &                      zph,rst,p,t,qv,qall,cdl,cdm,cdh,zph8s,      &
     &                      rh24,rh32,rh48,rh72,qsum,qsuml,qsumm,qsumh)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      character(len=3), intent(in) :: fproc
                       ! Control flag of processing type

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fpdz
                       ! Formal parameter of unique index of dz

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure

      real, intent(in) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

      real, intent(in) :: qall(0:ni+1,0:nj+1,1:nk)
                       ! Total water and ice mixing ratio

! Output variables

      real, intent(out) :: cdl(0:ni+1,0:nj+1)
                       ! Cloud cover in lower layer

      real, intent(out) :: cdm(0:ni+1,0:nj+1)
                       ! Cloud cover in middle layer

      real, intent(out) :: cdh(0:ni+1,0:nj+1)
                       ! Cloud cover in upper layer

! Internal shared variables

      integer cphopt   ! Option for cloud micro physics

      real dz          ! Grid distance in z direction

      real es0iv2      ! Inverse of es0 x 100

      real, intent(inout) :: zph8s(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates at scalar points

      real, intent(inout) :: rh24(0:ni+1,0:nj+1)
                       ! Relative humidity at 2400 [m]

      real, intent(inout) :: rh32(0:ni+1,0:nj+1)
                       ! Relative humidity at 3200 [m]

      real, intent(inout) :: rh48(0:ni+1,0:nj+1)
                       ! Relative humidity at 4800 [m]

      real, intent(inout) :: rh72(0:ni+1,0:nj+1)
                       ! Relative humidity at 7200 [m]

      real, intent(inout) :: qsum(0:ni+1,0:nj+1,1:nk)
                       ! Total water cover

      real, intent(inout) :: qsuml(0:ni+1,0:nj+1)
                       ! Total water cover in lower layer

      real, intent(inout) :: qsumm(0:ni+1,0:nj+1)
                       ! Total water cover in middle layer

      real, intent(inout) :: qsumh(0:ni+1,0:nj+1)
                       ! Total water cover in upper layer

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer irh      ! int(rhxx)

      real rhsfc       ! Relative humidity at surface

      real rha         ! Relative humidity above reference layer
      real rhb         ! Relative humidity below reference layer

      real dk          ! Weighting distance for interpolating

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpcphopt,cphopt)
      call getrname(fpdz,dz)

! -----

! Set the common used variable.

      es0iv2=100.e0/es0

! -----

!!! Calculate the cloud cover.

!$omp parallel default(shared) private(k)

! Set the no cloud cover in the case of dry air.

      if(fmois(1:3).eq.'dry') then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          cdl(i,j)=0.e0
          cdm(i,j)=0.e0
          cdh(i,j)=0.e0
        end do
        end do

!$omp end do

! -----

!! Estimate the cloud cover in the case of moist air.

      else if(fmois(1:5).eq.'moist') then

! Get the z physical coordinates at scalar points.

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            zph8s(i,j,k)=.5e0*(zph(i,j,k)+zph(i,j,k+1))
          end do
          end do

!$omp end do

        end do

! -----

! Estimate the cloud cover from the relative humidity.

        if(fproc(1:2).eq.'rh'.or.abs(cphopt).eq.0) then

!$omp do schedule(runtime) private(i,j,rhsfc)

          do j=1,nj-1
          do i=1,ni-1

            if(t(i,j,2).gt.tlow) then

              rhsfc=es0iv2*qv(i,j,2)*p(i,j,2)/(epsva+qv(i,j,2))         &
     &          *exp(17.269e0*(t(i,j,2)-t0)/(35.86e0-t(i,j,2)))

            else

              rhsfc=es0iv2*qv(i,j,2)*p(i,j,2)/(epsva+qv(i,j,2))         &
     &          *exp(21.875e0*(t(i,j,2)-t0)/(7.66e0-t(i,j,2)))

            end if

            rh24(i,j)=rhsfc
            rh32(i,j)=rhsfc
            rh48(i,j)=rhsfc
            rh72(i,j)=rhsfc

          end do
          end do

!$omp end do

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j,rha,rhb,dk)

            do j=1,nj-1
            do i=1,ni-1

              if(zph8s(i,j,k-1).lt.2400.e0                              &
     &          .and.zph8s(i,j,k).ge.2400.e0) then

                if(t(i,j,k).gt.tlow) then

                  rha=es0iv2*qv(i,j,k)*p(i,j,k)/(epsva+qv(i,j,k))       &
     &              *exp(17.269e0*(t(i,j,k)-t0)/(35.86e0-t(i,j,k)))

                else

                  rha=es0iv2*qv(i,j,k)*p(i,j,k)/(epsva+qv(i,j,k))       &
     &              *exp(21.875e0*(t(i,j,k)-t0)/(7.66e0-t(i,j,k)))

                end if

                if(t(i,j,k-1).gt.tlow) then

                  rhb=es0iv2*qv(i,j,k-1)*p(i,j,k-1)/(epsva+qv(i,j,k-1)) &
     &              *exp(17.269e0*(t(i,j,k-1)-t0)/(35.86e0-t(i,j,k-1)))

                else

                  rhb=es0iv2*qv(i,j,k-1)*p(i,j,k-1)/(epsva+qv(i,j,k-1)) &
     &              *exp(21.875e0*(t(i,j,k-1)-t0)/(7.66e0-t(i,j,k-1)))

                end if

                dk=(zph8s(i,j,k)-2400.e0)/(zph8s(i,j,k)-zph8s(i,j,k-1))

                rh24(i,j)=rha*(1.e0-dk)+rhb*dk

              end if

              if(zph8s(i,j,k-1).lt.3200.e0                              &
     &          .and.zph8s(i,j,k).ge.3200.e0) then

                if(t(i,j,k).gt.tlow) then

                  rha=es0iv2*qv(i,j,k)*p(i,j,k)/(epsva+qv(i,j,k))       &
     &              *exp(17.269e0*(t(i,j,k)-t0)/(35.86e0-t(i,j,k)))

                else

                  rha=es0iv2*qv(i,j,k)*p(i,j,k)/(epsva+qv(i,j,k))       &
     &              *exp(21.875e0*(t(i,j,k)-t0)/(7.66e0-t(i,j,k)))

                end if

                if(t(i,j,k-1).gt.tlow) then

                  rhb=es0iv2*qv(i,j,k-1)*p(i,j,k-1)/(epsva+qv(i,j,k-1)) &
     &              *exp(17.269e0*(t(i,j,k-1)-t0)/(35.86e0-t(i,j,k-1)))

                else

                  rhb=es0iv2*qv(i,j,k-1)*p(i,j,k-1)/(epsva+qv(i,j,k-1)) &
     &              *exp(21.875e0*(t(i,j,k-1)-t0)/(7.66e0-t(i,j,k-1)))

                end if

                dk=(zph8s(i,j,k)-3200.e0)/(zph8s(i,j,k)-zph8s(i,j,k-1))

                rh32(i,j)=rha*(1.e0-dk)+rhb*dk

              end if

              if(zph8s(i,j,k-1).lt.4800.e0                              &
     &          .and.zph8s(i,j,k).ge.4800.e0) then

                if(t(i,j,k).gt.tlow) then

                  rha=es0iv2*qv(i,j,k)*p(i,j,k)/(epsva+qv(i,j,k))       &
     &              *exp(17.269e0*(t(i,j,k)-t0)/(35.86e0-t(i,j,k)))

                else

                  rha=es0iv2*qv(i,j,k)*p(i,j,k)/(epsva+qv(i,j,k))       &
     &              *exp(21.875e0*(t(i,j,k)-t0)/(7.66e0-t(i,j,k)))

                end if

                if(t(i,j,k-1).gt.tlow) then

                  rhb=es0iv2*qv(i,j,k-1)*p(i,j,k-1)/(epsva+qv(i,j,k-1)) &
     &              *exp(17.269e0*(t(i,j,k-1)-t0)/(35.86e0-t(i,j,k-1)))

                else

                  rhb=es0iv2*qv(i,j,k-1)*p(i,j,k-1)/(epsva+qv(i,j,k-1)) &
     &              *exp(21.875e0*(t(i,j,k-1)-t0)/(7.66e0-t(i,j,k-1)))

                end if

                dk=(zph8s(i,j,k)-4800.e0)/(zph8s(i,j,k)-zph8s(i,j,k-1))

                rh48(i,j)=rha*(1.e0-dk)+rhb*dk

              end if

              if(zph8s(i,j,k-1).lt.7200.e0                              &
     &          .and.zph8s(i,j,k).ge.7200.e0) then

                if(t(i,j,k).gt.tlow) then

                  rha=es0iv2*qv(i,j,k)*p(i,j,k)/(epsva+qv(i,j,k))       &
     &              *exp(17.269e0*(t(i,j,k)-t0)/(35.86e0-t(i,j,k)))

                else

                  rha=es0iv2*qv(i,j,k)*p(i,j,k)/(epsva+qv(i,j,k))       &
     &              *exp(21.875e0*(t(i,j,k)-t0)/(7.66e0-t(i,j,k)))

                end if

                if(t(i,j,k-1).gt.tlow) then

                  rhb=es0iv2*qv(i,j,k-1)*p(i,j,k-1)/(epsva+qv(i,j,k-1)) &
     &              *exp(17.269e0*(t(i,j,k-1)-t0)/(35.86e0-t(i,j,k-1)))

                else

                  rhb=es0iv2*qv(i,j,k-1)*p(i,j,k-1)/(epsva+qv(i,j,k-1)) &
     &              *exp(21.875e0*(t(i,j,k-1)-t0)/(7.66e0-t(i,j,k-1)))

                end if

                dk=(zph8s(i,j,k)-7200.e0)/(zph8s(i,j,k)-zph8s(i,j,k-1))

                rh72(i,j)=rha*(1.e0-dk)+rhb*dk

              end if

            end do
            end do

!$omp end do

          end do

!$omp do schedule(runtime) private(i,j,irh,dk)

          do j=1,nj-1
          do i=1,ni-1

            irh=min(int(rh24(i,j)),100)

            dk=rh24(i,j)-aint(rh24(i,j))

            cdl(i,j)=((1.e0-dk)*rcdl(irh)+dk*rcdl(irh+1))               &
     &        *min(max(3200.e0-zph(i,j,2),0.e0),3200.e0)/3200.e0

            irh=min(int(rh32(i,j)),100)

            dk=rh32(i,j)-aint(rh32(i,j))

            cdm(i,j)=((1.e0-dk)*rcdm(irh)+dk*rcdm(irh+1))               &
     &        *min(max(4800.e0-zph(i,j,2),0.e0),1600.e0)/1600.e0

            irh=min(int(rh48(i,j)),100)

            dk=rh48(i,j)-aint(rh48(i,j))

            cdh(i,j)=((1.e0-dk)*rcdm(irh)+dk*rcdm(irh+1))               &
     &        *min(max(7200.e0-zph(i,j,2),0.e0),2400.e0)/2400.e0

            cdm(i,j)=.5e0*(cdm(i,j)+cdh(i,j))

            irh=min(int(rh72(i,j)),100)

            dk=rh72(i,j)-aint(rh72(i,j))

            cdh(i,j)=(1.e0-dk)*rcdh(irh)+dk*rcdh(irh+1)

          end do
          end do

!$omp end do

! -----

! Estimate the cloud cover from the hydrometeor mixing ratio.

        else if(fproc(1:3).eq.'mix'.and.abs(cphopt).ne.0) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              qsum(i,j,k)=rst(i,j,k)*qall(i,j,k)*dz
            end do
            end do

!$omp end do

          end do

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            qsuml(i,j)=0.e0
            qsumm(i,j)=0.e0
            qsumh(i,j)=0.e0
          end do
          end do

!$omp end do

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1

              if(zph8s(i,j,k).lt.2800.e0) then

                qsuml(i,j)=qsuml(i,j)+qsum(i,j,k)

              else if(zph8s(i,j,k).ge.2800.e0                           &
     &           .and.zph8s(i,j,k).lt.6000.e0) then

                qsumm(i,j)=qsumm(i,j)+qsum(i,j,k)

              else

                qsumh(i,j)=qsumh(i,j)+qsum(i,j,k)

              end if

            end do
            end do

!$omp end do

          end do

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            cdl(i,j)=min(tend7*(1.e0-exp(-15.e0*qsuml(i,j))),1.e0)
            cdm(i,j)=min(tend6*(1.e0-exp(-15.e0*qsumm(i,j))),1.e0)
            cdh(i,j)=min(tend3*(1.e0-exp(-15.e0*qsumh(i,j))),1.e0)
          end do
          end do

!$omp end do

        end if

! -----

      end if

!! -----

!$omp end parallel

!!! -----

      end subroutine s_cloudcov

!-----7--------------------------------------------------------------7--

      end module m_cloudcov
