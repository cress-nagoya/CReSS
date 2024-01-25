!***********************************************************************
      module m_bruntv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/11/24
!     Modification: 2000/01/17, 2000/03/08, 2000/04/18, 2000/07/05,
!                   2001/04/15, 2001/05/29, 2001/06/29, 2001/08/07,
!                   2001/10/18, 2001/12/11, 2002/01/15, 2002/04/02,
!                   2002/12/02, 2003/02/13, 2003/03/13, 2003/04/30,
!                   2003/05/19, 2003/11/28, 2003/12/12, 2004/02/01,
!                   2004/03/05, 2004/04/15, 2004/05/31, 2004/07/01,
!                   2004/08/01, 2004/09/01, 2004/09/10, 2005/04/04,
!                   2005/06/10, 2005/09/01, 2006/01/10, 2006/02/13,
!                   2006/04/03, 2006/10/20, 2007/01/20, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the Brunt-Vaisala frequency squared.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: bruntv, s_bruntv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface bruntv

        module procedure s_bruntv

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic exp
      intrinsic log

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_bruntv(fpcphopt,fpdziv,fpthresq,fmois,ni,nj,nk,jcb8w,&
     &                    pbr,ptbr,pp,ptp,qv,qall,nsq8w,pt,ptv,a,t)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fpdziv
                       ! Formal parameter of unique index of dziv

      integer, intent(in) :: fpthresq
                       ! Formal parameter of unique index of thresq

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation

      real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

      real, intent(in) :: qall(0:ni+1,0:nj+1,1:nk)
                       ! Total water and ice mixing ratio

! Output variable

      real, intent(out) :: nsq8w(0:ni+1,0:nj+1,1:nk)
                       ! Half value of Brunt Vaisala frequency squared
                       ! at w points

! Internal shared variables

      integer cphopt   ! Option for cloud micro physics

      integer nkm1     ! nk - 1

      real dziv        ! Inverse of dz

      real thresq      ! Minimum threshold value of mixing ratio

      real gdzv        ! g x dziv
      real gdzv05      ! 0.5 x g x dziv

      real rddvcp      ! rd / cp

      real cwmci       ! cw - ci

      real p0iv        ! 1.0 / p0

      real, intent(inout) :: pt(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature

      real, intent(inout) :: ptv(0:ni+1,0:nj+1,1:nk)
                       ! Virtual potential temperature

      real, intent(inout) :: a(0:ni+1,0:nj+1,1:nk)
                       ! Complicated coefficient table

      real, intent(inout) :: t(0:ni+1,0:nj+1)
                       ! Air temperature

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real lhcpt       ! Latent heat / (cp x t)

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpcphopt,cphopt)
      call getrname(fpdziv,dziv)
      call getrname(fpthresq,thresq)

! -----

! Set the common used variables.

      nkm1=nk-1

      gdzv=g*dziv
      gdzv05=.5e0*g*dziv

      rddvcp=rd/cp

      cwmci=cw-ci

      p0iv=1.e0/p0

! -----

!!! Calculate the Brunt-Vaisala frequency squared.

!$omp parallel default(shared) private(k)

! Get the potential temperature.

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          pt(i,j,k)=ptbr(i,j,k)+ptp(i,j,k)
        end do
        end do

!$omp end do

      end do

! -----

! Calculate the air stability for the dry air.

      if(fmois(1:3).eq.'dry') then

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            nsq8w(i,j,k)=gdzv*(pt(i,j,k)-pt(i,j,k-1))                   &
     &        /(jcb8w(i,j,k)*(ptbr(i,j,k-1)+ptbr(i,j,k)))
          end do
          end do

!$omp end do

        end do

! -----

!! Calculate the air stability for the moist air.

      else if(fmois(1:5).eq.'moist') then

! Get the virtual potential temperature.

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            ptv(i,j,k)=pt(i,j,k)*(1.e0+epsav*qv(i,j,k))/(1.e0+qv(i,j,k))
          end do
          end do

!$omp end do

        end do

! -----

! In the case of no cloud micro physics.

        if(abs(cphopt).eq.0) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              nsq8w(i,j,k)=gdzv*(ptv(i,j,k)-ptv(i,j,k-1))               &
     &          /(jcb8w(i,j,k)*(ptbr(i,j,k-1)+ptbr(i,j,k)))
            end do
            end do

!$omp end do

          end do

! -----

! In the case of performing cloud micro physics.

        else

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              t(i,j)=pt(i,j,k)                                          &
     &          *exp(rddvcp*log(p0iv*(pbr(i,j,k)+pp(i,j,k))))

              a(i,j,k)=lv0                                              &
     &          *exp((.167e0+3.67e-4*t(i,j))*log(t0/t(i,j)))

            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1

              if(t(i,j).le.tlow) then
                a(i,j,k)=a(i,j,k)+(lf0+cwmci*(t(i,j)-t0))
              end if

            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j,lhcpt)

            do j=1,nj-1
            do i=1,ni-1
              lhcpt=a(i,j,k)/(cp*t(i,j))

              pt(i,j,k)=pt(i,j,k)*exp(lhcpt*qv(i,j,k))

              a(i,j,k)=a(i,j,k)*qv(i,j,k)/(rd*t(i,j))
              a(i,j,k)=(1.e0+a(i,j,k))/(1.e0+epsva*lhcpt*a(i,j,k))
            end do
            end do

!$omp end do

          end do

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1

              if(qall(i,j,k).gt.thresq) then

                nsq8w(i,j,k)=gdzv05*((qall(i,j,k-1)-qall(i,j,k))        &
     &            +(pt(i,j,k)-pt(i,j,k-1))*(a(i,j,k-1)+a(i,j,k))        &
     &            /(ptbr(i,j,k-1)+ptbr(i,j,k)))/jcb8w(i,j,k)

              else

                nsq8w(i,j,k)=gdzv*(ptv(i,j,k)-ptv(i,j,k-1))             &
     &            /(jcb8w(i,j,k)*(ptbr(i,j,k-1)+ptbr(i,j,k)))

              end if

            end do
            end do

!$omp end do

          end do

        end if

! -----

      end if

!! -----

! Set the bottom and top boundary conditions.

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=1,ni-1
        nsq8w(i,j,1)=nsq8w(i,j,2)
        nsq8w(i,j,nk)=nsq8w(i,j,nkm1)
      end do
      end do

!$omp end do

! -----

!$omp end parallel

!!! -----

      end subroutine s_bruntv

!-----7--------------------------------------------------------------7--

      end module m_bruntv
