!***********************************************************************
      module m_inidisbw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/03/06
!     Modification: 2006/08/08, 2006/09/30, 2007/03/29, 2007/10/19,
!                   2008/01/11, 2008/05/02, 2008/08/25, 2009/01/30,
!                   2009/02/27, 2011/03/18, 2011/08/18, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the initial bin distribution of total water.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comphy
      use m_comtable
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: inidisbw, s_inidisbw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface inidisbw

        module procedure s_inidisbw

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic log

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_inidisbw(fpdziv,dtb,ni,nj,nk,nqw,nnw,jcb8w,rbv,pi,   &
     &                      w,t,brw,rbrw,ptptmp,qvtmp,qwtmp,ptp,qv,     &
     &                      mwbin,nwbin,c,d,d1,d2,d3,d4,cd1,ccd1,       &
     &                      gi,fi)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpdziv
                       ! Formal parameter of unique index of dziv

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

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: rbv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of base state density [cm^3/g]

      real, intent(in) :: pi(0:ni+1,0:nj+1,1:nk)
                       ! Exnar function

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

      real, intent(in) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

      real, intent(in) :: brw(1:nqw+1)
                       ! Radius of water bin boundary [cm]

      real, intent(in) :: rbrw(1:nqw+1,1:8)
                       ! Related parameters of brw

      real, intent(in) :: ptptmp(0:ni+1,0:nj+1,1:nk)
                       ! Temporary alternative
                       ! potential temperature perturbation

      real, intent(in) :: qvtmp(0:ni+1,0:nj+1,1:nk)
                       ! Temporary alternative
                       ! water vapor mixing ratio [g/g]

      real, intent(in) :: qwtmp(0:ni+1,0:nj+1,1:nk)
                       ! Temporary alternative
                       ! total water mixing ratio [g/g]

! Input and output variables

      real, intent(inout) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

      real, intent(inout) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio [g/g]

      real, intent(inout) :: mwbin(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Total water mass [g/cm^3]

      real, intent(inout) :: nwbin(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations [1/cm^3]

! Internal shared variables

      real dziv        ! Inverse of dz

      real qcciv       ! 1.0 / qccrit

      real cc40        ! 40.0 x cc
      real cc80        ! 80.0 x cc

      real cc43        ! 4.0 / 3.0 x cc

      real rwiv3       ! 1000.0 / rhow

      real dzvdt2      ! 0.5 x dziv x dtb

      real, intent(inout) :: c(0:ni+1,0:nj+1,1:nk)
                       ! Coefficient of 2nd order Gamma distributions

      real, intent(inout) :: d(0:ni+1,0:nj+1)
                       ! Coefficient of 2nd order Gamma distributions

      real, intent(inout) :: d1(0:ni+1,0:nj+1)
                       ! 1.0 / d

      real, intent(inout) :: d2(0:ni+1,0:nj+1)
                       ! 2.0 x d1 x d1

      real, intent(inout) :: d3(0:ni+1,0:nj+1)
                       ! 30.0 x d1 x d2

      real, intent(inout) :: d4(0:ni+1,0:nj+1)
                       ! 2.0 x d1 x d3

      real, intent(inout) :: cd1(0:ni+1,0:nj+1)
                       ! rwiv3 x c x d1

      real, intent(inout) :: ccd1(0:ni+1,0:nj+1)
                       ! cc43 x cd1

      real, intent(inout) :: gi(0:ni+1,0:nj+1,1:nqw+1)
                       ! Temporary array

      real, intent(inout) :: fi(0:ni+1,0:nj+1,1:nqw+1)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer n        ! Array index in water bin categories

      real cexp        ! exp(- d x brw)

      real nd          ! Total surturated number concentrations

      real lvcpi       ! Latent heat of evaporation / (cp x pi)

      real a           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getrname(fpdziv,dziv)

! -----

! Set the common used variable.

      qcciv=1.e0/qccrit

      cc40=40.e0*cc
      cc80=80.e0*cc

      cc43=fourd3*cc

      rwiv3=1000.e0/rhow

      dzvdt2=.5e0*dziv*dtb

! -----

!!! Set the initial bin distribution of total water, calculate the new
!!! potential temperature and water vapor mixing ratio and set the
!!! bottom boundary conditions.

!$omp parallel default(shared) private(k,n)

!! Set the initial bin distribution of total water.

      do k=2,nk-1

! Calculate the cloud condensation nuclei and the coefficient of
! 2nd order Gamma distributions.

!$omp do schedule(runtime) private(i,j,nd,a)

        do j=1,nj-1
        do i=1,ni-1

          if(qwtmp(i,j,k).gt.qccrit.or.qwtmp(i,j,k).lt.-10.e0) then

            if(qwtmp(i,j,k).gt.qccrit) then

              a=rbv(i,j,k)/qwtmp(i,j,k)

            else

              if(qv(i,j,k).gt.qccrit) then
                a=rbv(i,j,k)*qcciv
              else
                a=rbv(i,j,k)/qv(i,j,k)
              end if

            end if

            if(w(i,j,k).le..24e0) then

              nd=4710.e0*exp(1.19e0*log(w(i,j,k)+eps))*nccn(1)          &
     &          /(nccn(1)+(1090.e0*w(i,j,k)+33.2e0))

            else if(w(i,j,k).gt..24e0.and.w(i,j,k).le..5e0) then

              nd=(11700.e0*w(i,j,k)-1690.e0)*nccn(2)                    &
     &          /(nccn(2)+(10600.e0*w(i,j,k)-1480.e0))

            else if(w(i,j,k).gt..5e0.and.w(i,j,k).le.1.e0) then

              nd=4300.e0*exp(1.05e0*log(w(i,j,k)))*nccn(3)              &
     &          /(nccn(3)+(2760.e0*exp(.755e0*log(w(i,j,k)))))

            else if(w(i,j,k).gt.1.e0.and.w(i,j,k).le.3.e0) then

              nd=(7730.e0-15800.e0*exp(-1.08e0*w(i,j,k)))*nccn(4)       &
     &          /(nccn(4)+(6030.e0-24100.e0*exp(-1.87e0*w(i,j,k))))

            else

              nd=(1140.e0*w(i,j,k)-741.e0)*nccn(5)                      &
     &          /(nccn(5)+(909.e0*w(i,j,k)-56.2e0))

            end if

            c(i,j,k)=a*cc40*nd*nd

            d(i,j)=exp(oned3*log(a*cc80*nd))

            d1(i,j)=1.e0/d(i,j)
            d2(i,j)=2.e0*d1(i,j)*d1(i,j)
            d3(i,j)=30.e0*d1(i,j)*d2(i,j)
            d4(i,j)=2.e0*d1(i,j)*d3(i,j)

            cd1(i,j)=rwiv3*c(i,j,k)*d1(i,j)
            ccd1(i,j)=cc43*cd1(i,j)

          end if

        end do
        end do

!$omp end do

! -----

! Calculate the total mass and number concentrations for each bin.

        do n=1,nqw+1

!$omp do schedule(runtime) private(i,j,cexp)

          do j=1,nj-1
          do i=1,ni-1

            if(qwtmp(i,j,k).gt.qccrit.or.qwtmp(i,j,k).lt.-10.e0) then

              cexp=exp(-d(i,j)*brw(n))

              fi(i,j,n)=cexp*(d2(i,j)+rbrw(n,1)*d1(i,j)+rbrw(n,2))

              gi(i,j,n)=cexp*(d4(i,j)*(brw(n)+d1(i,j))+rbrw(n,2)*d3(i,j)&
     &          +rbrw(n,3)*d2(i,j)+rbrw(n,4)*d1(i,j)+rbrw(n,5))

            end if

          end do
          end do

!$omp end do

        end do

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1

          if(qwtmp(i,j,k).lt.-10.e0) then

            c(i,j,k)=0.e0

            d(i,j)=w(i,j,k)/jcb8w(i,j,k)*dzvdt2

          end if

        end do
        end do

!$omp end do

        do n=1,nqw

!$omp do schedule(runtime) private(i,j,a)

          do j=1,nj-1
          do i=1,ni-1

            if(qwtmp(i,j,k).gt.qccrit) then

             mwbin(i,j,k,n)=mwbin(i,j,k,n)                              &
     &         +ccd1(i,j)*(gi(i,j,n)-gi(i,j,n+1))

             nwbin(i,j,k,n)=nwbin(i,j,k,n)                              &
     &         +cd1(i,j)*(fi(i,j,n)-fi(i,j,n+1))

            else if(qwtmp(i,j,k).lt.-10.e0) then

             a=d(i,j)*(ccd1(i,j)*(gi(i,j,n)-gi(i,j,n+1))-mwbin(i,j,k,n))

             c(i,j,k)=c(i,j,k)+a

             mwbin(i,j,k,n)=mwbin(i,j,k,n)+a

             nwbin(i,j,k,n)=nwbin(i,j,k,n)                              &
     &         -d(i,j)*(nwbin(i,j,k,n)-cd1(i,j)*(fi(i,j,n)-fi(i,j,n+1)))

             a=d(i,j)*(ccd1(i,j)*(gi(i,j,n)-gi(i,j,n+1))-mwbin(i,j,k,n))

             c(i,j,k)=c(i,j,k)+a

             mwbin(i,j,k,n)=mwbin(i,j,k,n)+a

             nwbin(i,j,k,n)=nwbin(i,j,k,n)                              &
     &         -d(i,j)*(nwbin(i,j,k,n)-cd1(i,j)*(fi(i,j,n)-fi(i,j,n+1)))

            end if

          end do
          end do

!$omp end do

        end do

! -----

      end do

!! -----

!! Calculate the new potential temperature and water vapor mixing ratio
!! and set the bottom boundary conditions.

! Calculate the new potential temperature and water vapor mixing ratio.

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j,lvcpi,a)

        do j=1,nj-1
        do i=1,ni-1

          if(qwtmp(i,j,k).gt.qccrit) then

            ptp(i,j,k)=ptptmp(i,j,k)
            qv(i,j,k)=qvtmp(i,j,k)

          else if(qwtmp(i,j,k).lt.-10.e0) then

            lvcpi=lv0*exp((.167e0+3.67e-4*t(i,j,k))*log(t0/t(i,j,k)))
            lvcpi=lvcpi/(cp*pi(i,j,k))

            a=rbv(i,j,k)*c(i,j,k)

            ptp(i,j,k)=ptp(i,j,k)+a*lvcpi
            qv(i,j,k)=qv(i,j,k)-a

          end if

        end do
        end do

!$omp end do

      end do

! -----

! Set the bottom boundary conditions.

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=1,ni-1
        ptp(i,j,1)=ptp(i,j,2)
        qv(i,j,1)=qv(i,j,2)
      end do
      end do

!$omp end do

      do n=1,nqw

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          mwbin(i,j,1,n)=mwbin(i,j,2,n)
        end do
        end do

!$omp end do

      end do

      do n=1,nnw

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          nwbin(i,j,1,n)=nwbin(i,j,2,n)
        end do
        end do

!$omp end do

      end do

! -----

!! -----

!$omp end parallel

!!! -----

      end subroutine s_inidisbw

!-----7--------------------------------------------------------------7--

      end module m_inidisbw
