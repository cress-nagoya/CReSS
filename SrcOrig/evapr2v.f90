!***********************************************************************
      module m_evapr2v
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/11/01
!     Modification: 1999/11/19, 2000/01/17, 2000/03/08, 2000/04/18,
!                   2000/07/05, 2000/08/21, 2001/05/29, 2001/06/29,
!                   2001/10/18, 2001/11/20, 2002/01/15, 2002/04/02,
!                   2003/03/21, 2003/04/30, 2003/05/19, 2003/09/01,
!                   2003/12/12, 2004/03/22, 2004/04/15, 2004/09/01,
!                   2004/09/10, 2004/09/25, 2004/10/12, 2004/12/17,
!                   2005/01/31, 2006/02/13, 2006/04/03, 2007/10/19,
!                   2008/05/02, 2008/07/01, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the evaporation rate from the rain water to the water
!     vapor.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: evapr2v, s_evapr2v

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface evapr2v

        module procedure s_evapr2v

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
      subroutine s_evapr2v(fpthresq,dtb,ni,nj,nk,ptbr,rbr,pi,p,         &
     &                     ptpp,qvp,qrp,ptpf,qvf,qrf)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpthresq
                       ! Formal parameter of unique index of thresq

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: pi(0:ni+1,0:nj+1,1:nk)
                       ! Exner function

      real, intent(in) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure

      real, intent(in) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(in) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(in) :: qrp(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixing ratio at past

! Input and output variables

      real, intent(inout) :: ptpf(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at future

      real, intent(inout) :: qvf(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at future

      real, intent(inout) :: qrf(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixing ratio at future

! Internal shared variable

      real thresq      ! Minimum threshold value of mixing ratio

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real rbrqr       ! rbr x qr

      real t           ! Air temperature

      real esw         ! Saturation vapor pressure for water
      real qvsw        ! Saturation mixing ratio for water

      real vent        ! Ventilation factor

      real lvcpi       ! Latent heat of evapolation / (cp x pi)

      real evrv        ! Evaporation rate from rain water to water vapor

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getrname(fpthresq,thresq)

! -----

! Calculate the evaporation rate from the rain water to the water vapor.

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j,rbrqr,t,esw,qvsw,vent,lvcpi,evrv)

        do j=1,nj-1
        do i=1,ni-1

          if(qrp(i,j,k).gt.thresq) then

            t=(ptbr(i,j,k)+ptpp(i,j,k))*pi(i,j,k)

            esw=es0*exp(17.269e0*(t-t0)/(t-35.86e0))

            qvsw=epsva*esw/(p(i,j,k)-esw)

            if(qvp(i,j,k).lt.qvsw) then

              rbrqr=rbr(i,j,k)*qrp(i,j,k)

              vent=1.6e0+30.3922e0*exp(.2046e0*log(rbrqr))

              evrv=(1.e0-qvp(i,j,k)/qvsw)                               &
     &          *vent*exp(.525e0*log(rbrqr))*dtb                        &
     &          /((2.03e4+9.584e6/(p(i,j,k)*qvsw))*rbr(i,j,k))

              lvcpi=lv0*exp((.167e0+3.67e-4*t)*log(t0/t))/(cp*pi(i,j,k))

              if(qrf(i,j,k).gt.evrv) then

                ptpf(i,j,k)=ptpf(i,j,k)-evrv*lvcpi

                qvf(i,j,k)=qvf(i,j,k)+evrv
                qrf(i,j,k)=qrf(i,j,k)-evrv

              else

                ptpf(i,j,k)=ptpf(i,j,k)-qrf(i,j,k)*lvcpi

                qvf(i,j,k)=qvf(i,j,k)+qrf(i,j,k)
                qrf(i,j,k)=0.e0

              end if

            end if

          end if

        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_evapr2v

!-----7--------------------------------------------------------------7--

      end module m_evapr2v
