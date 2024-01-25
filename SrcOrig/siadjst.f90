!***********************************************************************
      module m_siadjst
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/05/31
!     Modification: 2004/06/10, 2004/09/01, 2004/09/10, 2004/09/25,
!                   2004/10/12, 2004/12/17, 2005/01/07, 2005/01/31,
!                   2005/04/04, 2005/10/05, 2006/02/13, 2006/09/30,
!                   2007/10/19, 2007/11/26, 2008/05/02, 2008/07/01,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform the saturation adjustment for ice.

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

      public :: siadjst, s_siadjst

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface siadjst

        module procedure s_siadjst

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
      subroutine s_siadjst(fpthresq,ni,nj,nk,ptbr,pi,p,ptp,qv,qi,nci)
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

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: pi(0:ni+1,0:nj+1,1:nk)
                       ! Exnar function

      real, intent(in) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure

! Input and output variables

      real, intent(inout) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

      real, intent(inout) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

      real, intent(inout) :: qi(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixing ratio

      real, intent(inout) :: nci(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud ice

! Internal shared variables

      real thresq      ! Minimum threshold value of mixing ratio

      real cwmci       ! cw - ci

      real mi0iv       ! Inverse of mi0

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real t           ! Air temperature

      real tcel        ! Ambient air temperature

      real esi         ! Saturation vapor pressure for ice
      real qvsi        ! Saturation mixing ratio for ice

      real lscpi       ! Latent heat of sublimation / (cp x pi)

      real dqi         ! Variation of cloud ice mixing ratio

      real a           ! Temporary variable
      real b           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getrname(fpthresq,thresq)

! -----

! Set the common used variables.

      cwmci=cw-ci

      mi0iv=1.e0/mi0

! -----

! Perform the saturation adjustment.

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j,t,tcel,esi,qvsi,lscpi,dqi,a,b)

        do j=1,nj-1
        do i=1,ni-1
          t=(ptbr(i,j,k)+ptp(i,j,k))*pi(i,j,k)

          if(t.le.tlow) then

            tcel=t-t0

            a=1.e0/(t-7.66e0)
            b=a*tcel

            esi=es0*exp(21.875e0*b)

            qvsi=epsva*esi/(p(i,j,k)-esi)

            if(qi(i,j,k).gt.thresq.or.qv(i,j,k).gt.qvsi) then

              lscpi=(lv0*exp((.167e0+3.67e-4*t)*log(t0/t))              &
     &          +(lf0+cwmci*tcel))/(cp*pi(i,j,k))

              dqi=(qvsi-qv(i,j,k))                                      &
     &          /(1.e0+21.875e0*a*(1.e0-b)*qvsi*lscpi*pi(i,j,k))

              if(qi(i,j,k).gt.dqi) then

                if(qi(i,j,k).gt.thresq) then

                  nci(i,j,k)=nci(i,j,k)-dqi*nci(i,j,k)/qi(i,j,k)

                else

                  nci(i,j,k)=nci(i,j,k)-dqi*mi0iv

                end if

                ptp(i,j,k)=ptp(i,j,k)-dqi*lscpi

                qv(i,j,k)=qv(i,j,k)+dqi
                qi(i,j,k)=qi(i,j,k)-dqi

              else

                nci(i,j,k)=0.e0

                ptp(i,j,k)=ptp(i,j,k)-qi(i,j,k)*lscpi

                qv(i,j,k)=qv(i,j,k)+qi(i,j,k)
                qi(i,j,k)=0.e0

              end if

            end if

            t=(ptbr(i,j,k)+ptp(i,j,k))*pi(i,j,k)

            if(t.le.tlow) then

              tcel=t-t0

              a=1.e0/(t-7.66e0)
              b=a*tcel

              esi=es0*exp(21.875e0*b)

              qvsi=epsva*esi/(p(i,j,k)-esi)

              if(qi(i,j,k).gt.thresq.or.qv(i,j,k).gt.qvsi) then

                lscpi=(lv0*exp((.167e0+3.67e-4*t)*log(t0/t))            &
     &            +(lf0+cwmci*tcel))/(cp*pi(i,j,k))

                dqi=(qvsi-qv(i,j,k))                                    &
     &            /(1.e0+21.875e0*a*(1.e0-b)*qvsi*lscpi*pi(i,j,k))

                if(qi(i,j,k).gt.dqi) then

                  if(qi(i,j,k).gt.thresq) then

                    nci(i,j,k)=nci(i,j,k)-dqi*nci(i,j,k)/qi(i,j,k)

                  else

                    nci(i,j,k)=nci(i,j,k)-dqi*mi0iv

                  end if

                  ptp(i,j,k)=ptp(i,j,k)-dqi*lscpi

                  qv(i,j,k)=qv(i,j,k)+dqi
                  qi(i,j,k)=qi(i,j,k)-dqi

                else

                  nci(i,j,k)=0.e0

                  ptp(i,j,k)=ptp(i,j,k)-qi(i,j,k)*lscpi

                  qv(i,j,k)=qv(i,j,k)+qi(i,j,k)
                  qi(i,j,k)=0.e0

                end if

              end if

            end if

          end if

        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_siadjst

!-----7--------------------------------------------------------------7--

      end module m_siadjst
