!***********************************************************************
      module m_swadjst
!***********************************************************************

!     Author      : Sakakibara Atsushi, Naito Daisuke
!     Date        : 1999/11/01
!     Modification: 1999/11/19, 1999/12/15, 2000/01/17, 2000/03/08,
!                   2000/04/18, 2000/07/05, 2000/08/21, 2001/04/15,
!                   2001/05/29, 2001/06/29, 2001/08/07, 2001/10/18,
!                   2003/02/13, 2003/03/28, 2003/04/30, 2003/05/19,
!                   2003/12/12, 2004/02/01, 2004/03/05, 2004/03/22,
!                   2004/04/15, 2004/05/31, 2004/06/10, 2004/09/01,
!                   2004/09/10, 2004/09/25, 2004/10/12, 2004/12/17,
!                   2006/02/13, 2006/09/30, 2007/10/19, 2007/11/26,
!                   2008/05/02, 2008/07/01, 2008/08/25, 2009/02/27,
!                   2011/03/29, 2011/04/06, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform the saturation adjustment for water.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comphy
      use m_getgamma
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: swadjst, s_swadjst

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface swadjst

        module procedure s_swadjst

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
      subroutine s_swadjst(fpcphopt,fpthresq,                           &
     &                     ni,nj,nk,ptbr,pi,w,p,ptp,qv,qc,ncc)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

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

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

      real, intent(in) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure

! Input and output variables

      real, intent(inout) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

      real, intent(inout) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

      real, intent(inout) :: qc(0:ni+1,0:nj+1,1:nk)
                       ! Cloud water mixing ratio

      real, intent(inout) :: ncc(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud water

! Internal shared variables

      integer cphopt   ! Option for cloud micro physics

      real thresq      ! Minimum threshold value of mixing ratio

      real mc0iv       ! Inverse of mc0

      real cdiaqc      ! Coefficient of mean diameter of cloud water

      real c1          ! 3.0 x 10^9
      real k1          ! 0.63

      real k1p2iv      ! 1.0 / (k1 + 2.0)

      real k105        ! 0.5 x k1
      real k10515      ! k105 + 1.5

      real beta        ! beta1 x beta2 / beta3

      real beta1       ! Value of Gamma function at k105
      real beta2       ! Value of Gamma function at 1.5
      real beta3       ! Value of Gamma function at k10515

      real coe1        ! g / rd
      real coe2        ! epsva x lv0 / cp
      real coe3        ! epsav / es0
      real coe4        ! coe2 x lv0 / rd
      real coe5        ! rv / (dv0 x es0)
      real coe6        ! lv0 / kpa
      real coe7        ! lv0 / kpa x lv0 / rv
      real coe8        ! 2.0 x cc x rhow
      real coe9        ! 2.0 x cc x k1 x c1 x beta

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real t           ! Air temperature

      real esw         ! Saturation vapor pressure for water
      real qvsw        ! Saturation mixing ratio for water

      real lvcpi       ! Latent heat of evaporation / (cp x pi)

      real dqc         ! Variation of cloud water mixing ratio

      real tiv         ! Inverse of air temperature

      real phi1        ! Temporary variable
      real phi2        ! Temporary variable

      real des         ! Temporary variable
      real nca         ! Temporary variable

      real a           ! Temporary variable
      real b           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpcphopt,cphopt)
      call getrname(fpthresq,thresq)

! -----

! Set the common used variables.

      mc0iv=1.e0/mc0
!MOD  mc0iv=.001e0/mc0

      cdiaqc=1.e0/(cc*rhow)

      c1=3.e9
      k1=.63e0

      k1p2iv=1.e0/(k1+2.e0)

      k105=.5e0*k1
      k10515=k105+1.5e0

      call getgamma(k105,beta1)
      call getgamma(1.5e0,beta2)
      call getgamma(k10515,beta3)

      beta=beta1*beta2/beta3

      coe1=g/rd
      coe2=epsva*lv0/cp
      coe3=epsav/es0
      coe4=coe2*lv0/rd
      coe5=rv/(dv0*es0)
      coe6=lv0/kpa
      coe7=lv0*lv0/(kpa*rv)
      coe8=2.e0*cc*rhow
      coe9=2.e0*cc*c1*k1*beta

! -----

!! Perform the saturation adjustment.

!$omp parallel default(shared) private(k)

! Perform calculating in the case the option abs(cphopt) is less than 3.

      if(abs(cphopt).le.3) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j,t,esw,qvsw,lvcpi,dqc,a,b)

          do j=1,nj-1
          do i=1,ni-1
            t=(ptbr(i,j,k)+ptp(i,j,k))*pi(i,j,k)

            a=1.e0/(t-35.86e0)
            b=a*(t-t0)

            esw=es0*exp(17.269e0*b)

            qvsw=epsva*esw/(p(i,j,k)-esw)

            if(qc(i,j,k).gt.thresq.or.qv(i,j,k).gt.qvsw) then

              lvcpi=lv0*exp((.167e0+3.67e-4*t)*log(t0/t))/(cp*pi(i,j,k))

              dqc=(qvsw-qv(i,j,k))                                      &
     &          /(1.e0+17.269e0*a*(1.e0-b)*qvsw*lvcpi*pi(i,j,k))

              if(qc(i,j,k).gt.dqc) then

                ptp(i,j,k)=ptp(i,j,k)-dqc*lvcpi

                qv(i,j,k)=qv(i,j,k)+dqc
                qc(i,j,k)=qc(i,j,k)-dqc

              else

                ptp(i,j,k)=ptp(i,j,k)-qc(i,j,k)*lvcpi

                qv(i,j,k)=qv(i,j,k)+qc(i,j,k)
                qc(i,j,k)=0.e0

              end if

              t=(ptbr(i,j,k)+ptp(i,j,k))*pi(i,j,k)

              a=1.e0/(t-35.86e0)
              b=a*(t-t0)

              esw=es0*exp(17.269e0*b)

              qvsw=epsva*esw/(p(i,j,k)-esw)

              if(qc(i,j,k).gt.thresq.or.qv(i,j,k).gt.qvsw) then

                lvcpi=lv0                                               &
     &            *exp((.167e0+3.67e-4*t)*log(t0/t))/(cp*pi(i,j,k))

                dqc=(qvsw-qv(i,j,k))                                    &
     &            /(1.e0+17.269e0*a*(1.e0-b)*qvsw*lvcpi*pi(i,j,k))

                if(qc(i,j,k).gt.dqc) then

                  ptp(i,j,k)=ptp(i,j,k)-dqc*lvcpi

                  qv(i,j,k)=qv(i,j,k)+dqc
                  qc(i,j,k)=qc(i,j,k)-dqc

                else

                  ptp(i,j,k)=ptp(i,j,k)-qc(i,j,k)*lvcpi

                  qv(i,j,k)=qv(i,j,k)+qc(i,j,k)
                  qc(i,j,k)=0.e0

                end if

              end if

            end if

          end do
          end do

!$omp end do

        end do

! -----

! Perform calculating in the case the option abs(cphopt) is equal to 4.

      else

        do k=1,nk-1

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,t,esw,qvsw,lvcpi,dqc,tiv,phi1,phi2,des,nca,a,b)

          do j=1,nj-1
          do i=1,ni-1
            t=(ptbr(i,j,k)+ptp(i,j,k))*pi(i,j,k)

            a=1.e0/(t-35.86e0)
            b=a*(t-t0)

            esw=es0*exp(17.269e0*b)

            qvsw=epsva*esw/(p(i,j,k)-esw)

            if(qc(i,j,k).gt.thresq.or.qv(i,j,k).gt.qvsw) then

              lvcpi=lv0*exp((.167e0+3.67e-4*t)*log(t0/t))/(cp*pi(i,j,k))

              dqc=(qvsw-qv(i,j,k))                                      &
     &          /(1.e0+17.269e0*a*(1.e0-b)*qvsw*lvcpi*pi(i,j,k))

              if(qc(i,j,k).gt.dqc) then

                tiv=1.e0/t

                phi1=coe1*tiv                                           &
     &            *(coe2*tiv-1.e0)*(rhow*coe5*t+tiv*(coe7*tiv-coe6))

                phi2=coe3*p(i,j,k)+coe4*tiv*tiv

                if(qc(i,j,k).gt.thresq) then

                  if(dqc.gt.0.e0) then

                    ncc(i,j,k)=ncc(i,j,k)-dqc*ncc(i,j,k)/qc(i,j,k)

                    ptp(i,j,k)=ptp(i,j,k)-dqc*lvcpi

                    qv(i,j,k)=qv(i,j,k)+dqc
                    qc(i,j,k)=qc(i,j,k)-dqc

                  else

                    if(w(i,j,k).gt.0.e0) then

                      des=phi1*w(i,j,k)/(coe8*phi2*ncc(i,j,k)           &
     &                  *exp(oned3*log(cdiaqc*qc(i,j,k)/ncc(i,j,k))))

                      nca=c1*exp(k1*log(des))

                      if(ncc(i,j,k).gt.nca) then
                        ncc(i,j,k)=nca
                      end if

                      ptp(i,j,k)=ptp(i,j,k)-dqc*lvcpi

                      qv(i,j,k)=qv(i,j,k)+dqc
                      qc(i,j,k)=qc(i,j,k)-dqc

                    end if

                  end if

                else

                  if(w(i,j,k).gt.0.e0) then

                    des=exp(k1p2iv                                      &
     &                *log(exp(1.5e0*log(phi1*w(i,j,k)))/(coe9*phi2)))

                    ncc(i,j,k)=ncc(i,j,k)+c1*exp(k1*log(des))

                    ptp(i,j,k)=ptp(i,j,k)-dqc*lvcpi

                    qv(i,j,k)=qv(i,j,k)+dqc
                    qc(i,j,k)=qc(i,j,k)-dqc

!ORIG             else

!ORIG               ncc(i,j,k)=ncc(i,j,k)-dqc*mc0iv

!ORIG               ptp(i,j,k)=ptp(i,j,k)-dqc*lvcpi

!ORIG               qv(i,j,k)=qv(i,j,k)+dqc
!ORIG               qc(i,j,k)=qc(i,j,k)-dqc

                  end if

                end if

              else

                ncc(i,j,k)=0.e0

                ptp(i,j,k)=ptp(i,j,k)-qc(i,j,k)*lvcpi

                qv(i,j,k)=qv(i,j,k)+qc(i,j,k)
                qc(i,j,k)=0.e0

              end if

              t=(ptbr(i,j,k)+ptp(i,j,k))*pi(i,j,k)

              a=1.e0/(t-35.86e0)
              b=a*(t-t0)

              esw=es0*exp(17.269e0*b)

              qvsw=epsva*esw/(p(i,j,k)-esw)

              if(qc(i,j,k).gt.thresq.or.qv(i,j,k).gt.qvsw) then

                lvcpi=lv0                                               &
     &            *exp((.167e0+3.67e-4*t)*log(t0/t))/(cp*pi(i,j,k))

                dqc=(qvsw-qv(i,j,k))                                    &
     &            /(1.e0+17.269e0*a*(1.e0-b)*qvsw*lvcpi*pi(i,j,k))

                if(qc(i,j,k).gt.dqc) then

                  tiv=1.e0/t

                  phi1=coe1*tiv                                         &
     &              *(coe2*tiv-1.e0)*(rhow*coe5*t+tiv*(coe7*tiv-coe6))

                  phi2=coe3*p(i,j,k)+coe4*tiv*tiv

                  if(qc(i,j,k).gt.thresq) then

                    if(dqc.gt.0.e0) then

                      ncc(i,j,k)=ncc(i,j,k)-dqc*ncc(i,j,k)/qc(i,j,k)

                      ptp(i,j,k)=ptp(i,j,k)-dqc*lvcpi

                      qv(i,j,k)=qv(i,j,k)+dqc
                      qc(i,j,k)=qc(i,j,k)-dqc

                    else

                      if(w(i,j,k).gt.0.e0) then

                        des=phi1*w(i,j,k)/(coe8*phi2*ncc(i,j,k)         &
     &                    *exp(oned3*log(cdiaqc*qc(i,j,k)/ncc(i,j,k))))

                        nca=c1*exp(k1*log(des))

                        if(ncc(i,j,k).gt.nca) then
                          ncc(i,j,k)=nca
                        end if

                        ptp(i,j,k)=ptp(i,j,k)-dqc*lvcpi

                        qv(i,j,k)=qv(i,j,k)+dqc
                        qc(i,j,k)=qc(i,j,k)-dqc

                      end if

                    end if

                  else

                    if(w(i,j,k).gt.0.e0) then

                      des=exp(k1p2iv                                    &
     &                  *log(exp(1.5e0*log(phi1*w(i,j,k)))/(coe9*phi2)))

                      ncc(i,j,k)=ncc(i,j,k)+c1*exp(k1*log(des))

                      ptp(i,j,k)=ptp(i,j,k)-dqc*lvcpi

                      qv(i,j,k)=qv(i,j,k)+dqc
                      qc(i,j,k)=qc(i,j,k)-dqc

!ORIG               else

!ORIG                 ncc(i,j,k)=ncc(i,j,k)-dqc*mc0iv

!ORIG                 ptp(i,j,k)=ptp(i,j,k)-dqc*lvcpi

!ORIG                 qv(i,j,k)=qv(i,j,k)+dqc
!ORIG                 qc(i,j,k)=qc(i,j,k)-dqc

                    end if

                  end if

                else

                  ncc(i,j,k)=0.e0

                  ptp(i,j,k)=ptp(i,j,k)-qc(i,j,k)*lvcpi

                  qv(i,j,k)=qv(i,j,k)+qc(i,j,k)
                  qc(i,j,k)=0.e0

                end if

              end if

            end if

          end do
          end do

!$omp end do

        end do

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_swadjst

!-----7--------------------------------------------------------------7--

      end module m_swadjst
