!***********************************************************************
      module m_aggregat
!***********************************************************************

!     Author      : Sakakibara Atsushi, Naito Daisuke
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2000/10/18, 2000/11/17, 2001/05/29,
!                   2001/06/29, 2001/10/18, 2001/11/20, 2002/01/15,
!                   2002/04/02, 2002/07/03, 2001/12/06, 2002/12/02,
!                   2003/04/30, 2003/05/19, 2003/12/12, 2004/03/22,
!                   2004/04/01, 2004/04/15, 2004/07/10, 2004/09/01,
!                   2004/09/10, 2004/09/25, 2004/12/17, 2005/01/07,
!                   2005/11/22, 2006/04/03, 2007/10/19, 2007/11/26,
!                   2008/04/17, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2011/03/18, 2011/03/29, 2011/04/06, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the aggregation rate for the cloud water, rain water,
!     cloud ice and snow.

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

      public :: aggregat, s_aggregat

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface aggregat

        module procedure s_aggregat

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
      subroutine s_aggregat(cphopt,dtb,thresq,ni,nj,nk,rbr,rbv,         &
     &                      qc,qr,qi,qs,ncc,ncr,nci,ncs,diaqc,diaqr,    &
     &                      agcn,agrn,agin,agsn)
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

      real, intent(in) :: qr(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixing ratio

      real, intent(in) :: qi(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixing ratio

      real, intent(in) :: qs(0:ni+1,0:nj+1,1:nk)
                       ! Snow mixing ratio

      real, intent(in) :: ncc(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud water

      real, intent(in) :: ncr(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of rain water

      real, intent(in) :: nci(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud ice

      real, intent(in) :: ncs(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of snow

      real, intent(in) :: diaqc(0:ni+1,0:nj+1,1:nk)
                       ! Mean diameter of cloud water

      real, intent(in) :: diaqr(0:ni+1,0:nj+1,1:nk)
                       ! Mean diameter of rain water

! Output variables

      real, intent(out) :: agcn(0:ni+1,0:nj+1,1:nk)
                       ! Aggregation rate for cloud water

      real, intent(out) :: agrn(0:ni+1,0:nj+1,1:nk)
                       ! Aggregation rate for rain water

      real, intent(out) :: agin(0:ni+1,0:nj+1,1:nk)
                       ! Aggregation rate for cloud ice

      real, intent(out) :: agsn(0:ni+1,0:nj+1,1:nk)
                       ! Aggregation rate for snow

! Internal shared variables

      real expo1       ! (2.0 + bus) / 3.0
      real expo2       ! (4.0 - bus) / 3.0

      real cagcn       ! Coefficient of aggregation rate for cloud water
      real cagrn1      ! Coefficient of aggregation rate for rain water
      real cagrn2      ! Coefficient of aggregation rate for rain water

      real cagin       ! Coefficient of aggregation rate for cloud ice
      real cagsn       ! Coefficient of aggregation rate for snow

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real ercol       ! Bulk collection efficiency for rain water

      real diaqc3      ! diaqc^3
!ORIG real diaqr3      ! diaqr^3

!-----7--------------------------------------------------------------7--

! Set the common used variables.

      expo1=oned3*(2.e0+bus)
      expo2=oned3*(4.e0-bus)

      cagcn=2.59e15*gf7*dtb

      cagrn1=2.59e15*gf7*dtb
      cagrn2=3.03e3*gf4*dtb

      cagin=r0*exp(3.e0*log(.5e0*aui*eii*xl/rhoi*dtb))

      cagsn=3.47222e-4*aus*ess*hfbus*dtb                                &
     &  *exp(oned3*(1.e0-bus)*log(cc))*exp(-oned3*(2.e0+bus)*log(rhos))

! -----

!!!! Calculate the aggregation rate.

!$omp parallel default(shared) private(k)

!!! In the case nk = 1.

      if(nk.eq.1) then

!! Perform calculating in the case the option abs(cphopt) is equal to 2.

        if(abs(cphopt).eq.2) then

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1

! Calculate the aggregation rate for cloud ice.

            if(qi(i,j,1).gt.thresq) then

              agin(i,j,1)=rbr(i,j,1)*nci(i,j,1)                         &
     &          *qi(i,j,1)*exp(oned3*log(cagin*rbv(i,j,1)))

            else

              agin(i,j,1)=0.e0

            end if

! -----

          end do
          end do

!$omp end do

!! -----

!! Perform calculating in the case the option abs(cphopt) is equal to 3.

        else if(abs(cphopt).eq.3) then

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1

! Calculate the aggregation rate for cloud ice.

            if(qi(i,j,1).gt.thresq) then

              agin(i,j,1)=rbr(i,j,1)*nci(i,j,1)                         &
     &          *qi(i,j,1)*exp(oned3*log(cagin*rbv(i,j,1)))

            else

              agin(i,j,1)=0.e0

            end if

! -----

! Calculate the aggregation rate for snow.

            if(qs(i,j,1).gt.thresq) then

              agsn(i,j,1)=cagsn*rbr(i,j,1)                              &
     &          *exp(expo1*log(qs(i,j,1)))*exp(expo2*log(ncs(i,j,1)))

            else

              agsn(i,j,1)=0.e0

            end if

! -----

          end do
          end do

!$omp end do

!! -----

!! Perform calculating in the case the option abs(cphopt) is equal to 4.

        else if(abs(cphopt).eq.4) then

!ORIG!$omp do schedule(runtime) private(i,j,diaqc3,diaqr3,ercol)
!$omp do schedule(runtime) private(i,j,diaqc3,ercol)

          do j=1,nj-1
          do i=1,ni-1

! Calculate the aggregation rate for cloud water.

            if(qc(i,j,1).gt.thresq) then

              diaqc3=diaqc(i,j,1)*diaqc(i,j,1)*diaqc(i,j,1)

              agcn(i,j,1)=cagcn                                         &
     &          *diaqc3*diaqc3*ncc(i,j,1)*ncc(i,j,1)*rbr(i,j,1)

            else

              agcn(i,j,1)=0.e0

            end if

! -----

! Calculate the aggregation rate for rain water.

            if(qr(i,j,1).gt.thresq) then

!ORIG         diaqr3=diaqr(i,j,1)*diaqr(i,j,1)*diaqr(i,j,1)

              if(diaqr(i,j,1).lt.6.e-4) then

                ercol=1.e0

              else if(diaqr(i,j,1).ge.6.e-4                             &
     &           .and.diaqr(i,j,1).lt.2.e-3) then

                ercol=exp(-2.5e3*(diaqr(i,j,1)-6.e-4))

              else

                ercol=0.e0

              end if

!ORIG         if(diaqr(i,j,1).lt.1.e-4) then

!ORIG           agrn(i,j,1)=cagrn1                                      &
!ORIG&            *diaqr3*diaqr3*ncr(i,j,1)*ncr(i,j,1)*rbr(i,j,1)

!ORIG         else

!ORIG           agrn(i,j,1)=cagrn2                                      &
!ORIG&            *diaqr3*ercol*ncr(i,j,1)*ncr(i,j,1)*rbr(i,j,1)

!ORIG         end if

              agrn(i,j,1)=cagrn2                                        &
     &          *diaqr(i,j,1)*diaqr(i,j,1)*diaqr(i,j,1)                 &
     &          *ercol*ncr(i,j,1)*ncr(i,j,1)*rbr(i,j,1)

            else

              agrn(i,j,1)=0.e0

            end if

! -----

! Calculate the aggregation rate for cloud ice.

            if(qi(i,j,1).gt.thresq) then

              agin(i,j,1)=rbr(i,j,1)*nci(i,j,1)                         &
     &          *qi(i,j,1)*exp(oned3*log(cagin*rbv(i,j,1)))

            else

              agin(i,j,1)=0.e0

            end if

! -----

! Calculate the aggregation rate for snow.

            if(qs(i,j,1).gt.thresq) then

              agsn(i,j,1)=cagsn*rbr(i,j,1)                              &
     &          *exp(expo1*log(qs(i,j,1)))*exp(expo2*log(ncs(i,j,1)))

            else

              agsn(i,j,1)=0.e0

            end if

! -----

          end do
          end do

!$omp end do

        end if

!! -----

!!! -----

!!! In the case nk > 1.

      else

!! Perform calculating in the case the option abs(cphopt) is equal to 2.

        if(abs(cphopt).eq.2) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1

! Calculate the aggregation rate for cloud ice.

              if(qi(i,j,k).gt.thresq) then

                agin(i,j,k)=rbr(i,j,k)*nci(i,j,k)                       &
     &            *qi(i,j,k)*exp(oned3*log(cagin*rbv(i,j,k)))

              else

                agin(i,j,k)=0.e0

              end if

! -----

            end do
            end do

!$omp end do

          end do

!! -----

!! Perform calculating in the case the option abs(cphopt) is equal to 3.

        else if(abs(cphopt).eq.3) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1

! Calculate the aggregation rate for cloud ice.

              if(qi(i,j,k).gt.thresq) then

                agin(i,j,k)=rbr(i,j,k)*nci(i,j,k)                       &
     &            *qi(i,j,k)*exp(oned3*log(cagin*rbv(i,j,k)))

              else

                agin(i,j,k)=0.e0

              end if

! -----

! Calculate the aggregation rate for snow.

              if(qs(i,j,k).gt.thresq) then

                agsn(i,j,k)=cagsn*rbr(i,j,k)                            &
     &            *exp(expo1*log(qs(i,j,k)))*exp(expo2*log(ncs(i,j,k)))

              else

                agsn(i,j,k)=0.e0

              end if

! -----

            end do
            end do

!$omp end do

          end do

!! -----

!! Perform calculating in the case the option abs(cphopt) is equal to 4.

        else if(abs(cphopt).eq.4) then

          do k=1,nk-1

!ORIG!$omp do schedule(runtime) private(i,j,diaqc3,diaqr3,ercol)
!$omp do schedule(runtime) private(i,j,diaqc3,ercol)

            do j=1,nj-1
            do i=1,ni-1

! Calculate the aggregation rate for cloud water.

              if(qc(i,j,k).gt.thresq) then

                diaqc3=diaqc(i,j,k)*diaqc(i,j,k)*diaqc(i,j,k)

                agcn(i,j,k)=cagcn                                       &
     &            *diaqc3*diaqc3*ncc(i,j,k)*ncc(i,j,k)*rbr(i,j,k)

              else

                agcn(i,j,k)=0.e0

              end if

! -----

! Calculate the aggregation rate for rain water.

              if(qr(i,j,k).gt.thresq) then

!ORIG           diaqr3=diaqr(i,j,k)*diaqr(i,j,k)*diaqr(i,j,k)

                if(diaqr(i,j,k).lt.6.e-4) then

                  ercol=1.e0

                else if(diaqr(i,j,k).ge.6.e-4                           &
     &             .and.diaqr(i,j,k).lt.2.e-3) then

                  ercol=exp(-2.5e3*(diaqr(i,j,k)-6.e-4))

                else

                  ercol=0.e0

                end if

!ORIG           if(diaqr(i,j,k).lt.1.e-4) then

!ORIG             agrn(i,j,k)=cagrn1                                    &
!ORIG&              *diaqr3*diaqr3*ncr(i,j,k)*ncr(i,j,k)*rbr(i,j,k)

!ORIG           else

!ORIG             agrn(i,j,k)=cagrn2                                    &
!ORIG&              *diaqr3*ercol*ncr(i,j,k)*ncr(i,j,k)*rbr(i,j,k)

!ORIG           end if

                agrn(i,j,k)=cagrn2                                      &
     &            *diaqr(i,j,k)*diaqr(i,j,k)*diaqr(i,j,k)               &
     &            *ercol*ncr(i,j,k)*ncr(i,j,k)*rbr(i,j,k)

              else

                agrn(i,j,k)=0.e0

              end if

! -----

! Calculate the aggregation rate for cloud ice.

              if(qi(i,j,k).gt.thresq) then

                agin(i,j,k)=rbr(i,j,k)*nci(i,j,k)                       &
     &            *qi(i,j,k)*exp(oned3*log(cagin*rbv(i,j,k)))

              else

                agin(i,j,k)=0.e0

              end if

! -----

! Calculate the aggregation rate for snow.

              if(qs(i,j,k).gt.thresq) then

                agsn(i,j,k)=cagsn*rbr(i,j,k)                            &
     &            *exp(expo1*log(qs(i,j,k)))*exp(expo2*log(ncs(i,j,k)))

              else

                agsn(i,j,k)=0.e0

              end if

! -----

            end do
            end do

!$omp end do

          end do

        end if

!! -----

!!! -----

      end if

!$omp end parallel

!!!! -----

      end subroutine s_aggregat

!-----7--------------------------------------------------------------7--

      end module m_aggregat
