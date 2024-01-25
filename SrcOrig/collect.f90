!***********************************************************************
      module m_collect
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2000/10/18, 2000/11/17, 2001/04/15,
!                   2001/10/18, 2001/11/20, 2001/12/11, 2002/01/15,
!                   2002/04/02, 2002/07/03, 2002/12/06, 2003/03/28,
!                   2003/04/30, 2003/05/19, 2003/10/31, 2003/11/05,
!                   2003/12/12, 2004/03/22, 2004/04/01, 2004/04/15,
!                   2004/06/10, 2004/07/10, 2004/08/01, 2004/09/01,
!                   2004/09/25, 2004/10/12, 2005/01/07, 2005/01/31,
!                   2005/04/04, 2006/01/10, 2006/03/06, 2006/04/03,
!                   2007/06/27, 2007/10/19, 2007/11/26, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2011/03/18, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the collection rate.

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

      public :: collect, s_collect

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface collect

        module procedure s_collect

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic exp
      intrinsic log
      intrinsic sqrt

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_collect(cphopt,dtb,thresq,ni,nj,nk,                  &
     &                     rbr,qc,qr,qi,qs,qg,ncr,ncs,ncg,              &
     &                     urq,usq,ugq,urn,usn,ugn,t,mu,mi,             &
     &                     diaqc,diaqr,diaqs,diaqg,clcr,clcs,clcg,      &
     &                     clri,clrs,clrg,clir,clis,clig,clsr,clsg,     &
     &                     clrin,clrsn,clsrn,clsgn,ecs)
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

      real, intent(in) :: qc(0:ni+1,0:nj+1,1:nk)
                       ! Cloud water mixing ratio

      real, intent(in) :: qr(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixing ratio

      real, intent(in) :: qi(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixing ratio

      real, intent(in) :: qs(0:ni+1,0:nj+1,1:nk)
                       ! Snow mixing ratio

      real, intent(in) :: qg(0:ni+1,0:nj+1,1:nk)
                       ! Graupel mixing ratio

      real, intent(in) :: ncr(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of rain water

      real, intent(in) :: ncs(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of snow

      real, intent(in) :: ncg(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of graupel

      real, intent(in) :: urq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of rain water

      real, intent(in) :: usq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of snow

      real, intent(in) :: ugq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of graupel

      real, intent(in) :: urn(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of rain water concentrations

      real, intent(in) :: usn(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of snow concentrations

      real, intent(in) :: ugn(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of graupel concentrations

      real, intent(in) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

      real, intent(in) :: mu(0:ni+1,0:nj+1,1:nk)
                       ! Viscosity of air

      real, intent(in) :: mi(0:ni+1,0:nj+1,1:nk)
                       ! Mean mass of cloud ice

      real, intent(in) :: diaqc(0:ni+1,0:nj+1,1:nk)
                       ! Mean diameter of cloud water

      real, intent(in) :: diaqr(0:ni+1,0:nj+1,1:nk)
                       ! Mean diameter of rain water

      real, intent(in) :: diaqs(0:ni+1,0:nj+1,1:nk)
                       ! Mean diameter of snow

      real, intent(in) :: diaqg(0:ni+1,0:nj+1,1:nk)
                       ! Mean diameter of graupel

! Output variables

      real, intent(out) :: clcr(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between cloud water and rain water

      real, intent(out) :: clcs(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between cloud water and snow

      real, intent(out) :: clcg(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between cloud water and graupel

      real, intent(out) :: clri(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between rain water and cloud ice

      real, intent(out) :: clrs(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! from rain water to snow

      real, intent(out) :: clrg(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between rain water and graupel

      real, intent(out) :: clir(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between rain water and cloud ice

      real, intent(out) :: clis(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between cloud ice and snow

      real, intent(out) :: clig(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between cloud ice and graupel

      real, intent(out) :: clsr(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! from snow to rain water

      real, intent(out) :: clsg(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between snow and graupel

      real, intent(out) :: clrin(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate for concentrations
                       ! between rain water and cloud ice

      real, intent(out) :: clrsn(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate for concentrations
                       ! between rain water and snow

      real, intent(out) :: clsrn(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate for concentrations
                       ! between rain water and snow

      real, intent(out) :: clsgn(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate for concentrations
                       ! between snow and graupel

      real, intent(out) :: ecs(0:ni+1,0:nj+1,1:nk)
                       ! Collection efficiency of snow for cloud water

! Internal shared variables

      real rwdv9       ! rhow / 9.0

      real bur2        ! 2.0 + bur
      real bus2        ! 2.0 + bus
      real bug2        ! 2.0 + bug

      real esgiv       ! Inverse of esgiv

      real cclcr       ! Coefficient of collection rate
      real cclcs       ! Coefficient of collection rate
      real cclcg       ! Coefficient of collection rate
      real cclri       ! Coefficient of collection rate
      real cclrs       ! Coefficient of collection rate
      real cclrg       ! Coefficient of collection rate
      real cclis       ! Coefficient of collection rate
      real cclig       ! Coefficient of collection rate
      real cclsr       ! Coefficient of collection rate
      real cclsg       ! Coefficient of collection rate
      real cclsg3      ! Coefficient of collection rate

      real cclrsn      ! Coefficient of collection rate
                       ! for concentrations

      real cclsgn      ! Coefficient of collection rate
                       ! for concentrations

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real cstk        ! rwdv9 x diaqc

      real r0rsq       ! (r0 x rbr)^0.5

      real qr2b        ! r0rsq x exp(bur2 x ln(diaqr))
      real qs2b        ! r0rsq x exp(bus2 x ln(diaqs))
      real qg2b        ! r0rsq x exp(bug2 x ln(diaqg))

      real diaqr2      ! diaqr^2
      real diaqs2      ! diaqs^2
      real diaqg2      ! diaqg^2

      real diaqr3      ! diaqr^3
      real diaqs3      ! diaqs^3

      real sink        ! Sink amount in collection

      real a           ! Temporary variable
      real b           ! Temporary variable
      real c           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Set the common used variables.

      rwdv9=oned9*rhow

      bur2=2.e0+bur
      bus2=2.e0+bus
      bug2=2.e0+bug

      esgiv=1.e0/esg

      cclcr=.25e0*aur*gf3bur*cc*dtb
      cclcs=.25e0*aus*gf3bus*cc*dtb
      cclcg=.25e0*aug*gf3bug*cc*dtb

      cclri=.25e0*aur*gf3bur*eir*cc*dtb
      cclrs=ers*cc*cc*rhow*dtb
      cclrg=erg*cc*cc*rhow*dtb

      cclis=.25e0*aus*gf3bus*eis*cc*dtb
      cclig=.25e0*aug*gf3bug*eig*cc*dtb

      cclsr=esr*cc*cc*rhos*dtb
      cclsg=esg*cc*cc*rhos*dtb
      cclsg3=cc*cc*rhos*dtb

      cclrsn=.5e0*ers*cc*dtb
      cclsgn=.5e0*esg*cc*dtb

! -----

!!!!! Calculate the collection rate.

!$omp parallel default(shared) private(k)

!!!! In the case nk = 1.

      if(nk.eq.1) then

!!! Perform calculating in the case the option abs(cphopt) is equal
!!! to 2.

        if(abs(cphopt).eq.2) then

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,cstk,r0rsq,qr2b,qs2b,qg2b)                         &
!$omp&   private(diaqr2,diaqs2,diaqg2,diaqr3,diaqs3,sink,a,b,c)

          do j=1,nj-1
          do i=1,ni-1

! Set the common used variables.

            r0rsq=sqrt(r0*rbr(i,j,1))

            if(qr(i,j,1).gt.thresq) then

              diaqr2=diaqr(i,j,1)*diaqr(i,j,1)
              diaqr3=diaqr(i,j,1)*diaqr2

              qr2b=r0rsq*exp(bur2*log(diaqr(i,j,1)))

            else

              diaqr2=0.e0
              diaqr3=0.e0

              qr2b=0.e0

            end if

            if(qs(i,j,1).gt.thresq) then

              diaqs2=diaqs(i,j,1)*diaqs(i,j,1)
              diaqs3=diaqs(i,j,1)*diaqs2

              qs2b=r0rsq*exp(bus2*log(diaqs(i,j,1)))

            else

              diaqs2=0.e0
              diaqs3=0.e0

              qs2b=0.e0

            end if

            if(qg(i,j,1).gt.thresq) then

              diaqg2=diaqg(i,j,1)*diaqg(i,j,1)

              qg2b=r0rsq*exp(bug2*log(diaqg(i,j,1)))

            else

              diaqg2=0.e0

              qg2b=0.e0

            end if

! -----

!! Calculate the collection rate between the cloud water and the rain
!! water, snow and graupel and between the rain water and the cloud ice,
!! the snow and the graupel.

            if(t(i,j,1).gt.tlow) then

! Calculate the collection rate between the cloud water and the rain
! water, snow and graupel.

              if(qc(i,j,1).gt.thresq) then

                cstk=rwdv9*diaqc(i,j,1)*diaqc(i,j,1)

                if(qr(i,j,1).gt.thresq) then

                  a=cstk*urq(i,j,1)/(diaqr(i,j,1)*mu(i,j,1))

                  b=a+.5e0

                  a=a*a/(b*b)

                  clcr(i,j,1)=cclcr*a*qr2b*ncr(i,j,1)*qc(i,j,1)

                else

                  clcr(i,j,1)=0.e0

                end if

                if(qs(i,j,1).gt.thresq) then

                  a=cstk*usq(i,j,1)/(diaqs(i,j,1)*mu(i,j,1))

                  b=a+.5e0

                  ecs(i,j,1)=a*a/(b*b)

                  clcs(i,j,1)=cclcs*qs2b*ecs(i,j,1)*ncs(i,j,1)*qc(i,j,1)

                else

                  clcs(i,j,1)=0.e0

                  ecs(i,j,1)=0.e0

                end if

                if(qg(i,j,1).gt.thresq) then

                  a=cstk*ugq(i,j,1)/(diaqg(i,j,1)*mu(i,j,1))

                  b=a+.5e0

                  a=a*a/(b*b)

                  clcg(i,j,1)=cclcg*a*qg2b*ncg(i,j,1)*qc(i,j,1)

                else

                  clcg(i,j,1)=0.e0

                end if

                sink=clcr(i,j,1)+clcs(i,j,1)+clcg(i,j,1)

                if(qc(i,j,1).lt.sink) then

                  a=qc(i,j,1)/sink

                  clcr(i,j,1)=clcr(i,j,1)*a
                  clcs(i,j,1)=clcs(i,j,1)*a
                  clcg(i,j,1)=clcg(i,j,1)*a

                end if

              else

                clcr(i,j,1)=0.e0
                clcs(i,j,1)=0.e0
                clcg(i,j,1)=0.e0

                ecs(i,j,1)=0.e0

              end if

! -----

! Calculate the collection rate between the rain water and the cloud
! ice, the snow and the graupel.

              if(qr(i,j,1).gt.thresq) then

                if(qi(i,j,1).gt.thresq) then

                  if(t(i,j,1).lt.t0) then

                    clri(i,j,1)=cclri*qr2b*ncr(i,j,1)*qi(i,j,1)
                    clir(i,j,1)=clri(i,j,1)

                  else

                    clri(i,j,1)=0.e0
                    clir(i,j,1)=0.e0

                  end if

                else

                  clri(i,j,1)=0.e0
                  clir(i,j,1)=0.e0

                end if

                if(qs(i,j,1).gt.thresq) then

                  a=2.e0*diaqr(i,j,1)*diaqs(i,j,1)

                  b=urq(i,j,1)-usq(i,j,1)

                  b=ncr(i,j,1)*ncs(i,j,1)*rbr(i,j,1)                    &
     &              *sqrt(b*b+.04e0*urq(i,j,1)*usq(i,j,1))

                  clrs(i,j,1)=b*cclrs*(5.e0*diaqr2+a+.5e0*diaqs2)*diaqr3
                  clsr(i,j,1)=b*cclsr*(5.e0*diaqs2+a+.5e0*diaqr2)*diaqs3

                else

                  clrs(i,j,1)=0.e0
                  clsr(i,j,1)=0.e0

                end if

                if(qg(i,j,1).gt.thresq) then

                  a=5.e0*diaqr2                                         &
     &              +2.e0*diaqr(i,j,1)*diaqg(i,j,1)+.5e0*diaqg2

                  b=urq(i,j,1)-ugq(i,j,1)

                  clrg(i,j,1)=a*cclrg*diaqr3*ncr(i,j,1)*ncg(i,j,1)      &
     &              *rbr(i,j,1)*sqrt(b*b+.04e0*urq(i,j,1)*ugq(i,j,1))

                else

                  clrg(i,j,1)=0.e0

                end if

                sink=clri(i,j,1)+clrs(i,j,1)+clrg(i,j,1)

                if(qr(i,j,1).lt.sink) then

                  a=qr(i,j,1)/sink

                  clri(i,j,1)=clri(i,j,1)*a
                  clrs(i,j,1)=clrs(i,j,1)*a
                  clrg(i,j,1)=clrg(i,j,1)*a

                end if

              else

                clri(i,j,1)=0.e0
                clrs(i,j,1)=0.e0
                clrg(i,j,1)=0.e0

                clir(i,j,1)=0.e0

                clsr(i,j,1)=0.e0

              end if

! -----

!! -----

! Fill in the array with 0 in the case the air temperature is lower than
! lowest super cooled point.

            else

              clcr(i,j,1)=0.e0
              clcs(i,j,1)=0.e0
              clcg(i,j,1)=0.e0

              clri(i,j,1)=0.e0
              clrs(i,j,1)=0.e0
              clrg(i,j,1)=0.e0

              clir(i,j,1)=0.e0

              clsr(i,j,1)=0.e0

              ecs(i,j,1)=0.e0

            end if

! -----

! Calculate the collection rate between the cloud ice and snow and
! graupel.

            if(qi(i,j,1).gt.thresq) then

              if(t(i,j,1).lt.t0) then

                if(qs(i,j,1).gt.thresq) then

                  clis(i,j,1)=cclis*qs2b*ncs(i,j,1)*qi(i,j,1)

                else

                  clis(i,j,1)=0.e0

                end if

                if(qg(i,j,1).gt.thresq) then

                  clig(i,j,1)=cclig*qg2b*ncg(i,j,1)*qi(i,j,1)

                else

                  clig(i,j,1)=0.e0

                end if

                sink=clir(i,j,1)+clis(i,j,1)+clig(i,j,1)

                if(qi(i,j,1).lt.sink) then

                  a=qi(i,j,1)/sink

                  clir(i,j,1)=clir(i,j,1)*a
                  clis(i,j,1)=clis(i,j,1)*a
                  clig(i,j,1)=clig(i,j,1)*a

                end if

              else

                clis(i,j,1)=0.e0
                clig(i,j,1)=0.e0

              end if

            else

              clis(i,j,1)=0.e0
              clig(i,j,1)=0.e0

            end if

! -----

! Calculate the collection rate between the snow and the graupel.

            if(qs(i,j,1).gt.thresq) then

              if(qg(i,j,1).gt.thresq) then

                if(t(i,j,1).lt.t0) then

                  a=5.e0*diaqs2                                         &
     &              +2.e0*diaqs(i,j,1)*diaqg(i,j,1)+.5e0*diaqg2

                  b=usq(i,j,1)-ugq(i,j,1)

                  clsg(i,j,1)=a*cclsg*diaqs3*ncs(i,j,1)*ncg(i,j,1)      &
     &              *rbr(i,j,1)*sqrt(b*b+.04e0*usq(i,j,1)*ugq(i,j,1))

                else

                  a=5.e0*diaqs2                                         &
     &              +2.e0*diaqs(i,j,1)*diaqg(i,j,1)+.5e0*diaqg2

                  b=usq(i,j,1)-ugq(i,j,1)

                  clsg(i,j,1)=a*cclsg3*diaqs3*ncs(i,j,1)*ncg(i,j,1)     &
     &              *rbr(i,j,1)*sqrt(b*b+.04e0*usq(i,j,1)*ugq(i,j,1))

                end if

              else

                clsg(i,j,1)=0.e0

              end if

              sink=clsr(i,j,1)+clsg(i,j,1)

              if(qs(i,j,1).lt.sink) then

                a=qs(i,j,1)/sink

                clsr(i,j,1)=clsr(i,j,1)*a
                clsg(i,j,1)=clsg(i,j,1)*a

              end if

            else

              clsg(i,j,1)=0.e0

            end if

! -----

          end do
          end do

!$omp end do

!!! -----

!!! Perform calculating in the case the option abs(cphopt) is greater
!!! than 2.

        else if(abs(cphopt).ge.3) then

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,cstk,r0rsq,qr2b,qs2b,qg2b)                         &
!$omp&   private(diaqr2,diaqs2,diaqg2,diaqr3,diaqs3,sink,a,b,c)

          do j=1,nj-1
          do i=1,ni-1

! Set the common used variables.

            r0rsq=sqrt(r0*rbr(i,j,1))

            if(qr(i,j,1).gt.thresq) then

              diaqr2=diaqr(i,j,1)*diaqr(i,j,1)
              diaqr3=diaqr(i,j,1)*diaqr2

              qr2b=r0rsq*exp(bur2*log(diaqr(i,j,1)))

            else

              diaqr2=0.e0
              diaqr3=0.e0

              qr2b=0.e0

            end if

            if(qs(i,j,1).gt.thresq) then

              diaqs2=diaqs(i,j,1)*diaqs(i,j,1)
              diaqs3=diaqs(i,j,1)*diaqs2

              qs2b=r0rsq*exp(bus2*log(diaqs(i,j,1)))

            else

              diaqs2=0.e0
              diaqs3=0.e0

              qs2b=0.e0

            end if

            if(qg(i,j,1).gt.thresq) then

              diaqg2=diaqg(i,j,1)*diaqg(i,j,1)

              qg2b=r0rsq*exp(bug2*log(diaqg(i,j,1)))

            else

              diaqg2=0.e0

              qg2b=0.e0

            end if

! -----

!! Calculate the collection rate between the cloud water and the rain
!! water, snow and graupel and between the rain water and the cloud ice,
!! the snow and the graupel.

            if(t(i,j,1).gt.tlow) then

! Calculate the collection rate between the cloud water and the rain
! water, snow and graupel.

              if(qc(i,j,1).gt.thresq) then

                cstk=rwdv9*diaqc(i,j,1)*diaqc(i,j,1)

                if(qr(i,j,1).gt.thresq) then

                  a=cstk*urq(i,j,1)/(diaqr(i,j,1)*mu(i,j,1))

                  b=a+.5e0

                  a=a*a/(b*b)

                  clcr(i,j,1)=cclcr*a*qr2b*ncr(i,j,1)*qc(i,j,1)

                else

                  clcr(i,j,1)=0.e0

                end if

                if(qs(i,j,1).gt.thresq) then

                  a=cstk*usq(i,j,1)/(diaqs(i,j,1)*mu(i,j,1))

                  b=a+.5e0

                  ecs(i,j,1)=a*a/(b*b)

                  clcs(i,j,1)=cclcs*qs2b*ecs(i,j,1)*ncs(i,j,1)*qc(i,j,1)

                else

                  clcs(i,j,1)=0.e0

                  ecs(i,j,1)=0.e0

                end if

                if(qg(i,j,1).gt.thresq) then

                  a=cstk*ugq(i,j,1)/(diaqg(i,j,1)*mu(i,j,1))

                  b=a+.5e0

                  a=a*a/(b*b)

                  clcg(i,j,1)=cclcg*a*qg2b*ncg(i,j,1)*qc(i,j,1)

                else

                  clcg(i,j,1)=0.e0

                end if

                sink=clcr(i,j,1)+clcs(i,j,1)+clcg(i,j,1)

                if(qc(i,j,1).lt.sink) then

                  a=qc(i,j,1)/sink

                  clcr(i,j,1)=clcr(i,j,1)*a
                  clcs(i,j,1)=clcs(i,j,1)*a
                  clcg(i,j,1)=clcg(i,j,1)*a

                end if

              else

                clcr(i,j,1)=0.e0
                clcs(i,j,1)=0.e0
                clcg(i,j,1)=0.e0

                ecs(i,j,1)=0.e0

              end if

! -----

! Calculate the collection rate between the rain water and the cloud
! ice, the snow and the graupel.

              if(qr(i,j,1).gt.thresq) then

                if(qi(i,j,1).gt.thresq) then

                  if(t(i,j,1).lt.t0) then

                    clri(i,j,1)=cclri*qr2b*ncr(i,j,1)*qi(i,j,1)
                    clir(i,j,1)=clri(i,j,1)

                    clrin(i,j,1)=clri(i,j,1)/mi(i,j,1)

                  else

                    clri(i,j,1)=0.e0
                    clir(i,j,1)=0.e0

                    clrin(i,j,1)=0.e0

                  end if

                else

                  clri(i,j,1)=0.e0
                  clir(i,j,1)=0.e0

                  clrin(i,j,1)=0.e0

                end if

                if(qs(i,j,1).gt.thresq) then

                  a=2.e0*diaqr(i,j,1)*diaqs(i,j,1)

                  b=ncr(i,j,1)*ncs(i,j,1)*rbr(i,j,1)

                  c=urq(i,j,1)-usq(i,j,1)

                  c=sqrt(c*c+.04e0*urq(i,j,1)*usq(i,j,1))

                  clrs(i,j,1)=b*c                                       &
     &              *cclrs*(5.e0*diaqr2+a+.5e0*diaqs2)*diaqr3

                  clsr(i,j,1)=b*c                                       &
     &              *cclsr*(5.e0*diaqs2+a+.5e0*diaqr2)*diaqs3

                  c=urn(i,j,1)-usn(i,j,1)

                  clrsn(i,j,1)=b*cclrsn*(diaqr2+.5e0*a+diaqs2)          &
     &              *sqrt(c*c+.04e0*urn(i,j,1)*usn(i,j,1))

                  clsrn(i,j,1)=clrsn(i,j,1)

                else

                  clrs(i,j,1)=0.e0
                  clsr(i,j,1)=0.e0

                  clrsn(i,j,1)=0.e0
                  clsrn(i,j,1)=0.e0

                end if

                if(qg(i,j,1).gt.thresq) then

                  a=5.e0*diaqr2                                         &
     &              +2.e0*diaqr(i,j,1)*diaqg(i,j,1)+.5e0*diaqg2

                  b=urq(i,j,1)-ugq(i,j,1)

                  clrg(i,j,1)=a*cclrg*diaqr3*ncr(i,j,1)*ncg(i,j,1)      &
     &              *rbr(i,j,1)*sqrt(b*b+.04e0*urq(i,j,1)*ugq(i,j,1))

                else

                  clrg(i,j,1)=0.e0

                end if

                sink=clri(i,j,1)+clrs(i,j,1)+clrg(i,j,1)

                if(qr(i,j,1).lt.sink) then

                  a=qr(i,j,1)/sink

                  clri(i,j,1)=clri(i,j,1)*a
                  clrs(i,j,1)=clrs(i,j,1)*a
                  clrg(i,j,1)=clrg(i,j,1)*a

                end if

                sink=clrin(i,j,1)+clrsn(i,j,1)

                if(ncr(i,j,1).lt.sink) then

                  a=ncr(i,j,1)/sink

                  clrin(i,j,1)=clrin(i,j,1)*a
                  clrsn(i,j,1)=clrsn(i,j,1)*a

                end if

              else

                clri(i,j,1)=0.e0
                clrs(i,j,1)=0.e0
                clrg(i,j,1)=0.e0

                clir(i,j,1)=0.e0

                clsr(i,j,1)=0.e0

                clrin(i,j,1)=0.e0
                clrsn(i,j,1)=0.e0
                clsrn(i,j,1)=0.e0

              end if

! -----

!! -----

! Fill in the array with 0 in the case the air temperature is lower than
! lowest super cooled point.

            else

              clcr(i,j,1)=0.e0
              clcs(i,j,1)=0.e0
              clcg(i,j,1)=0.e0

              clri(i,j,1)=0.e0
              clrs(i,j,1)=0.e0
              clrg(i,j,1)=0.e0

              clir(i,j,1)=0.e0

              clsr(i,j,1)=0.e0

              clrin(i,j,1)=0.e0
              clrsn(i,j,1)=0.e0
              clsrn(i,j,1)=0.e0

              ecs(i,j,1)=0.e0

            end if

! -----

! Calculate the collection rate between the cloud ice and snow and
! graupel.

            if(qi(i,j,1).gt.thresq) then

              if(t(i,j,1).lt.t0) then

                if(qs(i,j,1).gt.thresq) then

                  clis(i,j,1)=cclis*qs2b*ncs(i,j,1)*qi(i,j,1)

                else

                  clis(i,j,1)=0.e0

                end if

                if(qg(i,j,1).gt.thresq) then

                  clig(i,j,1)=cclig*qg2b*ncg(i,j,1)*qi(i,j,1)

                else

                  clig(i,j,1)=0.e0

                end if

                sink=clir(i,j,1)+clis(i,j,1)+clig(i,j,1)

                if(qi(i,j,1).lt.sink) then

                  a=qi(i,j,1)/sink

                  clir(i,j,1)=clir(i,j,1)*a
                  clis(i,j,1)=clis(i,j,1)*a
                  clig(i,j,1)=clig(i,j,1)*a

                end if

              else

                clis(i,j,1)=0.e0
                clig(i,j,1)=0.e0

              end if

            else

              clis(i,j,1)=0.e0
              clig(i,j,1)=0.e0

            end if

! -----

! Calculate the collection rate between the snow and the graupel.

            if(qs(i,j,1).gt.thresq) then

              if(qg(i,j,1).gt.thresq) then

                if(t(i,j,1).lt.t0) then

                  a=diaqs(i,j,1)*diaqg(i,j,1)

                  b=ncs(i,j,1)*ncg(i,j,1)*rbr(i,j,1)

                  c=usq(i,j,1)-ugq(i,j,1)

                  clsg(i,j,1)=b*cclsg*(5.e0*diaqs2+2.e0*a+.5e0*diaqg2)  &
     &              *diaqs3*sqrt(c*c+.04e0*usq(i,j,1)*ugq(i,j,1))

                  c=usn(i,j,1)-ugn(i,j,1)

                  clsgn(i,j,1)=b*cclsgn*(diaqs2+a+diaqg2)               &
     &              *sqrt(c*c+.04e0*usn(i,j,1)*ugn(i,j,1))

                else

                  a=diaqs(i,j,1)*diaqg(i,j,1)

                  b=esgiv*ncs(i,j,1)*ncg(i,j,1)*rbr(i,j,1)

                  c=usq(i,j,1)-ugq(i,j,1)

                  clsg(i,j,1)=b*cclsg*(5.e0*diaqs2+2.e0*a+.5e0*diaqg2)  &
     &              *diaqs3*sqrt(c*c+.04e0*usq(i,j,1)*ugq(i,j,1))

                  c=usn(i,j,1)-ugn(i,j,1)

                  clsgn(i,j,1)=b*cclsgn*(diaqs2+a+diaqg2)               &
     &              *sqrt(c*c+.04e0*usn(i,j,1)*ugn(i,j,1))

                end if

              else

                clsg(i,j,1)=0.e0
                clsgn(i,j,1)=0.e0

              end if

              sink=clsr(i,j,1)+clsg(i,j,1)

              if(qs(i,j,1).lt.sink) then

                a=qs(i,j,1)/sink

                clsr(i,j,1)=clsr(i,j,1)*a
                clsg(i,j,1)=clsg(i,j,1)*a

              end if

              sink=clsrn(i,j,1)+clsgn(i,j,1)

              if(ncs(i,j,1).lt.sink) then

                a=ncs(i,j,1)/sink

                clsrn(i,j,1)=clsrn(i,j,1)*a
                clsgn(i,j,1)=clsgn(i,j,1)*a

              end if

            else

              clsg(i,j,1)=0.e0
              clsgn(i,j,1)=0.e0

            end if

! -----

          end do
          end do

!$omp end do

        end if

!!! -----

!!!! -----

!!!! In the case nk > 1.

      else

!!! Perform calculating in the case the option abs(cphopt) is equal
!!! to 2.

        if(abs(cphopt).eq.2) then

          do k=1,nk-1

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,cstk,r0rsq,qr2b,qs2b,qg2b)                         &
!$omp&   private(diaqr2,diaqs2,diaqg2,diaqr3,diaqs3,sink,a,b,c)

            do j=1,nj-1
            do i=1,ni-1

! Set the common used variables.

              r0rsq=sqrt(r0*rbr(i,j,k))

              if(qr(i,j,k).gt.thresq) then

                diaqr2=diaqr(i,j,k)*diaqr(i,j,k)
                diaqr3=diaqr(i,j,k)*diaqr2

                qr2b=r0rsq*exp(bur2*log(diaqr(i,j,k)))

              else

                diaqr2=0.e0
                diaqr3=0.e0

                qr2b=0.e0

              end if

              if(qs(i,j,k).gt.thresq) then

                diaqs2=diaqs(i,j,k)*diaqs(i,j,k)
                diaqs3=diaqs(i,j,k)*diaqs2

                qs2b=r0rsq*exp(bus2*log(diaqs(i,j,k)))

              else

                diaqs2=0.e0
                diaqs3=0.e0

                qs2b=0.e0

              end if

              if(qg(i,j,k).gt.thresq) then

                diaqg2=diaqg(i,j,k)*diaqg(i,j,k)

                qg2b=r0rsq*exp(bug2*log(diaqg(i,j,k)))

              else

                diaqg2=0.e0

                qg2b=0.e0

              end if

! -----

!! Calculate the collection rate between the cloud water and the rain
!! water, snow and graupel and between the rain water and the cloud ice,
!! the snow and the graupel.

              if(t(i,j,k).gt.tlow) then

! Calculate the collection rate between the cloud water and the rain
! water, snow and graupel.

                if(qc(i,j,k).gt.thresq) then

                  cstk=rwdv9*diaqc(i,j,k)*diaqc(i,j,k)

                  if(qr(i,j,k).gt.thresq) then

                    a=cstk*urq(i,j,k)/(diaqr(i,j,k)*mu(i,j,k))

                    b=a+.5e0

                    a=a*a/(b*b)

                    clcr(i,j,k)=cclcr*a*qr2b*ncr(i,j,k)*qc(i,j,k)

                  else

                    clcr(i,j,k)=0.e0

                  end if

                  if(qs(i,j,k).gt.thresq) then

                    a=cstk*usq(i,j,k)/(diaqs(i,j,k)*mu(i,j,k))

                    b=a+.5e0

                    ecs(i,j,k)=a*a/(b*b)

                    clcs(i,j,k)=cclcs                                   &
     &                *qs2b*ecs(i,j,k)*ncs(i,j,k)*qc(i,j,k)

                  else

                    clcs(i,j,k)=0.e0

                    ecs(i,j,k)=0.e0

                  end if

                  if(qg(i,j,k).gt.thresq) then

                    a=cstk*ugq(i,j,k)/(diaqg(i,j,k)*mu(i,j,k))

                    b=a+.5e0

                    a=a*a/(b*b)

                    clcg(i,j,k)=cclcg*a*qg2b*ncg(i,j,k)*qc(i,j,k)

                  else

                    clcg(i,j,k)=0.e0

                  end if

                  sink=clcr(i,j,k)+clcs(i,j,k)+clcg(i,j,k)

                  if(qc(i,j,k).lt.sink) then

                    a=qc(i,j,k)/sink

                    clcr(i,j,k)=clcr(i,j,k)*a
                    clcs(i,j,k)=clcs(i,j,k)*a
                    clcg(i,j,k)=clcg(i,j,k)*a

                  end if

                else

                  clcr(i,j,k)=0.e0
                  clcs(i,j,k)=0.e0
                  clcg(i,j,k)=0.e0

                  ecs(i,j,k)=0.e0

                end if

! -----

! Calculate the collection rate between the rain water and the cloud
! ice, the snow and the graupel.

                if(qr(i,j,k).gt.thresq) then

                  if(qi(i,j,k).gt.thresq) then

                    if(t(i,j,k).lt.t0) then

                      clri(i,j,k)=cclri*qr2b*ncr(i,j,k)*qi(i,j,k)
                      clir(i,j,k)=clri(i,j,k)

                    else

                      clri(i,j,k)=0.e0
                      clir(i,j,k)=0.e0

                    end if

                  else

                    clri(i,j,k)=0.e0
                    clir(i,j,k)=0.e0

                  end if

                  if(qs(i,j,k).gt.thresq) then

                    a=2.e0*diaqr(i,j,k)*diaqs(i,j,k)

                    b=urq(i,j,k)-usq(i,j,k)

                    b=ncr(i,j,k)*ncs(i,j,k)*rbr(i,j,k)                  &
     &                *sqrt(b*b+.04e0*urq(i,j,k)*usq(i,j,k))

                    clrs(i,j,k)=b*cclrs                                 &
     &                *(5.e0*diaqr2+a+.5e0*diaqs2)*diaqr3

                    clsr(i,j,k)=b*cclsr                                 &
     &                *(5.e0*diaqs2+a+.5e0*diaqr2)*diaqs3

                  else

                    clrs(i,j,k)=0.e0
                    clsr(i,j,k)=0.e0

                  end if

                  if(qg(i,j,k).gt.thresq) then

                    a=5.e0*diaqr2                                       &
     &                +2.e0*diaqr(i,j,k)*diaqg(i,j,k)+.5e0*diaqg2

                    b=urq(i,j,k)-ugq(i,j,k)

                    clrg(i,j,k)=a*cclrg*diaqr3*ncr(i,j,k)*ncg(i,j,k)    &
     &                *rbr(i,j,k)*sqrt(b*b+.04e0*urq(i,j,k)*ugq(i,j,k))

                  else

                    clrg(i,j,k)=0.e0

                  end if

                  sink=clri(i,j,k)+clrs(i,j,k)+clrg(i,j,k)

                  if(qr(i,j,k).lt.sink) then

                    a=qr(i,j,k)/sink

                    clri(i,j,k)=clri(i,j,k)*a
                    clrs(i,j,k)=clrs(i,j,k)*a
                    clrg(i,j,k)=clrg(i,j,k)*a

                  end if

                else

                  clri(i,j,k)=0.e0
                  clrs(i,j,k)=0.e0
                  clrg(i,j,k)=0.e0

                  clir(i,j,k)=0.e0

                  clsr(i,j,k)=0.e0

                end if

! -----

!! -----

! Fill in the array with 0 in the case the air temperature is lower than
! lowest super cooled point.

              else

                clcr(i,j,k)=0.e0
                clcs(i,j,k)=0.e0
                clcg(i,j,k)=0.e0

                clri(i,j,k)=0.e0
                clrs(i,j,k)=0.e0
                clrg(i,j,k)=0.e0

                clir(i,j,k)=0.e0

                clsr(i,j,k)=0.e0

                ecs(i,j,k)=0.e0

              end if

! -----

! Calculate the collection rate between the cloud ice and snow and
! graupel.

              if(qi(i,j,k).gt.thresq) then

                if(t(i,j,k).lt.t0) then

                  if(qs(i,j,k).gt.thresq) then

                    clis(i,j,k)=cclis*qs2b*ncs(i,j,k)*qi(i,j,k)

                  else

                    clis(i,j,k)=0.e0

                  end if

                  if(qg(i,j,k).gt.thresq) then

                    clig(i,j,k)=cclig*qg2b*ncg(i,j,k)*qi(i,j,k)

                  else

                    clig(i,j,k)=0.e0

                  end if

                  sink=clir(i,j,k)+clis(i,j,k)+clig(i,j,k)

                  if(qi(i,j,k).lt.sink) then

                    a=qi(i,j,k)/sink

                    clir(i,j,k)=clir(i,j,k)*a
                    clis(i,j,k)=clis(i,j,k)*a
                    clig(i,j,k)=clig(i,j,k)*a

                  end if

                else

                  clis(i,j,k)=0.e0
                  clig(i,j,k)=0.e0

                end if

              else

                clis(i,j,k)=0.e0
                clig(i,j,k)=0.e0

              end if

! -----

! Calculate the collection rate between the snow and the graupel.

              if(qs(i,j,k).gt.thresq) then

                if(qg(i,j,k).gt.thresq) then

                  if(t(i,j,k).lt.t0) then

                    a=5.e0*diaqs2                                       &
     &                +2.e0*diaqs(i,j,k)*diaqg(i,j,k)+.5e0*diaqg2

                    b=usq(i,j,k)-ugq(i,j,k)

                    clsg(i,j,k)=a*cclsg*diaqs3*ncs(i,j,k)*ncg(i,j,k)    &
     &                *rbr(i,j,k)*sqrt(b*b+.04e0*usq(i,j,k)*ugq(i,j,k))

                  else

                    a=5.e0*diaqs2                                       &
     &                +2.e0*diaqs(i,j,k)*diaqg(i,j,k)+.5e0*diaqg2

                    b=usq(i,j,k)-ugq(i,j,k)

                    clsg(i,j,k)=a*cclsg3*diaqs3*ncs(i,j,k)*ncg(i,j,k)   &
     &                *rbr(i,j,k)*sqrt(b*b+.04e0*usq(i,j,k)*ugq(i,j,k))

                  end if

                else

                  clsg(i,j,k)=0.e0

                end if

                sink=clsr(i,j,k)+clsg(i,j,k)

                if(qs(i,j,k).lt.sink) then

                  a=qs(i,j,k)/sink

                  clsr(i,j,k)=clsr(i,j,k)*a
                  clsg(i,j,k)=clsg(i,j,k)*a

                end if

              else

                clsg(i,j,k)=0.e0

              end if

! -----

            end do
            end do

!$omp end do

          end do

!!! -----

!!! Perform calculating in the case the option abs(cphopt) is greater
!!! than 2.

        else if(abs(cphopt).ge.3) then

          do k=1,nk-1

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,cstk,r0rsq,qr2b,qs2b,qg2b)                         &
!$omp&   private(diaqr2,diaqs2,diaqg2,diaqr3,diaqs3,sink,a,b,c)

            do j=1,nj-1
            do i=1,ni-1

! Set the common used variables.

              r0rsq=sqrt(r0*rbr(i,j,k))

              if(qr(i,j,k).gt.thresq) then

                diaqr2=diaqr(i,j,k)*diaqr(i,j,k)
                diaqr3=diaqr(i,j,k)*diaqr2

                qr2b=r0rsq*exp(bur2*log(diaqr(i,j,k)))

              else

                diaqr2=0.e0
                diaqr3=0.e0

                qr2b=0.e0

              end if

              if(qs(i,j,k).gt.thresq) then

                diaqs2=diaqs(i,j,k)*diaqs(i,j,k)
                diaqs3=diaqs(i,j,k)*diaqs2

                qs2b=r0rsq*exp(bus2*log(diaqs(i,j,k)))

              else

                diaqs2=0.e0
                diaqs3=0.e0

                qs2b=0.e0

              end if

              if(qg(i,j,k).gt.thresq) then

                diaqg2=diaqg(i,j,k)*diaqg(i,j,k)

                qg2b=r0rsq*exp(bug2*log(diaqg(i,j,k)))

              else

                diaqg2=0.e0

                qg2b=0.e0

              end if

! -----

!! Calculate the collection rate between the cloud water and the rain
!! water, snow and graupel and between the rain water and the cloud ice,
!! the snow and the graupel.

              if(t(i,j,k).gt.tlow) then

! Calculate the collection rate between the cloud water and the rain
! water, snow and graupel.

                if(qc(i,j,k).gt.thresq) then

                  cstk=rwdv9*diaqc(i,j,k)*diaqc(i,j,k)

                  if(qr(i,j,k).gt.thresq) then

                    a=cstk*urq(i,j,k)/(diaqr(i,j,k)*mu(i,j,k))

                    b=a+.5e0

                    a=a*a/(b*b)

                    clcr(i,j,k)=cclcr*a*qr2b*ncr(i,j,k)*qc(i,j,k)

                  else

                    clcr(i,j,k)=0.e0

                  end if

                  if(qs(i,j,k).gt.thresq) then

                    a=cstk*usq(i,j,k)/(diaqs(i,j,k)*mu(i,j,k))

                    b=a+.5e0

                    ecs(i,j,k)=a*a/(b*b)

                    clcs(i,j,k)=cclcs                                   &
     &                *qs2b*ecs(i,j,k)*ncs(i,j,k)*qc(i,j,k)

                  else

                    clcs(i,j,k)=0.e0

                    ecs(i,j,k)=0.e0

                  end if

                  if(qg(i,j,k).gt.thresq) then

                    a=cstk*ugq(i,j,k)/(diaqg(i,j,k)*mu(i,j,k))

                    b=a+.5e0

                    a=a*a/(b*b)

                    clcg(i,j,k)=cclcg*a*qg2b*ncg(i,j,k)*qc(i,j,k)

                  else

                    clcg(i,j,k)=0.e0

                  end if

                  sink=clcr(i,j,k)+clcs(i,j,k)+clcg(i,j,k)

                  if(qc(i,j,k).lt.sink) then

                    a=qc(i,j,k)/sink

                    clcr(i,j,k)=clcr(i,j,k)*a
                    clcs(i,j,k)=clcs(i,j,k)*a
                    clcg(i,j,k)=clcg(i,j,k)*a

                  end if

                else

                  clcr(i,j,k)=0.e0
                  clcs(i,j,k)=0.e0
                  clcg(i,j,k)=0.e0

                  ecs(i,j,k)=0.e0

                end if

! -----

! Calculate the collection rate between the rain water and the cloud
! ice, the snow and the graupel.

                if(qr(i,j,k).gt.thresq) then

                  if(qi(i,j,k).gt.thresq) then

                    if(t(i,j,k).lt.t0) then

                      clri(i,j,k)=cclri*qr2b*ncr(i,j,k)*qi(i,j,k)
                      clir(i,j,k)=clri(i,j,k)

                      clrin(i,j,k)=clri(i,j,k)/mi(i,j,k)

                    else

                      clri(i,j,k)=0.e0
                      clir(i,j,k)=0.e0

                      clrin(i,j,k)=0.e0

                    end if

                  else

                    clri(i,j,k)=0.e0
                    clir(i,j,k)=0.e0

                    clrin(i,j,k)=0.e0

                  end if

                  if(qs(i,j,k).gt.thresq) then

                    a=2.e0*diaqr(i,j,k)*diaqs(i,j,k)

                    b=ncr(i,j,k)*ncs(i,j,k)*rbr(i,j,k)

                    c=urq(i,j,k)-usq(i,j,k)

                    c=sqrt(c*c+.04e0*urq(i,j,k)*usq(i,j,k))

                    clrs(i,j,k)=b*c                                     &
     &                *cclrs*(5.e0*diaqr2+a+.5e0*diaqs2)*diaqr3

                    clsr(i,j,k)=b*c                                     &
     &                *cclsr*(5.e0*diaqs2+a+.5e0*diaqr2)*diaqs3

                    c=urn(i,j,k)-usn(i,j,k)

                    clrsn(i,j,k)=b*cclrsn*(diaqr2+.5e0*a+diaqs2)        &
     &                *sqrt(c*c+.04e0*urn(i,j,k)*usn(i,j,k))

                    clsrn(i,j,k)=clrsn(i,j,k)

                  else

                    clrs(i,j,k)=0.e0
                    clsr(i,j,k)=0.e0

                    clrsn(i,j,k)=0.e0
                    clsrn(i,j,k)=0.e0

                  end if

                  if(qg(i,j,k).gt.thresq) then

                    a=5.e0*diaqr2                                       &
     &                +2.e0*diaqr(i,j,k)*diaqg(i,j,k)+.5e0*diaqg2

                    b=urq(i,j,k)-ugq(i,j,k)

                    clrg(i,j,k)=a*cclrg*diaqr3*ncr(i,j,k)*ncg(i,j,k)    &
     &                *rbr(i,j,k)*sqrt(b*b+.04e0*urq(i,j,k)*ugq(i,j,k))

                  else

                    clrg(i,j,k)=0.e0

                  end if

                  sink=clri(i,j,k)+clrs(i,j,k)+clrg(i,j,k)

                  if(qr(i,j,k).lt.sink) then

                    a=qr(i,j,k)/sink

                    clri(i,j,k)=clri(i,j,k)*a
                    clrs(i,j,k)=clrs(i,j,k)*a
                    clrg(i,j,k)=clrg(i,j,k)*a

                  end if

                  sink=clrin(i,j,k)+clrsn(i,j,k)

                  if(ncr(i,j,k).lt.sink) then

                    a=ncr(i,j,k)/sink

                    clrin(i,j,k)=clrin(i,j,k)*a
                    clrsn(i,j,k)=clrsn(i,j,k)*a

                  end if

                else

                  clri(i,j,k)=0.e0
                  clrs(i,j,k)=0.e0
                  clrg(i,j,k)=0.e0

                  clir(i,j,k)=0.e0

                  clsr(i,j,k)=0.e0

                  clrin(i,j,k)=0.e0
                  clrsn(i,j,k)=0.e0
                  clsrn(i,j,k)=0.e0

                end if

! -----

!! -----

! Fill in the array with 0 in the case the air temperature is lower than
! lowest super cooled point.

              else

                clcr(i,j,k)=0.e0
                clcs(i,j,k)=0.e0
                clcg(i,j,k)=0.e0

                clri(i,j,k)=0.e0
                clrs(i,j,k)=0.e0
                clrg(i,j,k)=0.e0

                clir(i,j,k)=0.e0

                clsr(i,j,k)=0.e0

                clrin(i,j,k)=0.e0
                clrsn(i,j,k)=0.e0
                clsrn(i,j,k)=0.e0

                ecs(i,j,k)=0.e0

              end if

! -----

! Calculate the collection rate between the cloud ice and snow and
! graupel.

              if(qi(i,j,k).gt.thresq) then

                if(t(i,j,k).lt.t0) then

                  if(qs(i,j,k).gt.thresq) then

                    clis(i,j,k)=cclis*qs2b*ncs(i,j,k)*qi(i,j,k)

                  else

                    clis(i,j,k)=0.e0

                  end if

                  if(qg(i,j,k).gt.thresq) then

                    clig(i,j,k)=cclig*qg2b*ncg(i,j,k)*qi(i,j,k)

                  else

                    clig(i,j,k)=0.e0

                  end if

                  sink=clir(i,j,k)+clis(i,j,k)+clig(i,j,k)

                  if(qi(i,j,k).lt.sink) then

                    a=qi(i,j,k)/sink

                    clir(i,j,k)=clir(i,j,k)*a
                    clis(i,j,k)=clis(i,j,k)*a
                    clig(i,j,k)=clig(i,j,k)*a

                  end if

                else

                  clis(i,j,k)=0.e0
                  clig(i,j,k)=0.e0

                end if

              else

                clis(i,j,k)=0.e0
                clig(i,j,k)=0.e0

              end if

! -----

! Calculate the collection rate between the snow and the graupel.

              if(qs(i,j,k).gt.thresq) then

                if(qg(i,j,k).gt.thresq) then

                  if(t(i,j,k).lt.t0) then

                    a=diaqs(i,j,k)*diaqg(i,j,k)

                    b=ncs(i,j,k)*ncg(i,j,k)*rbr(i,j,k)

                    c=usq(i,j,k)-ugq(i,j,k)

                    clsg(i,j,k)=b*cclsg*(5.e0*diaqs2+2.e0*a+.5e0*diaqg2)&
     &                *diaqs3*sqrt(c*c+.04e0*usq(i,j,k)*ugq(i,j,k))

                    c=usn(i,j,k)-ugn(i,j,k)

                    clsgn(i,j,k)=b*cclsgn*(diaqs2+a+diaqg2)             &
     &                *sqrt(c*c+.04e0*usn(i,j,k)*ugn(i,j,k))

                  else

                    a=diaqs(i,j,k)*diaqg(i,j,k)

                    b=esgiv*ncs(i,j,k)*ncg(i,j,k)*rbr(i,j,k)

                    c=usq(i,j,k)-ugq(i,j,k)

                    clsg(i,j,k)=b*cclsg*(5.e0*diaqs2+2.e0*a+.5e0*diaqg2)&
     &                *diaqs3*sqrt(c*c+.04e0*usq(i,j,k)*ugq(i,j,k))

                    c=usn(i,j,k)-ugn(i,j,k)

                    clsgn(i,j,k)=b*cclsgn*(diaqs2+a+diaqg2)             &
     &                *sqrt(c*c+.04e0*usn(i,j,k)*ugn(i,j,k))

                  end if

                else

                  clsg(i,j,k)=0.e0
                  clsgn(i,j,k)=0.e0

                end if

                sink=clsr(i,j,k)+clsg(i,j,k)

                if(qs(i,j,k).lt.sink) then

                  a=qs(i,j,k)/sink

                  clsr(i,j,k)=clsr(i,j,k)*a
                  clsg(i,j,k)=clsg(i,j,k)*a

                end if

                sink=clsrn(i,j,k)+clsgn(i,j,k)

                if(ncs(i,j,k).lt.sink) then

                  a=ncs(i,j,k)/sink

                  clsrn(i,j,k)=clsrn(i,j,k)*a
                  clsgn(i,j,k)=clsgn(i,j,k)*a

                end if

              else

                clsg(i,j,k)=0.e0
                clsgn(i,j,k)=0.e0

              end if

! -----

            end do
            end do

!$omp end do

          end do

        end if

!!! -----

      end if

!!!! -----

!$omp end parallel

!!!!! -----

      end subroutine s_collect

!-----7--------------------------------------------------------------7--

      end module m_collect
