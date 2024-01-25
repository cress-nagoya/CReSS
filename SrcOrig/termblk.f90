!***********************************************************************
      module m_termblk
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2000/10/27, 2000/11/17, 2001/04/15,
!                   2001/05/29, 2001/06/29, 2001/12/11, 2002/01/15,
!                   2002/04/02, 2002/12/02, 2003/03/28, 2003/04/30,
!                   2003/05/19, 2003/10/31, 2003/11/28, 2003/12/12,
!                   2004/03/05, 2004/03/22, 2004/04/01, 2004/04/15,
!                   2004/05/07, 2004/05/31, 2004/06/10, 2004/08/01,
!                   2004/08/20, 2004/09/01, 2004/09/10, 2004/09/25,
!                   2004/10/12, 2004/12/17, 2005/01/31, 2006/01/10,
!                   2006/04/03, 2006/09/30, 2007/10/19, 2007/11/26,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2011/01/14,
!                   2011/03/18, 2011/06/01, 2011/09/22, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set and calculate the terminal velocity.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comphy
      use m_getiname
      use m_getrname
      use m_temparam

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: termblk, s_termblk

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface termblk

        module procedure s_termblk

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
      subroutine s_termblk(fpcphopt,fphaiopt,fpthresq,ni,nj,nk,         &
     &                  rbv,qc,qr,qi,qs,qg,qh,ncc,ncr,nci,ncs,ncg,nch,  &
     &                  ucq,urq,uiq,usq,ugq,uhq,ucn,urn,uin,usn,ugn,uhn)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fphaiopt
                       ! Formal parameter of unique index of haiopt

      integer, intent(in) :: fpthresq
                       ! Formal parameter of unique index of thresq

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: rbv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of base state density

      real, intent(in) :: qc(0:ni+1,0:nj+1,1:nk)
                       ! Cloud water mixnig ratio

      real, intent(in) :: qr(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixnig ratio

      real, intent(in) :: qi(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixnig ratio

      real, intent(in) :: qs(0:ni+1,0:nj+1,1:nk)
                       ! Snow mixnig ratio

      real, intent(in) :: qg(0:ni+1,0:nj+1,1:nk)
                       ! Graupel mixnig ratio

      real, intent(in) :: qh(0:ni+1,0:nj+1,1:nk)
                       ! Hail mixnig ratio

      real, intent(in) :: ncc(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud water

      real, intent(in) :: ncr(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of rain water

      real, intent(in) :: nci(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud ice

      real, intent(in) :: ncs(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of snow

      real, intent(in) :: ncg(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of graupel

      real, intent(in) :: nch(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of hail

! Output variables

      real, intent(out) :: ucq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of cloud water

      real, intent(out) :: urq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of rain water

      real, intent(out) :: uiq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of cloud ice

      real, intent(out) :: usq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of snow

      real, intent(out) :: ugq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of graupel

      real, intent(out) :: uhq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of hail

      real, intent(out) :: ucn(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of cloud water concentrations

      real, intent(out) :: urn(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of rain water concentrations

      real, intent(out) :: uin(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of cloud ice concentrations

      real, intent(out) :: usn(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of snow concentrations

      real, intent(out) :: ugn(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of graupel concentrations

      real, intent(out) :: uhn(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of hail concentrations

! Internal shared variables

      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes

      real thresq      ! Minimum threshold value of mixing ratio

      real buc3        ! buc / 3.0
      real bur3        ! bur / 3.0
      real bui3        ! bui / 3.0
      real bus3        ! bus / 3.0
      real bug3        ! bug / 3.0
      real buh3        ! buh / 3.0

      real cucq        ! Coefficient of terminal velocity
                       ! of cloud water

      real curq        ! Coefficient of terminal velocity
                       ! of rain water

      real cusq        ! Coefficient of terminal velocity
                       ! of snow

      real cugq        ! Coefficient of terminal velocity
                       ! of graupel

      real cuhq        ! Coefficient of terminal velocity
                       ! of hail

      real cucn        ! Coefficient of terminal velocity
                       ! of cloud water concentrations

      real curn        ! Coefficient of terminal velocity
                       ! of rain water concentrations

      real cusn        ! Coefficient of terminal velocity
                       ! of snow concentrations

      real cugn        ! Coefficient of terminal velocity
                       ! of graupel concentrations

      real cuhn        ! Coefficient of terminal velocity
                       ! of hail concentrations

      real ccrw6       ! 6.0 / (cc x rhow)
      real ccri6       ! 6.0 / (cc x rhoi)

      real cdiaqc      ! Coefficient of mean diameter
                       ! of cloud water

      real cdiaqr      ! Coefficient of mean diameter
                       ! of rain water

      real cdiaqs      ! Coefficient of mean diameter
                       ! of snow

      real cdiaqg      ! Coefficient of mean diameter
                       ! of graupel

      real cdiaqh      ! Coefficient of mean diameter
                       ! of hail

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real lnr0r       ! ln(r0 x rbv)

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)
      call getrname(fpthresq,thresq)

! -----

! Set the common used variables.

      buc3=oned3*buc
      bur3=oned3*bur
      bui3=oned3*bui
      bus3=oned3*bus
      bug3=oned3*bug
      buh3=oned3*buh

      cucq=oned6*auc*gf4buc
      curq=oned6*aur*gf4bur
      cusq=oned6*aus*gf4bus
      cugq=oned6*aug*gf4bug
      cuhq=oned6*auh*gf4buh

      cucn=6.e0*gf1buc/gf4buc
      curn=6.e0*gf1bur/gf4bur
      cusn=6.e0*gf1bus/gf4bus
      cugn=6.e0*gf1bug/gf4bug
      cuhn=6.e0*gf1buh/gf4buh

      cdiaqc=1.e0/(cc*rhow)
      cdiaqr=1.e0/(cc*rhow)
      cdiaqs=1.e0/(cc*rhos)
      cdiaqg=1.e0/(cc*rhog)
      cdiaqh=1.e0/(cc*rhoh)

      ccrw6=6.e0/(cc*rhow)
      ccri6=6.e0/(cc*rhoi)

! -----

!! Set and calculate the terminal velocity.

!$omp parallel default(shared) private(k)

! Set the terminal velocity of the cloud water and cloud ice.

      if(flqcqi_opt.eq.1) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1

            if(qc(i,j,k).gt.thresq) then
              ucq(i,j,k)=ucqcst
              ucn(i,j,k)=ucncst

            else
              ucq(i,j,k)=0.e0
              ucn(i,j,k)=0.e0

            end if

            if(qi(i,j,k).gt.thresq) then
              uiq(i,j,k)=uiqcst
              uin(i,j,k)=uincst

            else
              uiq(i,j,k)=0.e0
              uin(i,j,k)=0.e0

            end if

          end do
          end do

!$omp end do

        end do

      else if(flqcqi_opt.eq.2) then

        if(abs(cphopt).le.3) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j,lnr0r)

            do j=1,nj-1
            do i=1,ni-1

              lnr0r=log(r0*rbv(i,j,k))

              if(qc(i,j,k).gt.thresq) then

                ucq(i,j,k)=auc*exp(guc*lnr0r)                           &
     &            *exp(buc3*log(ccrw6*qc(i,j,k)/ncc(i,j,k)))

                ucn(i,j,k)=ucq(i,j,k)

              else

                ucq(i,j,k)=0.e0
                ucn(i,j,k)=0.e0

              end if

              if(qi(i,j,k).gt.thresq) then

                uiq(i,j,k)=aui*exp(gui*lnr0r)                           &
     &            *exp(bui3*log(ccri6*qi(i,j,k)/nci(i,j,k)))

                uin(i,j,k)=uiq(i,j,k)

              else

                uiq(i,j,k)=0.e0
                uin(i,j,k)=0.e0

              end if

            end do
            end do

!$omp end do

          end do

        else

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j,lnr0r)

            do j=1,nj-1
            do i=1,ni-1

              lnr0r=log(r0*rbv(i,j,k))

              if(qc(i,j,k).gt.thresq) then

                ucq(i,j,k)=cucq*exp(guc*lnr0r)                          &
     &            *exp(buc3*log(cdiaqc*qc(i,j,k)/ncc(i,j,k)))

                ucn(i,j,k)=cucn*ucq(i,j,k)

              else

                ucq(i,j,k)=0.e0
                ucn(i,j,k)=0.e0

              end if

              if(qi(i,j,k).gt.thresq) then

                uiq(i,j,k)=aui*exp(gui*lnr0r)                           &
     &            *exp(bui3*log(ccri6*qi(i,j,k)/nci(i,j,k)))

                uin(i,j,k)=uiq(i,j,k)

              else

                uiq(i,j,k)=0.e0
                uin(i,j,k)=0.e0

              end if

            end do
            end do

!$omp end do

          end do

        end if

      end if

! -----

! Calculate the terminal velocity of the rain water, the snow and the
! graupel.

      if(haiopt.eq.0) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j,lnr0r)

          do j=1,nj-1
          do i=1,ni-1

            lnr0r=log(r0*rbv(i,j,k))

            if(qr(i,j,k).gt.thresq) then

              urq(i,j,k)=curq*exp(gur*lnr0r)                            &
     &          *exp(bur3*log(cdiaqr*qr(i,j,k)/ncr(i,j,k)))

              urn(i,j,k)=curn*urq(i,j,k)

            else

              urq(i,j,k)=0.e0
              urn(i,j,k)=0.e0

            end if

            if(qs(i,j,k).gt.thresq) then

              usq(i,j,k)=cusq*exp(gus*lnr0r)                            &
     &          *exp(bus3*log(cdiaqs*qs(i,j,k)/ncs(i,j,k)))

              usn(i,j,k)=cusn*usq(i,j,k)

            else

              usq(i,j,k)=0.e0
              usn(i,j,k)=0.e0

            end if

            if(qg(i,j,k).gt.thresq) then

              ugq(i,j,k)=cugq*exp(gug*lnr0r)                            &
     &          *exp(bug3*log(cdiaqg*qg(i,j,k)/ncg(i,j,k)))

              ugn(i,j,k)=cugn*ugq(i,j,k)

            else

              ugq(i,j,k)=0.e0
              ugn(i,j,k)=0.e0

            end if

          end do
          end do

!$omp end do

        end do

! -----

! Calculate the terminal velocity of the rain water, the snow, the
! graupel and the hail.

      else

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j,lnr0r)

          do j=1,nj-1
          do i=1,ni-1

            lnr0r=log(r0*rbv(i,j,k))

            if(qr(i,j,k).gt.thresq) then

              urq(i,j,k)=curq*exp(gur*lnr0r)                            &
     &          *exp(bur3*log(cdiaqr*qr(i,j,k)/ncr(i,j,k)))

              urn(i,j,k)=curn*urq(i,j,k)

            else

              urq(i,j,k)=0.e0
              urn(i,j,k)=0.e0

            end if

            if(qs(i,j,k).gt.thresq) then

              usq(i,j,k)=cusq*exp(gus*lnr0r)                            &
     &          *exp(bus3*log(cdiaqs*qs(i,j,k)/ncs(i,j,k)))

              usn(i,j,k)=cusn*usq(i,j,k)

            else

              usq(i,j,k)=0.e0
              usn(i,j,k)=0.e0

            end if

            if(qg(i,j,k).gt.thresq) then

              ugq(i,j,k)=cugq*exp(gug*lnr0r)                            &
     &          *exp(bug3*log(cdiaqg*qg(i,j,k)/ncg(i,j,k)))

              ugn(i,j,k)=cugn*ugq(i,j,k)

            else

              ugq(i,j,k)=0.e0
              ugn(i,j,k)=0.e0

            end if

            if(qh(i,j,k).gt.thresq) then

              uhq(i,j,k)=cuhq*exp(guh*lnr0r)                            &
     &          *exp(buh3*log(cdiaqh*qh(i,j,k)/nch(i,j,k)))

              uhn(i,j,k)=cuhn*uhq(i,j,k)

            else

              uhq(i,j,k)=0.e0
              uhn(i,j,k)=0.e0

            end if

          end do
          end do

!$omp end do

        end do

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_termblk

!-----7--------------------------------------------------------------7--

      end module m_termblk
