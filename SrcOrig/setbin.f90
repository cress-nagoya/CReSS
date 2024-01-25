!***********************************************************************
      module m_setbin
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/03/06
!     Modification: 2006/08/08, 2006/09/30, 2007/01/20, 2007/01/31,
!                   2007/04/11, 2007/10/19, 2008/01/11, 2008/05/02,
!                   2008/08/25, 2009/01/30, 2009/02/27, 2009/03/12,
!                   2009/11/13, 2011/03/18, 2011/08/18, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the bin parameters.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_combin
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

      public :: setbin, s_setbin

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface setbin

        module procedure s_setbin

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

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
      subroutine s_setbin(fpbbinw,fpsbinw,nqw)
!***********************************************************************

! Input variable

      integer, intent(in) :: fpbbinw
                       ! Formal parameter of unique index of bbinw

      integer, intent(in) :: fpsbinw
                       ! Formal parameter of unique index of sbinw

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

! Internal shared variables

      real bbinw       ! Base of exponential function
                       ! in mass components

      real sbinw       ! Exponent of exponential function
                       ! in mass components

      real clmb20      ! 125.7 x lmb20

      real dc0iv2      ! 2.0 / dc0

      real a           ! Parameter of bin resolution for mass
      real c           ! Parameter of bin resolution for radius

! Internal private variables

      integer n        ! Array index in bin categories

      integer nri      ! Array index in water bin radius
      integer nrj      ! Array index in water bin radius

      integer didx     ! Referenced index in rewbw
      integer cidx     ! Referenced index in rewbw

      real rcbw        ! Ratio of water bin radius

      real dd          ! Weighting coefficient for
                       ! interpolating refference coalescence efficiency

      real dc          ! Weighting coefficient for
                       ! interpolating refference coalescence efficiency

      real dd1m        ! 1.0 - dd
      real dc1m        ! 1.0 - dc

      real rbwij       ! rbw(nri) + rbw(nrj)

      real aa          ! Temporary array

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getrname(fpbbinw,bbinw)
      call getrname(fpsbinw,sbinw)

! -----

!!! Set the bin parameters.

! Set the common used variables.

      clmb20=125.7e0*lmb20

      dc0iv2=2.e0/dc0

      a=exp((1.e0/sbinw)*log(bbinw))
      c=exp((oned3/sbinw)*log(bbinw))

! -----

! Initialize the processed variables.

      brw(1)=rbwmin/sqrt(c)
      bmw(1,1)=.001e0*fourd3*cc*rhow*rbwmin*rbwmin*rbwmin/sqrt(a)

! -----

!! Calculate the bin parameters.

!$omp parallel default(shared)

! Set the common used variables and calculate the radius and mass at
! water bin boundary.

!$omp single private(n)

      do n=2,nqw+1

        brw(n)=c*brw(n-1)
        bmw(n,1)=a*bmw(n-1,1)

      end do

!$omp end single

! -----

! Calculate the other parameters for water bin boundary.

!$omp do schedule(runtime) private(n)

      do n=1,nqw+1

        rbrw(n,1)=2.e0*brw(n)

        rbrw(n,2)=brw(n)*brw(n)
        rbrw(n,3)=10.e0*brw(n)*rbrw(n,2)
        rbrw(n,4)=5.e0*rbrw(n,2)*rbrw(n,2)
        rbrw(n,5)=.1e0*brw(n)*rbrw(n,4)

        rbrw(n,6)=1.e0/brw(n)

        rbrw(n,7)=dc0iv2*rbrw(n,6)
        rbrw(n,8)=brw(n)/(brw(n)+2.16e-5)

      end do

!$omp end do

!$omp do schedule(runtime) private(n,aa)

      do n=1,nqw

        rbmw(n,1)=.5e0*(bmw(n+1,1)+bmw(n,1))
        rbmw(n,2)=oned3*(4.e0*rbmw(n,1)*rbmw(n,1)-bmw(n+1,1)*bmw(n,1))

        dbmw(n)=bmw(n+1,1)-bmw(n,1)

        aa=.01e0*dbmw(n)

        bmw(n+1,2)=bmw(n+1,1)-aa
        bmw(n,3)=bmw(n,1)+aa

      end do

!$omp end do

! -----

! Calculate the standard radius between adjacent water bins.

!$omp do schedule(runtime) private(n)

      do n=1,nqw
        rbw(n)=sqrt(brw(n)*brw(n+1))
      end do

!$omp end do

! -----

! Calculate the related parameters of standard radius between adjacent
! water bins.

!$omp do schedule(runtime) private(n)

      do n=1,nqw

        rrbw(n,1)=rbw(n)*rbw(n)
        rrbw(n,2)=rbw(n)*rrbw(n,1)

        rrbw(n,3)=1.e0/rbw(n)

        rrbw(n,4)=rbw(n)*(rbw(n)+clmb20)
        rrbw(n,5)=(1.e0+clmb20*rrbw(n,3))*rrbw(n,3)

      end do

!$omp end do

! -----

! Calculate the coalescence efficiency between optional water bins.

!$omp do schedule(runtime)                                              &
!$omp&   private(nri,nrj,n,didx,cidx,rcbw,dd,dc,dd1m,dc1m,rbwij)

      do nrj=1,nqw
      do nri=1,nqw

        n=max(nri,nrj)

        rcbw=rbw(min(nri,nrj))/rbw(n)

        if(rbw(n).lt.rrdbw(2)) then
          didx=2
        else if(rbw(n).ge.rrdbw(2).and.rbw(n).lt.rrdbw(3)) then
          didx=3
        else if(rbw(n).ge.rrdbw(3).and.rbw(n).lt.rrdbw(4)) then
          didx=4
        else if(rbw(n).ge.rrdbw(4).and.rbw(n).lt.rrdbw(5)) then
          didx=5
        else if(rbw(n).ge.rrdbw(5).and.rbw(n).lt.rrdbw(6)) then
          didx=6
        else if(rbw(n).ge.rrdbw(6).and.rbw(n).lt.rrdbw(7)) then
          didx=7
        else if(rbw(n).ge.rrdbw(7).and.rbw(n).lt.rrdbw(8)) then
          didx=8
        else if(rbw(n).ge.rrdbw(8).and.rbw(n).lt.rrdbw(9)) then
          didx=9
        else if(rbw(n).ge.rrdbw(9).and.rbw(n).lt.rrdbw(10)) then
          didx=10
        else if(rbw(n).ge.rrdbw(10).and.rbw(n).lt.rrdbw(11)) then
          didx=11
        else if(rbw(n).ge.rrdbw(11).and.rbw(n).lt.rrdbw(12)) then
          didx=12
        else if(rbw(n).ge.rrdbw(12).and.rbw(n).lt.rrdbw(13)) then
          didx=13
        else if(rbw(n).ge.rrdbw(13).and.rbw(n).lt.rrdbw(14)) then
          didx=14
        else
          didx=15
        end if

        if(rcbw.lt.rrcbw(2)) then
          cidx=2
        else if(rcbw.ge.rrcbw(2).and.rcbw.lt.rrcbw(3)) then
          cidx=3
        else if(rcbw.ge.rrcbw(3).and.rcbw.lt.rrcbw(4)) then
          cidx=4
        else if(rcbw.ge.rrcbw(4).and.rcbw.lt.rrcbw(5)) then
          cidx=5
        else if(rcbw.ge.rrcbw(5).and.rcbw.lt.rrcbw(6)) then
          cidx=6
        else if(rcbw.ge.rrcbw(6).and.rcbw.lt.rrcbw(7)) then
          cidx=7
        else if(rcbw.ge.rrcbw(7).and.rcbw.lt.rrcbw(8)) then
          cidx=8
        else if(rcbw.ge.rrcbw(8).and.rcbw.lt.rrcbw(9)) then
          cidx=9
        else if(rcbw.ge.rrcbw(9).and.rcbw.lt.rrcbw(10)) then
          cidx=10
        else if(rcbw.ge.rrcbw(10).and.rcbw.lt.rrcbw(11)) then
          cidx=11
        else if(rcbw.ge.rrcbw(11).and.rcbw.lt.rrcbw(12)) then
          cidx=12
        else if(rcbw.ge.rrcbw(12).and.rcbw.lt.rrcbw(13)) then
          cidx=13
        else if(rcbw.ge.rrcbw(13).and.rcbw.lt.rrcbw(14)) then
          cidx=14
        else if(rcbw.ge.rrcbw(14).and.rcbw.lt.rrcbw(15)) then
          cidx=15
        else if(rcbw.ge.rrcbw(15).and.rcbw.lt.rrcbw(16)) then
          cidx=16
        else if(rcbw.ge.rrcbw(16).and.rcbw.lt.rrcbw(17)) then
          cidx=17
        else if(rcbw.ge.rrcbw(17).and.rcbw.lt.rrcbw(18)) then
          cidx=18
        else if(rcbw.ge.rrcbw(18).and.rcbw.lt.rrcbw(19)) then
          cidx=19
        else if(rcbw.ge.rrcbw(19).and.rcbw.lt.rrcbw(20)) then
          cidx=20
        else
          cidx=21
        end if

        dd=min(1.e0,                                                    &
     &    max(0.e0,(rbw(n)-rrdbw(didx-1))/(rrdbw(didx)-rrdbw(didx-1))))

        dc=min(1.e0,                                                    &
     &    max(0.e0,(rcbw-rrcbw(cidx-1))/(rrcbw(cidx)-rrcbw(cidx-1))))

        dd1m=1.e0-dd
        dc1m=1.e0-dc

        ewbw(nri,nrj)                                                   &
     &    =dd*dc*rewbw(didx,cidx)+dd*dc1m*rewbw(didx,cidx-1)            &
     &    +dd1m*dc*rewbw(didx-1,cidx)+dd1m*dc1m*rewbw(didx-1,cidx-1)

        rbwij=rbw(nri)+rbw(nrj)

        ewbw(nri,nrj)=cc*rbwij*rbwij*ewbw(nri,nrj)

      end do
      end do

!$omp end do

! -----

!$omp end parallel

!! -----

!!! -----

      end subroutine s_setbin

!-----7--------------------------------------------------------------7--

      end module m_setbin
