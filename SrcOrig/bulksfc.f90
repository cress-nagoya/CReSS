!***********************************************************************
      module m_bulksfc
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/10/16
!     Modification: 2001/11/14, 2001/12/03, 2001/12/10, 2002/04/02,
!                   2002/12/02, 2003/04/30, 2003/05/19, 2003/07/15,
!                   2003/10/31, 2003/12/12, 2004/02/01, 2004/03/05,
!                   2004/04/01, 2004/04/15, 2004/05/07, 2004/08/01,
!                   2004/08/20, 2004/09/10, 2006/05/12, 2007/09/04,
!                   2007/10/19, 2008/05/02, 2008/07/01, 2008/08/25,
!                   2009/02/27, 2009/06/16, 2009/11/05, 2011/06/01,
!                   2011/09/22, 2011/11/10, 2011/12/17, 2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the bulk coefficients of surface flux.

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

      public :: bulksfc, s_bulksfc

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface bulksfc

        module procedure s_bulksfc

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic acos
      intrinsic atan
      intrinsic cos
      intrinsic exp
      intrinsic log
      intrinsic max
      intrinsic min
      intrinsic sqrt
      intrinsic tan

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_bulksfc(ni,nj,za,land,kai,z0m,z0h,rch,cm,ch)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

      real, intent(in) :: za(0:ni+1,0:nj+1)
                       ! z physical coordinates at lowest plane

      real, intent(in) :: kai(0:ni+1,0:nj+1)
                       ! Sea ice distribution

      real, intent(in) :: z0m(0:ni+1,0:nj+1)
                       ! Roughness length for velocity

      real, intent(in) :: z0h(0:ni+1,0:nj+1)
                       ! Roughness length for scalar

      real, intent(in) :: rch(0:ni+1,0:nj+1)
                       ! Richardson number

! Output variables

      real, intent(out) :: cm(0:ni+1,0:nj+1)
                       ! Bulk coefficient for velocity

      real, intent(out) :: ch(0:ni+1,0:nj+1)
                       ! Bulk coefficient for scalar

! Internal shared variables

      real cc05        ! 0.5 x cc

      real prgiv       ! 1.0 / prnumg
      real prwiv       ! 1.0 / prnumw

      real kp2         ! 2.0 x kappa
      real wkp2        ! 2.0 x wkappa

!ORIG real rkp         ! wkappa / kappa

      real kprg        ! kappa / prnumg
      real wkprw       ! wkappa / prnumw

      real kprg3       ! 10.0 x kappa / (3.0 x prnumg)
      real wkprw3      ! 10.0 x wkappa / (3.0 x prnumw)

      real icz0mv      ! 1.0 / icz0m
      real icz0hv      ! 1.0 / icz0h

      real rmg3v       ! 1.0 / (3.0 x rmg)
      real rms3v       ! 1.0 / (3.0 x rms)

      real cqg1        ! Constant of bulk coefficients for soil surface
      real cqg2        ! Constant of bulk coefficients for soil surface
      real cpg1        ! Constant of bulk coefficients for soil surface
      real cpg2        ! Constant of bulk coefficients for soil surface

      real cqs1        ! Constant of bulk coefficients for sea surface
      real cqs2        ! Constant of bulk coefficients for sea surface
      real cps1        ! Constant of bulk coefficients for sea surface
      real cps2        ! Constant of bulk coefficients for sea surface

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      real dz0m        ! za - z0m
      real dz0h        ! za - z0h

      real cmice       ! Bulk coefficient for velocity on ice surface
      real chice       ! Bulk coefficient for scalar on ice surface

      real a           ! Temporary variable
      real b           ! Temporary variable
      real c           ! Temporary variable
      real d           ! Temporary variable
      real e           ! Temporary variable
      real f           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Set the common used variables.

      cc05=.5e0*cc

      prgiv=1.e0/prnumg
      prwiv=1.e0/prnumw

      kp2=2.e0*kappa
      wkp2=2.e0*wkappa

!ORIG rkp=wkappa/kappa

      kprg=kappa/prnumg
      wkprw=wkappa/prnumw

      kprg3=tend3*kappa/prnumg
      wkprw3=tend3*wkappa/prnumw

      icz0mv=1.e0/icz0m
      icz0hv=1.e0/icz0h

      rmg3v=oned3/rmg
      rms3v=oned3/rms

      cqg1=oned9/(rmg*rmg)
      cqg2=rhg/(3.e0*rmg)
      cpg1=-oned27/(rmg*rmg*rmg)
      cpg2=(3.e0-rhg/rmg)/(6.e0*rmg)

      cqs1=oned9/(rms*rms)
      cqs2=rhs/(3.e0*rms)
      cps1=-oned27/(rms*rms*rms)
      cps2=(3.e0-rhs/rms)/(6.e0*rms)

! -----

!! Calculate the bulk coefficients of surface flux.

!$omp parallel default(shared)

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,dz0m,dz0h,cmice,chice,a,b,c,d,e,f)

      do j=1,nj-1
      do i=1,ni-1

! Set the common used variables.

        dz0m=za(i,j)-z0m(i,j)
        dz0h=za(i,j)-z0h(i,j)

        a=max(za(i,j)/z0m(i,j),1.01e0)
        b=max(za(i,j)/z0h(i,j),1.01e0)

! -----

! For the unstable case.

        if(rch(i,j).lt.0.e0) then

          if(land(i,j).lt.3) then

            a=log(a)
            b=log(b)

            c=rch(i,j)*prwiv
            c=c*c

            d=cqs1+cqs2*c
            e=cps1+cps2*c

            f=e*e-d*d*d

            c=(za(i,j)*dz0h*a*a)/(dz0m*dz0m*b)

            if(f.gt.0.e0) then

              f=exp(oned3*log(sqrt(f)+abs(e)))

              f=rms*c*(rms3v-(f+d/f))

            else

              f=sqrt(d)

              e=max(min(e/(d*f),1.e0),-1.e0)

              f=rms*c*(rms3v-2.e0*f*cos(oned3*acos(e)))

            end if

            f=sqrt(1.e0-min(f,1.e0))

            e=sqrt(f)

            c=2.e0*log(.5e0*(1.e0+e))                                   &
     &        +log(.5e0*(1.e0+f))-2.e0*atan(e)+cc05

            d=2.e0*log(.5e0*(1.e0+f))

            if(c.lt..5e0*a) then
              cm(i,j)=wkappa/(a-c)
            else
              cm(i,j)=wkp2/a
            end if

            if(d.lt..7e0*b) then
              ch(i,j)=wkprw/(b-d)
            else
              ch(i,j)=wkprw3/b
            end if

          else

            a=log(a)
            b=log(b)

            c=rch(i,j)*prgiv
            c=c*c

            d=cqg1+cqg2*c
            e=cpg1+cpg2*c

            f=e*e-d*d*d

            c=(za(i,j)*dz0h*a*a)/(dz0m*dz0m*b)

            if(f.gt.0.e0) then

              f=exp(oned3*log(sqrt(f)+abs(e)))

              f=rmg*c*(rmg3v-(f+d/f))

            else

              f=sqrt(d)

              e=max(min(e/(d*f),1.e0),-1.e0)

              f=rmg*c*(rmg3v-2.e0*f*cos(oned3*acos(e)))

            end if

            f=min(f,1.e0)

            e=sqrt(sqrt(1.e0-f))

            c=2.e0*log(.5e0*(1.e0+e))                                   &
     &        +log(.5e0*(1.e0+e*e))-2.e0*atan(e)+cc05

            f=sqrt(1.e0-.6e0*f)

            d=2.e0*log(.5e0*(1.e0+f))

            if(c.lt..5e0*a) then
              cm(i,j)=kappa/(a-c)
            else
              cm(i,j)=kp2/a
            end if

            if(d.lt..7e0*b) then
              ch(i,j)=kprg/(b-d)
            else
              ch(i,j)=kprg3/b
            end if

!ORIG       cm(i,j)=rkp*cm(i,j)
!ORIG       ch(i,j)=rkp*ch(i,j)

          end if

! -----

! For the stable case.

        else

          if(land(i,j).lt.3) then

            a=wkappa/log(a)
            b=wkappa/log(b)

            c=sqrt(1.e0+5.e0*rch(i,j))
            d=sqrt(1.e0+10.e0*rch(i,j)*c)

            cm(i,j)=a/d
            ch(i,j)=b*d/(prnumw*(1.e0+15.e0*rch(i,j)*c))

          else

            a=kappa/log(a)
            b=kappa/log(b)

            c=sqrt(1.e0+5.e0*rch(i,j))
            d=sqrt(1.e0+10.e0*rch(i,j)*c)

            cm(i,j)=a/d
            ch(i,j)=b*d/(prnumg*(1.e0+15.e0*rch(i,j)*c))

          end if

        end if

! -----

! Mix the bulk coefficients for the weighted average arrangement ice
! surface.

        if(land(i,j).eq.1) then

          dz0m=za(i,j)-icz0m
          dz0h=za(i,j)-icz0h

          a=max(za(i,j)*icz0mv,1.01e0)
          b=max(za(i,j)*icz0hv,1.01e0)

          if(rch(i,j).lt.0.e0) then

            a=log(a)
            b=log(b)

            c=rch(i,j)*prgiv
            c=c*c

            d=cqg1+cqg2*c
            e=cpg1+cpg2*c

            f=e*e-d*d*d

            c=(za(i,j)*dz0h*a*a)/(dz0m*dz0m*b)

            if(f.gt.0.e0) then

              f=exp(oned3*log(sqrt(f)+abs(e)))

              f=rmg*c*(rmg3v-(f+d/f))

            else

              f=sqrt(d)

              e=max(min(e/(d*f),1.e0),-1.e0)

              f=rmg*c*(rmg3v-2.e0*f*cos(oned3*acos(e)))

            end if

            f=min(f,1.e0)

            e=sqrt(sqrt(1.e0-f))

            c=2.e0*log(.5e0*(1.e0+e))                                   &
     &        +log(.5e0*(1.e0+e*e))-2.e0*atan(e)+cc05

            f=sqrt(1.e0-.6e0*f)

            d=2.e0*log(.5e0*(1.e0+f))

            if(c.lt..5e0*a) then
              cmice=kappa/(a-c)
            else
              cmice=kp2/a
            end if

            if(d.lt..7e0*b) then
              chice=kprg/(b-d)
            else
              chice=kprg3/b
            end if

          else

            a=kappa/log(a)
            b=kappa/log(b)

            c=sqrt(1.e0+5.e0*rch(i,j))
            d=sqrt(1.e0+10.e0*rch(i,j)*c)

            cmice=a/d
            chice=b*d/(prnumg*(1.e0+15.e0*rch(i,j)*c))

          end if

          a=1.e0-kai(i,j)

!ORIG     cm(i,j)=rkp*kai(i,j)*cmice+a*cm(i,j)
!ORIG     ch(i,j)=rkp*kai(i,j)*chice+a*ch(i,j)

          cm(i,j)=kai(i,j)*cmice+a*cm(i,j)
          ch(i,j)=kai(i,j)*chice+a*ch(i,j)

        end if

! -----

      end do
      end do

!$omp end do

!$omp end parallel

!! -----

      end subroutine s_bulksfc

!-----7--------------------------------------------------------------7--

      end module m_bulksfc
