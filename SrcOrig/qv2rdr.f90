!***********************************************************************
      module m_qv2rdr
!***********************************************************************

!     Author      : Mizutani Fumihiko, Sakakibara Atsushi
!     Date        : 2004/03/16
!     Modification: 2004/06/14, 2004/10/28, 2007/05/07, 2007/06/27,
!                   2007/07/30, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2009/03/31, 2009/08/20, 2009/11/13,
!                   2010/05/17, 2010/09/22, 2013/01/28, 2013/02/11,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform the analysis nudging to radar data of water vapor mixing
!     ratio.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy
      use m_getiname
      use m_temparam

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: qv2rdr, s_qv2rdr

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface qv2rdr

        module procedure s_qv2rdr

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic log
      intrinsic min

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_qv2rdr(fpngropt,ngrdmp,ni,nj,nk,zph,zph8s,etop,tlcl, &
     &                    pbr,ptbr,rst,ppp,ptpp,qvp,qvs,qvfrc,tsfc)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpngropt
                       ! Formal parameter of unique index of ngropt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: ngrdmp(1:2)
                       ! Analysis nudging dumping coefficient

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(in) :: zph8s(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates at scalar points

      real, intent(in) :: etop(0:ni+1,0:nj+1)
                       ! z physical coordinates at radar echo top

      real, intent(in) :: tlcl(0:ni+1,0:nj+1)
                       ! Air Temperature at Lifting Condensation Level

      real, intent(in) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jabobian

      real, intent(in) :: ppp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at past

      real, intent(in) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(in) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(in) :: qvs(0:ni+1,0:nj+1,1:nk)
                       ! Saturation mixing ratio

! Input and output variable

      real, intent(inout) :: qvfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in water vapor mixing ratio equation

! Internal shared variables

      integer ngropt   ! Option for analysis nudging to radar

      real rhr         ! 1.0 - rhqp

      real rddvcp      ! rd / cp

      real p0iv        ! 1.0 / p0

      real ngrqv       ! adjqv x ngrdmp

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real t           ! Air temperature
      real p           ! Pressure

      real ctop        ! Highest height of nudging

      real hvwei       ! Vertical weighted coefficient

      real rhndg       ! Relative humidity at current point
                       ! target to nudge

      real, intent(inout) :: tsfc(0:ni+1,0:nj+1)
                       ! Air temperature at lowest plane

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpngropt,ngropt)

! -----

! Set the common used variable.

      if(ngropt.eq.1.and.ngrdmp(1).gt.0.e0) then

        ngrqv=adjqv*ngrdmp(1)

      else if(ngropt.ge.2.and.ngrdmp(2).gt.0.e0) then

        ngrqv=adjqv*ngrdmp(2)

      else

        ngrqv=0.e0

      end if

      rhr=1.e0-rhqp

      rddvcp=rd/cp

      p0iv=1.e0/p0

! -----

! Perform the analysis nudging to radar data of water vapor mixing
! ratio.

!$omp parallel default(shared) private(k)

!$omp do schedule(runtime) private(i,j)

      do j=2,nj-2
      do i=2,ni-2
        tsfc(i,j)=(ptbr(i,j,2)+ptpp(i,j,2))                             &
     &    *exp(rddvcp*log(p0iv*(pbr(i,j,2)+ppp(i,j,2))))
      end do
      end do

!$omp end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j,ctop,hvwei,p,t,rhndg)

        do j=2,nj-2
        do i=2,ni-2

          ctop=min(etop(i,j),qvtop+zph(i,j,2))

!ORIG     if(zph8s(i,j,k).le.etop(i,j)) then
          if(zph8s(i,j,k).le.ctop) then

            p=pbr(i,j,k)+ppp(i,j,k)
            t=(ptbr(i,j,k)+ptpp(i,j,k))*exp(rddvcp*log(p0iv*p))

            if(t.le.tlcl(i,j)) then

              hvwei=1.e0-(zph8s(i,j,k)-zph(i,j,2))/(ctop-zph(i,j,2))

              qvfrc(i,j,k)=qvfrc(i,j,k)                                 &
     &          +ngrqv*hvwei*rst(i,j,k)*(qvs(i,j,k)-qvp(i,j,k))

            else

              rhndg=((tsfc(i,j)-t)/(tsfc(i,j)-tlcl(i,j)))*rhr+rhqp

              qvfrc(i,j,k)=qvfrc(i,j,k)                                 &
     &          +ngrqv*rst(i,j,k)*(rhndg*qvs(i,j,k)-qvp(i,j,k))

            end if

          end if

        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_qv2rdr

!-----7--------------------------------------------------------------7--

      end module m_qv2rdr
