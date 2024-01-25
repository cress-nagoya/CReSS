!***********************************************************************
      module m_chksat
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2008/07/01
!     Modification: 2008/08/25, 2009/02/27, 2009/11/13, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     check and avoid the super saturation mixing ratio.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy
      use m_getindx

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: chksat, s_chksat

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface chksat

        module procedure s_chksat

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
      subroutine s_chksat(fproc,xo,imin,imax,jmin,jmax,kmin,kmax,       &
     &                    pbr,ptbr,pp,ptp,qv)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fproc
                       ! Control flag of processing type

      character(len=3), intent(in) :: xo
                       ! Control flag of variable arrangement

      integer, intent(in) :: imin
                       ! Minimum array index in x direction

      integer, intent(in) :: imax
                       ! Maximum array index in x direction

      integer, intent(in) :: jmin
                       ! Minimum array index in y direction

      integer, intent(in) :: jmax
                       ! Maximum array index in y direction

      integer, intent(in) :: kmin
                       ! Minimum array index in z direction

      integer, intent(in) :: kmax
                       ! Maximum array index in z direction

      real, intent(in) :: pbr(imin:imax,jmin:jmax,kmin:kmax)
                       ! Base state pressure

      real, intent(in) :: ptbr(imin:imax,jmin:jmax,kmin:kmax)
                       ! Base state potential temperature

      real, intent(in) :: pp(imin:imax,jmin:jmax,kmin:kmax)
                       ! Pressure perturbation

      real, intent(in) :: ptp(imin:imax,jmin:jmax,kmin:kmax)
                       ! Potential temperature perturbation

! Input and output variable

      real, intent(inout) :: qv(imin:imax,jmin:jmax,kmin:kmax)
                       ! Water vapor mixing ratio

! Internal shared variables

      integer istr     ! Minimum do loops index in x direction
      integer iend     ! Maximum do loops index in x direction
      integer jstr     ! Minimum do loops index in y direction
      integer jend     ! Maximum do loops index in y direction
      integer kstr     ! Minimum do loops index in z direction
      integer kend     ! Maximum do loops index in z direction

      real rddvcp      ! rd / cp

      real p0iv        ! 1.0 / p0

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real t           ! Temperature

      real es          ! Saturation vapor pressure

      real qvs         ! Saturation mixing ratio

!-----7--------------------------------------------------------------7--

! Get the maximum and minimim indices of do loops.

      call getindx(xo,imin,imax,jmin,jmax,istr,iend,jstr,jend)

      if(xo(3:3).eq.'o') then
        kstr=1
        kend=kmax
      else if(xo(3:3).eq.'x') then
        kstr=2
        kend=kmax-2
      end if

! -----

! Set the common used variables.

      rddvcp=rd/cp

      p0iv=1.e0/p0

! -----

! Check and avoid the super saturation mixing ratio.

!$omp parallel default(shared) private(k)

      if(fproc(1:3).eq.'bar') then

        do k=kstr,kend

!$omp do schedule(runtime) private(i,j,t,es,qvs)

          do j=jstr,jend
          do i=istr,iend

            t=ptbr(i,j,k)*exp(rddvcp*log(p0iv*pbr(i,j,k)))

            if(t.gt.tlow) then

              es=es0*exp(17.269e0*(t-t0)/(t-35.86e0))

              qvs=epsva*es/(pbr(i,j,k)-es)

            else

              es=es0*exp(21.875e0*(t-t0)/(t-7.66e0))

              qvs=epsva*es/(pbr(i,j,k)-es)

            end if

            qv(i,j,k)=min(qv(i,j,k),qvs)

          end do
          end do

!$omp end do

        end do

      else if(fproc(1:5).eq.'total') then

        do k=kstr,kend

!$omp do schedule(runtime) private(i,j,t,es,qvs)

          do j=jstr,jend
          do i=istr,iend

            t=(ptbr(i,j,k)+ptp(i,j,k))                                  &
     &        *exp(rddvcp*log(p0iv*(pbr(i,j,k)+pp(i,j,k))))

            if(t.gt.tlow) then

              es=es0*exp(17.269e0*(t-t0)/(t-35.86e0))

              qvs=epsva*es/((pbr(i,j,k)+pp(i,j,k))-es)

            else

              es=es0*exp(21.875e0*(t-t0)/(t-7.66e0))

              qvs=epsva*es/((pbr(i,j,k)+pp(i,j,k))-es)

            end if

            qv(i,j,k)=min(qv(i,j,k),qvs)

          end do
          end do

!$omp end do

        end do

      end if

!$omp end parallel

! -----

      end subroutine s_chksat

!-----7--------------------------------------------------------------7--

      end module m_chksat
