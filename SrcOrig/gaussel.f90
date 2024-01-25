!***********************************************************************
      module m_gaussel
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/12/20
!     Modification: 2000/01/17, 2000/03/23, 2001/01/04, 2001/03/13,
!                   2001/04/15, 2001/05/29, 2001/10/18, 2002/04/02,
!                   2003/01/20, 2003/04/30, 2003/05/19, 2003/12/12,
!                   2004/08/20, 2004/09/01, 2004/09/10, 2005/01/31,
!                   2007/01/31, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     solve the tridiagonal equation with the Gauss elimination.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: gaussel, s_gaussel

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface gaussel

        module procedure s_gaussel

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic int
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_gaussel(fpimpopt,istr,iend,jstr,jend,kstr,kend,ni,nj,&
     &                     kmax,rr,ss,tt,ff,pv)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpimpopt
                       ! Formal parameter of unique index of impopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: kmax
                       ! Maximum array index in z direction

      integer, intent(in) :: istr
                       ! Minimum do loops index in x direction

      integer, intent(in) :: iend
                       ! Maximum do loops index in x direction

      integer, intent(in) :: jstr
                       ! Minimum do loops index in y direction

      integer, intent(in) :: jend
                       ! Maximum do loops index in y direction

      integer, intent(in) :: kstr
                       ! Minimum do loops index in z direction

      integer, intent(in) :: kend
                       ! Maximum do loops index in z direction

! Input and output variables

      real, intent(inout) :: rr(0:ni+1,0:nj+1,1:kmax)
                       ! Coefficient matrix

      real, intent(inout) :: ss(0:ni+1,0:nj+1,1:kmax)
                       ! Coefficient matrix

      real, intent(inout) :: tt(0:ni+1,0:nj+1,1:kmax)
                       ! Coefficient matrix

      real, intent(inout) :: ff(0:ni+1,0:nj+1,1:kmax)
                       ! Solved vector

! Internal shared variables

      integer impopt   ! Option for vertical implicit method

      integer kem1     ! kend - 1
      integer kem2     ! kend - 2

      real, intent(inout) :: pv(0:ni+1,0:nj+1,1:kmax)
                       ! List vector

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer kp       ! Pivoted array index at k
      integer kpm1     ! Pivoted array index at k - 1
      integer kpp1     ! Pivoted array index at k + 1

      real aa          ! Temporary variable

! Remark

!     ss: This variable is also temporary, because it is not used again.

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpimpopt,impopt)

! -----

! Set the common used variables.

      kem1=kend-1
      kem2=kend-2

! -----

!!! Solve the tridiagonal equation.

!$omp parallel default(shared) private(k)

!! Solve the tridiagonal equation with the Gauss elimination.

      if(impopt.eq.1) then

! Perform the forward eliminations.

        if(kstr.eq.kend) then

!$omp do schedule(runtime) private(i,j)

          do j=jstr,jend
          do i=istr,iend
            ff(i,j,kstr)=ff(i,j,kstr)/ss(i,j,kstr)
          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,j,aa)

          do j=jstr,jend
          do i=istr,iend
            aa=1.e0/ss(i,j,kstr)

            ss(i,j,kstr)=tt(i,j,kstr)*aa
            ff(i,j,kstr)=ff(i,j,kstr)*aa

          end do
          end do

!$omp end do

        end if

        if(kstr+2.le.kend) then

          do k=kstr+1,kend-1

!$omp do schedule(runtime) private(i,j,aa)

            do j=jstr,jend
            do i=istr,iend
              aa=1.e0/(ss(i,j,k)-rr(i,j,k)*ss(i,j,k-1))

              ss(i,j,k)=tt(i,j,k)*aa
              ff(i,j,k)=(ff(i,j,k)-rr(i,j,k)*ff(i,j,k-1))*aa

            end do
            end do

!$omp end do

          end do

        end if

        if(kstr+1.le.kend) then

!$omp do schedule(runtime) private(i,j,aa)

          do j=jstr,jend
          do i=istr,iend
            aa=1.e0/(ss(i,j,kend)-rr(i,j,kend)*ss(i,j,kem1))

            ff(i,j,kend)=(ff(i,j,kend)-rr(i,j,kend)*ff(i,j,kem1))*aa

          end do
          end do

!$omp end do

        end if

! -----

! Perform the back substitutions.

        if(kstr+1.le.kend) then

          do k=kend-1,kstr,-1

!$omp do schedule(runtime) private(i,j)

            do j=jstr,jend
            do i=istr,iend
              ff(i,j,k)=ff(i,j,k)-ss(i,j,k)*ff(i,j,k+1)
            end do
            end do

!$omp end do

          end do

        end if

! -----

!! -----

!! Solve the tridiagonal equation with the partial pivoting Gauss
!! elimination.

      else if(impopt.eq.2) then

! Perform the forward eliminations.

!$omp do schedule(runtime) private(i,j,aa)

        do j=2,nj-2
        do i=2,ni-2
          aa=1.e0/ss(i,j,kstr)

          ss(i,j,kstr)=tt(i,j,kstr)*aa
          ff(i,j,kstr)=ff(i,j,kstr)*aa

        end do
        end do

!$omp end do

!$omp do schedule(runtime)

        do k=kstr,kend
          pv(0,0,k)=real(k)+.1e0
        end do

!$omp end do

        do k=kstr,kend

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            pv(i,j,k)=pv(0,0,k)
          end do
          end do

!$omp end do

        end do

        do k=kstr+1,kend-2

!$omp do schedule(runtime) private(i,j,kp,kpm1,aa)

          do j=2,nj-2
          do i=2,ni-2

            if(abs(ss(i,j,int(pv(i,j,k)))).lt.abs(rr(i,j,k+1))) then

              aa=pv(i,j,k)

              pv(i,j,k)=pv(i,j,k+1)

              pv(i,j,k+1)=aa

            end if

            kp=int(pv(i,j,k))
            kpm1=int(pv(i,j,k-1))

            aa=1.e0/(ss(i,j,kp)-rr(i,j,kp)*ss(i,j,kpm1))

            ss(i,j,kp)=tt(i,j,kp)*aa
            ff(i,j,kp)=(ff(i,j,kp)-rr(i,j,kp)*ff(i,j,kpm1))*aa

          end do
          end do

!$omp end do

        end do

!$omp do schedule(runtime) private(i,j,kpm1,aa)

        do j=2,nj-2
        do i=2,ni-2
          kpm1=int(pv(i,j,kem2))

          aa=1.e0/(ss(i,j,kem1)-rr(i,j,kem1)*ss(i,j,kpm1))

          ss(i,j,kem1)=tt(i,j,kem1)*aa
          ff(i,j,kem1)=(ff(i,j,kem1)-rr(i,j,kem1)*ff(i,j,kpm1))*aa

          aa=1.e0/(ss(i,j,kend)-rr(i,j,kend)*ss(i,j,kem1))

          ff(i,j,kend)=(ff(i,j,kend)-rr(i,j,kend)*ff(i,j,kem1))*aa

        end do
        end do

!$omp end do

! -----

! Perform the back substitutions.

        do k=kend-1,kstr,-1

!$omp do schedule(runtime) private(i,j,kp,kpp1)

          do j=2,nj-2
          do i=2,ni-2
            kp=int(pv(i,j,k))
            kpp1=int(pv(i,j,k+1))

            ff(i,j,kp)=ff(i,j,kp)-ss(i,j,kp)*ff(i,j,kpp1)

          end do
          end do

!$omp end do

        end do

! -----

      end if

!! -----

!$omp end parallel

!!! -----

      end subroutine s_gaussel

!-----7--------------------------------------------------------------7--

      end module m_gaussel
