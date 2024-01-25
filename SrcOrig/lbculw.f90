!***********************************************************************
      module m_lbculw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/04/03
!     Modification: 2006/05/12, 2006/06/21, 2006/12/04, 2007/01/05,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the lateral boundary conditions for w components of velocity.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commpi
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: lbculw, s_lbculw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface lbculw

        module procedure s_lbculw

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_lbculw(fpwbc,fpebc,fpsbc,fpnbc,ni,nj,nk,w,           &
     &                    dxt36w,dxt36e,dyt36s,dyt36n)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpsbc
                       ! Formal parameter of unique index of sbc

      integer, intent(in) :: fpnbc
                       ! Formal parameter of unique index of nbc

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

! Input and output variables

      real, intent(inout) :: w(0:ni+1,0:nj+1,1:nk)
                       ! w components of velocity

      real, intent(inout) :: dxt36w(0:ni+1)
                       ! - dxiv^3 x dtb^3 / 6 and 0 on west boundary

      real, intent(inout) :: dxt36e(0:ni+1)
                       ! - dxiv^3 x dtb^3 / 6 and 0 on east boundary

      real, intent(inout) :: dyt36s(0:nj+1)
                       ! - dyiv^3 x dtb^3 / 6 and 0 on south boundary

      real, intent(inout) :: dyt36n(0:nj+1)
                       ! - dyiv^3 x dtb^3 / 6 and 0 on north boundary

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions
      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

      integer nim1     ! ni - 1
      integer nim2     ! ni - 2
      integer nim3     ! ni - 3

      integer njm1     ! nj - 1
      integer njm2     ! nj - 2
      integer njm3     ! nj - 3

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpsbc,sbc)
      call getiname(fpnbc,nbc)

! -----

! Set the common used variables.

      nim1=ni-1
      nim2=ni-2
      nim3=ni-3

      njm1=nj-1
      njm2=nj-2
      njm3=nj-3

! -----

! Set the coefficients with 0.

      if(abs(wbc).ne.1.and.abs(ebc).ne.1) then

        if(ebw.eq.1.and.isub.eq.0) then
          dxt36w(2)=0.e0
        end if

        if(ebe.eq.1.and.isub.eq.nisub-1) then
          dxt36e(ni-2)=0.e0
        end if

      end if

      if(abs(sbc).ne.1.and.abs(nbc).ne.1) then

        if(ebs.eq.1.and.jsub.eq.0) then
          dyt36s(2)=0.e0
        end if

        if(ebn.eq.1.and.jsub.eq.njsub-1) then
          dyt36n(nj-2)=0.e0
        end if

      end if

! -----

!! Set the lateral boundary conditions.

!$omp parallel default(shared) private(k)

! Set the west and east boundary conditions.

      if(abs(wbc).ne.1.and.abs(ebc).ne.1) then

        if(ebw.eq.1.and.isub.eq.0) then

          do k=1,nk

!$omp do schedule(runtime) private(j)

            do j=0,nj
              w(0,j,k)=w(3,j,k)-3.e0*(w(2,j,k)-w(1,j,k))
            end do

!$omp end do

          end do

        end if

        if(ebe.eq.1.and.isub.eq.nisub-1) then

          do k=1,nk

!$omp do schedule(runtime) private(j)

            do j=0,nj
              w(ni,j,k)=w(nim3,j,k)-3.e0*(w(nim2,j,k)-w(nim1,j,k))
            end do

!$omp end do

          end do

        end if

      end if

! -----

! Set the south and north boundary conditions.

      if(abs(sbc).ne.1.and.abs(nbc).ne.1) then

        if(ebs.eq.1.and.jsub.eq.0) then

          do k=1,nk

!$omp do schedule(runtime) private(i)

            do i=0,ni
              w(i,0,k)=w(i,3,k)-3.e0*(w(i,2,k)-w(i,1,k))
            end do

!$omp end do

          end do

        end if

        if(ebn.eq.1.and.jsub.eq.njsub-1) then

          do k=1,nk

!$omp do schedule(runtime) private(i)

            do i=0,ni
              w(i,nj,k)=w(i,njm3,k)-3.e0*(w(i,njm2,k)-w(i,njm1,k))
            end do

!$omp end do

          end do

        end if

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_lbculw

!-----7--------------------------------------------------------------7--

      end module m_lbculw
