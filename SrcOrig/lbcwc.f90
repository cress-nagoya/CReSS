!***********************************************************************
      module m_lbcwc
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/06/21
!     Modification: 2006/12/04, 2007/01/05, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the lateral boundary conditions for the zeta components of
!     contravariant velocity.

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

      public :: lbcwc, s_lbcwc

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface lbcwc

        module procedure s_lbcwc

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

!     none

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_lbcwc(fpwbc,fpebc,fpsbc,fpnbc,ni,nj,nk,wc)
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

! Input and output variable

      real, intent(inout) :: wc(0:ni+1,0:nj+1,1:nk)
                       ! zeta components of contravariant velocity

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions
      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

      integer nim1     ! ni - 1
      integer nim2     ! ni - 2
      integer njm1     ! nj - 1
      integer njm2     ! nj - 2

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
      njm1=nj-1
      njm2=nj-2

! -----

!! Set the lateral boundary conditions.

!$omp parallel default(shared) private(k)

! Set the west boundary conditions.

      if(ebw.eq.1.and.isub.eq.0) then

        if(wbc.eq.2) then

          do k=2,nk-1

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              wc(1,j,k)=wc(2,j,k)
            end do

!$omp end do

          end do

        else if(wbc.ge.3) then

          do k=2,nk-1

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              wc(1,j,k)=wc(2,j,k)
            end do

!$omp end do

          end do

        end if

      end if

! -----

! Set the east boundary conditions.

      if(ebe.eq.1.and.isub.eq.nisub-1) then

        if(ebc.eq.2) then

          do k=2,nk-1

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              wc(nim1,j,k)=wc(nim2,j,k)
            end do

!$omp end do

          end do

        else if(ebc.ge.3) then

          do k=2,nk-1

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              wc(nim1,j,k)=wc(nim2,j,k)
            end do

!$omp end do

          end do

        end if

      end if

! -----

! Set the south boundary conditions.

      if(ebs.eq.1.and.jsub.eq.0) then

        if(sbc.eq.2) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              wc(i,1,k)=wc(i,2,k)
            end do

!$omp end do

          end do

        else if(sbc.ge.3) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              wc(i,1,k)=wc(i,2,k)
            end do

!$omp end do

          end do

        end if

      end if

! -----

! Set the north boundary conditions.

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(nbc.eq.2) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              wc(i,njm1,k)=wc(i,njm2,k)
            end do

!$omp end do

          end do

        else if(nbc.ge.3) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              wc(i,njm1,k)=wc(i,njm2,k)
            end do

!$omp end do

          end do

        end if

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_lbcwc

!-----7--------------------------------------------------------------7--

      end module m_lbcwc
