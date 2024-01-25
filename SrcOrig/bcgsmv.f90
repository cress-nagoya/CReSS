!***********************************************************************
      module m_bcgsmv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/11/19
!     Modification: 2000/01/17, 2001/12/11, 2002/04/02, 2002/07/23,
!                   2003/04/30, 2003/05/19, 2003/12/12, 2004/03/05,
!                   2004/04/15, 2005/02/10, 2006/12/04, 2007/01/05,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the boundary conditions for the y components of velocity
!     in GPV data smoothing.

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

      public :: bcgsmv, s_bcgsmv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface bcgsmv

        module procedure s_bcgsmv

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
      subroutine s_bcgsmv(fpwbc,fpebc,fpsbc,fpnbc,ni,nj,nk,dfv)
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

      real, intent(inout) :: dfv(0:ni+1,0:nj+1,1:nk)
                       ! Diffusion term
                       ! of y components of velocity of GPV

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions
      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

      integer nim1     ! ni - 1
      integer nim2     ! ni - 2
      integer njm1     ! nj - 1
      integer nkm2     ! nk - 2
      integer nkm3     ! nk - 3

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
      nkm2=nk-2
      nkm3=nk-3

! -----

!! Set the boundary conditions.

!$omp parallel default(shared) private(k)

! Set the west boundary conditions.

      if(ebw.eq.1.and.isub.eq.0.and.abs(wbc).ne.1) then

        do k=2,nk-2

!$omp do schedule(runtime) private(j)

          do j=1,nj
            dfv(1,j,k)=dfv(2,j,k)
          end do

!$omp end do

        end do

      end if

! -----

! Set the east boundary conditions.

      if(ebe.eq.1.and.isub.eq.nisub-1.and.abs(ebc).ne.1) then

        do k=2,nk-2

!$omp do schedule(runtime) private(j)

          do j=1,nj
            dfv(nim1,j,k)=dfv(nim2,j,k)
          end do

!$omp end do

        end do

      end if

! -----

! Set the south boundary conditions.

      if(ebs.eq.1.and.jsub.eq.0.and.abs(sbc).ne.1) then

        do k=2,nk-2

!$omp do schedule(runtime) private(i)

          do i=1,ni-1
            dfv(i,1,k)=dfv(i,2,k)
          end do

!$omp end do

        end do

      end if

! -----

! Set the north boundary conditions.

      if(ebn.eq.1.and.jsub.eq.njsub-1.and.abs(nbc).ne.1) then

        do k=2,nk-2

!$omp do schedule(runtime) private(i)

          do i=1,ni-1
            dfv(i,nj,k)=dfv(i,njm1,k)
          end do

!$omp end do

        end do

      end if

! -----

! Set the bottom and top boundary conditions.

!$omp do schedule(runtime) private(i,j)

      do j=1,nj
      do i=1,ni-1
        dfv(i,j,2)=dfv(i,j,3)
        dfv(i,j,nkm2)=dfv(i,j,nkm3)
      end do
      end do

!$omp end do

! -----

!$omp end parallel

!! -----

      end subroutine s_bcgsmv

!-----7--------------------------------------------------------------7--

      end module m_bcgsmv
