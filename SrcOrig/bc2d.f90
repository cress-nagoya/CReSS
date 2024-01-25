!***********************************************************************
      module m_bc2d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/07/28,
!                   1999/08/03, 1999/08/18, 1999/08/23, 1999/09/01,
!                   1999/09/30, 1999/10/07, 1999/10/19, 2000/01/17,
!                   2001/12/11, 2002/04/02, 2003/04/30, 2003/05/19,
!                   2003/09/01, 2005/02/10, 2006/12/04, 2007/01/05,
!                   2007/01/31, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the boundary conditions for 2 dimensional variables.

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

      public :: bc2d, s_bc2d

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface bc2d

        module procedure s_bc2d

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
      subroutine s_bc2d(fpwbc,fpebc,fpsbc,fpnbc,ni,nj,var2d)
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

! Input and output variable

      real, intent(inout) :: var2d(0:ni+1,0:nj+1)
                       ! Optional 2 dimensional variable

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

!!! Set the lateral boundary conditions.

!! Set the boundary conditions at the opened sections.

!$omp parallel default(shared)

! Set the west boundary conditions.

      if(ebw.eq.1.and.isub.eq.0) then

        if(wbc.eq.2.or.wbc.eq.3) then

!$omp do schedule(runtime) private(j)

          do j=0,nj+1
            var2d(0,j)=var2d(2,j)
            var2d(1,j)=var2d(2,j)
          end do

!$omp end do

        else if(wbc.ge.4) then

!$omp do schedule(runtime) private(j)

          do j=0,nj+1
            var2d(0,j)=var2d(3,j)
            var2d(1,j)=var2d(3,j)
            var2d(2,j)=var2d(3,j)
          end do

!$omp end do

        end if

      end if

! -----

! Set the east boundary conditions.

      if(ebe.eq.1.and.isub.eq.nisub-1) then

        if(ebc.eq.2.or.ebc.eq.3) then

!$omp do schedule(runtime) private(j)

          do j=0,nj+1
            var2d(ni,j)=var2d(nim2,j)
            var2d(nim1,j)=var2d(nim2,j)
          end do

!$omp end do

        else if(ebc.ge.4) then

!$omp do schedule(runtime) private(j)

          do j=0,nj+1
            var2d(ni,j)=var2d(nim3,j)
            var2d(nim1,j)=var2d(nim3,j)
            var2d(nim2,j)=var2d(nim3,j)
          end do

!$omp end do

        end if

      end if

! -----

! Set the south boundary conditions.

      if(ebs.eq.1.and.jsub.eq.0) then

        if(sbc.eq.2.or.sbc.eq.3) then

!$omp do schedule(runtime) private(i)

          do i=0,ni+1
            var2d(i,0)=var2d(i,2)
            var2d(i,1)=var2d(i,2)
          end do

!$omp end do

        else if(sbc.ge.4) then

!$omp do schedule(runtime) private(i)

          do i=0,ni+1
            var2d(i,0)=var2d(i,3)
            var2d(i,1)=var2d(i,3)
            var2d(i,2)=var2d(i,3)
          end do

!$omp end do

        end if

      end if

! -----

! Set the north boundary conditions.

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(nbc.eq.2.or.nbc.eq.3) then

!$omp do schedule(runtime) private(i)

          do i=0,ni+1
            var2d(i,nj)=var2d(i,njm2)
            var2d(i,njm1)=var2d(i,njm2)
          end do

!$omp end do

        else if(nbc.ge.4) then

!$omp do schedule(runtime) private(i)

          do i=0,ni+1
            var2d(i,nj)=var2d(i,njm3)
            var2d(i,njm1)=var2d(i,njm3)
            var2d(i,njm2)=var2d(i,njm3)
          end do

!$omp end do

        end if

      end if

! -----

!$omp end parallel

!! -----

!! Set the boundary conditions at the closed corners.

! Set the boundary conditions at the south-west closed corner.

      if(ebsw.eq.1.and.isub.eq.0.and.jsub.eq.0) then

        if(wbc.eq.2.or.wbc.eq.3) then

          if(sbc.eq.2.or.sbc.eq.3) then

            var2d(0,0)=var2d(2,2)
            var2d(1,0)=var2d(2,2)
            var2d(2,0)=var2d(2,2)
            var2d(0,1)=var2d(2,2)
            var2d(1,1)=var2d(2,2)
            var2d(2,1)=var2d(2,2)
            var2d(0,2)=var2d(2,2)
            var2d(1,2)=var2d(2,2)

          else if(sbc.ge.4) then

            var2d(0,0)=var2d(2,3)
            var2d(1,0)=var2d(2,3)
            var2d(2,0)=var2d(2,3)
            var2d(0,1)=var2d(2,3)
            var2d(1,1)=var2d(2,3)
            var2d(2,1)=var2d(2,3)
            var2d(0,2)=var2d(2,3)
            var2d(1,2)=var2d(2,3)
            var2d(2,2)=var2d(2,3)
            var2d(0,3)=var2d(2,3)
            var2d(1,3)=var2d(2,3)

          end if

        else if(wbc.ge.4) then

          if(sbc.eq.2.or.sbc.eq.3) then

            var2d(0,0)=var2d(3,2)
            var2d(1,0)=var2d(3,2)
            var2d(2,0)=var2d(3,2)
            var2d(3,0)=var2d(3,2)
            var2d(0,1)=var2d(3,2)
            var2d(1,1)=var2d(3,2)
            var2d(2,1)=var2d(3,2)
            var2d(3,1)=var2d(3,2)
            var2d(0,2)=var2d(3,2)
            var2d(1,2)=var2d(3,2)
            var2d(2,2)=var2d(3,2)

          else if(sbc.ge.4) then

            var2d(0,0)=var2d(3,3)
            var2d(1,0)=var2d(3,3)
            var2d(2,0)=var2d(3,3)
            var2d(3,0)=var2d(3,3)
            var2d(0,1)=var2d(3,3)
            var2d(1,1)=var2d(3,3)
            var2d(2,1)=var2d(3,3)
            var2d(3,1)=var2d(3,3)
            var2d(0,2)=var2d(3,3)
            var2d(1,2)=var2d(3,3)
            var2d(2,2)=var2d(3,3)
            var2d(3,2)=var2d(3,3)
            var2d(0,3)=var2d(3,3)
            var2d(1,3)=var2d(3,3)
            var2d(2,3)=var2d(3,3)

          end if

        end if

      end if

! -----

! Set the boundary conditions at the south-east closed corner.

      if(ebse.eq.1.and.isub.eq.nisub-1.and.jsub.eq.0) then

        if(ebc.eq.2.or.ebc.eq.3) then

          if(sbc.eq.2.or.sbc.eq.3) then

            var2d(ni,0)=var2d(nim2,2)
            var2d(nim1,0)=var2d(nim2,2)
            var2d(nim2,0)=var2d(nim2,2)
            var2d(ni,1)=var2d(nim2,2)
            var2d(nim1,1)=var2d(nim2,2)
            var2d(nim2,1)=var2d(nim2,2)
            var2d(ni,2)=var2d(nim2,2)
            var2d(nim1,2)=var2d(nim2,2)

          else if(sbc.ge.4) then

            var2d(ni,0)=var2d(nim2,3)
            var2d(nim1,0)=var2d(nim2,3)
            var2d(nim2,0)=var2d(nim2,3)
            var2d(ni,1)=var2d(nim2,3)
            var2d(nim1,1)=var2d(nim2,3)
            var2d(nim2,1)=var2d(nim2,3)
            var2d(ni,2)=var2d(nim2,3)
            var2d(nim1,2)=var2d(nim2,3)
            var2d(nim2,2)=var2d(nim2,3)
            var2d(ni,3)=var2d(nim2,3)
            var2d(nim1,3)=var2d(nim2,3)

          end if

        else if(ebc.ge.4) then

          if(sbc.eq.2.or.sbc.eq.3) then

            var2d(ni,0)=var2d(nim3,2)
            var2d(nim1,0)=var2d(nim3,2)
            var2d(nim2,0)=var2d(nim3,2)
            var2d(nim3,0)=var2d(nim3,2)
            var2d(ni,1)=var2d(nim3,2)
            var2d(nim1,1)=var2d(nim3,2)
            var2d(nim2,1)=var2d(nim3,2)
            var2d(nim3,1)=var2d(nim3,2)
            var2d(ni,2)=var2d(nim3,2)
            var2d(nim1,2)=var2d(nim3,2)
            var2d(nim2,2)=var2d(nim3,2)

          else if(sbc.ge.4) then

            var2d(ni,0)=var2d(nim3,3)
            var2d(nim1,0)=var2d(nim3,3)
            var2d(nim2,0)=var2d(nim3,3)
            var2d(nim3,0)=var2d(nim3,3)
            var2d(ni,1)=var2d(nim3,3)
            var2d(nim1,1)=var2d(nim3,3)
            var2d(nim2,1)=var2d(nim3,3)
            var2d(nim3,1)=var2d(nim3,3)
            var2d(ni,2)=var2d(nim3,3)
            var2d(nim1,2)=var2d(nim3,3)
            var2d(nim2,2)=var2d(nim3,3)
            var2d(nim3,2)=var2d(nim3,3)
            var2d(ni,3)=var2d(nim3,3)
            var2d(nim1,3)=var2d(nim3,3)
            var2d(nim2,3)=var2d(nim3,3)

          end if

        end if

      end if

! -----

! Set the boundary conditions at the north-west closed corner.

      if(ebnw.eq.1.and.isub.eq.0.and.jsub.eq.njsub-1) then

        if(wbc.eq.2.or.wbc.eq.3) then

          if(nbc.eq.2.or.nbc.eq.3) then

            var2d(0,nj)=var2d(2,njm2)
            var2d(1,nj)=var2d(2,njm2)
            var2d(2,nj)=var2d(2,njm2)
            var2d(0,njm1)=var2d(2,njm2)
            var2d(1,njm1)=var2d(2,njm2)
            var2d(2,njm1)=var2d(2,njm2)
            var2d(0,njm2)=var2d(2,njm2)
            var2d(1,njm2)=var2d(2,njm2)

          else if(nbc.ge.4) then

            var2d(0,nj)=var2d(2,njm3)
            var2d(1,nj)=var2d(2,njm3)
            var2d(2,nj)=var2d(2,njm3)
            var2d(0,njm1)=var2d(2,njm3)
            var2d(1,njm1)=var2d(2,njm3)
            var2d(2,njm1)=var2d(2,njm3)
            var2d(0,njm2)=var2d(2,njm3)
            var2d(1,njm2)=var2d(2,njm3)
            var2d(2,njm2)=var2d(2,njm3)
            var2d(0,njm3)=var2d(2,njm3)
            var2d(1,njm3)=var2d(2,njm3)

          end if

        else if(wbc.ge.4) then

          if(nbc.eq.2.or.nbc.eq.3) then

            var2d(0,nj)=var2d(3,njm2)
            var2d(1,nj)=var2d(3,njm2)
            var2d(2,nj)=var2d(3,njm2)
            var2d(3,nj)=var2d(3,njm2)
            var2d(0,njm1)=var2d(3,njm2)
            var2d(1,njm1)=var2d(3,njm2)
            var2d(2,njm1)=var2d(3,njm2)
            var2d(3,njm1)=var2d(3,njm2)
            var2d(0,njm2)=var2d(3,njm2)
            var2d(1,njm2)=var2d(3,njm2)
            var2d(2,njm2)=var2d(3,njm2)

          else if(nbc.ge.4) then

            var2d(0,nj)=var2d(3,njm3)
            var2d(1,nj)=var2d(3,njm3)
            var2d(2,nj)=var2d(3,njm3)
            var2d(3,nj)=var2d(3,njm3)
            var2d(0,njm1)=var2d(3,njm3)
            var2d(1,njm1)=var2d(3,njm3)
            var2d(2,njm1)=var2d(3,njm3)
            var2d(3,njm1)=var2d(3,njm3)
            var2d(0,njm2)=var2d(3,njm3)
            var2d(1,njm2)=var2d(3,njm3)
            var2d(2,njm2)=var2d(3,njm3)
            var2d(3,njm2)=var2d(3,njm3)
            var2d(0,njm3)=var2d(3,njm3)
            var2d(1,njm3)=var2d(3,njm3)
            var2d(2,njm3)=var2d(3,njm3)

          end if

        end if

      end if

! -----

! Set the boundary conditions at the north-east closed corner.

      if(ebne.eq.1.and.isub.eq.nisub-1.and.jsub.eq.njsub-1) then

        if(ebc.eq.2.or.ebc.eq.3) then

          if(nbc.eq.2.or.nbc.eq.3) then

            var2d(ni,nj)=var2d(nim2,njm2)
            var2d(nim1,nj)=var2d(nim2,njm2)
            var2d(nim2,nj)=var2d(nim2,njm2)
            var2d(ni,njm1)=var2d(nim2,njm2)
            var2d(nim1,njm1)=var2d(nim2,njm2)
            var2d(nim2,njm1)=var2d(nim2,njm2)
            var2d(ni,njm2)=var2d(nim2,njm2)
            var2d(nim1,njm2)=var2d(nim2,njm2)

          else if(nbc.ge.4) then

            var2d(ni,nj)=var2d(nim2,njm3)
            var2d(nim1,nj)=var2d(nim2,njm3)
            var2d(nim2,nj)=var2d(nim2,njm3)
            var2d(ni,njm1)=var2d(nim2,njm3)
            var2d(nim1,njm1)=var2d(nim2,njm3)
            var2d(nim2,njm1)=var2d(nim2,njm3)
            var2d(ni,njm2)=var2d(nim2,njm3)
            var2d(nim1,njm2)=var2d(nim2,njm3)
            var2d(nim2,njm2)=var2d(nim2,njm3)
            var2d(ni,njm3)=var2d(nim2,njm3)
            var2d(nim1,njm3)=var2d(nim2,njm3)

          end if

        else if(ebc.ge.4) then

          if(nbc.eq.2.or.nbc.eq.3) then

            var2d(ni,nj)=var2d(nim3,njm2)
            var2d(nim1,nj)=var2d(nim3,njm2)
            var2d(nim2,nj)=var2d(nim3,njm2)
            var2d(nim3,nj)=var2d(nim3,njm2)
            var2d(ni,njm1)=var2d(nim3,njm2)
            var2d(nim1,njm1)=var2d(nim3,njm2)
            var2d(nim2,njm1)=var2d(nim3,njm2)
            var2d(nim3,njm1)=var2d(nim3,njm2)
            var2d(ni,njm2)=var2d(nim3,njm2)
            var2d(nim1,njm2)=var2d(nim3,njm2)
            var2d(nim2,njm2)=var2d(nim3,njm2)

          else if(nbc.ge.4) then

            var2d(ni,nj)=var2d(nim3,njm3)
            var2d(nim1,nj)=var2d(nim3,njm3)
            var2d(nim2,nj)=var2d(nim3,njm3)
            var2d(nim3,nj)=var2d(nim3,njm3)
            var2d(ni,njm1)=var2d(nim3,njm3)
            var2d(nim1,njm1)=var2d(nim3,njm3)
            var2d(nim2,njm1)=var2d(nim3,njm3)
            var2d(nim3,njm1)=var2d(nim3,njm3)
            var2d(ni,njm2)=var2d(nim3,njm3)
            var2d(nim1,njm2)=var2d(nim3,njm3)
            var2d(nim2,njm2)=var2d(nim3,njm3)
            var2d(nim3,njm2)=var2d(nim3,njm3)
            var2d(ni,njm3)=var2d(nim3,njm3)
            var2d(nim1,njm3)=var2d(nim3,njm3)
            var2d(nim2,njm3)=var2d(nim3,njm3)

          end if

        end if

      end if

! -----

!! -----

!!! -----

      end subroutine s_bc2d

!-----7--------------------------------------------------------------7--

      end module m_bc2d
