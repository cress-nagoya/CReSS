!***********************************************************************
      module m_lbcu
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/07/28,
!                   1999/08/03, 1999/08/18, 1999/08/23, 1999/09/06,
!                   1999/09/30, 1999/10/07, 1999/10/27, 1999/11/01,
!                   1999/12/06, 2000/01/17, 2000/03/17, 2000/03/23,
!                   2001/01/15, 2001/04/15, 2001/07/13, 2001/08/07,
!                   2001/09/13, 2001/12/11, 2002/04/02, 2002/06/06,
!                   2002/07/23, 2002/08/15, 2002/10/31, 2003/04/30,
!                   2003/05/19, 2003/06/27, 2003/11/05, 2003/11/28,
!                   2003/12/12, 2004/01/09, 2004/05/07, 2004/08/20,
!                   2005/01/31, 2006/06/21, 2006/12/04, 2007/01/05,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the lateral boundary conditions for the x components of
!     velocity.

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

      public :: lbcu, s_lbcu

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface lbcu

        module procedure s_lbcu

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
      subroutine s_lbcu(fpwbc,fpebc,fpsbc,fpnbc,fpadvopt,ni,nj,nk,uf)
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

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

! Input and output variable

      real, intent(inout) :: uf(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at future

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions
      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

      integer advopt   ! Option for advection scheme

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
      call getiname(fpadvopt,advopt)

! -----

! Set the common used variables.

      nim2=ni-2
      njm1=nj-1
      njm2=nj-2

! -----

!! Set the lateral boundary conditions.

!$omp parallel default(shared) private(k)

! Set the west boundary conditions.

      if(ebw.eq.1.and.isub.eq.0) then

        if(advopt.le.3) then

          if(wbc.eq.2) then

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                uf(1,j,k)=-uf(3,j,k)
              end do

!$omp end do

            end do

          else if(wbc.eq.3) then

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                uf(1,j,k)=uf(3,j,k)
              end do

!$omp end do

            end do

          end if

        else

          if(wbc.eq.2) then

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                uf(1,j,k)=-uf(3,j,k)
              end do

!$omp end do

            end do

          else if(wbc.ge.3) then

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                uf(1,j,k)=uf(3,j,k)
              end do

!$omp end do

            end do

          end if

        end if

      end if

! -----

! Set the east boundary conditions.

      if(ebe.eq.1.and.isub.eq.nisub-1) then

        if(advopt.le.3) then

          if(ebc.eq.2) then

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                uf(ni,j,k)=-uf(nim2,j,k)
              end do

!$omp end do

            end do

          else if(ebc.eq.3) then

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                uf(ni,j,k)=uf(nim2,j,k)
              end do

!$omp end do

            end do

          end if

        else

          if(ebc.eq.2) then

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                uf(ni,j,k)=-uf(nim2,j,k)
              end do

!$omp end do

            end do

          else if(ebc.ge.3) then

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                uf(ni,j,k)=uf(nim2,j,k)
              end do

!$omp end do

            end do

          end if

        end if

      end if

! -----

! Set the south boundary conditions.

      if(ebs.eq.1.and.jsub.eq.0) then

        if(advopt.le.3) then

          if(sbc.eq.2) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni
                uf(i,1,k)=uf(i,2,k)
              end do

!$omp end do

            end do

          else if(sbc.eq.3) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni
                uf(i,1,k)=uf(i,2,k)
              end do

!$omp end do

            end do

          end if

        else

          if(sbc.eq.2) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni
                uf(i,1,k)=uf(i,2,k)
              end do

!$omp end do

            end do

          else if(sbc.ge.3) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni
                uf(i,1,k)=uf(i,2,k)
              end do

!$omp end do

            end do

          end if

        end if

      end if

! -----

! Set the north boundary conditions.

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(advopt.le.3) then

          if(nbc.eq.2) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni
                uf(i,njm1,k)=uf(i,njm2,k)
              end do

!$omp end do

            end do

          else if(nbc.eq.3) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni
                uf(i,njm1,k)=uf(i,njm2,k)
              end do

!$omp end do

            end do

          end if

        else

          if(nbc.eq.2) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni
                uf(i,njm1,k)=uf(i,njm2,k)
              end do

!$omp end do

            end do

          else if(nbc.ge.3) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni
                uf(i,njm1,k)=uf(i,njm2,k)
              end do

!$omp end do

            end do

          end if

        end if

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_lbcu

!-----7--------------------------------------------------------------7--

      end module m_lbcu
