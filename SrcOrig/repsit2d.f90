!***********************************************************************
      module m_repsit2d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2005/07/01
!     Modification: 2005/08/05, 2006/09/21, 2006/12/04, 2007/01/05,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     reposition the restructed boundary variables form original restart
!     boundary variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commpi

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: repsit2d, s_repsit2d

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface repsit2d

        module procedure s_repsit2d

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
      subroutine s_repsit2d(xo,istrb,iendb,jstrb,jendb,                 &
     &                      dib,djb,ni,nj,nk,ni_rst,nj_rst,             &
     &                      var2dx,var2dy,var2dx_rst,var2dy_rst)
!***********************************************************************

! Input variables

      character(len=2), intent(in) :: xo
                       ! Control flag of variable arrangement

      integer, intent(in) :: istrb
                       ! Minimum do loops index in x direction
                       ! of lateral boundary

      integer, intent(in) :: iendb
                       ! Maximum do loops index in x direction
                       ! of lateral boundary

      integer, intent(in) :: jstrb
                       ! Minimum do loops index in y direction
                       ! of lateral boundary

      integer, intent(in) :: jendb
                       ! Maximum do loops index in y direction
                       ! of lateral boundary

      integer, intent(in) :: dib
                       ! Differential index to istrb

      integer, intent(in) :: djb
                       ! Differential index to jstrb

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: ni_rst
                       ! Restructed files dimension in x direction

      integer, intent(in) :: nj_rst
                       ! Restructed files dimension in y direction

      real, intent(in) :: var2dx(1:nj,1:nk,1:2)
                       ! Optional boundary variable
                       ! on west and east boundary

      real, intent(in) :: var2dy(1:ni,1:nk,1:2)
                       ! Optional boundary variable
                       ! on south and north boundary

! Input and output variables

      real, intent(inout) :: var2dx_rst(1:nj_rst,1:nk,1:2)
                       ! var2dx in restructed domain

      real, intent(inout) :: var2dy_rst(1:ni_rst,1:nk,1:2)
                       ! var2dx in restructed domain

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Reposition the restructed boundary variables form original restart
! boundary variables.

!$omp parallel default(shared) private(k)

      if(xo(1:2).eq.'ox') then

       if(ebw.eq.1.and.isub.eq.0.and.isub_rst.eq.0) then

         do k=1,nk

!$omp do schedule(runtime) private(j)

           do j=jstrb,jendb
             var2dx_rst(djb+j,k,1)=var2dx(j,k,1)
           end do

!$omp end do

         end do

       end if

       if(ebe.eq.1.and.isub.eq.nisub-1.and.isub_rst.eq.nisub_rst-1) then

         do k=1,nk

!$omp do schedule(runtime) private(j)

           do j=jstrb,jendb
             var2dx_rst(djb+j,k,2)=var2dx(j,k,2)
           end do

!$omp end do

         end do

       end if

       if(ebs.eq.1.and.jsub.eq.0.and.jsub_rst.eq.0) then

         do k=1,nk

!$omp do schedule(runtime) private(i)

           do i=istrb,iendb+1
             var2dy_rst(dib+i,k,1)=var2dy(i,k,1)
           end do

!$omp end do

         end do

       end if

       if(ebn.eq.1.and.jsub.eq.njsub-1.and.jsub_rst.eq.njsub_rst-1) then

         do k=1,nk

!$omp do schedule(runtime) private(i)

           do i=istrb,iendb+1
             var2dy_rst(dib+i,k,2)=var2dy(i,k,2)
           end do

!$omp end do

         end do

       end if

      else if(xo(1:2).eq.'xo') then

       if(ebw.eq.1.and.isub.eq.0.and.isub_rst.eq.0) then

         do k=1,nk

!$omp do schedule(runtime) private(j)

           do j=jstrb,jendb+1
             var2dx_rst(djb+j,k,1)=var2dx(j,k,1)
           end do

!$omp end do

         end do

       end if

       if(ebe.eq.1.and.isub.eq.nisub-1.and.isub_rst.eq.nisub_rst-1) then

         do k=1,nk

!$omp do schedule(runtime) private(j)

           do j=jstrb,jendb+1
             var2dx_rst(djb+j,k,2)=var2dx(j,k,2)
           end do

!$omp end do

         end do

       end if

       if(ebs.eq.1.and.jsub.eq.0.and.jsub_rst.eq.0) then

         do k=1,nk

!$omp do schedule(runtime) private(i)

           do i=istrb,iendb
             var2dy_rst(dib+i,k,1)=var2dy(i,k,1)
           end do

!$omp end do

         end do

       end if

       if(ebn.eq.1.and.jsub.eq.njsub-1.and.jsub_rst.eq.njsub_rst-1) then

         do k=1,nk

!$omp do schedule(runtime) private(i)

           do i=istrb,iendb
             var2dy_rst(dib+i,k,2)=var2dy(i,k,2)
           end do

!$omp end do

         end do

       end if

      else

       if(ebw.eq.1.and.isub.eq.0.and.isub_rst.eq.0) then

         do k=1,nk

!$omp do schedule(runtime) private(j)

           do j=jstrb,jendb
             var2dx_rst(djb+j,k,1)=var2dx(j,k,1)
           end do

!$omp end do

         end do

       end if

       if(ebe.eq.1.and.isub.eq.nisub-1.and.isub_rst.eq.nisub_rst-1) then

         do k=1,nk

!$omp do schedule(runtime) private(j)

           do j=jstrb,jendb
             var2dx_rst(djb+j,k,2)=var2dx(j,k,2)
           end do

!$omp end do

         end do

       end if

       if(ebs.eq.1.and.jsub.eq.0.and.jsub_rst.eq.0) then

         do k=1,nk

!$omp do schedule(runtime) private(i)

           do i=istrb,iendb
             var2dy_rst(dib+i,k,1)=var2dy(i,k,1)
           end do

!$omp end do

         end do

       end if

       if(ebn.eq.1.and.jsub.eq.njsub-1.and.jsub_rst.eq.njsub_rst-1) then

         do k=1,nk

!$omp do schedule(runtime) private(i)

           do i=istrb,iendb
             var2dy_rst(dib+i,k,2)=var2dy(i,k,2)
           end do

!$omp end do

         end do

       end if

      end if

!$omp end parallel

! -----

      end subroutine s_repsit2d

!-----7--------------------------------------------------------------7--

      end module m_repsit2d
