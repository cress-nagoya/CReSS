!***********************************************************************
      module m_getetop
!***********************************************************************

!     Author      : Mizutani Fumihiko, Sakakibara Atsushi
!     Date        : 2004/10/29
!     Modification: 2005/02/25, 2007/05/07, 2007/06/04, 2007/07/30,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2009/03/23, 2009/08/20, 2009/11/13, 2010/05/17,
!                   2011/09/22, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the z physical coordinates at scalar points and radar echo top
!     and total precipitation mixing ratio of radar data.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_getiname
      use m_temparam

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getetop, s_getetop

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getetop

        module procedure s_getetop

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic max

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_getetop(fpngropt,fphaiopt,ngrdmp,rtinc,ni,nj,nk,     &
     &                     nqw,nqi,zph,qwrdr,qwrtd,qirdr,qirtd,         &
     &                     zph8s,etop,qprdr)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpngropt
                       ! Formal parameter of unique index of ngropt

      integer, intent(in) :: fphaiopt
                       ! Formal parameter of unique index of haiopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

      integer, intent(in) :: nqi
                       ! Number of categories of ice hydrometeor

      real, intent(in) :: ngrdmp(1:2)
                       ! Analysis nudging damping coefficient for radar

      real, intent(in) :: rtinc(1:2)
                       ! Lapse of forecast time from radar data reading

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(in) :: qwrdr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor of radar data
                       ! at marked time

      real, intent(in) :: qwrtd(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Time tendency of
                       ! water hydrometeor of radar data

      real, intent(in) :: qirdr(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor of radar data
                       ! at marked time

      real, intent(in) :: qirtd(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Time tendency of
                       ! ice hydrometeor of radar data

! Output variables

      real, intent(out) :: zph8s(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates at scalar points

      real, intent(out) :: etop(0:ni+1,0:nj+1)
                       ! z physical coordinates at radar echo top

      real, intent(out) :: qprdr(0:ni+1,0:nj+1,1:nk)
                       ! Total precipitation mixing ratio of radar data
                       ! at marked time

! Internal shared variables

      integer ngropt   ! Option for analysis nudging to radar
      integer haiopt   ! Option for additional hail processes

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpngropt,ngropt)
      call getiname(fphaiopt,haiopt)

! -----

!! Get the z physical coordinates at scalar points and radar echo top
!! and total precipitation mixing ratio of radar data.

!$omp parallel default(shared) private(k)

! Initialize the array etop with undefined value.

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=1,ni-1
        etop(i,j)=lim36n
      end do
      end do

!$omp end do

! -----

! Perform calculation.

      if(ngropt.eq.1.and.ngrdmp(1).gt.0.e0) then

        if(haiopt.eq.0) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1

              if(qwrdr(i,j,k,2).gt.lim34n.and.                          &
     &           qwrtd(i,j,k,2).gt.lim34n.and.                          &
     &           qirdr(i,j,k,2).gt.lim34n.and.                          &
     &           qirtd(i,j,k,2).gt.lim34n.and.                          &
     &           qirdr(i,j,k,3).gt.lim34n.and.                          &
     &           qirtd(i,j,k,3).gt.lim34n) then

                zph8s(i,j,k)=.5e0*(zph(i,j,k)+zph(i,j,k+1))

                qprdr(i,j,k)=qwrdr(i,j,k,2)                             &
     &            +qirdr(i,j,k,2)+qirdr(i,j,k,3)

                qprdr(i,j,k)=qprdr(i,j,k)+rtinc(1)*(qwrtd(i,j,k,2)      &
     &            +qirtd(i,j,k,2)+qirtd(i,j,k,3))

                if(qprdr(i,j,k).gt.qpmin) then

                  etop(i,j)=max(etop(i,j),zph8s(i,j,k))

                end if

              else

                qprdr(i,j,k)=lim35n

              end if

            end do
            end do

!$omp end do

          end do

        else

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1

              if(qwrdr(i,j,k,2).gt.lim34n.and.                          &
     &           qwrtd(i,j,k,2).gt.lim34n.and.                          &
     &           qirdr(i,j,k,2).gt.lim34n.and.                          &
     &           qirtd(i,j,k,2).gt.lim34n.and.                          &
     &           qirdr(i,j,k,3).gt.lim34n.and.                          &
     &           qirtd(i,j,k,3).gt.lim34n.and.                          &
     &           qirdr(i,j,k,4).gt.lim34n.and.                          &
     &           qirtd(i,j,k,4).gt.lim34n) then

                zph8s(i,j,k)=.5e0*(zph(i,j,k)+zph(i,j,k+1))

                qprdr(i,j,k)=qwrdr(i,j,k,2)                             &
     &            +qirdr(i,j,k,2)+qirdr(i,j,k,3)+qirdr(i,j,k,4)

                qprdr(i,j,k)=qprdr(i,j,k)+rtinc(1)*(qwrtd(i,j,k,2)      &
     &            +qirtd(i,j,k,2)+qirtd(i,j,k,3)+qirtd(i,j,k,4))

                if(qprdr(i,j,k).gt.qpmin) then

                  etop(i,j)=max(etop(i,j),zph8s(i,j,k))

                end if

              else

                qprdr(i,j,k)=lim35n

              end if

            end do
            end do

!$omp end do

          end do

        end if

      else if(ngropt.ge.2.and.ngrdmp(2).gt.0.e0) then

        if(haiopt.eq.0) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1

              if(qwrdr(i,j,k,2).gt.lim34n.and.                          &
     &           qirdr(i,j,k,2).gt.lim34n.and.                          &
     &           qirdr(i,j,k,3).gt.lim34n) then

                zph8s(i,j,k)=.5e0*(zph(i,j,k)+zph(i,j,k+1))

                qprdr(i,j,k)=qwrdr(i,j,k,2)                             &
     &            +qirdr(i,j,k,2)+qirdr(i,j,k,3)

                if(qprdr(i,j,k).gt.qpmin) then

                  etop(i,j)=max(etop(i,j),zph8s(i,j,k))

                end if

              else

                qprdr(i,j,k)=lim35n

              end if

            end do
            end do

!$omp end do

          end do

        else

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1

              if(qwrdr(i,j,k,2).gt.lim34n.and.                          &
     &           qirdr(i,j,k,2).gt.lim34n.and.                          &
     &           qirdr(i,j,k,3).gt.lim34n.and.                          &
     &           qirdr(i,j,k,4).gt.lim34n) then

                zph8s(i,j,k)=.5e0*(zph(i,j,k)+zph(i,j,k+1))

                qprdr(i,j,k)=qwrdr(i,j,k,2)                             &
     &            +qirdr(i,j,k,2)+qirdr(i,j,k,3)+qirdr(i,j,k,4)

                if(qprdr(i,j,k).gt.qpmin) then

                  etop(i,j)=max(etop(i,j),zph8s(i,j,k))

                end if

              else

                qprdr(i,j,k)=lim35n

              end if

            end do
            end do

!$omp end do

          end do

        end if

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_getetop

!-----7--------------------------------------------------------------7--

      end module m_getetop
