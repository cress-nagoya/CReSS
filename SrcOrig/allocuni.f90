!***********************************************************************
      module m_allocuni
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/05/19
!     Modification: 2003/11/05, 2004/01/09, 2004/07/01, 2004/09/25,
!                   2005/02/10, 2006/01/10, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/01/31, 2007/04/11, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2008/10/10, 2009/01/05,
!                   2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     allocate the array for unite.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_comionum
      use m_commpi
      use m_comuni
      use m_cpondpe
      use m_destroy
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: allocuni, s_allocuni

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface allocuni

        module procedure s_allocuni

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_allocuni(fpuniopt_uni,ni,nj,nk,ni_uni,nj_uni,nio_uni)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpuniopt_uni
                       ! Formal parameter of unique index of uniopt_uni

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

! Output variables

      integer, intent(out) :: ni_uni
                       ! Model dimension of unite in x direction

      integer, intent(out) :: nj_uni
                       ! Model dimension of unite in y direction

      integer, intent(out) :: nio_uni
                       ! Maximum unit number of unite

! Internal shared variables

      integer uniopt_uni
                       ! Option for uniting process

      integer stat     ! Runtime status

      integer cstat    ! Runtime status at current allocate statement

! Internal private variables

      integer iio      ! Index of unit numbers table

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpuniopt_uni,uniopt_uni)

! -----

! Get the model dimension of unite.

      if(abs(mod(uniopt_uni,10)).le.4) then

        ni_uni=(ni-3)*nisub+3
        nj_uni=(nj-3)*njsub+3

      else if(abs(mod(uniopt_uni,10)).ge.5) then

        if(nio.le.nsub) then

          if(nio.le.nisub) then
            ni_uni=ni
            nj_uni=nj

          else
            ni_uni=(ni-3)*nisub+3
            nj_uni=nj

          end if

        else
          ni_uni=(ni-3)*nisub+3
          nj_uni=(nj-3)*njsub+3

        end if

      end if

      nio_uni=nio-1

! -----

!! Allocate the array for unite.

! Perform allocate.

      stat=0

      allocate(tmp1(1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp2(1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp3(1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp4(1:(nj-3)*njgrp*njsub),stat=cstat)

      stat=stat+abs(cstat)

      allocate(iodmp(0:nio_uni),stat=cstat)

      stat=stat+abs(cstat)

      allocate(var(2:ni_uni-2,2:nj_uni-2),stat=cstat)

      stat=stat+abs(cstat)

! -----

! If error occured, call the procedure destroy.

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('allocuni',8,'cont',5,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('allocuni',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

!! -----

! Fill in all array for the program unite with 0.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(k)

      do k=1,nk
        tmp1(k)=0.e0
        tmp2(k)=0.e0
        tmp3(k)=0.e0
      end do

!$omp end do

!$omp do schedule(runtime) private(j)

      do j=1,(nj-3)*njgrp*njsub
        tmp4(j)=0.e0
      end do

!$omp end do

!$omp do schedule(runtime) private(iio)

      do iio=0,nio_uni
        iodmp(iio)=0
      end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

      do j=2,nj_uni-2
      do i=2,ni_uni-2
        var(i,j)=0.e0
      end do
      end do

!$omp end do

!$omp end parallel

! -----

      end subroutine s_allocuni

!-----7--------------------------------------------------------------7--

      end module m_allocuni
