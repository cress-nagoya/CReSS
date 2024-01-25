!***********************************************************************
      module m_allocgrp
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/12/04
!     Modification: 2007/01/05, 2007/01/20, 2007/01/31, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/01/05, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     allocate the group domain arrangement table.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_comgrp
      use m_commpi
      use m_cpondpe
      use m_destroy
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: allocgrp, s_allocgrp

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface allocgrp

        module procedure s_allocgrp

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
      subroutine s_allocgrp(fpwbc,fpebc,fpsbc,fpnbc)
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

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions
      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

      integer nigrp1   ! nigrp + 1
      integer njgrp1   ! njgrp + 1

      integer igc      ! Array index in group domain in x direction
      integer jgc      ! Array index in group domain in y direction

      integer stat     ! Runtime status

      integer cstat    ! Runtime status at current allocate statement

! Internal private variables

      integer ijsc     ! Index of serial number table

      integer igc_sub  ! Substitute for igc
      integer jgc_sub  ! Substitute for jgc

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpsbc,sbc)
      call getiname(fpnbc,nbc)

! -----

! Set the common used variables.

      nigrp1=nigrp+1
      njgrp1=njgrp+1

! -----

!! Allocate the arrangement table.

! Perfom the allocate.

      stat=0

      allocate(chrgrp(1:nigrp,1:njgrp),stat=cstat)

      stat=stat+abs(cstat)

      allocate(grpxy(0:nigrp+1,0:njgrp+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(xgrp(0:nsrl-1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(ygrp(0:nsrl-1),stat=cstat)

      stat=stat+abs(cstat)

! -----

! If error occured, call the procedure destroy.

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('allocgrp',8,'cont',5,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('allocgrp',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

!! -----

!! Initialize the table and set boundary conditions.

! Initiaze the character table.

      do jgc=1,njgrp
      do igc=1,nigrp

        write(chrgrp(igc,jgc)(1:1),'(a1)') '-'

      end do
      end do

! -----

! Initialize the other table and set boundary conditions.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(igc_sub,jgc_sub)

      do jgc_sub=0,njgrp+1
      do igc_sub=0,nigrp+1
        grpxy(igc_sub,jgc_sub)=-1
      end do
      end do

!$omp end do

      if(abs(wbc).eq.1.and.abs(ebc).eq.1) then

!$omp do schedule(runtime) private(jgc_sub)

        do jgc_sub=0,njgrp+1
          grpxy(0,jgc_sub)=nsrl
          grpxy(nigrp1,jgc_sub)=nsrl
        end do

!$omp end do

      end if

      if(abs(sbc).eq.1.and.abs(nbc).eq.1) then

!$omp do schedule(runtime) private(igc_sub)

        do igc_sub=0,nigrp+1
          grpxy(igc_sub,0)=nsrl
          grpxy(igc_sub,njgrp1)=nsrl
        end do

!$omp end do

      end if

!$omp do schedule(runtime) private(ijsc)

      do ijsc=0,nsrl-1
        xgrp(ijsc)=-1
        ygrp(ijsc)=-1
      end do

!$omp end do

!$omp end parallel

! -----

!! -----

      end subroutine s_allocgrp

!-----7--------------------------------------------------------------7--

      end module m_allocgrp
