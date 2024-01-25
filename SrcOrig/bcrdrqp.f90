!***********************************************************************
      module m_bcrdrqp
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/07/30
!     Modification: 2008/05/02, 2008/08/25, 2009/02/27, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the boundary conditions for radar data.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bcycle
      use m_bcyclex
      use m_combuf
      use m_comindx
      use m_getbufgx
      use m_getbufgy
      use m_getbufsx
      use m_getbufsy
      use m_getiname
      use m_putbufgx
      use m_putbufgy
      use m_putbufsx
      use m_putbufsy
      use m_shiftgx
      use m_shiftgy
      use m_shiftsx
      use m_shiftsy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: bcrdrqp, s_bcrdrqp

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface bcrdrqp

        module procedure s_bcrdrqp

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
      subroutine s_bcrdrqp(fpwbc,fpebc,fpexbopt,ni,nj,nk,qprdr)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpexbopt
                       ! Formal parameter of unique index of exbopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

! Input and output variable

      real, intent(inout) :: qprdr(0:ni+1,0:nj+1,1:nk)
                       ! Precipitation mixing ratio of radar data
                       ! at marked time

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions

      integer exbopt   ! Option for external boundary forcing

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpexbopt,exbopt)

! -----

!! Set the lateral boundary conditions.

      if(exbopt.eq.0) then

! Exchange the value horizontally between sub domain.

        call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,qprdr,1,1,    &
     &                  sbuf)

        call s_shiftsx(idwbc,idebc,'bnd',nj,nk,1,sbuf,rbuf)

        call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,qprdr,1,1,    &
     &                  rbuf)

        call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,qprdr,1,1,    &
     &                  sbuf)

        call s_shiftsy(idsbc,idnbc,'bnd',ni,nk,1,sbuf,rbuf)

        call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,qprdr,1,1,    &
     &                  rbuf)

! -----

! Exchange the value horizontally betweein group domain.

        call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,qprdr,1,1,    &
     &                  sbuf)

        call s_shiftgx(idwbc,idebc,'bnd',nj,nk,1,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,qprdr,1,1,    &
     &                  rbuf)

        call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,qprdr,1,1,    &
     &                  sbuf)

        call s_shiftgy(idsbc,idnbc,'bnd',ni,nk,1,sbuf,rbuf)

        call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,qprdr,1,1,    &
     &                  rbuf)

! -----

! Set the periodic boundary conditions.

        call bcycle(idwbc,idebc,idsbc,idnbc,                            &
     &              2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,qprdr)

! -----

      end if

!! -----

!! Set the periodic boundary conditions in x direction.

      if(exbopt.ge.1.and.abs(wbc).eq.1.and.abs(ebc).eq.1) then

! Exchange the value horizontally betweein sub domain.

        call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,qprdr,1,1,    &
     &                  sbuf)

        call s_shiftsx(idwbc,idebc,'bnd',nj,nk,1,sbuf,rbuf)

        call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,qprdr,1,1,    &
     &                  rbuf)

! -----

! Exchange the value horizontally betweein group domain.

        call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,qprdr,1,1,    &
     &                  sbuf)

        call s_shiftgx(idwbc,idebc,'bnd',nj,nk,1,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,qprdr,1,1,    &
     &                  rbuf)

! -----

! Set the periodic boundary conditions in x direction.

        call bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,nk,qprdr)

! -----

      end if

!! -----

      end subroutine s_bcrdrqp

!-----7--------------------------------------------------------------7--

      end module m_bcrdrqp
