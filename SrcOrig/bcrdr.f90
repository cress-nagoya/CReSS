!***********************************************************************
      module m_bcrdr
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2005/02/10
!     Modification: 2006/09/21, 2006/09/30, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/07/30, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/02/13, 2013/03/27

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
      use m_getcname
      use m_getiname
      use m_inichar
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

      public :: bcrdr, s_bcrdr

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface bcrdr

        module procedure s_bcrdr

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
      subroutine s_bcrdr(fprdrvar,fpwbc,fpebc,fpexbopt,fpngropt,        &
     &                   ni,nj,nk,urdr,vrdr,wrdr,qprdr)
!***********************************************************************

! Input variables

      integer, intent(in) :: fprdrvar
                       ! Formal parameter of unique index of rdrvar

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpexbopt
                       ! Formal parameter of unique index of exbopt

      integer, intent(in) :: fpngropt
                       ! Formal parameter of unique index of ngropt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

! Input and output variables

      real, intent(inout) :: urdr(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity of radar data
                       ! at marked time

      real, intent(inout) :: vrdr(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity of radar data
                       ! at marked time

      real, intent(inout) :: wrdr(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity of radar data
                       ! at marked time

      real, intent(inout) :: qprdr(0:ni+1,0:nj+1,1:nk)
                       ! Precipitation mixing ratio of radar data
                       ! at marked time

! Internal shared variables

      character(len=108) rdrvar
                       ! Control flag of input radar data variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions

      integer exbopt   ! Option for external boundary forcing
      integer ngropt   ! Option for analysis nudging to radar

      integer ib       ! Exchanging variables number
      integer nb       ! Number of exchanging variables

      integer ni_sub   ! Substitute for ni
      integer nj_sub   ! Substitute for nj

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(rdrvar)

! -----

! Get the required namelist variables.

      call getcname(fprdrvar,rdrvar)
      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpexbopt,exbopt)
      call getiname(fpngropt,ngropt)

! -----

! Set the substituted variables.

      ni_sub=ni
      nj_sub=nj

! -----

! Count the number of variables.

      if(exbopt.eq.0                                                    &
     &  .or.(exbopt.ge.1.and.abs(wbc).eq.1.and.abs(ebc).eq.1)) then

        nb=0

        if(rdrvar(1:1).eq.'o') then
          nb=nb+1
        end if

        if(rdrvar(2:2).eq.'o') then
          nb=nb+1
        end if

        if(rdrvar(3:3).eq.'o') then
          nb=nb+1
        end if

        if(rdrvar(4:4).eq.'o'.and.ngropt.ge.2) then
          nb=nb+1
        end if

      end if

! -----

!!! Set the lateral boundary conditions.

      if(exbopt.eq.0) then

!! Exchange the value horizontally between sub domain.

! In x direction.

        ib=0

        if(rdrvar(1:1).eq.'o') then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,urdr,       &
     &                    ib,nb,sbuf)

        end if

        if(rdrvar(2:2).eq.'o') then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,vrdr,       &
     &                    ib,nb,sbuf)

        end if

        if(rdrvar(3:3).eq.'o') then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,wrdr,       &
     &                    ib,nb,sbuf)

        end if

        if(rdrvar(4:4).eq.'o'.and.ngropt.ge.2) then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,qprdr,      &
     &                    ib,nb,sbuf)

        end if

        call s_shiftsx(idwbc,idebc,'bnd',nj,nk,nb,sbuf,rbuf)

        ib=0

        if(rdrvar(1:1).eq.'o') then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,urdr,     &
     &                    ib,nb,rbuf)

        end if

        if(rdrvar(2:2).eq.'o') then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,vrdr,       &
     &                    ib,nb,rbuf)

        end if

        if(rdrvar(3:3).eq.'o') then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,wrdr,       &
     &                    ib,nb,rbuf)

        end if

        if(rdrvar(4:4).eq.'o'.and.ngropt.ge.2) then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,qprdr,      &
     &                    ib,nb,rbuf)

        end if

! -----

! In y direction.

        ib=0

        if(rdrvar(1:1).eq.'o') then

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,urdr,       &
     &                    ib,nb,sbuf)

        end if

        if(rdrvar(2:2).eq.'o') then

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',3,nj-2,ni,nj,nk,vrdr,       &
     &                    ib,nb,sbuf)

        end if

        if(rdrvar(3:3).eq.'o') then

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,wrdr,       &
     &                    ib,nb,sbuf)

        end if

        if(rdrvar(4:4).eq.'o'.and.ngropt.ge.2) then

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,qprdr,      &
     &                    ib,nb,sbuf)

        end if

        call s_shiftsy(idsbc,idnbc,'bnd',ni,nk,nb,sbuf,rbuf)

        ib=0

        if(rdrvar(1:1).eq.'o') then

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,urdr,       &
     &                    ib,nb,rbuf)

        end if

        if(rdrvar(2:2).eq.'o') then

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',1,nj_sub,ni,nj,nk,vrdr,     &
     &                    ib,nb,rbuf)

        end if

        if(rdrvar(3:3).eq.'o') then

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,wrdr,       &
     &                    ib,nb,rbuf)

        end if

        if(rdrvar(4:4).eq.'o'.and.ngropt.ge.2) then

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,qprdr,      &
     &                    ib,nb,rbuf)

        end if

! -----

!! -----

!! Exchange the value horizontally betweein group domain.

! In x direction.

        ib=0

        if(rdrvar(1:1).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,urdr,       &
     &                    ib,nb,sbuf)

        end if

        if(rdrvar(2:2).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,vrdr,       &
     &                    ib,nb,sbuf)

        end if

        if(rdrvar(3:3).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,wrdr,       &
     &                    ib,nb,sbuf)

        end if

        if(rdrvar(4:4).eq.'o'.and.ngropt.ge.2) then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,qprdr,      &
     &                    ib,nb,sbuf)

        end if

        call s_shiftgx(idwbc,idebc,'bnd',nj,nk,nb,sbuf,rbuf)

        ib=0

        if(rdrvar(1:1).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,urdr,     &
     &                    ib,nb,rbuf)

        end if

        if(rdrvar(2:2).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,vrdr,       &
     &                    ib,nb,rbuf)

        end if

        if(rdrvar(3:3).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,wrdr,       &
     &                    ib,nb,rbuf)

        end if

        if(rdrvar(4:4).eq.'o'.and.ngropt.ge.2) then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,qprdr,      &
     &                    ib,nb,rbuf)

        end if

! -----

! In y direction.

        ib=0

        if(rdrvar(1:1).eq.'o') then

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,urdr,       &
     &                    ib,nb,sbuf)

        end if

        if(rdrvar(2:2).eq.'o') then

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',3,nj-2,ni,nj,nk,vrdr,       &
     &                    ib,nb,sbuf)

        end if

        if(rdrvar(3:3).eq.'o') then

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,wrdr,       &
     &                    ib,nb,sbuf)

        end if

        if(rdrvar(4:4).eq.'o'.and.ngropt.ge.2) then

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,qprdr,      &
     &                    ib,nb,sbuf)

        end if

        call s_shiftgy(idsbc,idnbc,'bnd',ni,nk,nb,sbuf,rbuf)

        ib=0

        if(rdrvar(1:1).eq.'o') then

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,urdr,       &
     &                    ib,nb,rbuf)

        end if

        if(rdrvar(2:2).eq.'o') then

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',1,nj_sub,ni,nj,nk,vrdr,     &
     &                    ib,nb,rbuf)

        end if

        if(rdrvar(3:3).eq.'o') then

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,wrdr,       &
     &                    ib,nb,rbuf)

        end if

        if(rdrvar(4:4).eq.'o'.and.ngropt.ge.2) then

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,qprdr,      &
     &                    ib,nb,rbuf)

        end if

! -----

!! -----

! Set the periodic boundary conditions.

        if(rdrvar(1:1).eq.'o') then

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                3,1,ni-2,ni_sub,2,1,nj-2,nj-1,ni,nj,nk,urdr)

        end if

        if(rdrvar(2:2).eq.'o') then

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,3,1,nj-2,nj_sub,ni,nj,nk,vrdr)

        end if

        if(rdrvar(3:3).eq.'o') then

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,wrdr)

        end if

        if(rdrvar(4:4).eq.'o'.and.ngropt.ge.2) then

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,qprdr)

        end if

! -----

      end if

!!! -----

!! Set the periodic boundary conditions in x direction.

      if(exbopt.ge.1.and.abs(wbc).eq.1.and.abs(ebc).eq.1) then

! Exchange the value horizontally betweein sub domain.

        ib=0

        if(rdrvar(1:1).eq.'o') then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,urdr,       &
     &                    ib,nb,sbuf)

        end if

        if(rdrvar(2:2).eq.'o') then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,vrdr,       &
     &                    ib,nb,sbuf)

        end if

        if(rdrvar(3:3).eq.'o') then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,wrdr,       &
     &                    ib,nb,sbuf)

        end if

        if(rdrvar(4:4).eq.'o'.and.ngropt.ge.2) then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,qprdr,      &
     &                    ib,nb,sbuf)

        end if

        call s_shiftsx(idwbc,idebc,'bnd',nj,nk,nb,sbuf,rbuf)

        ib=0

        if(rdrvar(1:1).eq.'o') then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,urdr,     &
     &                    ib,nb,rbuf)

        end if

        if(rdrvar(2:2).eq.'o') then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,vrdr,       &
     &                    ib,nb,rbuf)

        end if

        if(rdrvar(3:3).eq.'o') then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,wrdr,       &
     &                    ib,nb,rbuf)

        end if

        if(rdrvar(4:4).eq.'o'.and.ngropt.ge.2) then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,qprdr,      &
     &                    ib,nb,rbuf)

        end if

! -----

! Exchange the value horizontally betweein group domain.

        ib=0

        if(rdrvar(1:1).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,urdr,       &
     &                    ib,nb,sbuf)

        end if

        if(rdrvar(2:2).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,vrdr,       &
     &                    ib,nb,sbuf)

        end if

        if(rdrvar(3:3).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,wrdr,       &
     &                    ib,nb,sbuf)

        end if

        if(rdrvar(4:4).eq.'o'.and.ngropt.ge.2) then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,qprdr,      &
     &                    ib,nb,sbuf)

        end if

        call s_shiftgx(idwbc,idebc,'bnd',nj,nk,nb,sbuf,rbuf)

        ib=0

        if(rdrvar(1:1).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,urdr,     &
     &                    ib,nb,rbuf)

        end if

        if(rdrvar(2:2).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,vrdr,       &
     &                    ib,nb,rbuf)

        end if

        if(rdrvar(3:3).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,wrdr,       &
     &                    ib,nb,rbuf)

        end if

        if(rdrvar(4:4).eq.'o'.and.ngropt.ge.2) then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,qprdr,      &
     &                    ib,nb,rbuf)

        end if

! -----

! Set the periodic boundary conditions in x direction.

        if(rdrvar(1:1).eq.'o') then

          call bcyclex(idwbc,idebc,3,1,ni-2,ni_sub,ni,nj,nk,urdr)

        end if

        if(rdrvar(2:2).eq.'o') then

          call bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,nk,vrdr)

        end if

        if(rdrvar(3:3).eq.'o') then

          call bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,nk,wrdr)

        end if

        if(rdrvar(4:4).eq.'o'.and.ngropt.ge.2) then

          call bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,nk,qprdr)

        end if

! -----

      end if

!! -----

      end subroutine s_bcrdr

!-----7--------------------------------------------------------------7--

      end module m_bcrdr
