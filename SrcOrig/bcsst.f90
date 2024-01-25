!***********************************************************************
      module m_bcsst
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/11/10
!     Modification: 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the boundary conditions for sea surface temperature data.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bc2d
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

      public :: bcsst, s_bcsst

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface bcsst

        module procedure s_bcsst

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
      subroutine s_bcsst(fpwbc,fpebc,fpexbopt,ni,nj,sst)
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

! Input and output variable

      real, intent(inout) :: sst(0:ni+1,0:nj+1)
                       ! Sea surface temperature

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

! Set the lateral boundary conditions.

      if(exbopt.eq.0) then

        call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,sst,1,1,sbuf)

        call s_shiftsx(idwbc,idebc,'bnd',nj,1,1,sbuf,rbuf)

        call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,sst,1,1,rbuf)

        call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,1,sst,1,1,sbuf)

        call s_shiftsy(idsbc,idnbc,'bnd',ni,1,1,sbuf,rbuf)

        call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,1,sst,1,1,rbuf)

        call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,sst,1,1,sbuf)

        call s_shiftgx(idwbc,idebc,'bnd',nj,1,1,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,sst,1,1,rbuf)

        call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,1,sst,1,1,sbuf)

        call s_shiftgy(idsbc,idnbc,'bnd',ni,1,1,sbuf,rbuf)

        call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,1,sst,1,1,rbuf)

        call s_bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,1,sst)

        call bc2d(idwbc,idebc,idsbc,idnbc,ni,nj,sst)

      end if

! -----

! Set the periodic boundary conditions in x direction.

      if(exbopt.ge.1.and.abs(wbc).eq.1.and.abs(ebc).eq.1) then

        call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,sst,1,1,sbuf)

        call s_shiftsx(idwbc,idebc,'bnd',nj,1,1,sbuf,rbuf)

        call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,sst,1,1,rbuf)

        call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,sst,1,1,sbuf)

        call s_shiftgx(idwbc,idebc,'bnd',nj,1,1,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,sst,1,1,rbuf)

        call s_bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,1,sst)

      end if

! -----

      end subroutine s_bcsst

!-----7--------------------------------------------------------------7--

      end module m_bcsst
