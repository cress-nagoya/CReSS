!***********************************************************************
      module m_bcasl
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/09/22
!     Modification: 2011/11/10, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the boundary conditions for aerosol data.

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

      public :: bcasl, s_bcasl

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface bcasl

        module procedure s_bcasl

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
      subroutine s_bcasl(fpwbc,fpebc,fpexbopt,ni,nj,nk,nqa,qagpv)
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

      integer, intent(in) :: nqa(0:4)
                       ! Number of types of aerosol

! Input and output variable

      real, intent(inout) :: qagpv(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Mixing ratio of aerosol data at marked time

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions

      integer exbopt   ! Option for external boundary forcing

      integer n        ! Array index in 4th direction

      integer ib       ! Exchanging variables number
      integer nb       ! Number of exchanging variables

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpexbopt,exbopt)

! -----

! Count the number of variables.

      if(exbopt.eq.0                                                    &
     &  .or.(exbopt.ge.1.and.abs(wbc).eq.1.and.abs(ebc).eq.1)) then

        nb=nqa(0)

      end if

! -----

!!! Set the lateral boundary conditions.

      if(exbopt.eq.0) then

!! Exchange the value between sub domain.

! In x direction.

        ib=0

        do n=1,nqa(0)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,            &
     &                    qagpv(0,0,1,n),ib,nb,sbuf)

        end do

        call s_shiftsx(idwbc,idebc,'bnd',nj,nk,nb,sbuf,rbuf)

        ib=0

        do n=1,nqa(0)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,            &
     &                    qagpv(0,0,1,n),ib,nb,rbuf)

        end do

! -----

! In y direction.

        ib=0

        do n=1,nqa(0)

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,            &
     &                    qagpv(0,0,1,n),ib,nb,sbuf)

        end do

        call s_shiftsy(idsbc,idnbc,'bnd',ni,nk,nb,sbuf,rbuf)

        ib=0

        do n=1,nqa(0)

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,            &
     &                    qagpv(0,0,1,n),ib,nb,rbuf)

        end do

! -----

!! -----

!! Exchange the value between group domain.

! In x direction.

        ib=0

        do n=1,nqa(0)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,            &
     &                    qagpv(0,0,1,n),ib,nb,sbuf)

        end do

        call s_shiftgx(idwbc,idebc,'bnd',nj,nk,nb,sbuf,rbuf)

        ib=0

        do n=1,nqa(0)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,            &
     &                    qagpv(0,0,1,n),ib,nb,rbuf)

        end do

! -----

! In y direction.

        ib=0

        do n=1,nqa(0)

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,            &
     &                    qagpv(0,0,1,n),ib,nb,sbuf)

        end do

        call s_shiftgy(idsbc,idnbc,'bnd',ni,nk,nb,sbuf,rbuf)

        ib=0

        do n=1,nqa(0)

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,            &
     &                    qagpv(0,0,1,n),ib,nb,rbuf)

        end do

! -----

!! -----

! Set the periodic boundary conditions.

        do n=1,nqa(0)

          call s_bcycle(idwbc,idebc,idsbc,idnbc,                        &
     &              2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,qagpv(0,0,1,n))

        end do

! -----

      end if

!!! -----

!! Set the periodic boundary conditions in x direction.

      if(exbopt.ge.1.and.abs(wbc).eq.1.and.abs(ebc).eq.1) then

! Exchange the value between sub domain.

        ib=0

        do n=1,nqa(0)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,            &
     &                    qagpv(0,0,1,n),ib,nb,sbuf)

        end do

        call s_shiftsx(idwbc,idebc,'bnd',nj,nk,nb,sbuf,rbuf)

        ib=0

        do n=1,nqa(0)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,            &
     &                    qagpv(0,0,1,n),ib,nb,rbuf)

        end do

! -----

! Exchange the value between group domain.

        ib=0

        do n=1,nqa(0)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,            &
     &                    qagpv(0,0,1,n),ib,nb,sbuf)

        end do

        call s_shiftgx(idwbc,idebc,'bnd',nj,nk,nb,sbuf,rbuf)

        ib=0

        do n=1,nqa(0)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,            &
     &                    qagpv(0,0,1,n),ib,nb,rbuf)

        end do

! -----

! Set the periodic boundary conditions.

        do n=1,nqa(0)

          call s_bcyclex(idwbc,idebc,                                   &
     &                   2,1,ni-2,ni-1,ni,nj,nk,qagpv(0,0,1,n))

        end do

! -----

      end if

!! -----

      end subroutine s_bcasl

!-----7--------------------------------------------------------------7--

      end module m_bcasl
