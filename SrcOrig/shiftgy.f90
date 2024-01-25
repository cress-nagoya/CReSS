!***********************************************************************
      module m_shiftgy
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/12/04
!     Modification: 2007/01/20, 2008/05/02, 2008/07/25, 2008/08/25,
!                   2009/01/05, 2009/02/27, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     exchange the value in y direction between group domain.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commpi
      use m_defmpi
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: shiftgy, s_shiftgy

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface shiftgy

        module procedure s_shiftgy

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
      subroutine s_shiftgy(fpsbc,fpnbc,fproc,ni,kmax,nb,sbufy,rbufy)
!***********************************************************************

! Input variables

      character(len=3), intent(in) :: fproc
                       ! Control flag of processing type

      integer, intent(in) :: fpsbc
                       ! Formal parameter of unique index of sbc

      integer, intent(in) :: fpnbc
                       ! Formal parameter of unique index of nbc

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: kmax
                       ! Maximum array index in z direction

      integer, intent(in) :: nb
                       ! Number of exchanging variables

      real, intent(in) :: sbufy(0:ni+1,1:kmax,1:2*nb)
                       ! Sending buffer in y direction

! Output variable

      real, intent(out) :: rbufy(0:ni+1,1:kmax,1:2*nb)
                       ! Receiving buffer in y direction

! Internal shared variables

      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

      integer dsts     ! South sending distnation
      integer dstn     ! North sending distnation

      integer srcs     ! South receiving source
      integer srcn     ! North receiving source

      integer statss   ! South sending request status
      integer statsn   ! North sending request status

      integer statrs   ! South receiving request status
      integer statrn   ! North receiving request status

      integer siz      ! Sending and receiving buffer size

      integer ierr     ! Error descriptor

      integer stat(mpi_status_size)
                       ! Runtime status table

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpsbc,sbc)
      call getiname(fpnbc,nbc)

! -----

!! Exchange the value in y direction.

      if(njgrp.ge.2.and.(sbc.ne.-1.and.nbc.ne.-1)) then

        if(fproc(1:3).eq.'all'.or.(sbc.eq.1.and.nbc.eq.1)) then

! Set the buffer size.

          siz=(ni+2)*kmax*nb

! -----

! Incliment the message tag.

          tag=tag+1

! -----

! Set the processor element number for sending distination and
! receiving source.

          if(fproc(1:3).eq.'bnd') then

            dsts=dsts_grp_bnd
            srcs=srcs_grp_bnd

            dstn=dstn_grp_bnd
            srcn=srcn_grp_bnd

          else

            dsts=dsts_grp
            srcs=srcs_grp

            dstn=dstn_grp
            srcn=srcn_grp

          end if

! -----

! Call the sending and receiving MPI function.

          call mpi_isend(sbufy(0,1,1),siz,mpi_real,                     &
     &                   dsts,tag,mpi_comm_cress,statss,ierr)

          call mpi_isend(sbufy(0,1,nb+1),siz,mpi_real,                  &
     &                   dstn,tag,mpi_comm_cress,statsn,ierr)

          call mpi_irecv(rbufy(0,1,nb+1),siz,mpi_real,                  &
     &                   srcs,tag,mpi_comm_cress,statrs,ierr)

          call mpi_irecv(rbufy(0,1,1),siz,mpi_real,                     &
     &                   srcn,tag,mpi_comm_cress,statrn,ierr)

! -----

! Call the waiting MPI function.

          call mpi_wait(statss,stat,ierr)
          call mpi_wait(statsn,stat,ierr)
          call mpi_wait(statrs,stat,ierr)
          call mpi_wait(statrn,stat,ierr)

! -----

        end if

      end if

!! -----

      end subroutine s_shiftgy

!-----7--------------------------------------------------------------7--

      end module m_shiftgy
