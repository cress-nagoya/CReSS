!***********************************************************************
      module m_shiftgx
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/12/04
!     Modification: 2007/01/20, 2008/05/02, 2008/07/25, 2008/08/25,
!                   2009/01/05, 2009/02/27, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     exchange the value in x direction between group domain.

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

      public :: shiftgx, s_shiftgx

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface shiftgx

        module procedure s_shiftgx

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
      subroutine s_shiftgx(fpwbc,fpebc,fproc,nj,kmax,nb,sbufx,rbufx)
!***********************************************************************

! Input variables

      character(len=3), intent(in) :: fproc
                       ! Control flag of processing type

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: kmax
                       ! Maximum array index in z direction

      integer, intent(in) :: nb
                       ! Number of exchanging variables

      real, intent(in) :: sbufx(0:nj+1,1:kmax,1:2*nb)
                       ! Sending buffer in x direction

! Output variable

      real, intent(out) :: rbufx(0:nj+1,1:kmax,1:2*nb)
                       ! Receiving buffer in x direction

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundaty conditions

      integer dstw     ! West sending distnation
      integer dste     ! East sending distnation

      integer srcw     ! West receiving source
      integer srce     ! East receiving source

      integer statsw   ! West sending request status
      integer statse   ! East sending request status

      integer statrw   ! West receiving request status
      integer statre   ! East receiving request status

      integer siz      ! Sending and receiving buffer size

      integer ierr     ! Error descriptor

      integer stat(mpi_status_size)
                       ! Runtime status table

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)

! -----

!! Exchange the value in x direction.

      if(nigrp.ge.2.and.(wbc.ne.-1.and.ebc.ne.-1)) then

        if(fproc(1:3).eq.'all'.or.(wbc.eq.1.and.ebc.eq.1)) then

! Set the buffer size.

          siz=(nj+2)*kmax*nb

! -----

! Incliment the message tag.

          tag=tag+1

! -----

! Set the processor element number for sending distination and
! receiving source.

          if(fproc(1:3).eq.'bnd') then

            dstw=dstw_grp_bnd
            srcw=srcw_grp_bnd

            dste=dste_grp_bnd
            srce=srce_grp_bnd

          else

            dstw=dstw_grp
            srcw=srcw_grp

            dste=dste_grp
            srce=srce_grp

          end if

! -----

! Call the sending and receiving MPI function.

          call mpi_isend(sbufx(0,1,1),siz,mpi_real,                     &
     &                   dstw,tag,mpi_comm_cress,statsw,ierr)

          call mpi_isend(sbufx(0,1,nb+1),siz,mpi_real,                  &
     &                   dste,tag,mpi_comm_cress,statse,ierr)

          call mpi_irecv(rbufx(0,1,nb+1),siz,mpi_real,                  &
     &                   srcw,tag,mpi_comm_cress,statrw,ierr)

          call mpi_irecv(rbufx(0,1,1),siz,mpi_real,                     &
     &                   srce,tag,mpi_comm_cress,statre,ierr)

! -----

! Call the waiting MPI function.

          call mpi_wait(statsw,stat,ierr)
          call mpi_wait(statse,stat,ierr)
          call mpi_wait(statrw,stat,ierr)
          call mpi_wait(statre,stat,ierr)

! -----

        end if

      end if

!! -----

      end subroutine s_shiftgx

!-----7--------------------------------------------------------------7--

      end module m_shiftgx
