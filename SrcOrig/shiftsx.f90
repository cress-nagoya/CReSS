!***********************************************************************
      module m_shiftsx
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/25
!     Modification: 1999/03/25, 1999/04/06, 1999/07/05, 1999/08/18,
!                   1999/08/23, 1999/10/07, 1999/10/12, 1999/11/01,
!                   1999/11/19, 2000/01/17, 2000/04/18, 2000/07/05,
!                   2001/05/29, 2002/04/02, 2002/06/06, 2003/04/30,
!                   2003/05/19, 2003/11/05, 2004/08/20, 2006/12/04,
!                   2007/01/05, 2007/01/20, 2008/05/02, 2008/07/25,
!                   2008/08/25, 2009/01/05, 2009/02/27, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     exchange the value in x direction between sub domain.

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

      public :: shiftsx, s_shiftsx

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface shiftsx

        module procedure s_shiftsx

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
      subroutine s_shiftsx(fpwbc,fpebc,fproc,nj,kmax,nb,sbufx,rbufx)
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

      if(nisub.ge.2) then

        if(fproc(1:3).eq.'all'.or.(wbc.eq.-1.and.ebc.eq.-1)) then

! Set the buffer size.

          siz=(nj+2)*kmax*nb

! -----

! Incliment the message tag.

          tag=tag+1

! -----

! Set the processor element number for sending distination and
! receiving source.

          if(fproc(1:3).eq.'bnd') then

            dstw=dstw_sub_bnd
            srcw=srcw_sub_bnd

            dste=dste_sub_bnd
            srce=srce_sub_bnd

          else

            dstw=dstw_sub
            srcw=srcw_sub

            dste=dste_sub
            srce=srce_sub

          end if

! -----

! Call the sending and receiving MPI function.

          call mpi_isend(sbufx(0,1,1),siz,mpi_real,                     &
     &                   dstw,tag,mpi_comm_cress_sub,statsw,ierr)

          call mpi_isend(sbufx(0,1,nb+1),siz,mpi_real,                  &
     &                   dste,tag,mpi_comm_cress_sub,statse,ierr)

          call mpi_irecv(rbufx(0,1,nb+1),siz,mpi_real,                  &
     &                   srcw,tag,mpi_comm_cress_sub,statrw,ierr)

          call mpi_irecv(rbufx(0,1,1),siz,mpi_real,                     &
     &                   srce,tag,mpi_comm_cress_sub,statre,ierr)

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

      end subroutine s_shiftsx

!-----7--------------------------------------------------------------7--

      end module m_shiftsx
