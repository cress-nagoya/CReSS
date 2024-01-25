!***********************************************************************
      module m_shiftsy
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/25
!     Modification: 1999/03/25, 1999/04/06, 1999/07/05, 1999/08/18,
!                   1999/08/23, 1999/10/07, 1999/10/12, 1999/11/01,
!                   2000/01/17, 2000/04/18, 2000/07/05, 2001/05/29,
!                   2002/04/02, 2002/06/06, 2003/04/30, 2003/05/19,
!                   2003/11/05, 2004/08/20, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2008/05/02, 2008/07/25, 2008/08/25,
!                   2009/01/05, 2009/02/27, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     exchange the value in y direction between sub domain.

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

      public :: shiftsy, s_shiftsy

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface shiftsy

        module procedure s_shiftsy

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
      subroutine s_shiftsy(fpsbc,fpnbc,fproc,ni,kmax,nb,sbufy,rbufy)
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

      if(njsub.ge.2) then

        if(fproc(1:3).eq.'all'.or.(sbc.eq.-1.and.nbc.eq.-1)) then

! Set the buffer size.

          siz=(ni+2)*kmax*nb

! -----

! Incliment the message tag.

          tag=tag+1

! -----

! Set the processor element number for sending distination and
! receiving source.

          if(fproc(1:3).eq.'bnd') then

            dsts=dsts_sub_bnd
            srcs=srcs_sub_bnd

            dstn=dstn_sub_bnd
            srcn=srcn_sub_bnd

          else

            dsts=dsts_sub
            srcs=srcs_sub

            dstn=dstn_sub
            srcn=srcn_sub

          end if

! -----

! Call the sending and receiving MPI function.

          call mpi_isend(sbufy(0,1,1),siz,mpi_real,                     &
     &                   dsts,tag,mpi_comm_cress_sub,statss,ierr)

          call mpi_isend(sbufy(0,1,nb+1),siz,mpi_real,                  &
     &                   dstn,tag,mpi_comm_cress_sub,statsn,ierr)

          call mpi_irecv(rbufy(0,1,nb+1),siz,mpi_real,                  &
     &                   srcs,tag,mpi_comm_cress_sub,statrs,ierr)

          call mpi_irecv(rbufy(0,1,1),siz,mpi_real,                     &
     &                   srcn,tag,mpi_comm_cress_sub,statrn,ierr)

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

      end subroutine s_shiftsy

!-----7--------------------------------------------------------------7--

      end module m_shiftsy
