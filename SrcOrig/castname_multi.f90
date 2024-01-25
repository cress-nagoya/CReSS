!***********************************************************************
      module m_castname
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/25
!     Modification: 1999/04/06, 2000/01/17, 2001/08/07, 2002/04/09,
!                   2003/05/19, 2003/07/15, 2003/11/05, 2004/01/09,
!                   2006/12/04, 2007/01/20, 2008/05/02, 2008/07/25,
!                   2008/08/25, 2008/10/10, 2009/01/30, 2009/02/27,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     broadcast the namelist variables from the root processor element
!     to the others.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commpi
      use m_comname
      use m_defmpi

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: castname, s_castname

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface castname

        module procedure s_castname

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
      subroutine s_castname
!***********************************************************************

! Internal shared variables

      integer n108     ! Sended character variable byte
      integer nin2     ! Sended integer variable byte

      integer ierr     ! Error descriptor

!-----7--------------------------------------------------------------7--

! Broadcast the namelist variables.

      n108=108*ncn
      nin2=nin+2

      call mpi_bcast(cname,n108,mpi_character,root,mpi_comm_cress,ierr)

      call mpi_bcast(iname,nin2,mpi_integer,root,mpi_comm_cress,ierr)

      call mpi_bcast(rname,nrn,mpi_real,root,mpi_comm_cress,ierr)

! -----

      end subroutine s_castname

!-----7--------------------------------------------------------------7--

      end module m_castname
