!***********************************************************************
      module m_inimpi
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/25
!     Modification: 1999/04/06, 1999/06/21, 1999/08/23, 1999/09/16,
!                   1999/09/30, 1999/10/12, 2000/01/05, 2000/01/17,
!                   2000/02/07, 2000/03/23, 2000/04/18, 2001/02/13,
!                   2001/05/29, 2001/06/06, 2002/04/02, 2002/06/18,
!                   2002/07/31, 2002/08/15, 2002/09/09, 2003/04/30,
!                   2003/05/19, 2003/11/05, 2004/05/31, 2004/08/20,
!                   2006/09/21, 2006/12/04, 2007/01/05, 2007/01/20,
!                   2007/04/11, 2008/05/02, 2008/07/25, 2008/08/25,
!                   2008/10/10, 2009/01/05, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     call the function to initialize the MPI processes and set the
!     parameters of parallelizing.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkstd
      use m_commpi
      use m_defmpi
      use m_outstd01

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: inimpi, s_inimpi

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface inimpi

        module procedure s_inimpi

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
      subroutine s_inimpi(pname,ncpn)
!***********************************************************************

! Input variables

      character(len=8), intent(in) :: pname
                       ! Running program name

      integer, intent(in) :: ncpn
                       ! Number of character of pname

! Internal shared variable

      integer ierr     ! Error descriptor

!-----7--------------------------------------------------------------7--

! Call the initial MPI function.

      call mpi_init(ierr)

! -----

! Initialize the parameters of communicators.

      mype=0

      npe=0

      mpi_comm_cress=0

      mpi_comm_cress_sub=0

! -----

! Get the new communicator.

      call mpi_comm_dup(mpi_comm_world,mpi_comm_cress,ierr)

! -----

! Get the number of entried processor elements.

      call mpi_comm_size(mpi_comm_cress,npe,ierr)

! -----

! Get my processor element number.

      call mpi_comm_rank(mpi_comm_cress,mype,ierr)

! -----

! Initialize the parameters of parallelizing.

      root=0

      tag=0

      if(npe.eq.1) then
        write(fpara(1:6),'(a6)') 'single'
      else
        write(fpara(1:6),'(a6)') 'multi '
      end if

      mysub=0

      nsub=0

      nisub=0
      njsub=0

      isub=0
      jsub=0

      mygrp=0

      ngrp=0

      nigrp=0
      njgrp=0

      igrp=0
      jgrp=0

      mysrl=0

      nsrl=0

      myred=0

      nred=0

      nired=0
      njred=0

      ired=0
      jred=0

      iwred=0
      jsred=0

      ebw=0
      ebe=0
      ebs=0
      ebn=0

      ebsw=0
      ebse=0
      ebnw=0
      ebne=0

      dstw_sub=0
      dste_sub=0
      dsts_sub=0
      dstn_sub=0

      srcw_sub=0
      srce_sub=0
      srcs_sub=0
      srcn_sub=0

      dstw_grp=0
      dste_grp=0
      dsts_grp=0
      dstn_grp=0

      srcw_grp=0
      srce_grp=0
      srcs_grp=0
      srcn_grp=0

      dstw_sub_bnd=0
      dste_sub_bnd=0
      dsts_sub_bnd=0
      dstn_sub_bnd=0

      srcw_sub_bnd=0
      srce_sub_bnd=0
      srcs_sub_bnd=0
      srcn_sub_bnd=0

      dstw_grp_bnd=0
      dste_grp_bnd=0
      dsts_grp_bnd=0
      dstn_grp_bnd=0

      srcw_grp_bnd=0
      srce_grp_bnd=0
      srcs_grp_bnd=0
      srcn_grp_bnd=0

      mysub_rst=0

      nsub_rst=0

      nisub_rst=0
      njsub_rst=0

      isub_rst=0
      jsub_rst=0

! -----

! Read in the fisrt message to standard i/o.

      if(mype.eq.root) then

        call outstd01(pname,ncpn)

      end if

      call chkstd(root)

! -----

      end subroutine s_inimpi

!-----7--------------------------------------------------------------7--

      end module m_inimpi
