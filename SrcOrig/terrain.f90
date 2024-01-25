!***********************************************************************
      program terrain
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/04/06
!     Modification: 1999/05/10, 1999/06/07, 1999/06/28, 1999/09/30,
!                   1999/10/22, 1999/11/01, 2000/01/17, 2000/04/18,
!                   2000/06/01, 2001/02/13, 2001/03/13, 2001/10/18,
!                   2002/01/07, 2002/04/09, 2002/06/18, 2002/07/03,
!                   2002/07/15, 2002/09/09, 2002/12/02, 2003/02/05,
!                   2003/05/19, 2003/07/15, 2003/11/05, 2003/12/12,
!                   2004/01/09, 2004/08/01, 2005/02/10, 2006/01/10,
!                   2006/09/30, 2006/12/04, 2007/01/05, 2007/01/20,
!                   2007/07/30, 2008/04/17, 2008/05/02, 2008/08/25,
!                   2008/10/10, 2009/01/30, 2009/02/27, 2011/08/18,
!                   2011/09/22, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! This is the main program for the pre processor terrain.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_allocgrp
      use m_allociot
      use m_alloctrn
      use m_chkkind
      use m_comindx
      use m_comtrn
      use m_defdim
      use m_endmpi
      use m_inidef
      use m_inimpi
      use m_ininame
      use m_rdgrp
      use m_setdim
      use m_setnltrn
      use m_trndrv

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

!     none

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Initialize the parameters of parallelizing.

      call inimpi('terrain ',7)

! -----

! Check the kind of the Fortran variables.

      call chkkind

! -----

! Initialize the model dimension and namelist variables.

      call inidef

! -----

! Initialize the table to archive namelist variables.

      call ininame

! -----

! Set the namelist variables.

      call setnltrn

! -----

! Set the model dimension variables.

      call setdim(idcphopt,idhaiopt,idaslopt)

! -----

! Allocate the table of unit numbers.

      call allociot

! -----

! Allocate the array for terrain.

      call alloctrn(ni,nj,nid_trn,njd_trn)

! -----

! Allocate the group domain arrangement table.

      call allocgrp(idwbc,idebc,idsbc,idnbc)

! -----

! Readout and check the group domain arrangement.

      call rdgrp(idexprim,idcrsdir,idncexp,idnccrs)

! -----

! Interpolate the terrain data.

      call trndrv(ni,nj,ri,rj,ht,tmp1,tmp2,tmp3,nid_trn,njd_trn,htdat)

! -----

! Finalize the MPI processes.

      call endmpi(0)

! -----

!-----7--------------------------------------------------------------7--

! Internal procedure

!     none

!-----7--------------------------------------------------------------7--

      end program terrain
