!***********************************************************************
      program unite
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/08/10
!     Modification: 2001/02/13, 2001/03/13, 2001/05/29, 2001/06/29,
!                   2002/01/15, 2002/06/18, 2002/09/02, 2002/09/09,
!                   2002/12/02, 2003/02/05, 2003/05/19, 2003/07/15,
!                   2003/11/05, 2003/12/12, 2004/01/09, 2004/07/01,
!                   2004/08/01, 2004/08/20, 2005/01/14, 2005/02/10,
!                   2006/01/10, 2006/09/30, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/04/11, 2007/06/27, 2007/07/30,
!                   2008/01/11, 2008/04/17, 2008/05/02, 2008/08/25,
!                   2008/10/10, 2009/01/30, 2009/02/27, 2011/08/18,
!                   2011/09/22, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! This is the main program for the post processor unite.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_allocgrp
      use m_allociot
      use m_allocuni
      use m_chkkind
      use m_comindx
      use m_comuni
      use m_defdim
      use m_endmpi
      use m_inidef
      use m_inimpi
      use m_ininame
      use m_rdgrp
      use m_setdim
      use m_setnluni
      use m_unidrv

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

      call inimpi('unite   ',5)

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

      call setnluni

! -----

! Set the model dimension variables.

      call setdim(idcphopt,idhaiopt,idaslopt)

! -----

! Allocate the table of unit numbers.

      call allociot

! -----

! Allocate the array for unite.

      call allocuni(iduniopt_uni,ni,nj,nk,ni_uni,nj_uni,nio_uni)

! -----

! Allocate the group domain arrangement table.

      call allocgrp(idwbc,idebc,idsbc,idnbc)

! -----

! Read out and check the group domain arrangement.

      call rdgrp(idexprim,idcrsdir,idncexp,idnccrs)

! -----

! Generate the united file from the dumped files.

      call unidrv(idfltyp_uni,idiniopt,iddmpmon,                        &
     &            iduniopt_uni,idugroup_uni,idflitv_uni,ni,nj,nk,       &
     &            tmp1,tmp2,tmp3,tmp4,ni_uni,nj_uni,var,nio_uni,iodmp)

! -----

! Finalize the MPI processes.

      call endmpi(0)

! -----

!-----7--------------------------------------------------------------7--

! Internal procedure

!     none

!-----7--------------------------------------------------------------7--

      end program unite
