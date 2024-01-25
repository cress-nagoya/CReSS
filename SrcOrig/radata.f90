!***********************************************************************
      program radata
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/09/09
!     Modification: 2002/12/02, 2003/02/05, 2003/05/19, 2003/07/15,
!                   2003/11/05, 2003/12/12, 2004/01/09, 2004/08/01,
!                   2005/02/10, 2005/08/05, 2006/01/10, 2006/09/30,
!                   2006/12/04, 2007/01/05, 2007/01/20, 2007/07/30,
!                   2008/04/17, 2008/05/02, 2008/08/25, 2008/10/10,
!                   2009/01/30, 2009/02/27, 2011/08/18, 2011/09/22,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! This is the main program for the pre processor radata.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_allocgrp
      use m_allociot
      use m_allocrdr
      use m_chkkind
      use m_comindx
      use m_comrdr
      use m_defdim
      use m_endmpi
      use m_inidef
      use m_inimpi
      use m_ininame
      use m_rdgrp
      use m_rdrdrv
      use m_setdays
      use m_setdim
      use m_setnlrdr

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

      call inimpi('radata  ',6)

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

! Set the referenced number of days.

      call setdays

! -----

! Set the namelist variables.

      call setnlrdr

! -----

! Set the model dimension variables.

      call setdim(idcphopt,idhaiopt,idaslopt)

! -----

! Allocate the table of unit numbers.

      call allociot

! -----

! Allocate the array for radata.

      call allocrdr(ni,nj,nk,nid_rdr,njd_rdr,nkd_rdr,km_rdr)

! -----

! Allocate the group domain arrangement table.

      call allocgrp(idwbc,idebc,idsbc,idnbc)

! -----

! Read out and check the group domain arrangement.

      call rdgrp(idexprim,idcrsdir,idncexp,idnccrs)

! -----

! Interpolate the radar data.

      call rdrdrv(idrdrvar,iddatype_rdr,ididate,idrotopt_rdr,idrdritv,  &
     &            idngrstr,ni,nj,nk,z,zph,u,v,w,qp,tmp1,tmp2,tmp3,      &
     &            tmp4,tmp5,tmp6,nid_rdr,njd_rdr,nkd_rdr,km_rdr,        &
     &            londat,zdat,udat,vdat,wdat,qpdat,dtmp1)

! -----

! Finalize the MPI processes.

      call endmpi(0)

! -----

!-----7--------------------------------------------------------------7--

! Internal procedure

!     none

!-----7--------------------------------------------------------------7--

      end program radata
