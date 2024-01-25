!***********************************************************************
      program surface
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/07/03
!     Modification: 2002/07/15, 2002/09/09, 2002/12/02, 2003/02/05,
!                   2003/05/19, 2003/07/15, 2003/11/05, 2003/12/12,
!                   2004/01/09, 2004/05/07, 2004/08/01, 2004/09/01,
!                   2005/01/14, 2005/02/10, 2006/01/10, 2006/09/30,
!                   2006/12/04, 2007/01/05, 2007/01/20, 2007/07/30,
!                   2008/04/17, 2008/05/02, 2008/08/25, 2008/10/10,
!                   2009/01/30, 2009/02/27, 2011/08/18, 2011/09/22,
!                   2011/11/10, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! This is the main program for the pre processor surface.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_allocgrp
      use m_allociot
      use m_allocsfc
      use m_chkkind
      use m_comindx
      use m_comsfc
      use m_defdim
      use m_endmpi
      use m_inidef
      use m_inimpi
      use m_ininame
      use m_rdgrp
      use m_setdays
      use m_setdim
      use m_setnlsfc
      use m_sfcdrv

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

      call inimpi('surface ',7)

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

      call setnlsfc

! -----

! Set the model dimension variables.

      call setdim(idcphopt,idhaiopt,idaslopt)

! -----

! Allocate the table of unit numbers.

      call allociot

! -----

! Allocate the array for surface.

      call allocsfc(ni,nj,nid_lnd,njd_lnd,nid_sst,njd_sst,              &
     &              nid_ice,njd_ice)

! -----

! Allocate the group domain arrangement table.

      call allocgrp(idwbc,idebc,idsbc,idnbc)

! -----

! Read out and check the group domain arrangement.

      call rdgrp(idexprim,idcrsdir,idncexp,idnccrs)

! -----

! Interpolate the surface data.

      call sfcdrv(idsfcdat,ididate,idintopt_lnd,idlnduse_lnd,           &
     &           idsstitv,idalbe_lnd,idbeta_lnd,idz0m_lnd,idz0h_lnd,    &
     &           idcap_lnd,idnuu_lnd,ni,nj,ri,rj,land,albe,beta,z0m,z0h,&
     &           cap,nuu,sst,kai,tmp1,tmp2,tmp3,nid_lnd,njd_lnd,landat, &
     &           ltmp1,nid_sst,njd_sst,sstdat,stmp1,nid_ice,njd_ice,    &
     &           icedat,itmp1)

! -----

! Finalize the MPI processes.

      call endmpi(0)

! -----

!-----7--------------------------------------------------------------7--

! Internal procedure

!     none

!-----7--------------------------------------------------------------7--

      end program surface
