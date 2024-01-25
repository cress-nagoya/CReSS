!***********************************************************************
      program asldata
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/09/22
!     Modification: 2011/11/10, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! This is the main program for the pre processor asldata.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_allocasl
      use m_allocgrp
      use m_allociot
      use m_asldrv
      use m_chkkind
      use m_comasl
      use m_comindx
      use m_defdim
      use m_endmpi
      use m_inidef
      use m_inimpi
      use m_ininame
      use m_rdgrp
      use m_setdays
      use m_setdim
      use m_setnlasl

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

      call inimpi('asldata ',7)

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

      call setnlasl

! -----

! Set the model dimension variables.

      call setdim(idcphopt,idhaiopt,idaslopt)

! -----

! Allocate the table of unit numbers.

      call allociot

! -----

! Allocate the array for asldata.

      call allocasl(ni,nj,nk,nqa,nid_asl,njd_asl,nkd_asl,km_asl)

! -----

! Allocate the group domain arrangement table.

      call allocgrp(idwbc,idebc,idsbc,idnbc)

! -----

! Read out and check the group domain arrangement.

      call rdgrp(idexprim,idcrsdir,idncexp,idnccrs)

! -----

! Interpolate the aerosol data.

      call asldrv(ididate,idaslitv,ni,nj,nk,nqa,                        &
     &            z,zph,qasl,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,             &
     &            nid_asl,njd_asl,nkd_asl,km_asl,zdat,qadat,dtmp1)

! -----

! Finalize the MPI processes.

      call endmpi(0)

! -----

!-----7--------------------------------------------------------------7--

! Internal procedure

!     none

!-----7--------------------------------------------------------------7--

      end program asldata
