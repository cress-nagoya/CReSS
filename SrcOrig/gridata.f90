!***********************************************************************
      program gridata
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/05/20
!     Modification: 1999/06/07, 1999/06/21, 1999/06/28, 1999/07/05,
!                   1999/08/03, 1999/09/01, 1999/09/30, 1999/10/22,
!                   1999/11/01, 1999/11/19, 1999/11/24, 2000/01/05,
!                   2000/01/17, 2000/04/18, 2000/12/18, 2001/01/15,
!                   2001/02/13, 2001/03/13, 2001/04/15, 2001/06/29,
!                   2001/10/18, 2002/01/07, 2002/04/09, 2002/06/18,
!                   2002/07/15, 2002/08/15, 2002/09/02, 2002/09/09,
!                   2002/12/02, 2003/02/05, 2003/05/19, 2003/07/15,
!                   2003/11/05, 2003/12/12, 2004/01/09, 2004/07/01,
!                   2004/08/01, 2005/02/10, 2005/08/05, 2006/01/10,
!                   2006/02/03, 2006/09/30, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/07/30, 2008/04/17, 2008/05/02,
!                   2008/08/19, 2008/08/25, 2008/10/10, 2008/12/11,
!                   2009/01/30, 2009/02/27, 2011/08/18, 2011/09/22,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! This is the main program for the pre processor gridata.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_allocgpv
      use m_allocgrp
      use m_allociot
      use m_chkkind
      use m_comgpv
      use m_comindx
      use m_defdim
      use m_endmpi
      use m_gpvdrv
      use m_inidef
      use m_inimpi
      use m_ininame
      use m_rdgrp
      use m_setdays
      use m_setdim
      use m_setnlgpv

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

      call inimpi('gridata ',7)

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

      call setnlgpv

! -----

! Set the model dimension variables.

      call setdim(idcphopt,idhaiopt,idaslopt)

! -----

! Allocate the table of unit numbers.

      call allociot

! -----

! Allocate the array for gridata.

      call allocgpv(ni,nj,nk,nid_gpv,njd_gpv,nkd_gpv,km_gpv)

! -----

! Allocate the group domain arrangement table.

      call allocgrp(idwbc,idebc,idsbc,idnbc)

! -----

! Read out and check the group domain arrangement.

      call rdgrp(idexprim,idcrsdir,idncexp,idnccrs)

! -----

! Interpolate the grid data.

      call gpvdrv(idgpvvar,ididate,iddatype_gpv,idtrnopt,idexbopt,      &
     &            idrotopt_gpv,idgpvitv,ni,nj,nk,land,z,zph,ubr,vbr,    &
     &            pbr,ptbr,qvbr,u,v,w,pp,ptp,qv,qc,qr,qi,qs,qg,qh,      &
     &            z1d,u1d,v1d,p1d,pt1d,qv1d,tmp1,tmp2,tmp3,tmp4,tmp5,   &
     &            tmp6,tmp7,nid_gpv,njd_gpv,nkd_gpv,km_gpv,londat,htdat,&
     &            zdat,udat,vdat,wdat,pdat,ptdat,qvdat,qcdat,qrdat,     &
     &            qidat,qsdat,qgdat,qhdat,pbdat,ptbdat,                 &
     &            dtmp1,dtmp2)

! -----

! Finalize the MPI processes.

      call endmpi(0)

! -----

!-----7--------------------------------------------------------------7--

! Internal procedure

!     none

!-----7--------------------------------------------------------------7--

      end program gridata
