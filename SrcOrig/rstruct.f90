!***********************************************************************
      program rstruct
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2005/07/01
!     Modification: 2005/08/05, 2006/01/10, 2006/02/13, 2006/04/03,
!                   2006/07/21, 2006/09/30, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/05/07, 2007/07/30, 2008/04/17,
!                   2008/05/02, 2008/08/25, 2008/10/10, 2009/01/30,
!                   2009/02/27, 2011/08/18, 2011/09/22, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! This is the main program for the post processor rstruct.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_allocgrp
      use m_allociot
      use m_allocrst
      use m_chkkind
      use m_comindx
      use m_comrst
      use m_defdim
      use m_endmpi
      use m_inidef
      use m_inimpi
      use m_ininame
      use m_rdgrp
      use m_rstdrv
      use m_setdim
      use m_setnlrst

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

      call inimpi('rstruct ',7)

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

      call setnlrst

! -----

! Set the model dimension variables.

      call setdim(idcphopt,idhaiopt,idaslopt)

! -----

! Allocate the table of unit numbers.

      call allociot

! -----

! Allocate the array for rstruct.

      call allocrst(iddmpvar,idsavmem,idwbc,idebc,idsbc,idnbc,          &
     &              idsfcopt,idadvopt,idcphopt,idqcgopt,idaslopt,       &
     &              idtrkopt,idtubopt,iddiaopt,ni,nj,nk,nqw,nnw,nqi,nni,&
     &              nqa,nund,ni_rst,nj_rst)

! -----

! Allocate the group domain arrangement table.

      call allocgrp(idwbc,idebc,idsbc,idnbc)

! -----

! Read out and check the group domain arrangement.

      call rdgrp(idexprim,idcrsdir,idncexp,idnccrs)

! -----

! Restruct the restart files.

      call rstdrv(idrmopt_rst,idflitv_rst,                              &
     &            ni,nj,nk,nqw,nnw,nqi,nni,nqa,nund,ni_rst,nj_rst,      &
     &            ubr,vbr,pbr,ptbr,qvbr,u,up,v,vp,w,wp,pp,ppp,ptp,ptpp, &
     &            qv,qvp,qwtr,qwtrp,nwtr,nwtrp,qice,qicep,nice,nicep,   &
     &            qcwtr,qcwtrp,qcice,qcicep,qasl,qaslp,qt,qtp,tke,tkep, &
     &            ucpx,ucpy,vcpx,vcpy,wcpx,wcpy,pcpx,pcpy,ptcpx,ptcpy,  &
     &            qvcpx,qvcpy,qwcpx,qwcpy,nwcpx,nwcpy,                  &
     &            qicpx,qicpy,nicpx,nicpy,qcwcpx,qcwcpy,qcicpx,qcicpy,  &
     &            qacpx,qacpy,qtcpx,qtcpy,tkecpx,tkecpy,maxvl,          &
     &            prwtr,price,pdia,z0m,z0h,tund,tundp,ubr_rst,vbr_rst,  &
     &            pbr_rst,ptbr_rst,qvbr_rst,u_rst,up_rst,v_rst,vp_rst,  &
     &            w_rst,wp_rst,pp_rst,ppp_rst,ptp_rst,ptpp_rst,         &
     &            qv_rst,qvp_rst,qwtr_rst,qwtrp_rst,nwtr_rst,nwtrp_rst, &
     &            qice_rst,qicep_rst,nice_rst,nicep_rst,                &
     &            qcwtr_rst,qcwtrp_rst,qcice_rst,qcicep_rst,            &
     &            qasl_rst,qaslp_rst,qt_rst,qtp_rst,tke_rst,tkep_rst,   &
     &            ucpx_rst,ucpy_rst,vcpx_rst,vcpy_rst,                  &
     &            wcpx_rst,wcpy_rst,pcpx_rst,pcpy_rst,                  &
     &            ptcpx_rst,ptcpy_rst,qvcpx_rst,qvcpy_rst,              &
     &            qwcpx_rst,qwcpy_rst,nwcpx_rst,nwcpy_rst,              &
     &            qicpx_rst,qicpy_rst,nicpx_rst,nicpy_rst,              &
     &            qcwcpx_rst,qcwcpy_rst,qcicpx_rst,qcicpy_rst,          &
     &            qacpx_rst,qacpy_rst,qtcpx_rst,qtcpy_rst,              &
     &            tkecpx_rst,tkecpy_rst,maxvl_rst,prwtr_rst,price_rst,  &
     &            pdia_rst,z0m_rst,z0h_rst,tund_rst,tundp_rst)

! -----

! Finalize the MPI processes.

      call endmpi(0)

! -----

!-----7--------------------------------------------------------------7--

! Internal procedure

!     none

!-----7--------------------------------------------------------------7--

      end program rstruct
