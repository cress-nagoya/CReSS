!***********************************************************************
      program solver
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/05/10
!     Modification: 1999/05/20, 1999/06/07, 1999/06/21, 1999/07/05,
!                   1999/07/23, 1999/07/28, 1999/08/03, 1999/08/09,
!                   1999/08/23, 1999/09/06, 1999/09/30, 1999/10/12,
!                   1999/10/22, 1999/11/01, 1999/11/19, 1999/11/24,
!                   1999/12/06, 2000/01/05, 2000/01/17, 2000/03/08,
!                   2000/03/17, 2000/04/18, 2000/06/01, 2000/07/05,
!                   2000/12/18, 2001/02/13, 2001/03/13, 2001/04/15,
!                   2001/06/29, 2001/07/13, 2001/08/07, 2001/09/13,
!                   2001/10/18, 2001/11/20, 2001/12/11, 2002/01/07,
!                   2002/01/15, 2002/04/02, 2002/04/09, 2002/06/18,
!                   2002/07/03, 2002/07/15, 2002/07/23, 2002/08/15,
!                   2002/09/09, 2002/10/15, 2002/10/31, 2002/12/02,
!                   2003/01/04, 2003/02/05, 2003/03/21, 2003/05/19,
!                   2003/07/15, 2003/11/05, 2003/11/28, 2003/12/12,
!                   2004/01/09, 2004/02/01, 2004/04/15, 2004/05/31,
!                   2004/06/10, 2004/07/01, 2004/08/01, 2004/08/20,
!                   2004/09/25, 2005/01/14, 2005/02/10, 2005/04/04,
!                   2006/01/10, 2006/02/13, 2006/04/03, 2006/06/21,
!                   2006/07/21, 2006/09/30, 2006/11/06, 2006/11/27,
!                   2006/12/04, 2007/01/05, 2007/01/20, 2007/01/31,
!                   2007/05/07, 2007/07/30, 2008/01/11, 2008/04/17,
!                   2008/05/02, 2008/07/01, 2008/08/25, 2008/10/10,
!                   2008/12/11, 2009/01/30, 2009/02/27, 2011/08/18,
!                   2011/09/22, 2011/11/10, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! This is the main program for the solver processer solver.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_allocbuf
      use m_allociot
      use m_allocmph
      use m_allocslv
      use m_chkbyte
      use m_chkkind
      use m_comindx
      use m_comslv
      use m_currpe
      use m_defdim
      use m_endmpi
      use m_inidef
      use m_inimod
      use m_inimpi
      use m_inivar
      use m_rdgrp
      use m_setdim
      use m_setgrid
      use m_sethalo
      use m_setnlslv
      use m_setref
      use m_slvdrv

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

      call inimpi('solver  ',6)

! -----

! Check the kind of the Fortran variables.

      call chkkind

! -----

! Initialize the model dimension and namelist variables.

      call inidef

! -----

! Initialize the common used module variables.

      call inimod

! -----

! Set the referenced table.

      call setref

! -----

! Set the namelist variables.

      call setnlslv

! -----

! Set the model dimension variables.

      call setdim(idcphopt,idhaiopt,idaslopt)

! -----

! Allocate the table of unit numbers.

      call allociot

! -----

! Allocate the array for solver.

      call allocslv(iddmpvar,idsavmem,idwbc,idebc,idsbc,idnbc,idnggopt, &
     &            idexbopt,idlspopt,idvspopt,idngropt,idiniopt,idsfcopt,&
     &            idadvopt,idcphopt,idqcgopt,idaslopt,idtrkopt,idtubopt,&
     &            ni,nj,nk,nqw,nnw,nqi,nni,km,nqa,nund,nlev,            &
     &            nid_rdr,njd_rdr,nkd_rdr,km_rdr)

! -----

! Allocate the array for cloud physics.

      call allocmph(idsavmem,idcphopt,idhaiopt,idaslopt,ni,nj,nk,nqw)

! -----

! Allocate the communication buffer.

      call allocbuf(idwbc,idebc,idsbc,idnbc,idgwmopt,idadvopt,idsmtopt, &
     &              idcphopt,idqcgopt,idaslopt,idtrkopt,idtubopt,       &
     &              ni,nj,nk,nqw,nnw,nqi,nni,nqa)

! -----

! Read out and check the group domain arrangement.

      call rdgrp(idexprim,idcrsdir,idncexp,idnccrs)

! -----

! Set the parameters of parallelizing.

      call currpe('all     ',3,'unset')

! -----

! Set the parameters of sendding and receiving.

      call sethalo(idwbc,idebc,idsbc,idnbc)

! -----

! Check the endian of running machine.

      call chkbyte

! -----

! Set the model grid.

      call s_setgrid(idwbc,idebc,idtrnopt,idexbopt,idmfcopt,idcoropt,   &
     &               area,ni,nj,nk,zph,zsth,lat,lon,j31,j32,jcb,        &
     &               jcb8u,jcb8v,jcb8w,mf,mf8u,mf8v,rmf,rmf8u,rmf8v,    &
     &               fc,tmp1,tmp2,tmp3)

! -----

! Set the initial conditions.

      call s_inivar(idiniopt,idmovopt,                                  &
     &              idlspopt,idvspopt,idsfcopt,idcphopt,idaslopt,       &
     &              fmois,ksp0,ni,nj,nk,nqw,nnw,nqi,nni,nqa,nund,       &
     &              zph,lat,lon,j31,j32,jcb,jcb8w,mf,fc,ubr,vbr,        &
     &              pbr,ptbr,qvbr,rbr,rst,rst8u,rst8v,rst8w,            &
     &              rcsq,u,up,v,vp,w,wp,pp,ppp,ptp,ptpp,qv,qvp,         &
     &              qwtr,qwtrp,nwtr,nwtrp,qice,qicep,nice,nicep,        &
     &              qcwtr,qcwtrp,qcice,qcicep,qasl,qaslp,               &
     &              qt,qtp,tke,tkep,rbcx,rbcy,rbcxy,rbct,               &
     &              ucpx,ucpy,vcpx,vcpy,wcpx,wcpy,pcpx,pcpy,            &
     &              ptcpx,ptcpy,qvcpx,qvcpy,qwcpx,qwcpy,                &
     &              nwcpx,nwcpy,qicpx,qicpy,nicpx,nicpy,                &
     &              qcwcpx,qcwcpy,qcicpx,qcicpy,qacpx,qacpy,            &
     &              qtcpx,qtcpy,tkecpx,tkecpy,maxvl,prwtr,price,        &
     &              pdia,land,albe,beta,z0m,z0h,cap,nuu,kai,            &
     &              tund,tundp,tmp1,tmp2,tmp3,nlev,z1d,u1d,v1d,         &
     &              p1d,pt1d,qv1d,ltmp1,ltmp2,ltmp3)

! -----

! Start the time integration.

      call slvdrv(ididate,idnggopt,idexbopt,idlspopt,idvspopt,          &
     &            idngropt,idsfcopt,idmasopt,idadvopt,idaslopt,         &
     &            idtubopt,idresopt,idmxnopt,fmois,ksp0,area,           &
     &            ni,nj,nk,nqw,nnw,nqi,nni,km,nqa,nund,zph,             &
     &            zsth,lat,lon,j31,j32,jcb,jcb8u,jcb8v,jcb8w,           &
     &            mf,mf8u,mf8v,rmf,rmf8u,rmf8v,fc,ubr,vbr,              &
     &            pbr,ptbr,qvbr,rbr,rst,rst8u,rst8v,rst8w,rcsq,         &
     &            rbcx,rbcy,rbcxy,rbct,land,albe,beta,cap,nuu,          &
     &            kai,u,up,v,vp,w,wp,pp,ppp,ptp,ptpp,qv,qvp,            &
     &            qwtr,qwtrp,nwtr,nwtrp,qice,qicep,nice,nicep,          &
     &            qcwtr,qcwtrp,qcice,qcicep,qasl,qaslp,qt,qtp,          &
     &            tke,tkep,ucpx,ucpy,vcpx,vcpy,wcpx,wcpy,               &
     &            pcpx,pcpy,ptcpx,ptcpy,qvcpx,qvcpy,qwcpx,qwcpy,        &
     &            nwcpx,nwcpy,qicpx,qicpy,nicpx,nicpy,                  &
     &            qcwcpx,qcwcpy,qcicpx,qcicpy,qacpx,qacpy,              &
     &            qtcpx,qtcpy,tkecpx,tkecpy,ugpv,utd,vgpv,vtd,          &
     &            wgpv,wtd,ppgpv,pptd,ptpgpv,ptptd,qvgpv,qvtd,          &
     &            qwgpv,qwtd,qigpv,qitd,qagpv,qatd,urdr,vrdr,           &
     &            wrdr,qwrdr,qwrtd,qirdr,qirtd,sst,sstd,                &
     &            maxvl,prwtr,price,qall,qallp,pdia,z0m,z0h,            &
     &            tund,tundp,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,             &
     &            tmp7,tmp8,tmp9,tmp10,tmp11,tmp12,tmp13,               &
     &            tmp14,tmp15,tmp16,tmp17,tmp18,tmp19,                  &
     &            qwtmp,nwtmp,qitmp,nitmp,qcwtmp,qcitmp,qatmp,          &
     &            qttmp,tketmp,tutmp,nid_rdr,njd_rdr,nkd_rdr,           &
     &            km_rdr,lon_rdr,z_rdr,u_rdr,v_rdr,w_rdr,qp_rdr,        &
     &            tmp1_rdr)

! -----

! Finalize the MPI processes.

      call endmpi(0)

! -----

!-----7--------------------------------------------------------------7--

! Internal procedure

!     none

!-----7--------------------------------------------------------------7--

      end program solver
