!***********************************************************************
      module m_sethalo
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/12/04
!     Modification: 2007/01/05, 2007/01/20, 2008/05/02, 2008/07/25,
!                   2008/08/25, 2009/01/05, 2009/01/30, 2009/02/27,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the parameters of sendding and receiving.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_comgrp
      use m_comindx
      use m_commpi
      use m_comname
      use m_cpondpe
      use m_defmpi
      use m_destroy
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: sethalo, s_sethalo

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface sethalo

        module procedure s_sethalo

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_sethalo(fpwbc,fpebc,fpsbc,fpnbc)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpsbc
                       ! Formal parameter of unique index of sbc

      integer, intent(in) :: fpnbc
                       ! Formal parameter of unique index of nbc

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions
      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

      integer holes    ! Number of holes
                       ! in group domain arrangement to southward
      integer holen    ! Number of holes
                       ! in group domain arrangement to northward

      integer mytmp    ! Temporary used my sub domain number
                       ! in group domain

      integer ntmp     ! Temporary used number of sub domain
                       ! in group domain

      integer stat     ! Runtime status

      integer ierr     ! Error descriptor

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpsbc,sbc)
      call getiname(fpnbc,nbc)

! -----

!! Set the sendding destination and receiving source.

! For the sub domain.

      dstw_sub=mysub-1
      dste_sub=mysub+1
      srcw_sub=mysub+1
      srce_sub=mysub-1

      if(dstw_sub.lt.jsub*nisub) then
        dstw_sub=dstw_sub+nisub
      end if

      if(dste_sub.ge.(jsub+1)*nisub) then
        dste_sub=dste_sub-nisub
      end if

      if(srcw_sub.ge.(jsub+1)*nisub) then
        srcw_sub=srcw_sub-nisub
      end if

      if(srce_sub.lt.jsub*nisub) then
        srce_sub=srce_sub+nisub
      end if

      dsts_sub=mysub-nisub
      dstn_sub=mysub+nisub
      srcs_sub=mysub+nisub
      srcn_sub=mysub-nisub

      if(dsts_sub.lt.0) then
        dsts_sub=dsts_sub+nsub
      end if

      if(dstn_sub.ge.nsub) then
        dstn_sub=dstn_sub-nsub
      end if

      if(srcs_sub.ge.nsub) then
        srcs_sub=srcs_sub-nsub
      end if

      if(srcn_sub.lt.0) then
        srcn_sub=srcn_sub+nsub
      end if

      dstw_sub_bnd=dstw_sub
      dste_sub_bnd=dste_sub
      dsts_sub_bnd=dsts_sub
      dstn_sub_bnd=dstn_sub

      srcw_sub_bnd=srcw_sub
      srce_sub_bnd=srce_sub
      srcs_sub_bnd=srcs_sub
      srcn_sub_bnd=srcn_sub

! -----

! For the group domain.

      holes=(mygrp-mysrl)-((nigrp*(jgrp-1)+igrp)-grpxy(igrp+1,jgrp))
      holen=((nigrp*(jgrp+1)+igrp)-grpxy(igrp+1,jgrp+2))-(mygrp-mysrl)

      if(sbc.eq.1.and.nbc.eq.1) then
        holes=0
        holen=0
      end if

      dstw_grp=(mygrp-1)*nsub+(jsub+1)*nisub-1
      dste_grp=(mygrp+1)*nsub+jsub*nisub
      srcw_grp=(mygrp+1)*nsub+jsub*nisub
      srce_grp=(mygrp-1)*nsub+(jsub+1)*nisub-1

      if(dstw_grp.lt.jgrp*nigrp*nsub) then
        dstw_grp=dstw_grp+nigrp*nsub
      end if

      if(dste_grp.ge.(jgrp+1)*nigrp*nsub) then
        dste_grp=dste_grp-nigrp*nsub
      end if

      if(srcw_grp.ge.(jgrp+1)*nigrp*nsub) then
        srcw_grp=srcw_grp-nigrp*nsub
      end if

      if(srce_grp.lt.jgrp*nigrp*nsub) then
        srce_grp=srce_grp+nigrp*nsub
      end if

      dstw_grp=dstw_grp-(mygrp-mysrl)*nsub
      dste_grp=dste_grp-(mygrp-mysrl)*nsub
      srcw_grp=srcw_grp-(mygrp-mysrl)*nsub
      srce_grp=srce_grp-(mygrp-mysrl)*nsub

      dsts_grp=(mygrp-nigrp)*nsub+nisub*(njsub-1)+isub
      dstn_grp=(mygrp+nigrp)*nsub+isub
      srcs_grp=(mygrp+nigrp)*nsub+isub
      srcn_grp=(mygrp-nigrp)*nsub+nisub*(njsub-1)+isub

      if(dsts_grp.lt.0) then
        dsts_grp=dsts_grp+nsub*ngrp
      end if

      if(dstn_grp.ge.nsub*ngrp) then
        dstn_grp=dstn_grp-nsub*ngrp
      end if

      if(srcs_grp.ge.nsub*ngrp) then
        srcs_grp=srcs_grp-nsub*ngrp
      end if

      if(srcn_grp.lt.0) then
        srcn_grp=srcn_grp+nsub*ngrp
      end if

      dsts_grp=dsts_grp-(mygrp-mysrl-holes)*nsub
      dstn_grp=dstn_grp-(mygrp-mysrl+holen)*nsub
      srcs_grp=srcs_grp-(mygrp-mysrl+holen)*nsub
      srcn_grp=srcn_grp-(mygrp-mysrl-holes)*nsub

      dstw_grp_bnd=dstw_grp
      dste_grp_bnd=dste_grp
      dsts_grp_bnd=dsts_grp
      dstn_grp_bnd=dstn_grp

      srcw_grp_bnd=srcw_grp
      srce_grp_bnd=srce_grp
      srcs_grp_bnd=srcs_grp
      srcn_grp_bnd=srcn_grp

! -----

! Set the additional parameters.

      if(isub.ne.0) then
        dstw_grp=mpi_proc_null
        srce_grp=mpi_proc_null

        dstw_grp_bnd=mpi_proc_null
        srce_grp_bnd=mpi_proc_null

      end if

      if(isub.ne.nisub-1) then
        dste_grp=mpi_proc_null
        srcw_grp=mpi_proc_null

        dste_grp_bnd=mpi_proc_null
        srcw_grp_bnd=mpi_proc_null

      end if

      if(jsub.ne.0) then
        dsts_grp=mpi_proc_null
        srcn_grp=mpi_proc_null

        dsts_grp_bnd=mpi_proc_null
        srcn_grp_bnd=mpi_proc_null

      end if

      if(jsub.ne.njsub-1) then
        dstn_grp=mpi_proc_null
        srcs_grp=mpi_proc_null

        dstn_grp_bnd=mpi_proc_null
        srcs_grp_bnd=mpi_proc_null

      end if

      if(wbc.eq.-1) then

        if(isub.ne.0) then
          dstw_sub_bnd=mpi_proc_null
        end if

        if(isub.ne.nisub-1) then
          srcw_sub_bnd=mpi_proc_null
        end if

      else if(wbc.ne.-1) then

        if(isub.eq.0) then
          dstw_sub=mpi_proc_null
        end if

        if(isub.eq.nisub-1) then
          srcw_sub=mpi_proc_null
        end if

      end if

      if(ebc.eq.-1) then

        if(isub.ne.nisub-1) then
          dste_sub_bnd=mpi_proc_null
        end if

        if(isub.ne.0) then
          srce_sub_bnd=mpi_proc_null
        end if

      else if(ebc.ne.-1) then

        if(isub.eq.nisub-1) then
          dste_sub=mpi_proc_null
        end if

        if(isub.eq.0) then
          srce_sub=mpi_proc_null
        end if

      end if

      if(sbc.eq.-1) then

        if(jsub.ne.0) then
          dsts_sub_bnd=mpi_proc_null
        end if

        if(jsub.ne.njsub-1) then
          srcs_sub_bnd=mpi_proc_null
        end if

      else if(sbc.ne.-1) then

        if(jsub.eq.0) then
          dsts_sub=mpi_proc_null
        end if

        if(jsub.eq.njsub-1) then
          srcs_sub=mpi_proc_null
        end if

      end if

      if(nbc.eq.-1) then

        if(jsub.ne.njsub-1) then
          dstn_sub_bnd=mpi_proc_null
        end if

        if(jsub.ne.0) then
          srcn_sub_bnd=mpi_proc_null
        end if

      else if(nbc.ne.-1) then

        if(jsub.eq.njsub-1) then
          dstn_sub=mpi_proc_null
        end if

        if(jsub.eq.0) then
          srcn_sub=mpi_proc_null
        end if

      end if

      if(wbc.eq.1) then

        if(igrp.ne.0) then
          dstw_grp_bnd=mpi_proc_null
        end if

        if(igrp.ne.nigrp-1) then
          srcw_grp_bnd=mpi_proc_null
        end if

      else if(wbc.ne.1) then

        if(ebw.eq.1) then
          dstw_grp=mpi_proc_null
        end if

        if(ebe.eq.1) then
          srcw_grp=mpi_proc_null
        end if

      end if

      if(ebc.eq.1) then

        if(igrp.ne.nigrp-1) then
          dste_grp_bnd=mpi_proc_null
        end if

        if(igrp.ne.0) then
          srce_grp_bnd=mpi_proc_null
        end if

      else if(ebc.ne.1) then

        if(ebe.eq.1) then
          dste_grp=mpi_proc_null
        end if

        if(ebw.eq.1) then
          srce_grp=mpi_proc_null
        end if

      end if

      if(sbc.eq.1) then

        if(jgrp.ne.0) then
          dsts_grp_bnd=mpi_proc_null
        end if

        if(jgrp.ne.njgrp-1) then
          srcs_grp_bnd=mpi_proc_null
        end if

      else if(sbc.ne.1) then

        if(ebs.eq.1) then
          dsts_grp=mpi_proc_null
        end if

        if(ebn.eq.1) then
          srcs_grp=mpi_proc_null
        end if

      end if

      if(nbc.eq.1) then

        if(jgrp.ne.njgrp-1) then
          dstn_grp_bnd=mpi_proc_null
        end if

        if(jgrp.ne.0) then
          srcn_grp_bnd=mpi_proc_null
        end if

      else if(nbc.ne.1) then

        if(ebn.eq.1) then
          dstn_grp=mpi_proc_null
        end if

        if(ebs.eq.1) then
          srcn_grp=mpi_proc_null
        end if

      end if

! -----

!! -----

!! Get the splited new communicator.

! Initialize the runtime status, stat.

      stat=0

! -----

! Get the new communicator.

      call mpi_comm_split(mpi_comm_cress,mysrl,mype,mpi_comm_cress_sub, &
     &                    ierr)

! -----

! Get the number of entried processor elements.

      call mpi_comm_size(mpi_comm_cress_sub,ntmp,ierr)

      if(nsub.ne.ntmp) then
        stat=stat+1
      end if

! -----

! Get my processor element number.

      call mpi_comm_rank(mpi_comm_cress_sub,mytmp,ierr)

      if(mysub.ne.mytmp) then
        stat=stat+1
      end if

! -----

! If error occured, call the procedure destroy.

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('sethalo ',7,'cont',7,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('sethalo ',7,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

!! -----

! Convert the added or subtracted indices on boundary.

      if(ebw.eq.0.or.isub.ne.0) then
        iname(idiwest)=0
      else

        if(abs(wbc).eq.1) then
          iname(idiwest)=0
        else
          iname(idiwest)=1
        end if

      end if

      if(ebe.eq.0.or.isub.ne.nisub-1) then
        iname(idieast)=0
      else

        if(abs(ebc).eq.1) then
          iname(idieast)=0
        else
          iname(idieast)=1
        end if

      end if

      if(ebs.eq.0.or.jsub.ne.0) then
        iname(idjsouth)=0
      else

        if(abs(sbc).eq.1) then
          iname(idjsouth)=0
        else
          iname(idjsouth)=1
        end if

      end if

      if(ebn.eq.0.or.jsub.ne.njsub-1) then
        iname(idjnorth)=0
      else

        if(abs(nbc).eq.1) then
          iname(idjnorth)=0
        else
          iname(idjnorth)=1
        end if

      end if

! -----

      end subroutine s_sethalo

!-----7--------------------------------------------------------------7--

      end module m_sethalo
