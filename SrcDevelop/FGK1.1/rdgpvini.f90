!***********************************************************************
      module m_rdgpvini
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/05/10
!     Modification: 1999/05/20, 1999/06/07, 1999/06/14, 1999/06/21,
!                   1999/07/05, 1999/08/03, 1999/08/23, 1999/09/30,
!                   1999/11/01, 1999/11/19, 2000/01/17, 2000/04/18,
!                   2001/01/15, 2001/02/13, 2001/03/13, 2001/04/15,
!                   2001/05/29, 2001/06/29, 2002/04/02, 2002/04/09,
!                   2002/06/18, 2002/07/18, 2002/08/15, 2002/12/02,
!                   2003/03/28, 2003/04/30, 2003/05/19, 2003/06/27,
!                   2003/11/28, 2003/12/12, 2004/01/09, 2004/03/05,
!                   2004/04/15, 2004/05/31, 2004/08/20, 2004/09/01,
!                   2004/09/10, 2004/10/12, 2005/01/14, 2005/02/10,
!                   2006/01/10, 2006/02/03, 2006/09/21, 2006/12/04,
!                   2007/01/05, 2007/01/20, 2007/05/07, 2007/05/14,
!                   2007/05/21, 2007/06/27, 2007/08/24, 2008/05/02,
!                   2008/07/01, 2008/08/25, 2008/10/10, 2008/12/11,
!                   2009/01/30, 2009/02/27, 2011/08/18, 2011/09/22,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     open and read out the data from the interpolated GPV data file at
!     forecast start time.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_adjstq
      use m_bcini
      use m_chkerr
      use m_chkfile
      use m_chkmoist
      use m_chksat
      use m_chkstd
      use m_comindx
      use m_comkind
      use m_commpi
      use m_comname
      use m_cpondpe
      use m_destroy
      use m_getcname
      use m_getiname
      use m_getrname
      use m_getunit
      use m_inichar
      use m_outstd03
      use m_putunit
      use m_setcst3d
      use m_smoogpv

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: rdgpvini, s_rdgpvini

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rdgpvini

        module procedure s_rdgpvini

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic int

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_rdgpvini(fpexprim,fpcrsdir,fpgpvvar,fpncexp,fpnccrs, &
     &                      fpwlngth,fpgsmopt,fpcphopt,fphaiopt,fpstime,&
     &                      fmois,ni,nj,nk,nqw,nqi,ubr,vbr,pbr,ptbr,    &
     &                      qvbr,u,v,w,pp,ptp,qv,qwtr,qice,tmp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpexprim
                       ! Formal parameter of unique index of exprim

      integer, intent(in) :: fpcrsdir
                       ! Formal parameter of unique index of crsdir

      integer, intent(in) :: fpgpvvar
                       ! Formal parameter of unique index of gpvvar

      integer, intent(in) :: fpncexp
                       ! Formal parameter of unique index of ncexp

      integer, intent(in) :: fpnccrs
                       ! Formal parameter of unique index of nccrs

      integer, intent(in) :: fpwlngth
                       ! Formal parameter of unique index of wlngth

      integer, intent(in) :: fpgsmopt
                       ! Formal parameter of unique index of gsmopt

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fphaiopt
                       ! Formal parameter of unique index of haiopt

      integer, intent(in) :: fpstime
                       ! Formal parameter of unique index of stime

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

      integer, intent(in) :: nqi
                       ! Number of categories of ice hydrometeor

! Output variables

      character(len=5), intent(out) :: fmois
                       ! Control flag of air moisture

      real, intent(out) :: ubr(0:ni+1,0:nj+1,1:nk)
                       ! Base state x components of velocity

      real, intent(out) :: vbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state y components of velocity

      real, intent(out) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(out) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(out) :: qvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state water vapor mixing ratio

      real, intent(out) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(out) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(out) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

      real, intent(out) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation

      real, intent(out) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

      real, intent(out) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

      real, intent(out) :: qwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor

      real, intent(out) :: qice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor

! Internal shared variables

      character(len=108) exprim
                       ! Optional run name

      character(len=108) crsdir
                       ! User specified directory for CReSS files

      character(len=108) gpvvar
                       ! Control flag of input GPV data variables

      character(len=108) gpvfl
                       ! Opened file name

      integer ncexp    ! Number of character of exprim
      integer nccrs    ! Number of character of crsdir

      integer wlngth   ! Word length of direct access file

      integer gsmopt   ! Option for GPV data smoothing
      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes

      integer in       ! Namelist table index

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer ncfl     ! Number of character of GPV data file

      integer recgpv   ! Current record number of GPV data file

      integer iogpv    ! Unit number of GPV data file

      integer siz      ! Record length of GPV data file

      integer stat     ! Runtime status

      integer broot    ! Broadcasting root

      real stime       ! Forecast start time

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

!-----7--------------------------------------------------------------7--

! Initialize the character variables.

      call inichar(exprim)
      call inichar(crsdir)
      call inichar(gpvvar)

! -----

! Get the required namelist variables.

      call getcname(fpexprim,exprim)
      call getcname(fpcrsdir,crsdir)
      call getcname(fpgpvvar,gpvvar)
      call getiname(fpncexp,ncexp)
      call getiname(fpnccrs,nccrs)
      call getiname(fpwlngth,wlngth)
      call getiname(fpgsmopt,gsmopt)
      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)
      call getrname(fpstime,stime)

! -----

!! Open and read out the data from the interpolated GPV data checking
!! file at forecast start time.

! Initialize the character variable.

      if(mype.eq.root) then

        call inichar(gpvfl)

      end if

! -----

! Get the unit number.

      if(mype.eq.root) then

        call getunit(iogpv)

      end if

! -----

! Open the interpolated GPV data checking file at forecast start time.

      if(mype.eq.root) then

        gpvfl(1:ncexp)=exprim(1:ncexp)

        write(gpvfl(ncexp+1:ncexp+13),'(a13)') 'gpv.check.txt'

        open(iogpv,iostat=stat,err=100,                                 &
     &       file=crsdir(1:nccrs)//gpvfl(1:ncexp+13),                   &
     &       status='old',access='sequential',form='formatted',         &
     &       blank='null',position='rewind',action='read')

      else

        stat=0

      end if

  100 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdgpvini',8,'cont',1,'              ',14,iogpv, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdgpvini',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('rdgpvini',8,gpvfl,ncexp+13,iogpv,1,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read out the data from the interpolated GPV data checking file at
! forecast start time.

      if(mype.eq.root) then

        read(iogpv,'(a)',iostat=stat,end=110,err=110)                   &
     &      (rcname(in),in=1,ncn)

        read(iogpv,*,iostat=stat,end=110,err=110)                       &
     &      (riname(in),in=1,nin)

        read(iogpv,*,iostat=stat,end=110,err=110)                       &
     &      (rrname(in),in=1,nrn)

      else

        stat=0

      end if

  110 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdgpvini',8,'cont',3,'              ',14,iogpv, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdgpvini',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('rdgpvini',8,gpvfl,108,iogpv,3,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Check the interpolated GPV data checking file at forecast start time.

      if(mype.eq.root) then

        call chkfile('gpv',stat,ncn,nin,nrn,                            &
     &               cname,iname,rname,rcname,riname,rrname)

      else

        stat=0

      end if

      call chkerr(stat)

      if(stat.lt.0) then

        call destroy('chkfile ',7,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

! Close the interpolated GPV data checking file at forecast start time.

      if(mype.eq.root) then

        close(iogpv,iostat=stat,err=120,status='keep')

      else

        stat=0

      end if

  120 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdgpvini',8,'cont',2,'              ',14,iogpv, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdgpvini',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('rdgpvini',8,gpvfl,108,iogpv,2,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit number.

      if(mype.eq.root) then

        call putunit(iogpv)

      end if

! -----

!! -----

!! Open and read out the data from the interpolated GPV data file at
!! forecast start time.

! Initialize the character variable.

      call inichar(gpvfl)

! -----

! Get the unit number.

      call getunit(iogpv)

! -----

! Open the interpolated GPV data file at forecast start time.

      siz=(ni+2)*(nj+2)*nk*wlngth

      gpvfl(1:ncexp)=exprim(1:ncexp)

      if(ngrp.eq.1) then

        ncfl=ncexp+22

        write(gpvfl(ncexp+1:ncexp+22),'(a3,i8.8,a3,i4.4,a4)')           &
     &          'gpv',int(stime+.1e0),'.pe',mysub,'.bin'

      else

        ncfl=ncexp+31

        write(gpvfl(ncexp+1:ncexp+31),'(a3,i8.8,2(a4,i4.4),a4)')        &
     &     'gpv',int(stime+.1e0),'.grp',mygrp,'-sub',mysub,'.bin'

      end if

      open(iogpv,iostat=stat,err=130,                                   &
     &     file=crsdir(1:nccrs)//gpvfl(1:ncfl),                         &
     &     status='old',access='direct',form='unformatted',             &
     &     recl=siz,action='read')

  130 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdgpvini',8,'cont',1,'              ',14,iogpv, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdgpvini',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        if(fpara(1:5).eq.'multi') then

          if(ngrp.eq.1) then

            write(gpvfl(ncexp+15:ncexp+18),'(a4)') 'XXXX'

          else

            write(gpvfl(ncexp+16:ncexp+19),'(a4)') 'XXXX'
            write(gpvfl(ncexp+24:ncexp+27),'(a4)') 'YYYY'

          end if

        end if

        call outstd03('rdgpvini',8,gpvfl,ncfl,iogpv,1,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read out the data from the interpolated GPV data file at forecast
! start time.

      if(gpvvar(2:2).eq.'o') then

        recgpv=9

        read(iogpv,rec=1,iostat=stat,err=140)                           &
     &      (((ubr(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

        read(iogpv,rec=2,iostat=stat,err=140)                           &
     &      (((vbr(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

        read(iogpv,rec=3,iostat=stat,err=140)                           &
     &      (((pbr(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

        read(iogpv,rec=4,iostat=stat,err=140)                           &
     &      (((ptbr(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

        read(iogpv,rec=5,iostat=stat,err=140)                           &
     &      (((qvbr(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

        read(iogpv,rec=6,iostat=stat,err=140)                           &
     &      (((u(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

        read(iogpv,rec=7,iostat=stat,err=140)                           &
     &      (((v(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

        read(iogpv,rec=8,iostat=stat,err=140)                           &
     &      (((pp(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

        read(iogpv,rec=9,iostat=stat,err=140)                           &
     &      (((ptp(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      else

        recgpv=8

        read(iogpv,rec=1,iostat=stat,err=140)                           &
     &      (((ubr(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

        read(iogpv,rec=2,iostat=stat,err=140)                           &
     &      (((vbr(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

        read(iogpv,rec=3,iostat=stat,err=140)                           &
     &      (((pbr(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

        read(iogpv,rec=4,iostat=stat,err=140)                           &
     &      (((ptbr(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

        read(iogpv,rec=5,iostat=stat,err=140)                           &
     &      (((u(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

        read(iogpv,rec=6,iostat=stat,err=140)                           &
     &      (((v(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

        read(iogpv,rec=7,iostat=stat,err=140)                           &
     &      (((pp(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

        read(iogpv,rec=8,iostat=stat,err=140)                           &
     &      (((ptp(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      end if

      if(gpvvar(1:1).eq.'o') then
        recgpv=recgpv+1

        read(iogpv,rec=recgpv,iostat=stat,err=140)                      &
     &      (((w(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      end if

      if(gpvvar(2:2).eq.'o') then
        recgpv=recgpv+1

        read(iogpv,rec=recgpv,iostat=stat,err=140)                      &
     &      (((qv(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      end if

      if(abs(cphopt).lt.10) then

        if(abs(cphopt).ge.1) then

          if(gpvvar(3:3).eq.'o') then
            recgpv=recgpv+1

            read(iogpv,rec=recgpv,iostat=stat,err=140)                  &
     &          (((qwtr(i,j,k,1),i=0,ni+1),j=0,nj+1),k=1,nk)

          end if

          if(gpvvar(4:4).eq.'o') then
            recgpv=recgpv+1

            read(iogpv,rec=recgpv,iostat=stat,err=140)                  &
     &          (((qwtr(i,j,k,2),i=0,ni+1),j=0,nj+1),k=1,nk)

          end if

        end if

        if(abs(cphopt).ge.2) then

          if(gpvvar(5:5).eq.'o') then
            recgpv=recgpv+1

            read(iogpv,rec=recgpv,iostat=stat,err=140)                  &
     &          (((qice(i,j,k,1),i=0,ni+1),j=0,nj+1),k=1,nk)

          end if

          if(gpvvar(6:6).eq.'o') then
            recgpv=recgpv+1

            read(iogpv,rec=recgpv,iostat=stat,err=140)                  &
     &          (((qice(i,j,k,2),i=0,ni+1),j=0,nj+1),k=1,nk)

          end if

          if(gpvvar(7:7).eq.'o') then
            recgpv=recgpv+1

            read(iogpv,rec=recgpv,iostat=stat,err=140)                  &
     &          (((qice(i,j,k,3),i=0,ni+1),j=0,nj+1),k=1,nk)

          end if

          if(gpvvar(8:8).eq.'o') then
            recgpv=recgpv+1

            if(haiopt.eq.0) then

              read(iogpv,rec=recgpv,iostat=stat,err=140)                &
     &            (((tmp1(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

              do k=1,nk
              do j=0,nj+1
              do i=0,ni+1
                qice(i,j,k,3)=qice(i,j,k,3)+tmp1(i,j,k)
              end do
              end do
              end do

            else

              read(iogpv,rec=recgpv,iostat=stat,err=140)                &
     &            (((qice(i,j,k,4),i=0,ni+1),j=0,nj+1),k=1,nk)

            end if

          end if

        end if

      end if

  140 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdgpvini',8,'cont',3,'              ',14,iogpv, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdgpvini',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('rdgpvini',8,gpvfl,108,iogpv,3,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Close the interpolated GPV data file at forecast start time.

      close(iogpv,iostat=stat,err=150,status='keep')

  150 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdgpvini',8,'cont',2,'              ',14,iogpv, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdgpvini',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('rdgpvini',8,gpvfl,108,iogpv,2,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit number.

      call putunit(iogpv)

! -----

!! -----

! In the case of no data, force the base state water vapor mixing ratio
! to 0.

      if(gpvvar(2:2).eq.'x') then

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qvbr)

      end if

! -----

! Check and avoid the super saturation mixing ratio.

      if(gpvvar(2:2).eq.'o') then

        call chksat('total','xxx',0,ni+1,0,nj+1,1,nk,pbr,ptbr,pp,ptp,qv)

      end if

! -----

! Set the boundary conditions.

      call bcini(idgpvvar,idwbc,idebc,idexbopt,idcphopt,idhaiopt,       &
     &           ni,nj,nk,nqw,nqi,ubr,vbr,pbr,ptbr,qvbr,                &
     &           u,v,w,pp,ptp,qv,qwtr,qice)

! -----

! Smooth the interpolated GPV data.

      if(gsmopt.eq.1) then

        call smoogpv(idgpvvar,idcphopt,idhaiopt,idgsmcnt,               &
     &               ni,nj,nk,nqw,nqi,u,v,w,pp,ptp,qv,qwtr,qice,tmp1)

      end if

! -----

! Force mixing ratio more than user specified value.

      call adjstq(idcphopt,idhaiopt,ni,nj,nk,nqw,nqi,qv,qwtr,qice)

! -----

! Check the air moisture.

      call chkmoist(fmois,ni,nj,nk,qvbr)

! -----

      end subroutine s_rdgpvini

!-----7--------------------------------------------------------------7--

      end module m_rdgpvini
