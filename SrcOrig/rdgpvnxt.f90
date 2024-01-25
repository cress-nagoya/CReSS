!***********************************************************************
      module m_rdgpvnxt
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/03/13
!     Modification: 2001/04/15, 2001/05/29, 2001/06/06, 2001/06/29,
!                   2002/04/02, 2002/04/09, 2002/06/18, 2002/07/15,
!                   2002/08/15, 2002/09/09, 2002/12/02, 2003/03/28,
!                   2003/04/30, 2003/05/19, 2003/06/27, 2003/11/05,
!                   2003/11/28, 2003/12/12, 2004/01/09, 2004/03/05,
!                   2004/04/15, 2004/05/31, 2004/06/10, 2004/07/01,
!                   2004/08/01, 2004/08/20, 2004/09/01, 2004/09/10,
!                   2004/10/12, 2005/01/14, 2005/02/10, 2006/01/10,
!                   2006/02/03, 2006/08/18, 2006/09/21, 2006/12/04,
!                   2007/01/05, 2007/01/20, 2007/05/07, 2007/05/14,
!                   2007/05/21, 2007/06/27, 2007/07/30, 2007/08/24,
!                   2007/09/25, 2008/05/02, 2008/07/01, 2008/07/25,
!                   2008/08/25, 2008/10/10, 2008/12/11, 2009/01/30,
!                   2009/02/27, 2009/03/23, 2011/08/18, 2011/09/22,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read out the data from the interpolated GPV data file.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_adjstq
      use m_bcgpv
      use m_chkerr
      use m_chkfile
      use m_chksat
      use m_chkstd
      use m_comindx
      use m_comkind
      use m_commpi
      use m_comname
      use m_comsave
      use m_cpondpe
      use m_destroy
      use m_getcname
      use m_getiname
      use m_getrname
      use m_getunit
      use m_inichar
      use m_outstd03
      use m_putunit
      use m_setgpv
      use m_smoogpv

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: rdgpvnxt, s_rdgpvnxt

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rdgpvnxt

        module procedure s_rdgpvnxt

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic int
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_rdgpvnxt(fpexprim,fpcrsdir,fpgpvvar,fpncexp,fpnccrs, &
     &                      fpwlngth,fpgsmopt,fpcphopt,fphaiopt,fpetime,&
     &                      fpgpvitv,fgpv,ctime,ftime,gtinc,ni,nj,nk,   &
     &                      nqw,nqi,pbr,ptbr,ugpv,utd,vgpv,vtd,wgpv,wtd,&
     &                      ppgpv,pptd,ptpgpv,ptptd,qvgpv,qvtd,         &
     &                      qwgpv,qwtd,qigpv,qitd,tmp1)
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

      integer, intent(in) :: fpetime
                       ! Formal parameter of unique index of etime

      integer, intent(in) :: fpgpvitv
                       ! Formal parameter of unique index of gpvitv

      integer, intent(in) :: fgpv
                       ! Descriptor to put into motion for GPV nudging

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

      integer(kind=i8), intent(in) :: ftime
                       ! Model forecast time at 1 step future

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

      real, intent(in) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

! Input and output variables

      real, intent(inout) :: ugpv(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity of GPV data
                       ! at marked time

      real, intent(inout) :: utd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! x components of velocity of GPV data

      real, intent(inout) :: vgpv(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity of GPV data
                       ! at marked time

      real, intent(inout) :: vtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! y components of velocity of GPV data

      real, intent(inout) :: wgpv(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity of GPV data
                       ! at marked time

      real, intent(inout) :: wtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! z components of velocity of GPV data

      real, intent(inout) :: ppgpv(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation of GPV data
                       ! at marked time

      real, intent(inout) :: pptd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! pressure perturbation of GPV data

      real, intent(inout) :: ptpgpv(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation of GPV data
                       ! at marked time

      real, intent(inout) :: ptptd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! potential temperature perturbation of GPV data

      real, intent(inout) :: qvgpv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio of GPV data
                       ! at marked time

      real, intent(inout) :: qvtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! water vapor mixing ratio of GPV data

      real, intent(inout) :: qwgpv(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor of GPV data at marked time

      real, intent(inout) :: qwtd(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Time tendency of water hydrometeor of GPV data

      real, intent(inout) :: qigpv(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor of GPV data at marked time

      real, intent(inout) :: qitd(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Time tendency of ice hydrometeor of GPV data

! Output variable

      real, intent(out) :: gtinc
                       ! Lapse of forecast time from GPV data reading

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

      integer ird      ! Index of count to read out

      integer ncfl     ! Number of character of GPV data file

      integer recgpv   ! Current record number of GPV data file

      integer iogpv    ! Unit number of GPV data file

      integer siz      ! Record length of GPV data file

      integer stat     ! Runtime status

      integer broot    ! Broadcasting root

      integer dategr   ! GPV data file date

      integer gpv01    ! int(gpvitv + 0.1)

      integer(kind=i8) gpv103
                       ! 1000 x int(gpvitv + 0.1)

      integer(kind=i8) crtime
                       ! Current forecast time to read out

      real etime       ! Forecast stop time

      real gpvitv      ! Time interval of GPV data

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getrname(fpetime,etime)
      call getrname(fpgpvitv,gpvitv)

! -----

! Set the common used variables.

      gpv01=int(gpvitv+.1e0)

      gpv103=1000_i8*int(gpvitv+.1e0,i8)

      crtime=gpv103*(ctime/gpv103)

      gtinc=.001e0*real(ctime-crtime)

      if(fgpv.eq.0) then

        nxtgpv=crtime/1000_i8

      end if

! -----

!!! Open and read out the data from the interpolated GPV data file when
!!! the current forecast time reaches marked time.

      if(ctime/1000_i8.ge.nxtgpv                                        &
     &  .and.ctime.lt.1000_i8*int(etime+.1e0,i8)) then

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

! -----

!! Open and read out the data from the interpolated GPV data checking
!! file.

        if(fgpv.eq.0) then

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

! Open the interpolated GPV data checking file.

          if(mype.eq.root) then

            gpvfl(1:ncexp)=exprim(1:ncexp)

            write(gpvfl(ncexp+1:ncexp+13),'(a13)') 'gpv.check.txt'

            open(iogpv,iostat=stat,err=100,                             &
     &           file=crsdir(1:nccrs)//gpvfl(1:ncexp+13),               &
     &           status='old',access='sequential',form='formatted',     &
     &           blank='null',position='rewind',action='read')

          else

            stat=0

          end if

  100     call chkerr(stat)

          if(stat.lt.0) then

            if(mype.eq.-stat-1) then

              call destroy('rdgpvnxt',8,'cont',1,'              ',14,   &
     &                     iogpv,stat)

            end if

            call cpondpe

            call destroy('rdgpvnxt',8,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          if(mype.eq.root) then

            call outstd03('rdgpvnxt',8,gpvfl,ncexp+13,iogpv,1,1,ftime)

          end if

          broot=stat-1

          call chkstd(broot)

! -----

! Read out the data from the interpolated GPV data checking file.

          if(mype.eq.root) then

            read(iogpv,'(a)',iostat=stat,end=110,err=110)               &
     &          (rcname(in),in=1,ncn)

            read(iogpv,*,iostat=stat,end=110,err=110)                   &
     &          (riname(in),in=1,nin)

            read(iogpv,*,iostat=stat,end=110,err=110)                   &
     &          (rrname(in),in=1,nrn)

          else

            stat=0

          end if

  110     call chkerr(stat)

          if(stat.lt.0) then

            if(mype.eq.-stat-1) then

              call destroy('rdgpvnxt',8,'cont',3,'              ',14,   &
     &                     iogpv,stat)

            end if

            call cpondpe

            call destroy('rdgpvnxt',8,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          if(mype.eq.root) then

            call outstd03('rdgpvnxt',8,gpvfl,108,iogpv,3,1,ftime)

          end if

          broot=stat-1

          call chkstd(broot)

! -----

! Check the interpolated GPV data checking file.

          if(mype.eq.root) then

            call chkfile('gpv',stat,ncn,nin,nrn,                        &
     &                   cname,iname,rname,rcname,riname,rrname)

          else

            stat=0

          end if

          call chkerr(stat)

          if(stat.lt.0) then

            call destroy('chkfile ',7,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

! -----

! Close the interpolated GPV data checking file.

          if(mype.eq.root) then

            close(iogpv,iostat=stat,err=120,status='keep')

          else

            stat=0

          end if

  120     call chkerr(stat)

          if(stat.lt.0) then

            if(mype.eq.-stat-1) then

              call destroy('rdgpvnxt',8,'cont',2,'              ',14,   &
     &                     iogpv,stat)

            end if

            call cpondpe

            call destroy('rdgpvnxt',8,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          if(mype.eq.root) then

            call outstd03('rdgpvnxt',8,gpvfl,108,iogpv,2,1,ftime)

          end if

          broot=stat-1

          call chkstd(broot)

! -----

! Return the unit number.

          if(mype.eq.root) then

            call putunit(iogpv)

          end if

! -----

        end if

!! -----

!! Open and read out the data from the interpolated GPV data file.

        siz=(ni+2)*(nj+2)*nk*wlngth

        do ird=2,1,-1

! Initialize the character variable.

          call inichar(gpvfl)

! -----

! Get the unit number.

          call getunit(iogpv)

! -----

! Open the interpolated GPV data file.

          dategr=crtime/1000_i8
          dategr=dategr+abs(ird-2)*gpv01

          gpvfl(1:ncexp)=exprim(1:ncexp)

          if(ngrp.eq.1) then

            ncfl=ncexp+22

            write(gpvfl(ncexp+1:ncexp+22),'(a3,i8.8,a3,i4.4,a4)')       &
     &                       'gpv',dategr,'.pe',mysub,'.bin'

          else

            ncfl=ncexp+31

            write(gpvfl(ncexp+1:ncexp+31),'(a3,i8.8,2(a4,i4.4),a4)')    &
     &                  'gpv',dategr,'.grp',mygrp,'-sub',mysub,'.bin'

          end if

          open(iogpv,iostat=stat,err=130,                               &
     &         file=crsdir(1:nccrs)//gpvfl(1:ncfl),                     &
     &         status='old',access='direct',form='unformatted',         &
     &         recl=siz,action='read')

  130     call chkerr(stat)

          if(stat.lt.0) then

            if(mype.eq.-stat-1) then

              call destroy('rdgpvnxt',8,'cont',1,'              ',14,   &
     &                     iogpv,stat)

            end if

            call cpondpe

            call destroy('rdgpvnxt',8,'stop',1001,'              ',14,  &
     &                   101,stat)

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

            call outstd03('rdgpvnxt',8,gpvfl,ncfl,iogpv,1,1,ftime)

          end if

          broot=stat-1

          call chkstd(broot)

! -----

! Read out the data from the interpolated GPV data file.

          if(dategr.eq.0) then

            if(gpvvar(2:2).eq.'o') then

              recgpv=9

              read(iogpv,rec=6,iostat=stat,err=140)                     &
     &            (((utd(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

              read(iogpv,rec=7,iostat=stat,err=140)                     &
     &            (((vtd(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

              read(iogpv,rec=8,iostat=stat,err=140)                     &
     &            (((pptd(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

              read(iogpv,rec=9,iostat=stat,err=140)                     &
     &            (((ptptd(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

            else

              recgpv=8

              read(iogpv,rec=5,iostat=stat,err=140)                     &
     &            (((utd(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

              read(iogpv,rec=6,iostat=stat,err=140)                     &
     &            (((vtd(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

              read(iogpv,rec=7,iostat=stat,err=140)                     &
     &            (((pptd(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

              read(iogpv,rec=8,iostat=stat,err=140)                     &
     &            (((ptptd(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

            end if

          else

            recgpv=4

            read(iogpv,rec=1,iostat=stat,err=140)                       &
     &          (((utd(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

            read(iogpv,rec=2,iostat=stat,err=140)                       &
     &          (((vtd(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

            read(iogpv,rec=3,iostat=stat,err=140)                       &
     &          (((pptd(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

            read(iogpv,rec=4,iostat=stat,err=140)                       &
     &          (((ptptd(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

          end if

          if(gpvvar(1:1).eq.'o') then
            recgpv=recgpv+1

            read(iogpv,rec=recgpv,iostat=stat,err=140)                  &
     &          (((wtd(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

          end if

          if(gpvvar(2:2).eq.'o') then
            recgpv=recgpv+1

            read(iogpv,rec=recgpv,iostat=stat,err=140)                  &
     &          (((qvtd(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

          end if

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              if(gpvvar(3:3).eq.'o') then
                recgpv=recgpv+1

                read(iogpv,rec=recgpv,iostat=stat,err=140)              &
     &              (((qwtd(i,j,k,1),i=0,ni+1),j=0,nj+1),k=1,nk)

              end if

              if(gpvvar(4:4).eq.'o') then
                recgpv=recgpv+1

                read(iogpv,rec=recgpv,iostat=stat,err=140)              &
     &              (((qwtd(i,j,k,2),i=0,ni+1),j=0,nj+1),k=1,nk)

              end if

            end if

            if(abs(cphopt).ge.2) then

              if(gpvvar(5:5).eq.'o') then
                recgpv=recgpv+1

                read(iogpv,rec=recgpv,iostat=stat,err=140)              &
     &              (((qitd(i,j,k,1),i=0,ni+1),j=0,nj+1),k=1,nk)

              end if

              if(gpvvar(6:6).eq.'o') then
                recgpv=recgpv+1

                read(iogpv,rec=recgpv,iostat=stat,err=140)              &
     &              (((qitd(i,j,k,2),i=0,ni+1),j=0,nj+1),k=1,nk)

              end if

              if(gpvvar(7:7).eq.'o') then
                recgpv=recgpv+1

                read(iogpv,rec=recgpv,iostat=stat,err=140)              &
     &              (((qitd(i,j,k,3),i=0,ni+1),j=0,nj+1),k=1,nk)

              end if

              if(gpvvar(8:8).eq.'o') then
                recgpv=recgpv+1

                if(haiopt.eq.0) then

                  read(iogpv,rec=recgpv,iostat=stat,err=140)            &
     &                (((tmp1(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

                  do k=1,nk
                  do j=0,nj+1
                  do i=0,ni+1
                    qitd(i,j,k,3)=qitd(i,j,k,3)+tmp1(i,j,k)
                  end do
                  end do
                  end do

                else

                  read(iogpv,rec=recgpv,iostat=stat,err=140)            &
     &                (((qitd(i,j,k,4),i=0,ni+1),j=0,nj+1),k=1,nk)

                end if

              end if

            end if

          end if

  140     call chkerr(stat)

          if(stat.lt.0) then

            if(mype.eq.-stat-1) then

              call destroy('rdgpvnxt',8,'cont',3,'              ',14,   &
     &                     iogpv,stat)

            end if

            call cpondpe

            call destroy('rdgpvnxt',8,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          if(mype.eq.stat-1) then

            call outstd03('rdgpvnxt',8,gpvfl,108,iogpv,3,1,ftime)

          end if

          broot=stat-1

          call chkstd(broot)

! -----

! Close the interpolated GPV data file.

          close(iogpv,iostat=stat,err=150,status='keep')

  150     call chkerr(stat)

          if(stat.lt.0) then

            if(mype.eq.-stat-1) then

              call destroy('rdgpvnxt',8,'cont',2,'              ',14,   &
     &                     iogpv,stat)

            end if

            call cpondpe

            call destroy('rdgpvnxt',8,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          if(mype.eq.stat-1) then

            call outstd03('rdgpvnxt',8,gpvfl,108,iogpv,2,1,ftime)

          end if

          broot=stat-1

          call chkstd(broot)

! -----

! Return the unit number.

          call putunit(iogpv)

! -----

! Check and avoid the super saturation mixing ratio.

          if(gpvvar(2:2).eq.'o') then

            call chksat('total','xxx',0,ni+1,0,nj+1,1,nk,pbr,ptbr,      &
     &                  pptd,ptptd,qvtd)

          end if

! -----

! Set the boundary conditions.

          call bcgpv(idgpvvar,idwbc,idebc,idexbopt,idcphopt,idhaiopt,   &
     &               ni,nj,nk,nqw,nqi,utd,vtd,wtd,pptd,ptptd,qvtd,      &
     &               qwtd,qitd)

! -----

! Smooth the interpolated GPV data.

          if(gsmopt.eq.1) then

            call smoogpv(idgpvvar,idcphopt,idhaiopt,idgsmcnt,           &
     &                   ni,nj,nk,nqw,nqi,utd,vtd,wtd,pptd,ptptd,       &
     &                   qvtd,qwtd,qitd,tmp1)

          end if

! -----

! Force mixing ratio more than user specified value.

          call adjstq(idcphopt,idhaiopt,ni,nj,nk,nqw,nqi,qvtd,qwtd,qitd)

! -----

! Set the time tendency of variables at current marked time and the
! variables at initial or restart time.

          call setgpv(idgpvvar,idcphopt,idhaiopt,idgpvitv,ird,ni,nj,nk, &
     &                nqw,nqi,ugpv,utd,vgpv,vtd,wgpv,wtd,ppgpv,pptd,    &
     &                ptpgpv,ptptd,qvgpv,qvtd,qwgpv,qwtd,qigpv,qitd)

! -----

        end do

!! -----

! Calculate the next time to read out.

        nxtgpv=nxtgpv+gpv01

! -----

      end if

!!! -----

      end subroutine s_rdgpvnxt

!-----7--------------------------------------------------------------7--

      end module m_rdgpvnxt
