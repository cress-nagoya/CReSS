!***********************************************************************
      module m_rdrdrnxt
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/07/30
!     Modification: 2007/08/24, 2007/09/25, 2008/05/02, 2008/07/25,
!                   2008/08/25, 2008/10/10, 2009/01/05, 2009/01/30,
!                   2009/02/27, 2009/03/23, 2009/11/13, 2011/08/18,
!                   2011/09/22, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read out the data from the interpolated radar data file.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bcrdrqp
      use m_chkerr
      use m_chkfile
      use m_chkstd
      use m_comindx
      use m_comkind
      use m_commpi
      use m_comname
      use m_comsave
      use m_cpondpe
      use m_destroy
      use m_distrqp
      use m_getcname
      use m_getiname
      use m_getrname
      use m_getunit
      use m_inichar
      use m_outstd03
      use m_putunit
      use m_setrdrqp

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: rdrdrnxt, s_rdrdrnxt

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rdrdrnxt

        module procedure s_rdrdrnxt

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic int
      intrinsic max
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_rdrdrnxt(fpexprim,fpcrsdir,fprdrvar,                 &
     &                      fpncexp,fpnccrs,fpwlngth,fpcphopt,          &
     &                      fprdritv,fpngrstr,fpngrend,frdr,ctime,ftime,&
     &                      rtinc,ni,nj,nk,nqw,nqi,rbr,qwtr,qice,       &
     &                      qwrdr,qwrtd,qirdr,qirtd,qprdr)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpexprim
                       ! Formal parameter of unique index of exprim

      integer, intent(in) :: fpcrsdir
                       ! Formal parameter of unique index of crsdir

      integer, intent(in) :: fprdrvar
                       ! Formal parameter of unique index of rdrvar

      integer, intent(in) :: fpncexp
                       ! Formal parameter of unique index of ncexp

      integer, intent(in) :: fpnccrs
                       ! Formal parameter of unique index of nccrs

      integer, intent(in) :: fpwlngth
                       ! Formal parameter of unique index of wlngth

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fprdritv
                       ! Formal parameter of unique index of rdritv

      integer, intent(in) :: fpngrstr
                       ! Formal parameter of unique index of ngrstr

      integer, intent(in) :: fpngrend
                       ! Formal parameter of unique index of ngrend

      integer, intent(in) :: frdr(1:2)
                       ! Descriptor to put into motion for radar nudging

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

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: qwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor

      real, intent(in) :: qice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor

! Input and output variables

      real, intent(inout) :: qwrdr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor of radar data
                       ! at marked time

      real, intent(inout) :: qwrtd(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Time tendency of
                       ! water hydrometeor of radar data

      real, intent(inout) :: qirdr(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor of radar data
                       ! at marked time

      real, intent(inout) :: qirtd(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Time tendency of
                       ! ice hydrometeor of radar data

! Output variable

      real, intent(out) :: rtinc(1:2)
                       ! Lapse of forecast time from radar data reading

! Internal shared variables

      character(len=108) exprim
                       ! Optional run name

      character(len=108) crsdir
                       ! User specified directory for CReSS files

      character(len=108) rdrvar
                       ! Control flag of input radar data variables

      character(len=108) rdrfl
                       ! Opened file name

      integer ncexp    ! Number of character of exprim
      integer nccrs    ! Number of character of crsdir

      integer wlngth   ! Word length of direct access file

      integer cphopt   ! Option for cloud micro physics

      integer in       ! Namelist table index

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer ird      ! Index of count to read out

      integer ncfl     ! Number of character of radar data file

      integer recrdr   ! Current record number of radar data file

      integer iordr    ! Unit number of radar data file

      integer siz      ! Record length of radar data file

      integer stat     ! Runtime status

      integer broot    ! Broadcasting root

      integer daterd   ! Radar data file date

      integer rdr01    ! int(rdritv + 0.1)

      integer(kind=i8) str01
                       ! int(ngrstr + 0.1)

      integer(kind=i8) rdr103
                       ! 1000 x int(rdritv + 0.1)

      integer(kind=i8) str103
                       ! 1000 x int(ngrstr + 0.1)

      integer(kind=i8) end103
                       ! 1000 x int(ngrend + 0.1)

      integer(kind=i8) crtime
                       ! Current forecast time to read out

      real rdritv      ! Time interval of radar data

      real ngrstr      ! Analysis nudging start time to radar data
      real ngrend      ! Analysis nudging end time to radar data

      real, intent(inout) :: qprdr(0:ni+1,0:nj+1,1:nk)
                       ! Precipitation mixing ratio of radar data
                       ! at marked time

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(rdrvar)

! -----

! Get the required namelist variables.

      call getcname(fprdrvar,rdrvar)
      call getrname(fprdritv,rdritv)
      call getrname(fpngrstr,ngrstr)
      call getrname(fpngrend,ngrend)

! -----

! Set the common used variables.

      rdr01=int(rdritv+.1e0)
      str01=int(ngrstr+.1e0)

      rdr103=1000_i8*int(rdritv+.1e0,i8)
      str103=1000_i8*int(ngrstr+.1e0,i8)
      end103=1000_i8*int(ngrend+.1e0,i8)

      crtime=rdr103*((ctime-str103)/rdr103)+str103

      rtinc(1)=.001e0*real(ctime-crtime)

      if(frdr(1).eq.0) then

        nxtrdr(1)=max(crtime/1000_i8,str01)

      end if

! -----

!!! Open and read out the data from the interpolated radar data file
!!! when the current forecast time reaches marked time.

      if(rdrvar(4:4).eq.'o') then

        if(ctime/1000_i8.ge.nxtrdr(1).and.ctime.lt.end103) then

! Initialize the character variables.

          call inichar(exprim)
          call inichar(crsdir)

! -----

! Get the required namelist variables.

          call getcname(fpexprim,exprim)
          call getcname(fpcrsdir,crsdir)
          call getiname(fpncexp,ncexp)
          call getiname(fpnccrs,nccrs)
          call getiname(fpwlngth,wlngth)
          call getiname(fpcphopt,cphopt)

! -----

!! Open and read out the data from the interpolated radar data checking
!! file.

          if(frdr(1).eq.0.or.nxtrdr(1).eq.str01) then

! Initialize the character variable.

            if(mype.eq.root) then

              call inichar(rdrfl)

            end if

! -----

! Get the unit number.

            if(mype.eq.root) then

              call getunit(iordr)

            end if

! -----

! Open the interpolated radar data checking file.

            if(mype.eq.root) then

              rdrfl(1:ncexp)=exprim(1:ncexp)

              write(rdrfl(ncexp+1:ncexp+13),'(a13)') 'rdr.check.txt'

              open(iordr,iostat=stat,err=100,                           &
     &             file=crsdir(1:nccrs)//rdrfl(1:ncexp+13),             &
     &             status='old',access='sequential',form='formatted',   &
     &             blank='null',position='rewind',action='read')

            else

              stat=0

            end if

  100       call chkerr(stat)

            if(stat.lt.0) then

              if(mype.eq.-stat-1) then

                call destroy('rdrdrnxt',8,'cont',1,'              ',14, &
     &                       iordr,stat)

              end if

              call cpondpe

              call destroy('rdrdrnxt',8,'stop',1001,'              ',14,&
     &                     101,stat)

            end if

            if(mype.eq.root) then

              call outstd03('rdrdrnxt',8,rdrfl,ncexp+13,iordr,1,1,ftime)

            end if

            broot=stat-1

            call chkstd(broot)

! -----

! Read out the data from the interpolated radar data checking file.

            if(mype.eq.root) then

              read(iordr,'(a)',iostat=stat,end=110,err=110)             &
     &            (rcname(in),in=1,ncn)

              read(iordr,*,iostat=stat,end=110,err=110)                 &
     &            (riname(in),in=1,nin)

              read(iordr,*,iostat=stat,end=110,err=110)                 &
     &            (rrname(in),in=1,nrn)

            else

              stat=0

            end if

  110       call chkerr(stat)

            if(stat.lt.0) then

              if(mype.eq.-stat-1) then

                call destroy('rdrdrnxt',8,'cont',3,'              ',14, &
     &                       iordr,stat)

              end if

              call cpondpe

              call destroy('rdrdrnxt',8,'stop',1001,'              ',14,&
     &                     101,stat)

            end if

            if(mype.eq.root) then

              call outstd03('rdrdrnxt',8,rdrfl,108,iordr,3,1,ftime)

            end if

            broot=stat-1

            call chkstd(broot)

! -----

! Check the interpolated radar data checking file.

            if(mype.eq.root) then

              call chkfile('rdr',stat,ncn,nin,nrn,                      &
     &                     cname,iname,rname,rcname,riname,rrname)

            else

              stat=0

            end if

            call chkerr(stat)

            if(stat.lt.0) then

              call destroy('chkfile ',7,'stop',1001,'              ',14,&
     &                     101,stat)

            end if

! -----

! Close the interpolated radar data checking file.

            if(mype.eq.root) then

              close(iordr,iostat=stat,err=120,status='keep')

            else

              stat=0

            end if

  120       call chkerr(stat)

            if(stat.lt.0) then

              if(mype.eq.-stat-1) then

                call destroy('rdrdrnxt',8,'cont',2,'              ',14, &
     &                       iordr,stat)

              end if

              call cpondpe

              call destroy('rdrdrnxt',8,'stop',1001,'              ',14,&
     &                     101,stat)

            end if

            if(mype.eq.root) then

              call outstd03('rdrdrnxt',8,rdrfl,108,iordr,2,1,ftime)

            end if

            broot=stat-1

            call chkstd(broot)

! -----

! Return the unit number.

            if(mype.eq.root) then

              call putunit(iordr)

            end if

! -----

          end if

!! -----

!! Open and read out the data from the interpolated radar data file.

          siz=(ni+2)*(nj+2)*nk*wlngth

          do ird=2,1,-1

! Initialize the character variable.

            call inichar(rdrfl)

! -----

! Get the unit number.

            call getunit(iordr)

! -----

! Open the interpolated radar data file.

            daterd=crtime/1000_i8
            daterd=daterd+abs(ird-2)*rdr01

            rdrfl(1:ncexp)=exprim(1:ncexp)

            if(ngrp.eq.1) then

              ncfl=ncexp+22

              write(rdrfl(ncexp+1:ncexp+22),'(a3,i8.8,a3,i4.4,a4)')     &
     &                         'rdr',daterd,'.pe',mysub,'.bin'

            else

              ncfl=ncexp+31

              write(rdrfl(ncexp+1:ncexp+31),'(a3,i8.8,2(a4,i4.4),a4)')  &
     &            'rdr',daterd,'.grp',mygrp,'-sub',mysub,'.bin'

            end if

            open(iordr,iostat=stat,err=130,                             &
     &           file=crsdir(1:nccrs)//rdrfl(1:ncfl),                   &
     &           status='old',access='direct',form='unformatted',       &
     &           recl=siz,action='read')

  130       call chkerr(stat)

            if(stat.lt.0) then

              if(mype.eq.-stat-1) then

                call destroy('rdrdrnxt',8,'cont',1,'              ',14, &
     &                       iordr,stat)

              end if

              call cpondpe

              call destroy('rdrdrnxt',8,'stop',1001,'              ',14,&
     &                     101,stat)

            end if

            if(mype.eq.stat-1) then

              if(fpara(1:5).eq.'multi') then

                if(ngrp.eq.1) then

                  write(rdrfl(ncexp+15:ncexp+18),'(a4)') 'XXXX'

                else

                  write(rdrfl(ncexp+16:ncexp+19),'(a4)') 'XXXX'
                  write(rdrfl(ncexp+24:ncexp+27),'(a4)') 'YYYY'

                end if

              end if

              call outstd03('rdrdrnxt',8,rdrfl,ncfl,iordr,1,1,ftime)

            end if

            broot=stat-1

            call chkstd(broot)

! -----

! Read out the data from the interpolated radar data file.

            recrdr=0

            if(rdrvar(1:1).eq.'o') then
              recrdr=recrdr+1
            end if

            if(rdrvar(2:2).eq.'o') then
              recrdr=recrdr+1
            end if

            if(rdrvar(3:3).eq.'o') then
              recrdr=recrdr+1
            end if

            recrdr=recrdr+1

            read(iordr,rec=recrdr,iostat=stat,err=140)                  &
     &          (((qprdr(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

  140       call chkerr(stat)

            if(stat.lt.0) then

              if(mype.eq.-stat-1) then

                call destroy('rdrdrnxt',8,'cont',3,'              ',14, &
     &                       iordr,stat)

              end if

              call cpondpe

              call destroy('rdrdrnxt',8,'stop',1001,'              ',14,&
     &                     101,stat)

            end if

            if(mype.eq.stat-1) then

              call outstd03('rdrdrnxt',8,rdrfl,108,iordr,3,1,ftime)

            end if

            broot=stat-1

            call chkstd(broot)

! -----

! Close the interpolated radar data file.

            close(iordr,iostat=stat,err=150,status='keep')

  150       call chkerr(stat)

            if(stat.lt.0) then

              if(mype.eq.-stat-1) then

                call destroy('rdrdrnxt',8,'cont',2,'              ',14, &
     &                       iordr,stat)

              end if

              call cpondpe

              call destroy('rdrdrnxt',8,'stop',1001,'              ',14,&
     &                     101,stat)

            end if

            if(mype.eq.stat-1) then

              call outstd03('rdrdrnxt',8,rdrfl,108,iordr,2,1,ftime)

            end if

            broot=stat-1

            call chkstd(broot)

! -----

! Return the unit number.

            call putunit(iordr)

! -----

! Set the bondary conditions for interpolated radar data.

            call bcrdrqp(idwbc,idebc,idexbopt,ni,nj,nk,qprdr)

! -----

! Distribute the observed precipitation to the rain, snow and graupel
! mixing ratio.

            if(abs(cphopt).lt.10) then

              call distrqp(iddatype_rdr,idcphopt,idhaiopt,idthresq,     &
     &                 ni,nj,nk,nqw,nqi,rbr,qwtr,qice,qprdr,qwrtd,qirtd)

            end if

! -----

! Set the time tendency of variables at current marked time and the
! variables at initial or restart time.

            if(abs(cphopt).lt.10) then

              call setrdrqp(idcphopt,idhaiopt,idrdritv,ird,             &
     &                      ni,nj,nk,nqw,nqi,qwrdr,qwrtd,qirdr,qirtd)

            end if

! -----

          end do

!! -----

! Calculate the next time to read out.

          nxtrdr(1)=nxtrdr(1)+rdr01

! -----

        end if

      end if

!!! -----

      end subroutine s_rdrdrnxt

!-----7--------------------------------------------------------------7--

      end module m_rdrdrnxt
