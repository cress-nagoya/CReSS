!***********************************************************************
      module m_outctl
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/04/11
!     Modification: 2007/05/07, 2007/05/14, 2007/06/21, 2007/07/30,
!                   2007/08/24, 2007/10/19, 2007/11/13, 2007/11/26,
!                   2008/01/11, 2008/04/17, 2008/05/02, 2008/07/25,
!                   2008/08/25, 2008/10/10, 2009/01/30, 2009/02/27,
!                   2011/09/22, 2013/01/28, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     create the GrADS control file.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_chkfile
      use m_chkstd
      use m_comindx
      use m_comkind
      use m_commpi
      use m_comname
      use m_copy1d
      use m_cpondpe
      use m_destroy
      use m_getcname
      use m_getiname
      use m_getrname
      use m_getunit
      use m_getz
      use m_inichar
      use m_outstd03
      use m_paractl
      use m_putunit
      use m_stretch

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: outctl, s_outctl

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outctl

        module procedure s_outctl

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic int
      intrinsic mod
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_outctl(fpexprim,fpcrsdir,fpncexp,fpnccrs,            &
     &                    fpmpopt,fpnspol,fpsthopt,fpdmplev,            &
     &                    fpuniopt_uni,fpdx,fpdy,fpdz,fptlat1,fptlat2,  &
     &                    fptlon,fpflitv_uni,fproc,nstp0,nstp1,ni,nj,nk,&
     &                    z1d,zsth,dzsth,mlat)
!***********************************************************************

! Input variables

      character(len=3), intent(in) :: fproc
                       ! Control flag of processing type

      integer, intent(in) :: fpexprim
                       ! Formal parameter of unique index of exprim

      integer, intent(in) :: fpcrsdir
                       ! Formal parameter of unique index of crsdir

      integer, intent(in) :: fpncexp
                       ! Formal parameter of unique index of ncexp

      integer, intent(in) :: fpnccrs
                       ! Formal parameter of unique index of nccrs

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpnspol
                       ! Formal parameter of unique index of nspol

      integer, intent(in) :: fpsthopt
                       ! Formal parameter of unique index of sthopt

      integer, intent(in) :: fpdmplev
                       ! Formal parameter of unique index of dmplev

      integer, intent(in) :: fpuniopt_uni
                       ! Formal parameter of unique index of uniopt_uni

      integer, intent(in) :: fpdx
                       ! Formal parameter of unique index of dx

      integer, intent(in) :: fpdy
                       ! Formal parameter of unique index of dy

      integer, intent(in) :: fpdz
                       ! Formal parameter of unique index of dz

      integer, intent(in) :: fptlat1
                       ! Formal parameter of unique index of tlat1

      integer, intent(in) :: fptlat2
                       ! Formal parameter of unique index of tlat2

      integer, intent(in) :: fptlon
                       ! Formal parameter of unique index of tlon

      integer, intent(in) :: fpflitv_uni
                       ! Formal parameter of unique index of flitv_uni

      integer(kind=i8), intent(in) :: nstp0
                       ! Start index of main do loop

      integer(kind=i8), intent(in) :: nstp1
                       ! End index of main do loop

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

! Internal shared variables

      character(len=108) exprim
                       ! Optional run name

      character(len=108) crsdir
                       ! User specified directory for CReSS files

      character(len=75) vcap
                       ! Caption for dumped variable

      character(len=108) chkfl
                       ! Opened checking file name

      character(len=108) ctlfl
                       ! Opened control file name

      integer ncexp    ! Number of character of exprim
      integer nccrs    ! Number of character of crsdir

      integer mpopt    ! Option for map projection
      integer nspol    ! Option for projected region
      integer sthopt   ! Option for vertical grid stretching

      integer dmplev   ! Option for z coordinates of dumped variables

      integer uniopt_uni
                       ! Option for uniting process

      integer nym3     ! (nj - 3) x njgrp x njsub
      integer nyrm3    ! (nj - 3) x njsub x (jsred + njred)

      integer nkm2     ! nk - 2

      integer j        ! Array index in y drection
      integer k        ! Array index in z drection

      integer icnt     ! Optional used counter

      integer in       ! Namelist table index

      integer iflitv   ! int(flitv_uni + 0.1)

      integer varcnt   ! United variables count

      integer ncvc     ! Number of character of vcap

      integer ncchk    ! Number of character of checking file name
      integer ncctl    ! Number of character of control file name

      integer iochk    ! Unit number of checking file
      integer ioctl    ! Unit number of control file

      integer stat     ! Runtime status

      integer broot    ! Broadcasting root

      real dx          ! Grid distance in x direction
      real dy          ! Grid distance in y direction
      real dz          ! Grid distance in z direction

      real tlat1       ! True latitude 1
      real tlat2       ! True latitude 2
      real tlon        ! True longitude

      real flitv_uni   ! Time interval of processed file

      real latsw       ! Latitude at south-west corner
      real lonsw       ! Longitude at south-west corner

      real latmin      ! Minimum latitude in model domain
      real lonmin      ! Minimum longitude in model domain

      real dlat        ! Grid distance in latitude direction
      real dlon        ! Grid distance in longitude direction

      real ipole       ! Reference real index at pole in x direction
      real jpole       ! Reference real index at pole in y direction

      real, intent(inout) :: z1d(1:nk)
                       ! Constant height at scalar points

      real, intent(inout) :: zsth(1:nk)
                       ! 1 dimensional stretched z coordinates

      real, intent(inout) :: dzsth(1:nk)
                       ! Distance of
                       ! 1 dimensional stretched z coordinates

      real, intent(inout) :: mlat(1:(nj-3)*njgrp*njsub)
                       ! Latitude with Mercator projection

! Internal private variable

      integer k_sub    ! Substitute for k

!-----7--------------------------------------------------------------7--

! Initialize the character variables.

      call inichar(exprim)
      call inichar(crsdir)

! -----

! Get the required namelist variables.

      call getcname(fpexprim,exprim)
      call getcname(fpcrsdir,crsdir)
      call getiname(fpncexp,ncexp)
      call getiname(fpnccrs,nccrs)
      call getiname(fpmpopt,mpopt)
      call getiname(fpnspol,nspol)
      call getiname(fpsthopt,sthopt)
      call getiname(fpdmplev,dmplev)
      call getiname(fpuniopt_uni,uniopt_uni)
      call getrname(fpdx,dx)
      call getrname(fpdy,dy)
      call getrname(fpdz,dz)
      call getrname(fptlat1,tlat1)
      call getrname(fptlat2,tlat2)
      call getrname(fptlon,tlon)
      call getrname(fpflitv_uni,flitv_uni)

! -----

! Set the common used variables.

      nym3=(nj-3)*njgrp*njsub
      nyrm3=(nj-3)*njsub*(jsred+njred)

      nkm2=nk-2

! -----

! Set the parameters of GrADS control file.

      if(uniopt_uni.ge.-10) then

        call paractl(idmpopt,idnspol,                                   &
     &             iduniopt_uni,iddx,iddy,iddxiv,iddyiv,ni,nj,          &
     &             latsw,lonsw,latmin,lonmin,dlat,dlon,ipole,jpole,mlat)

      end if

! -----

! Calculate the 1 dimensional stretched z coordinates.

      if(fproc(1:3).eq.'dmp') then

        if(mype.eq.root) then

          if(mod(dmplev,10).eq.3) then

            call getz(iddz,idzsfc,nk,z1d)

            if(sthopt.eq.0) then

              call copy1d(1,nk,z1d,zsth)

            else

              call stretch(idsthopt,idzsfc,iddzmin,idlayer1,idlayer2,   &
     &                     z1d(nk-2),z1d(nk-1),nk,zsth,dzsth)

            end if

          end if

        end if

      end if

! -----

! Calculate the constant height at scalar points.

!$omp parallel default(shared)

      if(fproc(1:3).eq.'dmp') then

        if(mype.eq.root) then

          if(mod(dmplev,10).eq.3) then

!$omp do schedule(runtime) private(k_sub)

            do k_sub=2,nk-2
              z1d(k_sub)=.5e0*(zsth(k_sub)+zsth(k_sub+1))
            end do

!$omp end do

          end if

        end if

      end if

!$omp end parallel

! -----

!! Create the GrADS control file.

! Initialize the character variables.

      if(mype.eq.root) then

        call inichar(chkfl)
        call inichar(ctlfl)

      end if

! -----

! Get the unit numbers.

      if(mype.eq.root) then

        call getunit(iochk)
        call getunit(ioctl)

      end if

! -----

! Open the dumped data checking file.

      if(mype.eq.root) then

        chkfl(1:ncexp)=exprim(1:ncexp)

        if(fproc(1:3).eq.'dmp') then

          ncchk=ncexp+13

          write(chkfl(ncexp+1:ncexp+13),'(a13)') 'dmp.check.txt'

          open(iochk,iostat=stat,err=100,                               &
     &         file=crsdir(1:nccrs)//chkfl(1:ncexp+13),                 &
     &         status='old',access='sequential',form='formatted',       &
     &         blank='null',position='rewind',action='read')

        else if(fproc(1:3).eq.'mon') then

          ncchk=ncexp+13

          write(chkfl(ncexp+1:ncexp+13),'(a13)') 'mon.check.txt'

          open(iochk,iostat=stat,err=100,                               &
     &         file=crsdir(1:nccrs)//chkfl(1:ncexp+13),                 &
     &         status='old',access='sequential',form='formatted',       &
     &         blank='null',position='rewind',action='read')

        else if(fproc(1:3).eq.'geo') then

          ncchk=ncexp+19

          write(chkfl(ncexp+1:ncexp+19),'(a19)') 'geography.check.txt'

          open(iochk,iostat=stat,err=100,                               &
     &         file=crsdir(1:nccrs)//chkfl(1:ncexp+19),                 &
     &         status='old',access='sequential',form='formatted',       &
     &         blank='null',position='rewind',action='read')

        end if

      else

        stat=0

      end if

  100 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outctl  ',6,'cont',1,'              ',14,iochk, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outctl  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('outctl  ',6,chkfl,ncchk,iochk,1,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read out the data from the dumped data checking file and count the
! united variables.

      varcnt=0

      if(mype.eq.root) then

        read(iochk,'(a)',iostat=stat,end=110,err=110)                   &
     &      (rcname(in),in=1,ncn)

        read(iochk,*,iostat=stat,end=110,err=110)                       &
     &      (riname(in),in=1,nin)

        read(iochk,*,iostat=stat,end=110,err=110)                       &
     &      (rrname(in),in=1,nrn)

        read(iochk,'(a75,i3)',iostat=stat,end=110,err=110)              &
     &       vcap(1:75),ncvc

        do_varsec: do

          varcnt=varcnt+1

          read(iochk,'(a75,i3)',iostat=stat,end=120,err=120)            &
     &         vcap(1:75),ncvc

  120     if(stat.ne.0) then

            varcnt=varcnt-1

            if(varcnt.gt.0) then
              stat=0
            end if

            exit do_varsec

          end if

        end do do_varsec

        if(stat.eq.0) then

          rewind(iochk,iostat=stat,err=110)

          read(iochk,'(a)',iostat=stat,end=110,err=110)                 &
     &        (rcname(in),in=1,ncn)

          read(iochk,*,iostat=stat,end=110,err=110)                     &
     &        (riname(in),in=1,nin)

          read(iochk,*,iostat=stat,end=110,err=110)                     &
     &        (rrname(in),in=1,nrn)

          read(iochk,'(a75,i3)',iostat=stat,end=110,err=110)            &
     &         vcap(1:75),ncvc

        end if

      else

        stat=0

      end if

  110 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outctl  ',6,'cont',3,'              ',14,iochk, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outctl  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('outctl  ',6,chkfl,108,iochk,3,0,0_i8)

        call chkfile('ctl',stat,ncn,nin,nrn,                            &
     &               cname,iname,rname,rcname,riname,rrname)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Open the GrADS control file.

      if(mype.eq.root) then

        ctlfl(1:ncexp)=exprim(1:ncexp)

        if(fproc(1:3).eq.'dmp') then

          ncctl=ncexp+7

          write(ctlfl(ncexp+1:ncexp+7),'(a7)') 'dmp.ctl'

          open(ioctl,iostat=stat,err=130,                               &
     &         file=crsdir(1:nccrs)//ctlfl(1:ncexp+7),                  &
     &         status='new',access='sequential',form='formatted',       &
     &         blank='null',action='write')

        else if(fproc(1:3).eq.'mon') then

          ncctl=ncexp+7

          write(ctlfl(ncexp+1:ncexp+7),'(a7)') 'mon.ctl'

          open(ioctl,iostat=stat,err=130,                               &
     &         file=crsdir(1:nccrs)//ctlfl(1:ncexp+7),                  &
     &         status='new',access='sequential',form='formatted',       &
     &         blank='null',action='write')

        else if(fproc(1:3).eq.'geo') then

          ncctl=ncexp+13

          write(ctlfl(ncexp+1:ncexp+13),'(a13)') 'geography.ctl'

          open(ioctl,iostat=stat,err=130,                               &
     &         file=crsdir(1:nccrs)//ctlfl(1:ncexp+13),                 &
     &         status='new',access='sequential',form='formatted',       &
     &         blank='null',action='write')

        end if

      else

        stat=0

      end if

  130 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outctl  ',6,'warn',1,'              ',14,ioctl, &
     &                 stat)

        end if

        call cpondpe

        go to 140

      end if

      if(mype.eq.root) then

        call outstd03('outctl  ',6,ctlfl,ncctl,ioctl,1,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read in the massage to the GrADS control file.

      if(mype.eq.root) then

        if(fproc(1:3).eq.'dmp') then

          iflitv=int(flitv_uni+.1e0)

          if(mod(iflitv,10).ne.0) then

            write(ioctl,'(a6,a,a21)',iostat=stat,err=150)               &
     &           'dset ^',exprim(1:ncexp),'dmp0000%y4.united.bin'

          else if(mod(iflitv,100).ne.0) then

            write(ioctl,'(a6,a,a21)',iostat=stat,err=150)               &
     &           'dset ^',exprim(1:ncexp),'dmp000%y40.united.bin'

            iflitv=iflitv/10

          else if(mod(iflitv,1000).ne.0) then

            write(ioctl,'(a6,a,a21)',iostat=stat,err=150)               &
     &           'dset ^',exprim(1:ncexp),'dmp00%y400.united.bin'

            iflitv=iflitv/100

          else if(mod(iflitv,10000).ne.0) then

            write(ioctl,'(a6,a,a21)',iostat=stat,err=150)               &
     &           'dset ^',exprim(1:ncexp),'dmp0%y4000.united.bin'

            iflitv=iflitv/1000

          else

            write(ioctl,'(a6,a,a21)',iostat=stat,err=150)               &
     &           'dset ^',exprim(1:ncexp),'dmp%y40000.united.bin'

            iflitv=iflitv/10000

          end if

        else if(fproc(1:3).eq.'mon') then

          iflitv=int(flitv_uni+.1e0)

          if(mod(iflitv,10).ne.0) then

            write(ioctl,'(a6,a,a21)',iostat=stat,err=150)               &
     &           'dset ^',exprim(1:ncexp),'mon0000%y4.united.bin'

          else if(mod(iflitv,100).ne.0) then

            write(ioctl,'(a6,a,a21)',iostat=stat,err=150)               &
     &           'dset ^',exprim(1:ncexp),'mon000%y40.united.bin'

            iflitv=iflitv/10

          else if(mod(iflitv,1000).ne.0) then

            write(ioctl,'(a6,a,a21)',iostat=stat,err=150)               &
     &           'dset ^',exprim(1:ncexp),'mon00%y400.united.bin'

            iflitv=iflitv/100

          else if(mod(iflitv,10000).ne.0) then

            write(ioctl,'(a6,a,a21)',iostat=stat,err=150)               &
     &           'dset ^',exprim(1:ncexp),'mon0%y4000.united.bin'

            iflitv=iflitv/1000

          else

            write(ioctl,'(a6,a,a21)',iostat=stat,err=150)               &
     &           'dset ^',exprim(1:ncexp),'mon%y40000.united.bin'

            iflitv=iflitv/10000

          end if

        else if(fproc(1:3).eq.'geo') then

          write(ioctl,'(a6,a,a20)',iostat=stat,err=150)                 &
     &         'dset ^',exprim(1:ncexp),'geography.united.bin'

        end if

        write(ioctl,'(a)',iostat=stat,err=150) vcap(1:ncvc)

        write(ioctl,'(a12)',iostat=stat,err=150) 'undef -1.e35'

        if(uniopt_uni.lt.-10) then

          write(ioctl,'(a5,i7,a34)',iostat=stat,err=150) 'xdef ',       &
     &         (ni-3)*nired*nisub,' linear      1.00000     1.0000000'

          write(ioctl,'(a5,i7,a34)',iostat=stat,err=150) 'ydef ',       &
     &         (nj-3)*njred*njsub,' linear      1.00000     1.0000000'

        else

          if(mpopt.eq.1) then

            if(abs(uniopt_uni).eq.2.or.abs(uniopt_uni).eq.4) then

              write(ioctl,'(a5,2i7)',iostat=stat,err=150,advance='no')  &
     &              'pdef ',(ni-3)*nired*nisub,(nj-3)*njred*njsub

            else

              write(ioctl,'(a5,2i7)',iostat=stat,err=150,advance='no')  &
     &              'pdef ',(ni-3)*nigrp*nisub,(nj-3)*njgrp*njsub

            end if

            if(nspol.eq.1) then

              write(ioctl,'(a7,2f11.3,f9.3,f10.6)',iostat=stat,err=150) &
     &                    '    nps',ipole,jpole,tlon,.001e0*dx

            else

              write(ioctl,'(a7,2f11.3,f9.3,f10.6)',iostat=stat,err=150) &
     &                    '    sps',ipole,jpole,tlon,.001e0*dx

            end if

            if(abs(uniopt_uni).eq.2.or.abs(uniopt_uni).eq.4) then

              write(ioctl,'(a5,i7,a9,f11.5,f14.7)',iostat=stat,err=150) &
     &              'xdef ',(ni-3)*nired*nisub+2,' linear  ',lonmin,dlon

              write(ioctl,'(a5,i7,a9,f11.5,f14.7)',iostat=stat,err=150) &
     &              'ydef ',(nj-3)*njred*njsub+2,' linear  ',latmin,dlat

            else

              write(ioctl,'(a5,i7,a9,f11.5,f14.7)',iostat=stat,err=150) &
     &              'xdef ',(ni-3)*nigrp*nisub+2,' linear  ',lonmin,dlon

              write(ioctl,'(a5,i7,a9,f11.5,f14.7)',iostat=stat,err=150) &
     &              'ydef ',(nj-3)*njgrp*njsub+2,' linear  ',latmin,dlat

            end if

          else if(mpopt.eq.2) then

            if(abs(uniopt_uni).eq.2.or.abs(uniopt_uni).eq.4) then

              write(ioctl,'(a5,2i7)',iostat=stat,err=150,advance='no')  &
     &              'pdef ',(ni-3)*nired*nisub,(nj-3)*njred*njsub

            else

              write(ioctl,'(a5,2i7)',iostat=stat,err=150,advance='no')  &
     &              'pdef ',(ni-3)*nigrp*nisub,(nj-3)*njgrp*njsub

            end if

            if(dmplev.lt.10) then

!ORIG         write(ioctl,'(a7)',                                       &
!ORIG&              iostat=stat,err=150,advance='no') '   lccr'

              write(ioctl,'(a7)',                                       &
     &              iostat=stat,err=150,advance='no') '    lcc'

            else

              write(ioctl,'(a7)',                                       &
     &              iostat=stat,err=150,advance='no') '    lcc'

            end if

            write(ioctl,'(2f11.5)',                                     &
     &            iostat=stat,err=150,advance='no') latsw,lonsw

            write(ioctl,'(a12,2f8.3,f9.3,2f10.3)',iostat=stat,err=150)  &
     &                  ' 1.000 1.000',tlat1,tlat2,tlon,dx,dy

            if(abs(uniopt_uni).eq.2.or.abs(uniopt_uni).eq.4) then

              write(ioctl,'(a5,i7,a9,f11.5,f14.7)',iostat=stat,err=150) &
     &              'xdef ',(ni-3)*nired*nisub+2,' linear  ',lonmin,dlon

              write(ioctl,'(a5,i7,a9,f11.5,f14.7)',iostat=stat,err=150) &
     &              'ydef ',(nj-3)*njred*njsub+2,' linear  ',latmin,dlat

            else

              write(ioctl,'(a5,i7,a9,f11.5,f14.7)',iostat=stat,err=150) &
     &              'xdef ',(ni-3)*nigrp*nisub+2,' linear  ',lonmin,dlon

              write(ioctl,'(a5,i7,a9,f11.5,f14.7)',iostat=stat,err=150) &
     &              'ydef ',(nj-3)*njgrp*njsub+2,' linear  ',latmin,dlat

            end if

          else if(mpopt.eq.3.or.mpopt.eq.13) then

           if(abs(uniopt_uni).eq.2.or.abs(uniopt_uni).eq.4) then

            write(ioctl,'(a5,i7,a9,f11.5,f14.7)',iostat=stat,err=150)   &
     &            'xdef ',(ni-3)*nired*nisub,' linear  ',lonmin,dlon

            write(ioctl,'(a5,i7,a7)',iostat=stat,err=150)               &
     &            'ydef ',(nj-3)*njred*njsub,' levels'

            icnt=0

            do j=(nj-3)*njsub*jsred+1,(nj-3)*njsub*(jsred+njred)

              icnt=icnt+1

              if(icnt.eq.1) then

                write(ioctl,'(a1)',                                     &
     &                iostat=stat,err=150,advance='no') ' '

              end if

              if(j.eq.nyrm3.or.icnt.eq.7) then

                write(ioctl,'(f11.3)',iostat=stat,err=150) mlat(j)

                icnt=0

              else

                write(ioctl,'(f11.3)',                                  &
     &                iostat=stat,err=150,advance='no') mlat(j)

              end if

            end do

           else

            write(ioctl,'(a5,i7,a9,f11.5,f14.7)',iostat=stat,err=150)   &
     &            'xdef ',(ni-3)*nigrp*nisub,' linear  ',lonmin,dlon

            write(ioctl,'(a5,i7,a7)',iostat=stat,err=150)               &
     &            'ydef ',(nj-3)*njgrp*njsub,' levels'

            icnt=0

            do j=1,(nj-3)*njgrp*njsub

              icnt=icnt+1

              if(icnt.eq.1) then

                write(ioctl,'(a1)',                                     &
     &                iostat=stat,err=150,advance='no') ' '

              end if

              if(j.eq.nym3.or.icnt.eq.7) then

                write(ioctl,'(f11.3)',iostat=stat,err=150) mlat(j)

                icnt=0

              else

                write(ioctl,'(f11.3)',                                  &
     &                iostat=stat,err=150,advance='no') mlat(j)

              end if

            end do

           end if

          else

           if(abs(uniopt_uni).eq.2.or.abs(uniopt_uni).eq.4) then

            write(ioctl,'(a5,i7,a9,f11.5,f14.7)',iostat=stat,err=150)   &
     &            'xdef ',(ni-3)*nired*nisub,' linear  ',lonmin,dlon

            write(ioctl,'(a5,i7,a9,f11.5,f14.7)',iostat=stat,err=150)   &
     &            'ydef ',(nj-3)*njred*njsub,' linear  ',latmin,dlat

           else

            write(ioctl,'(a5,i7,a9,f11.5,f14.7)',iostat=stat,err=150)   &
     &            'xdef ',(ni-3)*nigrp*nisub,' linear  ',lonmin,dlon

            write(ioctl,'(a5,i7,a9,f11.5,f14.7)',iostat=stat,err=150)   &
     &            'ydef ',(nj-3)*njgrp*njsub,' linear  ',latmin,dlat

           end if

          end if

        end if

        if(fproc(1:3).eq.'dmp') then

          if(uniopt_uni.lt.-20) then

            write(ioctl,'(a5,i7,a34)',iostat=stat,err=150)              &
     &            'zdef ',nk-3,' linear      1.00000     1.0000000'

          else

            if(mod(dmplev,10).eq.2) then

              write(ioctl,'(a5,i7,a9,f11.5,f14.7)',iostat=stat,err=150) &
     &              'zdef ',nk-3,' linear  ',.5e0*dz,dz

            else if(mod(dmplev,10).eq.3) then

              write(ioctl,'(a5,i7,a7)',iostat=stat,err=150)             &
     &              'zdef ',nk-3,' levels'

              icnt=0

              do k=2,nk-2

                icnt=icnt+1

                if(icnt.eq.1) then

                  write(ioctl,'(a1)',                                   &
     &                  iostat=stat,err=150,advance='no') ' '

                end if

                if(k.eq.nkm2.or.icnt.eq.7) then

                  write(ioctl,'(f11.3)',iostat=stat,err=150) z1d(k)

                  icnt=0

                else

                  write(ioctl,'(f11.3)',                                &
     &                  iostat=stat,err=150,advance='no') z1d(k)

                end if

              end do

            else

              write(ioctl,'(a5,i7,a34)',iostat=stat,err=150)            &
     &              'zdef ',nk-3,' linear      1.00000     1.0000000'

            end if

          end if

          write(ioctl,'(a5,i7,a27,i5,a2)',iostat=stat,err=150) 'tdef ', &
     &        nstp1-nstp0+1_i8,' linear     00:00Z01Jan0000',iflitv,'yr'

        else if(fproc(1:3).eq.'mon') then

          write(ioctl,'(a46)',iostat=stat,err=150)                      &
     &                'zdef       1 linear      1.00000     1.0000000'

          write(ioctl,'(a5,i7,a27,i5,a2)',iostat=stat,err=150) 'tdef ', &
     &        nstp1-nstp0+1_i8,' linear     00:00Z01Jan0000',iflitv,'yr'

        else if(fproc(1:3).eq.'geo') then

          write(ioctl,'(a46)',iostat=stat,err=150)                      &
     &                'zdef       1 linear      1.00000     1.0000000'

          write(ioctl,'(a46)',iostat=stat,err=150)                      &
     &                'tdef       1 linear     00:00Z01Jan0000    1yr'

        end if

        write(ioctl,'(a5,i7)',iostat=stat,err=150) 'vars ',varcnt

        do icnt=1,varcnt

          read(iochk,'(a75,i3)') vcap(1:75),ncvc

          write(ioctl,'(a)',iostat=stat,err=150) vcap(1:ncvc)

        end do

        write(ioctl,'(a7)',iostat=stat,err=150) 'endvars'

      else

        stat=0

      end if

  150 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outctl  ',6,'cont',3,'              ',14,ioctl, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outctl  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('outctl  ',6,ctlfl,108,ioctl,4,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Close the GrADS control file.

      if(mype.eq.root) then

        close(ioctl,iostat=stat,err=160,status='keep')

      else

        stat=0

      end if

  160 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outctl  ',6,'cont',2,'              ',14,ioctl, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outctl  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('outctl  ',6,ctlfl,108,ioctl,2,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Close the dumped data checking file.

      if(mype.eq.root) then

        close(iochk,iostat=stat,err=170,status='keep')

      else

        stat=0

      end if

  170 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outctl  ',6,'cont',2,'              ',14,iochk, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outctl  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('outctl  ',6,chkfl,108,iochk,2,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit numbers.

  140 if(mype.eq.root) then

        call putunit(ioctl)
        call putunit(iochk)

      end if

! -----

!! -----

      end subroutine s_outctl

!-----7--------------------------------------------------------------7--

      end module m_outctl
