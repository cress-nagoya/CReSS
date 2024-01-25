!***********************************************************************
      module m_outgeo
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/03/13
!     Modification: 2001/04/15, 2001/05/29, 2001/11/20, 2002/04/02,
!                   2002/06/18, 2002/07/03, 2002/07/15, 2003/01/04,
!                   2003/03/28, 2003/04/30, 2003/05/19, 2003/06/27,
!                   2003/12/12, 2004/04/15, 2004/05/31, 2004/06/10,
!                   2004/08/20, 2004/08/31, 2004/09/25, 2005/01/14,
!                   2005/02/10, 2006/09/21, 2006/11/06, 2006/12/04,
!                   2007/01/05, 2007/01/20, 2007/04/11, 2007/05/14,
!                   2007/08/24, 2008/01/11, 2008/04/17, 2008/05/02,
!                   2008/08/25, 2008/10/10, 2009/01/30, 2009/02/27,
!                   2011/08/09, 2011/09/22, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     open and read in data to the geography file.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_chkopen
      use m_chkstd
      use m_comcapt
      use m_comdmp
      use m_comkind
      use m_commpi
      use m_comname
      use m_cpondpe
      use m_destroy
      use m_getcname
      use m_getiname
      use m_getunit
      use m_inichar
      use m_outcap
      use m_outstd03
      use m_outstd10
      use m_putunit

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: outgeo, s_outgeo

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outgeo

        module procedure s_outgeo

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_outgeo(fpexprim,fpcrsdir,fpncexp,fpnccrs,fpwlngth,   &
     &                    fpmfcopt,fpcoropt,fpsfcopt,fpdmpfmt,          &
     &                    ni,nj,ht,lat,lon,mf,fc,land)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpexprim
                       ! Formal parameter of unique index of exprim

      integer, intent(in) :: fpcrsdir
                       ! Formal parameter of unique index of crsdir

      integer, intent(in) :: fpncexp
                       ! Formal parameter of unique index of ncexp

      integer, intent(in) :: fpnccrs
                       ! Formal parameter of unique index of nccrs

      integer, intent(in) :: fpwlngth
                       ! Formal parameter of unique index of wlngth

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

      integer, intent(in) :: fpcoropt
                       ! Formal parameter of unique index of coropt

      integer, intent(in) :: fpsfcopt
                       ! Formal parameter of unique index of sfcopt

      integer, intent(in) :: fpdmpfmt
                       ! Formal parameter of unique index of dmpfmt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

      real, intent(in) :: ht(0:ni+1,0:nj+1)
                       ! Terrain height

      real, intent(in) :: lat(0:ni+1,0:nj+1)
                       ! Latitude

      real, intent(in) :: lon(0:ni+1,0:nj+1)
                       ! Longitude

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: fc(0:ni+1,0:nj+1,1:2)
                       ! 0.25 x Coriolis parameters

! Internal shared variables

      character(len=108) exprim
                       ! Optional run name

      character(len=108) crsdir
                       ! User specified directory for CReSS files

      character(len=108) geofl
                       ! Opened file name

      integer ncexp    ! Number of character of exprim
      integer nccrs    ! Number of character of crsdir

      integer wlngth   ! Word length of direct access file

      integer mfcopt   ! Option for map scale factor
      integer coropt   ! Option for coriolis force
      integer sfcopt   ! Option for surface physics

      integer dmpfmt   ! Option for dumped file format

      integer in       ! Namelist table index

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      integer ncfl     ! Number of character of geography file name

      integer ncend    ! Number of output character for endian option

      integer extpe    ! Control flag of file name extension

      integer recgeo   ! Current record number of geography file

      integer iogeo    ! Unit number of geography file

      integer siz      ! Record length of geography file

      integer stat     ! Runtime status

      integer broot    ! Broadcasting root

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
      call getiname(fpwlngth,wlngth)
      call getiname(fpmfcopt,mfcopt)
      call getiname(fpcoropt,coropt)
      call getiname(fpsfcopt,sfcopt)
      call getiname(fpdmpfmt,dmpfmt)

! -----

!! Open and read in the data to the geography data checking file.

! Initialize the character variable.

      if(mype.eq.root) then

        call inichar(geofl)

      end if

! -----

! Get the unit number.

      if(mype.eq.root) then

        call getunit(iogeo)

      end if

! -----

! Open the geography data checking file.

      if(mype.eq.root) then

        geofl(1:ncexp)=exprim(1:ncexp)

        write(geofl(ncexp+1:ncexp+19),'(a19)') 'geography.check.txt'

        open(iogeo,iostat=stat,err=100,                                 &
     &       file=crsdir(1:nccrs)//geofl(1:ncexp+19),                   &
     &       status='new',access='sequential',form='formatted',         &
     &       blank='null',action='write')

  100   if(stat.eq.0) then

          ncfl=ncexp+19

        else

          ncfl=ncexp+24

          write(geofl(ncexp+20:ncexp+24),'(a5)') '.swap'

          open(iogeo,iostat=stat,err=110,                               &
     &         file=crsdir(1:nccrs)//geofl(1:ncexp+24),                 &
     &         status='new',access='sequential',form='formatted',       &
     &         blank='null',action='write')

        end if

      else

        stat=0

      end if

  110 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outgeo  ',6,'cont',1,'              ',14,iogeo, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outgeo  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('outgeo  ',6,geofl,ncfl,iogeo,1,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read in the data to the geography data checking file.

      if(mype.eq.root) then

        write(iogeo,'(a)',iostat=stat,err=120)                          &
     &       (cname(in),in=1,ncn)

        write(iogeo,*,iostat=stat,err=120)                              &
     &       (iname(in),in=1,nin)

        write(iogeo,*,iostat=stat,err=120)                              &
     &       (rname(in),in=1,nrn)

        if(border(1:7).eq.'unknown') then

          ncend=16

          write(iogeo,'(a30,a30,a15,i3)',iostat=stat,err=120)           &
     &                'options template              ',                 &
     &                '                              ',                 &
     &                '               ',ncend

        else if(border(1:3).eq.'big') then

          ncend=27

          write(iogeo,'(a30,a30,a15,i3)',iostat=stat,err=120)           &
     &                'options template big_endian   ',                 &
     &                '                              ',                 &
     &                '               ',ncend

        else if(border(1:6).eq.'little') then

          ncend=30

          write(iogeo,'(a30,a30,a15,i3)',iostat=stat,err=120)           &
     &                'options template little_endian',                 &
     &                '                              ',                 &
     &                '               ',ncend

        end if

      else

        stat=0

      end if

  120 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outgeo  ',6,'cont',3,'              ',14,iogeo, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outgeo  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      call outcap('geo','ht    ',capt(66),ncpt(66),iogeo,3)
      call outcap('geo','alat  ',capt(67),ncpt(67),iogeo,3)
      call outcap('geo','alon  ',capt(68),ncpt(68),iogeo,3)

      if(mfcopt.eq.1) then

        call outcap('geo','map   ',capt(69),ncpt(69),iogeo,3)

      end if

      if(coropt.ge.1) then

        call outcap('geo','fs    ',capt(70),ncpt(70),iogeo,3)

        if(coropt.eq.2) then

          call outcap('geo','fc    ',capt(71),ncpt(71),iogeo,3)

        end if

      end if

      if(sfcopt.ge.1) then

        call outcap('geo','land  ',capt(72),ncpt(72),iogeo,3)

      end if

      if(mype.eq.root) then

        call outstd03('outgeo  ',6,geofl,108,iogeo,4,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Close the geography data checking file.

      if(mype.eq.root) then

        close(iogeo,iostat=stat,err=130,status='keep')

      else

        stat=0

      end if

  130 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outgeo  ',6,'cont',2,'              ',14,iogeo, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outgeo  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('outgeo  ',6,geofl,108,iogeo,2,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit number.

      if(mype.eq.root) then

        call putunit(iogeo)

      end if

! -----

!! -----

!! Open and read in data to the geography file.

! Initialize the character variable.

      call inichar(geofl)

! -----

! Get the unit number.

      call getunit(iogeo)

! -----

! Open the geography file.

      geofl(1:ncexp)=exprim(1:ncexp)

      if(ngrp.eq.1) then

        ncfl=ncexp+20

        write(geofl(ncexp+1:ncexp+16),'(a12,i4.4)')                     &
     &                   'geography.pe',mysub

      else

        ncfl=ncexp+29

        write(geofl(ncexp+1:ncexp+25),'(a13,i4.4,a4,i4.4)')             &
     &                  'geography.grp',mygrp,'-sub',mysub

      end if

      if(dmpfmt.eq.1) then

        write(geofl(ncfl-3:ncfl),'(a4)') '.txt'

        open(iogeo,iostat=stat,err=140,                                 &
     &       file=crsdir(1:nccrs)//geofl(1:ncfl),                       &
     &       status='new',access='sequential',form='formatted',         &
     &       blank='null',action='write')

  140   if(stat.eq.0) then

          extpe=0

        else

          extpe=1

          ncfl=ncfl+5

          write(geofl(ncfl-4:ncfl),'(a5)') '.swap'

          open(iogeo,iostat=stat,err=150,                               &
     &         file=crsdir(1:nccrs)//geofl(1:ncfl),                     &
     &         status='new',access='sequential',form='formatted',       &
     &         blank='null',action='write')

        end if

      else if(dmpfmt.eq.2) then

        siz=(ni-3)*(nj-3)*wlngth

        write(geofl(ncfl-3:ncfl),'(a4)') '.bin'

        open(iogeo,iostat=stat,err=160,                                 &
     &       file=crsdir(1:nccrs)//geofl(1:ncfl),                       &
     &       status='new',access='direct',form='unformatted',           &
     &       recl=siz,action='write')

  160   if(stat.eq.0) then

          extpe=0

        else

          extpe=1

          ncfl=ncfl+5

          write(geofl(ncfl-4:ncfl),'(a5)') '.swap'

          open(iogeo,iostat=stat,err=150,                               &
     &         file=crsdir(1:nccrs)//geofl(1:ncfl),                     &
     &         status='new',access='direct',form='unformatted',         &
     &         recl=siz,action='write')

        end if

      end if

  150 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outgeo  ',6,'cont',1,'              ',14,iogeo, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outgeo  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      call chkopen(extpe,iogeo)

      if(mype.eq.stat-1) then

        if(fpara(1:5).eq.'multi') then

          if(ngrp.eq.1) then

            write(geofl(ncexp+13:ncexp+16),'(a4)') 'XXXX'

          else

            write(geofl(ncexp+14:ncexp+17),'(a4)') 'XXXX'
            write(geofl(ncexp+22:ncexp+25),'(a4)') 'YYYY'

          end if

        end if

        call outstd03('outgeo  ',6,geofl,ncfl,iogeo,1,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read in data to the geography file.

      if(dmpfmt.eq.1) then

        write(iogeo,*,iostat=stat,err=170)                              &
     &              ((ht(i,j),i=2,ni-2),j=2,nj-2)

        write(iogeo,*,iostat=stat,err=170)                              &
     &              ((lat(i,j),i=2,ni-2),j=2,nj-2)

        write(iogeo,*,iostat=stat,err=170)                              &
     &              ((lon(i,j),i=2,ni-2),j=2,nj-2)

        if(mfcopt.eq.1) then

          write(iogeo,*,iostat=stat,err=170)                            &
     &                ((mf(i,j),i=2,ni-2),j=2,nj-2)

        end if

        if(coropt.ge.1) then

          write(iogeo,*,iostat=stat,err=170)                            &
     &                ((4.e0*fc(i,j,1),i=2,ni-2),j=2,nj-2)

          if(coropt.eq.2) then

            write(iogeo,*,iostat=stat,err=170)                          &
     &                  ((4.e0*fc(i,j,2),i=2,ni-2),j=2,nj-2)

          end if

        end if

        if(sfcopt.ge.1) then

          write(iogeo,*,iostat=stat,err=170)                            &
     &                ((real(land(i,j)),i=2,ni-2),j=2,nj-2)

        end if

      else if(dmpfmt.eq.2) then

        recgeo=3

        write(iogeo,rec=1,iostat=stat,err=170)                          &
     &       ((ht(i,j),i=2,ni-2),j=2,nj-2)

        write(iogeo,rec=2,iostat=stat,err=170)                          &
     &       ((lat(i,j),i=2,ni-2),j=2,nj-2)

        write(iogeo,rec=3,iostat=stat,err=170)                          &
     &       ((lon(i,j),i=2,ni-2),j=2,nj-2)

        if(mfcopt.eq.1) then

          recgeo=recgeo+1

          write(iogeo,rec=recgeo,iostat=stat,err=170)                   &
     &         ((mf(i,j),i=2,ni-2),j=2,nj-2)

        end if

        if(coropt.ge.1) then

          recgeo=recgeo+1

          write(iogeo,rec=recgeo,iostat=stat,err=170)                   &
     &         ((4.e0*fc(i,j,1),i=2,ni-2),j=2,nj-2)

          if(coropt.eq.2) then

            recgeo=recgeo+1

            write(iogeo,rec=recgeo,iostat=stat,err=170)                 &
     &           ((4.e0*fc(i,j,2),i=2,ni-2),j=2,nj-2)

          end if

        end if

        if(sfcopt.ge.1) then

          recgeo=recgeo+1

          write(iogeo,rec=recgeo,iostat=stat,err=170)                   &
     &         ((real(land(i,j)),i=2,ni-2),j=2,nj-2)

        end if

      end if

  170 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outgeo  ',6,'cont',4,'              ',14,iogeo, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outgeo  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd10('ht    ',2,capt(66),ncpt(66),1)
        call outstd10('alat  ',4,capt(67),ncpt(67),2)
        call outstd10('alon  ',4,capt(68),ncpt(68),3)

        if(mfcopt.eq.1) then

          call outstd10('map   ',3,capt(69),ncpt(69),4)

        end if

        if(coropt.ge.1) then

          call outstd10('fs    ',2,capt(70),ncpt(70),5)

          if(coropt.eq.2) then

            call outstd10('fc    ',2,capt(71),ncpt(71),6)

          end if

        end if

        if(sfcopt.ge.1) then

          call outstd10('rland ',5,capt(72),ncpt(72),7)

        end if

      end if

      call cpondpe

! -----

! Close the geography file.

      close(iogeo,iostat=stat,err=180,status='keep')

  180 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outgeo  ',6,'cont',2,'              ',14,iogeo, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outgeo  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('outgeo  ',6,geofl,108,iogeo,2,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit number.

      call putunit(iogeo)

! -----

!! -----

      end subroutine s_outgeo

!-----7--------------------------------------------------------------7--

      end module m_outgeo
