!***********************************************************************
      module m_rdsfc
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/07/03
!     Modification: 2002/07/15, 2003/03/28, 2003/04/30, 2003/05/19,
!                   2003/06/27, 2003/07/15, 2003/12/12, 2004/01/09,
!                   2004/05/31, 2004/08/01, 2004/08/20, 2004/09/01,
!                   2005/01/14, 2005/02/10, 2006/09/21, 2006/12/04,
!                   2007/01/05, 2007/01/20, 2007/08/24, 2008/05/02,
!                   2008/08/25, 2008/10/10, 2009/01/30, 2009/02/27,
!                   2011/11/10, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     open and read out the data from the interpolated surface file.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkfile
      use m_chkerr
      use m_chkstd
      use m_comkind
      use m_commpi
      use m_comname
      use m_cpondpe
      use m_destroy
      use m_getcname
      use m_getiname
      use m_getunit
      use m_inichar
      use m_outstd03
      use m_putunit

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: rdsfc, s_rdsfc

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rdsfc

        module procedure s_rdsfc

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

!     none

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_rdsfc(fpexprim,fpcrsdir,fpsfcdat,fpncexp,fpnccrs,    &
     &                   fpwlngth,ni,nj,land,albe,beta,z0m,z0h,         &
     &                   cap,nuu,kai)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpexprim
                       ! Formal parameter of unique index of exprim

      integer, intent(in) :: fpcrsdir
                       ! Formal parameter of unique index of crsdir

      integer, intent(in) :: fpsfcdat
                       ! Formal parameter of unique index of sfcdat

      integer, intent(in) :: fpncexp
                       ! Formal parameter of unique index of ncexp

      integer, intent(in) :: fpnccrs
                       ! Formal parameter of unique index of nccrs

      integer, intent(in) :: fpwlngth
                       ! Formal parameter of unique index of wlngth

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

! Output variables

      integer, intent(out) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

      real, intent(out) :: albe(0:ni+1,0:nj+1)
                       ! Albedo

      real, intent(out) :: beta(0:ni+1,0:nj+1)
                       ! Evapotranspiration efficiency

      real, intent(out) :: z0m(0:ni+1,0:nj+1)
                       ! Roughness length for velocity

      real, intent(out) :: z0h(0:ni+1,0:nj+1)
                       ! Roughness length for scalar

      real, intent(out) :: cap(0:ni+1,0:nj+1)
                       ! Thermal capacity

      real, intent(out) :: nuu(0:ni+1,0:nj+1)
                       ! Thermal diffusivity

      real, intent(out) :: kai(0:ni+1,0:nj+1)
                       ! Sea ice distribution

! Internal shared variables

      character(len=108) exprim
                       ! Optional run name

      character(len=108) crsdir
                       ! User specified directory for CReSS files

      character(len=108) sfcdat
                       ! Control flag of input surface data type

      character(len=108) sfcfl
                       ! Opened file name

      integer ncexp    ! Number of character of exprim
      integer nccrs    ! Number of character of crsdir

      integer wlngth   ! Word length of direct access file

      integer in       ! Namelist table index

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      integer ncfl     ! Number of character of surface file

      integer recsfc   ! Current record number of surface file

      integer iosfc    ! Unit number of surface file

      integer siz      ! Record length of surface file

      integer stat     ! Runtime status

      integer broot    ! Broadcasting root

!-----7--------------------------------------------------------------7--

! Initialize the character variables.

      call inichar(exprim)
      call inichar(crsdir)
      call inichar(sfcdat)

! -----

! Get the required namelist variables.

      call getcname(fpexprim,exprim)
      call getcname(fpcrsdir,crsdir)
      call getcname(fpsfcdat,sfcdat)
      call getiname(fpncexp,ncexp)
      call getiname(fpnccrs,nccrs)
      call getiname(fpwlngth,wlngth)

! -----

!! Open and read out the data from the interpolated surface checking
!! file.

! Initialize the character variable.

      if(mype.eq.root) then

        call inichar(sfcfl)

      end if

! -----

! Get the unit number.

      if(mype.eq.root) then

        call getunit(iosfc)

      end if

! -----

! Open the interpolated surface checking file.

      if(mype.eq.root) then

        sfcfl(1:ncexp)=exprim(1:ncexp)

        write(sfcfl(ncexp+1:ncexp+17),'(a17)') 'surface.check.txt'

        open(iosfc,iostat=stat,err=100,                                 &
     &       file=crsdir(1:nccrs)//sfcfl(1:ncexp+17),                   &
     &       status='old',access='sequential',form='formatted',         &
     &       blank='null',position='rewind',action='read')

      else

        stat=0

      end if

  100 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdsfc   ',5,'cont',1,'              ',14,iosfc, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdsfc   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

       call outstd03('rdsfc   ',5,sfcfl,ncexp+17,iosfc,1,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read out the data from the interpolated surface checking file.

      if(mype.eq.root) then

        read(iosfc,'(a)',iostat=stat,end=110,err=110)                   &
     &      (rcname(in),in=1,ncn)

        read(iosfc,*,iostat=stat,end=110,err=110)                       &
     &      (riname(in),in=1,nin)

        read(iosfc,*,iostat=stat,end=110,err=110)                       &
     &      (rrname(in),in=1,nrn)

      else

        stat=0

      end if

  110 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdsfc   ',5,'cont',3,'              ',14,iosfc, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdsfc   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('rdsfc   ',5,sfcfl,108,iosfc,3,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Check the interpolated surface checking file.

      if(mype.eq.root) then

        call chkfile('sfc',stat,ncn,nin,nrn,                            &
     &               cname,iname,rname,rcname,riname,rrname)

      else

        stat=0

      end if

      call chkerr(stat)

      if(stat.lt.0) then

        call destroy('chkfile ',7,'stop',1001,'              ',14,10,   &
     &               stat)

      end if

! -----

! Close the interpolated surface checking file.

      if(mype.eq.root) then

        close(iosfc,iostat=stat,err=120,status='keep')

      else

        stat=0

      end if

  120 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdsfc   ',5,'cont',2,'              ',14,iosfc, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdsfc   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('rdsfc   ',5,sfcfl,108,iosfc,2,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit number.

      if(mype.eq.root) then

        call putunit(iosfc)

      end if

! -----

!! -----

!! Open and read out the data from the interpolated surface file.

! Initialize the character variable.

      call inichar(sfcfl)

! -----

! Get the unit number.

      call getunit(iosfc)

! -----

! Open the interpolated surface file.

      siz=(ni+2)*(nj+2)*wlngth

      sfcfl(1:ncexp)=exprim(1:ncexp)

      if(ngrp.eq.1) then

        ncfl=ncexp+18

        write(sfcfl(ncexp+1:ncexp+18),'(a10,i4.4,a4)')                  &
     &                     'surface.pe',mysub,'.bin'

      else

        ncfl=ncexp+27

        write(sfcfl(ncexp+1:ncexp+27),'(a11,i4.4,a4,i4.4,a4)')          &
     &                    'surface.grp',mygrp,'-sub',mysub,'.bin'

      end if

      open(iosfc,iostat=stat,err=130,                                   &
     &     file=crsdir(1:nccrs)//sfcfl(1:ncfl),                         &
     &     status='old',access='direct',form='unformatted',             &
     &     recl=siz,action='read')

  130 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdsfc   ',5,'cont',1,'              ',14,iosfc, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdsfc   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        if(fpara(1:5).eq.'multi') then

          if(ngrp.eq.1) then

            write(sfcfl(ncexp+11:ncexp+14),'(a4)') 'XXXX'

          else

            write(sfcfl(ncexp+12:ncexp+15),'(a4)') 'XXXX'
            write(sfcfl(ncexp+20:ncexp+23),'(a4)') 'YYYY'

          end if

        end if

        call outstd03('rdsfc   ',5,sfcfl,ncfl,iosfc,1,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read out the data from the interpolated surface file.

      recsfc=0

      if(sfcdat(1:1).eq.'o') then

        recsfc=recsfc+1

        read(iosfc,rec=recsfc,iostat=stat,err=140)                      &
     &      ((land(i,j),i=0,ni+1),j=0,nj+1)

        recsfc=recsfc+1

        read(iosfc,rec=recsfc,iostat=stat,err=140)                      &
     &      ((albe(i,j),i=0,ni+1),j=0,nj+1)

        recsfc=recsfc+1

        read(iosfc,rec=recsfc,iostat=stat,err=140)                      &
     &      ((beta(i,j),i=0,ni+1),j=0,nj+1)

        recsfc=recsfc+1

        read(iosfc,rec=recsfc,iostat=stat,err=140)                      &
     &      ((z0m(i,j),i=0,ni+1),j=0,nj+1)

        recsfc=recsfc+1

        read(iosfc,rec=recsfc,iostat=stat,err=140)                      &
     &      ((z0h(i,j),i=0,ni+1),j=0,nj+1)

        recsfc=recsfc+1

        read(iosfc,rec=recsfc,iostat=stat,err=140)                      &
     &      ((cap(i,j),i=0,ni+1),j=0,nj+1)

        recsfc=recsfc+1

        read(iosfc,rec=recsfc,iostat=stat,err=140)                      &
     &      ((nuu(i,j),i=0,ni+1),j=0,nj+1)

      end if

      if(sfcdat(3:3).eq.'o') then

        recsfc=recsfc+1

        read(iosfc,rec=recsfc,iostat=stat,err=140)                      &
     &      ((kai(i,j),i=0,ni+1),j=0,nj+1)

      end if

  140 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdsfc   ',5,'cont',3,'              ',14,iosfc, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdsfc   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('rdsfc   ',5,sfcfl,108,iosfc,3,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Close the interpolated surface file.

      close(iosfc,iostat=stat,err=150,status='keep')

  150 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdsfc   ',5,'cont',2,'              ',14,iosfc, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdsfc   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('rdsfc   ',5,sfcfl,108,iosfc,2,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit number.

      call putunit(iosfc)

! -----

!! -----

      end subroutine s_rdsfc

!-----7--------------------------------------------------------------7--

      end module m_rdsfc
