!***********************************************************************
      module m_outrst
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2005/07/01
!     Modification: 2005/08/05, 2005/10/05, 2006/01/10, 2006/02/13,
!                   2006/04/03, 2006/09/21, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/05/07, 2007/05/14, 2007/07/30,
!                   2007/08/24, 2007/11/26, 2008/03/12, 2008/05/02,
!                   2008/07/01, 2008/08/25, 2008/10/10, 2008/12/11,
!                   2009/01/30, 2009/02/27, 2009/12/05, 2011/08/18,
!                   2011/09/22, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     open and read in the data to the restructed restart file.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_chkstd
      use m_comindx
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

      public :: outrst, s_outrst

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outrst

        module procedure s_outrst

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_outrst(fpexprim,fpcrsdir,fpdmpvar,fpncexp,fpnccrs,   &
     &               fpwbc,fpebc,fpsbc,fpnbc,fpsfcopt,fpadvopt,fpcphopt,&
     &               fpqcgopt,fpaslopt,fptrkopt,fptubopt,fpdiaopt,      &
     &               fmois,ctime,nk,nqw,nnw,nqi,nni,nqa,nund,           &
     &               ni_rst,nj_rst,ubr_rst,vbr_rst,pbr_rst,ptbr_rst,    &
     &               qvbr_rst,u_rst,up_rst,v_rst,vp_rst,w_rst,wp_rst,   &
     &               pp_rst,ppp_rst,ptp_rst,ptpp_rst,qv_rst,qvp_rst,    &
     &               qwtr_rst,qwtrp_rst,nwtr_rst,nwtrp_rst,             &
     &               qice_rst,qicep_rst,nice_rst,nicep_rst,             &
     &               qcwtr_rst,qcwtrp_rst,qcice_rst,qcicep_rst,         &
     &               qasl_rst,qaslp_rst,qt_rst,qtp_rst,tke_rst,tkep_rst,&
     &               ucpx_rst,ucpy_rst,vcpx_rst,vcpy_rst,               &
     &               wcpx_rst,wcpy_rst,pcpx_rst,pcpy_rst,               &
     &               ptcpx_rst,ptcpy_rst,qvcpx_rst,qvcpy_rst,           &
     &               qwcpx_rst,qwcpy_rst,nwcpx_rst,nwcpy_rst,           &
     &               qicpx_rst,qicpy_rst,nicpx_rst,nicpy_rst,           &
     &               qcwcpx_rst,qcwcpy_rst,qcicpx_rst,qcicpy_rst,       &
     &               qacpx_rst,qacpy_rst,qtcpx_rst,qtcpy_rst,           &
     &               tkecpx_rst,tkecpy_rst,maxvl_rst,                   &
     &               prwtr_rst,price_rst,pdia_rst,z0m_rst,z0h_rst,      &
     &               tund_rst,tundp_rst)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: fpexprim
                       ! Formal parameter of unique index of exprim

      integer, intent(in) :: fpcrsdir
                       ! Formal parameter of unique index of crsdir

      integer, intent(in) :: fpdmpvar
                       ! Formal parameter of unique index of dmpvar

      integer, intent(in) :: fpncexp
                       ! Formal parameter of unique index of ncexp

      integer, intent(in) :: fpnccrs
                       ! Formal parameter of unique index of nccrs

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpsbc
                       ! Formal parameter of unique index of sbc

      integer, intent(in) :: fpnbc
                       ! Formal parameter of unique index of nbc

      integer, intent(in) :: fpsfcopt
                       ! Formal parameter of unique index of sfcopt

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fpqcgopt
                       ! Formal parameter of unique index of qcgopt

      integer, intent(in) :: fpaslopt
                       ! Formal parameter of unique index of aslopt

      integer, intent(in) :: fptrkopt
                       ! Formal parameter of unique index of trkopt

      integer, intent(in) :: fptubopt
                       ! Formal parameter of unique index of tubopt

      integer, intent(in) :: fpdiaopt
                       ! Formal parameter of unique index of diaopt

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

      integer, intent(in) :: nnw
                       ! Number of categories of water concentrations

      integer, intent(in) :: nqi
                       ! Number of categories of ice hydrometeor

      integer, intent(in) :: nni
                       ! Number of categories of ice concentrations

      integer, intent(in) :: nqa(0:4)
                       ! Number of types of aerosol

      integer, intent(in) :: nund
                       ! Number of soil and sea layers

      integer, intent(in) :: ni_rst
                       ! Restructed files dimension in x direction

      integer, intent(in) :: nj_rst
                       ! Restructed files dimension in y direction

      real, intent(in) :: ubr_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! ubr in restructed domain

      real, intent(in) :: vbr_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! vbr in restructed domain

      real, intent(in) :: pbr_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! pbr in restructed domain

      real, intent(in) :: ptbr_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! ptbr in restructed domain

      real, intent(in) :: qvbr_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! qvbr in restructed domain

      real, intent(in) :: u_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! u in restructed domain

      real, intent(in) :: up_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! up in restructed domain

      real, intent(in) :: v_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! v in restructed domain

      real, intent(in) :: vp_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! vp in restructed domain

      real, intent(in) :: w_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! w in restructed domain

      real, intent(in) :: wp_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! wp in restructed domain

      real, intent(in) :: pp_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! pp in restructed domain

      real, intent(in) :: ppp_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! ppp in restructed domain

      real, intent(in) :: ptp_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! ptp in restructed domain

      real, intent(in) :: ptpp_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! ptpp in restructed domain

      real, intent(in) :: qv_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! qv in restructed domain

      real, intent(in) :: qvp_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! qvp in restructed domain

      real, intent(in) :: qwtr_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqw)
                       ! qwtr in restructed domain

      real, intent(in) :: qwtrp_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqw)
                       ! qwtrp in restructed domain

      real, intent(in) :: nwtr_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nnw)
                       ! nwtr in restructed domain

      real, intent(in) :: nwtrp_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nnw)
                       ! nwtrp in restructed domain

      real, intent(in) :: qice_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqi)
                       ! qice in restructed domain

      real, intent(in) :: qicep_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqi)
                       ! qicep in restructed domain

      real, intent(in) :: nice_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nni)
                       ! nice in restructed domain

      real, intent(in) :: nicep_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nni)
                       ! nicep in restructed domain

      real, intent(in) :: qcwtr_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqw)
                       ! qcwtr in restructed domain

      real, intent(in) :: qcwtrp_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqw)
                       ! qcwtrp in restructed domain

      real, intent(in) :: qcice_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqi)
                       ! qcice in restructed domain

      real, intent(in) :: qcicep_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqi)
                       ! qcicep in restructed domain

      real, intent(in) :: qasl_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqa(0))
                       ! qasl in restructed domain

      real, intent(in) :: qaslp_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqa(0))
                       ! qaslp in restructed domain

      real, intent(in) :: qt_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! qt in restructed domain

      real, intent(in) :: qtp_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! qtp in restructed domain

      real, intent(in) :: tke_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! tke in restructed domain

      real, intent(in) :: tkep_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! tkep in restructed domain

      real, intent(in) :: ucpx_rst(1:nj_rst,1:nk,1:2)
                       ! ucpx in restructed domain

      real, intent(in) :: ucpy_rst(1:ni_rst,1:nk,1:2)
                       ! ucpy in restructed domain

      real, intent(in) :: vcpx_rst(1:nj_rst,1:nk,1:2)
                       ! vcpx in restructed domain

      real, intent(in) :: vcpy_rst(1:ni_rst,1:nk,1:2)
                       ! vcpy in restructed domain

      real, intent(in) :: wcpx_rst(1:nj_rst,1:nk,1:2)
                       ! wcpx in restructed domain

      real, intent(in) :: wcpy_rst(1:ni_rst,1:nk,1:2)
                       ! wcpy in restructed domain

      real, intent(in) :: pcpx_rst(1:nj_rst,1:nk,1:2)
                       ! pcpx in restructed domain

      real, intent(in) :: pcpy_rst(1:ni_rst,1:nk,1:2)
                       ! pcpy in restructed domain

      real, intent(in) :: ptcpx_rst(1:nj_rst,1:nk,1:2)
                       ! ptcpx in restructed domain

      real, intent(in) :: ptcpy_rst(1:ni_rst,1:nk,1:2)
                       ! ptcpy in restructed domain

      real, intent(in) :: qvcpx_rst(1:nj_rst,1:nk,1:2)
                       ! qvcpx in restructed domain

      real, intent(in) :: qvcpy_rst(1:ni_rst,1:nk,1:2)
                       ! qvcpy in restructed domain

      real, intent(in) :: qwcpx_rst(1:nj_rst,1:nk,1:2,1:nqw)
                       ! qwcpx in restructed domain

      real, intent(in) :: qwcpy_rst(1:ni_rst,1:nk,1:2,1:nqw)
                       ! qwcpy in restructed domain

      real, intent(in) :: nwcpx_rst(1:nj_rst,1:nk,1:2,1:nnw)
                       ! nwcpx in restructed domain

      real, intent(in) :: nwcpy_rst(1:ni_rst,1:nk,1:2,1:nnw)
                       ! nwcpy in restructed domain

      real, intent(in) :: qicpx_rst(1:nj_rst,1:nk,1:2,1:nqi)
                       ! qicpx in restructed domain

      real, intent(in) :: qicpy_rst(1:ni_rst,1:nk,1:2,1:nqi)
                       ! qicpy in restructed domain

      real, intent(in) :: nicpx_rst(1:nj_rst,1:nk,1:2,1:nni)
                       ! nicpx in restructed domain

      real, intent(in) :: nicpy_rst(1:ni_rst,1:nk,1:2,1:nni)
                       ! nicpy in restructed domain

      real, intent(in) :: qcwcpx_rst(1:nj_rst,1:nk,1:2,1:nqw)
                       ! qcwcpx in restructed domain

      real, intent(in) :: qcwcpy_rst(1:ni_rst,1:nk,1:2,1:nqw)
                       ! qcwcpy in restructed domain

      real, intent(in) :: qcicpx_rst(1:nj_rst,1:nk,1:2,1:nqi)
                       ! qcicpx in restructed domain

      real, intent(in) :: qcicpy_rst(1:ni_rst,1:nk,1:2,1:nqi)
                       ! qcicpy in restructed domain

      real, intent(in) :: qacpx_rst(1:nj_rst,1:nk,1:2,1:nqa(0))
                       ! qacpx in restructed domain

      real, intent(in) :: qacpy_rst(1:ni_rst,1:nk,1:2,1:nqa(0))
                       ! qacpy in restructed domain

      real, intent(in) :: qtcpx_rst(1:nj_rst,1:nk,1:2)
                       ! qtcpx in restructed domain

      real, intent(in) :: qtcpy_rst(1:ni_rst,1:nk,1:2)
                       ! qtcpy in restructed domain

      real, intent(in) :: tkecpx_rst(1:nj_rst,1:nk,1:2)
                       ! tkecpx in restructed domain

      real, intent(in) :: tkecpy_rst(1:ni_rst,1:nk,1:2)
                       ! tkecpy in restructed domain

      real, intent(in) :: maxvl_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! maxvl in restructed domain

      real, intent(in) :: prwtr_rst(0:ni_rst+1,0:nj_rst+1,1:2,1:nqw)
                       ! prwtr in restructed domain

      real, intent(in) :: price_rst(0:ni_rst+1,0:nj_rst+1,1:2,1:nqi)
                       ! price in restructed domain

      real, intent(in) :: z0m_rst(0:ni_rst+1,0:nj_rst+1)
                       ! z0m in restructed domain

      real, intent(in) :: z0h_rst(0:ni_rst+1,0:nj_rst+1)
                       ! z0h in restructed domain

      real, intent(in) :: pdia_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! pdia in restructed domain

      real, intent(in) :: tund_rst(0:ni_rst+1,0:nj_rst+1,1:nund)
                       ! tund in restructed domain

      real, intent(in) :: tundp_rst(0:ni_rst+1,0:nj_rst+1,1:nund)
                       ! tundp in restructed domain

! Internal shared variables

      character(len=108) exprim
                       ! Optional run name

      character(len=108) crsdir
                       ! User specified directory for CReSS files

      character(len=108) dmpvar
                       ! Control flag of dumped variables

      character(len=108) rstfl
                       ! Opened file name

      integer ncexp    ! Number of character of exprim
      integer nccrs    ! Number of character of crsdir

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions
      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

      integer sfcopt   ! Option for surface physics
      integer advopt   ! Option for advection scheme
      integer cphopt   ! Option for cloud micro physics
      integer qcgopt   ! Option for charging distribution
      integer aslopt   ! Option for aerosol processes
      integer trkopt   ! Option for mixing ratio tracking
      integer tubopt   ! Option for turbulent mixing
      integer diaopt   ! Option for diabatic calculation

      integer in       ! Namelist table index

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer n        ! Array index in 4th direction

      integer ncfl     ! Number of character of restructed restart file

      integer iorst    ! Unit number of restructed restart file

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

! -----

!! Open and read in the data to the restructed restart checking file.

! Initialize the character variable.

      if(mype.eq.root) then

        call inichar(rstfl)

      end if

! -----

! Get the unit number.

      if(mype.eq.root) then

        call getunit(iorst)

      end if

! -----

! Open the restructed restart checking file.

      if(mype.eq.root) then

        rstfl(1:ncexp)=exprim(1:ncexp)

        ncfl=ncexp+26

        write(rstfl(ncexp+1:ncexp+26),'(a3,i8.8,a15)')                  &
     &                            'res',ctime/1000_i8,'.check.txt.swap'

        open(iorst,iostat=stat,err=100,                                 &
     &       file=crsdir(1:nccrs)//rstfl(1:ncfl),                       &
     &       status='new',access='sequential',form='formatted',         &
     &       blank='null',action='write')

      else

        stat=0

      end if

  100 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outrst  ',6,'cont',1,'              ',14,iorst, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outrst  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('outrst  ',6,rstfl,ncfl,iorst,1,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Reset the namelist table for parallelizing.

      if(mype.eq.root) then

        riname(idnumpe)=nsrl*nsub_rst

        riname(idxsub)=nisub_rst
        riname(idysub)=njsub_rst

      end if

! -----

! Read in the data to the restructed restart checking file.

      if(mype.eq.root) then

        write(iorst,'(a)',iostat=stat,err=110)                          &
     &       (rcname(in),in=1,ncn)

        write(iorst,*,iostat=stat,err=110)                              &
     &       (riname(in),in=1,nin)

        write(iorst,*,iostat=stat,err=110)                              &
     &       (rrname(in),in=1,nrn)

      else

        stat=0

      end if

  110 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outrst  ',6,'cont',4,'              ',14,iorst, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outrst  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('outrst  ',6,rstfl,108,iorst,4,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Close the restructed restart checking file.

      if(mype.eq.root) then

        close(iorst,iostat=stat,err=120,status='keep')

      else

        stat=0

      end if

  120 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outrst  ',6,'cont',2,'              ',14,iorst, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outrst  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('outrst  ',6,rstfl,108,iorst,2,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit number.

      if(mype.eq.root) then

        call putunit(iorst)

      end if

! -----

!! -----

!! Open and read in the data to the restructed restart file.

! Initialize the character variable.

      call inichar(dmpvar)

! -----

! Get the required namelist variables.

      call getcname(fpdmpvar,dmpvar)
      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpsbc,sbc)
      call getiname(fpnbc,nbc)
      call getiname(fpsfcopt,sfcopt)
      call getiname(fpadvopt,advopt)
      call getiname(fpcphopt,cphopt)
      call getiname(fpqcgopt,qcgopt)
      call getiname(fpaslopt,aslopt)
      call getiname(fptrkopt,trkopt)
      call getiname(fptubopt,tubopt)
      call getiname(fpdiaopt,diaopt)

! -----

! Initialize the character variable.

      call inichar(rstfl)

! -----

! Get the unit number.

      call getunit(iorst)

! -----

! Open the restructed restart file.

      rstfl(1:ncexp)=exprim(1:ncexp)

      if(ngrp.eq.1) then

        ncfl=ncexp+27

        write(rstfl(ncexp+1:ncexp+27),'(a3,i8.8,a3,i4.4,a9)')           &
     &            'res',ctime/1000_i8,'.pe',mysub_rst,'.bin.swap'

      else

        ncfl=ncexp+36

        write(rstfl(ncexp+1:ncexp+36),'(a3,i8.8,2(a4,i4.4),a9)')        &
     &     'res',ctime/1000_i8,'.grp',mygrp,'-sub',mysub_rst,'.bin.swap'

      end if

      open(iorst,iostat=stat,err=130,                                   &
     &     file=crsdir(1:nccrs)//rstfl(1:ncfl),                         &
     &     status='new',access='sequential',form='unformatted',         &
     &     position='rewind',action='write')

  130 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outrst  ',6,'cont',1,'              ',14,iorst, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outrst  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        if(fpara(1:5).eq.'multi') then

          if(ngrp.eq.1) then

            write(rstfl(ncexp+15:ncexp+18),'(a4)') 'XXXX'

          else

            write(rstfl(ncexp+16:ncexp+19),'(a4)') 'XXXX'
            write(rstfl(ncexp+24:ncexp+27),'(a4)') 'YYYY'

          end if

        end if

        call outstd03('outrst  ',6,rstfl,ncfl,iorst,1,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read in the data to the restructed restart file.

      if(advopt.le.3.and.sfcopt.ge.1) then

        write(iorst,iostat=stat,err=140)                                &
     &       (((tund_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nund)

      end if

      if(sfcopt.ge.1) then

        write(iorst,iostat=stat,err=140)                                &
     &       (((tundp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nund)

        write(iorst,iostat=stat,err=140)                                &
     &       ((z0m_rst(i,j),i=0,ni_rst+1),j=0,nj_rst+1)

        write(iorst,iostat=stat,err=140)                                &
     &       ((z0h_rst(i,j),i=0,ni_rst+1),j=0,nj_rst+1)

      end if

      write(iorst,iostat=stat,err=140) fmois(1:5)

      write(iorst,iostat=stat,err=140)                                  &
     &     (((ubr_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      write(iorst,iostat=stat,err=140)                                  &
     &     (((vbr_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      write(iorst,iostat=stat,err=140)                                  &
     &     (((pbr_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      write(iorst,iostat=stat,err=140)                                  &
     &     (((ptbr_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      write(iorst,iostat=stat,err=140)                                  &
     &     (((qvbr_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      if(advopt.le.3) then

        write(iorst,iostat=stat,err=140)                                &
     &       (((u_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

        write(iorst,iostat=stat,err=140)                                &
     &       (((v_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

        write(iorst,iostat=stat,err=140)                                &
     &       (((w_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

        write(iorst,iostat=stat,err=140)                                &
     &       (((pp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

        write(iorst,iostat=stat,err=140)                                &
     &       (((ptp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

      write(iorst,iostat=stat,err=140)                                  &
     &     (((up_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      write(iorst,iostat=stat,err=140)                                  &
     &     (((vp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      write(iorst,iostat=stat,err=140)                                  &
     &     (((wp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      write(iorst,iostat=stat,err=140)                                  &
     &     (((ppp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      write(iorst,iostat=stat,err=140)                                  &
     &     (((ptpp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

        write(iorst,iostat=stat,err=140)                                &
     &       (((ucpx_rst(j,k,i),j=1,nj_rst),k=1,nk),i=1,2)

        write(iorst,iostat=stat,err=140)                                &
     &       (((ucpy_rst(i,k,j),i=1,ni_rst),k=1,nk),j=1,2)

        write(iorst,iostat=stat,err=140)                                &
     &       (((vcpx_rst(j,k,i),j=1,nj_rst),k=1,nk),i=1,2)

        write(iorst,iostat=stat,err=140)                                &
     &       (((vcpy_rst(i,k,j),i=1,ni_rst),k=1,nk),j=1,2)

        write(iorst,iostat=stat,err=140)                                &
     &       (((wcpx_rst(j,k,i),j=1,nj_rst),k=1,nk),i=1,2)

        write(iorst,iostat=stat,err=140)                                &
     &       (((wcpy_rst(i,k,j),i=1,ni_rst),k=1,nk),j=1,2)

        write(iorst,iostat=stat,err=140)                                &
     &       (((pcpx_rst(j,k,i),j=1,nj_rst),k=1,nk),i=1,2)

        write(iorst,iostat=stat,err=140)                                &
     &       (((pcpy_rst(i,k,j),i=1,ni_rst),k=1,nk),j=1,2)

        write(iorst,iostat=stat,err=140)                                &
     &       (((ptcpx_rst(j,k,i),j=1,nj_rst),k=1,nk),i=1,2)

        write(iorst,iostat=stat,err=140)                                &
     &       (((ptcpy_rst(i,k,j),i=1,ni_rst),k=1,nk),j=1,2)

      end if

      if(fmois(1:5).eq.'moist'.and.advopt.le.3) then

        write(iorst,iostat=stat,err=140)                                &
     &       (((qv_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

      if(fmois(1:5).eq.'moist') then

        write(iorst,iostat=stat,err=140)                                &
     &       (((qvp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4)) then

        write(iorst,iostat=stat,err=140)                                &
     &       (((qvcpx_rst(j,k,i),j=1,nj_rst),k=1,nk),i=1,2)

        write(iorst,iostat=stat,err=140)                                &
     &       (((qvcpy_rst(i,k,j),i=1,ni_rst),k=1,nk),j=1,2)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (advopt.le.3.and.abs(cphopt).ge.1)) then

        write(iorst,iostat=stat,err=140)                                &
     &       ((((qwtr_rst(i,j,k,n),                                     &
     &               i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqw)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (advopt.le.3.and.(abs(cphopt).ge.4.and.abs(cphopt).lt.20))) then

        write(iorst,iostat=stat,err=140)                                &
     &       ((((nwtr_rst(i,j,k,n),                                     &
     &               i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nnw)

      end if

      if(fmois(1:5).eq.'moist'.and.abs(cphopt).ge.1) then

        write(iorst,iostat=stat,err=140)                                &
     &       ((((qwtrp_rst(i,j,k,n),                                    &
     &                i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqw)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (abs(cphopt).ge.4.and.abs(cphopt).lt.20)) then

        write(iorst,iostat=stat,err=140)                                &
     &       ((((nwtrp_rst(i,j,k,n),                                    &
     &                i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nnw)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  ((abs(cphopt).ge.1.and.abs(cphopt).lt.20).and.                  &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        write(iorst,iostat=stat,err=140)                                &
     &       ((((qwcpx_rst(j,k,i,n),j=1,nj_rst),k=1,nk),i=1,2),n=1,nqw)

        write(iorst,iostat=stat,err=140)                                &
     &       ((((qwcpy_rst(i,k,j,n),i=1,ni_rst),k=1,nk),j=1,2),n=1,nqw)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  ((abs(cphopt).ge.4.and.abs(cphopt).lt.20).and.                  &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        write(iorst,iostat=stat,err=140)                                &
     &       ((((nwcpx_rst(j,k,i,n),j=1,nj_rst),k=1,nk),i=1,2),n=1,nnw)

        write(iorst,iostat=stat,err=140)                                &
     &       ((((nwcpy_rst(i,k,j,n),i=1,ni_rst),k=1,nk),j=1,2),n=1,nnw)

      end if

      if(fmois(1:5).eq.'moist'.and.(advopt.le.3.and.                    &
     &  (abs(cphopt).ge.2.and.mod(abs(cphopt),10).ne.1))) then

        write(iorst,iostat=stat,err=140)                                &
     &       ((((qice_rst(i,j,k,n),                                     &
     &               i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqi)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (advopt.le.3.and.(abs(cphopt).ge.2.and.                         &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11))) then

        if(abs(cphopt).eq.2) then

          write(iorst,iostat=stat,err=140)                              &
     &         (((nice_rst(i,j,k,1),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

        else

          write(iorst,iostat=stat,err=140)                              &
     &         ((((nice_rst(i,j,k,n),                                   &
     &                 i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nni)

        end if

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (abs(cphopt).ge.2.and.mod(abs(cphopt),10).ne.1)) then

        write(iorst,iostat=stat,err=140)                                &
     &       ((((qicep_rst(i,j,k,n),                                    &
     &                i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqi)

      end if

      if(fmois(1:5).eq.'moist'.and.(abs(cphopt).ge.2.and.               &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11)) then

        if(abs(cphopt).eq.2) then

          write(iorst,iostat=stat,err=140)                              &
     &         (((nicep_rst(i,j,k,1),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

        else

          write(iorst,iostat=stat,err=140)                              &
     &         ((((nicep_rst(i,j,k,n),                                  &
     &                  i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nni)

        end if

      end if

      if(fmois(1:5).eq.'moist'.and.((abs(cphopt).ge.2.and.              &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11).and.                  &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        write(iorst,iostat=stat,err=140)                                &
     &       ((((qicpx_rst(j,k,i,n),j=1,nj_rst),k=1,nk),i=1,2),n=1,nqi)

        write(iorst,iostat=stat,err=140)                                &
     &       ((((qicpy_rst(i,k,j,n),i=1,ni_rst),k=1,nk),j=1,2),n=1,nqi)

        if(abs(cphopt).eq.2) then

          write(iorst,iostat=stat,err=140)                              &
     &         (((nicpx_rst(j,k,i,1),j=1,nj_rst),k=1,nk),i=1,2)

          write(iorst,iostat=stat,err=140)                              &
     &         (((nicpy_rst(i,k,j,1),i=1,ni_rst),k=1,nk),j=1,2)

        else

          write(iorst,iostat=stat,err=140)                              &
     &        ((((nicpx_rst(j,k,i,n),j=1,nj_rst),k=1,nk),i=1,2),n=1,nni)

          write(iorst,iostat=stat,err=140)                              &
     &        ((((nicpy_rst(i,k,j,n),i=1,ni_rst),k=1,nk),j=1,2),n=1,nni)

        end if

      end if

      if(fmois(1:5).eq.'moist'.and.abs(cphopt).ge.1) then

        if(abs(cphopt).lt.10) then

          write(iorst,iostat=stat,err=140)                              &
     &        (((prwtr_rst(i,j,1,n),i=0,ni_rst+1),j=0,nj_rst+1),n=1,nqw)

          write(iorst,iostat=stat,err=140)                              &
     &        (((prwtr_rst(i,j,2,n),i=0,ni_rst+1),j=0,nj_rst+1),n=1,nqw)

        else

          write(iorst,iostat=stat,err=140)                              &
     &         ((prwtr_rst(i,j,1,1),i=0,ni_rst+1),j=0,nj_rst+1)

          write(iorst,iostat=stat,err=140)                              &
     &         ((prwtr_rst(i,j,2,1),i=0,ni_rst+1),j=0,nj_rst+1)

        end if

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (abs(cphopt).ge.2.and.mod(abs(cphopt),10).ne.1)) then

        if(abs(cphopt).lt.10) then

          write(iorst,iostat=stat,err=140)                              &
     &        (((price_rst(i,j,1,n),i=0,ni_rst+1),j=0,nj_rst+1),n=1,nqi)

          write(iorst,iostat=stat,err=140)                              &
     &        (((price_rst(i,j,2,n),i=0,ni_rst+1),j=0,nj_rst+1),n=1,nqi)

        else

          write(iorst,iostat=stat,err=140)                              &
     &         ((price_rst(i,j,1,1),i=0,ni_rst+1),j=0,nj_rst+1)

          write(iorst,iostat=stat,err=140)                              &
     &         ((price_rst(i,j,2,1),i=0,ni_rst+1),j=0,nj_rst+1)

        end if

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (advopt.le.3.and.cphopt.lt.0.and.qcgopt.eq.2)) then

        write(iorst,iostat=stat,err=140)                                &
     &       ((((qcwtr_rst(i,j,k,n),                                    &
     &                i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqw)

      end if

      if(fmois(1:5).eq.'moist'.and.(advopt.le.3.and.cphopt.lt.0)) then

        write(iorst,iostat=stat,err=140)                                &
     &       ((((qcice_rst(i,j,k,n),                                    &
     &                i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqi)

      end if

      if(fmois(1:5).eq.'moist'.and.(cphopt.lt.0.and.qcgopt.eq.2)) then

        write(iorst,iostat=stat,err=140)                                &
     &       ((((qcwtrp_rst(i,j,k,n),                                   &
     &                 i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqw)

      end if

      if(fmois(1:5).eq.'moist'.and.cphopt.lt.0) then

        write(iorst,iostat=stat,err=140)                                &
     &       ((((qcicep_rst(i,j,k,n),                                   &
     &                 i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqi)

      end if

      if(fmois(1:5).eq.'moist'.and.(cphopt.lt.0.and.qcgopt.eq.2.and.    &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        write(iorst,iostat=stat,err=140)                                &
     &       ((((qcwcpx_rst(j,k,i,n),j=1,nj_rst),k=1,nk),i=1,2),n=1,nqw)

        write(iorst,iostat=stat,err=140)                                &
     &       ((((qcwcpy_rst(i,k,j,n),i=1,ni_rst),k=1,nk),j=1,2),n=1,nqw)

      end if

      if(fmois(1:5).eq.'moist'.and.(cphopt.lt.0.and.                    &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        write(iorst,iostat=stat,err=140)                                &
     &       ((((qcicpx_rst(j,k,i,n),j=1,nj_rst),k=1,nk),i=1,2),n=1,nqi)

        write(iorst,iostat=stat,err=140)                                &
     &       ((((qcicpy_rst(i,k,j,n),i=1,ni_rst),k=1,nk),j=1,2),n=1,nqi)

      end if

      if(advopt.le.3.and.aslopt.ge.1) then

        write(iorst,iostat=stat,err=140)                                &
     &       ((((qasl_rst(i,j,k,n),                                     &
     &               i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqa(0))

      end if

      if(aslopt.ge.1) then

        write(iorst,iostat=stat,err=140)                                &
     &       ((((qaslp_rst(i,j,k,n),                                    &
     &                i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqa(0))

      end if

      if(aslopt.ge.1.and.                                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4)) then

        write(iorst,iostat=stat,err=140)                                &
     &    ((((qacpx_rst(j,k,i,n),j=1,nj_rst),k=1,nk),i=1,2),n=1,nqa(0))

        write(iorst,iostat=stat,err=140)                                &
     &    ((((qacpy_rst(i,k,j,n),i=1,ni_rst),k=1,nk),j=1,2),n=1,nqa(0))

      end if

      if(advopt.le.3.and.trkopt.ge.1) then

        write(iorst,iostat=stat,err=140)                                &
     &       (((qt_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

      if(trkopt.ge.1) then

        write(iorst,iostat=stat,err=140)                                &
     &       (((qtp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

      if(trkopt.ge.1.and.                                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4)) then

        write(iorst,iostat=stat,err=140)                                &
     &       (((qtcpx_rst(j,k,i),j=1,nj_rst),k=1,nk),i=1,2)

        write(iorst,iostat=stat,err=140)                                &
     &       (((qtcpy_rst(i,k,j),i=1,ni_rst),k=1,nk),j=1,2)

      end if

      if(advopt.le.3.and.tubopt.ge.2) then

        write(iorst,iostat=stat,err=140)                                &
     &       (((tke_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

      if(tubopt.ge.2) then

        write(iorst,iostat=stat,err=140)                                &
     &       (((tkep_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

      if(tubopt.ge.2.and.                                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4)) then

        write(iorst,iostat=stat,err=140)                                &
     &       (((tkecpx_rst(j,k,i),j=1,nj_rst),k=1,nk),i=1,2)

        write(iorst,iostat=stat,err=140)                                &
     &       (((tkecpy_rst(i,k,j),i=1,ni_rst),k=1,nk),j=1,2)

      end if

      if(tubopt.ge.2.and.                                               &
     &  (dmpvar(12:12).eq.'+'.or.dmpvar(12:12).eq.'-')) then

        write(iorst,iostat=stat,err=140)                                &
     &       (((maxvl_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

      if(diaopt.eq.1) then

        write(iorst,iostat=stat,err=140)                                &
     &       (((pdia_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

  140 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outrst  ',6,'cont',4,'              ',14,iorst, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outrst  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('outrst  ',6,rstfl,108,iorst,4,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Close the restructed restart file.

      close(iorst,iostat=stat,err=150,status='keep')

  150 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outrst  ',6,'cont',2,'              ',14,iorst, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outrst  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('outrst  ',6,rstfl,108,iorst,2,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit number.

      call putunit(iorst)

! -----

!! -----

      end subroutine s_outrst

!-----7--------------------------------------------------------------7--

      end module m_outrst
