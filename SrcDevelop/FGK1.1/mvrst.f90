!***********************************************************************
      module m_mvrst
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2005/08/05
!     Modification: 2005/10/05, 2006/01/10, 2006/02/13, 2006/04/03,
!                   2006/07/21, 2006/09/21, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/05/07, 2007/05/14, 2007/07/30,
!                   2007/08/24, 2007/11/26, 2008/03/12, 2008/05/02,
!                   2008/07/01, 2008/08/25, 2008/10/10, 2008/12/11,
!                   2009/01/30, 2009/02/27, 2011/08/18, 2011/09/22,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     rename restructed restart files.

!-----7--------------------------------------------------------------7--

! Module reference

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

      public :: mvrst, s_mvrst

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface mvrst

        module procedure s_mvrst

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
      subroutine s_mvrst(fpexprim,fpcrsdir,fpdmpvar,fpncexp,fpnccrs,    &
     &               fpwbc,fpebc,fpsbc,fpnbc,fpsfcopt,fpadvopt,fpcphopt,&
     &               fpqcgopt,fpaslopt,fptrkopt,fptubopt,fpdiaopt,      &
     &               ctime,fmois,nk,nqw,nnw,nqi,nni,nqa,nund,           &
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

! Internal shared variables

      character(len=108) exprim
                       ! Optional run name

      character(len=108) crsdir
                       ! User specified directory for CReSS files

      character(len=108) dmpvar
                       ! Control flag of dumped variables

      character(len=108) resfl
                       ! Opened file name

      character(len=108) rstfl
                       ! Opened file name

      character(len=5), intent(inout) :: fmois
                       ! Control flag of air moisture

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

      integer ncres    ! Number of character of restart file name
      integer ncrst    ! Number of character of restructed file name

      integer iores    ! Unit number of restart file
      integer iorst    ! Unit number of restructed restart file

      integer stat     ! Runtime status

      integer broot    ! Broadcasting root

      real, intent(inout) :: ubr_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! ubr in restructed domain

      real, intent(inout) :: vbr_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! vbr in restructed domain

      real, intent(inout) :: pbr_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! pbr in restructed domain

      real, intent(inout) :: ptbr_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! ptbr in restructed domain

      real, intent(inout) :: qvbr_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! qvbr in restructed domain

      real, intent(inout) :: u_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! u in restructed domain

      real, intent(inout) :: up_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! up in restructed domain

      real, intent(inout) :: v_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! v in restructed domain

      real, intent(inout) :: vp_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! vp in restructed domain

      real, intent(inout) :: w_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! w in restructed domain

      real, intent(inout) :: wp_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! wp in restructed domain

      real, intent(inout) :: pp_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! pp in restructed domain

      real, intent(inout) :: ppp_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! ppp in restructed domain

      real, intent(inout) :: ptp_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! ptp in restructed domain

      real, intent(inout) :: ptpp_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! ptpp in restructed domain

      real, intent(inout) :: qv_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! qv in restructed domain

      real, intent(inout) :: qvp_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! qvp in restructed domain

      real, intent(inout) :: qwtr_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqw)
                       ! qwtr in restructed domain

      real, intent(inout) :: qwtrp_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqw)
                       ! qwtrp in restructed domain

      real, intent(inout) :: nwtr_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nnw)
                       ! nwtr in restructed domain

      real, intent(inout) :: nwtrp_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nnw)
                       ! nwtrp in restructed domain

      real, intent(inout) :: qice_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqi)
                       ! qice in restructed domain

      real, intent(inout) :: qicep_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqi)
                       ! qicep in restructed domain

      real, intent(inout) :: nice_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nni)
                       ! nice in restructed domain

      real, intent(inout) :: nicep_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nni)
                       ! nicep in restructed domain

      real, intent(inout) ::                                            &
     &                   qcwtr_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqw)
                       ! qcwtr in restructed domain

      real, intent(inout) ::                                            &
     &                   qcwtrp_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqw)
                       ! qcwtrp in restructed domain

      real, intent(inout) ::                                            &
     &                   qcice_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqi)
                       ! qcice in restructed domain

      real, intent(inout) ::                                            &
     &                   qcicep_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqi)
                       ! qcicep in restructed domain

      real, intent(inout) ::                                            &
     &                   qasl_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqa(0))
                       ! qasl in restructed domain

      real, intent(inout) ::                                            &
     &                   qaslp_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqa(0))
                       ! qaslp in restructed domain

      real, intent(inout) :: qt_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! qt in restructed domain

      real, intent(inout) :: qtp_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! qtp in restructed domain

      real, intent(inout) :: tke_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! tke in restructed domain

      real, intent(inout) :: tkep_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! tkep in restructed domain

      real, intent(inout) :: ucpx_rst(1:nj_rst,1:nk,1:2)
                       ! ucpx in restructed domain

      real, intent(inout) :: ucpy_rst(1:ni_rst,1:nk,1:2)
                       ! ucpy in restructed domain

      real, intent(inout) :: vcpx_rst(1:nj_rst,1:nk,1:2)
                       ! vcpx in restructed domain

      real, intent(inout) :: vcpy_rst(1:ni_rst,1:nk,1:2)
                       ! vcpy in restructed domain

      real, intent(inout) :: wcpx_rst(1:nj_rst,1:nk,1:2)
                       ! wcpx in restructed domain

      real, intent(inout) :: wcpy_rst(1:ni_rst,1:nk,1:2)
                       ! wcpy in restructed domain

      real, intent(inout) :: pcpx_rst(1:nj_rst,1:nk,1:2)
                       ! pcpx in restructed domain

      real, intent(inout) :: pcpy_rst(1:ni_rst,1:nk,1:2)
                       ! pcpy in restructed domain

      real, intent(inout) :: ptcpx_rst(1:nj_rst,1:nk,1:2)
                       ! ptcpx in restructed domain

      real, intent(inout) :: ptcpy_rst(1:ni_rst,1:nk,1:2)
                       ! ptcpy in restructed domain

      real, intent(inout) :: qvcpx_rst(1:nj_rst,1:nk,1:2)
                       ! qvcpx in restructed domain

      real, intent(inout) :: qvcpy_rst(1:ni_rst,1:nk,1:2)
                       ! qvcpy in restructed domain

      real, intent(inout) :: qwcpx_rst(1:nj_rst,1:nk,1:2,1:nqw)
                       ! qwcpx in restructed domain

      real, intent(inout) :: qwcpy_rst(1:ni_rst,1:nk,1:2,1:nqw)
                       ! qwcpy in restructed domain

      real, intent(inout) :: nwcpx_rst(1:nj_rst,1:nk,1:2,1:nnw)
                       ! nwcpx in restructed domain

      real, intent(inout) :: nwcpy_rst(1:ni_rst,1:nk,1:2,1:nnw)
                       ! nwcpy in restructed domain

      real, intent(inout) :: qicpx_rst(1:nj_rst,1:nk,1:2,1:nqi)
                       ! qicpx in restructed domain

      real, intent(inout) :: qicpy_rst(1:ni_rst,1:nk,1:2,1:nqi)
                       ! qicpy in restructed domain

      real, intent(inout) :: nicpx_rst(1:nj_rst,1:nk,1:2,1:nni)
                       ! nicpx in restructed domain

      real, intent(inout) :: nicpy_rst(1:ni_rst,1:nk,1:2,1:nni)
                       ! nicpy in restructed domain

      real, intent(inout) :: qcwcpx_rst(1:nj_rst,1:nk,1:2,1:nqw)
                       ! qcwcpx in restructed domain

      real, intent(inout) :: qcwcpy_rst(1:ni_rst,1:nk,1:2,1:nqw)
                       ! qcwcpy in restructed domain

      real, intent(inout) :: qcicpx_rst(1:nj_rst,1:nk,1:2,1:nqi)
                       ! qcicpx in restructed domain

      real, intent(inout) :: qcicpy_rst(1:ni_rst,1:nk,1:2,1:nqi)
                       ! qcicpy in restructed domain

      real, intent(inout) :: qacpx_rst(1:nj_rst,1:nk,1:2,1:nqa(0))
                       ! qacpx in restructed domain

      real, intent(inout) :: qacpy_rst(1:ni_rst,1:nk,1:2,1:nqa(0))
                       ! qacpy in restructed domain

      real, intent(inout) :: qtcpx_rst(1:nj_rst,1:nk,1:2)
                       ! qtcpx in restructed domain

      real, intent(inout) :: qtcpy_rst(1:ni_rst,1:nk,1:2)
                       ! qtcpy in restructed domain

      real, intent(inout) :: tkecpx_rst(1:nj_rst,1:nk,1:2)
                       ! tkecpx in restructed domain

      real, intent(inout) :: tkecpy_rst(1:ni_rst,1:nk,1:2)
                       ! tkecpy in restructed domain

      real, intent(inout) :: maxvl_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! maxvl in restructed domain

      real, intent(inout) :: prwtr_rst(0:ni_rst+1,0:nj_rst+1,1:2,1:nqw)
                       ! prwtr in restructed domain

      real, intent(inout) :: price_rst(0:ni_rst+1,0:nj_rst+1,1:2,1:nqi)
                       ! price in restructed domain

      real, intent(inout) :: pdia_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! pdia in restructed domain

      real, intent(inout) :: z0m_rst(0:ni_rst+1,0:nj_rst+1)
                       ! z0m in restructed domain

      real, intent(inout) :: z0h_rst(0:ni_rst+1,0:nj_rst+1)
                       ! z0h in restructed domain

      real, intent(inout) :: tund_rst(0:ni_rst+1,0:nj_rst+1,1:nund)
                       ! tund in restructed domain

      real, intent(inout) :: tundp_rst(0:ni_rst+1,0:nj_rst+1,1:nund)
                       ! tundp in restructed domain

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

!! Rename the restructed restart checking file.

! Initialize the character variables.

      if(mype.eq.root) then

        call inichar(rstfl)
        call inichar(resfl)

      end if

! -----

! Get the unit numbers.

      if(mype.eq.root) then

        call getunit(iorst)
        call getunit(iores)

      end if

! -----

! Open the restructed restart checking file.

      if(mype.eq.root) then

        rstfl(1:ncexp)=exprim(1:ncexp)

        ncrst=ncexp+26

        write(rstfl(ncexp+1:ncexp+26),'(a3,i8.8,a15)')                  &
     &                            'res',ctime/1000_i8,'.check.txt.swap'

        open(iorst,iostat=stat,err=100,                                 &
     &       file=crsdir(1:nccrs)//rstfl(1:ncrst),                      &
     &       status='old',access='sequential',form='formatted',         &
     &       blank='null',action='readwrite')

      else

        stat=0

      end if

  100 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('mvrst   ',5,'cont',1,'              ',14,iorst, &
     &                 stat)

        end if

        call cpondpe

        call destroy('mvrst   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('mvrst   ',5,rstfl,ncrst,iorst,1,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Open the new restart checking file.

      if(mype.eq.root) then

        resfl(1:ncexp)=exprim(1:ncexp)

        ncres=ncexp+21

        write(resfl(ncexp+1:ncexp+21),'(a3,i8.8,a10)')                  &
     &                                'res',ctime/1000_i8,'.check.txt'

        open(iores,iostat=stat,err=110,                                 &
     &       file=crsdir(1:nccrs)//resfl(1:ncres),                      &
     &       status='new',access='sequential',form='formatted',         &
     &       blank='null',action='write')

      else

        stat=0

      end if

  110 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('mvrst   ',5,'cont',1,'              ',14,iores, &
     &                 stat)

        end if

        call cpondpe

        call destroy('mvrst   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('mvrst   ',5,resfl,ncres,iores,1,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read out the data from the restructed restart checking file.

      if(mype.eq.root) then

        read(iorst,'(a)',iostat=stat,end=120,err=120)                   &
     &      (rcname(in),in=1,ncn)

        read(iorst,*,iostat=stat,end=120,err=120)                       &
     &      (riname(in),in=1,nin)

        read(iorst,*,iostat=stat,end=120,err=120)                       &
     &      (rrname(in),in=1,nrn)

      else

        stat=0

      end if

  120 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('mvrst   ',5,'cont',3,'              ',14,iorst, &
     &                 stat)

        end if

        call cpondpe

        call destroy('mvrst   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('mvrst   ',5,rstfl,108,iorst,3,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read in the data to the new restart checking file.

      if(mype.eq.root) then

        write(iores,'(a)',iostat=stat,err=130)                          &
     &       (rcname(in),in=1,ncn)

        write(iores,*,iostat=stat,err=130)                              &
     &       (riname(in),in=1,nin)

        write(iores,*,iostat=stat,err=130)                              &
     &       (rrname(in),in=1,nrn)

      else

        stat=0

      end if

  130 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('mvrst   ',5,'cont',4,'              ',14,iores, &
     &                 stat)

        end if

        call cpondpe

        call destroy('mvrst   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('mvrst   ',5,resfl,108,iores,4,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Close the restructed restart checking file.

      if(mype.eq.root) then

        close(iorst,iostat=stat,err=140,status='delete')

      else

        stat=0

      end if

  140 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('mvrst   ',5,'cont',2,'              ',14,iorst, &
     &                 stat)

        end if

        call cpondpe

        call destroy('mvrst   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('mvrst   ',5,rstfl,108,iorst,2,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Close the new restart checking file.

      if(mype.eq.root) then

        close(iores,iostat=stat,err=150,status='keep')

      else

        stat=0

      end if

  150 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('mvrst   ',5,'cont',2,'              ',14,iores, &
     &                 stat)

        end if

        call cpondpe

        call destroy('mvrst   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('mvrst   ',5,resfl,108,iores,2,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit numbers.

      if(mype.eq.root) then

        call putunit(iores)
        call putunit(iorst)

      end if

! -----

!! -----

!! Rename restructed restart files.

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

! Initialize the character variables.

      call inichar(rstfl)
      call inichar(resfl)

! -----

! Get the unit numbers.

      call getunit(iorst)
      call getunit(iores)

! -----

! Open the restructed restart file.

      rstfl(1:ncexp)=exprim(1:ncexp)

      if(ngrp.eq.1) then

        ncrst=ncexp+27

        write(rstfl(ncexp+1:ncexp+27),'(a3,i8.8,a3,i4.4,a9)')           &
     &            'res',ctime/1000_i8,'.pe',mysub_rst,'.bin.swap'

      else

        ncrst=ncexp+36

        write(rstfl(ncexp+1:ncexp+36),'(a3,i8.8,2(a4,i4.4),a9)')        &
     &     'res',ctime/1000_i8,'.grp',mygrp,'-sub',mysub_rst,'.bin.swap'

      end if

      open(iorst,iostat=stat,err=160,                                   &
     &     file=crsdir(1:nccrs)//rstfl(1:ncrst),                        &
     &     status='old',access='sequential',form='unformatted',         &
     &     position='rewind',action='readwrite')

  160 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('mvrst   ',5,'cont',1,'              ',14,iorst, &
     &                 stat)

        end if

        call cpondpe

        call destroy('mvrst   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('mvrst   ',5,rstfl,ncrst,iorst,1,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Open the new restart file.

      resfl(1:ncexp)=exprim(1:ncexp)

      if(ngrp.eq.1) then

        ncres=ncexp+22

        write(resfl(ncexp+1:ncexp+22),'(a3,i8.8,a3,i4.4,a4)')           &
     &            'res',ctime/1000_i8,'.pe',mysub_rst,'.bin'

      else

        ncres=ncexp+31

        write(resfl(ncexp+1:ncexp+31),'(a3,i8.8,2(a4,i4.4),a4)')        &
     &       'res',ctime/1000_i8,'.grp',mygrp,'-sub',mysub_rst,'.bin'

      end if

      open(iores,iostat=stat,err=170,                                   &
     &     file=crsdir(1:nccrs)//resfl(1:ncres),                        &
     &     status='new',access='sequential',form='unformatted',         &
     &     position='rewind',action='write')

  170 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('mvrst   ',5,'cont',1,'              ',14,iores, &
     &                 stat)

        end if

        call cpondpe

        call destroy('mvrst   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('mvrst   ',5,resfl,ncres,iores,1,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read out the data from the restructed restart file.

      if(advopt.le.3.and.sfcopt.ge.1) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((tund_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nund)

      end if

      if(sfcopt.ge.1) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((tundp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nund)

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      ((z0m_rst(i,j),i=0,ni_rst+1),j=0,nj_rst+1)

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      ((z0h_rst(i,j),i=0,ni_rst+1),j=0,nj_rst+1)

      end if

      read(iorst,iostat=stat,end=180,err=180) fmois(1:5)

      read(iorst,iostat=stat,end=180,err=180)                           &
     &    (((ubr_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      read(iorst,iostat=stat,end=180,err=180)                           &
     &    (((vbr_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      read(iorst,iostat=stat,end=180,err=180)                           &
     &    (((pbr_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      read(iorst,iostat=stat,end=180,err=180)                           &
     &    (((ptbr_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      read(iorst,iostat=stat,end=180,err=180)                           &
     &    (((qvbr_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      if(advopt.le.3) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((u_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((v_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((w_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((pp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((ptp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

      read(iorst,iostat=stat,end=180,err=180)                           &
     &    (((up_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      read(iorst,iostat=stat,end=180,err=180)                           &
     &    (((vp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      read(iorst,iostat=stat,end=180,err=180)                           &
     &    (((wp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      read(iorst,iostat=stat,end=180,err=180)                           &
     &    (((ppp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      read(iorst,iostat=stat,end=180,err=180)                           &
     &    (((ptpp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((ucpx_rst(j,k,i),j=1,nj_rst),k=1,nk),i=1,2)

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((ucpy_rst(i,k,j),i=1,ni_rst),k=1,nk),j=1,2)

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((vcpx_rst(j,k,i),j=1,nj_rst),k=1,nk),i=1,2)

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((vcpy_rst(i,k,j),i=1,ni_rst),k=1,nk),j=1,2)

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((wcpx_rst(j,k,i),j=1,nj_rst),k=1,nk),i=1,2)

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((wcpy_rst(i,k,j),i=1,ni_rst),k=1,nk),j=1,2)

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((pcpx_rst(j,k,i),j=1,nj_rst),k=1,nk),i=1,2)

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((pcpy_rst(i,k,j),i=1,ni_rst),k=1,nk),j=1,2)

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((ptcpx_rst(j,k,i),j=1,nj_rst),k=1,nk),i=1,2)

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((ptcpy_rst(i,k,j),i=1,ni_rst),k=1,nk),j=1,2)

      end if

      if(fmois(1:5).eq.'moist'.and.advopt.le.3) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((qv_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

      if(fmois(1:5).eq.'moist') then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((qvp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4)) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((qvcpx_rst(j,k,i),j=1,nj_rst),k=1,nk),i=1,2)

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((qvcpy_rst(i,k,j),i=1,ni_rst),k=1,nk),j=1,2)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (advopt.le.3.and.abs(cphopt).ge.1)) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      ((((qwtr_rst(i,j,k,n),                                      &
     &              i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqw)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (advopt.le.3.and.(abs(cphopt).ge.4.and.abs(cphopt).lt.20))) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      ((((nwtr_rst(i,j,k,n),                                      &
     &              i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nnw)

      end if

      if(fmois(1:5).eq.'moist'.and.abs(cphopt).ge.1) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      ((((qwtrp_rst(i,j,k,n),                                     &
     &               i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqw)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (abs(cphopt).ge.4.and.abs(cphopt).lt.20)) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      ((((nwtrp_rst(i,j,k,n),                                     &
     &               i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nnw)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  ((abs(cphopt).ge.1.and.abs(cphopt).lt.20).and.                  &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      ((((qwcpx_rst(j,k,i,n),j=1,nj_rst),k=1,nk),i=1,2),n=1,nqw)

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      ((((qwcpy_rst(i,k,j,n),i=1,ni_rst),k=1,nk),j=1,2),n=1,nqw)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  ((abs(cphopt).ge.4.and.abs(cphopt).lt.20).and.                  &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      ((((nwcpx_rst(j,k,i,n),j=1,nj_rst),k=1,nk),i=1,2),n=1,nnw)

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      ((((nwcpy_rst(i,k,j,n),i=1,ni_rst),k=1,nk),j=1,2),n=1,nnw)

      end if

      if(fmois(1:5).eq.'moist'.and.(advopt.le.3.and.                    &
     &  (abs(cphopt).ge.2.and.mod(abs(cphopt),10).ne.1))) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      ((((qice_rst(i,j,k,n),                                      &
     &              i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqi)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (advopt.le.3.and.(abs(cphopt).ge.2.and.                         &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11))) then

        if(abs(cphopt).eq.2) then

          read(iorst,iostat=stat,end=180,err=180)                       &
     &        (((nice_rst(i,j,k,1),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

        else

          read(iorst,iostat=stat,end=180,err=180)                       &
     &        ((((nice_rst(i,j,k,n),                                    &
     &                i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nni)

        end if

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (abs(cphopt).ge.2.and.mod(abs(cphopt),10).ne.1)) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      ((((qicep_rst(i,j,k,n),                                     &
     &               i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqi)

      end if

      if(fmois(1:5).eq.'moist'.and.(abs(cphopt).ge.2.and.               &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11)) then

        if(abs(cphopt).eq.2) then

          read(iorst,iostat=stat,end=180,err=180)                       &
     &        (((nicep_rst(i,j,k,1),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

        else

          read(iorst,iostat=stat,end=180,err=180)                       &
     &        ((((nicep_rst(i,j,k,n),                                   &
     &                 i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nni)

        end if

      end if

      if(fmois(1:5).eq.'moist'.and.((abs(cphopt).ge.2.and.              &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11).and.                  &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      ((((qicpx_rst(j,k,i,n),j=1,nj_rst),k=1,nk),i=1,2),n=1,nqi)

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      ((((qicpy_rst(i,k,j,n),i=1,ni_rst),k=1,nk),j=1,2),n=1,nqi)

        if(abs(cphopt).eq.2) then

          read(iorst,iostat=stat,end=180,err=180)                       &
     &        (((nicpx_rst(j,k,i,1),j=1,nj_rst),k=1,nk),i=1,2)

          read(iorst,iostat=stat,end=180,err=180)                       &
     &        (((nicpy_rst(i,k,j,1),i=1,ni_rst),k=1,nk),j=1,2)

        else

          read(iorst,iostat=stat,end=180,err=180)                       &
     &        ((((nicpx_rst(j,k,i,n),j=1,nj_rst),k=1,nk),i=1,2),n=1,nni)

          read(iorst,iostat=stat,end=180,err=180)                       &
     &        ((((nicpy_rst(i,k,j,n),i=1,ni_rst),k=1,nk),j=1,2),n=1,nni)

        end if

      end if

      if(fmois(1:5).eq.'moist'.and.abs(cphopt).ge.1) then

        if(abs(cphopt).lt.10) then

          read(iorst,iostat=stat,end=180,err=180)                       &
     &        (((prwtr_rst(i,j,1,n),i=0,ni_rst+1),j=0,nj_rst+1),n=1,nqw)

          read(iorst,iostat=stat,end=180,err=180)                       &
     &        (((prwtr_rst(i,j,2,n),i=0,ni_rst+1),j=0,nj_rst+1),n=1,nqw)

        else

          read(iorst,iostat=stat,end=180,err=180)                       &
     &        ((prwtr_rst(i,j,1,1),i=0,ni_rst+1),j=0,nj_rst+1)

          read(iorst,iostat=stat,end=180,err=180)                       &
     &        ((prwtr_rst(i,j,2,1),i=0,ni_rst+1),j=0,nj_rst+1)

        end if

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (abs(cphopt).ge.2.and.mod(abs(cphopt),10).ne.1)) then

        if(abs(cphopt).lt.10) then

          read(iorst,iostat=stat,end=180,err=180)                       &
     &        (((price_rst(i,j,1,n),i=0,ni_rst+1),j=0,nj_rst+1),n=1,nqi)

          read(iorst,iostat=stat,end=180,err=180)                       &
     &        (((price_rst(i,j,2,n),i=0,ni_rst+1),j=0,nj_rst+1),n=1,nqi)

        else

          read(iorst,iostat=stat,end=180,err=180)                       &
     &        ((price_rst(i,j,1,1),i=0,ni_rst+1),j=0,nj_rst+1)

          read(iorst,iostat=stat,end=180,err=180)                       &
     &        ((price_rst(i,j,2,1),i=0,ni_rst+1),j=0,nj_rst+1)

        end if

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (advopt.le.3.and.cphopt.lt.0.and.qcgopt.eq.2)) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      ((((qcwtr_rst(i,j,k,n),                                     &
     &               i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqw)

      end if

      if(fmois(1:5).eq.'moist'.and.(advopt.le.3.and.cphopt.lt.0)) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      ((((qcice_rst(i,j,k,n),                                     &
     &               i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqi)

      end if

      if(fmois(1:5).eq.'moist'.and.(cphopt.lt.0.and.qcgopt.eq.2)) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      ((((qcwtrp_rst(i,j,k,n),                                    &
     &                i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqw)

      end if

      if(fmois(1:5).eq.'moist'.and.cphopt.lt.0) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      ((((qcicep_rst(i,j,k,n),                                    &
     &                i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqi)

      end if

      if(fmois(1:5).eq.'moist'.and.(cphopt.lt.0.and.qcgopt.eq.2.and.    &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      ((((qcwcpx_rst(j,k,i,n),j=1,nj_rst),k=1,nk),i=1,2),n=1,nqw)

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      ((((qcwcpy_rst(i,k,j,n),i=1,ni_rst),k=1,nk),j=1,2),n=1,nqw)

      end if

      if(fmois(1:5).eq.'moist'.and.(cphopt.lt.0.and.                    &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      ((((qcicpx_rst(j,k,i,n),j=1,nj_rst),k=1,nk),i=1,2),n=1,nqi)

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      ((((qcicpy_rst(i,k,j,n),i=1,ni_rst),k=1,nk),j=1,2),n=1,nqi)

      end if

      if(advopt.le.3.and.aslopt.ge.1) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      ((((qasl_rst(i,j,k,n),                                      &
     &              i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqa(0))

      end if

      if(aslopt.ge.1) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      ((((qaslp_rst(i,j,k,n),                                     &
     &               i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqa(0))

      end if

      if(aslopt.ge.1.and.                                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4)) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &     ((((qacpx_rst(j,k,i,n),j=1,nj_rst),k=1,nk),i=1,2),n=1,nqa(0))

        read(iorst,iostat=stat,end=180,err=180)                         &
     &     ((((qacpy_rst(i,k,j,n),i=1,ni_rst),k=1,nk),j=1,2),n=1,nqa(0))

      end if

      if(advopt.le.3.and.trkopt.ge.1) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((qt_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

      if(trkopt.ge.1) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((qtp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

      if(trkopt.ge.1.and.                                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4)) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((qtcpx_rst(j,k,i),j=1,nj_rst),k=1,nk),i=1,2)

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((qtcpy_rst(i,k,j),i=1,ni_rst),k=1,nk),j=1,2)

      end if

      if(advopt.le.3.and.tubopt.ge.2) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((tke_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

      if(tubopt.ge.2) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((tkep_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

      if(tubopt.ge.2.and.                                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4)) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((tkecpx_rst(j,k,i),j=1,nj_rst),k=1,nk),i=1,2)

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((tkecpy_rst(i,k,j),i=1,ni_rst),k=1,nk),j=1,2)

      end if

      if(tubopt.ge.2.and.                                               &
     &  (dmpvar(12:12).eq.'+'.or.dmpvar(12:12).eq.'-')) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((maxvl_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

      if(diaopt.eq.1) then

        read(iorst,iostat=stat,end=180,err=180)                         &
     &      (((pdia_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

  180 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('mvrst   ',5,'cont',3,'              ',14,iorst, &
     &                 stat)

        end if

        call cpondpe

        call destroy('mvrst   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('mvrst   ',5,rstfl,108,iorst,3,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read in the data to the new restart file.

      if(advopt.le.3.and.sfcopt.ge.1) then

        write(iores,iostat=stat,err=190)                                &
     &       (((tund_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nund)

      end if

      if(sfcopt.ge.1) then

        write(iores,iostat=stat,err=190)                                &
     &       (((tundp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nund)

        write(iores,iostat=stat,err=190)                                &
     &       ((z0m_rst(i,j),i=0,ni_rst+1),j=0,nj_rst+1)

        write(iores,iostat=stat,err=190)                                &
     &       ((z0h_rst(i,j),i=0,ni_rst+1),j=0,nj_rst+1)

      end if

      write(iores,iostat=stat,err=190) fmois(1:5)

      write(iores,iostat=stat,err=190)                                  &
     &     (((ubr_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      write(iores,iostat=stat,err=190)                                  &
     &     (((vbr_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      write(iores,iostat=stat,err=190)                                  &
     &     (((pbr_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      write(iores,iostat=stat,err=190)                                  &
     &     (((ptbr_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      write(iores,iostat=stat,err=190)                                  &
     &     (((qvbr_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      if(advopt.le.3) then

        write(iores,iostat=stat,err=190)                                &
     &       (((u_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

        write(iores,iostat=stat,err=190)                                &
     &       (((v_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

        write(iores,iostat=stat,err=190)                                &
     &       (((w_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

        write(iores,iostat=stat,err=190)                                &
     &       (((pp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

        write(iores,iostat=stat,err=190)                                &
     &       (((ptp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

      write(iores,iostat=stat,err=190)                                  &
     &     (((up_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      write(iores,iostat=stat,err=190)                                  &
     &     (((vp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      write(iores,iostat=stat,err=190)                                  &
     &     (((wp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      write(iores,iostat=stat,err=190)                                  &
     &     (((ppp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      write(iores,iostat=stat,err=190)                                  &
     &     (((ptpp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

        write(iores,iostat=stat,err=190)                                &
     &       (((ucpx_rst(j,k,i),j=1,nj_rst),k=1,nk),i=1,2)

        write(iores,iostat=stat,err=190)                                &
     &       (((ucpy_rst(i,k,j),i=1,ni_rst),k=1,nk),j=1,2)

        write(iores,iostat=stat,err=190)                                &
     &       (((vcpx_rst(j,k,i),j=1,nj_rst),k=1,nk),i=1,2)

        write(iores,iostat=stat,err=190)                                &
     &       (((vcpy_rst(i,k,j),i=1,ni_rst),k=1,nk),j=1,2)

        write(iores,iostat=stat,err=190)                                &
     &       (((wcpx_rst(j,k,i),j=1,nj_rst),k=1,nk),i=1,2)

        write(iores,iostat=stat,err=190)                                &
     &       (((wcpy_rst(i,k,j),i=1,ni_rst),k=1,nk),j=1,2)

        write(iores,iostat=stat,err=190)                                &
     &       (((pcpx_rst(j,k,i),j=1,nj_rst),k=1,nk),i=1,2)

        write(iores,iostat=stat,err=190)                                &
     &       (((pcpy_rst(i,k,j),i=1,ni_rst),k=1,nk),j=1,2)

        write(iores,iostat=stat,err=190)                                &
     &       (((ptcpx_rst(j,k,i),j=1,nj_rst),k=1,nk),i=1,2)

        write(iores,iostat=stat,err=190)                                &
     &       (((ptcpy_rst(i,k,j),i=1,ni_rst),k=1,nk),j=1,2)

      end if

      if(fmois(1:5).eq.'moist'.and.advopt.le.3) then

        write(iores,iostat=stat,err=190)                                &
     &       (((qv_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

      if(fmois(1:5).eq.'moist') then

        write(iores,iostat=stat,err=190)                                &
     &       (((qvp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4)) then

        write(iores,iostat=stat,err=190)                                &
     &       (((qvcpx_rst(j,k,i),j=1,nj_rst),k=1,nk),i=1,2)

        write(iores,iostat=stat,err=190)                                &
     &       (((qvcpy_rst(i,k,j),i=1,ni_rst),k=1,nk),j=1,2)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (advopt.le.3.and.abs(cphopt).ge.1)) then

        write(iores,iostat=stat,err=190)                                &
     &       ((((qwtr_rst(i,j,k,n),                                     &
     &               i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqw)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (advopt.le.3.and.(abs(cphopt).ge.4.and.abs(cphopt).lt.20))) then

        write(iores,iostat=stat,err=190)                                &
     &       ((((nwtr_rst(i,j,k,n),                                     &
     &               i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nnw)

      end if

      if(fmois(1:5).eq.'moist'.and.abs(cphopt).ge.1) then

        write(iores,iostat=stat,err=190)                                &
     &       ((((qwtrp_rst(i,j,k,n),                                    &
     &                i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqw)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (abs(cphopt).ge.4.and.abs(cphopt).lt.20)) then

        write(iores,iostat=stat,err=190)                                &
     &       ((((nwtrp_rst(i,j,k,n),                                    &
     &                i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nnw)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  ((abs(cphopt).ge.1.and.abs(cphopt).lt.20).and.                  &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        write(iores,iostat=stat,err=190)                                &
     &       ((((qwcpx_rst(j,k,i,n),j=1,nj_rst),k=1,nk),i=1,2),n=1,nqw)

        write(iores,iostat=stat,err=190)                                &
     &       ((((qwcpy_rst(i,k,j,n),i=1,ni_rst),k=1,nk),j=1,2),n=1,nqw)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  ((abs(cphopt).ge.4.and.abs(cphopt).lt.20).and.                  &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        write(iores,iostat=stat,err=190)                                &
     &       ((((nwcpx_rst(j,k,i,n),j=1,nj_rst),k=1,nk),i=1,2),n=1,nnw)

        write(iores,iostat=stat,err=190)                                &
     &       ((((nwcpy_rst(i,k,j,n),i=1,ni_rst),k=1,nk),j=1,2),n=1,nnw)

      end if

      if(fmois(1:5).eq.'moist'.and.(advopt.le.3.and.                    &
     &  (abs(cphopt).ge.2.and.mod(abs(cphopt),10).ne.1))) then

        write(iores,iostat=stat,err=190)                                &
     &       ((((qice_rst(i,j,k,n),                                     &
     &               i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqi)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (advopt.le.3.and.(abs(cphopt).ge.2.and.                         &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11))) then

        if(abs(cphopt).eq.2) then

          write(iores,iostat=stat,err=190)                              &
     &         (((nice_rst(i,j,k,1),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

        else

          write(iores,iostat=stat,err=190)                              &
     &         ((((nice_rst(i,j,k,n),                                   &
     &                 i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nni)

        end if

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (abs(cphopt).ge.2.and.mod(abs(cphopt),10).ne.1)) then

        write(iores,iostat=stat,err=190)                                &
     &       ((((qicep_rst(i,j,k,n),                                    &
     &                i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqi)

      end if

      if(fmois(1:5).eq.'moist'.and.(abs(cphopt).ge.2.and.               &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11)) then

        if(abs(cphopt).eq.2) then

          write(iores,iostat=stat,err=190)                              &
     &         (((nicep_rst(i,j,k,1),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

        else

          write(iores,iostat=stat,err=190)                              &
     &         ((((nicep_rst(i,j,k,n),                                  &
     &                  i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nni)

        end if

      end if

      if(fmois(1:5).eq.'moist'.and.((abs(cphopt).ge.2.and.              &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11).and.                  &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        write(iores,iostat=stat,err=190)                                &
     &       ((((qicpx_rst(j,k,i,n),j=1,nj_rst),k=1,nk),i=1,2),n=1,nqi)

        write(iores,iostat=stat,err=190)                                &
     &       ((((qicpy_rst(i,k,j,n),i=1,ni_rst),k=1,nk),j=1,2),n=1,nqi)

        if(abs(cphopt).eq.2) then

          write(iores,iostat=stat,err=190)                              &
     &        (((nicpx_rst(j,k,i,1),j=1,nj_rst),k=1,nk),i=1,2)

          write(iores,iostat=stat,err=190)                              &
     &        (((nicpy_rst(i,k,j,1),i=1,ni_rst),k=1,nk),j=1,2)

        else

          write(iores,iostat=stat,err=190)                              &
     &        ((((nicpx_rst(j,k,i,n),j=1,nj_rst),k=1,nk),i=1,2),n=1,nni)

          write(iores,iostat=stat,err=190)                              &
     &        ((((nicpy_rst(i,k,j,n),i=1,ni_rst),k=1,nk),j=1,2),n=1,nni)

        end if

      end if

      if(fmois(1:5).eq.'moist'.and.abs(cphopt).ge.1) then

        if(abs(cphopt).lt.10) then

          write(iores,iostat=stat,err=190)                              &
     &        (((prwtr_rst(i,j,1,n),i=0,ni_rst+1),j=0,nj_rst+1),n=1,nqw)

          write(iores,iostat=stat,err=190)                              &
     &        (((prwtr_rst(i,j,2,n),i=0,ni_rst+1),j=0,nj_rst+1),n=1,nqw)

        else

          write(iores,iostat=stat,err=190)                              &
     &         ((prwtr_rst(i,j,1,1),i=0,ni_rst+1),j=0,nj_rst+1)

          write(iores,iostat=stat,err=190)                              &
     &         ((prwtr_rst(i,j,2,1),i=0,ni_rst+1),j=0,nj_rst+1)

        end if

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (abs(cphopt).ge.2.and.mod(abs(cphopt),10).ne.1)) then

        if(abs(cphopt).lt.10) then

          write(iores,iostat=stat,err=190)                              &
     &        (((price_rst(i,j,1,n),i=0,ni_rst+1),j=0,nj_rst+1),n=1,nqi)

          write(iores,iostat=stat,err=190)                              &
     &        (((price_rst(i,j,2,n),i=0,ni_rst+1),j=0,nj_rst+1),n=1,nqi)

        else

          write(iores,iostat=stat,err=190)                              &
     &         ((price_rst(i,j,1,1),i=0,ni_rst+1),j=0,nj_rst+1)

          write(iores,iostat=stat,err=190)                              &
     &         ((price_rst(i,j,2,1),i=0,ni_rst+1),j=0,nj_rst+1)

        end if

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (advopt.le.3.and.cphopt.lt.0.and.qcgopt.eq.2)) then

        write(iores,iostat=stat,err=190)                                &
     &       ((((qcwtr_rst(i,j,k,n),                                    &
     &                i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqw)

      end if

      if(fmois(1:5).eq.'moist'.and.(advopt.le.3.and.cphopt.lt.0)) then

        write(iores,iostat=stat,err=190)                                &
     &       ((((qcice_rst(i,j,k,n),                                    &
     &                i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqi)

      end if

      if(fmois(1:5).eq.'moist'.and.(cphopt.lt.0.and.qcgopt.eq.2)) then

        write(iores,iostat=stat,err=190)                                &
     &       ((((qcwtrp_rst(i,j,k,n),                                   &
     &                 i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqw)

      end if

      if(fmois(1:5).eq.'moist'.and.cphopt.lt.0) then

        write(iores,iostat=stat,err=190)                                &
     &       ((((qcicep_rst(i,j,k,n),                                   &
     &                 i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqi)

      end if

      if(fmois(1:5).eq.'moist'.and.(cphopt.lt.0.and.qcgopt.eq.2.and.    &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        write(iores,iostat=stat,err=190)                                &
     &       ((((qcwcpx_rst(j,k,i,n),j=1,nj_rst),k=1,nk),i=1,2),n=1,nqw)

        write(iores,iostat=stat,err=190)                                &
     &       ((((qcwcpy_rst(i,k,j,n),i=1,ni_rst),k=1,nk),j=1,2),n=1,nqw)

      end if

      if(fmois(1:5).eq.'moist'.and.(cphopt.lt.0.and.                    &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        write(iores,iostat=stat,err=190)                                &
     &       ((((qcicpx_rst(j,k,i,n),j=1,nj_rst),k=1,nk),i=1,2),n=1,nqi)

        write(iores,iostat=stat,err=190)                                &
     &       ((((qcicpy_rst(i,k,j,n),i=1,ni_rst),k=1,nk),j=1,2),n=1,nqi)

      end if

      if(advopt.le.3.and.aslopt.ge.1) then

        write(iores,iostat=stat,err=190)                                &
     &       ((((qasl_rst(i,j,k,n),                                     &
     &               i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqa(0))

      end if

      if(aslopt.ge.1) then

        write(iores,iostat=stat,err=190)                                &
     &       ((((qaslp_rst(i,j,k,n),                                    &
     &                i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk),n=1,nqa(0))

      end if

      if(aslopt.ge.1.and.                                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4)) then

        write(iores,iostat=stat,err=190)                                &
     &    ((((qacpx_rst(j,k,i,n),j=1,nj_rst),k=1,nk),i=1,2),n=1,nqa(0))

        write(iores,iostat=stat,err=190)                                &
     &    ((((qacpy_rst(i,k,j,n),i=1,ni_rst),k=1,nk),j=1,2),n=1,nqa(0))

      end if

      if(advopt.le.3.and.trkopt.ge.1) then

        write(iores,iostat=stat,err=190)                                &
     &       (((qt_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

      if(trkopt.ge.1) then

        write(iores,iostat=stat,err=190)                                &
     &       (((qtp_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

      if(trkopt.ge.1.and.                                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4)) then

        write(iores,iostat=stat,err=190)                                &
     &       (((qtcpx_rst(j,k,i),j=1,nj_rst),k=1,nk),i=1,2)

        write(iores,iostat=stat,err=190)                                &
     &       (((qtcpy_rst(i,k,j),i=1,ni_rst),k=1,nk),j=1,2)

      end if

      if(advopt.le.3.and.tubopt.ge.2) then

        write(iores,iostat=stat,err=190)                                &
     &       (((tke_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

      if(tubopt.ge.2) then

        write(iores,iostat=stat,err=190)                                &
     &       (((tkep_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

      if(tubopt.ge.2.and.                                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4)) then

        write(iores,iostat=stat,err=190)                                &
     &       (((tkecpx_rst(j,k,i),j=1,nj_rst),k=1,nk),i=1,2)

        write(iores,iostat=stat,err=190)                                &
     &       (((tkecpy_rst(i,k,j),i=1,ni_rst),k=1,nk),j=1,2)

      end if

      if(tubopt.ge.2.and.                                               &
     &  (dmpvar(12:12).eq.'+'.or.dmpvar(12:12).eq.'-')) then

        write(iores,iostat=stat,err=190)                                &
     &       (((maxvl_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

      if(diaopt.eq.1) then

        write(iores,iostat=stat,err=190)                                &
     &       (((pdia_rst(i,j,k),i=0,ni_rst+1),j=0,nj_rst+1),k=1,nk)

      end if

  190 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('mvrst   ',5,'cont',4,'              ',14,iores, &
     &                 stat)

        end if

        call cpondpe

        call destroy('mvrst   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('mvrst   ',5,resfl,108,iores,4,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Close the restructed restart file.

      close(iorst,iostat=stat,err=200,status='delete')

  200 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('mvrst   ',5,'cont',2,'              ',14,iorst, &
     &                 stat)

        end if

        call cpondpe

        call destroy('mvrst   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('mvrst   ',5,rstfl,108,iorst,2,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Close the new restart file.

      close(iores,iostat=stat,err=210,status='keep')

  210 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('mvrst   ',5,'cont',2,'              ',14,iores, &
     &                 stat)

        end if

        call cpondpe

        call destroy('mvrst   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('mvrst   ',5,resfl,108,iores,2,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit numbers.

      call putunit(iorst)
      call putunit(iores)

! -----

!! -----

      end subroutine s_mvrst

!-----7--------------------------------------------------------------7--

      end module m_mvrst
