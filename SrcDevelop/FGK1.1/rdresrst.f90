!***********************************************************************
      module m_rdresrst
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/01/25, 1999/03/25, 1999/04/06, 1999/05/20,
!                   1999/06/14, 1999/06/21, 1999/07/05, 1999/08/23,
!                   1999/09/30, 1999/11/01, 1999/11/19, 2000/01/05,
!                   2000/01/17, 2000/03/23, 2000/04/18, 2000/06/01,
!                   2001/02/13, 2001/04/15, 2001/05/29, 2001/06/29,
!                   2001/07/13, 2001/08/07, 2001/09/13, 2001/10/18,
!                   2001/11/14, 2002/01/15, 2002/04/02, 2002/06/18,
!                   2002/07/03, 2002/07/15, 2002/07/23, 2002/08/15,
!                   2002/12/02, 2003/03/13, 2003/03/28, 2003/04/30,
!                   2003/05/19, 2003/06/27, 2003/11/05, 2003/11/28,
!                   2003/12/12, 2004/01/09, 2004/02/01, 2004/04/15,
!                   2004/05/31, 2004/07/01, 2004/08/01, 2004/08/20,
!                   2004/09/01, 2005/01/14, 2005/02/10, 2005/04/04,
!                   2005/08/05, 2005/10/05, 2006/01/10, 2006/02/13,
!                   2006/04/03, 2006/07/21, 2006/09/21, 2006/12/04,
!                   2007/01/05, 2007/01/20, 2007/05/07, 2007/07/30,
!                   2007/08/24, 2007/11/26, 2008/01/11, 2008/03/12,
!                   2008/05/02, 2008/07/01, 2008/08/25, 2008/10/10,
!                   2008/12/11, 2009/01/30, 2009/02/27, 2011/08/18,
!                   2011/09/22, 2011/11/10, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     open and read out the data from the restart file.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_chkfile
      use m_chkstd
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

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: rdresrst, s_rdresrst

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rdresrst

        module procedure s_rdresrst

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
      subroutine s_rdresrst(fpexprim,fpcrsdir,fpdmpvar,fpncexp,fpnccrs, &
     &                      fpwbc,fpebc,fpsbc,fpnbc,fpsfcopt,fpadvopt,  &
     &                      fpcphopt,fpqcgopt,fpaslopt,fptrkopt,        &
     &                      fptubopt,fpdiaopt,ctime,fmois,ni,nj,nk,     &
     &                      nqw,nnw,nqi,nni,nqa,nund,ubr,vbr,pbr,ptbr,  &
     &                      qvbr,u,up,v,vp,w,wp,pp,ppp,ptp,ptpp,qv,qvp, &
     &                      qwtr,qwtrp,nwtr,nwtrp,qice,qicep,nice,nicep,&
     &                      qcwtr,qcwtrp,qcice,qcicep,qasl,qaslp,qt,qtp,&
     &                      tke,tkep,ucpx,ucpy,vcpx,vcpy,wcpx,wcpy,     &
     &                      pcpx,pcpy,ptcpx,ptcpy,qvcpx,qvcpy,          &
     &                      qwcpx,qwcpy,nwcpx,nwcpy,qicpx,qicpy,        &
     &                      nicpx,nicpy,qcwcpx,qcwcpy,qcicpx,qcicpy,    &
     &                      qacpx,qacpy,qtcpx,qtcpy,tkecpx,tkecpy,maxvl,&
     &                      prwtr,price,pdia,z0m,z0h,tund,tundp)
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

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

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
                       ! x components of velocity at present

      real, intent(out) :: up(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at past

      real, intent(out) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at present

      real, intent(out) :: vp(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at past

      real, intent(out) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at present

      real, intent(out) :: wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

      real, intent(out) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at present

      real, intent(out) :: ppp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at past

      real, intent(out) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at present

      real, intent(out) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(out) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at present

      real, intent(out) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(out) :: qwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at present

      real, intent(out) :: qwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at past

      real, intent(out) :: nwtr(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at present

      real, intent(out) :: nwtrp(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at past

      real, intent(out) :: qice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at present

      real, intent(out) :: qicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at past

      real, intent(out) :: nice(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at present

      real, intent(out) :: nicep(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at past

      real, intent(out) :: qcwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at present

      real, intent(out) :: qcwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at past

      real, intent(out) :: qcice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at present

      real, intent(out) :: qcicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at past

      real, intent(out) :: qasl(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at present

      real, intent(out) :: qaslp(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at past

      real, intent(out) :: qt(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at present

      real, intent(out) :: qtp(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at past

      real, intent(out) :: tke(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at present

      real, intent(out) :: tkep(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at past

      real, intent(out) :: ucpx(1:nj,1:nk,1:2)
                       ! Phase speed of x components of velocity
                       ! on west and east boundary

      real, intent(out) :: ucpy(1:ni,1:nk,1:2)
                       ! Phase speed of x components of velocity
                       ! on south and north boundary

      real, intent(out) :: vcpx(1:nj,1:nk,1:2)
                       ! Phase speed of y components of velocity
                       ! on west and east boundary

      real, intent(out) :: vcpy(1:ni,1:nk,1:2)
                       ! Phase speed of y components of velocity
                       ! on south and north boundary

      real, intent(out) :: wcpx(1:nj,1:nk,1:2)
                       ! Phase speed of z components of velocity
                       ! on west and east boundary

      real, intent(out) :: wcpy(1:ni,1:nk,1:2)
                       ! Phase speed of z components of velocity
                       ! on south and north boundary

      real, intent(out) :: pcpx(1:nj,1:nk,1:2)
                       ! Phase speed of pressure
                       ! on west and east boundary

      real, intent(out) :: pcpy(1:ni,1:nk,1:2)
                       ! Phase speed of pressure
                       ! on south and north boundary

      real, intent(out) :: ptcpx(1:nj,1:nk,1:2)
                       ! Phase speed of potential temperature
                       ! on west and east boundary

      real, intent(out) :: ptcpy(1:ni,1:nk,1:2)
                       ! Phase speed of potential temperature
                       ! on south and north boundary

      real, intent(out) :: qvcpx(1:nj,1:nk,1:2)
                       ! Phase speed of water vapor mixing ratio
                       ! on west and east boundary

      real, intent(out) :: qvcpy(1:ni,1:nk,1:2)
                       ! Phase speed of water vapor mixing ratio
                       ! on south and north boundary

      real, intent(out) :: qwcpx(1:nj,1:nk,1:2,1:nqw)
                       ! Phase speed of water hydrometeor
                       ! on west and east boundary

      real, intent(out) :: qwcpy(1:ni,1:nk,1:2,1:nqw)
                       ! Phase speed of water hydrometeor
                       ! on south and north boundary

      real, intent(out) :: nwcpx(1:nj,1:nk,1:2,1:nnw)
                       ! Phase speed of water concentrations
                       ! on west and east boundary

      real, intent(out) :: nwcpy(1:ni,1:nk,1:2,1:nnw)
                       ! Phase speed of water concentrations
                       ! on south and north boundary

      real, intent(out) :: qicpx(1:nj,1:nk,1:2,1:nqi)
                       ! Phase speed of ice hydrometeor
                       ! on west and east boundary

      real, intent(out) :: qicpy(1:ni,1:nk,1:2,1:nqi)
                       ! Phase speed of ice hydrometeor
                       ! on south and north boundary

      real, intent(out) :: nicpx(1:nj,1:nk,1:2,1:nni)
                       ! Phase speed of ice concentrations
                       ! on west and east boundary

      real, intent(out) :: nicpy(1:ni,1:nk,1:2,1:nni)
                       ! Phase speed of ice concentrations
                       ! on south and north boundary

      real, intent(out) :: qcwcpx(1:nj,1:nk,1:2,1:nqw)
                       ! Phase speed of charging distribution for water
                       ! on west and east boundary

      real, intent(out) :: qcwcpy(1:ni,1:nk,1:2,1:nqw)
                       ! Phase speed of charging distribution for water
                       ! on south and north boundary

      real, intent(out) :: qcicpx(1:nj,1:nk,1:2,1:nqi)
                       ! Phase speed of charging distribution for ice
                       ! on west and east boundary

      real, intent(out) :: qcicpy(1:ni,1:nk,1:2,1:nqi)
                       ! Phase speed of charging distribution for ice
                       ! on south and north boundary

      real, intent(out) :: qacpx(1:nj,1:nk,1:2,1:nqa(0))
                       ! Phase speed of aerosol mixing ratio
                       ! on west and east boundary

      real, intent(out) :: qacpy(1:ni,1:nk,1:2,1:nqa(0))
                       ! Phase speed of aerosol mixing ratio
                       ! on south and north boundary

      real, intent(out) :: qtcpx(1:nj,1:nk,1:2)
                       ! Phase speed of tracer mixing ratio
                       ! on west and east boundary

      real, intent(out) :: qtcpy(1:ni,1:nk,1:2)
                       ! Phase speed of tracer mixing ratio
                       ! on south and north boundary

      real, intent(out) :: tkecpx(1:nj,1:nk,1:2)
                       ! Phase speed of turbulent kinetic energy
                       ! on west and east boundary

      real, intent(out) :: tkecpy(1:ni,1:nk,1:2)
                       ! Phase speed of turbulent kinetic energy
                       ! on south and north boundary

      real, intent(out) :: maxvl(0:ni+1,0:nj+1,1:nk)
                       ! Maximum instantaneous wind velocity

      real, intent(out) :: prwtr(0:ni+1,0:nj+1,1:2,1:nqw)
                       ! Precipitation and accumulation for water

      real, intent(out) :: price(0:ni+1,0:nj+1,1:2,1:nqi)
                       ! Precipitation and accumulation for ice

      real, intent(out) :: pdia(0:ni+1,0:nj+1,1:nk)
                       ! Diabatic value

      real, intent(out) :: z0m(0:ni+1,0:nj+1)
                       ! Roughness length for velocity

      real, intent(out) :: z0h(0:ni+1,0:nj+1)
                       ! Roughness length for scalar

      real, intent(out) :: tund(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at present

      real, intent(out) :: tundp(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at past

! Internal shared variables

      character(len=108) exprim
                       ! Optional run name

      character(len=108) crsdir
                       ! User specified directory for CReSS files

      character(len=108) dmpvar
                       ! Control flag of dumped variables

      character(len=108) resfl
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

      integer ncfl     ! Number of character of restart file

      integer iores    ! Unit number of restart file

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

!! Open and read out the data from the restart checking file.

! Initialize the character variable.

      if(mype.eq.root) then

        call inichar(resfl)

      end if

! -----

! Get the unit number.

      if(mype.eq.root) then

        call getunit(iores)

      end if

! -----

! Open the restart checking file.

      if(mype.eq.root) then

        resfl(1:ncexp)=exprim(1:ncexp)

        write(resfl(ncexp+1:ncexp+21),'(a3,i8.8,a10)')                  &
     &                            'res',ctime/1000_i8,'.check.txt'

        open(iores,iostat=stat,err=100,                                 &
     &       file=crsdir(1:nccrs)//resfl(1:ncexp+21),                   &
     &       status='old',access='sequential',form='formatted',         &
     &       blank='null',position='rewind',action='read')

      else

        stat=0

      end if

  100 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdresrst',8,'cont',1,'              ',14,iores, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdresrst',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('rdresrst',8,resfl,ncexp+21,iores,1,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read out the data from the restart checking file.

      if(mype.eq.root) then

        read(iores,'(a)',iostat=stat,end=110,err=110)                   &
     &      (rcname(in),in=1,ncn)

        read(iores,*,iostat=stat,end=110,err=110)                       &
     &      (riname(in),in=1,nin)

        read(iores,*,iostat=stat,end=110,err=110)                       &
     &      (rrname(in),in=1,nrn)

      else

        stat=0

      end if

  110 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdresrst',8,'cont',3,'              ',14,iores, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdresrst',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('rdresrst',8,resfl,108,iores,3,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Check the restart checking file.

      if(mype.eq.root) then

        call chkfile('rst',stat,ncn,nin,nrn,                            &
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

! Close the restart checking file.

      if(mype.eq.root) then

        close(iores,iostat=stat,err=120,status='keep')

      else

        stat=0

      end if

  120 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdresrst',8,'cont',2,'              ',14,iores, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdresrst',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('rdresrst',8,resfl,108,iores,2,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit number.

      if(mype.eq.root) then

        call putunit(iores)

      end if

! -----

!! -----

!! Open and read out the data from the restart file.

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

      call inichar(resfl)

! -----

! Get the unit number.

      call getunit(iores)

! -----

! Open the restart file.

      resfl(1:ncexp)=exprim(1:ncexp)

      if(ngrp.eq.1) then

        ncfl=ncexp+22

        write(resfl(ncexp+1:ncexp+22),'(a3,i8.8,a3,i4.4,a4)')           &
     &            'res',ctime/1000_i8,'.pe',mysub,'.bin'

      else

        ncfl=ncexp+31

        write(resfl(ncexp+1:ncexp+31),'(a3,i8.8,2(a4,i4.4),a4)')        &
     &            'res',ctime/1000_i8,'.grp',mygrp,'-sub',mysub,'.bin'

      end if

      open(iores,iostat=stat,err=130,                                   &
     &     file=crsdir(1:nccrs)//resfl(1:ncfl),                         &
     &     status='old',access='sequential',form='unformatted',         &
     &     position='rewind',action='read')

  130 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdresrst',8,'cont',1,'              ',14,iores, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdresrst',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        if(fpara(1:5).eq.'multi') then

          if(ngrp.eq.1) then

            write(resfl(ncexp+15:ncexp+18),'(a4)') 'XXXX'

          else

            write(resfl(ncexp+16:ncexp+19),'(a4)') 'XXXX'
            write(resfl(ncexp+24:ncexp+27),'(a4)') 'YYYY'

          end if

        end if

        call outstd03('rdresrst',8,resfl,ncfl,iores,1,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read out the data from the restart file.

      if(advopt.le.3.and.sfcopt.ge.1) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((tund(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nund)

      end if

      if(sfcopt.ge.1) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((tundp(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nund)

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((z0m(i,j),i=0,ni+1),j=0,nj+1)

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((z0h(i,j),i=0,ni+1),j=0,nj+1)

      end if

      read(iores,iostat=stat,end=140,err=140) fmois(1:5)

      read(iores,iostat=stat,end=140,err=140)                           &
     &    (((ubr(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      read(iores,iostat=stat,end=140,err=140)                           &
     &    (((vbr(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      read(iores,iostat=stat,end=140,err=140)                           &
     &    (((pbr(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      read(iores,iostat=stat,end=140,err=140)                           &
     &    (((ptbr(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      read(iores,iostat=stat,end=140,err=140)                           &
     &    (((qvbr(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      if(advopt.le.3) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((u(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((v(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((w(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((pp(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((ptp(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      end if

      read(iores,iostat=stat,end=140,err=140)                           &
     &    (((up(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      read(iores,iostat=stat,end=140,err=140)                           &
     &    (((vp(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      read(iores,iostat=stat,end=140,err=140)                           &
     &    (((wp(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      read(iores,iostat=stat,end=140,err=140)                           &
     &    (((ppp(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      read(iores,iostat=stat,end=140,err=140)                           &
     &    (((ptpp(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((ucpx(j,k,i),j=1,nj),k=1,nk),i=1,2)

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((ucpy(i,k,j),i=1,ni),k=1,nk),j=1,2)

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((vcpx(j,k,i),j=1,nj),k=1,nk),i=1,2)

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((vcpy(i,k,j),i=1,ni),k=1,nk),j=1,2)

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((wcpx(j,k,i),j=1,nj),k=1,nk),i=1,2)

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((wcpy(i,k,j),i=1,ni),k=1,nk),j=1,2)

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((pcpx(j,k,i),j=1,nj),k=1,nk),i=1,2)

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((pcpy(i,k,j),i=1,ni),k=1,nk),j=1,2)

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((ptcpx(j,k,i),j=1,nj),k=1,nk),i=1,2)

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((ptcpy(i,k,j),i=1,ni),k=1,nk),j=1,2)

      end if

      if(fmois(1:5).eq.'moist'.and.advopt.le.3) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((qv(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      end if

      if(fmois(1:5).eq.'moist') then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((qvp(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4)) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((qvcpx(j,k,i),j=1,nj),k=1,nk),i=1,2)

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((qvcpy(i,k,j),i=1,ni),k=1,nk),j=1,2)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (advopt.le.3.and.abs(cphopt).ge.1)) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((((qwtr(i,j,k,n),i=0,ni+1),j=0,nj+1),k=1,nk),n=1,nqw)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (advopt.le.3.and.(abs(cphopt).ge.4.and.abs(cphopt).lt.20))) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((((nwtr(i,j,k,n),i=0,ni+1),j=0,nj+1),k=1,nk),n=1,nnw)

      end if

      if(fmois(1:5).eq.'moist'.and.abs(cphopt).ge.1) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((((qwtrp(i,j,k,n),i=0,ni+1),j=0,nj+1),k=1,nk),n=1,nqw)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (abs(cphopt).ge.4.and.abs(cphopt).lt.20)) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((((nwtrp(i,j,k,n),i=0,ni+1),j=0,nj+1),k=1,nk),n=1,nnw)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  ((abs(cphopt).ge.1.and.abs(cphopt).lt.20).and.                  &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((((qwcpx(j,k,i,n),j=1,nj),k=1,nk),i=1,2),n=1,nqw)

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((((qwcpy(i,k,j,n),i=1,ni),k=1,nk),j=1,2),n=1,nqw)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  ((abs(cphopt).ge.4.and.abs(cphopt).lt.20).and.                  &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((((nwcpx(j,k,i,n),j=1,nj),k=1,nk),i=1,2),n=1,nnw)

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((((nwcpy(i,k,j,n),i=1,ni),k=1,nk),j=1,2),n=1,nnw)

      end if

      if(fmois(1:5).eq.'moist'.and.(advopt.le.3.and.                    &
     &  (abs(cphopt).ge.2.and.mod(abs(cphopt),10).ne.1))) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((((qice(i,j,k,n),i=0,ni+1),j=0,nj+1),k=1,nk),n=1,nqi)

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (advopt.le.3.and.(abs(cphopt).ge.2.and.                         &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11))) then

        if(abs(cphopt).eq.2) then

          read(iores,iostat=stat,end=140,err=140)                       &
     &        (((nice(i,j,k,1),i=0,ni+1),j=0,nj+1),k=1,nk)

        else

          read(iores,iostat=stat,end=140,err=140)                       &
     &        ((((nice(i,j,k,n),i=0,ni+1),j=0,nj+1),k=1,nk),n=1,nni)

        end if

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (abs(cphopt).ge.2.and.mod(abs(cphopt),10).ne.1)) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((((qicep(i,j,k,n),i=0,ni+1),j=0,nj+1),k=1,nk),n=1,nqi)

      end if

      if(fmois(1:5).eq.'moist'.and.(abs(cphopt).ge.2.and.               &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11)) then

        if(abs(cphopt).eq.2) then

          read(iores,iostat=stat,end=140,err=140)                       &
     &        (((nicep(i,j,k,1),i=0,ni+1),j=0,nj+1),k=1,nk)

        else

          read(iores,iostat=stat,end=140,err=140)                       &
     &        ((((nicep(i,j,k,n),i=0,ni+1),j=0,nj+1),k=1,nk),n=1,nni)

        end if

      end if

      if(fmois(1:5).eq.'moist'.and.((abs(cphopt).ge.2.and.              &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11).and.                  &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((((qicpx(j,k,i,n),j=1,nj),k=1,nk),i=1,2),n=1,nqi)

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((((qicpy(i,k,j,n),i=1,ni),k=1,nk),j=1,2),n=1,nqi)

        if(abs(cphopt).eq.2) then

          read(iores,iostat=stat,end=140,err=140)                       &
     &        (((nicpx(j,k,i,1),j=1,nj),k=1,nk),i=1,2)

          read(iores,iostat=stat,end=140,err=140)                       &
     &        (((nicpy(i,k,j,1),i=1,ni),k=1,nk),j=1,2)

        else

          read(iores,iostat=stat,end=140,err=140)                       &
     &        ((((nicpx(j,k,i,n),j=1,nj),k=1,nk),i=1,2),n=1,nni)

          read(iores,iostat=stat,end=140,err=140)                       &
     &        ((((nicpy(i,k,j,n),i=1,ni),k=1,nk),j=1,2),n=1,nni)

        end if

      end if

      if(fmois(1:5).eq.'moist'.and.abs(cphopt).ge.1) then

        if(abs(cphopt).lt.10) then

          read(iores,iostat=stat,end=140,err=140)                       &
     &        (((prwtr(i,j,1,n),i=0,ni+1),j=0,nj+1),n=1,nqw)

          read(iores,iostat=stat,end=140,err=140)                       &
     &        (((prwtr(i,j,2,n),i=0,ni+1),j=0,nj+1),n=1,nqw)

        else

          read(iores,iostat=stat,end=140,err=140)                       &
     &        ((prwtr(i,j,1,1),i=0,ni+1),j=0,nj+1)

          read(iores,iostat=stat,end=140,err=140)                       &
     &        ((prwtr(i,j,2,1),i=0,ni+1),j=0,nj+1)

        end if

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (abs(cphopt).ge.2.and.mod(abs(cphopt),10).ne.1)) then

        if(abs(cphopt).lt.10) then

          read(iores,iostat=stat,end=140,err=140)                       &
     &        (((price(i,j,1,n),i=0,ni+1),j=0,nj+1),n=1,nqi)

          read(iores,iostat=stat,end=140,err=140)                       &
     &        (((price(i,j,2,n),i=0,ni+1),j=0,nj+1),n=1,nqi)

        else

          read(iores,iostat=stat,end=140,err=140)                       &
     &        ((price(i,j,1,1),i=0,ni+1),j=0,nj+1)

          read(iores,iostat=stat,end=140,err=140)                       &
     &        ((price(i,j,2,1),i=0,ni+1),j=0,nj+1)

        end if

      end if

      if(fmois(1:5).eq.'moist'.and.                                     &
     &  (advopt.le.3.and.cphopt.lt.0.and.qcgopt.eq.2)) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((((qcwtr(i,j,k,n),i=0,ni+1),j=0,nj+1),k=1,nk),n=1,nqw)

      end if

      if(fmois(1:5).eq.'moist'.and.(advopt.le.3.and.cphopt.lt.0)) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((((qcice(i,j,k,n),i=0,ni+1),j=0,nj+1),k=1,nk),n=1,nqi)

      end if

      if(fmois(1:5).eq.'moist'.and.(cphopt.lt.0.and.qcgopt.eq.2)) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((((qcwtrp(i,j,k,n),i=0,ni+1),j=0,nj+1),k=1,nk),n=1,nqw)

      end if

      if(fmois(1:5).eq.'moist'.and.cphopt.lt.0) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((((qcicep(i,j,k,n),i=0,ni+1),j=0,nj+1),k=1,nk),n=1,nqi)

      end if

      if(fmois(1:5).eq.'moist'.and.(cphopt.lt.0.and.qcgopt.eq.2.and.    &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((((qcwcpx(j,k,i,n),j=1,nj),k=1,nk),i=1,2),n=1,nqw)

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((((qcwcpy(i,k,j,n),i=1,ni),k=1,nk),j=1,2),n=1,nqw)

      end if

      if(fmois(1:5).eq.'moist'.and.(cphopt.lt.0.and.                    &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((((qcicpx(j,k,i,n),j=1,nj),k=1,nk),i=1,2),n=1,nqi)

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((((qcicpy(i,k,j,n),i=1,ni),k=1,nk),j=1,2),n=1,nqi)

      end if

      if(advopt.le.3.and.aslopt.ge.1) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((((qasl(i,j,k,n),i=0,ni+1),j=0,nj+1),k=1,nk),n=1,nqa(0))

      end if

      if(aslopt.ge.1) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((((qaslp(i,j,k,n),i=0,ni+1),j=0,nj+1),k=1,nk),n=1,nqa(0))

      end if

      if(aslopt.ge.1.and.                                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4)) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((((qacpx(j,k,i,n),j=1,nj),k=1,nk),i=1,2),n=1,nqa(0))

        read(iores,iostat=stat,end=140,err=140)                         &
     &      ((((qacpy(i,k,j,n),i=1,ni),k=1,nk),j=1,2),n=1,nqa(0))

      end if

      if(advopt.le.3.and.trkopt.ge.1) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((qt(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      end if

      if(trkopt.ge.1) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((qtp(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      end if

      if(trkopt.ge.1.and.                                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4)) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((qtcpx(j,k,i),j=1,nj),k=1,nk),i=1,2)

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((qtcpy(i,k,j),i=1,ni),k=1,nk),j=1,2)

      end if

      if(advopt.le.3.and.tubopt.ge.2) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((tke(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      end if

      if(tubopt.ge.2) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((tkep(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      end if

      if(tubopt.ge.2.and.                                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4)) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((tkecpx(j,k,i),j=1,nj),k=1,nk),i=1,2)

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((tkecpy(i,k,j),i=1,ni),k=1,nk),j=1,2)

      end if

      if(tubopt.ge.2.and.                                               &
     &  (dmpvar(12:12).eq.'+'.or.dmpvar(12:12).eq.'-')) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((maxvl(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      end if

      if(diaopt.eq.1) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((pdia(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      end if

  140 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdresrst',8,'cont',3,'              ',14,iores, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdresrst',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('rdresrst',8,resfl,108,iores,3,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Close the restart file.

      close(iores,iostat=stat,err=150,status='keep')

  150 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdresrst',8,'cont',2,'              ',14,iores, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdresrst',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('rdresrst',8,resfl,108,iores,2,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit number.

      call putunit(iores)

! -----

!! -----

      end subroutine s_rdresrst

!-----7--------------------------------------------------------------7--

      end module m_rdresrst
