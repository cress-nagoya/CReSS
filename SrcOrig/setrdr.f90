!***********************************************************************
      module m_setrdr
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/09/09
!     Modification: 2002/10/31, 2002/12/02, 2003/03/28, 2003/04/30,
!                   2003/05/19, 2003/09/01, 2003/11/05, 2003/12/12,
!                   2004/03/05, 2004/04/15, 2004/05/07, 2004/05/31,
!                   2004/06/10, 2004/08/20, 2004/09/01, 2005/02/10,
!                   2005/04/04, 2005/08/05, 2006/01/10, 2006/08/18,
!                   2006/09/21, 2006/09/30, 2006/11/06, 2006/12/04,
!                   2007/01/20, 2007/07/30, 2008/04/17, 2008/05/02,
!                   2008/07/25, 2008/08/25, 2008/10/10, 2009/02/27,
!                   2009/11/13, 2011/08/18, 2011/09/22, 2013/02/13,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the inferior procedures to interpolate the radar data.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bcrdr
      use m_castvar
      use m_chkerr
      use m_chkrdr
      use m_chkstd
      use m_comindx
      use m_comkind
      use m_commpi
      use m_comsave
      use m_dbz2kg
      use m_destroy
      use m_distrqp
      use m_getcname
      use m_getdate
      use m_getiname
      use m_getrname
      use m_getxy
      use m_getz
      use m_inichar
      use m_intrprdr
      use m_rdradar
      use m_rotuvm2s
      use m_rotuvs2m
      use m_setproj
      use m_vintrdr
      use m_xy2ll

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: setrdr, s_setrdr

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface setrdr

        module procedure s_setrdr

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic int
      intrinsic max
      intrinsic mod
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_setrdr(fprdrvar,fpdatype_rdr,fpidate,fpcphopt,       &
     &                 fprotopt_rdr,fprdritv,fpngrstr,fpngrend,fpngraff,&
     &                 ctime,ftime,frdr,rtinc,ni,nj,nk,nqw,nqi,zph,     &
     &                 rbr,qwtr,qice,urdr,vrdr,wrdr,qwrdr,qirdr,        &
     &                 z,qprdr,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,           &
     &                 nid_rdr,njd_rdr,nkd_rdr,km_rdr,lon_rdr,          &
     &                 z_rdr,u_rdr,v_rdr,w_rdr,qp_rdr,tmp1_rdr)
!***********************************************************************

! Input variables

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

      integer(kind=i8), intent(in) :: ftime
                       ! Model forecast time at 1 step future

      integer, intent(in) :: fprdrvar
                       ! Formal parameter of unique index of rdrvar

      integer, intent(in) :: fpdatype_rdr
                       ! Formal parameter of unique index of datype_rdr

      integer, intent(in) :: fpidate
                       ! Formal parameter of unique index of idate

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fprotopt_rdr
                       ! Formal parameter of unique index of rotopt_rdr

      integer, intent(in) :: fprdritv
                       ! Formal parameter of unique index of rdritv

      integer, intent(in) :: fpngrstr
                       ! Formal parameter of unique index of ngrstr

      integer, intent(in) :: fpngrend
                       ! Formal parameter of unique index of ngrend

      integer, intent(in) :: fpngraff
                       ! Formal parameter of unique index of ngraff

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

      integer, intent(in) :: nid_rdr
                       ! Radar data dimension in x direction

      integer, intent(in) :: njd_rdr
                       ! Radar data dimension in y direction

      integer, intent(in) :: nkd_rdr
                       ! Radar data dimension in z direction

      integer, intent(in) :: km_rdr
                       ! Dimension of max(nk, nkd_rdr)

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: qwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor

      real, intent(in) :: qice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor

! Input and output variable

      integer, intent(inout) :: frdr(1:2)
                       ! Descriptor to put into motion
                       ! for radar nudging

! Output variables

      real, intent(out) :: rtinc(1:2)
                       ! Lapse of forecast time from radar data reading

      real, intent(out) :: urdr(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity of radar data
                       ! at marked time

      real, intent(out) :: vrdr(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity of radar data
                       ! at marked time

      real, intent(out) :: wrdr(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity of radar data
                       ! at marked time

      real, intent(out) :: qwrdr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor of radar data
                       ! at marked time

      real, intent(out) :: qirdr(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor of radar data
                       ! at marked time

! Internal shared variables

      character(len=108) rdrvar
                       ! Control flag of input radar data variables

      character(len=108) datype_rdr
                       ! Control flag of radar data type

      character(len=108) idate
                       ! Forecast start date
                       ! with Gregorian calendar, yyyymmddhhmm

      character(len=12) nxtdat
                       ! Forecast date to read out
                       ! with Gregorian calendar, yyyymmddhhmm

      integer cphopt   ! Option for cloud micro physics

      integer rotopt_rdr
                       ! Option for rotation of wind direction
                       ! in radar data

      integer(kind=i8) str01
                       ! int(ngrstr + 0.1)

      integer(kind=i8) nxtmin
                       ! 60000 x (nxtrdr(2) / 60)

      integer(kind=i8) nxtsec
                       ! mod(nxtrdr(2), 60)

      integer(kind=i8) crtime
                       ! Current forecast time to read out

      integer(kind=i8) rdr103
                       ! 1000 x int(rdritv + 0.1)

      integer(kind=i8) str103
                       ! 1000 x int(ngrstr + 0.1)

      integer(kind=i8) end103
                       ! 1000 x int(ngrend + 0.1)

      integer(kind=i8) aff103
                       ! 1000 x int(ngraff + 0.1)

      integer stat     ! Runtime status in this routine

      integer broot    ! Broadcasting root

      integer nid_rdr_sub
                       ! Substitute for nid_rdr

      integer njd_rdr_sub
                       ! Substitute for njd_rdr

      integer nkd_rdr_sub
                       ! Substitute for nkd_rdr

      real rdritv      ! Time interval of radar data file

      real ngrstr      ! Analysis nudging start time to radar data
      real ngrend      ! Analysis nudging end time to radar data
      real ngraff      ! Analysis nudging affected time to radar data

      real x0          ! x origin of model grid
      real y0          ! y origin of model grid

      real cpj(1:7)    ! Map projection parameters of model grid

      real x0rdr       ! x origin of radar data grid
      real y0rdr       ! y origin of radar data grid

      real cpjrdr(1:7) ! Map projection parameters of radar data grid

      real, intent(inout) :: z(1:nk)
                       ! zeta coordinates

      real, intent(inout) :: qprdr(0:ni+1,0:nj+1,1:nk)
                       ! Precipitation mixing ratio of radar data
                       ! at marked time

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp4(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp5(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp6(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: lon_rdr(1:nid_rdr,1:njd_rdr)
                       ! Longitude in radar data

      real, intent(inout) :: z_rdr(1:nid_rdr,1:njd_rdr,1:nkd_rdr)
                       ! z physical coordinates in radar data

      real, intent(inout) :: u_rdr(1:nid_rdr,1:njd_rdr,1:km_rdr)
                       ! x components of velocity in radar data

      real, intent(inout) :: v_rdr(1:nid_rdr,1:njd_rdr,1:km_rdr)
                       ! y components of velocity in radar data

      real, intent(inout) :: w_rdr(1:nid_rdr,1:njd_rdr,1:km_rdr)
                       ! z components of velocity in radar data

      real, intent(inout) :: qp_rdr(1:nid_rdr,1:njd_rdr,1:km_rdr)
                       ! Precipitation mixing ratio in radar data

      real, intent(inout) :: tmp1_rdr(1:nid_rdr,1:njd_rdr,1:nk)
                       ! Temporary array

!-----7--------------------------------------------------------------7--

! Initialize character variables.

      call inichar(rdrvar)
      call inichar(datype_rdr)
      call inichar(idate)

! -----

! Get the required namelist variables.

      call getcname(fprdrvar,rdrvar)
      call getcname(fpdatype_rdr,datype_rdr)
      call getcname(fpidate,idate)
      call getiname(fpcphopt,cphopt)
      call getiname(fprotopt_rdr,rotopt_rdr)
      call getrname(fprdritv,rdritv)
      call getrname(fpngrstr,ngrstr)
      call getrname(fpngrend,ngrend)
      call getrname(fpngraff,ngraff)

! -----

! Set the common used variables.

      str01=int(ngrstr+.1e0)

      rdr103=1000_i8*int(rdritv+.1e0,i8)
      str103=1000_i8*int(ngrstr+.1e0,i8)
      end103=1000_i8*int(ngrend+.1e0,i8)
      aff103=1000_i8*int(ngraff+.1e0,i8)

      crtime=rdr103*((ctime+aff103-str103)/rdr103)+str103

      if(crtime.gt.end103) then
        crtime=crtime-rdr103
      end if

      rtinc(2)=abs(.001e0*real(ctime-crtime))

      if(frdr(2).eq.0) then

        if(ctime-aff103.gt.crtime) then
          nxtrdr(2)=max((crtime+rdr103)/1000_i8,str01)
        else
          nxtrdr(2)=max(crtime/1000_i8,str01)
        end if

      end if

! -----

! Set the substituted variables.

      nid_rdr_sub=nid_rdr
      njd_rdr_sub=njd_rdr
      nkd_rdr_sub=nkd_rdr

! -----

!!! Interpolate the radar data.

      if(frdr(2).ge.0.and.                                              &
     &  (ctime+aff103)/1000_i8.ge.nxtrdr(2).and.ctime.lt.end103) then

! Calculate the forecast date to read out.

        nxtmin=60000_i8*(nxtrdr(2)/60_i8)
        nxtsec=mod(nxtrdr(2),60)

        call getdate(idate,nxtmin,nxtdat)

! -----

! Read out the data from the radar data file.

        stat=0

        if(mype.eq.root) then

          call rdradar(iddatdir,idrdrvar,idncdat,idwlngth,              &
     &                 nxtdat,ftime,nxtsec,stat,nid_rdr,njd_rdr,nkd_rdr,&
     &                 z_rdr,u_rdr,v_rdr,w_rdr,qp_rdr)

        end if

        call chkerr(stat)

        if(stat.lt.0) then

          call destroy('rdradar ',7,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        call chkstd(root)

! -----

! Check the radar data variables.

        stat=0

        if(mype.eq.root) then

          call chkrdr(idrdrvar,iddatype_rdr,ftime,stat,                 &
     &           nid_rdr,njd_rdr,nkd_rdr,z_rdr,u_rdr,v_rdr,w_rdr,qp_rdr)

        end if

        call chkerr(stat)

        if(stat.lt.0) then

          call destroy('chkrdr  ',6,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        call chkstd(root)

! -----

! Broadcast the radar data variables.

        broot=stat-1

        call chkstd(broot)

        if(stat.gt.0) then

          stat=0

          call castvar(1,nid_rdr,1,njd_rdr,1,nkd_rdr,z_rdr)

          if(rdrvar(1:1).eq.'o') then

            call castvar(1,nid_rdr,1,njd_rdr,1,nkd_rdr,u_rdr)

          end if

          if(rdrvar(2:2).eq.'o') then

            call castvar(1,nid_rdr,1,njd_rdr,1,nkd_rdr,v_rdr)

          end if

          if(rdrvar(3:3).eq.'o') then

            call castvar(1,nid_rdr,1,njd_rdr,1,nkd_rdr,w_rdr)

          end if

          if(rdrvar(4:4).eq.'o') then

            call castvar(1,nid_rdr,1,njd_rdr,1,nkd_rdr,qp_rdr)

          end if

        end if

! -----

!! Perform interpolating in the case the runtime status is equal to 0.

        if(stat.eq.0) then

! Get the zeta coordinates at the model grid points.

          call getz(iddz,idzsfc,nk,z)

! -----

! Set the map projection parameters of the model grid.

          call setproj(idmpopt,idnspol,iddx,iddy,idulat,idulon,idriu,   &
     &                 idrju,idtlat1,idtlat2,idtlon,'solver  ',6,       &
     &                 x0,y0,cpj)

! -----

! Set the map projection parameters of the data grid.

          call setproj(idmpopt_rdr,idnspol_rdr,iddx_rdr,iddy_rdr,       &
     &                 idulat_rdr,idulon_rdr,idriu_rdr,idrju_rdr,       &
     &                 idtlat1_rdr,idtlat2_rdr,idtlon_rdr,              &
     &                 'radata  ',6,x0rdr,y0rdr,cpjrdr)

! -----

! Rotate the x and the y components of velocity to the model grid.

          if(rdrvar(1:1).eq.'o'.and.rdrvar(2:2).eq.'o') then

            call s_getxy(iddx_rdr,iddy_rdr,'oo',1,nid_rdr,1,njd_rdr,    &
     &                   tmp1_rdr(1,1,1),tmp1_rdr(1,1,2))

            call s_xy2ll(idmpopt_rdr,idnspol_rdr,idtlon_rdr,            &
     &                   'radata  ',6,'oo',x0rdr,y0rdr,cpjrdr,          &
     &                   1,nid_rdr,1,njd_rdr,tmp1_rdr(1,1,1),           &
     &                   tmp1_rdr(1,1,2),tmp1_rdr(1,1,3),lon_rdr)

            if(rotopt_rdr.eq.1) then

              call rotuvm2s(idmpopt_rdr,idnspol_rdr,idtlon_rdr,         &
     &                      1,nid_rdr,1,njd_rdr,1,nkd_rdr,cpjrdr,       &
     &                      1,nid_rdr_sub,1,njd_rdr_sub,1,nkd_rdr_sub,  &
     &                      lon_rdr,u_rdr,v_rdr)

            end if

            call rotuvs2m(idmpopt,idnspol,idtlon,cpj,                   &
     &                    nid_rdr,njd_rdr,nkd_rdr,lon_rdr,u_rdr,v_rdr)

          end if

! -----

! Convert the reflection intensity to the precipitation mixing ratio.

          if(rdrvar(4:4).eq.'o'.and.datype_rdr(1:1).eq.'r') then

            call dbz2kg(idrdrcoe_rdr,idrdrexp_rdr,                      &
     &                  nid_rdr,njd_rdr,nkd_rdr,qp_rdr)

          end if

! -----

! Interpolate the radar data to the flat plane vertically.

          call vintrdr(idrdrvar,nid_rdr,njd_rdr,nkd_rdr,z_rdr,nk,z,     &
     &                 tmp1_rdr,km_rdr,u_rdr,v_rdr,w_rdr,qp_rdr)

! -----

! Interpolate the radar data to the model grid points horizontally.

          call intrprdr(idrdrvar,x0,y0,cpj,x0rdr,y0rdr,cpjrdr,ni,nj,nk, &
     &                  z,zph,urdr,vrdr,wrdr,qprdr,tmp1,tmp2,tmp3,tmp4, &
     &                  tmp5,tmp6,nid_rdr,njd_rdr,u_rdr,v_rdr,w_rdr,    &
     &                  qp_rdr)

! -----

! Set the bondary conditions for radar data.

          call bcrdr(idrdrvar,idwbc,idebc,idexbopt,idngropt,ni,nj,nk,   &
     &               urdr,vrdr,wrdr,qprdr)

! -----

! Distribute the observed precipitation to the rain, snow and graupel
! mixing ratio.

          if(rdrvar(4:4).eq.'o'.and.abs(cphopt).lt.10) then

            call distrqp(iddatype_rdr,idcphopt,idhaiopt,idthresq,       &
     &                 ni,nj,nk,nqw,nqi,rbr,qwtr,qice,qprdr,qwrdr,qirdr)

          end if

! -----

! Calculate the next time to read out.

          nxtrdr(2)=nxtrdr(2)+int(rdritv+.1e0)

! -----

!! -----

! Set the error frdr.

        else

          frdr(2)=-1

        end if

! -----

      end if

!!! -----

      end subroutine s_setrdr

!-----7--------------------------------------------------------------7--

      end module m_setrdr
