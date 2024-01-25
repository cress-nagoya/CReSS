!***********************************************************************
      module m_rdrdrv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/09/09
!     Modification: 2002/10/31, 2002/12/02, 2003/04/30, 2003/05/19,
!                   2003/09/01, 2003/12/12, 2004/03/05, 2004/04/15,
!                   2004/05/07, 2004/05/31, 2004/06/10, 2004/08/20,
!                   2004/09/01, 2005/02/10, 2005/04/04, 2005/08/05,
!                   2006/09/21, 2006/09/30, 2006/11/06, 2007/01/20,
!                   2007/07/30, 2008/04/17, 2008/05/02, 2008/08/25,
!                   2008/10/10, 2009/01/05, 2009/02/27, 2013/02/13,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the inferior procedures to interpolate the radar data.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_chkrdr
      use m_chkstd
      use m_comindx
      use m_comkind
      use m_commpi
      use m_currpe
      use m_dbz2kg
      use m_destroy
      use m_getcname
      use m_getdate
      use m_getiname
      use m_getrname
      use m_getxy
      use m_getz
      use m_getzph
      use m_inichar
      use m_intrprdr
      use m_outrdr
      use m_outstd05
      use m_rdradar
      use m_rdrstep
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

      public :: rdrdrv, s_rdrdrv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rdrdrv

        module procedure s_rdrdrv

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic int
      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_rdrdrv(fprdrvar,fpdatype_rdr,fpidate,                &
     &                    fprotopt_rdr,fprdritv,fpngrstr,ni,nj,nk,      &
     &                    z,zph,u,v,w,qp,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, &
     &                    nid_rdr,njd_rdr,nkd_rdr,km_rdr,londat,zdat,   &
     &                    udat,vdat,wdat,qpdat,dtmp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: fprdrvar
                       ! Formal parameter of unique index of rdrvar

      integer, intent(in) :: fpdatype_rdr
                       ! Formal parameter of unique index of datype_rdr

      integer, intent(in) :: fpidate
                       ! Formal parameter of unique index of idate

      integer, intent(in) :: fprotopt_rdr
                       ! Formal parameter of unique index of rotopt_rdr

      integer, intent(in) :: fprdritv
                       ! Formal parameter of unique index of rdritv

      integer, intent(in) :: fpngrstr
                       ! Formal parameter of unique index of ngrstr

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nid_rdr
                       ! Radar data dimension in x direction

      integer, intent(in) :: njd_rdr
                       ! Radar data dimension in y direction

      integer, intent(in) :: nkd_rdr
                       ! Radar data dimension in z direction

      integer, intent(in) :: km_rdr
                       ! Dimension of max(nk, nkd_rdr)

! Internal shared variables

      character(len=108) rdrvar
                       ! Control flag of input radar data variables

      character(len=108) datype_rdr
                       ! Control flag of radar data type

      character(len=108) idate
                       ! Forecast start date
                       ! with Gregorian calendar, yyyymmddhhmm

      character(len=12) cdate
                       ! Current forecast date
                       ! with Gregorian calendar, yyyymmddhhmm

      integer rotopt_rdr
                       ! Option for
                       ! rotation of wind direction in radar data

      integer(kind=i8) it
                       ! Index of main do loop

      integer(kind=i8) nstp0
                       ! Start index of main do loop

      integer(kind=i8) nstp1
                       ! End index of main do loop

      integer(kind=i8) ctime
                       ! Model current forecast time

      integer(kind=i8) cmin
                       ! 60000 x (ctime / 60000)

      integer(kind=i8) csec
                       ! mod(ctime, 60000)

      integer(kind=i8) rdr103
                       ! 1000 x int(rdritv + 0.1)

      integer(kind=i8) str103
                       ! 1000 x int(ngrstr + 0.1)

      integer stat     ! Runtime status

      integer nid_rdr_sub
                       ! Substitute for nid_rdr

      integer njd_rdr_sub
                       ! Substitute for njd_rdr

      integer nkd_rdr_sub
                       ! Substitute for nkd_rdr

      real rdritv      ! Time interval of radar data file

      real ngrstr      ! Analysis nudging start time to radar data

      real x0          ! x origin of model grid
      real y0          ! y origin of model grid

      real cpj(1:7)    ! Map projection parameters of model grid

      real x0rdr       ! x origin of radar data grid
      real y0rdr       ! y origin of radar data grid

      real cpjrdr(1:7) ! Map projection parameters of radar data grid

      real, intent(inout) :: z(1:nk)
                       ! zeta coordinates

      real, intent(inout) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(inout) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(inout) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(inout) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

      real, intent(inout) :: qp(0:ni+1,0:nj+1,1:nk)
                       ! Precipitation mixing ratio

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

      real, intent(inout) :: londat(1:nid_rdr,1:njd_rdr)
                       ! Longitude in data

      real, intent(inout) :: zdat(1:nid_rdr,1:njd_rdr,1:nkd_rdr)
                       ! z physical coordinates in data

      real, intent(inout) :: udat(1:nid_rdr,1:njd_rdr,1:km_rdr)
                       ! x components of velocity in data

      real, intent(inout) :: vdat(1:nid_rdr,1:njd_rdr,1:km_rdr)
                       ! y components of velocity in data

      real, intent(inout) :: wdat(1:nid_rdr,1:njd_rdr,1:km_rdr)
                       ! z components of velocity in data

      real, intent(inout) :: qpdat(1:nid_rdr,1:njd_rdr,1:km_rdr)
                       ! Precipitation mixing ratio in data

      real, intent(inout) :: dtmp1(1:nid_rdr,1:njd_rdr,1:nk)
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
      call getiname(fprotopt_rdr,rotopt_rdr)
      call getrname(fprdritv,rdritv)
      call getrname(fpngrstr,ngrstr)

! -----

! Set the substituted variables.

      nid_rdr_sub=nid_rdr
      njd_rdr_sub=njd_rdr
      nkd_rdr_sub=nkd_rdr

! -----

! Read in the message to standard i/o.

      call outstd05(0)

! -----

! Set the map projection parameters of the model grid.

      call setproj(idmpopt,idnspol,iddx,iddy,idulat,idulon,idriu,idrju, &
     &             idtlat1,idtlat2,idtlon,'solver  ',6,x0,y0,cpj)

! -----

! Set the map projection parameters of the data grid.

      call setproj(idmpopt_rdr,idnspol_rdr,iddx_rdr,iddy_rdr,idulat_rdr,&
     &           idulon_rdr,idriu_rdr,idrju_rdr,idtlat1_rdr,idtlat2_rdr,&
     &           idtlon_rdr,'radata  ',6,x0rdr,y0rdr,cpjrdr)

! -----

!! Get the x and the y coordinates and calculate the latitude and the
!! longitude at the data grid points.

      if(rdrvar(1:1).eq.'o'.and.rdrvar(2:2).eq.'o') then

! Get the x and the y coordinates at the data grid points.

        call s_getxy(iddx_rdr,iddy_rdr,'oo',1,nid_rdr,1,njd_rdr,        &
     &               dtmp1(1,1,1),dtmp1(1,1,2))

! -----

! Calculate the latitude and the longitude at the data grid points.

        call s_xy2ll(idmpopt_rdr,idnspol_rdr,idtlon_rdr,'radata  ',6,   &
     &               'oo',x0rdr,y0rdr,cpjrdr,1,nid_rdr,1,njd_rdr,       &
     &               dtmp1(1,1,1),dtmp1(1,1,2),dtmp1(1,1,3),londat)

! -----

      end if

!! -----

! Get the zeta coordinates at the model grid points.

      call getz(iddz,idzsfc,nk,z)

! -----

! Calculate the number of steps of the main do loop.

      call rdrstep(idngropt,idrdritv,idngrstr,idngrend,idstime,idetime, &
     &             nstp0,nstp1)

! -----

!!! Create the model input files.

! Set the common used variables.

      rdr103=1000_i8*int(rdritv+.1e0,i8)
      str103=1000_i8*int(ngrstr+.1e0,i8)

! -----

      do it=nstp0,nstp1

! Calculate the current forecast date.

        ctime=rdr103*(it-1_i8)+str103

        cmin=60000_i8*(ctime/60000_i8)
        csec=mod(ctime,60000_i8)

        call getdate(idate,cmin,cdate)

! -----

! Read out the data from the radar data file.

        stat=0

        call rdradar(iddatdir,idrdrvar,idncdat,idwlngth,                &
     &               cdate,ctime,csec,stat,nid_rdr,njd_rdr,nkd_rdr,     &
     &               zdat,udat,vdat,wdat,qpdat)

        call chkerr(stat)

        if(stat.lt.0) then

          call destroy('rdradar ',7,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        call chkstd(root)

! -----

! Check the radar data variables.

        stat=0

        call chkrdr(idrdrvar,iddatype_rdr,ctime,stat,                   &
     &              nid_rdr,njd_rdr,nkd_rdr,zdat,udat,vdat,wdat,qpdat)

        call chkerr(stat)

        if(stat.lt.0) then

          call destroy('chkrdr  ',6,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        call chkstd(root)

! -----

! Rotate the x and the y components of velocity to the model grid.

        if(rdrvar(1:1).eq.'o'.and.rdrvar(2:2).eq.'o') then

          if(rotopt_rdr.eq.1) then

           call rotuvm2s(idmpopt_rdr,idnspol_rdr,idtlon_rdr,            &
     &                   1,nid_rdr,1,njd_rdr,1,nkd_rdr,cpjrdr,          &
     &                   1,nid_rdr_sub,1,njd_rdr_sub,1,nkd_rdr_sub,     &
     &                   londat,udat,vdat)

          end if

          call rotuvs2m(idmpopt,idnspol,idtlon,cpj,                     &
     &                  nid_rdr,njd_rdr,nkd_rdr,londat,udat,vdat)

        end if

! -----

! Convert the reflection intensity to the precipitation mixing ratio.

        if(rdrvar(4:4).eq.'o'.and.datype_rdr(1:1).eq.'r') then

          call dbz2kg(idrdrcoe_rdr,idrdrexp_rdr,nid_rdr,njd_rdr,nkd_rdr,&
     &                qpdat)

        end if

! -----

! Interpolate the radar data to the flat plane vertically.

        call vintrdr(idrdrvar,nid_rdr,njd_rdr,nkd_rdr,zdat,nk,z,dtmp1,  &
     &               km_rdr,udat,vdat,wdat,qpdat)

! -----

!! Create the model input file for each processor element.

        do mype=0,npe-1

! Calculate the current processor element number.

         call currpe('all     ',3,'unset')

! -----

! Calculate the z physical coordinates.

         call s_getzph(idtrnopt,idexbopt,ni,nj,nk,z,zph,tmp1,tmp5,tmp6)

! -----

! Interpolate the radar data to the model grid points horizontally.

         call intrprdr(idrdrvar,x0,y0,cpj,x0rdr,y0rdr,cpjrdr,ni,nj,nk,  &
     &                 z,zph,u,v,w,qp,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,    &
     &                 nid_rdr,njd_rdr,udat,vdat,wdat,qpdat)

! -----

! Read in the data to the interpolated radar file.

         call outrdr(idexprim,idcrsdir,idrdrvar,idncexp,idnccrs,        &
     &               idwlngth,it,nstp0,ctime,ni,nj,nk,u,v,w,qp)

! -----

        end do

!! -----

! Read in the message to standard i/o.

        call outstd05(0)

! -----

      end do

!!! -----

      end subroutine s_rdrdrv

!-----7--------------------------------------------------------------7--

      end module m_rdrdrv
