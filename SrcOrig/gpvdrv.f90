!***********************************************************************
      module m_gpvdrv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/04/18
!     Modification: 2000/12/18, 2001/01/15, 2001/03/13, 2001/04/15,
!                   2001/05/29, 2001/06/29, 2002/01/07, 2002/04/02,
!                   2002/04/09, 2002/06/18, 2002/07/03, 2002/07/15,
!                   2002/08/15, 2002/09/02, 2002/09/09, 2002/10/31,
!                   2003/01/20, 2003/04/30, 2003/05/19, 2003/09/01,
!                   2003/12/12, 2004/01/09, 2004/03/05, 2004/04/15,
!                   2004/05/07, 2004/05/31, 2004/06/10, 2004/07/01,
!                   2004/08/01, 2004/08/20, 2004/09/01, 2004/09/10,
!                   2005/02/10, 2005/04/04, 2005/08/05, 2006/02/03,
!                   2006/08/18, 2006/09/21, 2006/09/30, 2006/11/06,
!                   2006/11/27, 2007/01/20, 2007/03/23, 2007/06/27,
!                   2007/07/30, 2007/09/04, 2008/05/02, 2008/07/01,
!                   2008/08/19, 2008/08/25, 2008/10/10, 2008/12/11,
!                   2009/02/27, 2009/03/31, 2011/09/22, 2013/02/13,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the inferior procedures to interpolate the grid data.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_chkgpv
      use m_chkstd
      use m_comindx
      use m_comkind
      use m_commpi
      use m_currpe
      use m_destroy
      use m_get1d
      use m_getbase
      use m_getcname
      use m_getdate
      use m_getiname
      use m_getnews
      use m_getrname
      use m_getxy
      use m_getz
      use m_getz11
      use m_getzph
      use m_gpvstep
      use m_inichar
      use m_intrpgpv
      use m_outgpv
      use m_outstd05
      use m_pc2kg
      use m_rdgpv
      use m_rotuvm2s
      use m_rotuvs2m
      use m_setproj
      use m_sparprt
      use m_t2pt
      use m_trndamp
      use m_vintbase
      use m_vintgpv
      use m_xy2ll

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: gpvdrv, s_gpvdrv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface gpvdrv

        module procedure s_gpvdrv

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
      subroutine s_gpvdrv(fpgpvvar,fpidate,fpdatype_gpv,                &
     &                    fptrnopt,fpexbopt,fprotopt_gpv,fpgpvitv,      &
     &                    ni,nj,nk,land,z,zph,ubr,vbr,pbr,ptbr,qvbr,    &
     &                    u,v,w,pp,ptp,qv,qc,qr,qi,qs,qg,qh,            &
     &                    z1d,u1d,v1d,p1d,pt1d,qv1d,tmp1,tmp2,tmp3,tmp4,&
     &                    tmp5,tmp6,tmp7,nid_gpv,njd_gpv,nkd_gpv,km_gpv,&
     &                    londat,htdat,zdat,udat,vdat,wdat,pdat,ptdat,  &
     &                    qvdat,qcdat,qrdat,qidat,qsdat,qgdat,qhdat,    &
     &                    pbdat,ptbdat,dtmp1,dtmp2)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpgpvvar
                       ! Formal parameter of unique index of gpvvar

      integer, intent(in) :: fpidate
                       ! Formal parameter of unique index of idate

      integer, intent(in) :: fpdatype_gpv
                       ! Formal parameter of unique index of datype_gpv

      integer, intent(in) :: fptrnopt
                       ! Formal parameter of unique index of trnopt

      integer, intent(in) :: fpexbopt
                       ! Formal parameter of unique index of exbopt

      integer, intent(in) :: fprotopt_gpv
                       ! Formal parameter of unique index of rotopt_gpv

      integer, intent(in) :: fpgpvitv
                       ! Formal parameter of unique index of gpvitv

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nid_gpv
                       ! GPV data dimension in x direction

      integer, intent(in) :: njd_gpv
                       ! GPV data dimension in y direction

      integer, intent(in) :: nkd_gpv
                       ! GPV data dimension in z direction

      integer, intent(in) :: km_gpv
                       ! Dimension of max(nk, nkd_gpv)

! Internal shared variables

      character(len=108) gpvvar
                       ! Control flag of input GPV data variables

      character(len=108) datype_gpv
                       ! Control flag of GPV data type

      character(len=108) idate
                       ! Forecast start date
                       ! with Gregorian calendar, yyyymmddhhmm

      character(len=12) cdate
                       ! Current forecast date
                       ! with Gregorian calendar, yyyymmddhhmm

      integer trnopt   ! Option for terrain height setting
      integer exbopt   ! Option for external boundary forcing

      integer rotopt_gpv
                       ! Option for
                       ! rotation of wind direction in GPV data

      integer idstr    ! Minimum index
                       ! of model grid in data region in x direction

      integer idend    ! Maximum index
                       ! of model grid in data region in x direction

      integer jdstr    ! Minimum index
                       ! of model grid in data region in y direction

      integer jdend    ! Maximum index
                       ! of model grid in data region in y direction

      integer kref     ! Reference index

      integer(kind=i8) it
                       ! Index of main do loop

      integer(kind=i8) nstp0
                       ! Start index of main do loop

      integer(kind=i8) nstp1
                       ! End index of main do loop

      integer(kind=i8) ctime
                       ! Model current forecast time

      integer(kind=i8) gpv103
                       ! 1000 x int(gpvitv + 0.1)

      integer stat     ! Runtime status

      integer nid_gpv_sub
                       ! Substitute for nid_gpv

      integer njd_gpv_sub
                       ! Substitute for njd_gpv

      integer nkd_gpv_sub
                       ! Substitute for nkd_gpv

      integer, intent(inout) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

      real gpvitv      ! Time interval of GPV data file

      real x0          ! x origin of model grid
      real y0          ! y origin of model grid

      real cpj(1:7)    ! Map projection parameters of model grid

      real x0gpv       ! x origin of GPV data grid
      real y0gpv       ! y origin of GPV data grid

      real cpjgpv(1:7) ! Map projection parameters of GPV data grid

      real, intent(inout) :: z(1:nk)
                       ! zeta coordinates

      real, intent(inout) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(inout) :: ubr(0:ni+1,0:nj+1,1:nk)
                       ! Base state x components of velocity

      real, intent(inout) :: vbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state y components of velocity

      real, intent(inout) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(inout) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(inout) :: qvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state water vapor mixing ratio

      real, intent(inout) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(inout) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(inout) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

      real, intent(inout) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation

      real, intent(inout) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

      real, intent(inout) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

      real, intent(inout) :: qc(0:ni+1,0:nj+1,1:nk)
                       ! Cloud water mixing ratio

      real, intent(inout) :: qr(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixing ratio

      real, intent(inout) :: qi(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixing ratio

      real, intent(inout) :: qs(0:ni+1,0:nj+1,1:nk)
                       ! Snow mixing ratio

      real, intent(inout) :: qg(0:ni+1,0:nj+1,1:nk)
                       ! Graupel mixing ratio

      real, intent(inout) :: qh(0:ni+1,0:nj+1,1:nk)
                       ! Hail mixing ratio

      real, intent(inout) :: z1d(0:4*nk-3)
                       ! Horizontally averaged z coordinates

      real, intent(inout) :: u1d(0:4*nk-3)
                       ! Horizontally averaged x components of velocity

      real, intent(inout) :: v1d(0:4*nk-3)
                       ! Horizontally averaged y components of velocity

      real, intent(inout) :: p1d(0:4*nk-3)
                       ! Horizontally averaged pressure

      real, intent(inout) :: pt1d(0:4*nk-3)
                       ! Horizontally averaged potential temperature

      real, intent(inout) :: qv1d(0:4*nk-3)
                       ! Horizontally averaged water vapor mixing ratio

      real, intent(inout) :: tmp1(1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp4(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp5(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp6(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp7(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: londat(1:nid_gpv,1:njd_gpv)
                       ! Longitude in data

      real, intent(inout) :: htdat(1:nid_gpv,1:njd_gpv)
                       ! Terrain height in data

      real, intent(inout) :: zdat(1:nid_gpv,1:njd_gpv,1:nkd_gpv)
                       ! z physical coordinates in data

      real, intent(inout) :: udat(1:nid_gpv,1:njd_gpv,1:km_gpv)
                       ! x components of velocity in data

      real, intent(inout) :: vdat(1:nid_gpv,1:njd_gpv,1:km_gpv)
                       ! y components of velocity in data

      real, intent(inout) :: wdat(1:nid_gpv,1:njd_gpv,1:km_gpv)
                       ! z components of velocity in data

      real, intent(inout) :: pdat(1:nid_gpv,1:njd_gpv,1:km_gpv)
                       ! Pressure in data

      real, intent(inout) :: ptdat(1:nid_gpv,1:njd_gpv,1:km_gpv)
                       ! Potential temperature in data

      real, intent(inout) :: qvdat(1:nid_gpv,1:njd_gpv,1:km_gpv)
                       ! Water vapor mixing ratio in data

      real, intent(inout) :: qcdat(1:nid_gpv,1:njd_gpv,1:km_gpv)
                       ! Cloud water mixing ratio in data

      real, intent(inout) :: qrdat(1:nid_gpv,1:njd_gpv,1:km_gpv)
                       ! Rain water mixing ratio in data

      real, intent(inout) :: qidat(1:nid_gpv,1:njd_gpv,1:km_gpv)
                       ! Cloud ice mixing ratio in data

      real, intent(inout) :: qsdat(1:nid_gpv,1:njd_gpv,1:km_gpv)
                       ! Snow mixing ratio in data

      real, intent(inout) :: qgdat(1:nid_gpv,1:njd_gpv,1:km_gpv)
                       ! Graupel mixing ratio in data

      real, intent(inout) :: qhdat(1:nid_gpv,1:njd_gpv,1:km_gpv)
                       ! Hail mixing ratio in data

      real, intent(inout) :: pbdat(1:nid_gpv,1:njd_gpv,1:km_gpv)
                       ! Base state pressure in data

      real, intent(inout) :: ptbdat(1:nid_gpv,1:njd_gpv,1:km_gpv)
                       ! Base state potential temperature in data

      real, intent(inout) :: dtmp1(1:nid_gpv,1:njd_gpv)
                       ! Temporary array

      real, intent(inout) :: dtmp2(1:nid_gpv,1:njd_gpv,1:km_gpv)
                       ! Temporary array

!-----7--------------------------------------------------------------7--

! Initialize character variables.

      call inichar(gpvvar)
      call inichar(datype_gpv)
      call inichar(idate)

! -----

! Get the required namelist variables.

      call getcname(fpgpvvar,gpvvar)
      call getcname(fpidate,idate)
      call getcname(fpdatype_gpv,datype_gpv)
      call getiname(fptrnopt,trnopt)
      call getiname(fpexbopt,exbopt)
      call getiname(fprotopt_gpv,rotopt_gpv)
      call getrname(fpgpvitv,gpvitv)

! -----

! Set the substituted variables.

      nid_gpv_sub=nid_gpv
      njd_gpv_sub=njd_gpv
      nkd_gpv_sub=nkd_gpv

! -----

! Read in the message to standard i/o.

      call outstd05(0)

! -----

! Set the map projection parameters of the model grid.

      call setproj(idmpopt,idnspol,iddx,iddy,idulat,idulon,idriu,idrju, &
     &             idtlat1,idtlat2,idtlon,'solver  ',6,x0,y0,cpj)

! -----

! Set the map projection parameters of the data grid.

      call setproj(idmpopt_gpv,idnspol_gpv,iddx_gpv,iddy_gpv,idulat_gpv,&
     &           idulon_gpv,idriu_gpv,idrju_gpv,idtlat1_gpv,idtlat2_gpv,&
     &           idtlon_gpv,'gridata ',7,x0gpv,y0gpv,cpjgpv)

! -----

! Get the x and the y coordinates at the data grid points.

      call s_getxy(iddx_gpv,iddy_gpv,'oo',1,nid_gpv,1,njd_gpv,          &
     &             dtmp2(1,1,1),dtmp2(1,1,2))

! -----

! Calculate the latitude and the longitude at the data grid points.

      call s_xy2ll(idmpopt_gpv,idnspol_gpv,idtlon_gpv,'gridata ',7,     &
     &             'oo',x0gpv,y0gpv,cpjgpv,1,nid_gpv,1,njd_gpv,         &
     &             dtmp2(1,1,1),dtmp2(1,1,2),dtmp1,londat)

! -----

! Get the mininum and maximum data indices to create the base state
! variables.

      do mype=0,npe-1

        call currpe('all     ',3,'unset')

        call s_getnews(nid_gpv,njd_gpv,idstr,idend,jdstr,jdend,         &
     &                 x0,y0,cpj,x0gpv,y0gpv,cpjgpv,ni,nj,              &
     &                 tmp2,tmp3,tmp4,tmp5,tmp6,tmp7)

      end do

! -----

! Get the zeta coordinates at the model grid points.

      call getz11(iddz,idzsfc,nk,z)

! -----

! Calculate the number of steps of the main do loop.

      call gpvstep(idnggopt,idexbopt,idlspopt,idvspopt,idgpvitv,        &
     &             idstime,idetime,nstp0,nstp1)

! -----

!! Get the base state variables in the case of restart running.

      if(nstp0.ge.2_i8) then

! Calculate the current forecast date.

        ctime=0_i8

        call getdate(idate,ctime,cdate)

! -----

! Read out the data from the GPV data file.

        stat=0

        call rdgpv(iddatdir,idgpvvar,idncdat,idwlngth,idexbopt,         &
     &             cdate,ctime,stat,nid_gpv,njd_gpv,nkd_gpv,            &
     &             htdat,zdat,udat,vdat,wdat,pdat,ptdat,qvdat,          &
     &             qcdat,qrdat,qidat,qsdat,qgdat,qhdat)

        call chkerr(stat)

        if(stat.lt.0) then

          call destroy('rdgpv   ',5,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        call chkstd(root)

! -----

! Check the GPV data variables.

        stat=0

        call chkgpv(idgpvvar,iddatype_gpv,idexbopt,ctime,stat,          &
     &             nid_gpv,njd_gpv,nkd_gpv,htdat,zdat,udat,vdat,wdat,   &
     &             pdat,ptdat,qvdat,qcdat,qrdat,qidat,qsdat,qgdat,qhdat)

        call chkerr(stat)

        if(stat.lt.0) then

          call destroy('chkgpv  ',6,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        call chkstd(root)

! -----

! Rotate the x and the y components of velocity to the model grid.

        if(rotopt_gpv.eq.1) then

          call rotuvm2s(idmpopt_gpv,idnspol_gpv,idtlon_gpv,             &
     &                  1,nid_gpv,1,njd_gpv,1,nkd_gpv,cpjgpv,           &
     &                  1,nid_gpv_sub,1,njd_gpv_sub,1,nkd_gpv_sub,      &
     &                  londat,udat,vdat)

        end if

        call rotuvs2m(idmpopt,idnspol,idtlon,cpj,                       &
     &                nid_gpv,njd_gpv,nkd_gpv,londat,udat,vdat)

! -----

! Convert the temperature to the potential temperature.

        if(datype_gpv(1:1).eq.'t') then

          call t2pt(nid_gpv,njd_gpv,nkd_gpv,pdat,ptdat)

        end if

! -----

! Convert the relative humidity to the water vapor mixing ratio.

        if(gpvvar(2:2).eq.'o'.and.datype_gpv(2:2).eq.'r') then

          call pc2kg(nid_gpv,njd_gpv,nkd_gpv,pdat,ptdat,qvdat)

        end if

! -----

! Interpolating the GPV data to the flat plane to create the base state
! variables.

        call vintbase(idgpvvar,idrefsfc_gpv,idstr,idend,jdstr,jdend,    &
     &                kref,nid_gpv,njd_gpv,nkd_gpv,zdat,dtmp1,nk,z,     &
     &                km_gpv,pdat,ptdat,udat,vdat,qvdat,pbdat,ptbdat,   &
     &                dtmp2)

! -----

! Get the horizontally averaged variables.

        call get1d(idgpvvar,idetrvar_gpv,'all   ',                      &
     &             idstr,idend,jdstr,jdend,kref,nid_gpv,njd_gpv,nkd_gpv,&
     &             zdat,nk,z,udat,vdat,pbdat,ptbdat,qvdat,              &
     &             z1d,u1d,v1d,p1d,pt1d,qv1d,tmp1)

! -----

! Read in the message to standard i/o.

        call outstd05(0)

! -----

      end if

!! -----

!!! Create the model input files.

! Set the common used variable.

      gpv103=1000_i8*int(gpvitv+.1e0,i8)

! -----

      do it=nstp0,nstp1

! Calculate the current forecast date.

        ctime=gpv103*(it-1_i8)

        call getdate(idate,ctime,cdate)

! -----

! Read out the data from the GPV data file.

        stat=0

        call rdgpv(iddatdir,idgpvvar,idncdat,idwlngth,idexbopt,         &
     &             cdate,ctime,stat,nid_gpv,njd_gpv,nkd_gpv,            &
     &             htdat,zdat,udat,vdat,wdat,pdat,ptdat,qvdat,          &
     &             qcdat,qrdat,qidat,qsdat,qgdat,qhdat)

        call chkerr(stat)

        if(stat.lt.0) then

          call destroy('rdgpv   ',5,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        call chkstd(root)

! -----

! Check the GPV data variables.

        stat=0

        call chkgpv(idgpvvar,iddatype_gpv,idexbopt,ctime,stat,          &
     &             nid_gpv,njd_gpv,nkd_gpv,htdat,zdat,udat,vdat,wdat,   &
     &             pdat,ptdat,qvdat,qcdat,qrdat,qidat,qsdat,qgdat,qhdat)

        call chkerr(stat)

        if(stat.lt.0) then

          call destroy('chkgpv  ',6,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        call chkstd(root)

! -----

! Rotate the x and the y components of velocity to the model grid.

        if(rotopt_gpv.eq.1) then

          call rotuvm2s(idmpopt_gpv,idnspol_gpv,idtlon_gpv,             &
     &                  1,nid_gpv,1,njd_gpv,1,nkd_gpv,cpjgpv,           &
     &                  1,nid_gpv_sub,1,njd_gpv_sub,1,nkd_gpv_sub,      &
     &                  londat,udat,vdat)

        end if

        call rotuvs2m(idmpopt,idnspol,idtlon,cpj,                       &
     &                nid_gpv,njd_gpv,nkd_gpv,londat,udat,vdat)

! -----

! Convert the temperature to the potential temperature.

        if(datype_gpv(1:1).eq.'t') then

          call t2pt(nid_gpv,njd_gpv,nkd_gpv,pdat,ptdat)

        end if

! -----

! Convert the relative humidity to the water vapor mixing ratio.

        if(gpvvar(2:2).eq.'o'.and.datype_gpv(2:2).eq.'r') then

          call pc2kg(nid_gpv,njd_gpv,nkd_gpv,pdat,ptdat,qvdat)

        end if

! -----

! Interpolating the GPV data to the flat plane to create the base state
! variables.

        if(it.eq.1_i8) then

          call vintbase(idgpvvar,idrefsfc_gpv,idstr,idend,jdstr,jdend,  &
     &                  kref,nid_gpv,njd_gpv,nkd_gpv,zdat,dtmp1,nk,z,   &
     &                  km_gpv,pdat,ptdat,udat,vdat,qvdat,pbdat,ptbdat, &
     &                  dtmp2)

        end if

! -----

! Get the horizontally averaged variables.

        if(it.eq.1_i8) then

          call get1d(idgpvvar,idetrvar_gpv,'all   ',                    &
     &             idstr,idend,jdstr,jdend,kref,nid_gpv,njd_gpv,nkd_gpv,&
     &             zdat,nk,z,udat,vdat,pbdat,ptbdat,qvdat,              &
     &             z1d,u1d,v1d,p1d,pt1d,qv1d,tmp1)

        end if

! -----

! Separate the pressure and the potential temperature to the base state
! and the perturbation value.

        call sparprt(nid_gpv,njd_gpv,nkd_gpv,zdat,pdat,ptdat,           &
     &               pbdat,ptbdat,nk,z1d,p1d,pt1d)

! -----

! Interpolate the GPV data to the flat plane vertically.

        call vintgpv(idgpvvar,idrefsfc_gpv,it,nid_gpv,njd_gpv,nkd_gpv,  &
     &               zdat,dtmp1,nk,z,z1d,p1d,pt1d,dtmp2,pbdat,ptbdat,   &
     &               km_gpv,udat,vdat,wdat,pdat,ptdat,qvdat,            &
     &               qcdat,qrdat,qidat,qsdat,qgdat,qhdat)

! -----

! Calculate the damped terrain height.

        if(it.eq.nstp0.and.trnopt.eq.2.and.mod(exbopt,10).eq.2) then

          do mype=0,npe-1

            call currpe('all     ',3,'unset')

            call s_trndamp(x0,y0,cpj,x0gpv,y0gpv,cpjgpv,ni,nj,land,zph, &
     &                     tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,               &
     &                     nid_gpv,njd_gpv,htdat)

          end do

        end if

! -----

!! Create the model input file for each processor element.

        do mype=0,npe-1

! Calculate the current processor element number.

         call currpe('all     ',3,'unset')

! -----

! Calculate the z physical coordinates.

         call getz(iddz,idzsfc,nk,tmp1)

         call s_getzph(idtrnopt,idexbopt,ni,nj,nk,tmp1,zph,             &
     &                 tmp2,tmp6,tmp7)

! -----

! Extract the base state variables for the model grid from the
! interpolated GPV data.

         if(it.eq.1_i8) then

           call getbase(idgpvvar,ni,nj,nk,zph,z1d,u1d,v1d,              &
     &                  p1d,pt1d,qv1d,ubr,vbr,pbr,ptbr,qvbr,tmp6)

         end if

! -----

! Interpolate the GPV data to the model grid points horizontally.

         call intrpgpv(idgpvvar,x0,y0,cpj,x0gpv,y0gpv,cpjgpv,           &
     &                 ni,nj,nk,z,zph,u,v,w,pp,ptp,qv,qc,qr,qi,qs,qg,qh,&
     &                 tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,nid_gpv,njd_gpv,   &
     &                 udat,vdat,wdat,pdat,ptdat,qvdat,qcdat,qrdat,     &
     &                 qidat,qsdat,qgdat,qhdat)

! -----

! Read in the data to the interpolated GPV file.

         call outgpv(idexprim,idcrsdir,idgpvvar,idncexp,idnccrs,        &
     &               idwlngth,it,nstp0,ctime,ni,nj,nk,ubr,vbr,          &
     &               pbr,ptbr,qvbr,u,v,w,pp,ptp,qv,qc,qr,qi,qs,qg,qh)

! -----

        end do

!! -----

! Read in the message to standard i/o.

        call outstd05(0)

! -----

      end do

!!! -----

      end subroutine s_gpvdrv

!-----7--------------------------------------------------------------7--

      end module m_gpvdrv
