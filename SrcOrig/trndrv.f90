!***********************************************************************
      module m_trndrv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/04/18
!     Modification: 2000/06/01, 2000/12/18, 2001/05/29, 2002/06/18,
!                   2002/07/03, 2002/09/09, 2003/04/30, 2003/05/19,
!                   2003/06/27, 2003/12/12, 2004/04/15, 2004/05/07,
!                   2004/05/31, 2004/06/10, 2005/02/10, 2005/04/04,
!                   2006/01/10, 2006/11/06, 2007/01/20, 2007/08/24,
!                   2008/05/02, 2008/08/25, 2008/10/10, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the inferior procedures to interpolate the terrain data.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_chkmxn
      use m_chkstd
      use m_comindx
      use m_comkind
      use m_commpi
      use m_currpe
      use m_destroy
      use m_getrij
      use m_hint2d
      use m_outstd05
      use m_outstd12
      use m_outtrn
      use m_rdheight
      use m_setproj

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: trndrv, s_trndrv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface trndrv

        module procedure s_trndrv

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
      subroutine s_trndrv(ni,nj,ri,rj,ht,tmp1,tmp2,tmp3,nid_trn,njd_trn,&
     &                    htdat)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nid_trn
                       ! Terrain data dimension in x direction

      integer, intent(in) :: njd_trn
                       ! Terrain data dimension in y direction

! Internal shared variables

      integer stat     ! Runtime status

      real x0          ! x origin of model grid
      real y0          ! y origin of model grid

      real cpj(1:7)    ! Map projection parameters of model grid

      real x0trn       ! x origin of terrain data grid
      real y0trn       ! y origin of terrain data grid

      real cpjtrn(1:7) ! Map projection parameters of terrain data grid

      real, intent(inout) :: ri(0:ni+1,0:nj+1)
                       ! Real indices in data region in x direction

      real, intent(inout) :: rj(0:ni+1,0:nj+1)
                       ! Real indices in data region in y direction

      real, intent(inout) :: ht(0:ni+1,0:nj+1)
                       ! Terrain height

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: htdat(1:nid_trn,1:njd_trn)
                       ! Terrain height in data

! Remark

!     ht: This variable is also temporary.

!-----7--------------------------------------------------------------7--

! Set the map projection parameters of the model grid.

      call setproj(idmpopt,idnspol,iddx,iddy,idulat,idulon,idriu,idrju, &
     &             idtlat1,idtlat2,idtlon,'solver  ',6,x0,y0,cpj)

! -----

! Set the map projection parameters of the data grid.

      call setproj(idmpopt_trn,idnspol_trn,iddx_trn,iddy_trn,idulat_trn,&
     &           idulon_trn,idriu_trn,idrju_trn,idtlat1_trn,idtlat2_trn,&
     &           idtlon_trn,'terrain ',7,x0trn,y0trn,cpjtrn)

! -----

! Read out the data from the terrain data file.

      stat=0

      call rdheight(iddatdir,idncdat,idwlngth,stat,nid_trn,njd_trn,     &
     &              htdat)

      call chkerr(stat)

      if(stat.lt.0) then

        call destroy('rdheight',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      call chkstd(root)

! -----

! Check the terrain data file.

      stat=0

      call outstd12(0,'    ',4,'       ',0_i8,0,0,1,0.e0,0,0,1,0.e0)

      call s_chkmxn('ht  ',2,'[m]    ','all  ',0_i8,stat,               &
     &            -1000.e0,1000000.e0,0.e0,0.e0,nid_trn,njd_trn,1,htdat)

      call chkerr(stat)

      if(stat.lt.0) then

        call destroy('chkmxn  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      call chkstd(root)

! -----

!! Create the model input file for each processor element.

      do mype=0,npe-1

! Get the current processor element number.

        call currpe('all     ',3,'unset')

! -----

! Calculate the real indices at the model grid points in the data
! region.

        call getrij(idmpopt_trn,idnspol_trn,idtlon_trn,                 &
     &              iddxiv_trn,iddyiv_trn,'terrain ',7,'xx',            &
     &              x0,y0,cpj,x0trn,y0trn,cpjtrn,ni,nj,ri,rj,ht,        &
     &              tmp1,tmp2,tmp3)

! -----

! Interpolate the terrain data to the model grid horizontally.

        call hint2d(idmpopt_trn,idintopt_trn,'xx','cal ',ni,nj,ri,rj,   &
     &              tmp1,tmp2,ht,nid_trn,njd_trn,htdat)

! -----

! Read in the data to the interpolated terrain file.

        call outtrn(idexprim,idcrsdir,idncexp,idnccrs,idwlngth,         &
     &              'terrain     ',7,0,ni,nj,ht)

! -----

      end do

!! -----

! Read in the message to standard i/o.

      call outstd05(0)

! -----

      end subroutine s_trndrv

!-----7--------------------------------------------------------------7--

      end module m_trndrv
