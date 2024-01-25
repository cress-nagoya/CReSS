!***********************************************************************
      module m_phycood
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/03/25, 1999/04/06, 1999/07/05, 1999/08/23,
!                   1999/10/12, 1999/11/01, 1999/11/19, 2000/01/17,
!                   2000/07/05, 2000/12/18, 2001/01/15, 2001/04/15,
!                   2001/05/29, 2001/10/17, 2001/12/11, 2002/04/02,
!                   2002/09/09, 2003/04/30, 2003/05/19, 2003/12/12,
!                   2004/08/01, 2004/09/10, 2005/02/10, 2005/08/05,
!                   2006/01/10, 2007/01/20, 2007/01/31, 2007/05/14,
!                   2007/10/19, 2008/05/02, 2008/06/09, 2008/08/25,
!                   2009/01/05, 2009/02/27, 2009/11/13, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the 3 dimensional z physical coordinates.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_comindx
      use m_commath
      use m_commpi
      use m_copy1d
      use m_cpondpe
      use m_destroy
      use m_getiname
      use m_getrname
      use m_stretch

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: phycood, s_phycood

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface phycood

        module procedure s_phycood

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic max
      intrinsic min

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_phycood(fpsthopt,fpzsfc,fpzflat,ni,nj,nk,z,ht,       &
     &                     zph,zsth,dzsth)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpsthopt
                       ! Formal parameter of unique index of sthopt

      integer, intent(in) :: fpzsfc
                       ! Formal parameter of unique index of zsfc

      integer, intent(in) :: fpzflat
                       ! Formal parameter of unique index of zflat

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: z(1:nk)
                       ! zeta coordinates

      real, intent(in) :: ht(0:ni+1,0:nj+1)
                       ! Terrain height

! Output variables

      real, intent(out) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(out) :: zsth(1:nk)
                       ! 1 dimensional stretched z coordinates

! Internal shared variables

      integer sthopt   ! Option for vertical grid stretching

      integer nkm1     ! nk - 1
      integer nkm2     ! nk - 2

      integer kflat    ! Index of lowest flat level

      integer stat     ! Runtime status

      real zsfc        ! Sea surface terrain height

      real zflat       ! Lowest flat level

      real zflat0      ! Re-calculated lowest flat level

      real htmax       ! Highest terrain height

      real htuiv       ! Inverse of distance of surface to flat level

      real htuivz      ! htuiv x zflat0

      real, intent(inout) :: dzsth(1:nk)
                       ! Distance of
                       ! 1 dimensional stretched z coordinates

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpsthopt,sthopt)
      call getrname(fpzsfc,zsfc)
      call getrname(fpzflat,zflat)

! -----

! Copy the z coordinates to the array zsth in the case the vertical
! stretching is not applied.

      if(sthopt.eq.0) then

        call copy1d(1,nk,z,zsth)

! -----

! Calculate the 1 dimensional physical coordinates.

      else

        call stretch(idsthopt,idzsfc,iddzmin,idlayer1,idlayer2,z(nk-2), &
     &               z(nk-1),nk,zsth,dzsth)

      end if

! -----

!! Check the highest mountain height.

! Initialize the processed variable, htmax.

      htmax=lim36n

! -----

! Get the highest mountain height.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j) reduction(max: htmax)

      do j=0,nj
      do i=0,ni
        htmax=max(ht(i,j),htmax)
      end do
      end do

!$omp end do

!$omp end parallel

! -----

! If error occured, call the procedure destroy.

      stat=0

      if(zflat.lt.1.5e0*htmax) then
        stat=1
      end if

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('phycood ',7,'cont',7,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('phycood ',7,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

!! -----

!! Calculate the lowest flat level.

! Initialize the processed variable, kflat.

      kflat=nk

! -----

! Get the index of lowest flat level.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(k) reduction(min: kflat)

      do k=2,nk-1

        if(zsth(k).gt.zflat) then

          kflat=min(k,kflat)

        end if

      end do

!$omp end do

!$omp end parallel

! -----

! Finally get the lowest flat level.

      zflat0=min(zsth(kflat-1),z(nk-1))

      if(zflat0.le.zsfc) then
        zflat0=z(nk-1)
      end if

! -----

!! -----

!! Get the 3 dimensional z physical coordinates, reset the 1 dimensional
!! stretched z coordinates and set the bottom and top boundary
!! conditions.

! Set the common used variables.

      nkm1=nk-1
      nkm2=nk-2

      htuiv=1.e0/(zflat0-zsfc)
      htuivz=1.e0/(zflat0-zsfc)*zflat0

! -----

! Get the z physical coordinates and reset the stretched z coordinates.

!$omp parallel default(shared) private(k)

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=0,nj
        do i=0,ni

          if(zsth(k).gt.zflat0) then
            zph(i,j,k)=zsth(k)

          else
            zph(i,j,k)=htuiv*(zflat0-ht(i,j))*(zsth(k)-zsfc)+ht(i,j)

          end if

        end do
        end do

!$omp end do

      end do

!$omp do schedule(runtime) private(i,j)

      do j=0,nj
      do i=0,ni
        zph(i,j,2)=ht(i,j)
        zph(i,j,1)=2.e0*zph(i,j,2)-zph(i,j,3)
        zph(i,j,nk)=2.e0*zph(i,j,nkm1)-zph(i,j,nkm2)
      end do
      end do

!$omp end do

!$omp do schedule(runtime)

      do k=2,nk-1

        if(zsth(k).le.zflat0) then
          zsth(k)=htuivz*(zsth(k)-zsfc)
        end if

      end do

!$omp end do

!$omp end parallel

! -----

! Set the bottom and top boundary conditions.

      zsth(2)=0.e0
      zsth(1)=-zsth(3)
      zsth(nk)=2.e0*zsth(nk-1)-zsth(nk-2)

! -----

!! -----

      end subroutine s_phycood

!-----7--------------------------------------------------------------7--

      end module m_phycood
