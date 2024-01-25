!***********************************************************************
      module m_rbcv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/07/28,
!                   1999/08/03, 1999/08/18, 1999/08/23, 1999/09/06,
!                   1999/09/30, 1999/10/07, 1999/10/27, 1999/11/01,
!                   1999/12/06, 2000/01/17, 2000/03/17, 2000/03/23,
!                   2001/01/15, 2001/04/15, 2001/04/25, 2001/07/13,
!                   2001/08/07, 2001/09/13, 2001/12/11, 2002/04/02,
!                   2002/06/06, 2002/07/23, 2002/08/15, 2002/10/31,
!                   2003/04/30, 2003/05/19, 2003/06/27, 2003/11/05,
!                   2003/11/28, 2003/12/12, 2004/01/09, 2004/05/07,
!                   2004/08/20, 2005/01/31, 2006/09/21, 2006/12/04,
!                   2007/01/05, 2007/03/10, 2007/05/07, 2007/05/21,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2008/12/11,
!                   2009/02/27, 2009/03/23, 2011/09/22, 2013/01/28,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the radiative lateral boundary conditions for the y components
!     of velocity.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commpi
      use m_getcname
      use m_getiname
      use m_getrname
      use m_inichar

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: rbcv, s_rbcv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rbcv

        module procedure s_rbcv

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic mod
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_rbcv(fplbcvar,fpwbc,fpebc,fpsbc,fpnbc,fpnggopt,      &
     &                  fplspopt,fpvspopt,fplbnews,fplbnorm,            &
     &                  isstp,dts,gtinc,ni,nj,nk,vbr,                   &
     &                  vcpx,vcpy,vgpv,vtd,v)
!***********************************************************************

! Input variables

      integer, intent(in) :: fplbcvar
                       ! Formal parameter of unique index of lbcvar

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpsbc
                       ! Formal parameter of unique index of sbc

      integer, intent(in) :: fpnbc
                       ! Formal parameter of unique index of nbc

      integer, intent(in) :: fpnggopt
                       ! Formal parameter of unique index of nggopt

      integer, intent(in) :: fplspopt
                       ! Formal parameter of unique index of lspopt

      integer, intent(in) :: fpvspopt
                       ! Formal parameter of unique index of vspopt

      integer, intent(in) :: fplbnews
                       ! Formal parameter of unique index of lbnews

      integer, intent(in) :: fplbnorm
                       ! Formal parameter of unique index of lbnorm

      integer, intent(in) :: isstp
                       ! Index of small time steps integration

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dts
                       ! Small time steps interval

      real, intent(in) :: gtinc
                       ! Lapse of forecast time from GPV data reading

      real, intent(in) :: vbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state y components of velocity

      real, intent(in) :: vcpx(1:nj,1:nk,1:2)
                       ! Phase speed of y components of velocity
                       ! on west and east boundary

      real, intent(in) :: vcpy(1:ni,1:nk,1:2)
                       ! Phase speed of y components of velocity
                       ! on south and north boundary

      real, intent(in) :: vgpv(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: vtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! y components of velocity of GPV data

! Input and output variable

      real, intent(inout) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

! Internal shared variables

      character(len=108) lbcvar
                       ! Control flag of
                       ! lateral boundary forced variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions
      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

      integer nggopt   ! Option for analysis nudging to GPV
      integer lspopt   ! Option for lateral sponge damping
      integer vspopt   ! Option for vertical sponge damping

      integer nim1     ! ni - 1
      integer nim2     ! ni - 2
      integer njm1     ! nj - 1

      real lbnews      ! Boundary damping coefficient

      real lbnorm      ! Boundary damping coefficient
                       ! for u and v in normal

      real tdmpdt      ! lbnews x dts
      real ndmpdt      ! lbnorm x dts

      real tpdt        ! gtinc + real(isstp - 1) x dts

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(lbcvar)

! -----

! Get the required namelist variables.

      call getcname(fplbcvar,lbcvar)
      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpsbc,sbc)
      call getiname(fpnbc,nbc)
      call getiname(fpnggopt,nggopt)
      call getiname(fplspopt,lspopt)
      call getiname(fpvspopt,vspopt)
      call getrname(fplbnews,lbnews)
      call getrname(fplbnorm,lbnorm)

! -----

! Set the common used variables.

      nim1=ni-1
      nim2=ni-2
      njm1=nj-1

      if(lbcvar(2:2).eq.'o') then
        tdmpdt=lbnews*dts
        ndmpdt=lbnorm*dts
      else
        tdmpdt=0.e0
        ndmpdt=0.e0
      end if

      tpdt=gtinc+real(isstp-1)*dts

! -----

!! Set the radiative lateral boundary conditions.

!$omp parallel default(shared)

! Set the south boundary conditions.

      if(ebs.eq.1.and.jsub.eq.0) then

        if(sbc.ge.4) then

          if(nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

!$omp do schedule(runtime) private(i,k)

            do k=2,nk-2
            do i=1,ni-1
              v(i,1,k)=v(i,1,k)-vcpy(i,k,1)*(v(i,2,k)-v(i,1,k))         &
     &          -ndmpdt*(v(i,1,k)-(vgpv(i,1,k)+vtd(i,1,k)*tpdt))
            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(i,k)

            do k=2,nk-2
            do i=1,ni-1
              v(i,1,k)=v(i,1,k)-vcpy(i,k,1)*(v(i,2,k)-v(i,1,k))         &
     &          -ndmpdt*(v(i,1,k)-vbr(i,1,k))
            end do
            end do

!$omp end do

          end if

        end if

      end if

! -----

! Set the north boundary conditions.

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(nbc.ge.4) then

          if(nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

!$omp do schedule(runtime) private(i,k)

            do k=2,nk-2
            do i=1,ni-1
              v(i,nj,k)=v(i,nj,k)+vcpy(i,k,2)*(v(i,njm1,k)-v(i,nj,k))   &
     &          -ndmpdt*(v(i,nj,k)-(vgpv(i,nj,k)+vtd(i,nj,k)*tpdt))
            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(i,k)

            do k=2,nk-2
            do i=1,ni-1
              v(i,nj,k)=v(i,nj,k)+vcpy(i,k,2)*(v(i,njm1,k)-v(i,nj,k))   &
     &          -ndmpdt*(v(i,nj,k)-vbr(i,nj,k))
            end do
            end do

!$omp end do

          end if

        end if

      end if

! -----

! Set the west boundary conditions.

      if(ebw.eq.1.and.isub.eq.0) then

        if(wbc.ge.4) then

          if(nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

!$omp do schedule(runtime) private(j,k)

            do k=2,nk-2
            do j=2,nj-1
              v(1,j,k)=v(1,j,k)-vcpx(j,k,1)*(v(2,j,k)-v(1,j,k))         &
     &          -tdmpdt*(v(1,j,k)-(vgpv(1,j,k)+vtd(1,j,k)*tpdt))
            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(j,k)

            do k=2,nk-2
            do j=2,nj-1
              v(1,j,k)=v(1,j,k)-vcpx(j,k,1)*(v(2,j,k)-v(1,j,k))         &
     &          -tdmpdt*(v(1,j,k)-vbr(1,j,k))
            end do
            end do

!$omp end do

          end if

        end if

      end if

! -----

! Set the east boundary conditions.

      if(ebe.eq.1.and.isub.eq.nisub-1) then

        if(ebc.ge.4) then

          if(nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

!$omp do schedule(runtime) private(j,k)

            do k=2,nk-2
            do j=2,nj-1
              v(nim1,j,k)=v(nim1,j,k)                                   &
     &          +vcpx(j,k,2)*(v(nim2,j,k)-v(nim1,j,k))-tdmpdt           &
     &          *(v(nim1,j,k)-(vgpv(nim1,j,k)+vtd(nim1,j,k)*tpdt))
            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(j,k)

            do k=2,nk-2
            do j=2,nj-1
              v(nim1,j,k)=v(nim1,j,k)                                   &
     &          +vcpx(j,k,2)*(v(nim2,j,k)-v(nim1,j,k))                  &
     &          -tdmpdt*(v(nim1,j,k)-vbr(nim1,j,k))
            end do
            end do

!$omp end do

          end if

        end if

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_rbcv

!-----7--------------------------------------------------------------7--

      end module m_rbcv
