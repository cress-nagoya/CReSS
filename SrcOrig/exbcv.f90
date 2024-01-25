!***********************************************************************
      module m_exbcv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/06/07
!     Modification: 1999/07/05, 1999/08/03, 1999/08/09, 1999/09/06,
!                   1999/09/30, 1999/11/01, 2000/01/17, 2000/02/02,
!                   2000/04/18, 2001/01/15, 2001/03/13, 2001/06/06,
!                   2001/06/29, 2001/07/13, 2001/08/07, 2001/12/11,
!                   2002/04/02, 2002/06/06, 2002/07/23, 2002/08/15,
!                   2002/10/31, 2003/04/30, 2003/05/19, 2003/06/27,
!                   2003/11/05, 2003/11/28, 2003/12/12, 2004/05/07,
!                   2004/08/01, 2004/08/20, 2005/01/31, 2005/02/10,
!                   2006/09/21, 2006/12/04, 2007/01/05, 2007/05/07,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2008/12/11,
!                   2009/02/27, 2009/03/23, 2011/09/22, 2013/01/28,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     force the lateral boundary value to the external boundary value
!     for the y components of velocity.

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

      public :: exbcv, s_exbcv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface exbcv

        module procedure s_exbcv

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_exbcv(fpexbvar,fpwbc,fpebc,fpexnews,fpexnorm,        &
     &                   isstp,dts,gtinc,ni,nj,nk,vcpx,vcpy,vgpv,vtd,v)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpexbvar
                       ! Formal parameter of unique index of exbvar

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpexnews
                       ! Formal parameter of unique index of exnews

      integer, intent(in) :: fpexnorm
                       ! Formal parameter of unique index of exnorm

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

      character(len=108) exbvar
                       ! Control flag of
                       ! extrenal boundary forced variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions

      integer nim1     ! ni - 1
      integer nim2     ! ni - 2
      integer njm1     ! nj - 1

      real exnews      ! Boundary damping coefficient

      real exnorm      ! Boundary damping coefficient
                       ! for u and v in normal

      real tdmpdt      ! exnews x dts
      real ndmpdt      ! exnorm x dts

      real tpdt        ! gtinc + real(isstp - 1) x dts

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real vb1         ! Temporary variable
      real vb2         ! Temporary variable

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(exbvar)

! -----

! Get the required namelist variables.

      call getcname(fpexbvar,exbvar)
      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getrname(fpexnews,exnews)
      call getrname(fpexnorm,exnorm)

! -----

! Set the common used variables.

      nim1=ni-1
      nim2=ni-2
      njm1=nj-1

      tdmpdt=exnews*dts
      ndmpdt=exnorm*dts

      tpdt=gtinc+real(isstp-1)*dts

! -----

!! Force the lateral boundary value to the external boundary value.

!$omp parallel default(shared)

! Force the south boundary value to the external boundary value.

      if(ebs.eq.1.and.jsub.eq.0) then

        if(exbvar(2:2).eq.'-') then

!$omp do schedule(runtime) private(i,k,vb1,vb2)

          do k=2,nk-2
          do i=1,ni-1
            vb1=vgpv(i,1,k)+vtd(i,1,k)*tpdt
            vb2=vgpv(i,2,k)+vtd(i,2,k)*tpdt

            v(i,1,k)=v(i,1,k)+vtd(i,1,k)*dts                            &
     &        -vcpy(i,k,1)*((v(i,2,k)-v(i,1,k))-(vb2-vb1))              &
     &        -ndmpdt*(v(i,1,k)-vb1)

          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,k)

          do k=2,nk-2
          do i=2,ni-2
            v(i,1,k)=v(i,1,k)+vtd(i,1,k)*dts
          end do
          end do

!$omp end do

        end if

      end if

! -----

! Force the north boundary value to the external boundary value.

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(exbvar(2:2).eq.'-') then

!$omp do schedule(runtime) private(i,k,vb1,vb2)

          do k=2,nk-2
          do i=1,ni-1
            vb1=vgpv(i,nj,k)+vtd(i,nj,k)*tpdt
            vb2=vgpv(i,njm1,k)+vtd(i,njm1,k)*tpdt

            v(i,nj,k)=v(i,nj,k)+vtd(i,nj,k)*dts                         &
     &        +vcpy(i,k,2)*((v(i,njm1,k)-v(i,nj,k))-(vb2-vb1))          &
     &        -ndmpdt*(v(i,nj,k)-vb1)

          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,k)

          do k=2,nk-2
          do i=2,ni-2
            v(i,nj,k)=v(i,nj,k)+vtd(i,nj,k)*dts
          end do
          end do

!$omp end do

        end if

      end if

! -----

! Force the west boundary value to the external boundary value.

      if(ebw.eq.1.and.isub.eq.0) then

        if(abs(wbc).ne.1) then

          if(exbvar(2:2).eq.'-'.or.exbvar(2:2).eq.'+') then

!$omp do schedule(runtime) private(j,k,vb1,vb2)

            do k=2,nk-2
            do j=2,nj-1
              vb1=vgpv(1,j,k)+vtd(1,j,k)*tpdt
              vb2=vgpv(2,j,k)+vtd(2,j,k)*tpdt

              v(1,j,k)=v(1,j,k)+vtd(1,j,k)*dts                          &
     &          -vcpx(j,k,1)*((v(2,j,k)-v(1,j,k))-(vb2-vb1))            &
     &          -tdmpdt*(v(1,j,k)-vb1)

            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(j,k)

            do k=2,nk-2
            do j=1,nj
              v(1,j,k)=v(1,j,k)+vtd(1,j,k)*dts
            end do
            end do

!$omp end do

          end if

        end if

      end if

! -----

! Force the east boundary value to the external boundary value.

      if(ebe.eq.1.and.isub.eq.nisub-1) then

        if(abs(ebc).ne.1) then

          if(exbvar(2:2).eq.'-'.or.exbvar(2:2).eq.'+') then

!$omp do schedule(runtime) private(j,k,vb1,vb2)

            do k=2,nk-2
            do j=2,nj-1
              vb1=vgpv(nim1,j,k)+vtd(nim1,j,k)*tpdt
              vb2=vgpv(nim2,j,k)+vtd(nim2,j,k)*tpdt

              v(nim1,j,k)=v(nim1,j,k)+vtd(nim1,j,k)*dts                 &
     &          +vcpx(j,k,2)*((v(nim2,j,k)-v(nim1,j,k))-(vb2-vb1))      &
     &          -tdmpdt*(v(nim1,j,k)-vb1)

            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(j,k)

            do k=2,nk-2
            do j=1,nj
              v(nim1,j,k)=v(nim1,j,k)+vtd(nim1,j,k)*dts
            end do
            end do

!$omp end do

          end if

        end if

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_exbcv

!-----7--------------------------------------------------------------7--

      end module m_exbcv
