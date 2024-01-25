!***********************************************************************
      module m_exbcw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/06/07
!     Modification: 1999/07/05, 1999/08/03, 1999/08/09, 1999/09/30,
!                   1999/11/01, 2000/01/17, 2000/02/02, 2000/04/18,
!                   2001/03/13, 2001/06/06, 2001/06/29, 2001/07/13,
!                   2001/08/07, 2001/12/11, 2002/04/02, 2002/06/06,
!                   2002/07/23, 2002/08/15, 2002/10/31, 2003/03/28,
!                   2003/04/30, 2003/05/19, 2003/06/27, 2003/11/05,
!                   2003/11/28, 2003/12/12, 2004/05/07, 2004/08/20,
!                   2005/01/31, 2005/02/10, 2006/09/21, 2006/12/04,
!                   2007/01/05, 2007/01/31, 2007/05/07, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2008/12/11, 2009/02/27,
!                   2009/03/23, 2011/09/22, 2013/01/28, 2013/02/13,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     force the lateral boundary value to the external boundary value
!     for the z components of velocity.

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

      public :: exbcw, s_exbcw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface exbcw

        module procedure s_exbcw

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
      subroutine s_exbcw(fpexbvar,fpwbc,fpebc,fpexnews,                 &
     &                   isstp,dts,gtinc,ni,nj,nk,wcpx,wcpy,wgpv,wtd,w)
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

      real, intent(in) :: wcpx(1:nj,1:nk,1:2)
                       ! Phase speed of z components of velocity
                       ! on west and east boundary

      real, intent(in) :: wcpy(1:ni,1:nk,1:2)
                       ! Phase speed of z components of velocity
                       ! on south and north boundary

      real, intent(in) :: wgpv(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: wtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! z components of velocity of GPV data

! Input and output variable

      real, intent(inout) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

! Internal shared variables

      character(len=108) exbvar
                       ! Control flag of
                       ! extrenal boundary forced variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions

      integer nim1     ! ni - 1
      integer nim2     ! ni - 2
      integer njm1     ! nj - 1
      integer njm2     ! nj - 2

      real exnews      ! Boundary damping coefficient

      real dmpdt       ! exnews x dts

      real tpdt        ! gtinc + real(isstp - 1) x dts

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real wb1         ! Temporary variable
      real wb2         ! Temporary variable
      real wb2i        ! Temporary variable
      real wb2j        ! Temporary variable

      real radwe       ! Temporary variable
      real radsn       ! Temporary variable

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(exbvar)

! -----

! Get the required namelist variables.

      call getcname(fpexbvar,exbvar)
      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getrname(fpexnews,exnews)

! -----

! Set the common used variables.

      nim1=ni-1
      nim2=ni-2
      njm1=nj-1
      njm2=nj-2

      dmpdt=exnews*dts

      tpdt=gtinc+real(isstp-1)*dts

! -----

!! Force the lateral boundary value to the external boundary value.

!$omp parallel default(shared)

! Force the boundary value to the external boundary value at the four
! corners.

      if(abs(wbc).ne.1.and.abs(ebc).ne.1) then

        if(exbvar(3:3).eq.'-') then

          if(ebs.eq.1.and.jsub.eq.0) then

            if(ebw.eq.1.and.isub.eq.0) then

!$omp do schedule(runtime) private(k,wb1,wb2i,wb2j,radwe,radsn)

              do k=2,nk-1
                wb1=wgpv(1,1,k)+wtd(1,1,k)*tpdt
                wb2i=wgpv(2,1,k)+wtd(2,1,k)*tpdt
                wb2j=wgpv(1,2,k)+wtd(1,2,k)*tpdt

                radwe=wcpx(1,k,1)*((w(2,1,k)-w(1,1,k))-(wb2i-wb1))
                radsn=wcpy(1,k,1)*((w(1,2,k)-w(1,1,k))-(wb2j-wb1))

                w(1,1,k)=w(1,1,k)+wtd(1,1,k)*dts                        &
     &            -(radwe+radsn)-dmpdt*(w(1,1,k)-wb1)

              end do

!$omp end do

            end if

            if(ebe.eq.1.and.isub.eq.nisub-1) then

!$omp do schedule(runtime) private(k,wb1,wb2i,wb2j,radwe,radsn)

              do k=2,nk-1
                wb1=wgpv(nim1,1,k)+wtd(nim1,1,k)*tpdt
                wb2i=wgpv(nim2,1,k)+wtd(nim2,1,k)*tpdt
                wb2j=wgpv(nim1,2,k)+wtd(nim1,2,k)*tpdt

                radwe=wcpx(1,k,2)                                       &
     &            *((w(nim2,1,k)-w(nim1,1,k))-(wb2i-wb1))

                radsn=wcpy(nim1,k,1)                                    &
     &            *((w(nim1,2,k)-w(nim1,1,k))-(wb2j-wb1))

                w(nim1,1,k)=w(nim1,1,k)+wtd(nim1,1,k)*dts               &
     &            +(radwe-radsn)-dmpdt*(w(nim1,1,k)-wb1)

              end do

!$omp end do

            end if

          end if

          if(ebn.eq.1.and.jsub.eq.njsub-1) then

            if(ebw.eq.1.and.isub.eq.0) then

!$omp do schedule(runtime) private(k,wb1,wb2i,wb2j,radwe,radsn)

              do k=2,nk-1
                wb1=wgpv(1,njm1,k)+wtd(1,njm1,k)*tpdt
                wb2i=wgpv(2,njm1,k)+wtd(2,njm1,k)*tpdt
                wb2j=wgpv(1,njm2,k)+wtd(1,njm2,k)*tpdt

                radwe=wcpx(njm1,k,1)                                    &
     &            *((w(2,njm1,k)-w(1,njm1,k))-(wb2i-wb1))

                radsn=wcpy(1,k,2)                                       &
     &            *((w(1,njm2,k)-w(1,njm1,k))-(wb2j-wb1))

                w(1,njm1,k)=w(1,njm1,k)+wtd(1,njm1,k)*dts               &
     &            -(radwe-radsn)-dmpdt*(w(1,njm1,k)-wb1)

              end do

!$omp end do

            end if

            if(ebe.eq.1.and.isub.eq.nisub-1) then

!$omp do schedule(runtime) private(k,wb1,wb2i,wb2j,radwe,radsn)

              do k=2,nk-1
                wb1=wgpv(nim1,njm1,k)+wtd(nim1,njm1,k)*tpdt
                wb2i=wgpv(nim2,njm1,k)+wtd(nim2,njm1,k)*tpdt
                wb2j=wgpv(nim1,njm2,k)+wtd(nim1,njm2,k)*tpdt

                radwe=wcpx(njm1,k,2)                                    &
     &            *((w(nim2,njm1,k)-w(nim1,njm1,k))-(wb2i-wb1))

                radsn=wcpy(nim1,k,2)                                    &
     &            *((w(nim1,njm2,k)-w(nim1,njm1,k))-(wb2j-wb1))

                w(nim1,njm1,k)=w(nim1,njm1,k)+wtd(nim1,njm1,k)*dts      &
     &            +(radwe+radsn)-dmpdt*(w(nim1,njm1,k)-wb1)

              end do

!$omp end do

            end if

          end if

        end if

      end if

! -----

! Force the west boundary value to the external boundary value.

      if(ebw.eq.1.and.isub.eq.0) then

        if(abs(wbc).ne.1) then

          if(exbvar(3:3).eq.'-') then

!$omp do schedule(runtime) private(j,k,wb1,wb2)

            do k=2,nk-1
            do j=2,nj-2
              wb1=wgpv(1,j,k)+wtd(1,j,k)*tpdt
              wb2=wgpv(2,j,k)+wtd(2,j,k)*tpdt

              w(1,j,k)=w(1,j,k)+wtd(1,j,k)*dts                          &
     &          -wcpx(j,k,1)*((w(2,j,k)-w(1,j,k))-(wb2-wb1))            &
     &          -dmpdt*(w(1,j,k)-wb1)

            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(j,k)

            do k=2,nk-1
            do j=2,nj-2
              w(1,j,k)=w(1,j,k)+wtd(1,j,k)*dts
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

          if(exbvar(3:3).eq.'-') then

!$omp do schedule(runtime) private(j,k,wb1,wb2)

            do k=2,nk-1
            do j=2,nj-2
              wb1=wgpv(nim1,j,k)+wtd(nim1,j,k)*tpdt
              wb2=wgpv(nim2,j,k)+wtd(nim2,j,k)*tpdt

              w(nim1,j,k)=w(nim1,j,k)+wtd(nim1,j,k)*dts                 &
     &          +wcpx(j,k,2)*((w(nim2,j,k)-w(nim1,j,k))-(wb2-wb1))      &
     &          -dmpdt*(w(nim1,j,k)-wb1)

            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(j,k)

            do k=2,nk-1
            do j=2,nj-2
              w(nim1,j,k)=w(nim1,j,k)+wtd(nim1,j,k)*dts
            end do
            end do

!$omp end do

          end if

        end if

      end if

! -----

! Force the south boundary value to the external boundary value.

      if(ebs.eq.1.and.jsub.eq.0) then

        if(exbvar(3:3).eq.'-') then

!$omp do schedule(runtime) private(i,k,wb1,wb2)

          do k=2,nk-1
          do i=2,ni-2
            wb1=wgpv(i,1,k)+wtd(i,1,k)*tpdt
            wb2=wgpv(i,2,k)+wtd(i,2,k)*tpdt

            w(i,1,k)=w(i,1,k)+wtd(i,1,k)*dts                            &
     &        -wcpy(i,k,1)*((w(i,2,k)-w(i,1,k))-(wb2-wb1))              &
     &        -dmpdt*(w(i,1,k)-wb1)

          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,k)

          do k=2,nk-1
          do i=1,ni-1
            w(i,1,k)=w(i,1,k)+wtd(i,1,k)*dts
          end do
          end do

!$omp end do

        end if

      end if

! -----

! Force the north boundary value to the external boundary value.

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(exbvar(3:3).eq.'-') then

!$omp do schedule(runtime) private(i,k,wb1,wb2)

          do k=2,nk-1
          do i=2,ni-2
            wb1=wgpv(i,njm1,k)+wtd(i,njm1,k)*tpdt
            wb2=wgpv(i,njm2,k)+wtd(i,njm2,k)*tpdt

            w(i,njm1,k)=w(i,njm1,k)+wtd(i,njm1,k)*dts                   &
     &        +wcpy(i,k,2)*((w(i,njm2,k)-w(i,njm1,k))-(wb2-wb1))        &
     &        -dmpdt*(w(i,njm1,k)-wb1)

          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,k)

          do k=2,nk-1
          do i=1,ni-1
            w(i,njm1,k)=w(i,njm1,k)+wtd(i,njm1,k)*dts
          end do
          end do

!$omp end do

        end if

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_exbcw

!-----7--------------------------------------------------------------7--

      end module m_exbcw
