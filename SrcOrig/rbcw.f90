!***********************************************************************
      module m_rbcw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/07/28,
!                   1999/08/03, 1999/08/18, 1999/08/23, 1999/09/06,
!                   1999/09/30, 1999/10/07, 1999/11/01, 1999/12/06,
!                   2000/01/17, 2000/03/17, 2000/03/23, 2001/04/15,
!                   2001/07/13, 2001/08/07, 2001/09/13, 2001/12/11,
!                   2002/04/02, 2002/06/06, 2002/07/23, 2002/08/15,
!                   2002/10/31, 2003/03/28, 2003/04/30, 2003/05/19,
!                   2003/06/27, 2003/11/05, 2003/11/28, 2003/12/12,
!                   2004/01/09, 2004/05/07, 2004/08/20, 2005/01/31,
!                   2006/09/21, 2006/12/04, 2007/01/05, 2007/01/31,
!                   2007/03/10, 2007/05/07, 2007/05/21, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2008/12/11, 2009/02/27,
!                   2009/03/23, 2011/09/22, 2013/01/28, 2013/02/13,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the radiative lateral boundary conditions for the z components
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

      public :: rbcw, s_rbcw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rbcw

        module procedure s_rbcw

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
      subroutine s_rbcw(fpgpvvar,fplbcvar,fpwbc,fpebc,fpsbc,fpnbc,      &
     &                  fpnggopt,fplspopt,fpvspopt,fplbnews,isstp,      &
     &                  dts,gtinc,ni,nj,nk,wcpx,wcpy,wgpv,wtd,w)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpgpvvar
                       ! Formal parameter of unique index of gpvvar

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

      character(len=108) gpvvar
                       ! Control flag of input GPV data variables

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
      integer njm2     ! nj - 2

      real lbnews      ! Boundary damping coefficient

      real dmpdt       ! lbnews x dts

      real tpdt        ! gtinc + real(isstp - 1) x dts

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real radwe       ! Temporary variable
      real radsn       ! Temporary variable

!-----7--------------------------------------------------------------7--

! Initialize the character variables.

      call inichar(gpvvar)
      call inichar(lbcvar)

! -----

! Get the required namelist variables.

      call getcname(fpgpvvar,gpvvar)
      call getcname(fplbcvar,lbcvar)
      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpsbc,sbc)
      call getiname(fpnbc,nbc)
      call getiname(fpnggopt,nggopt)
      call getiname(fplspopt,lspopt)
      call getiname(fpvspopt,vspopt)
      call getrname(fplbnews,lbnews)

! -----

! Set the common used variables.

      nim1=ni-1
      nim2=ni-2
      njm1=nj-1
      njm2=nj-2

      if(lbcvar(3:3).eq.'o') then
        dmpdt=lbnews*dts
      else
        dmpdt=0.e0
      end if

      tpdt=gtinc+real(isstp-1)*dts

! -----

!! Set the radiative lateral boundary conditions.

!$omp parallel default(shared)

! Set the boundary conditions at the four corners.

      if(ebs.eq.1.and.jsub.eq.0) then

        if(ebw.eq.1.and.isub.eq.0.and.wbc.ge.4.and.sbc.ge.4) then

          if(gpvvar(1:1).eq.'o'.and.                                    &
     &      (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(k,radwe,radsn)

            do k=2,nk-1
              radwe=(w(2,1,k)-w(1,1,k))*wcpx(1,k,1)
              radsn=(w(1,2,k)-w(1,1,k))*wcpy(1,k,1)

              w(1,1,k)=w(1,1,k)-(radwe+radsn)                           &
     &          -dmpdt*(w(1,1,k)-(wgpv(1,1,k)+wtd(1,1,k)*tpdt))

            end do

!$omp end do

          else

!$omp do schedule(runtime) private(k,radwe,radsn)

            do k=2,nk-1
              radwe=(w(2,1,k)-w(1,1,k))*wcpx(1,k,1)
              radsn=(w(1,2,k)-w(1,1,k))*wcpy(1,k,1)

              w(1,1,k)=w(1,1,k)-(radwe+radsn)-dmpdt*w(1,1,k)

            end do

!$omp end do

          end if

        end if

        if(ebe.eq.1.and.isub.eq.nisub-1.and.ebc.ge.4.and.sbc.ge.4) then

          if(gpvvar(1:1).eq.'o'.and.                                    &
     &      (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(k,radwe,radsn)

            do k=2,nk-1
              radwe=(w(nim2,1,k)-w(nim1,1,k))*wcpx(1,k,2)
              radsn=(w(nim1,2,k)-w(nim1,1,k))*wcpy(nim1,k,1)

              w(nim1,1,k)=w(nim1,1,k)+(radwe-radsn)                     &
     &          -dmpdt*(w(nim1,1,k)-(wgpv(nim1,1,k)+wtd(nim1,1,k)*tpdt))

            end do

!$omp end do

          else

!$omp do schedule(runtime) private(k,radwe,radsn)

            do k=2,nk-1
              radwe=(w(nim2,1,k)-w(nim1,1,k))*wcpx(1,k,2)
              radsn=(w(nim1,2,k)-w(nim1,1,k))*wcpy(nim1,k,1)

              w(nim1,1,k)=w(nim1,1,k)+(radwe-radsn)-dmpdt*w(nim1,1,k)

            end do

!$omp end do

          end if

        end if

      end if

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(ebw.eq.1.and.isub.eq.0.and.wbc.ge.4.and.nbc.ge.4) then

          if(gpvvar(1:1).eq.'o'.and.                                    &
     &      (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(k,radwe,radsn)

            do k=2,nk-1
              radwe=(w(2,njm1,k)-w(1,njm1,k))*wcpx(njm1,k,1)
              radsn=(w(1,njm2,k)-w(1,njm1,k))*wcpy(1,k,2)

              w(1,njm1,k)=w(1,njm1,k)-(radwe-radsn)                     &
     &         -dmpdt*(w(1,njm1,k)-(wgpv(1,njm1,k)+wtd(1,njm1,k)*tpdt))

            end do

!$omp end do

          else

!$omp do schedule(runtime) private(k,radwe,radsn)

            do k=2,nk-1
              radwe=(w(2,njm1,k)-w(1,njm1,k))*wcpx(njm1,k,1)
              radsn=(w(1,njm2,k)-w(1,njm1,k))*wcpy(1,k,2)

              w(1,njm1,k)=w(1,njm1,k)-(radwe-radsn)-dmpdt*w(1,njm1,k)

            end do

!$omp end do

          end if

        end if

        if(ebe.eq.1.and.isub.eq.nisub-1.and.ebc.ge.4.and.nbc.ge.4) then

          if(gpvvar(1:1).eq.'o'.and.                                    &
     &      (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(k,radwe,radsn)

            do k=2,nk-1
              radwe=(w(nim2,njm1,k)-w(nim1,njm1,k))*wcpx(njm1,k,2)
              radsn=(w(nim1,njm2,k)-w(nim1,njm1,k))*wcpy(nim1,k,2)

              w(nim1,njm1,k)=w(nim1,njm1,k)                             &
     &          +(radwe+radsn)-dmpdt*(w(nim1,njm1,k)                    &
     &          -(wgpv(nim1,njm1,k)+wtd(nim1,njm1,k)*tpdt))

            end do

!$omp end do

          else

!$omp do schedule(runtime) private(k,radwe,radsn)

            do k=2,nk-1
              radwe=(w(nim2,njm1,k)-w(nim1,njm1,k))*wcpx(njm1,k,2)
              radsn=(w(nim1,njm2,k)-w(nim1,njm1,k))*wcpy(nim1,k,2)

              w(nim1,njm1,k)=w(nim1,njm1,k)+(radwe+radsn)               &
     &          -dmpdt*w(nim1,njm1,k)

            end do

!$omp end do

          end if

        end if

      end if

! -----

! Set the west boundary conditions.

      if(ebw.eq.1.and.isub.eq.0) then

        if(wbc.ge.4) then

          if(gpvvar(1:1).eq.'o'.and.                                    &
     &      (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(j,k)

            do k=2,nk-1
            do j=2,nj-2
              w(1,j,k)=w(1,j,k)-wcpx(j,k,1)*(w(2,j,k)-w(1,j,k))         &
     &          -dmpdt*(w(1,j,k)-(wgpv(1,j,k)+wtd(1,j,k)*tpdt))
            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(j,k)

            do k=2,nk-1
            do j=2,nj-2
              w(1,j,k)=w(1,j,k)                                         &
     &          -wcpx(j,k,1)*(w(2,j,k)-w(1,j,k))-dmpdt*w(1,j,k)
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

          if(gpvvar(1:1).eq.'o'.and.                                    &
     &      (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(j,k)

            do k=2,nk-1
            do j=2,nj-2
              w(nim1,j,k)=w(nim1,j,k)                                   &
     &          +wcpx(j,k,2)*(w(nim2,j,k)-w(nim1,j,k))                  &
     &          -dmpdt*(w(nim1,j,k)-(wgpv(nim1,j,k)+wtd(nim1,j,k)*tpdt))
            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(j,k)

            do k=2,nk-1
            do j=2,nj-2
              w(nim1,j,k)=w(nim1,j,k)                                   &
     &          +wcpx(j,k,2)*(w(nim2,j,k)-w(nim1,j,k))-dmpdt*w(nim1,j,k)
            end do
            end do

!$omp end do

          end if

        end if

      end if

! -----

! Set the south boundary conditions.

      if(ebs.eq.1.and.jsub.eq.0) then

        if(sbc.ge.4) then

          if(gpvvar(1:1).eq.'o'.and.                                    &
     &      (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(i,k)

            do k=2,nk-1
            do i=2,ni-2
              w(i,1,k)=w(i,1,k)-wcpy(i,k,1)*(w(i,2,k)-w(i,1,k))         &
     &          -dmpdt*(w(i,1,k)-(wgpv(i,1,k)+wtd(i,1,k)*tpdt))
            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(i,k)

            do k=2,nk-1
            do i=2,ni-2
              w(i,1,k)=w(i,1,k)                                         &
     &          -wcpy(i,k,1)*(w(i,2,k)-w(i,1,k))-dmpdt*w(i,1,k)
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

          if(gpvvar(1:1).eq.'o'.and.                                    &
     &      (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(i,k)

            do k=2,nk-1
            do i=2,ni-2
              w(i,njm1,k)=w(i,njm1,k)                                   &
     &          +wcpy(i,k,2)*(w(i,njm2,k)-w(i,njm1,k))                  &
     &          -dmpdt*(w(i,njm1,k)-(wgpv(i,njm1,k)+wtd(i,njm1,k)*tpdt))
            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(i,k)

            do k=2,nk-1
            do i=2,ni-2
              w(i,njm1,k)=w(i,njm1,k)                                   &
     &          +wcpy(i,k,2)*(w(i,njm2,k)-w(i,njm1,k))-dmpdt*w(i,njm1,k)
            end do
            end do

!$omp end do

          end if

        end if

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_rbcw

!-----7--------------------------------------------------------------7--

      end module m_rbcw
