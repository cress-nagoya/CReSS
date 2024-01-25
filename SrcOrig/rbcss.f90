!***********************************************************************
      module m_rbcss
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/07/28,
!                   1999/08/03, 1999/08/18, 1999/08/23, 1999/09/01,
!                   1999/09/06, 1999/09/30, 1999/10/07, 1999/11/01,
!                   1999/12/06, 2000/01/17, 2000/03/17, 2001/03/23,
!                   2001/04/15, 2001/07/13, 2001/08/07, 2001/11/20,
!                   2001/12/11, 2002/04/02, 2002/06/06, 2002/07/23,
!                   2002/08/15, 2002/10/31, 2003/03/28, 2003/04/30,
!                   2003/05/19, 2003/06/27, 2003/11/05, 2003/11/28,
!                   2003/12/12, 2004/01/09, 2004/04/15, 2004/05/07,
!                   2004/08/20, 2005/01/31, 2006/09/21, 2006/12/04,
!                   2007/01/05, 2007/01/31, 2007/03/10, 2007/05/07,
!                   2007/05/21, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2008/12/11, 2009/02/27, 2009/03/23, 2011/09/22,
!                   2013/01/28, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the radiative lateral boundary conditions for optional scalar
!     variable.

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

      public :: rbcss, s_rbcss

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rbcss

        module procedure s_rbcss

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
      subroutine s_rbcss(fplbcvar,fpwbc,fpebc,fpsbc,fpnbc,fpnggopt,     &
     &                   fplspopt,fpvspopt,fplbnews,apl,isstp,dts,      &
     &                   gtinc,ni,nj,nk,scpx,scpy,sgpv,std,s)
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

      integer, intent(in) :: apl
                       ! Pointer of lbcvar

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

      real, intent(in) :: scpx(1:nj,1:nk,1:2)
                       ! Phase speed of optional scalar variable
                       ! on west and east boundary

      real, intent(in) :: scpy(1:ni,1:nk,1:2)
                       ! Phase speed of optional scalar variable
                       ! on south and north boundary

      real, intent(in) :: sgpv(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable of GPV data
                       ! at marked time

      real, intent(in) :: std(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! optional scalar variable of GPV data

! Input and output variable

      real, intent(inout) :: s(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable

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

! -----

! Set the common used variables.

      nim1=ni-1
      nim2=ni-2
      njm1=nj-1
      njm2=nj-2

      if(lbcvar(apl:apl).eq.'o') then
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

          if(nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

!$omp do schedule(runtime) private(k,radwe,radsn)

            do k=2,nk-2
              radwe=(s(2,1,k)-s(1,1,k))*scpx(1,k,1)
              radsn=(s(1,2,k)-s(1,1,k))*scpy(1,k,1)

              s(1,1,k)=s(1,1,k)-(radwe+radsn)                           &
     &          -dmpdt*(s(1,1,k)-(sgpv(1,1,k)+std(1,1,k)*tpdt))

            end do

!$omp end do

          else

!$omp do schedule(runtime) private(k,radwe,radsn)

            do k=2,nk-2
              radwe=(s(2,1,k)-s(1,1,k))*scpx(1,k,1)
              radsn=(s(1,2,k)-s(1,1,k))*scpy(1,k,1)

              s(1,1,k)=s(1,1,k)-(radwe+radsn)-dmpdt*s(1,1,k)

            end do

!$omp end do

          end if

        end if

        if(ebe.eq.1.and.isub.eq.nisub-1.and.ebc.ge.4.and.sbc.ge.4) then

          if(nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

!$omp do schedule(runtime) private(k,radwe,radsn)

            do k=2,nk-2
              radwe=(s(nim2,1,k)-s(nim1,1,k))*scpx(1,k,2)
              radsn=(s(nim1,2,k)-s(nim1,1,k))*scpy(nim1,k,1)

              s(nim1,1,k)=s(nim1,1,k)+(radwe-radsn)                     &
     &          -dmpdt*(s(nim1,1,k)-(sgpv(nim1,1,k)+std(nim1,1,k)*tpdt))

            end do

!$omp end do

          else

!$omp do schedule(runtime) private(k,radwe,radsn)

            do k=2,nk-2
              radwe=(s(nim2,1,k)-s(nim1,1,k))*scpx(1,k,2)
              radsn=(s(nim1,2,k)-s(nim1,1,k))*scpy(nim1,k,1)

              s(nim1,1,k)=s(nim1,1,k)+(radwe-radsn)-dmpdt*s(nim1,1,k)

            end do

!$omp end do

          end if

        end if

      end if

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(ebw.eq.1.and.isub.eq.0.and.wbc.ge.4.and.nbc.ge.4) then

          if(nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

!$omp do schedule(runtime) private(k,radwe,radsn)

            do k=2,nk-2
              radwe=(s(2,njm1,k)-s(1,njm1,k))*scpx(njm1,k,1)
              radsn=(s(1,njm2,k)-s(1,njm1,k))*scpy(1,k,2)

              s(1,njm1,k)=s(1,njm1,k)-(radwe-radsn)                     &
     &          -dmpdt*(s(1,njm1,k)-(sgpv(1,njm1,k)+std(1,njm1,k)*tpdt))

            end do

!$omp end do

          else

!$omp do schedule(runtime) private(k,radwe,radsn)

            do k=2,nk-2
              radwe=(s(2,njm1,k)-s(1,njm1,k))*scpx(njm1,k,1)
              radsn=(s(1,njm2,k)-s(1,njm1,k))*scpy(1,k,2)

              s(1,njm1,k)=s(1,njm1,k)-(radwe-radsn)-dmpdt*s(1,njm1,k)

            end do

!$omp end do

          end if

        end if

        if(ebe.eq.1.and.isub.eq.nisub-1.and.ebc.ge.4.and.nbc.ge.4) then

          if(nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

!$omp do schedule(runtime) private(k,radwe,radsn)

            do k=2,nk-2
              radwe=(s(nim2,njm1,k)-s(nim1,njm1,k))*scpx(njm1,k,2)
              radsn=(s(nim1,njm2,k)-s(nim1,njm1,k))*scpy(nim1,k,2)

              s(nim1,njm1,k)=s(nim1,njm1,k)                             &
     &          +(radwe+radsn)-dmpdt*(s(nim1,njm1,k)                    &
     &          -(sgpv(nim1,njm1,k)+std(nim1,njm1,k)*tpdt))

            end do

!$omp end do

          else

!$omp do schedule(runtime) private(k,radwe,radsn)

            do k=2,nk-2
              radwe=(s(nim2,njm1,k)-s(nim1,njm1,k))*scpx(njm1,k,2)
              radsn=(s(nim1,njm2,k)-s(nim1,njm1,k))*scpy(nim1,k,2)

              s(nim1,njm1,k)=s(nim1,njm1,k)+(radwe+radsn)               &
     &          -dmpdt*s(nim1,njm1,k)

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
            do j=2,nj-2
              s(1,j,k)=s(1,j,k)-scpx(j,k,1)*(s(2,j,k)-s(1,j,k))         &
     &          -dmpdt*(s(1,j,k)-(sgpv(1,j,k)+std(1,j,k)*tpdt))
            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(j,k)

            do k=2,nk-2
            do j=2,nj-2
              s(1,j,k)=s(1,j,k)                                         &
     &          -scpx(j,k,1)*(s(2,j,k)-s(1,j,k))-dmpdt*s(1,j,k)
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
            do j=2,nj-2
              s(nim1,j,k)=s(nim1,j,k)                                   &
     &          +scpx(j,k,2)*(s(nim2,j,k)-s(nim1,j,k))                  &
     &          -dmpdt*(s(nim1,j,k)-(sgpv(nim1,j,k)+std(nim1,j,k)*tpdt))
            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(j,k)

            do k=2,nk-2
            do j=2,nj-2
              s(nim1,j,k)=s(nim1,j,k)                                   &
     &          +scpx(j,k,2)*(s(nim2,j,k)-s(nim1,j,k))-dmpdt*s(nim1,j,k)
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

          if(nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

!$omp do schedule(runtime) private(i,k)

            do k=2,nk-2
            do i=2,ni-2
              s(i,1,k)=s(i,1,k)-scpy(i,k,1)*(s(i,2,k)-s(i,1,k))         &
     &          -dmpdt*(s(i,1,k)-(sgpv(i,1,k)+std(i,1,k)*tpdt))
            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(i,k)

            do k=2,nk-2
            do i=2,ni-2
              s(i,1,k)=s(i,1,k)                                         &
     &          -scpy(i,k,1)*(s(i,2,k)-s(i,1,k))-dmpdt*s(i,1,k)
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
            do i=2,ni-2
              s(i,njm1,k)=s(i,njm1,k)                                   &
     &          +scpy(i,k,2)*(s(i,njm2,k)-s(i,njm1,k))                  &
     &          -dmpdt*(s(i,njm1,k)-(sgpv(i,njm1,k)+std(i,njm1,k)*tpdt))
            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(i,k)

            do k=2,nk-2
            do i=2,ni-2
              s(i,njm1,k)=s(i,njm1,k)                                   &
     &          +scpy(i,k,2)*(s(i,njm2,k)-s(i,njm1,k))-dmpdt*s(i,njm1,k)
            end do
            end do

!$omp end do

          end if

        end if

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_rbcss

!-----7--------------------------------------------------------------7--

      end module m_rbcss
