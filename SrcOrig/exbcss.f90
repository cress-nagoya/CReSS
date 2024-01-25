!***********************************************************************
      module m_exbcss
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/06/07
!     Modification: 1999/07/05, 1999/08/03, 1999/08/09, 1999/09/30,
!                   1999/11/01, 2000/01/17, 2000/02/02, 2000/04/18,
!                   2001/01/15, 2001/03/13, 2001/04/15, 2001/05/29,
!                   2001/06/06, 2001/06/29, 2001/07/13, 2001/08/07,
!                   2001/11/20, 2001/12/11, 2002/04/02, 2002/06/06,
!                   2002/07/23, 2002/08/15, 2002/10/31, 2003/03/28,
!                   2003/04/30, 2003/05/19, 2003/06/27, 2003/11/05,
!                   2003/11/28, 2003/12/12, 2004/04/15, 2004/08/20,
!                   2005/01/31, 2005/02/10, 2006/09/21, 2006/12/04,
!                   2007/01/05, 2007/01/31, 2007/05/07, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2008/12/11, 2009/02/27,
!                   2009/03/23, 2011/09/22, 2013/01/28, 2013/02/13,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     force the lateral boundary value to the external boundary value
!     for scalar variables.

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

      public :: exbcss, s_exbcss

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface exbcss

        module procedure s_exbcss

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
      subroutine s_exbcss(fpexbvar,fpwbc,fpebc,fpexnews,ape,            &
     &                    isstp,dts,gtinc,ni,nj,nk,scpx,scpy,sgpv,std,s)
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

      integer, intent(in) :: ape
                       ! Pointer of exbvar

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

      real sb1         ! Temporary variable
      real sb2         ! Temporary variable
      real sb2i        ! Temporary variable
      real sb2j        ! Temporary variable

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

        if(exbvar(ape:ape).eq.'-') then

          if(ebs.eq.1.and.jsub.eq.0) then

            if(ebw.eq.1.and.isub.eq.0) then

!$omp do schedule(runtime) private(k,sb1,sb2i,sb2j,radwe,radsn)

              do k=2,nk-2
                sb1=sgpv(1,1,k)+std(1,1,k)*tpdt
                sb2i=sgpv(2,1,k)+std(2,1,k)*tpdt
                sb2j=sgpv(1,2,k)+std(1,2,k)*tpdt

                radwe=scpx(1,k,1)*((s(2,1,k)-s(1,1,k))-(sb2i-sb1))
                radsn=scpy(1,k,1)*((s(1,2,k)-s(1,1,k))-(sb2j-sb1))

                s(1,1,k)=s(1,1,k)+std(1,1,k)*dts                        &
     &            -(radwe+radsn)-dmpdt*(s(1,1,k)-sb1)

              end do

!$omp end do

            end if

            if(ebe.eq.1.and.isub.eq.nisub-1) then

!$omp do schedule(runtime) private(k,sb1,sb2i,sb2j,radwe,radsn)

              do k=2,nk-2
                sb1=sgpv(nim1,1,k)+std(nim1,1,k)*tpdt
                sb2i=sgpv(nim2,1,k)+std(nim2,1,k)*tpdt
                sb2j=sgpv(nim1,2,k)+std(nim1,2,k)*tpdt

                radwe=scpx(1,k,2)                                       &
     &            *((s(nim2,1,k)-s(nim1,1,k))-(sb2i-sb1))

                radsn=scpy(nim1,k,1)                                    &
     &            *((s(nim1,2,k)-s(nim1,1,k))-(sb2j-sb1))

                s(nim1,1,k)=s(nim1,1,k)+std(nim1,1,k)*dts               &
     &            +(radwe-radsn)-dmpdt*(s(nim1,1,k)-sb1)

              end do

!$omp end do

            end if

          end if

          if(ebn.eq.1.and.jsub.eq.njsub-1) then

            if(ebw.eq.1.and.isub.eq.0) then

!$omp do schedule(runtime) private(k,sb1,sb2i,sb2j,radwe,radsn)

              do k=2,nk-2
                sb1=sgpv(1,njm1,k)+std(1,njm1,k)*tpdt
                sb2i=sgpv(2,njm1,k)+std(2,njm1,k)*tpdt
                sb2j=sgpv(1,njm2,k)+std(1,njm2,k)*tpdt

                radwe=scpx(njm1,k,1)                                    &
     &            *((s(2,njm1,k)-s(1,njm1,k))-(sb2i-sb1))

                radsn=scpy(1,k,2)                                       &
     &            *((s(1,njm2,k)-s(1,njm1,k))-(sb2j-sb1))

                s(1,njm1,k)=s(1,njm1,k)+std(1,njm1,k)*dts               &
     &            -(radwe-radsn)-dmpdt*(s(1,njm1,k)-sb1)

              end do

!$omp end do

            end if

            if(ebe.eq.1.and.isub.eq.nisub-1) then

!$omp do schedule(runtime) private(k,sb1,sb2i,sb2j,radwe,radsn)

              do k=2,nk-2
                sb1=sgpv(nim1,njm1,k)+std(nim1,njm1,k)*tpdt
                sb2i=sgpv(nim2,njm1,k)+std(nim2,njm1,k)*tpdt
                sb2j=sgpv(nim1,njm2,k)+std(nim1,njm2,k)*tpdt

                radwe=scpx(njm1,k,2)                                    &
     &            *((s(nim2,njm1,k)-s(nim1,njm1,k))-(sb2i-sb1))

                radsn=scpy(nim1,k,2)                                    &
     &            *((s(nim1,njm2,k)-s(nim1,njm1,k))-(sb2j-sb1))

                s(nim1,njm1,k)=s(nim1,njm1,k)+std(nim1,njm1,k)*dts      &
     &            +(radwe+radsn)-dmpdt*(s(nim1,njm1,k)-sb1)

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

          if(exbvar(ape:ape).eq.'-') then

!$omp do schedule(runtime) private(j,k,sb1,sb2)

            do k=2,nk-2
            do j=2,nj-2
              sb1=sgpv(1,j,k)+std(1,j,k)*tpdt
              sb2=sgpv(2,j,k)+std(2,j,k)*tpdt

              s(1,j,k)=s(1,j,k)+std(1,j,k)*dts                          &
     &          -scpx(j,k,1)*((s(2,j,k)-s(1,j,k))-(sb2-sb1))            &
     &          -dmpdt*(s(1,j,k)-sb1)

            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(j,k)

            do k=2,nk-2
            do j=2,nj-2
              s(1,j,k)=s(1,j,k)+std(1,j,k)*dts
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

          if(exbvar(ape:ape).eq.'-') then

!$omp do schedule(runtime) private(j,k,sb1,sb2)

            do k=2,nk-2
            do j=2,nj-2
              sb1=sgpv(nim1,j,k)+std(nim1,j,k)*tpdt
              sb2=sgpv(nim2,j,k)+std(nim2,j,k)*tpdt

              s(nim1,j,k)=s(nim1,j,k)+std(nim1,j,k)*dts                 &
     &          +scpx(j,k,2)*((s(nim2,j,k)-s(nim1,j,k))-(sb2-sb1))      &
     &          -dmpdt*(s(nim1,j,k)-sb1)

            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(j,k)

            do k=2,nk-2
            do j=2,nj-2
              s(nim1,j,k)=s(nim1,j,k)+std(nim1,j,k)*dts
            end do
            end do

!$omp end do

          end if

        end if

      end if

! -----

! Force the south boundary value to the external boundary value.

      if(ebs.eq.1.and.jsub.eq.0) then

        if(exbvar(ape:ape).eq.'-') then

!$omp do schedule(runtime) private(i,k,sb1,sb2)

          do k=2,nk-2
          do i=2,ni-2
            sb1=sgpv(i,1,k)+std(i,1,k)*tpdt
            sb2=sgpv(i,2,k)+std(i,2,k)*tpdt

            s(i,1,k)=s(i,1,k)+std(i,1,k)*dts                            &
     &        -scpy(i,k,1)*((s(i,2,k)-s(i,1,k))-(sb2-sb1))              &
     &        -dmpdt*(s(i,1,k)-sb1)

          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,k)

          do k=2,nk-2
          do i=1,ni-1
            s(i,1,k)=s(i,1,k)+std(i,1,k)*dts
          end do
          end do

!$omp end do

        end if

      end if

! -----

! Force the north boundary value to the external boundary value.

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(exbvar(ape:ape).eq.'-') then

!$omp do schedule(runtime) private(i,k,sb1,sb2)

          do k=2,nk-2
          do i=2,ni-2
            sb1=sgpv(i,njm1,k)+std(i,njm1,k)*tpdt
            sb2=sgpv(i,njm2,k)+std(i,njm2,k)*tpdt

            s(i,njm1,k)=s(i,njm1,k)+std(i,njm1,k)*dts                   &
     &        +scpy(i,k,2)*((s(i,njm2,k)-s(i,njm1,k))-(sb2-sb1))        &
     &        -dmpdt*(s(i,njm1,k)-sb1)

          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,k)

          do k=2,nk-2
          do i=1,ni-1
            s(i,njm1,k)=s(i,njm1,k)+std(i,njm1,k)*dts
          end do
          end do

!$omp end do

        end if

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_exbcss

!-----7--------------------------------------------------------------7--

      end module m_exbcss
