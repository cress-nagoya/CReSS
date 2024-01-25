!***********************************************************************
      module m_lspuvw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/06/07
!     Modification: 1999/07/05, 1999/08/03, 1999/09/06, 1999/09/30,
!                   1999/10/12, 1999/11/01, 2000/01/17, 2000/02/02,
!                   2000/02/07, 2000/04/18, 2001/01/15, 2001/03/13,
!                   2001/04/15, 2001/05/29, 2001/06/06, 2001/06/29,
!                   2001/07/13, 2001/08/07, 2002/04/02, 2002/07/23,
!                   2002/08/15, 2002/12/11, 2003/01/04, 2003/04/30,
!                   2003/05/19, 2003/10/10, 2003/12/12, 2004/05/07,
!                   2004/05/31, 2004/08/01, 2004/08/20, 2006/09/21,
!                   2007/05/07, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2008/12/11, 2009/02/27, 2009/03/23, 2011/09/22,
!                   2013/01/28, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the lateral sponge damping for the velocity.

!-----7--------------------------------------------------------------7--

! Module reference

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

      public :: lspuvw, s_lspuvw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface lspuvw

        module procedure s_lspuvw

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_lspuvw(fpgpvvar,fplspvar,fplspopt,                   &
     &                    fpwdnews,fpwdnorm,fplsnews,fplsnorm,fplspsmt, &
     &                    gtinc,ni,nj,nk,ubr,vbr,rst8u,rst8v,rst8w,     &
     &                    up,vp,wp,rbcx,rbcy,rbcxy,ugpv,utd,vgpv,vtd,   &
     &                    wgpv,wtd,ufrc,vfrc,wfrc,rbc8uv,tmp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpgpvvar
                       ! Formal parameter of unique index of gpvvar

      integer, intent(in) :: fplspvar
                       ! Formal parameter of unique index of lspvar

      integer, intent(in) :: fplspopt
                       ! Formal parameter of unique index of lspopt

      integer, intent(in) :: fpwdnews
                       ! Formal parameter of unique index of wdnews

      integer, intent(in) :: fpwdnorm
                       ! Formal parameter of unique index of wdnorm

      integer, intent(in) :: fplsnews
                       ! Formal parameter of unique index of lsnews

      integer, intent(in) :: fplsnorm
                       ! Formal parameter of unique index of lsnorm

      integer, intent(in) :: fplspsmt
                       ! Formal parameter of unique index of lspsmt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: gtinc
                       ! Lapse of forecast time from GPV data reading

      real, intent(in) :: ubr(0:ni+1,0:nj+1,1:nk)
                       ! Base state x components of velocity

      real, intent(in) :: vbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state y components of velocity

      real, intent(in) :: rst8u(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at u points

      real, intent(in) :: rst8v(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at v points

      real, intent(in) :: rst8w(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at w points

      real, intent(in) :: up(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at past

      real, intent(in) :: vp(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at past

      real, intent(in) :: wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

      real, intent(in) :: rbcx(1:ni)
                       ! Relaxed lateral sponge damping coefficients
                       ! in x direction

      real, intent(in) :: rbcy(1:nj)
                       ! Relaxed lateral sponge damping coefficients
                       ! in y direction

      real, intent(in) :: rbcxy(1:ni,1:nj)
                       ! Relaxed lateral sponge damping coefficients

      real, intent(in) :: ugpv(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: utd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! x components of velocity of GPV data

      real, intent(in) :: vgpv(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: vtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! y components of velocity of GPV data

      real, intent(in) :: wgpv(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: wtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! z components of velocity of GPV data

! Input and output variables

      real, intent(inout) :: ufrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in u equation

      real, intent(inout) :: vfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in v equation

      real, intent(inout) :: wfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in w equation

! Internal shared variables

      character(len=108) gpvvar
                       ! Control flag of input GPV data variables

      character(len=108) lspvar
                       ! Control flag of lateral sponge damped variables

      integer lspopt   ! Option for lateral sponge damping

      integer wdnews   ! Lateral sponge damping thickness

      integer wdnorm   ! Lateral sponge damping thickness
                       ! for u and v in normal

      real lsnews      ! Lateral sponge damping coefficient

      real lsnorm      ! Lateral sponge damping coefficient
                       ! for u and v in normal

      real lspsmt      ! Lateral sponge smoothing coefficient

      real, intent(inout) :: rbc8uv(0:ni+1,0:nj+1)
                       ! Relaxed lateral sponge damping coefficients
                       ! at u or v points

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real a           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Initialize the character variables.

      call inichar(gpvvar)
      call inichar(lspvar)

! -----

! Get the required namelist variables.

      call getcname(fpgpvvar,gpvvar)
      call getcname(fplspvar,lspvar)
      call getiname(fplspopt,lspopt)
      call getiname(fpwdnews,wdnews)
      call getiname(fpwdnorm,wdnorm)
      call getrname(fplsnews,lsnews)
      call getrname(fplsnorm,lsnorm)
      call getrname(fplspsmt,lspsmt)

! -----

!!! Calculate the lateral sponge damping for the velocity.

!$omp parallel default(shared) private(k)

!! For the x components of velocity.

      if((wdnews.ge.1.or.wdnorm.ge.1).and.lspvar(1:1).ne.'x') then

! Damp to the GPV data.

        if(mod(lspopt,10).eq.1) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni
              tmp1(i,j,k)=rst8u(i,j,k)                                  &
     &          *(up(i,j,k)-(ugpv(i,j,k)+utd(i,j,k)*gtinc))
            end do
            end do

!$omp end do

          end do

! -----

! Damp to the base state value.

        else

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni
              tmp1(i,j,k)=rst8u(i,j,k)*(up(i,j,k)-ubr(i,j,k))
            end do
            end do

!$omp end do

          end do

        end if

! -----

! Finally get the lateral sponge damping term.

        if(lspvar(1:1).eq.'o') then

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-1
            rbc8uv(i,j)=.5e0*(rbcxy(i-1,j)+rbcxy(i,j))
          end do
          end do

!$omp end do

          if(lspopt.lt.10) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-1
                ufrc(i,j,k)=ufrc(i,j,k)-lsnews*rbc8uv(i,j)*tmp1(i,j,k)
              end do
              end do

!$omp end do

            end do

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j,a)

              do j=2,nj-2
              do i=2,ni-1
                a=2.e0*tmp1(i,j,k)

                ufrc(i,j,k)=ufrc(i,j,k)-rbc8uv(i,j)*(lsnews*tmp1(i,j,k) &
     &            -lspsmt*(((tmp1(i+1,j,k)+tmp1(i-1,j,k))-a)            &
     &            +((tmp1(i,j+1,k)+tmp1(i,j-1,k))-a)))

              end do
              end do

!$omp end do

            end do

          end if

        else if(lspvar(1:1).eq.'+') then

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-1
            rbc8uv(i,j)=.5e0*(rbcxy(i-1,j)+rbcxy(i,j))
          end do
          end do

!$omp end do

          if(lspopt.lt.10) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-1
                ufrc(i,j,k)=ufrc(i,j,k)                                 &
     &            -(lsnews*rbc8uv(i,j)+lsnorm*rbcx(i))*tmp1(i,j,k)
              end do
              end do

!$omp end do

            end do

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j,a)

              do j=2,nj-2
              do i=2,ni-1
                a=2.e0*tmp1(i,j,k)

                ufrc(i,j,k)=ufrc(i,j,k)                                 &
     &            -(lsnews*rbc8uv(i,j)+lsnorm*rbcx(i))*tmp1(i,j,k)

                ufrc(i,j,k)=ufrc(i,j,k)+lspsmt*(rbc8uv(i,j)+rbcx(i))    &
     &            *(((tmp1(i+1,j,k)+tmp1(i-1,j,k))-a)                   &
     &            +((tmp1(i,j+1,k)+tmp1(i,j-1,k))-a))

              end do
              end do

!$omp end do

            end do

          end if

        else if(lspvar(1:1).eq.'-') then

          if(lspopt.lt.10) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-1
                ufrc(i,j,k)=ufrc(i,j,k)-lsnorm*rbcx(i)*tmp1(i,j,k)
              end do
              end do

!$omp end do

            end do

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j,a)

              do j=2,nj-2
              do i=2,ni-1
                a=2.e0*tmp1(i,j,k)

                ufrc(i,j,k)=ufrc(i,j,k)-rbcx(i)*(lsnorm*tmp1(i,j,k)     &
     &            -lspsmt*(((tmp1(i+1,j,k)+tmp1(i-1,j,k))-a)            &
     &            +((tmp1(i,j+1,k)+tmp1(i,j-1,k))-a)))

              end do
              end do

!$omp end do

            end do

          end if

        end if

! -----

      end if

!! -----

!! For the y components of velocity.

      if((wdnews.ge.1.or.wdnorm.ge.1).and.lspvar(2:2).ne.'x') then

! Damp to the GPV data.

        if(mod(lspopt,10).eq.1) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=1,nj
            do i=1,ni-1
              tmp1(i,j,k)=rst8v(i,j,k)                                  &
     &          *(vp(i,j,k)-(vgpv(i,j,k)+vtd(i,j,k)*gtinc))
            end do
            end do

!$omp end do

          end do

! -----

! Damp to the base state value.

        else

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=1,nj
            do i=1,ni-1
              tmp1(i,j,k)=rst8v(i,j,k)*(vp(i,j,k)-vbr(i,j,k))
            end do
            end do

!$omp end do

          end do

        end if

! -----

! Finally get the lateral sponge damping term.

        if(lspvar(2:2).eq.'o') then

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=2,ni-2
            rbc8uv(i,j)=.5e0*(rbcxy(i,j-1)+rbcxy(i,j))
          end do
          end do

!$omp end do

          if(lspopt.lt.10) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-1
              do i=2,ni-2
                vfrc(i,j,k)=vfrc(i,j,k)-lsnews*rbc8uv(i,j)*tmp1(i,j,k)
              end do
              end do

!$omp end do

            end do

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j,a)

              do j=2,nj-1
              do i=2,ni-2
                a=2.e0*tmp1(i,j,k)

                vfrc(i,j,k)=vfrc(i,j,k)-rbc8uv(i,j)*(lsnews*tmp1(i,j,k) &
     &            -lspsmt*(((tmp1(i+1,j,k)+tmp1(i-1,j,k))-a)            &
     &            +((tmp1(i,j+1,k)+tmp1(i,j-1,k))-a)))

              end do
              end do

!$omp end do

            end do

          end if

        else if(lspvar(2:2).eq.'+') then

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=2,ni-2
            rbc8uv(i,j)=.5e0*(rbcxy(i,j-1)+rbcxy(i,j))
          end do
          end do

!$omp end do

          if(lspopt.lt.10) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-1
              do i=2,ni-2
                vfrc(i,j,k)=vfrc(i,j,k)                                 &
     &            -(lsnews*rbc8uv(i,j)+lsnorm*rbcy(j))*tmp1(i,j,k)
              end do
              end do

!$omp end do

            end do

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j,a)

              do j=2,nj-1
              do i=2,ni-2
                a=2.e0*tmp1(i,j,k)

                vfrc(i,j,k)=vfrc(i,j,k)                                 &
     &            -(lsnews*rbc8uv(i,j)+lsnorm*rbcy(j))*tmp1(i,j,k)

                vfrc(i,j,k)=vfrc(i,j,k)+lspsmt*(rbc8uv(i,j)+rbcy(j))    &
     &            *(((tmp1(i+1,j,k)+tmp1(i-1,j,k))-a)                   &
     &            +((tmp1(i,j+1,k)+tmp1(i,j-1,k))-a))

              end do
              end do

!$omp end do

            end do

          end if

        else if(lspvar(2:2).eq.'-') then

          if(lspopt.lt.10) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-1
              do i=2,ni-2
                vfrc(i,j,k)=vfrc(i,j,k)-lsnorm*rbcy(j)*tmp1(i,j,k)
              end do
              end do

!$omp end do

            end do

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j,a)

              do j=2,nj-1
              do i=2,ni-2
                a=2.e0*tmp1(i,j,k)

                vfrc(i,j,k)=vfrc(i,j,k)-rbcy(j)*(lsnorm*tmp1(i,j,k)     &
     &            -lspsmt*(((tmp1(i+1,j,k)+tmp1(i-1,j,k))-a)            &
     &            +((tmp1(i,j+1,k)+tmp1(i,j-1,k))-a)))

              end do
              end do

!$omp end do

            end do

          end if

        end if

! -----

      end if

!! -----

!! For the z components of velocity.

      if(wdnews.ge.1.and.lspvar(3:3).eq.'o') then

! Damp to the GPV data.

        if(gpvvar(1:1).eq.'o'.and.mod(lspopt,10).eq.1) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              tmp1(i,j,k)=rst8w(i,j,k)                                  &
     &          *(wp(i,j,k)-(wgpv(i,j,k)+wtd(i,j,k)*gtinc))
            end do
            end do

!$omp end do

          end do

! -----

! Damp to the 0.

        else

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              tmp1(i,j,k)=rst8w(i,j,k)*wp(i,j,k)
            end do
            end do

!$omp end do

          end do

        end if

! -----

! Finally get the lateral sponge damping term.

        if(lspopt.lt.10) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-2
              wfrc(i,j,k)=wfrc(i,j,k)-lsnews*rbcxy(i,j)*tmp1(i,j,k)
            end do
            end do

!$omp end do

          end do

        else

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j,a)

            do j=2,nj-2
            do i=2,ni-2
              a=2.e0*tmp1(i,j,k)

              wfrc(i,j,k)=wfrc(i,j,k)-rbcxy(i,j)*(lsnews*tmp1(i,j,k)    &
     &          -lspsmt*(((tmp1(i+1,j,k)+tmp1(i-1,j,k))-a)              &
     &          +((tmp1(i,j+1,k)+tmp1(i,j-1,k))-a)))

            end do
            end do

!$omp end do

          end do

        end if

! -----

      end if

!! -----

!$omp end parallel

!!! -----

      end subroutine s_lspuvw

!-----7--------------------------------------------------------------7--

      end module m_lspuvw
