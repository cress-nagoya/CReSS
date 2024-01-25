!***********************************************************************
      module m_lspdmp
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/06/07
!     Modification: 1999/07/05, 1999/08/03, 1999/08/09, 1999/09/30,
!                   1999/11/01, 2000/01/17, 2000/02/02, 2000/07/05,
!                   2001/01/15, 2001/04/15, 2001/05/29, 2001/06/06,
!                   2001/07/13, 2001/08/07, 2002/04/02, 2002/08/15,
!                   2002/10/31, 2002/12/11, 2003/04/30, 2003/05/19,
!                   2003/09/01, 2003/10/10, 2003/12/12, 2004/05/31,
!                   2004/08/01, 2004/08/20, 2005/02/10, 2006/12/04,
!                   2007/01/05, 2007/01/31, 2007/10/10, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2009/08/20,
!                   2009/11/13, 2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the relaxed lateral sponge damping coefficients.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_commpi
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: lspdmp, s_lspdmp

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface lspdmp

        module procedure s_lspdmp

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic cos
      intrinsic max
      intrinsic min
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_lspdmp(fpwbc,fpebc,fpsbc,fpnbc,fpwdnews,fpwdnorm,    &
     &                    ni,nj,rbcx,rbcy,rbcxy)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpsbc
                       ! Formal parameter of unique index of sbc

      integer, intent(in) :: fpnbc
                       ! Formal parameter of unique index of nbc

      integer, intent(in) :: fpwdnews
                       ! Formal parameter of unique index of wdnews

      integer, intent(in) :: fpwdnorm
                       ! Formal parameter of unique index of wdnorm

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in x direction

! Output variables

      real, intent(out) :: rbcx(1:ni)
                       ! Relaxed lateral sponge damping coefficients
                       ! in x direction

      real, intent(out) :: rbcy(1:nj)
                       ! Relaxed lateral sponge damping coefficients
                       ! in y direction

      real, intent(out) :: rbcxy(1:ni,1:nj)
                       ! Relaxed lateral sponge damping coefficients

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions
      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

      integer wdnews   ! Lateral sponge damping thickness

      integer wdnorm   ! Lateral sponge damping thickness
                       ! for u and v in normal

      integer nx       ! Composite model dimension in x direction
      integer ny       ! Composite model dimension in y direction

      integer igs      ! Start index in group domain in x direction
      integer jgs      ! Start index in group domain in y direction

      integer wd2      ! wdnews + 2
      integer wd1x     ! nx - wdnews - 1
      integer wd1y     ! ny - wdnews - 1

      integer wd23     ! 2 x wdnews + 3
      integer wd21x    ! 2 x (nx - wdnews) - 1
      integer wd21y    ! 2 x (ny - wdnews) - 1

      integer wd2n     ! wdnorm + 2
      integer wd1xn    ! nx - wdnorm - 1
      integer wd1yn    ! ny - wdnorm - 1

      integer wd23n    ! 2 x wdnorm + 3
      integer wd21xn   ! 2 x (nx - wdnorm) - 1
      integer wd21yn   ! 2 x (ny - wdnorm) - 1

      real wdiv5       ! 0.5 / real(wdnews)
      real wdiv5n      ! 0.5 / real(wdnorm)

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      integer igc      ! Current index in group domain in x direction
      integer jgc      ! Current index in group domain in y direction

      real rbcx8s      ! Relaxed lateral sponge damping coefficient
                       ! in x direction at scalar points

      real rbcy8s      ! Relaxed lateral sponge damping coefficient
                       ! in y direction at scalar points

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpsbc,sbc)
      call getiname(fpnbc,nbc)
      call getiname(fpwdnews,wdnews)
      call getiname(fpwdnorm,wdnorm)

! -----

! Set the common used variables.

      nx=(ni-3)*nisub+3
      ny=(nj-3)*njsub+3

      igs=(ni-3)*isub
      jgs=(nj-3)*jsub

      wd2=wdnews+2
      wd1x=nx-wdnews-1
      wd1y=ny-wdnews-1

      wd23=2*wdnews+3
      wd21x=2*(nx-wdnews)-1
      wd21y=2*(ny-wdnews)-1

      wd2n=wdnorm+2
      wd1xn=nx-wdnorm-1
      wd1yn=ny-wdnorm-1

      wd23n=2*wdnorm+3
      wd21xn=2*(nx-wdnorm)-1
      wd21yn=2*(ny-wdnorm)-1

      if(wdnews.ge.1) then
        wdiv5=.5e0/real(wdnews)
      else
        wdiv5=0.e0
      end if

      if(wdnorm.ge.1) then
        wdiv5n=.5e0/real(wdnorm)
      else
        wdiv5n=0.e0
      end if

! -----

!!!! Calculate the relaxed lateral sponge damping coefficients.

!!! Calculate the relaxed lateral sponge damping coefficients in x and y
!!! direction.

!$omp parallel default(shared)

!! Set the positive coefficients for damping case.

      if(wdnews.ge.1) then

! Calculate the linear relaxed lateral sponge damping coefficients in x
! direction.

        if(abs(wbc).eq.1.and.abs(ebc).eq.1) then

!$omp do schedule(runtime) private(i)

          do i=1,ni
            rbcx(i)=0.e0
          end do

!$omp end do

        else

          if(ebw.eq.1) then

!$omp do schedule(runtime) private(i,igc)

            do i=1,ni
              igc=igs+i

              if(igc.le.wd2) then
                rbcx(i)=wdiv5*real(wd23-2*igc)
              else
                rbcx(i)=0.e0
              end if

            end do

!$omp end do

          else

!$omp do schedule(runtime) private(i)

            do i=1,ni
              rbcx(i)=0.e0
            end do

!$omp end do

          end if

          if(ebe.eq.1) then

!$omp do schedule(runtime) private(i,igc)

            do i=1,ni
              igc=igs+i

              if(igc.ge.wd1x) then
                rbcx(i)=rbcx(i)+wdiv5*real(2*igc-wd21x)
              end if

            end do

!$omp end do

          end if

        end if

! -----

! Calculate the linear relaxed lateral sponge damping coefficients in y
! direction.

        if(abs(sbc).eq.1.and.abs(nbc).eq.1) then

!$omp do schedule(runtime) private(j)

          do j=1,nj
            rbcy(j)=0.e0
          end do

!$omp end do

        else

          if(ebs.eq.1) then

!$omp do schedule(runtime) private(j,jgc)

            do j=1,nj
              jgc=jgs+j

              if(jgc.le.wd2) then
                rbcy(j)=wdiv5*real(wd23-2*jgc)
              else
                rbcy(j)=0.e0
              end if

            end do

!$omp end do

          else

!$omp do schedule(runtime) private(j)

            do j=1,nj
              rbcy(j)=0.e0
            end do

!$omp end do

          end if

          if(ebn.eq.1) then

!$omp do schedule(runtime) private(j,jgc)

            do j=1,nj
              jgc=jgs+j

              if(jgc.ge.wd1y) then
                rbcy(j)=rbcy(j)+wdiv5*real(2*jgc-wd21y)
              end if

            end do

!$omp end do

          end if

        end if

! -----

! Get the relaxed lateral sponge damping coefficients horizontally.

!$omp do schedule(runtime) private(i,j,rbcx8s,rbcy8s)

        do j=1,nj-1
        do i=1,ni-1
          rbcx8s=.5e0*(rbcx(i)+rbcx(i+1))
          rbcy8s=.5e0*(rbcy(j)+rbcy(j+1))

          rbcxy(i,j)=max(rbcx8s,rbcy8s,0.e0)

        end do
        end do

!$omp end do

! -----

! Calculate the relaxed lateral sponge damping coefficients at the
! closed corners.

        if(ebsw.eq.1) then

!$omp do schedule(runtime) private(i,igc)

          do i=1,ni
            igc=igs+i

            if(igc.le.wd2) then
              rbcx(i)=wdiv5*real(wd23-2*igc)
            else
              rbcx(i)=0.e0
            end if

          end do

!$omp end do

!$omp do schedule(runtime) private(j,jgc)

          do j=1,nj
            jgc=jgs+j

            if(jgc.le.wd2) then
              rbcy(j)=wdiv5*real(wd23-2*jgc)
            else
              rbcy(j)=0.e0
            end if

          end do

!$omp end do

!$omp do schedule(runtime) private(i,j,rbcx8s,rbcy8s)

          do j=1,nj-1
          do i=1,ni-1
            rbcx8s=.5e0*(rbcx(i)+rbcx(i+1))
            rbcy8s=.5e0*(rbcy(j)+rbcy(j+1))

            rbcxy(i,j)=rbcxy(i,j)+min(rbcx8s,rbcy8s)

          end do
          end do

!$omp end do

        end if

        if(ebse.eq.1) then

!$omp do schedule(runtime) private(i,igc)

          do i=1,ni
            igc=igs+i

            if(igc.ge.wd1x) then
              rbcx(i)=wdiv5*real(2*igc-wd21x)
            else
              rbcx(i)=0.e0
            end if

          end do

!$omp end do

!$omp do schedule(runtime) private(j,jgc)

          do j=1,nj
            jgc=jgs+j

            if(jgc.le.wd2) then
              rbcy(j)=wdiv5*real(wd23-2*jgc)
            else
              rbcy(j)=0.e0
            end if

          end do

!$omp end do

!$omp do schedule(runtime) private(i,j,rbcx8s,rbcy8s)

          do j=1,nj-1
          do i=1,ni-1
            rbcx8s=.5e0*(rbcx(i)+rbcx(i+1))
            rbcy8s=.5e0*(rbcy(j)+rbcy(j+1))

            rbcxy(i,j)=rbcxy(i,j)+min(rbcx8s,rbcy8s)

          end do
          end do

!$omp end do

        end if

        if(ebnw.eq.1) then

!$omp do schedule(runtime) private(i,igc)

          do i=1,ni
            igc=igs+i

            if(igc.le.wd2) then
              rbcx(i)=wdiv5*real(wd23-2*igc)
            else
              rbcx(i)=0.e0
            end if

          end do

!$omp end do

!$omp do schedule(runtime) private(j,jgc)

          do j=1,nj
            jgc=jgs+j

            if(jgc.ge.wd1y) then
              rbcy(j)=wdiv5*real(2*jgc-wd21y)
            else
              rbcy(j)=0.e0
            end if

          end do

!$omp end do

!$omp do schedule(runtime) private(i,j,rbcx8s,rbcy8s)

          do j=1,nj-1
          do i=1,ni-1
            rbcx8s=.5e0*(rbcx(i)+rbcx(i+1))
            rbcy8s=.5e0*(rbcy(j)+rbcy(j+1))

            rbcxy(i,j)=rbcxy(i,j)+min(rbcx8s,rbcy8s)

          end do
          end do

!$omp end do

        end if

        if(ebne.eq.1) then

!$omp do schedule(runtime) private(i,igc)

          do i=1,ni
            igc=igs+i

            if(igc.ge.wd1x) then
              rbcx(i)=wdiv5*real(2*igc-wd21x)
            else
              rbcx(i)=0.e0
            end if

          end do

!$omp end do

!$omp do schedule(runtime) private(j,jgc)

          do j=1,nj
            jgc=jgs+j

            if(jgc.ge.wd1y) then
              rbcy(j)=wdiv5*real(2*jgc-wd21y)
            else
              rbcy(j)=0.e0
            end if

          end do

!$omp end do

!$omp do schedule(runtime) private(i,j,rbcx8s,rbcy8s)

          do j=1,nj-1
          do i=1,ni-1
            rbcx8s=.5e0*(rbcx(i)+rbcx(i+1))
            rbcy8s=.5e0*(rbcy(j)+rbcy(j+1))

            rbcxy(i,j)=rbcxy(i,j)+min(rbcx8s,rbcy8s)

          end do
          end do

!$omp end do

        end if

! -----

! Finally apply the cosine function.

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          rbcxy(i,j)=.5e0*(1.e0-cos(cc*max(rbcxy(i,j),0.e0)))
        end do
        end do

!$omp end do

! -----

!! -----

! Set the coefficients with 0 for no damping case.

      else

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          rbcxy(i,j)=0.e0
        end do
        end do

!$omp end do

      end if

! -----

!$omp end parallel

!!! -----

!!! Calculate the relaxed lateral sponge damping coefficients for normal
!!! direction to lateral boundary.

!$omp parallel default(shared)

!! Set the positive coefficients for damping case.

      if(wdnorm.ge.1) then

! Get the relaxed lateral sponge damping coefficients with cosine
! function in x direction.

        if(abs(wbc).eq.1.and.abs(ebc).eq.1) then

!$omp do schedule(runtime) private(i)

          do i=1,ni
            rbcx(i)=0.e0
          end do

!$omp end do

        else

          if(ebw.eq.1) then

!$omp do schedule(runtime) private(i,igc)

            do i=1,ni
              igc=igs+i

              if(igc.le.wd2n) then
                rbcx(i)=wdiv5n*real(wd23n-2*igc)
              else
                rbcx(i)=0.e0
              end if

            end do

!$omp end do

          else

!$omp do schedule(runtime) private(i)

            do i=1,ni
              rbcx(i)=0.e0
            end do

!$omp end do

          end if

          if(ebe.eq.1) then

!$omp do schedule(runtime) private(i,igc)

            do i=1,ni
              igc=igs+i

              if(igc.ge.wd1xn) then
                rbcx(i)=rbcx(i)+wdiv5n*real(2*igc-wd21xn)
              end if

            end do

!$omp end do

          end if

!$omp do schedule(runtime) private(i)

          do i=1,ni
            rbcx(i)=.5e0*(1.e0-cos(cc*max(rbcx(i),0.e0)))
          end do

!$omp end do

        end if

! -----

! Get the relaxed lateral sponge damping coefficients with cosine
! function in y direction.

        if(abs(sbc).eq.1.and.abs(nbc).eq.1) then

!$omp do schedule(runtime) private(j)

          do j=1,nj
            rbcy(j)=0.e0
          end do

!$omp end do

        else

          if(ebs.eq.1) then

!$omp do schedule(runtime) private(j,jgc)

            do j=1,nj
              jgc=jgs+j

              if(jgc.le.wd2n) then
                rbcy(j)=wdiv5n*real(wd23n-2*jgc)
              else
                rbcy(j)=0.e0
              end if

            end do

!$omp end do

          else

!$omp do schedule(runtime) private(j)

            do j=1,nj
              rbcy(j)=0.e0
            end do

!$omp end do

          end if

          if(ebn.eq.1) then

!$omp do schedule(runtime) private(j,jgc)

            do j=1,nj
              jgc=jgs+j

              if(jgc.ge.wd1yn) then
                rbcy(j)=rbcy(j)+wdiv5n*real(2*jgc-wd21yn)
              end if

            end do

!$omp end do

          end if

!$omp do schedule(runtime) private(j)

          do j=1,nj
            rbcy(j)=.5e0*(1.e0-cos(cc*max(rbcy(j),0.e0)))
          end do

!$omp end do

        end if

! -----

!! -----

! Set the coefficients with 0 for no damping case.

      else

!$omp do schedule(runtime) private(i)

        do i=1,ni
          rbcx(i)=0.e0
        end do

!$omp end do

!$omp do schedule(runtime) private(j)

        do j=1,nj
          rbcy(j)=0.e0
        end do

!$omp end do

      end if

! -----

!$omp end parallel

!!! -----

!!!! -----

      end subroutine s_lspdmp

!-----7--------------------------------------------------------------7--

      end module m_lspdmp
