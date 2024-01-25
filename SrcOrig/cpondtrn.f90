!***********************************************************************
      module m_cpondtrn
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/06/07
!     Modification: 1999/06/28, 1999/08/03, 1999/08/09, 1999/09/30,
!                   1999/10/12, 2000/01/17, 2000/02/02, 2000/07/05,
!                   2001/01/15, 2001/04/15, 2001/05/29, 2001/06/06,
!                   2002/04/02, 2002/08/15, 2002/08/27, 2003/04/30,
!                   2003/05/19, 2003/10/10, 2003/11/28, 2004/05/31,
!                   2006/12/04, 2007/01/05, 2007/01/31, 2007/10/19,
!                   2008/05/02, 2008/08/19, 2008/08/25, 2009/02/27,
!                   2009/11/13, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     correspond and damp the model height to the external data height.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commpi
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: cpondtrn, s_cpondtrn

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface cpondtrn

        module procedure s_cpondtrn

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic max
      intrinsic min
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_cpondtrn(fpexbwid,ni,nj,land,htref,ht,dfx,dfy,dfxy)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpexbwid
                       ! Formal parameter of unique index of exbwid

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

      real, intent(in) :: htref(0:ni+1,0:nj+1)
                       ! Horizontally interpolated height

! Input and output variable

      real, intent(inout) :: ht(0:ni+1,0:nj+1)
                       ! Terrain height

! Internal shared variables

      integer exbwid   ! Corresponded thickness
                       ! between model and data height

      integer nx       ! Composite model dimension in x direction
      integer ny       ! Composite model dimension in y direction

      integer igs      ! Start index in group domain in x direction
      integer jgs      ! Start index in group domain in y direction

      integer wd1      ! exbwid+1

      integer wd1x     ! nx - exbwid - 1
      integer wd1y     ! ny - exbwid - 1

      real wdiv        ! 1.0 / real(exbwid)

      real, intent(inout) :: dfx(0:ni+1)
                       ! Interpolating ratio in x direction

      real, intent(inout) :: dfy(0:nj+1)
                       ! Interpolating ratio in y direction

      real, intent(inout) :: dfxy(0:ni+1,0:nj+1)
                       ! Interpolating ratio in x and y direction

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      integer igc      ! Current index in group domain in x direction
      integer jgc      ! Current index in group domain in y direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpexbwid,exbwid)

! -----

! Set the common used variables.

      nx=(ni-3)*nisub+3
      ny=(nj-3)*njsub+3

      igs=(ni-3)*isub
      jgs=(nj-3)*jsub

      wd1=exbwid+1

      wd1x=nx-exbwid-1
      wd1y=ny-exbwid-1

      wdiv=1.e0/real(exbwid)

! -----

! Initialize the processed variables.

      if(exbwid.ge.1) then

        if(ebw.eq.1) then

          if(isub.eq.0) then
            dfx(0)=1.e0
          end if

        end if

        if(ebe.eq.1) then

          if(isub.eq.nisub-1) then
            dfx(ni)=1.e0
          end if

        end if

        if(ebs.eq.1) then

          if(jsub.eq.0) then
            dfy(0)=1.e0
          end if

        end if

        if(ebn.eq.1) then

          if(jsub.eq.njsub-1) then
            dfy(nj)=1.e0
          end if

        end if

        if(ebsw.eq.1) then

          if(isub.eq.0) then
            dfx(0)=1.e0
          end if

          if(jsub.eq.0) then
            dfy(0)=1.e0
          end if

        end if

        if(ebse.eq.1) then

          if(isub.eq.nisub-1) then
            dfx(ni)=1.e0
          end if

          if(jsub.eq.0) then
            dfy(0)=1.e0
          end if

        end if

        if(ebnw.eq.1) then

          if(isub.eq.0) then
            dfx(0)=1.e0
          end if

          if(jsub.eq.njsub-1) then
            dfy(nj)=1.e0
          end if

        end if

        if(ebne.eq.1) then

          if(isub.eq.nisub-1) then
            dfx(ni)=1.e0
          end if

          if(jsub.eq.njsub-1) then
            dfy(nj)=1.e0
          end if

        end if

      end if

! -----

!! Correspond and damp the model height to the external data height.

!$omp parallel default(shared)

      if(exbwid.ge.1) then

! Calculate the interpolating ratio in x direction.

        if(ebw.eq.1) then

!$omp do schedule(runtime) private(i,igc)

          do i=0,ni
            igc=igs+i

            if(igc.le.wd1) then
              dfx(i)=wdiv*real(wd1-igc)
            else
              dfx(i)=0.e0
            end if

          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i)

          do i=0,ni
            dfx(i)=0.e0
          end do

!$omp end do

        end if

        if(ebe.eq.1) then

!$omp do schedule(runtime) private(i,igc)

          do i=0,ni
            igc=igs+i

            if(igc.ge.wd1x) then
              dfx(i)=dfx(i)+wdiv*real(igc-wd1x)
            end if

          end do

!$omp end do

        end if

! -----

! Calculate the interpolating ratio in y direction.

        if(ebs.eq.1) then

!$omp do schedule(runtime) private(j,jgc)

          do j=0,nj
            jgc=jgs+j

            if(jgc.le.wd1) then
              dfy(j)=wdiv*real(wd1-jgc)
            else
              dfy(j)=0.e0
            end if

          end do

!$omp end do

        else

!$omp do schedule(runtime) private(j)

          do j=0,nj
            dfy(j)=0.e0
          end do

!$omp end do

        end if

        if(ebn.eq.1) then

!$omp do schedule(runtime) private(j,jgc)

          do j=0,nj
            jgc=jgs+j

            if(jgc.ge.wd1y) then
              dfy(j)=dfy(j)+wdiv*real(jgc-wd1y)
            end if

          end do

!$omp end do

        end if

! -----

! Get the interpolating ratio in x and y direction.

!$omp do schedule(runtime) private(i,j)

        do j=0,nj
        do i=0,ni
          dfxy(i,j)=max(dfx(i),dfy(j))
        end do
        end do

!$omp end do

! -----

! Calculate the interpolating ratio at the closed corners.

        if(ebsw.eq.1) then

!$omp do schedule(runtime) private(i,igc)

          do i=0,ni
            igc=igs+i

            if(igc.le.wd1) then
              dfx(i)=wdiv*real(wd1-igc)
            else
              dfx(i)=0.e0
            end if

          end do

!$omp end do

!$omp do schedule(runtime) private(j,jgc)

          do j=0,nj
            jgc=jgs+j

            if(jgc.le.wd1) then
              dfy(j)=wdiv*real(wd1-jgc)
            else
              dfy(j)=0.e0
            end if

          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=0,nj
          do i=0,ni
            dfxy(i,j)=dfxy(i,j)+min(dfx(i),dfy(j))
          end do
          end do

!$omp end do

        end if

        if(ebse.eq.1) then

!$omp do schedule(runtime) private(i,igc)

          do i=0,ni
            igc=igs+i

            if(igc.ge.wd1x) then
              dfx(i)=wdiv*real(igc-wd1x)
            else
              dfx(i)=0.e0
            end if

          end do

!$omp end do

!$omp do schedule(runtime) private(j,jgc)

          do j=0,nj
            jgc=jgs+j

            if(jgc.le.wd1) then
              dfy(j)=wdiv*real(wd1-jgc)
            else
              dfy(j)=0.e0
            end if

          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=0,nj
          do i=0,ni
            dfxy(i,j)=dfxy(i,j)+min(dfx(i),dfy(j))
          end do
          end do

!$omp end do

        end if

        if(ebnw.eq.1) then

!$omp do schedule(runtime) private(i,igc)

          do i=0,ni
            igc=igs+i

            if(igc.le.wd1) then
              dfx(i)=wdiv*real(wd1-igc)
            else
              dfx(i)=0.e0
            end if

          end do

!$omp end do

!$omp do schedule(runtime) private(j,jgc)

          do j=0,nj
            jgc=jgs+j

            if(jgc.ge.wd1y) then
              dfy(j)=wdiv*real(jgc-wd1y)
            else
              dfy(j)=0.e0
            end if

          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=0,nj
          do i=0,ni
            dfxy(i,j)=dfxy(i,j)+min(dfx(i),dfy(j))
          end do
          end do

!$omp end do

        end if

        if(ebne.eq.1) then

!$omp do schedule(runtime) private(i,igc)

          do i=0,ni
            igc=igs+i

            if(igc.ge.wd1x) then
              dfx(i)=wdiv*real(igc-wd1x)
            else
              dfx(i)=0.e0
            end if

          end do

!$omp end do

!$omp do schedule(runtime) private(j,jgc)

          do j=0,nj
            jgc=jgs+j

            if(jgc.ge.wd1y) then
              dfy(j)=wdiv*real(jgc-wd1y)
            else
              dfy(j)=0.e0
            end if

          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=0,nj
          do i=0,ni
            dfxy(i,j)=dfxy(i,j)+min(dfx(i),dfy(j))
          end do
          end do

!$omp end do

        end if

! -----

! Finally correspond and damp the model height to the external data
! height.

!$omp do schedule(runtime) private(i,j)

        do j=0,nj
        do i=0,ni

          if(land(i,j).ge.3) then

            ht(i,j)=(1.e0-dfxy(i,j))*ht(i,j)+dfxy(i,j)*htref(i,j)

          end if

        end do
        end do

!$omp end do

! -----

      end if

!$omp end parallel

!! -----

      end subroutine s_cpondtrn

!-----7--------------------------------------------------------------7--

      end module m_cpondtrn
