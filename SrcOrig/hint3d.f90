!***********************************************************************
      module m_hint3d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/03/24
!     Modification: 1999/03/25, 1999/04/06, 1999/05/10, 1999/05/20,
!                   1999/06/28, 1999/07/05, 1999/08/18, 1999/09/30,
!                   1999/10/12, 1999/11/01, 2000/01/17, 2001/04/15,
!                   2001/05/29, 2001/11/20, 2002/04/02, 2002/06/18,
!                   2002/07/03, 2002/07/15, 2003/04/30, 2003/05/19,
!                   2003/10/31, 2003/12/12, 2004/05/07, 2004/08/01,
!                   2004/08/20, 2004/09/10, 2005/02/10, 2006/09/21,
!                   2007/01/20, 2007/05/14, 2007/10/19, 2008/05/02,
!                   2008/06/09, 2008/08/19, 2008/08/25, 2009/01/05,
!                   2009/02/27, 2009/11/13, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     interpolate the data variables to the model grid horizontally.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_commpi
      use m_cpondpe
      use m_destroy
      use m_getiname
      use m_getindx

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: hint3d, s_hint3d

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface hint3d

        module procedure s_hint3d

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic aint
      intrinsic floor
      intrinsic int
      intrinsic max
      intrinsic min
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_hint3d(fpmpopt,fpintopt,xo,fproc,                    &
     &                    ni,nj,nk,ri,rj,di,dj,var,nid,njd,varef)
!***********************************************************************

! Input variables

      character(len=2), intent(in) :: xo
                       ! Control flag of variable arrangement

      character(len=4), intent(in) :: fproc
                       ! Control flag of processing type

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpintopt
                       ! Formal parameter of unique index of intopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nid
                       ! Data dimension in x direction

      integer, intent(in) :: njd
                       ! Data dimension in y direction

      real, intent(in) :: ri(0:ni+1,0:nj+1)
                       ! Real indices in data region in x direction

      real, intent(in) :: rj(0:ni+1,0:nj+1)
                       ! Real indices in data region in y direction

      real, intent(in) :: varef(1:nid,1:njd,1:nk)
                       ! Optional vertically interpolated variable

! Input and output variables

      real, intent(inout) :: di(0:ni+1,0:nj+1)
                       ! Distance between data and model grid points
                       ! in x direction

      real, intent(inout) :: dj(0:ni+1,0:nj+1)
                       ! Distance between data and model grid points
                       ! in y direction

! Output variable

      real, intent(out) :: var(0:ni+1,0:nj+1,1:nk)
                       ! Optional interpolated variable

! Internal shared variables

      integer mpopt    ! Option for map projection
      integer intopt   ! Option for interpolating method

      integer stat     ! Runtime status

      integer istr     ! Minimum do loops index in x direction
      integer iend     ! Maximum do loops index in x direction
      integer jstr     ! Minimum do loops index in y direction
      integer jend     ! Maximum do loops index in y direction

      integer idmin    ! Minimum index of model grid in data region
                       ! in x direction

      integer idmax    ! Maximum index of model grid in data region
                       ! in x direction

      integer jdmin    ! Minimum index of model grid in data region
                       ! in y direction

      integer jdmax    ! Maximum index of model grid in data region
                       ! in y direction

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer id       ! Data array index in x direction
      integer jd       ! Data array index in y direction

      integer idm1     ! id - 1
      integer idp1     ! id + 1
      integer idp2     ! id + 2

      real di1         ! 1.0 - di
      real dj1         ! 1.0 - dj

      real di2         ! di x di
      real di12        ! di1 x di1

      real xint1       ! Temporary variable
      real xint2       ! Temporary variable
      real xint3       ! Temporary variable
      real xint4       ! Temporary variable

      real a           ! Temporary variable
      real b           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpmpopt,mpopt)
      call getiname(fpintopt,intopt)

! -----

! Get the maximum and minimim indices of do loops.

      call getindx(xo,0,ni+1,0,nj+1,istr,iend,jstr,jend)

! -----

!! Calculate the distance between data grid points and model grid
!! points.

      if(fproc(1:3).eq.'cal') then

! Initialize the processed variables.

        idmin=nid
        idmax=1

        jdmin=njd
        jdmax=1

! -----

! Get the distance between data grid points and model grid points.

!$omp parallel default(shared)

        if(mpopt.lt.10) then

!$omp do schedule(runtime) private(i,j,id,jd)                           &
!$omp&   reduction(min: idmin,jdmin) reduction(max: idmax,jdmax)

          do j=jstr,jend
          do i=istr,iend
            id=int(ri(i,j))
            jd=int(rj(i,j))

            idmin=min(id,idmin)
            idmax=max(id,idmax)

            jdmin=min(jd,jdmin)
            jdmax=max(jd,jdmax)

            di(i,j)=ri(i,j)-aint(ri(i,j))
            dj(i,j)=rj(i,j)-aint(rj(i,j))

          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,j,jd)                              &
!$omp&   reduction(min: jdmin) reduction(max: jdmax)

          do j=jstr,jend
          do i=istr,iend
            jd=int(rj(i,j))

            jdmin=min(jd,jdmin)
            jdmax=max(jd,jdmax)

            di(i,j)=ri(i,j)-real(floor(ri(i,j)))
            dj(i,j)=rj(i,j)-aint(rj(i,j))

          end do
          end do

!$omp end do

        end if

!$omp end parallel

! -----

! If error occured, call the procedure destroy.

        stat=0

        if(intopt.eq.1) then

          if(mpopt.lt.10) then

            if(idmin.lt.1.or.idmax+1.gt.nid                             &
     &        .or.jdmin.lt.1.or.jdmax+1.gt.njd) then

              stat=1

            end if

          else

            if(jdmin.lt.1.or.jdmax+1.gt.njd) then

              stat=1

            end if

          end if

        end if

        if(intopt.eq.2) then

          if(mpopt.lt.10) then

            if(idmin-1.lt.1.or.idmax+2.gt.nid                           &
     &        .or.jdmin-1.lt.1.or.jdmax+2.gt.njd) then

              stat=1

            end if

          else

            if(jdmin-1.lt.1.or.jdmax+2.gt.njd) then

              stat=1

            end if

          end if

        end if

        call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('hint3d  ',6,'cont',7,'              ',14,101, &
     &                   stat)

          end if

          call cpondpe

          call destroy('hint3d  ',6,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

! -----

      end if

!! -----

!! Interpolate the data variables to the model grid horizontally.

!$omp parallel default(shared) private(k)

! Interpolate the data variable to the model grid horizontally with the
! linear interpolation.

      if(intopt.eq.1) then

        if(mpopt.lt.10) then

          do k=1,nk

!$omp do schedule(runtime) private(i,j,id,jd,di1,xint1,xint2)

            do j=jstr,jend
            do i=istr,iend
              id=int(ri(i,j))
              jd=int(rj(i,j))

              di1=1.e0-di(i,j)

              xint1=di1*varef(id,jd,k)+di(i,j)*varef(id+1,jd,k)
              xint2=di1*varef(id,jd+1,k)+di(i,j)*varef(id+1,jd+1,k)

              var(i,j,k)=(1.e0-dj(i,j))*xint1+dj(i,j)*xint2

            end do
            end do

!$omp end do

          end do

        else

          do k=1,nk

!$omp do schedule(runtime) private(i,j,id,jd,idp1,di1,xint1,xint2)

            do j=jstr,jend
            do i=istr,iend
              id=floor(ri(i,j))
              jd=int(rj(i,j))

              idp1=id+1

              if(id.lt.1) then
                id=id+nid
              end if

              if(id.gt.nid) then
                id=id-nid
              end if

              if(idp1.lt.1) then
                idp1=idp1+nid
              end if

              if(idp1.gt.nid) then
                idp1=idp1-nid
              end if

              di1=1.e0-di(i,j)

              xint1=di1*varef(id,jd,k)+di(i,j)*varef(idp1,jd,k)
              xint2=di1*varef(id,jd+1,k)+di(i,j)*varef(idp1,jd+1,k)

              var(i,j,k)=(1.e0-dj(i,j))*xint1+dj(i,j)*xint2

            end do
            end do

!$omp end do

          end do

        end if

! -----

! Interpolate the data variable to the model grid horizontally with the
! parabolic interpolation.

      else if(intopt.eq.2) then

        if(mpopt.lt.10) then

          do k=1,nk

!$omp do schedule(runtime) private(i,j,id,jd)                           &
!$omp&   private(di1,dj1,di2,di12,xint1,xint2,xint3,xint4,a,b)

            do j=jstr,jend
            do i=istr,iend
              id=int(ri(i,j))
              jd=int(rj(i,j))

              di1=1.e0-di(i,j)
              dj1=1.e0-dj(i,j)

              di2=di(i,j)*di(i,j)
              di12=di1*di1

              a=(.5e0*(varef(id-1,jd-1,k)+varef(id+1,jd-1,k))           &
     &          -varef(id,jd-1,k))*di2+.5e0*(varef(id+1,jd-1,k)         &
     &          -varef(id-1,jd-1,k))*di(i,j)+varef(id,jd-1,k)

              b=(.5e0*(varef(id,jd-1,k)+varef(id+2,jd-1,k))             &
     &          -varef(id+1,jd-1,k))*di12+.5e0*(varef(id,jd-1,k)        &
     &          -varef(id+2,jd-1,k))*di1+varef(id+1,jd-1,k)

              xint1=di1*a+di(i,j)*b

              a=(.5e0*(varef(id-1,jd,k)+varef(id+1,jd,k))               &
     &          -varef(id,jd,k))*di2+.5e0*(varef(id+1,jd,k)             &
     &          -varef(id-1,jd,k))*di(i,j)+varef(id,jd,k)

              b=(.5e0*(varef(id,jd,k)+varef(id+2,jd,k))                 &
     &          -varef(id+1,jd,k))*di12+.5e0*(varef(id,jd,k)            &
     &          -varef(id+2,jd,k))*di1+varef(id+1,jd,k)

              xint2=di1*a+di(i,j)*b

              a=(.5e0*(varef(id-1,jd+1,k)+varef(id+1,jd+1,k))           &
     &          -varef(id,jd+1,k))*di2+.5e0*(varef(id+1,jd+1,k)         &
     &          -varef(id-1,jd+1,k))*di(i,j)+varef(id,jd+1,k)

              b=(.5e0*(varef(id,jd+1,k)+varef(id+2,jd+1,k))             &
     &          -varef(id+1,jd+1,k))*di12+.5e0*(varef(id,jd+1,k)        &
     &          -varef(id+2,jd+1,k))*di1+varef(id+1,jd+1,k)

              xint3=di1*a+di(i,j)*b

              a=(.5e0*(varef(id-1,jd+2,k)+varef(id+1,jd+2,k))           &
     &          -varef(id,jd+2,k))*di2+.5e0*(varef(id+1,jd+2,k)         &
     &          -varef(id-1,jd+2,k))*di(i,j)+varef(id,jd+2,k)

              b=(.5e0*(varef(id,jd+2,k)+varef(id+2,jd+2,k))             &
     &          -varef(id+1,jd+2,k))*di12+.5e0*(varef(id,jd+2,k)        &
     &          -varef(id+2,jd+2,k))*di1+varef(id+1,jd+2,k)

              xint4=di1*a+di(i,j)*b

              a=(.5e0*(xint1+xint3)-xint2)*dj(i,j)*dj(i,j)              &
     &          +.5e0*(xint3-xint1)*dj(i,j)+xint2

              b=(.5e0*(xint2+xint4)-xint3)*dj1*dj1                      &
     &          +.5e0*(xint2-xint4)*dj1+xint3

              var(i,j,k)=dj1*a+dj(i,j)*b

            end do
            end do

!$omp end do

          end do

        else

          do k=1,nk

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,id,jd,idm1,idp1,idp2)                              &
!$omp&   private(di1,dj1,di2,di12,xint1,xint2,xint3,xint4,a,b)

            do j=jstr,jend
            do i=istr,iend
              id=floor(ri(i,j))
              jd=int(rj(i,j))

              idm1=id-1
              idp1=id+1
              idp2=id+2

              if(id.lt.1) then
                id=id+nid
              end if

              if(id.gt.nid) then
                id=id-nid
              end if

              if(idm1.lt.1) then
                idm1=idm1+nid
              end if

              if(idm1.gt.nid) then
                idm1=idm1-nid
              end if

              if(idp1.lt.1) then
                idp1=idp1+nid
              end if

              if(idp1.gt.nid) then
                idp1=idp1-nid
              end if

              if(idp2.lt.1) then
                idp2=idp2+nid
              end if

              if(idp2.gt.nid) then
                idp2=idp2-nid
              end if

              di1=1.e0-di(i,j)
              dj1=1.e0-dj(i,j)

              di2=di(i,j)*di(i,j)
              di12=di1*di1

              a=(.5e0*(varef(idm1,jd-1,k)+varef(idp1,jd-1,k))           &
     &          -varef(id,jd-1,k))*di2+.5e0*(varef(idp1,jd-1,k)         &
     &          -varef(idm1,jd-1,k))*di(i,j)+varef(id,jd-1,k)

              b=(.5e0*(varef(id,jd-1,k)+varef(idp2,jd-1,k))             &
     &          -varef(idp1,jd-1,k))*di12+.5e0*(varef(id,jd-1,k)        &
     &          -varef(idp2,jd-1,k))*di1+varef(idp1,jd-1,k)

              xint1=di1*a+di(i,j)*b

              a=(.5e0*(varef(idm1,jd,k)+varef(idp1,jd,k))               &
     &          -varef(id,jd,k))*di2+.5e0*(varef(idp1,jd,k)             &
     &          -varef(idm1,jd,k))*di(i,j)+varef(id,jd,k)

              b=(.5e0*(varef(id,jd,k)+varef(idp2,jd,k))                 &
     &          -varef(idp1,jd,k))*di12+.5e0*(varef(id,jd,k)            &
     &          -varef(idp2,jd,k))*di1+varef(idp1,jd,k)

              xint2=di1*a+di(i,j)*b

              a=(.5e0*(varef(idm1,jd+1,k)+varef(idp1,jd+1,k))           &
     &          -varef(id,jd+1,k))*di2+.5e0*(varef(idp1,jd+1,k)         &
     &          -varef(idm1,jd+1,k))*di(i,j)+varef(id,jd+1,k)

              b=(.5e0*(varef(id,jd+1,k)+varef(idp2,jd+1,k))             &
     &          -varef(idp1,jd+1,k))*di12+.5e0*(varef(id,jd+1,k)        &
     &          -varef(idp2,jd+1,k))*di1+varef(idp1,jd+1,k)

              xint3=di1*a+di(i,j)*b

              a=(.5e0*(varef(idm1,jd+2,k)+varef(idp1,jd+2,k))           &
     &          -varef(id,jd+2,k))*di2+.5e0*(varef(idp1,jd+2,k)         &
     &          -varef(idm1,jd+2,k))*di(i,j)+varef(id,jd+2,k)

              b=(.5e0*(varef(id,jd+2,k)+varef(idp2,jd+2,k))             &
     &          -varef(idp1,jd+2,k))*di12+.5e0*(varef(id,jd+2,k)        &
     &          -varef(idp2,jd+2,k))*di1+varef(idp1,jd+2,k)

              xint4=di1*a+di(i,j)*b

              a=(.5e0*(xint1+xint3)-xint2)*dj(i,j)*dj(i,j)              &
     &          +.5e0*(xint3-xint1)*dj(i,j)+xint2

              b=(.5e0*(xint2+xint4)-xint3)*dj1*dj1                      &
     &          +.5e0*(xint2-xint4)*dj1+xint3

              var(i,j,k)=dj1*a+dj(i,j)*b

            end do
            end do

!$omp end do

          end do

        end if

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_hint3d

!-----7--------------------------------------------------------------7--

      end module m_hint3d
