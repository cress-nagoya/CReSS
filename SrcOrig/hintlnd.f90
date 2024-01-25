!***********************************************************************
      module m_hintlnd
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/07/03
!     Modification: 2003/04/30, 2003/05/19, 2004/09/10, 2005/02/10,
!                   2007/01/20, 2007/05/14, 2007/10/19, 2008/05/02,
!                   2008/08/19, 2008/08/25, 2009/01/05, 2009/02/27,
!                   2009/11/13, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     interpolate the land use data to the model grid.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_commpi
      use m_cpondpe
      use m_destroy
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: hintlnd, s_hintlnd

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface hintlnd

        module procedure s_hintlnd

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic int
      intrinsic max
      intrinsic min
      intrinsic nint

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_hintlnd(fpmpopt,ni,nj,ri,rj,land,nid,njd,landat)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nid
                       ! Data dimension in x direction

      integer, intent(in) :: njd
                       ! Data dimension in y direction

      integer, intent(in) :: landat(1:nid,1:njd)
                       ! Land use in data

      real, intent(in) :: ri(0:ni+1,0:nj+1)
                       ! Real indices in data region in x direction

      real, intent(in) :: rj(0:ni+1,0:nj+1)
                       ! Real indices in data region in y direction

! Output variable

      integer, intent(out) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

! Internal shared variables

      integer mpopt    ! Option for map projection

      integer stat     ! Runtime status

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

      integer id       ! Data array index in x direction
      integer jd       ! Data array index in y direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpmpopt,mpopt)

! -----

! Initialize the processed variables.

      idmin=nid
      idmax=1

      jdmin=njd
      jdmax=1

! -----

! Get the required indices at the four corners in data grid.

!$omp parallel default(shared)

      if(mpopt.lt.10) then

!$omp do schedule(runtime) private(i,j,id,jd)                           &
!$omp&   reduction(min: idmin,jdmin) reduction(max: idmax,jdmax)

        do j=0,nj
        do i=0,ni
          id=int(ri(i,j))
          jd=int(rj(i,j))

          idmin=min(id,idmin)
          idmax=max(id,idmax)

          jdmin=min(jd,jdmin)
          jdmax=max(jd,jdmax)

        end do
        end do

!$omp end do

      else

!$omp do schedule(runtime) private(i,j,jd)                              &
!$omp&   reduction(min: jdmin) reduction(max: jdmax)

        do j=0,nj
        do i=0,ni
          jd=int(rj(i,j))

          jdmin=min(jd,jdmin)
          jdmax=max(jd,jdmax)

        end do
        end do

!$omp end do

      end if

!$omp end parallel

! -----

! If error occured, call the procedure destroy.

      stat=0

      if(mpopt.lt.10) then

        if(idmin.lt.1.or.idmax+1.gt.nid                                 &
     &    .or.jdmin.lt.1.or.jdmax+1.gt.njd) then

          stat=1

        end if

      else

        if(jdmin.lt.1.or.jdmax+1.gt.njd) then

          stat=1

        end if

      end if

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('hintlnd ',7,'cont',7,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('hintlnd ',7,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

! Interpolate the land use data to the model grid.

!$omp parallel default(shared)

      if(mpopt.lt.10) then

!$omp do schedule(runtime) private(i,j,id,jd)

        do j=0,nj
        do i=0,ni
          id=nint(ri(i,j))
          jd=nint(rj(i,j))

          land(i,j)=landat(id,jd)

        end do
        end do

!$omp end do

      else

!$omp do schedule(runtime) private(i,j,id,jd)

        do j=0,nj
        do i=0,ni
          id=nint(ri(i,j))
          jd=nint(rj(i,j))

          if(id.lt.1) then
            id=id+nid
          end if

          if(id.gt.nid) then
            id=id-nid
          end if

          land(i,j)=landat(id,jd)

        end do
        end do

!$omp end do

      end if

!$omp end parallel

! -----

      end subroutine s_hintlnd

!-----7--------------------------------------------------------------7--

      end module m_hintlnd
