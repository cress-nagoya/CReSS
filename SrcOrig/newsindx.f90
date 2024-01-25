!***********************************************************************
      module m_newsindx
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/07/01
!     Modification: 2004/08/01, 2004/09/10, 2005/02/10, 2006/09/30,
!                   2007/01/20, 2007/01/31, 2007/10/19, 2008/05/02,
!                   2008/06/09, 2008/08/25, 2008/12/11, 2009/01/05,
!                   2009/02/27, 2009/11/13, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the mininum and maximum data indices to create the base state
!     variables.

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

      public :: newsindx, s_newsindx

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface newsindx

        module procedure s_newsindx

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic floor
      intrinsic int
      intrinsic max
      intrinsic min

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_newsindx(fpmpopt,fpmpopt_gpv,                        &
     &                      nid,njd,idstr,idend,jdstr,jdend,ni,nj,ri,rj)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmpopt_gpv
                       ! Formal parameter of unique index of mpopt_gpv

      integer, intent(in) :: nid
                       ! Data dimension in x direction

      integer, intent(in) :: njd
                       ! Data dimension in y direction

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      real, intent(in) :: ri(0:ni+1,0:nj+1)
                       ! Real indices in data region in x direction

      real, intent(in) :: rj(0:ni+1,0:nj+1)
                       ! Real indices in data region in y direction

! Input and output variables

      integer, intent(inout) :: idstr
                       ! Minimum index of model grid in data region
                       ! in x direction

      integer, intent(inout) :: idend
                       ! Maximum index of model grid in data region
                       ! in x direction

      integer, intent(inout) :: jdstr
                       ! Minimum index of model grid in data region
                       ! in y direction

      integer, intent(inout) :: jdend
                       ! Maximum index of model grid in data region
                       ! in y direction

! Internal shared variables

      integer mpopt    ! Option for map projection

      integer mpopt_gpv
                       ! Option for map projection for GPV data

      integer stat     ! Runtime status

      integer cidstr   ! Minimum index of model grid in data region
                       ! in x direction of current process

      integer cidend   ! Maximum index of model grid in data region
                       ! in x direction of current process

      integer cjdstr   ! Minimum index of model grid in data region
                       ! in y direction of current process

      integer cjdend   ! Maximum index of model grid in data region
                       ! in y direction of current process

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      integer id       ! Data array index in x direction
      integer jd       ! Data array index in y direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpmpopt,mpopt)
      call getiname(fpmpopt_gpv,mpopt_gpv)

! -----

! Initialize the processed variables.

      cidstr=2*nid
      cidend=-nid

      cjdstr=njd
      cjdend=1

! -----

! Get the mininum and maximum data indices to create the base state
! variables.

!$omp parallel default(shared)

      if(mpopt.lt.10) then

!$omp do schedule(runtime) private(i,j,id,jd)                           &
!$omp&   reduction(min: cidstr,cjdstr) reduction(max: cidend,cjdend)

        do j=0,nj
        do i=0,ni
          id=floor(ri(i,j))
          jd=int(rj(i,j))

          cidstr=min(id,cidstr)
          cidend=max(id,cidend)

          cjdstr=min(jd,cjdstr)
          cjdend=max(jd,cjdend)

        end do
        end do

!$omp end do

      else

!$omp do schedule(runtime) private(i,j,jd)                              &
!$omp&   reduction(min: cjdstr) reduction(max: cjdend)

        do j=0,nj
        do i=0,ni
          jd=int(rj(i,j))

          cjdstr=min(jd,cjdstr)
          cjdend=max(jd,cjdend)

        end do
        end do

!$omp end do

      end if

!$omp end parallel

! -----

! Finally reset the mininum and maximum data indices to create the base
! state variables.

      if(mpopt.lt.10) then

        cidend=cidend+1
        cjdend=cjdend+1

        if(mpopt_gpv.ge.10) then

          if(cidstr.lt.1) then
            cidstr=cidstr+nid
          end if

          if(cidstr.gt.nid) then
            cidstr=cidstr-nid
          end if

          if(cidend.lt.1) then
            cidend=cidend+nid
          end if

          if(cidend.gt.nid) then
            cidend=cidend-nid
          end if

        end if

        idstr=min(cidstr,idstr)
        idend=max(cidend,idend)

        jdstr=min(cjdstr,jdstr)
        jdend=max(cjdend,jdend)

      else

        idstr=1
        idend=nid

        jdstr=min(cjdstr,jdstr)
        jdend=max(cjdend+1,jdend)

      end if

! -----

! If error occured, call the procedure destroy.

      stat=0

      if(mpopt_gpv.lt.10) then

        if(idstr.lt.1.or.idend.gt.nid                                   &
     &    .or.jdstr.lt.1.or.jdend.gt.njd) then

          stat=1

        end if

      else

        if(jdstr.lt.1.or.jdend.gt.njd) then

          stat=1

        end if

      end if

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('newsindx',8,'cont',7,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('newsindx',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

      end subroutine s_newsindx

!-----7--------------------------------------------------------------7--

      end module m_newsindx
