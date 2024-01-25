!***********************************************************************
      module m_undefsst
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/07/26
!     Modification: 2002/09/09, 2002/11/11, 2003/04/30, 2003/05/19,
!                   2003/07/15, 2003/11/28, 2004/04/10, 2004/08/01,
!                   2004/09/01, 2004/09/10, 2004/10/12, 2005/02/10,
!                   2007/01/20, 2007/01/31, 2007/05/14, 2007/10/19,
!                   2008/05/02, 2008/06/09, 2008/08/25, 2009/01/05,
!                   2009/02/27, 2009/11/13, 2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     convert the undefined value to the interpolated value in the sea
!     surface temperature data.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_destroy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: undefsst, s_undefsst

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface undefsst

        module procedure s_undefsst

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic aint
      intrinsic max
      intrinsic min

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_undefsst(stat,nid,njd,sstdat,und)
!***********************************************************************

! Input variables

      integer, intent(in) :: nid
                       ! Data dimension in x direction

      integer, intent(in) :: njd
                       ! Data dimension in y direction

! Input and output variable

      real, intent(inout) :: sstdat(1:nid,1:njd)
                       ! Sea suface temperature in data

! Output variable

      integer, intent(out) :: stat
                       ! Runtime status

! Internal shared variables

      integer nidm1    ! nid - 1
      integer njdm1    ! njd - 1

      real rstat       ! Real runtime status

      real t268v       ! Inverse of lowest sea surface temperature

      real minsst      ! Minimum value of sea surface temperature data
      real maxsst      ! Maximum value of sea surface temperature data

      real, intent(inout) :: und(1:nid,1:njd)
                       ! Control flag to specify undefined value point

! Internal private variables

      integer id       ! Data array index in x direction
      integer jd       ! Data array index in y direction

      real weiw        ! Weighting parameter
      real weie        ! Weighting parameter
      real weis        ! Weighting parameter
      real wein        ! Weighting parameter

      real weiws       ! Weighting parameter
      real weiwn       ! Weighting parameter
      real weies       ! Weighting parameter
      real weien       ! Weighting parameter

      real sumwei      ! weiw + weie + weis + wein
                       ! + weiws + weiwn + weies + weien

!-----7--------------------------------------------------------------7--

! Set the common used variables.

      nidm1=nid-1
      njdm1=njd-1

      t268v=1.e0/268.16e0

! -----

! Initialize the processed variable, rstat.

      rstat=0.e0

! -----

! Check the undefined value.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(id,jd) reduction(+: rstat)

      do jd=1,njd
      do id=1,nid

        if(sstdat(id,jd).ge.268.16e0.and.sstdat(id,jd).le.323.16e0) then
          rstat=rstat+1.e0
        end if

      end do
      end do

!$omp end do

!$omp end parallel

! -----

! If error occured, call the procedure destroy.

      if(rstat.gt.0.e0) then
        stat=0
      else
        stat=1
      end if

      if(stat.ne.0) then

        call destroy('undefsst',8,'cont',7,'              ',14,101,stat)

        return

      end if

! -----

!!! Convert the undefined value to the interpolated value.

      iterate: do

! Initialize the processed variables.

        minsst=lim36
        maxsst=lim36n

! -----

!! Check and convert the undefined value.

!$omp parallel default(shared)

! Check the undefined value.

!$omp do schedule(runtime) private(id,jd)

        do jd=1,njd
        do id=1,nid

          if(sstdat(id,jd).lt.268.16e0                                  &
     &      .or.sstdat(id,jd).gt.323.16e0) then

            und(id,jd)=.1e0

          else

            und(id,jd)=sstdat(id,jd)

          end if

        end do
        end do

!$omp end do

! -----

! Convert the undefined value.

!$omp do schedule(runtime) private(id,jd)                               &
!$omp&   private(weiw,weie,weis,wein,weiws,weiwn,weies,weien,sumwei)

        do jd=2,njd-1
        do id=2,nid-1

          if(und(id,jd).lt.1.e0) then

            if(und(id-1,jd).gt.1.e0.or.und(id+1,jd).gt.1.e0.or.         &
     &         und(id,jd-1).gt.1.e0.or.und(id,jd+1).gt.1.e0.or.         &
     &         und(id-1,jd-1).gt.1.e0.or.und(id-1,jd+1).gt.1.e0.or.     &
     &         und(id+1,jd-1).gt.1.e0.or.und(id+1,jd+1).gt.1.e0) then

              weiw=aint(und(id-1,jd)*t268v+.1e0)
              weie=aint(und(id+1,jd)*t268v+.1e0)
              weis=aint(und(id,jd-1)*t268v+.1e0)
              wein=aint(und(id,jd+1)*t268v+.1e0)

              weiws=aint(und(id-1,jd-1)*t268v+.1e0)
              weiwn=aint(und(id-1,jd+1)*t268v+.1e0)
              weies=aint(und(id+1,jd-1)*t268v+.1e0)
              weien=aint(und(id+1,jd+1)*t268v+.1e0)

              sumwei=weiw+weie+weis+wein+weiws+weiwn+weies+weien

              sstdat(id,jd)=weiw*und(id-1,jd)+weie*und(id+1,jd)         &
     &          +weis*und(id,jd-1)+wein*und(id,jd+1)

              sstdat(id,jd)=(sstdat(id,jd)                              &
     &          +weiws*und(id-1,jd-1)+weiwn*und(id-1,jd+1)              &
     &          +weies*und(id+1,jd-1)+weien*und(id+1,jd+1))/sumwei

            end if

          end if

        end do
        end do

!$omp end do

! -----

! Check the convergence of the iteration.

!$omp do schedule(runtime) private(id,jd)                               &
!$omp&   reduction(min: minsst) reduction(max: maxsst)

        do jd=2,njd-1
        do id=2,nid-1
          minsst=min(sstdat(id,jd),minsst)
          maxsst=max(sstdat(id,jd),maxsst)
        end do
        end do

!$omp end do

! -----

!$omp end parallel

!! -----

! Exit the iteration.

        if(minsst.ge.268.16e0.and.maxsst.le.323.16e0) then

          exit iterate

        end if

! -----

      end do iterate

!!! -----

! Set the boundary conditions.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(jd)

      do jd=1,njd
        sstdat(1,jd)=sstdat(2,jd)
        sstdat(nid,jd)=sstdat(nidm1,jd)
      end do

!$omp end do

!$omp do schedule(runtime) private(id)

      do id=1,nid
        sstdat(id,1)=sstdat(id,2)
        sstdat(id,njd)=sstdat(id,njdm1)
      end do

!$omp end do

!$omp end parallel

! -----

      end subroutine s_undefsst

!-----7--------------------------------------------------------------7--

      end module m_undefsst
