!***********************************************************************
      module m_gettrn
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/03/25, 1999/04/06, 1999/06/07, 1999/06/21,
!                   1999/08/18, 1999/08/23, 1999/09/01, 1999/09/06,
!                   1999/09/30, 1999/10/12, 2000/01/17, 2000/04/18,
!                   2000/12/18, 2001/01/15, 2001/05/29, 2001/11/20,
!                   2002/04/02, 2002/06/18, 2002/07/15, 2002/10/31,
!                   2003/04/30, 2003/05/19, 2003/10/31, 2005/08/05,
!                   2005/12/13, 2006/09/21, 2007/01/05, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2008/10/10, 2009/02/27,
!                   2009/11/13, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the terrain height.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_getiname
      use m_getrname
      use m_rdtrn

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: gettrn, s_gettrn

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface gettrn

        module procedure s_gettrn

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic max

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_gettrn(fptrnopt,fpzsfc,fpmnthgh,fpmntwx,fpmntwy,     &
     &                    fpmntcx,fpmntcy,dvname,ncdvn,fmsg,            &
     &                    ni,nj,xs,ys,ht)
!***********************************************************************

! Input variables

      character(len=12), intent(in) :: dvname
                       ! Optional data variable name

      integer, intent(in) :: fptrnopt
                       ! Formal parameter of unique index of trnopt

      integer, intent(in) :: fpzsfc
                       ! Formal parameter of unique index of zsfc

      integer, intent(in) :: fpmnthgh
                       ! Formal parameter of unique index of mnthgh

      integer, intent(in) :: fpmntwx
                       ! Formal parameter of unique index of mntwx

      integer, intent(in) :: fpmntwy
                       ! Formal parameter of unique index of mntwy

      integer, intent(in) :: fpmntcx
                       ! Formal parameter of unique index of mntcx

      integer, intent(in) :: fpmntcy
                       ! Formal parameter of unique index of mntcy

      integer, intent(in) :: ncdvn
                       ! Number of character of dvname

      integer, intent(in) :: fmsg
                       ! Control flag of message type in outstd03

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      real, intent(in) :: xs(0:ni+1)
                       ! x coordinates at scalar points

      real, intent(in) :: ys(0:nj+1)
                       ! y coordinates at scalar points

! Output variable

      real, intent(out) :: ht(0:ni+1,0:nj+1)
                       ! Terrain height

! Internal shared variables

      integer trnopt   ! Option for terrain height setting

      real zsfc        ! Sea surface terrain height

      real mnthgh(1:2) ! Flat or bell shaped mountain height
                       ! and base level height

      real mntwx       ! Bell shaped mountain width in x direction
      real mntwy       ! Bell shaped mountain width in y direction

      real mntcx       ! Center in x coordinates of
                       ! bell shaped mountain

      real mntcy       ! Center in y coordinates of
                       ! bell shaped mountain

      real wxiv        ! 1.0 / mntwx
      real wyiv        ! 1.0 / mntwy

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      real a           ! Temporary variable
      real b           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fptrnopt,trnopt)

! -----

!! Fill in the array ht with the surface height in the case the flat
!! terrain is applied.

      if(trnopt.eq.0) then

! Get the required namelist variables.

        call getrname(fpzsfc,zsfc)
        call getrname(fpmnthgh,mnthgh(1))
        call getrname(fpmnthgh+1,mnthgh(2))

! -----

! Set the flat terrain.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j)

        do j=0,nj
        do i=0,ni
          ht(i,j)=max(mnthgh(1)+mnthgh(2),zsfc)
        end do
        end do

!$omp end do

!$omp end parallel

! -----

!! -----

!! Set the bell shaped mountain.

      else if(trnopt.eq.1) then

! Get the required namelist variables.

        call getrname(fpzsfc,zsfc)
        call getrname(fpmnthgh,mnthgh(1))
        call getrname(fpmnthgh+1,mnthgh(2))
        call getrname(fpmntwx,mntwx)
        call getrname(fpmntwy,mntwy)
        call getrname(fpmntcx,mntcx)
        call getrname(fpmntcy,mntcy)

! -----

! Set the common used variables.

        wxiv=1.e0/mntwx
        wyiv=1.e0/mntwy

! -----

! Set the bell shaped mountain.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j,a,b)

        do j=0,nj
        do i=0,ni
          a=wxiv*(xs(i)-mntcx)
          b=wyiv*(ys(j)-mntcy)

          ht(i,j)=max(mnthgh(1)/(1.e0+(a*a+b*b))+mnthgh(2),zsfc)

        end do
        end do

!$omp end do

!$omp end parallel

! -----

!! -----

! Read out the data from the terrain file.

      else if(trnopt.eq.2) then

        call rdtrn(idexprim,idcrsdir,idncexp,idnccrs,idwlngth,          &
     &             dvname,ncdvn,fmsg,ni,nj,ht)

      end if

! -----

      end subroutine s_gettrn

!-----7--------------------------------------------------------------7--

      end module m_gettrn
