!***********************************************************************
      module m_outmxn
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/03/23
!     Modification: 2000/04/18, 2000/06/01, 2000/07/05, 2001/01/15,
!                   2001/02/13, 2001/04/15, 2001/05/29, 2001/06/29,
!                   2002/04/02, 2002/07/15, 2002/12/02, 2003/03/28,
!                   2003/04/30, 2003/05/19, 2003/11/28, 2003/12/12,
!                   2004/04/15, 2004/05/31, 2004/06/10, 2004/09/01,
!                   2005/02/10, 2006/01/10, 2006/02/13, 2006/07/21,
!                   2006/09/11, 2006/09/21, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/07/30, 2007/11/26, 2008/05/02,
!                   2008/06/09, 2008/08/25, 2008/10/10, 2009/02/27,
!                   2009/11/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read in the maximum and minimum value of optional prognostic
!     variable to the standard i/o.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkstd
      use m_comkind
      use m_commath
      use m_commpi
      use m_getmxn
      use m_outstd12

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: outmxn, s_outmxn

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outmxn

        module procedure s_outmxn

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic max
      intrinsic min
      intrinsic sign

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_outmxn(vname,ncvn,cunit,ctime,ni,nj,nk,              &
     &                    istr,iend,jstr,jend,kstr,kend,outcnt,var)
!***********************************************************************

! Input variables

      character(len=4), intent(in) :: vname
                       ! Optional variable name

      character(len=7), intent(in) :: cunit
                       ! Optional variable unit

      integer, intent(in) :: ncvn
                       ! Number of character of vname

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: istr
                       ! Minimum do loops index in x direction

      integer, intent(in) :: iend
                       ! Maximum do loops index in x direction

      integer, intent(in) :: jstr
                       ! Minimum do loops index in y direction

      integer, intent(in) :: jend
                       ! Maximum do loops index in y direction

      integer, intent(in) :: kstr
                       ! Minimum do loops index in z direction

      integer, intent(in) :: kend
                       ! Maximum do loops index in z direction

      real, intent(in) :: var(0:ni+1,0:nj+1,1:nk)
                       ! Optional variable

! Input and output variable

      integer, intent(inout) :: outcnt
                       ! Counter of output variables

! Internal shared variables

      integer ies      ! Start index in entire domain in x direction
      integer jes      ! Start index in entire domain in y direction

      integer maxi     ! Index of maximum value point in x direction
      integer maxj     ! Index of maxinum value point in y direction
      integer maxk     ! Index of maxinum value point in z direction

      integer mini     ! Index of mininum value point in x direction
      integer minj     ! Index of mininum value point in y direction
      integer mink     ! Index of mininum value point in z direction

      real maxvl       ! Maximum value
      real minvl       ! Minimum value

      real maxeps      ! maxvl + sign(eps, maxvl)
      real mineps      ! minvl + sign(eps, minvl)

      real chkeps      ! Very small constant
                       ! for real variables comparison

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer ipies    ! i + start index
                       ! in entire domain in x direction
      integer jpjes    ! j + start index
                       ! in entire domain in y direction

      real cvl         ! Temporary variable

!-----7--------------------------------------------------------------7--

!! Calculate and read in the maximum and minimum value of optional
!! prognostic variable.

! Increse the outcnt.

      outcnt=outcnt+1

! -----

! Calculate the start indices in entire domain.

      ies=(ni-3)*(nisub*igrp+isub)
      jes=(nj-3)*(njsub*jgrp+jsub)

! -----

! Initialize the processed variables.

      maxi=0
      maxj=0
      maxk=0

      mini=(ni-3)*nigrp*nisub+4
      minj=(nj-3)*njgrp*njsub+4
      mink=nk+1

      maxvl=lim36n
      minvl=lim36

      maxeps=lim36n
      mineps=lim36

! -----

! Set the common used variable.

      chkeps=1.e-5

! -----

! Get the maximum and minimum value of optional variable.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j,k,cvl)                           &
!$omp&   reduction(max: maxvl,maxeps) reduction(min: minvl,mineps)

      do k=kstr,kend
      do j=jstr,jend
      do i=istr,iend

        cvl=var(i,j,k)+sign(eps,var(i,j,k))

        maxvl=max(var(i,j,k),maxvl)
        minvl=min(var(i,j,k),minvl)

        maxeps=max(cvl,maxeps)
        mineps=min(cvl,mineps)

      end do
      end do
      end do

!$omp end do

!$omp end parallel

! -----

! Reset the maximum and minimum value of optional variable.

      maxeps=1.e0/maxeps
      mineps=1.e0/mineps

! -----

! Get the indices of maximum and minimum value of optional variable.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j,k,ipies,jpjes,cvl)               &
!$omp&   reduction(max: maxi,maxj,maxk) reduction(min: mini,minj,mink)

      do k=kstr,kend
      do j=jstr,jend
      do i=istr,iend

        ipies=i+ies
        jpjes=j+jes

        cvl=var(i,j,k)+sign(eps,var(i,j,k))

        if(abs(cvl*maxeps-1.e0).lt.chkeps) then

          maxi=max(ipies,maxi)
          maxj=max(jpjes,maxj)
          maxk=max(k,maxk)

        end if

        if(abs(cvl*mineps-1.e0).lt.chkeps) then

          mini=min(ipies,mini)
          minj=min(jpjes,minj)
          mink=min(k,mink)

        end if

      end do
      end do
      end do

!$omp end do

!$omp end parallel

! -----

! Gather the maximum or minimum value and indices from the other
! processor elements to the root.

      call getmxn('max',ni,nj,nk,maxi,maxj,maxk,maxvl)
      call getmxn('min',ni,nj,nk,mini,minj,mink,minvl)

! -----

! Read in the message to the standard i/o.

      if(mype.eq.root) then

        call outstd12(2,vname,ncvn,cunit,ctime,                         &
     &                maxi,maxj,maxk,maxvl,mini,minj,mink,minvl)

      end if

      call chkstd(root)

! -----

!! -----

      end subroutine s_outmxn

!-----7--------------------------------------------------------------7--

      end module m_outmxn
