!***********************************************************************
      module m_getmxn
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/03/23
!     Modification: 2001/02/13, 2002/04/02, 2002/10/31, 2003/04/30,
!                   2003/05/19, 2003/11/05, 2004/05/31, 2004/08/20,
!                   2006/12/04, 2007/01/05, 2007/01/20, 2008/05/02,
!                   2008/06/09, 2008/07/25, 2008/08/25, 2008/10/10,
!                   2009/02/27, 2009/11/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the maximum and minimum value of optional prognostic variable
!     in the group domain.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_combuf
      use m_commath
      use m_commpi
      use m_defmpi

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getmxn, s_getmxn

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getmxn

        module procedure s_getmxn

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
      subroutine s_getmxn(fproc,ni,nj,nk,mxni,mxnj,mxnk,mxnvl)
!***********************************************************************

! Input variables

      character(len=3), intent(in) :: fproc
                       ! Control flag of processing type

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

! Input and output variables

      integer, intent(inout) :: mxni
                       ! Index of maximum or minimum value point
                       ! in x direction

      integer, intent(inout) :: mxnj
                       ! Index of maximum or minimum value point
                       ! in y direction

      integer, intent(inout) :: mxnk
                       ! Index of maximum or minimum value point
                       ! in z direction

      real, intent(inout) :: mxnvl
                       ! Maximum or minimum value

! Internal shared variables

      integer maxi     ! Index of maximum value point in x direction
      integer maxj     ! Index of maxinum value point in y direction
      integer maxk     ! Index of maxinum value point in z direction

      integer mini     ! Index of mininum value point in x direction
      integer minj     ! Index of mininum value point in y direction
      integer mink     ! Index of mininum value point in z direction

      integer ierr     ! Error descriptor

      integer bufijk(1:3)
                       ! Sending buffer of mxni, mxnj and mxnk

      real maxvl       ! Maximum value
      real minvl       ! Minimum value

      real maxeps      ! maxvl + sign(eps, maxvl)
      real mineps      ! minvl + sign(eps, minvl)

      real chkeps      ! Very small constant
                       ! for real variables comparison

! Internal private variables

      integer ijpe     ! Common used array index

      real cvl         ! Temporary variable

!-----7--------------------------------------------------------------7--

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

! Gather the maximum or minimum value and their indices from the other
! processor elements to the root.

      bufijk(1)=mxni
      bufijk(2)=mxnj
      bufijk(3)=mxnk

      call mpi_gather(mxnvl,1,mpi_real,mxnbuf,1,mpi_real,root,          &
     &                mpi_comm_cress,ierr)

      call mpi_gather(bufijk,3,mpi_integer,idxbuf,3,mpi_integer,root,   &
     &                mpi_comm_cress,ierr)

! -----

!!! Get the maximum or minimum value and their indices.

      if(npe.ge.2) then

!! Get the maximum value and its indices.

        if(fproc(1:3).eq.'max') then

! Get the maximum value.

!$omp parallel default(shared)

!$omp do schedule(runtime)                                              &
!$omp&   private(ijpe,cvl) reduction(max: maxvl,maxeps)

          do ijpe=0,npe-1

            cvl=mxnbuf(ijpe)+sign(eps,mxnbuf(ijpe))

            maxvl=max(mxnbuf(ijpe),maxvl)

            maxeps=max(cvl,maxeps)

          end do

!$omp end do

!$omp end parallel

! -----

! Reset the maximum value.

          maxeps=1.e0/maxeps

! -----

! Get the indices of maximum value.

!$omp parallel default(shared)

!$omp do schedule(runtime)                                              &
!$omp&   private(ijpe,cvl) reduction(max: maxi,maxj,maxk)

          do ijpe=0,npe-1

            cvl=mxnbuf(ijpe)+sign(eps,mxnbuf(ijpe))

            if(abs(cvl*maxeps-1.e0).lt.chkeps) then

              maxi=max(idxbuf(1,ijpe),maxi)
              maxj=max(idxbuf(2,ijpe),maxj)
              maxk=max(idxbuf(3,ijpe),maxk)

            end if

          end do

!$omp end do

!$omp end parallel

! -----

!! -----

!! Get the minimum value and its indices.

        else if(fproc(1:3).eq.'min') then

! Get the minimum value.

!$omp parallel default(shared)

!$omp do schedule(runtime)                                              &
!$omp&   private(ijpe,cvl) reduction(min: minvl,mineps)

          do ijpe=0,npe-1

            cvl=mxnbuf(ijpe)+sign(eps,mxnbuf(ijpe))

            minvl=min(mxnbuf(ijpe),minvl)

            mineps=min(cvl,mineps)

          end do

!$omp end do

!$omp end parallel

! -----

! Reset the minimum value.

          mineps=1.e0/mineps

! -----

! Get the indices of minimum value.

!$omp parallel default(shared)

!$omp do schedule(runtime)                                              &
!$omp&   private(ijpe,cvl) reduction(min: mini,minj,mink)

          do ijpe=0,npe-1

            cvl=mxnbuf(ijpe)+sign(eps,mxnbuf(ijpe))

            if(abs(cvl*mineps-1.e0).lt.chkeps) then

              mini=min(idxbuf(1,ijpe),mini)
              minj=min(idxbuf(2,ijpe),minj)
              mink=min(idxbuf(3,ijpe),mink)

            end if

          end do

!$omp end do

!$omp end parallel

! -----

        end if

!! -----

      end if

!!! -----

! Finally set the output variables.

      if(npe.ge.2) then

        if(fproc(1:3).eq.'max') then

          mxnvl=maxvl

          mxni=maxi
          mxnj=maxj
          mxnk=maxk

        else if(fproc(1:3).eq.'min') then

          mxnvl=minvl

          mxni=mini
          mxnj=minj
          mxnk=mink

        end if

      end if

! -----

      end subroutine s_getmxn

!-----7--------------------------------------------------------------7--

      end module m_getmxn
