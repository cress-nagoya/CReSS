!***********************************************************************
      module m_chkmxn
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/04/10
!     Modification: 2004/04/15, 2004/05/31, 2004/06/10, 2005/02/10,
!                   2006/01/10, 2006/09/21, 2007/01/05, 2007/01/20,
!                   2007/07/30, 2008/05/02, 2008/06/09, 2008/08/25,
!                   2008/10/10, 2009/01/05, 2009/02/27, 2009/11/05

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     check the maximum and minimum value of optional data.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
      use m_commath
      use m_destroy
      use m_outstd12

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: chkmxn, s_chkmxn

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface chkmxn

        module procedure s_chkmxn

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
      subroutine s_chkmxn(vname,ncvn,cunit,fproc,ctime,stat,            &
     &                    limmin,limmax,undmin,undmax,                  &
     &                    nid,njd,nkd,vardat)
!***********************************************************************

! Input variables

      character(len=4), intent(in) :: vname
                       ! Optional variable name

      character(len=7), intent(in) :: cunit
                       ! Optional variable unit

      character(len=5), intent(in) :: fproc
                       ! Control flag of processing type

      integer, intent(in) :: ncvn
                       ! Number of character of vname

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

      integer, intent(in) :: nid
                       ! Data dimension in x direction

      integer, intent(in) :: njd
                       ! Data dimension in y direction

      integer, intent(in) :: nkd
                       ! Data dimension in z direction

      real, intent(in) :: limmin
                       ! Minimum limit value

      real, intent(in) :: limmax
                       ! Maximum limit value

      real, intent(in) :: undmin
                       ! Minimum undefined limit value

      real, intent(in) :: undmax
                       ! Maximum undefined limit value

      real, intent(in) :: vardat(1:nid,1:njd,1:nkd)
                       ! Optional variable in data

! Output variable

      integer, intent(out) :: stat
                       ! Runtime status

! Internal shared variables

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

      integer id       ! Array index in x direction
      integer jd       ! Array index in y direction
      integer kd       ! Array index in z direction

      real cvl         ! Temporary variable

!-----7--------------------------------------------------------------7--

! Initialize the processed variables.

      maxi=0
      maxj=0
      maxk=0

      mini=nid+1
      minj=njd+1
      mink=nkd+1

      maxvl=lim36n
      minvl=lim36

      maxeps=lim36n
      mineps=lim36

! -----

! Set the common used variable.

      chkeps=1.e-5

! -----

!! Get the maximum and minimum value of optional data and their indices.

! Get the maximum and minimum value of optional data.

!$omp parallel default(shared)

      if(fproc(1:3).eq.'all') then

!$omp do schedule(runtime) private(id,jd,kd,cvl)                        &
!$omp&   reduction(max: maxvl,maxeps) reduction(min: minvl,mineps)

        do kd=1,nkd
        do jd=1,njd
        do id=1,nid

          cvl=vardat(id,jd,kd)+sign(eps,vardat(id,jd,kd))

          maxvl=max(vardat(id,jd,kd),maxvl)
          minvl=min(vardat(id,jd,kd),minvl)

          maxeps=max(cvl,maxeps)
          mineps=min(cvl,mineps)

        end do
        end do
        end do

!$omp end do

      else

!$omp do schedule(runtime) private(id,jd,kd,cvl)                        &
!$omp&   reduction(max: maxvl,maxeps) reduction(min: minvl,mineps)

        do kd=1,nkd
        do jd=1,njd
        do id=1,nid

          if(vardat(id,jd,kd).ge.undmin                                 &
     &      .and.vardat(id,jd,kd).le.undmax) then

            cvl=vardat(id,jd,kd)+sign(eps,vardat(id,jd,kd))

            maxvl=max(vardat(id,jd,kd),maxvl)
            minvl=min(vardat(id,jd,kd),minvl)

            maxeps=max(cvl,maxeps)
            mineps=min(cvl,mineps)

          end if

        end do
        end do
        end do

!$omp end do

      end if

!$omp end parallel

! -----

! Reset the maximum and minimum value of optional data.

      maxeps=1.e0/maxeps
      mineps=1.e0/mineps

! -----

! Get the indices of maximum and minimum value of optional data.

!$omp parallel default(shared)

      if(fproc(1:3).eq.'all') then

!$omp do schedule(runtime) private(id,jd,kd,cvl)                        &
!$omp&   reduction(max: maxi,maxj,maxk) reduction(min: mini,minj,mink)

        do kd=1,nkd
        do jd=1,njd
        do id=1,nid

          cvl=vardat(id,jd,kd)+sign(eps,vardat(id,jd,kd))

          if(abs(cvl*maxeps-1.e0).lt.chkeps) then

            maxi=max(id,maxi)
            maxj=max(jd,maxj)
            maxk=max(kd,maxk)

          end if

          if(abs(cvl*mineps-1.e0).lt.chkeps) then

            mini=min(id,mini)
            minj=min(jd,minj)
            mink=min(kd,mink)

          end if

        end do
        end do
        end do

!$omp end do

      else

!$omp do schedule(runtime) private(id,jd,kd,cvl)                        &
!$omp&   reduction(max: maxi,maxj,maxk) reduction(min: mini,minj,mink)

        do kd=1,nkd
        do jd=1,njd
        do id=1,nid

          if(vardat(id,jd,kd).ge.undmin                                 &
     &      .and.vardat(id,jd,kd).le.undmax) then

            cvl=vardat(id,jd,kd)+sign(eps,vardat(id,jd,kd))

            if(abs(cvl*maxeps-1.e0).lt.chkeps) then

              maxi=max(id,maxi)
              maxj=max(jd,maxj)
              maxk=max(kd,maxk)

            end if

            if(abs(cvl*mineps-1.e0).lt.chkeps) then

              mini=min(id,mini)
              minj=min(jd,minj)
              mink=min(kd,mink)

            end if

          end if

        end do
        end do
        end do

!$omp end do

      end if

!$omp end parallel

! -----

!! -----

! Read in the maximum and minimum value of optional data to the standard
! i/o.

      call outstd12(2,vname,ncvn,cunit,ctime,                           &
     &              maxi,maxj,maxk,maxvl,mini,minj,mink,minvl)

! -----

! Check the maximum and minimum value of optional data.

      if(maxvl.gt.limmax.or.minvl.lt.limmin) then

        stat=1

        call destroy('chkmxn  ',6,'cont',7,'              ',14,101,stat)

      else

        stat=0

      end if

! -----

      end subroutine s_chkmxn

!-----7--------------------------------------------------------------7--

      end module m_chkmxn
