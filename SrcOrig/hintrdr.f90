!***********************************************************************
      module m_hintrdr
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/09/09
!     Modification: 2002/12/02, 2003/04/30, 2003/05/19, 2003/12/12,
!                   2004/08/01, 2004/08/20, 2005/02/10, 2006/09/21,
!                   2007/01/20, 2007/10/19, 2008/05/02, 2008/06/09,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     interpolate the data variables to the model grid horizontally.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_getiname
      use m_getindx

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: hintrdr, s_hintrdr

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface hintrdr

        module procedure s_hintrdr

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic aint
      intrinsic int
      intrinsic floor
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_hintrdr(fpmpopt,xo,fproc,ni,nj,nk,ri,rj,di,dj,var,   &
     &                     nid,njd,varef)
!***********************************************************************

! Input variables

      character(len=2), intent(in) :: xo
                       ! Control flag of variable arrangement

      character(len=4), intent(in) :: fproc
                       ! Control flag of processing type

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

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

      integer istr     ! Minimum do loops index in x direction
      integer iend     ! Maximum do loops index in x direction
      integer jstr     ! Minimum do loops index in y direction
      integer jend     ! Maximum do loops index in y direction

      integer nidm1    ! nid - 1
      integer njdm1    ! njd - 1

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer id       ! Data array index in x direction
      integer jd       ! Data array index in y direction

      integer idp1     ! id + 1

      real di1         ! 1.0 - di

      real xint1       ! Temporary variable
      real xint2       ! Temporary variable

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpmpopt,mpopt)

! -----

! Get the maximum and minimim indices of do loops.

      call getindx(xo,1,ni,1,nj,istr,iend,jstr,jend)

! -----

! Set the common used variables.

      nidm1=nid-1
      njdm1=njd-1

! -----

!! Interpolate the data variables to the model grid horizontally.

!$omp parallel default(shared) private(k)

! Calculate the distance between data grid points and model grid
! points.

      if(fproc(1:3).eq.'cal') then

        if(mpopt.lt.10) then

!$omp do schedule(runtime) private(i,j,id,jd)

          do j=jstr,jend
          do i=istr,iend
            id=int(ri(i,j))
            jd=int(rj(i,j))

            if((id.lt.1.or.id.gt.nidm1)                                 &
     &        .or.(jd.lt.1.or.jd.gt.njdm1)) then

              di(i,j)=lim35n
              dj(i,j)=lim35n

            else

              di(i,j)=ri(i,j)-aint(ri(i,j))
              dj(i,j)=rj(i,j)-aint(rj(i,j))

            end if

          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,j,jd)

          do j=jstr,jend
          do i=istr,iend
            jd=int(rj(i,j))

            if(jd.lt.1.or.jd.gt.njdm1) then

              di(i,j)=lim35n
              dj(i,j)=lim35n

            else

              di(i,j)=ri(i,j)-real(floor(ri(i,j)))
              dj(i,j)=rj(i,j)-aint(rj(i,j))

            end if

          end do
          end do

!$omp end do

        end if

      end if

! -----

! Perform interpolation.

      if(mpopt.lt.10) then

        do k=1,nk

!$omp do schedule(runtime) private(i,j,id,jd,di1,xint1,xint2)

          do j=jstr,jend
          do i=istr,iend

            if(di(i,j).gt.lim34n.and.dj(i,j).gt.lim34n) then

              id=int(ri(i,j))
              jd=int(rj(i,j))

              if(varef(id,jd,k).gt.lim34n                               &
     &          .and.varef(id+1,jd,k).gt.lim34n                         &
     &          .and.varef(id,jd+1,k).gt.lim34n                         &
     &          .and.varef(id+1,jd+1,k).gt.lim34n) then

                di1=1.e0-di(i,j)

                xint1=di1*varef(id,jd,k)+di(i,j)*varef(id+1,jd,k)
                xint2=di1*varef(id,jd+1,k)+di(i,j)*varef(id+1,jd+1,k)

                var(i,j,k)=(1.e0-dj(i,j))*xint1+dj(i,j)*xint2

              else

                var(i,j,k)=lim35n

              end if

            else

              var(i,j,k)=lim35n

            end if

          end do
          end do

!$omp end do

        end do

      else

        do k=1,nk

!$omp do schedule(runtime) private(i,j,id,jd,idp1,di1,xint1,xint2)

          do j=jstr,jend
          do i=istr,iend

            if(di(i,j).gt.lim34n.and.dj(i,j).gt.lim34n) then

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

              if(varef(id,jd,k).gt.lim34n                               &
     &          .and.varef(idp1,jd,k).gt.lim34n                         &
     &          .and.varef(id,jd+1,k).gt.lim34n                         &
     &          .and.varef(idp1,jd+1,k).gt.lim34n) then

                di1=1.e0-di(i,j)

                xint1=di1*varef(id,jd,k)+di(i,j)*varef(idp1,jd,k)
                xint2=di1*varef(id,jd+1,k)+di(i,j)*varef(idp1,jd+1,k)

                var(i,j,k)=(1.e0-dj(i,j))*xint1+dj(i,j)*xint2

              else

                var(i,j,k)=lim35n

              end if

            else

              var(i,j,k)=lim35n

            end if

          end do
          end do

!$omp end do

        end do

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_hintrdr

!-----7--------------------------------------------------------------7--

      end module m_hintrdr
