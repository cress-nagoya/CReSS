!***********************************************************************
      module m_smoo4p
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/09/16
!     Modification: 1999/09/30, 1999/10/12, 1999/11/01, 1999/11/19,
!                   1999/11/24, 2000/01/17, 2000/12/18, 2001/04/15,
!                   2001/06/06, 2001/06/29, 2001/11/20, 2002/04/02,
!                   2002/07/23, 2003/01/04, 2003/04/30, 2003/05/19,
!                   2003/07/15, 2003/11/05, 2003/11/28, 2003/12/12,
!                   2004/01/23, 2004/03/05, 2004/04/15, 2004/08/20,
!                   2004/09/10, 2005/04/04, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform the 4th order smoothing for the pressure.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: smoo4p, s_smoo4p

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface smoo4p

        module procedure s_smoo4p

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_smoo4p(fpsmtopt,fpiwest,fpieast,fpjsouth,fpjnorth,   &
     &                    fpsmhcoe,fpsmvcoe,ni,nj,nk,pp,pfrc,pp2,       &
     &                    tmp1,tmp2,tmp3)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpsmtopt
                       ! Formal parameter of unique index of smtopt

      integer, intent(in) :: fpiwest
                       ! Formal parameter of unique index of iwest

      integer, intent(in) :: fpieast
                       ! Formal parameter of unique index of ieast

      integer, intent(in) :: fpjsouth
                       ! Formal parameter of unique index of jsouth

      integer, intent(in) :: fpjnorth
                       ! Formal parameter of unique index of jnorth

      integer, intent(in) :: fpsmhcoe
                       ! Formal parameter of unique index of smhcoe

      integer, intent(in) :: fpsmvcoe
                       ! Formal parameter of unique index of smvcoe

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation

! Input and output variable

      real, intent(inout) :: pfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in pressure equation

! Internal shared variables

      integer smtopt   ! Option for numerical smoothing

      integer iwest    ! Added index on west boundary
      integer jsouth   ! Added index on south boundary

      integer ieast    ! Subtracted index on east boundary
      integer jnorth   ! Subtracted index on north boundary

      integer nkm1     ! nk - 1
      integer nkm2     ! nk - 2

      real smhcoe      ! Horizontal smoothig coefficient
      real smvcoe      ! Vertical smoothig coefficient

      real, intent(inout) :: pp2(0:ni+1,0:nj+1,1:nk)
                       ! 2.0 x pp

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpsmtopt,smtopt)
      call getiname(fpiwest,iwest)
      call getiname(fpieast,ieast)
      call getiname(fpjsouth,jsouth)
      call getiname(fpjnorth,jnorth)
      call getrname(fpsmhcoe,smhcoe)
      call getrname(fpsmvcoe,smvcoe)

! -----

! Set the common used variables.

      nkm1=nk-1
      nkm2=nk-2

! -----

! Calculate the 4th order pressure numerical smoothing.

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=jsouth,nj-jnorth
        do i=iwest,ni-ieast
          pp2(i,j,k)=2.e0*pp(i,j,k)
        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=1+iwest,ni-1-ieast
          tmp1(i,j,k)=(pp(i+1,j,k)+pp(i-1,j,k))-pp2(i,j,k)
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=1+jsouth,nj-1-jnorth
        do i=2,ni-2
          tmp2(i,j,k)=(pp(i,j+1,k)+pp(i,j-1,k))-pp2(i,j,k)
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          tmp3(i,j,k)=(pp(i,j,k+1)+pp(i,j,k-1))-pp2(i,j,k)
        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          pfrc(i,j,k)=pfrc(i,j,k)                                       &
     &      +smhcoe*(tmp1(i,j,k)+tmp2(i,j,k))+smvcoe*tmp3(i,j,k)
        end do
        end do

!$omp end do

      end do

      if(mod(smtopt,10).eq.2) then

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2+jsouth,nj-2-jnorth
          do i=2+iwest,ni-2-ieast
            pfrc(i,j,k)=pfrc(i,j,k)                                     &
     &        +smhcoe*((tmp1(i,j,k)-(tmp1(i+1,j,k)+tmp1(i-1,j,k)))      &
     &        +(tmp2(i,j,k)-(tmp2(i,j+1,k)+tmp2(i,j-1,k))))
          end do
          end do

!$omp end do

        end do

      else

!$omp do schedule(runtime) private(i,j)

        do j=2+jsouth,nj-2-jnorth
        do i=2+iwest,ni-2-ieast
          tmp3(i,j,1)=tmp3(i,j,2)
          tmp3(i,j,nkm1)=tmp3(i,j,nkm2)
        end do
        end do

!$omp end do

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2+jsouth,nj-2-jnorth
          do i=2+iwest,ni-2-ieast
            pfrc(i,j,k)=pfrc(i,j,k)                                     &
     &        +(smvcoe*(tmp3(i,j,k)-(tmp3(i,j,k+1)+tmp3(i,j,k-1)))      &
     &        +smhcoe*((tmp1(i,j,k)-(tmp1(i+1,j,k)+tmp1(i-1,j,k)))      &
     &        +(tmp2(i,j,k)-(tmp2(i,j+1,k)+tmp2(i,j-1,k)))))
          end do
          end do

!$omp end do

        end do

      end if

!$omp end parallel

! -----

      end subroutine s_smoo4p

!-----7--------------------------------------------------------------7--

      end module m_smoo4p
