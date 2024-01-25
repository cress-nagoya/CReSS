!***********************************************************************
      module m_estimsfc
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/07/03
!     Modification: 2003/04/30, 2003/05/19, 2004/05/07, 2004/09/10,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2008/10/10,
!                   2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     re-estimate the surface value on the water, ice and snow surface.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: estimsfc, s_estimsfc

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface estimsfc

        module procedure s_estimsfc

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

!     none

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_estimsfc(fpnumctg_lnd,lnduse_lnd,sfcvar_lnd,         &
     &                      ni,nj,land,sfcvar)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpnumctg_lnd
                       ! Formal parameter of unique index of numctg_lnd

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: land(0:ni+1,0:nj+1)
                       ! Land use on surface

      integer, intent(in) :: lnduse_lnd(1:100)
                       ! User specified land use namelist table

      real, intent(in) :: sfcvar_lnd(1:100)
                       ! User specified surface namelist table

! Output variable

      real, intent(out) :: sfcvar(0:ni+1,0:nj+1)
                       ! Optional surface variable

! Internal shared variable

      integer numctg_lnd
                       ! Number of land use data categories

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      integer ic       ! Index of user specified land use namelist table

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpnumctg_lnd,numctg_lnd)

! -----

! Reestimate the surface value on the water, ice and snow surface.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j,ic)

      do ic=1,numctg_lnd

        do j=1,nj-1
        do i=1,ni-1

          if(land(i,j).eq.lnduse_lnd(ic)) then

            if(land(i,j).lt.10) then

              sfcvar(i,j)=sfcvar_lnd(ic)

            end if

          end if

        end do
        end do

      end do

!$omp end do

!$omp end parallel

! -----

      end subroutine s_estimsfc

!-----7--------------------------------------------------------------7--

      end module m_estimsfc
