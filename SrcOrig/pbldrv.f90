!***********************************************************************
      module m_pbldrv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/10/16
!     Modification: 2001/11/20, 2001/12/10, 2002/01/15, 2002/04/02,
!                   2002/08/15, 2002/12/27, 2003/01/20, 2003/04/30,
!                   2003/05/19, 2003/11/05, 2003/12/12, 2004/04/01,
!                   2004/05/07, 2004/09/10, 2007/01/20, 2007/06/21,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2011/09/22, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the inferior procedures for the virtical diffusion in the
!     planetary boundary layer.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_comphy
      use m_eddypbl
      use m_getiname
      use m_pblptv
      use m_pblqv
      use m_pblu
      use m_pblv

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: pbldrv, s_pbldrv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface pbldrv

        module procedure s_pbldrv

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
      subroutine s_pbldrv(fplevpbl,fmois,dtb,ni,nj,nk,zph,jcb8w,        &
     &                    ptbr,rbr,rst,rst8u,rst8v,ce,ct,cq,qvsfc,      &
     &                    ptv,u,v,ptp,qv,kms,khs,tmp1,tmp2,tmp3,tmp4)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: fplevpbl
                       ! Formal parameter of unique index of levpbl

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: rst8u(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at u points

      real, intent(in) :: rst8v(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at v points

      real, intent(in) :: ce(0:ni+1,0:nj+1)
                       ! Exchange coefficient of surface momentum flux

      real, intent(in) :: ct(0:ni+1,0:nj+1)
                       ! Exchange coefficient of surface heat flux

      real, intent(in) :: cq(0:ni+1,0:nj+1)
                       ! Exchange coefficient of surface moisture flux

      real, intent(in) :: qvsfc(0:ni+1,0:nj+1)
                       ! Water vapor mixing ratio on surface

! Input and output variables

      real, intent(inout) :: ptv(0:ni+1,0:nj+1,1:nk)
                       ! Virtual potential temperature

      real, intent(inout) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(inout) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(inout) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

      real, intent(inout) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

! Internal shared variables

      integer levpbl   ! Number of planetary boundary layer

      real, intent(inout) :: kms(0:ni+1,0:nj+1,1:nk)
                       ! Eddy viscosity in planetaty boundary layer

      real, intent(inout) :: khs(0:ni+1,0:nj+1,1:nk)
                       ! Eddy diffusivity in planetaty boundary layer

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp4(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fplevpbl,levpbl)

! -----

! Calculate the eddy viscosity and diffusivity in planetaty boundary
! layer.

      if(levpbl.ge.2) then

        call eddypbl(idlevpbl,iddz,ni,nj,nk,zph,jcb8w,u,v,ptv,kms,khs,  &
     &               tmp1)

      end if

! -----

! Calculate the virtical diffusion for the virtual potential
! temperature.

      call pblptv(idlevpbl,iddziv,dtb,ni,nj,nk,jcb8w,rbr,rst,ct,khs,ptv,&
     &            tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the virtical diffusion for the water vapor mixing ratio.

      if(fmois(1:5).eq.'moist') then

        call pblqv(idlevpbl,iddziv,dtb,ni,nj,nk,jcb8w,rbr,rst,qvsfc,    &
     &             cq,khs,qv,tmp1,tmp2,tmp3,tmp4)

      end if

! -----

! Calculate the virtical diffusion for the x components of velocity.

      call pblu(idlevpbl,iddziv,dtb,ni,nj,nk,jcb8w,rbr,rst8u,ce,kms,u,  &
     &          tmp1,tmp2,tmp3,tmp4,khs)

! -----

! Calculate the virtical diffusion for the y components of velocity.

      call pblv(idlevpbl,iddziv,dtb,ni,nj,nk,jcb8w,rbr,rst8v,ce,kms,v,  &
     &          tmp1,tmp2,tmp3,tmp4,khs)

! -----

! Finally convert the virtual potential temperature to the potential
! temperature.

!$omp parallel default(shared) private(k)

      if(fmois(1:3).eq.'dry') then

        do k=1,levpbl+1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            ptp(i,j,k)=ptv(i,j,k)-ptbr(i,j,k)
          end do
          end do

!$omp end do

        end do

      else

        do k=1,levpbl+1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            ptp(i,j,k)=ptv(i,j,k)                                       &
     &        *(1.e0+qv(i,j,k))/(1.e0+epsav*qv(i,j,k))-ptbr(i,j,k)
          end do
          end do

!$omp end do

        end do

      end if

!$omp end parallel

! -----

      end subroutine s_pbldrv

!-----7--------------------------------------------------------------7--

      end module m_pbldrv
