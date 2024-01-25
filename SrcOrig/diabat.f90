!***********************************************************************
      module m_diabat
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/11/24
!     Modification: 1999/12/15, 2000/01/17, 2000/04/18, 2000/06/01,
!                   2001/01/15, 2001/05/29, 2001/11/20, 2002/04/02,
!                   2002/12/02, 2003/01/04, 2003/04/30, 2003/05/19,
!                   2003/08/01, 2003/11/28, 2003/12/12, 2004/09/10,
!                   2005/08/05, 2005/10/05, 2005/12/13, 2006/01/10,
!                   2006/02/13, 2006/04/03, 2006/06/21, 2006/11/06,
!                   2007/01/20, 2007/10/19, 2008/05/02, 2008/06/09,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the diabatic in the pressure equation.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_advp
      use m_comindx
      use m_comphy
      use m_getiname
      use m_phy2cnt
      use m_smoo2p
      use m_smoo4p

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: diabat, s_diabat

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface diabat

        module procedure s_diabat

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_diabat(fpsmtopt,fpcphopt,                            &
     &                    fpiwest,fpieast,fpjsouth,fpjnorth,fmois,      &
     &                    dtb,ni,nj,nk,j31,j32,jcb,jcb8u,jcb8v,jcb8w,   &
     &                    mf,mf8u,mf8v,ptbr,rcsq,u,v,w,ptp,ptpp,ptpf,   &
     &                    qv,qvp,qvf,qall,qallp,qallf,pdia,wc,          &
     &                    ptv,ptvp,ptvf,tmp1,tmp2,tmp3,tmp4,tmp5)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: fpsmtopt
                       ! Formal parameter of unique index of smtopt

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fpiwest
                       ! Formal parameter of unique index of iwest

      integer, intent(in) :: fpieast
                       ! Formal parameter of unique index of ieast

      integer, intent(in) :: fpjsouth
                       ! Formal parameter of unique index of jsouth

      integer, intent(in) :: fpjnorth
                       ! Formal parameter of unique index of jnorth

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in x direction

      integer, intent(in) :: nk
                       ! Model dimension in x direction

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: j31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of Jacobian

      real, intent(in) :: j32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of Jacobian

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: jcb8u(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at u points

      real, intent(in) :: jcb8v(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at v points

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: mf8u(0:ni+1,0:nj+1)
                       ! Map scale factors at u points

      real, intent(in) :: mf8v(0:ni+1,0:nj+1)
                       ! Map scale factors at v points

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: rcsq(0:ni+1,0:nj+1,1:nk)
                       ! rbr x sound wave speed squared

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at present

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at present

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at present

      real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at present

      real, intent(in) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(in) :: ptpf(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at future

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at present

      real, intent(in) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(in) :: qvf(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at future

! Input and output variables

      real, intent(inout) :: qall(0:ni+1,0:nj+1,1:nk)
                       ! Total water and ice mixing ratio at present

      real, intent(inout) :: qallp(0:ni+1,0:nj+1,1:nk)
                       ! Total water and ice mixing ratio at past

      real, intent(inout) :: qallf(0:ni+1,0:nj+1,1:nk)
                       ! Total water and ice mixing ratio at future

! Output variable

      real, intent(out) :: pdia(0:ni+1,0:nj+1,1:nk)
                       ! Diabatic value

! Internal shared variables

      integer smtopt   ! Option for numerical smoothing
      integer cphopt   ! Option for cloud micro physics

      integer iwest    ! Added index on west boundary
      integer jsouth   ! Added index on south boundary

      integer ieast    ! Subtracted index on east boundary
      integer jnorth   ! Subtracted index on north boundary

      real dtb2iv      ! 0.5 / dtb

      real, intent(inout) :: wc(0:ni+1,0:nj+1,1:nk)
                       ! zeta components of contravariant velocity

      real, intent(inout) :: ptv(0:ni+1,0:nj+1,1:nk)
                       ! Virtual potential temperature at present

      real, intent(inout) :: ptvp(0:ni+1,0:nj+1,1:nk)
                       ! Virtual potential temperature at past

      real, intent(inout) :: ptvf(0:ni+1,0:nj+1,1:nk)
                       ! Virtual potential temperature at future

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp4(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp5(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

! Remark

!     qall,qallp,qallf: These variables are also temporary, because they
!                       are not used again.

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpsmtopt,smtopt)
      call getiname(fpcphopt,cphopt)
      call getiname(fpiwest,iwest)
      call getiname(fpieast,ieast)
      call getiname(fpjsouth,jsouth)
      call getiname(fpjnorth,jnorth)

! -----

! Set the common used variable.

      dtb2iv=.5e0/dtb

! -----

!! Calculate the diabatic in the pressure equation.

! Cauculate the virtual potential temperature.

!$omp parallel default(shared) private(k)

      if(fmois(1:3).eq.'dry') then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=jsouth,nj-jnorth
          do i=iwest,ni-ieast
            ptv(i,j,k)=ptbr(i,j,k)+ptp(i,j,k)
            ptvp(i,j,k)=ptbr(i,j,k)+ptpp(i,j,k)
            ptvf(i,j,k)=ptbr(i,j,k)+ptpf(i,j,k)
          end do
          end do

!$omp end do

        end do

      else if(fmois(1:5).eq.'moist') then

        if(abs(cphopt).eq.0) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=jsouth,nj-jnorth
            do i=iwest,ni-ieast

              ptv(i,j,k)=(ptbr(i,j,k)+ptp(i,j,k))                       &
     &          *(1.e0+epsav*qv(i,j,k))/(1.e0+qv(i,j,k))

              ptvp(i,j,k)=(ptbr(i,j,k)+ptpp(i,j,k))                     &
     &          *(1.e0+epsav*qvp(i,j,k))/(1.e0+qvp(i,j,k))

              ptvf(i,j,k)=(ptbr(i,j,k)+ptpf(i,j,k))                     &
     &          *(1.e0+epsav*qvf(i,j,k))/(1.e0+qvf(i,j,k))

            end do
            end do

!$omp end do

          end do

        else

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=jsouth,nj-jnorth
            do i=iwest,ni-ieast

              ptv(i,j,k)=(ptbr(i,j,k)+ptp(i,j,k))                       &
     &          *(1.e0+epsav*qv(i,j,k))/(1.e0+qv(i,j,k)+qall(i,j,k))

              ptvp(i,j,k)=(ptbr(i,j,k)+ptpp(i,j,k))                     &
     &          *(1.e0+epsav*qvp(i,j,k))/(1.e0+qvp(i,j,k)+qallp(i,j,k))

              ptvf(i,j,k)=(ptbr(i,j,k)+ptpf(i,j,k))                     &
     &          *(1.e0+epsav*qvf(i,j,k))/(1.e0+qvf(i,j,k)+qallf(i,j,k))

            end do
            end do

!$omp end do

          end do

        end if

      end if

!$omp end parallel

! -----

! Calculate the zeta components of contravariant velocity.

      call s_phy2cnt(idsthopt,idtrnopt,idmpopt,idmfcopt,idoneopt,       &
     &               ni,nj,nk,j31,j32,jcb8w,mf,u,v,w,wc,tmp1,tmp2,tmp3)

! -----

! Calculate the advection term.

      call advp(idadvopt,idmpopt,idmfcopt,idzeropt,                     &
     &          idiwest,idieast,idjsouth,idjnorth,iddxiv,iddyiv,iddziv, &
     &          ni,nj,nk,mf8u,mf8v,jcb8u,jcb8v,jcb8w,u,v,wc,ptv,pdia,   &
     &          tmp1,tmp2,tmp3,tmp4,tmp5,qall,qallp,qallf)

! -----

! Calculate the 2nd order smoothing.

      if(mod(smtopt,10).eq.1) then

        call smoo2p(idsmhcoe,idsmvcoe,ni,nj,nk,ptvp,pdia)

! -----

! Calculate the 4th order smoothing.

      else if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

        call smoo4p(idsmtopt,idiwest,idieast,idjsouth,idjnorth,         &
     &              idsmhcoe,idsmvcoe,ni,nj,nk,ptvp,pdia,               &
     &              tmp1,tmp2,tmp3,tmp4)

      end if

! -----

!! Add the time tendency to the diabatic term and get the diabatic
!! value.

!$omp parallel default(shared) private(k)

! Add the time tendency to the diabatic term.

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          pdia(i,j,k)                                                   &
     &      =(ptvf(i,j,k)-ptvp(i,j,k))*jcb(i,j,k)*dtb2iv-pdia(i,j,k)
        end do
        end do

!$omp end do

      end do

! -----

! Finally get the diabatic value.

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          pdia(i,j,k)=rcsq(i,j,k)*pdia(i,j,k)/ptv(i,j,k)
        end do
        end do

!$omp end do

      end do

! -----

!$omp end parallel

!! -----

      end subroutine s_diabat

!-----7--------------------------------------------------------------7--

      end module m_diabat
