!***********************************************************************
      module m_timeflt
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/10/12, 1999/11/01,
!                   2000/01/17, 2000/04/18, 2000/06/01, 2002/04/02,
!                   2002/12/02, 2003/04/30, 2003/05/19, 2003/11/28,
!                   2003/12/12, 2004/02/01, 2004/05/31, 2004/08/20,
!                   2004/09/10, 2005/04/04, 2005/10/05, 2006/01/10,
!                   2006/02/13, 2006/07/21, 2007/01/20, 2007/10/19,
!                   2007/11/26, 2008/05/02, 2008/07/01, 2008/08/25,
!                   2009/01/30, 2009/02/27, 2011/08/18, 2011/09/22,
!                   2011/11/10, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform the Asselin time filter for the prognostic variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: timeflt, s_timeflt

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface timeflt

        module procedure s_timeflt

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_timeflt(fpsfcopt,fpcphopt,fphaiopt,fpqcgopt,         &
     &                  fpaslopt,fptrkopt,fptubopt,fpfilcoe,fmois,      &
     &                  dtsoil,ni,nj,nk,nqw,nnw,nqi,nni,nqa,nund,land,  &
     &                  up,uf,vp,vf,wp,wf,ppp,ppf,ptpp,ptpf,qvp,qvf,    &
     &                  qwtrp,qwtrf,nwtrp,nwtrf,qicep,qicef,nicep,nicef,&
     &                  qcwtrp,qcwtrf,qcicep,qcicef,qaslp,qaslf,qtp,qtf,&
     &                  tkep,tkef,tundp,tundf,u,v,w,pp,ptp,qv,          &
     &                  qwtr,nwtr,qice,nice,qcwtr,qcice,                &
     &                  qasl,qt,tke,tund)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: fpsfcopt
                       ! Formal parameter of unique index of sfcopt

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fphaiopt
                       ! Formal parameter of unique index of haiopt

      integer, intent(in) :: fpqcgopt
                       ! Formal parameter of unique index of qcgopt

      integer, intent(in) :: fpaslopt
                       ! Formal parameter of unique index of aslopt

      integer, intent(in) :: fptrkopt
                       ! Formal parameter of unique index of trkopt

      integer, intent(in) :: fptubopt
                       ! Formal parameter of unique index of tubopt

      integer, intent(in) :: fpfilcoe
                       ! Formal parameter of unique index of filcoe

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

      integer, intent(in) :: nnw
                       ! Number of categories of water concentrations

      integer, intent(in) :: nqi
                       ! Number of categories of ice hydrometeor

      integer, intent(in) :: nni
                       ! Number of categories of ice concentrations

      integer, intent(in) :: nqa(0:4)
                       ! Number of types of aerosol

      integer, intent(in) :: nund
                       ! Number of soil and sea layers

      integer, intent(in) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

      real, intent(in) :: dtsoil
                       ! Time interval of soil temperature calculation

      real, intent(in) :: up(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at past

      real, intent(in) :: uf(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at future

      real, intent(in) :: vp(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at past

      real, intent(in) :: vf(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at future

      real, intent(in) :: wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

      real, intent(in) :: wf(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at future

      real, intent(in) :: ppp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at past

      real, intent(in) :: ppf(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at future

      real, intent(in) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(in) :: ptpf(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at future

      real, intent(in) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(in) :: qvf(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at future

      real, intent(in) :: qwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at past

      real, intent(in) :: qwtrf(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at future

      real, intent(in) :: nwtrp(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at past

      real, intent(in) :: nwtrf(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at future

      real, intent(in) :: qicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at past

      real, intent(in) :: qicef(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at future

      real, intent(in) :: nicep(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at past

      real, intent(in) :: nicef(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at future

      real, intent(in) :: qcwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at past

      real, intent(in) :: qcwtrf(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at future

      real, intent(in) :: qcicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at past

      real, intent(in) :: qcicef(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at future

      real, intent(in) :: qaslp(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at past

      real, intent(in) :: qaslf(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at future

      real, intent(in) :: qtp(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at past

      real, intent(in) :: qtf(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at future

      real, intent(in) :: tkep(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at past

      real, intent(in) :: tkef(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at future

      real, intent(in) :: tundp(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at past

      real, intent(in) :: tundf(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at future

! Input and output variables

      real, intent(inout) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at present

      real, intent(inout) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at present

      real, intent(inout) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at present

      real, intent(inout) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at present

      real, intent(inout) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at present

      real, intent(inout) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at present

      real, intent(inout) :: qwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at present

      real, intent(inout) :: nwtr(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at present

      real, intent(inout) :: qice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at present

      real, intent(inout) :: nice(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at present

      real, intent(inout) :: qcwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at present

      real, intent(inout) :: qcice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at present

      real, intent(inout) :: qasl(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at present

      real, intent(inout) :: qt(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at present

      real, intent(inout) :: tke(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at present

      real, intent(inout) :: tund(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at present

! Internal shared variables

      integer sfcopt   ! Option for surface physics
      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes
      integer qcgopt   ! Option for charging distribution
      integer aslopt   ! Option for aerosol processes
      integer trkopt   ! Option for mixing ratio tracking
      integer tubopt   ! Option for turbulent mixing

      real filcoe      ! Coefficient of Asselin time filter

      real fc2         ! 2.0 x filcoe

      real m1fc2       ! 1.0 - 2.0 x filcoe
      real m1fc4       ! 1.0 - 4.0 x filcoe

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer n        ! Array index in 4th direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpsfcopt,sfcopt)
      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)
      call getiname(fpqcgopt,qcgopt)
      call getiname(fpaslopt,aslopt)
      call getiname(fptrkopt,trkopt)
      call getiname(fptubopt,tubopt)
      call getrname(fpfilcoe,filcoe)

! -----

! Set the common used variables.

      fc2=2.e0*filcoe

      m1fc2=1.e0-2.e0*filcoe
      m1fc4=1.e0-4.e0*filcoe

! -----

!!!! Perform the Asselin time filter.

!$omp parallel default(shared) private(k,n)

! Perform the Asselin time filter for the velocity.

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni
          u(i,j,k)=m1fc2*u(i,j,k)+filcoe*(uf(i,j,k)+up(i,j,k))
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=1,nj
        do i=1,ni-1
          v(i,j,k)=m1fc2*v(i,j,k)+filcoe*(vf(i,j,k)+vp(i,j,k))
        end do
        end do

!$omp end do

      end do

      do k=1,nk

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          w(i,j,k)=m1fc2*w(i,j,k)+filcoe*(wf(i,j,k)+wp(i,j,k))
        end do
        end do

!$omp end do

      end do

! -----

! Perform the Asselin time filter for the pressure and potential
! temperature.

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          pp(i,j,k)=m1fc2*pp(i,j,k)+filcoe*(ppf(i,j,k)+ppp(i,j,k))
          ptp(i,j,k)=m1fc2*ptp(i,j,k)+filcoe*(ptpf(i,j,k)+ptpp(i,j,k))
        end do
        end do

!$omp end do

      end do

! -----

!!! Perform the Asselin time filter for the hydrometeor.

! Perform the Asselin time filter for the water vapor mixing ratio.

      if(fmois(1:5).eq.'moist') then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            qv(i,j,k)=m1fc2*qv(i,j,k)+filcoe*(qvf(i,j,k)+qvp(i,j,k))
          end do
          end do

!$omp end do

        end do

! -----

!! For the bulk categories.

        if(abs(cphopt).lt.10) then

! Perform the Asselin time filter for the water hydrometeor.

          if(abs(cphopt).ge.1) then

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1
                qwtr(i,j,k,1)=m1fc2*qwtr(i,j,k,1)                       &
     &            +filcoe*(qwtrf(i,j,k,1)+qwtrp(i,j,k,1))

                qwtr(i,j,k,2)=m1fc2*qwtr(i,j,k,2)                       &
     &            +filcoe*(qwtrf(i,j,k,2)+qwtrp(i,j,k,2))

              end do
              end do

!$omp end do

            end do

          end if

! -----

! Perform the Asselin time filter for the water concentrations.

          if(abs(cphopt).eq.4) then

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1
                nwtr(i,j,k,1)=m1fc2*nwtr(i,j,k,1)                       &
     &            +filcoe*(nwtrf(i,j,k,1)+nwtrp(i,j,k,1))

                nwtr(i,j,k,2)=m1fc2*nwtr(i,j,k,2)                       &
     &            +filcoe*(nwtrf(i,j,k,2)+nwtrp(i,j,k,2))

              end do
              end do

!$omp end do

            end do

          end if

! -----

! Perform the Asselin time filter for the ice hydrometeor.

          if(abs(cphopt).ge.2) then

            if(haiopt.eq.0) then

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=1,nj-1
                do i=1,ni-1
                  qice(i,j,k,1)=m1fc2*qice(i,j,k,1)                     &
     &              +filcoe*(qicef(i,j,k,1)+qicep(i,j,k,1))

                  qice(i,j,k,2)=m1fc2*qice(i,j,k,2)                     &
     &              +filcoe*(qicef(i,j,k,2)+qicep(i,j,k,2))

                  qice(i,j,k,3)=m1fc2*qice(i,j,k,3)                     &
     &              +filcoe*(qicef(i,j,k,3)+qicep(i,j,k,3))

                end do
                end do

!$omp end do

              end do

            else

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=1,nj-1
                do i=1,ni-1
                  qice(i,j,k,1)=m1fc2*qice(i,j,k,1)                     &
     &              +filcoe*(qicef(i,j,k,1)+qicep(i,j,k,1))

                  qice(i,j,k,2)=m1fc2*qice(i,j,k,2)                     &
     &              +filcoe*(qicef(i,j,k,2)+qicep(i,j,k,2))

                  qice(i,j,k,3)=m1fc2*qice(i,j,k,3)                     &
     &              +filcoe*(qicef(i,j,k,3)+qicep(i,j,k,3))

                  qice(i,j,k,4)=m1fc2*qice(i,j,k,4)                     &
     &              +filcoe*(qicef(i,j,k,4)+qicep(i,j,k,4))

                end do
                end do

!$omp end do

              end do

            end if

          end if

! -----

! Perform the Asselin time filter for the ice concentrations.

          if(abs(cphopt).eq.2) then

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1
                nice(i,j,k,1)=m1fc2*nice(i,j,k,1)                       &
     &            +filcoe*(nicef(i,j,k,1)+nicep(i,j,k,1))
              end do
              end do

!$omp end do

            end do

          else if(abs(cphopt).ge.3) then

            if(haiopt.eq.0) then

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=1,nj-1
                do i=1,ni-1
                  nice(i,j,k,1)=m1fc2*nice(i,j,k,1)                     &
     &              +filcoe*(nicef(i,j,k,1)+nicep(i,j,k,1))

                  nice(i,j,k,2)=m1fc2*nice(i,j,k,2)                     &
     &              +filcoe*(nicef(i,j,k,2)+nicep(i,j,k,2))

                  nice(i,j,k,3)=m1fc2*nice(i,j,k,3)                     &
     &              +filcoe*(nicef(i,j,k,3)+nicep(i,j,k,3))

                end do
                end do

!$omp end do

              end do

            else

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=1,nj-1
                do i=1,ni-1
                  nice(i,j,k,1)=m1fc2*nice(i,j,k,1)                     &
     &              +filcoe*(nicef(i,j,k,1)+nicep(i,j,k,1))

                  nice(i,j,k,2)=m1fc2*nice(i,j,k,2)                     &
     &              +filcoe*(nicef(i,j,k,2)+nicep(i,j,k,2))

                  nice(i,j,k,3)=m1fc2*nice(i,j,k,3)                     &
     &              +filcoe*(nicef(i,j,k,3)+nicep(i,j,k,3))

                  nice(i,j,k,4)=m1fc2*nice(i,j,k,4)                     &
     &              +filcoe*(nicef(i,j,k,4)+nicep(i,j,k,4))

                end do
                end do

!$omp end do

              end do

            end if

          end if

! -----

! Perform the Asselin time filter for the charging distribution.

          if(cphopt.lt.0) then

            if(qcgopt.eq.2) then

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=1,nj-1
                do i=1,ni-1
                  qcwtr(i,j,k,1)=m1fc2*qcwtr(i,j,k,1)                   &
     &              +filcoe*(qcwtrf(i,j,k,1)+qcwtrp(i,j,k,1))

                  qcwtr(i,j,k,2)=m1fc2*qcwtr(i,j,k,2)                   &
     &              +filcoe*(qcwtrf(i,j,k,2)+qcwtrp(i,j,k,2))

                end do
                end do

!$omp end do

              end do

            end if

            if(haiopt.eq.0) then

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=1,nj-1
                do i=1,ni-1
                  qcice(i,j,k,1)=m1fc2*qcice(i,j,k,1)                   &
     &              +filcoe*(qcicef(i,j,k,1)+qcicep(i,j,k,1))

                  qcice(i,j,k,2)=m1fc2*qcice(i,j,k,2)                   &
     &              +filcoe*(qcicef(i,j,k,2)+qcicep(i,j,k,2))

                  qcice(i,j,k,3)=m1fc2*qcice(i,j,k,3)                   &
     &              +filcoe*(qcicef(i,j,k,3)+qcicep(i,j,k,3))

                end do
                end do

!$omp end do

              end do

            else

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=1,nj-1
                do i=1,ni-1
                  qcice(i,j,k,1)=m1fc2*qcice(i,j,k,1)                   &
     &              +filcoe*(qcicef(i,j,k,1)+qcicep(i,j,k,1))

                  qcice(i,j,k,2)=m1fc2*qcice(i,j,k,2)                   &
     &              +filcoe*(qcicef(i,j,k,2)+qcicep(i,j,k,2))

                  qcice(i,j,k,3)=m1fc2*qcice(i,j,k,3)                   &
     &              +filcoe*(qcicef(i,j,k,3)+qcicep(i,j,k,3))

                  qcice(i,j,k,4)=m1fc2*qcice(i,j,k,4)                   &
     &              +filcoe*(qcicef(i,j,k,4)+qcicep(i,j,k,4))

                end do
                end do

!$omp end do

              end do

            end if

          end if

! -----

!! -----

!! For the bin categories.

        else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

! Perform the Asselin time filter for the water hydrometeor.

          if(abs(cphopt).ge.11) then

            do n=1,nqw

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=1,nj-1
                do i=1,ni-1
                  qwtr(i,j,k,n)=m1fc2*qwtr(i,j,k,n)                     &
     &              +filcoe*(qwtrf(i,j,k,n)+qwtrp(i,j,k,n))
                end do
                end do

!$omp end do

              end do

            end do

            do n=1,nnw

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=1,nj-1
                do i=1,ni-1
                  nwtr(i,j,k,n)=m1fc2*nwtr(i,j,k,n)                     &
     &              +filcoe*(nwtrf(i,j,k,n)+nwtrp(i,j,k,n))
                end do
                end do

!$omp end do

              end do

            end do

          end if

! -----

! Perform the Asselin time filter for the ice hydrometeor.

          if(abs(cphopt).eq.12) then

            do n=1,nqi

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=1,nj-1
                do i=1,ni-1
                  qice(i,j,k,n)=m1fc2*qice(i,j,k,n)                     &
     &              +filcoe*(qicef(i,j,k,n)+qicep(i,j,k,n))
                end do
                end do

!$omp end do

              end do

            end do

            do n=1,nni

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=1,nj-1
                do i=1,ni-1
                  nice(i,j,k,n)=m1fc2*nice(i,j,k,n)                     &
     &              +filcoe*(nicef(i,j,k,n)+nicep(i,j,k,n))
                end do
                end do

!$omp end do

              end do

            end do

          end if

! -----

        end if

!! -----

      end if

!!! -----

! Perform the Asselin time filter for the aerosol.

      if(aslopt.ge.1) then

        do n=1,nqa(0)

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              qasl(i,j,k,n)=m1fc2*qasl(i,j,k,n)                         &
     &          +filcoe*(qaslf(i,j,k,n)+qaslp(i,j,k,n))
            end do
            end do

!$omp end do

          end do

        end do

      end if

! -----

! Perform the Asselin time filter for the tracer.

      if(trkopt.ge.1) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            qt(i,j,k)=m1fc2*qt(i,j,k)+filcoe*(qtf(i,j,k)+qtp(i,j,k))
          end do
          end do

!$omp end do

        end do

      end if

! -----

! Perform the Asselin time filter for the turbulent kinetic energy.

      if(tubopt.ge.2) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            tke(i,j,k)=m1fc4*tke(i,j,k)+fc2*(tkef(i,j,k)+tkep(i,j,k))
          end do
          end do

!$omp end do

        end do

      end if

! -----

! Perform the Asselin time filter for the soil and sea temperature.

      if(sfcopt.ge.1) then

        if(dtsoil.gt.0.e0) then

          if(sfcopt.eq.1.or.sfcopt.eq.11) then

            do k=1,nund

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1

                if(land(i,j).lt.3.or.land(i,j).ge.10) then

                  tund(i,j,k)=m1fc2*tund(i,j,k)                         &
     &              +filcoe*(tundf(i,j,k)+tundp(i,j,k))

                end if

              end do
              end do

!$omp end do

            end do

          else if(sfcopt.eq.2.or.sfcopt.eq.3.or.sfcopt.ge.12) then

            do k=1,nund

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1

                if(land(i,j).ge.10) then

                  tund(i,j,k)=m1fc2*tund(i,j,k)                         &
     &              +filcoe*(tundf(i,j,k)+tundp(i,j,k))

                end if

              end do
              end do

!$omp end do

            end do

          end if

        end if

      end if

! -----

!$omp end parallel

!!!! -----

      end subroutine s_timeflt

!-----7--------------------------------------------------------------7--

      end module m_timeflt
