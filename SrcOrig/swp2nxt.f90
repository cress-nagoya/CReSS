!***********************************************************************
      module m_swp2nxt
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/04/06, 1999/08/03, 1999/08/23,
!                   1999/10/12, 1999/11/01, 2000/01/17, 2000/04/18,
!                   2000/06/01, 2002/04/02, 2002/12/02, 2003/04/30,
!                   2003/05/19, 2003/11/28, 2003/12/12, 2004/05/31,
!                   2004/08/20, 2005/01/31, 2005/10/05, 2006/01/10,
!                   2006/02/13, 2006/04/03, 2006/07/21, 2007/01/20,
!                   2007/10/19, 2007/11/26, 2008/05/02, 2008/07/01,
!                   2008/08/25, 2009/01/30, 2009/02/27, 2011/08/18,
!                   2011/09/22, 2011/11/10, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     swap the prognostic variables to the next time step.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: swp2nxt, s_swp2nxt

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface swp2nxt

        module procedure s_swp2nxt

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
      subroutine s_swp2nxt(fpsfcopt,fpadvopt,fpcphopt,fphaiopt,         &
     &                    fpqcgopt,fpaslopt,fptrkopt,fptubopt,          &
     &                    fpiwest,fpieast,fpjsouth,fpjnorth,fmois,      &
     &                    dtsoil,ni,nj,nk,nqw,nnw,nqi,nni,nqa,nund,     &
     &                    uf,vf,wf,ppf,ptpf,qvf,qwtrf,nwtrf,qicef,nicef,&
     &                    qcwtrf,qcicef,qaslf,qtf,tkef,tundf,u,v,w,     &
     &                    pp,ptp,qv,qwtr,nwtr,qice,nice,qcwtr,qcice,    &
     &                    qasl,qt,tke,tund,up,vp,wp,ppp,ptpp,qvp,       &
     &                    qwtrp,nwtrp,qicep,nicep,qcwtrp,qcicep,        &
     &                    qaslp,qtp,tkep,tundp)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: fpsfcopt
                       ! Formal parameter of unique index of sfcopt

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

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

      real, intent(in) :: dtsoil
                       ! Time interval of soil temperature calculation

      real, intent(in) :: uf(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at future

      real, intent(in) :: vf(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at future

      real, intent(in) :: wf(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at future

      real, intent(in) :: ppf(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at future

      real, intent(in) :: ptpf(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at future

      real, intent(in) :: qvf(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at future

      real, intent(in) :: qwtrf(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at future

      real, intent(in) :: nwtrf(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at future

      real, intent(in) :: qicef(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at future

      real, intent(in) :: nicef(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at future

      real, intent(in) :: qcwtrf(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at future

      real, intent(in) :: qcicef(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at future

      real, intent(in) :: qaslf(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at future

      real, intent(in) :: qtf(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at future

      real, intent(in) :: tkef(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at future

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

! Output variables

      real, intent(out) :: up(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at past

      real, intent(out) :: vp(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at past

      real, intent(out) :: wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

      real, intent(out) :: ppp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at past

      real, intent(out) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(out) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(out) :: qwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at past

      real, intent(out) :: nwtrp(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at past

      real, intent(out) :: qicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at past

      real, intent(out) :: nicep(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at past

      real, intent(out) :: qcwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at past

      real, intent(out) :: qcicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at past

      real, intent(out) :: qaslp(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at past

      real, intent(out) :: qtp(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at past

      real, intent(out) :: tkep(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at past

      real, intent(out) :: tundp(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at past

! Internal shared variables

      integer sfcopt   ! Option for surface physics
      integer advopt   ! Option for advection scheme
      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes
      integer qcgopt   ! Option for charging distribution
      integer aslopt   ! Option for aerosol processes
      integer trkopt   ! Option for mixing ratio tracking
      integer tubopt   ! Option for turbulent mixing

      integer iwest    ! Added index on west boundary
      integer ieast    ! Subtracted index on east boundary
      integer jsouth   ! Added index on south boundary
      integer jnorth   ! Subtracted index on north boundary

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer n        ! Array index in 4th direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpsfcopt,sfcopt)
      call getiname(fpadvopt,advopt)
      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)
      call getiname(fpqcgopt,qcgopt)
      call getiname(fpaslopt,aslopt)
      call getiname(fptrkopt,trkopt)
      call getiname(fptubopt,tubopt)
      call getiname(fpiwest,iwest)
      call getiname(fpieast,ieast)
      call getiname(fpjsouth,jsouth)
      call getiname(fpjnorth,jnorth)

! -----

!!!!! Swap the prognostic variables to the next time step.

!$omp parallel default(shared) private(k,n)

!!!! Swap the prognostic variables to the next time step in the case the
!!!! centered advection scheme is performed.

      if(advopt.le.3) then

! Swap the velocity prognostic variables to the next time step.

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=jsouth,nj-jnorth
          do i=iwest,ni+1-ieast
            up(i,j,k)=u(i,j,k)
            u(i,j,k)=uf(i,j,k)
          end do
          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=jsouth,nj+1-jnorth
          do i=iwest,ni-ieast
            vp(i,j,k)=v(i,j,k)
            v(i,j,k)=vf(i,j,k)
          end do
          end do

!$omp end do

        end do

        do k=1,nk

!$omp do schedule(runtime) private(i,j)

          do j=jsouth,nj-jnorth
          do i=iwest,ni-ieast
            wp(i,j,k)=w(i,j,k)
            w(i,j,k)=wf(i,j,k)
          end do
          end do

!$omp end do

        end do

! -----

! Swap the pressure and potential temperature perturbation to the next
! time step.

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=jsouth,nj-jnorth
          do i=iwest,ni-ieast
            ppp(i,j,k)=pp(i,j,k)
            ptpp(i,j,k)=ptp(i,j,k)

            pp(i,j,k)=ppf(i,j,k)
            ptp(i,j,k)=ptpf(i,j,k)

          end do
          end do

!$omp end do

        end do

! -----

!!! Swap the hydrometeor to the next time step.

        if(fmois(1:5).eq.'moist') then

! Swap the water vapor mixing raito to the next time step.

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=jsouth,nj-jnorth
            do i=iwest,ni-ieast
              qvp(i,j,k)=qv(i,j,k)
              qv(i,j,k)=qvf(i,j,k)
            end do
            end do

!$omp end do

          end do

! -----

!! For the bulk categories.

          if(abs(cphopt).lt.10) then

! Swap the water hydrometeor to the next time step.

            if(abs(cphopt).ge.1) then

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=jsouth,nj-jnorth
                do i=iwest,ni-ieast
                  qwtrp(i,j,k,1)=qwtr(i,j,k,1)
                  qwtrp(i,j,k,2)=qwtr(i,j,k,2)

                  qwtr(i,j,k,1)=qwtrf(i,j,k,1)
                  qwtr(i,j,k,2)=qwtrf(i,j,k,2)

                end do
                end do

!$omp end do

              end do

            end if

! -----

! Swap the water concentrations to the next time step.

            if(abs(cphopt).eq.4) then

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=jsouth,nj-jnorth
                do i=iwest,ni-ieast
                  nwtrp(i,j,k,1)=nwtr(i,j,k,1)
                  nwtrp(i,j,k,2)=nwtr(i,j,k,2)

                  nwtr(i,j,k,1)=nwtrf(i,j,k,1)
                  nwtr(i,j,k,2)=nwtrf(i,j,k,2)

                end do
                end do

!$omp end do

              end do

            end if

! -----

! Swap the ice hydrometeor to the next time step.

            if(abs(cphopt).ge.2) then

              if(haiopt.eq.0) then

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                  do j=jsouth,nj-jnorth
                  do i=iwest,ni-ieast
                    qicep(i,j,k,1)=qice(i,j,k,1)
                    qicep(i,j,k,2)=qice(i,j,k,2)
                    qicep(i,j,k,3)=qice(i,j,k,3)

                    qice(i,j,k,1)=qicef(i,j,k,1)
                    qice(i,j,k,2)=qicef(i,j,k,2)
                    qice(i,j,k,3)=qicef(i,j,k,3)

                  end do
                  end do

!$omp end do

                end do

              else

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                  do j=jsouth,nj-jnorth
                  do i=iwest,ni-ieast
                    qicep(i,j,k,1)=qice(i,j,k,1)
                    qicep(i,j,k,2)=qice(i,j,k,2)
                    qicep(i,j,k,3)=qice(i,j,k,3)
                    qicep(i,j,k,4)=qice(i,j,k,4)

                    qice(i,j,k,1)=qicef(i,j,k,1)
                    qice(i,j,k,2)=qicef(i,j,k,2)
                    qice(i,j,k,3)=qicef(i,j,k,3)
                    qice(i,j,k,4)=qicef(i,j,k,4)

                  end do
                  end do

!$omp end do

                end do

              end if

            end if

! -----

! Swap the ice concentrations to the next time step.

            if(abs(cphopt).eq.2) then

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=jsouth,nj-jnorth
                do i=iwest,ni-ieast
                  nicep(i,j,k,1)=nice(i,j,k,1)
                  nice(i,j,k,1)=nicef(i,j,k,1)
                end do
                end do

!$omp end do

              end do

            else if(abs(cphopt).ge.3) then

              if(haiopt.eq.0) then

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                  do j=jsouth,nj-jnorth
                  do i=iwest,ni-ieast
                    nicep(i,j,k,1)=nice(i,j,k,1)
                    nicep(i,j,k,2)=nice(i,j,k,2)
                    nicep(i,j,k,3)=nice(i,j,k,3)

                    nice(i,j,k,1)=nicef(i,j,k,1)
                    nice(i,j,k,2)=nicef(i,j,k,2)
                    nice(i,j,k,3)=nicef(i,j,k,3)

                  end do
                  end do

!$omp end do

                end do

              else

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                  do j=jsouth,nj-jnorth
                  do i=iwest,ni-ieast
                    nicep(i,j,k,1)=nice(i,j,k,1)
                    nicep(i,j,k,2)=nice(i,j,k,2)
                    nicep(i,j,k,3)=nice(i,j,k,3)
                    nicep(i,j,k,4)=nice(i,j,k,4)

                    nice(i,j,k,1)=nicef(i,j,k,1)
                    nice(i,j,k,2)=nicef(i,j,k,2)
                    nice(i,j,k,3)=nicef(i,j,k,3)
                    nice(i,j,k,4)=nicef(i,j,k,4)

                  end do
                  end do

!$omp end do

                end do

              end if

            end if

! -----

! Swap the charging distributions to the next time step.

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                  do j=jsouth,nj-jnorth
                  do i=iwest,ni-ieast
                    qcwtrp(i,j,k,1)=qcwtr(i,j,k,1)
                    qcwtrp(i,j,k,2)=qcwtr(i,j,k,2)

                    qcwtr(i,j,k,1)=qcwtrf(i,j,k,1)
                    qcwtr(i,j,k,2)=qcwtrf(i,j,k,2)

                  end do
                  end do

!$omp end do

                end do

              end if

              if(haiopt.eq.0) then

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                  do j=jsouth,nj-jnorth
                  do i=iwest,ni-ieast
                    qcicep(i,j,k,1)=qcice(i,j,k,1)
                    qcicep(i,j,k,2)=qcice(i,j,k,2)
                    qcicep(i,j,k,3)=qcice(i,j,k,3)

                    qcice(i,j,k,1)=qcicef(i,j,k,1)
                    qcice(i,j,k,2)=qcicef(i,j,k,2)
                    qcice(i,j,k,3)=qcicef(i,j,k,3)

                  end do
                  end do

!$omp end do

                end do

              else

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                  do j=jsouth,nj-jnorth
                  do i=iwest,ni-ieast
                    qcicep(i,j,k,1)=qcice(i,j,k,1)
                    qcicep(i,j,k,2)=qcice(i,j,k,2)
                    qcicep(i,j,k,3)=qcice(i,j,k,3)
                    qcicep(i,j,k,4)=qcice(i,j,k,4)

                    qcice(i,j,k,1)=qcicef(i,j,k,1)
                    qcice(i,j,k,2)=qcicef(i,j,k,2)
                    qcice(i,j,k,3)=qcicef(i,j,k,3)
                    qcice(i,j,k,4)=qcicef(i,j,k,4)

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

! Swap the water hydrometeor to the next time step.

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                  do j=jsouth,nj-jnorth
                  do i=iwest,ni-ieast
                    qwtrp(i,j,k,n)=qwtr(i,j,k,n)
                    qwtr(i,j,k,n)=qwtrf(i,j,k,n)
                  end do
                  end do

!$omp end do

                end do

              end do

              do n=1,nnw

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                  do j=jsouth,nj-jnorth
                  do i=iwest,ni-ieast
                    nwtrp(i,j,k,n)=nwtr(i,j,k,n)
                    nwtr(i,j,k,n)=nwtrf(i,j,k,n)
                  end do
                  end do

!$omp end do

                end do

              end do

            end if

! -----

! Swap the ice hydrometeor to the next time step.

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                  do j=jsouth,nj-jnorth
                  do i=iwest,ni-ieast
                    qicep(i,j,k,n)=qice(i,j,k,n)
                    qice(i,j,k,n)=qicef(i,j,k,n)
                  end do
                  end do

!$omp end do

                end do

              end do

              do n=1,nni

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                  do j=jsouth,nj-jnorth
                  do i=iwest,ni-ieast
                    nicep(i,j,k,n)=nice(i,j,k,n)
                    nice(i,j,k,n)=nicef(i,j,k,n)
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

! Swap the aerosol to the next time step.

        if(aslopt.ge.1) then

          do n=1,nqa(0)

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=jsouth,nj-jnorth
              do i=iwest,ni-ieast
                qaslp(i,j,k,n)=qasl(i,j,k,n)
                qasl(i,j,k,n)=qaslf(i,j,k,n)
              end do
              end do

!$omp end do

            end do

          end do

        end if

! -----

! Swap the tracer to the next time step.

        if(trkopt.ge.1) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=jsouth,nj-jnorth
            do i=iwest,ni-ieast
              qtp(i,j,k)=qt(i,j,k)
              qt(i,j,k)=qtf(i,j,k)
            end do
            end do

!$omp end do

          end do

        end if

! -----

! Swap the turbulent kinetic energy to the next time step.

        if(tubopt.ge.2) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=jsouth,nj-jnorth
            do i=iwest,ni-ieast
              tkep(i,j,k)=tke(i,j,k)
              tke(i,j,k)=tkef(i,j,k)
            end do
            end do

!$omp end do

          end do

        end if

! -----

! Swap the soil and sea temperature to the next time step.

        if(sfcopt.ge.1) then

          if(dtsoil.gt.0.e0) then

            do k=1,nund

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1
                tundp(i,j,k)=tund(i,j,k)
                tund(i,j,k)=tundf(i,j,k)
              end do
              end do

!$omp end do

            end do

          end if

        end if

! -----

!!!! ----

!!!! Swap the prognostic variables to the next time step in the case the
!!!! Cubic Lagrange advection scheme is performed.

      else

! Swap the velocity prognostic variables to the next time step.

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=jsouth,nj-jnorth
          do i=iwest,ni+1-ieast
            up(i,j,k)=uf(i,j,k)
          end do
          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=jsouth,nj+1-jnorth
          do i=iwest,ni-ieast
            vp(i,j,k)=vf(i,j,k)
          end do
          end do

!$omp end do

        end do

        do k=1,nk

!$omp do schedule(runtime) private(i,j)

          do j=jsouth,nj-jnorth
          do i=iwest,ni-ieast
            wp(i,j,k)=wf(i,j,k)
          end do
          end do

!$omp end do

        end do

! -----

! Swap the pressure and potential temperature perturbation to the next
! time step.

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=jsouth,nj-jnorth
          do i=iwest,ni-ieast
            ppp(i,j,k)=ppf(i,j,k)
            ptpp(i,j,k)=ptpf(i,j,k)
          end do
          end do

!$omp end do

        end do

! -----

!!! Swap the hydrometeor to the next time step.

        if(fmois(1:5).eq.'moist') then

! Swap the water vapor mixing raito to the next time step.

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=jsouth,nj-jnorth
            do i=iwest,ni-ieast
              qvp(i,j,k)=qvf(i,j,k)
            end do
            end do

!$omp end do

          end do

! -----

!! For the bulk categories.

          if(abs(cphopt).lt.10) then

! Swap the water hydrometeor to the next time step.

            if(abs(cphopt).ge.1) then

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=jsouth,nj-jnorth
                do i=iwest,ni-ieast
                  qwtrp(i,j,k,1)=qwtrf(i,j,k,1)
                  qwtrp(i,j,k,2)=qwtrf(i,j,k,2)
                end do
                end do

!$omp end do

              end do

            end if

! -----

! Swap the water concentrations to the next time step.

            if(abs(cphopt).eq.4) then

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=jsouth,nj-jnorth
                do i=iwest,ni-ieast
                  nwtrp(i,j,k,1)=nwtrf(i,j,k,1)
                  nwtrp(i,j,k,2)=nwtrf(i,j,k,2)
                end do
                end do

!$omp end do

              end do

            end if

! -----

! Swap the ice hydrometeor to the next time step.

            if(abs(cphopt).ge.2) then

              if(haiopt.eq.0) then

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                  do j=jsouth,nj-jnorth
                  do i=iwest,ni-ieast
                    qicep(i,j,k,1)=qicef(i,j,k,1)
                    qicep(i,j,k,2)=qicef(i,j,k,2)
                    qicep(i,j,k,3)=qicef(i,j,k,3)
                  end do
                  end do

!$omp end do

                end do

              else

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                  do j=jsouth,nj-jnorth
                  do i=iwest,ni-ieast
                    qicep(i,j,k,1)=qicef(i,j,k,1)
                    qicep(i,j,k,2)=qicef(i,j,k,2)
                    qicep(i,j,k,3)=qicef(i,j,k,3)
                    qicep(i,j,k,4)=qicef(i,j,k,4)
                  end do
                  end do

!$omp end do

                end do

              end if

            end if

! -----

! Swap the ice concentrations to the next time step.

            if(abs(cphopt).eq.2) then

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=jsouth,nj-jnorth
                do i=iwest,ni-ieast
                  nicep(i,j,k,1)=nicef(i,j,k,1)
                end do
                end do

!$omp end do

              end do

            else if(abs(cphopt).ge.3) then

              if(haiopt.eq.0) then

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                  do j=jsouth,nj-jnorth
                  do i=iwest,ni-ieast
                    nicep(i,j,k,1)=nicef(i,j,k,1)
                    nicep(i,j,k,2)=nicef(i,j,k,2)
                    nicep(i,j,k,3)=nicef(i,j,k,3)
                  end do
                  end do

!$omp end do

                end do

              else

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                  do j=jsouth,nj-jnorth
                  do i=iwest,ni-ieast
                    nicep(i,j,k,1)=nicef(i,j,k,1)
                    nicep(i,j,k,2)=nicef(i,j,k,2)
                    nicep(i,j,k,3)=nicef(i,j,k,3)
                    nicep(i,j,k,4)=nicef(i,j,k,4)
                  end do
                  end do

!$omp end do

                end do

              end if

            end if

! -----

! Swap the charging distributions to the next time step.

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                  do j=jsouth,nj-jnorth
                  do i=iwest,ni-ieast
                    qcwtrp(i,j,k,1)=qcwtrf(i,j,k,1)
                    qcwtrp(i,j,k,2)=qcwtrf(i,j,k,2)
                  end do
                  end do

!$omp end do

                end do

              end if

              if(haiopt.eq.0) then

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                  do j=jsouth,nj-jnorth
                  do i=iwest,ni-ieast
                    qcicep(i,j,k,1)=qcicef(i,j,k,1)
                    qcicep(i,j,k,2)=qcicef(i,j,k,2)
                    qcicep(i,j,k,3)=qcicef(i,j,k,3)
                  end do
                  end do

!$omp end do

                end do

              else

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                  do j=jsouth,nj-jnorth
                  do i=iwest,ni-ieast
                    qcicep(i,j,k,1)=qcicef(i,j,k,1)
                    qcicep(i,j,k,2)=qcicef(i,j,k,2)
                    qcicep(i,j,k,3)=qcicef(i,j,k,3)
                    qcicep(i,j,k,4)=qcicef(i,j,k,4)
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

! Swap the water hydrometeor to the next time step.

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                  do j=jsouth,nj-jnorth
                  do i=iwest,ni-ieast
                    qwtrp(i,j,k,n)=qwtrf(i,j,k,n)
                  end do
                  end do

!$omp end do

                end do

              end do

              do n=1,nnw

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                  do j=jsouth,nj-jnorth
                  do i=iwest,ni-ieast
                    nwtrp(i,j,k,n)=nwtrf(i,j,k,n)
                  end do
                  end do

!$omp end do

                end do

              end do

            end if

! -----

! Swap the ice hydrometeor to the next time step.

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                  do j=jsouth,nj-jnorth
                  do i=iwest,ni-ieast
                    qicep(i,j,k,n)=qicef(i,j,k,n)
                  end do
                  end do

!$omp end do

                end do

              end do

              do n=1,nni

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                  do j=jsouth,nj-jnorth
                  do i=iwest,ni-ieast
                    nicep(i,j,k,n)=nicef(i,j,k,n)
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

! Swap the aerosol to the next time step.

        if(aslopt.ge.1) then

          do n=1,nqa(0)

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=jsouth,nj-jnorth
              do i=iwest,ni-ieast
                qaslp(i,j,k,n)=qaslf(i,j,k,n)
              end do
              end do

!$omp end do

            end do

          end do

        end if

! -----

! Swap the tracer to the next time step.

        if(trkopt.ge.1) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=jsouth,nj-jnorth
            do i=iwest,ni-ieast
              qtp(i,j,k)=qtf(i,j,k)
            end do
            end do

!$omp end do

          end do

        end if

! -----

! Swap the turbulent kinetic energy to the next time step.

        if(tubopt.ge.2) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=jsouth,nj-jnorth
            do i=iwest,ni-ieast
              tkep(i,j,k)=tkef(i,j,k)
            end do
            end do

!$omp end do

          end do

        end if

! -----

! Swap the soil and sea temperature to the next time step.

        if(sfcopt.ge.1) then

          if(dtsoil.gt.0.e0) then

            do k=1,nund

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1
                tundp(i,j,k)=tundf(i,j,k)
              end do
              end do

!$omp end do

            end do

          end if

        end if

! -----

      end if

!!!! ----

!$omp end parallel

!!!!! -----

      end subroutine s_swp2nxt

!-----7--------------------------------------------------------------7--

      end module m_swp2nxt
