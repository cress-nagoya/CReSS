!***********************************************************************
      module m_steps
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/01/25, 1999/03/25, 1999/06/07,
!                   1999/07/05, 1999/07/21, 1999/08/03, 1999/08/09,
!                   1999/08/18, 1999/08/23, 1999/09/06, 1999/09/30,
!                   1999/10/07, 1999/10/12, 1999/11/01, 1999/11/19,
!                   1999/12/06, 2000/01/17, 2000/02/02, 2000/04/18,
!                   2000/12/18, 2001/01/15, 2001/03/13, 2001/05/29,
!                   2001/06/06, 2001/07/13, 2001/08/07, 2001/09/13,
!                   2001/11/20, 2002/04/02, 2002/06/06, 2002/07/23,
!                   2002/08/15, 2002/10/31, 2003/03/21, 2003/04/30,
!                   2003/05/19, 2003/06/27, 2003/11/05, 2003/11/28,
!                   2003/12/12, 2004/01/09, 2004/04/01, 2004/04/15,
!                   2004/05/07, 2004/05/31, 2004/08/20, 2004/09/10,
!                   2005/01/07, 2005/02/10, 2005/10/05, 2006/01/10,
!                   2006/02/13, 2006/04/03, 2006/05/12, 2006/07/21,
!                   2006/09/21, 2006/09/30, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/05/07, 2007/10/19, 2007/11/26,
!                   2008/05/02, 2008/08/25, 2008/12/11, 2009/01/30,
!                   2009/02/27, 2009/03/23, 2009/11/13, 2010/04/01,
!                   2011/08/18, 2011/09/22, 2011/11/10, 2013/01/28,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     solve the scalar variables to the next time step.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bc4news
      use m_bcycle
      use m_combuf
      use m_comindx
      use m_exbcpt
      use m_exbcq
      use m_getbufgx
      use m_getbufgy
      use m_getbufsx
      use m_getbufsy
      use m_getcname
      use m_getiname
      use m_inichar
      use m_lbcs
      use m_putbufgx
      use m_putbufgy
      use m_putbufsx
      use m_putbufsy
      use m_rbcpt
      use m_rbcq
      use m_rbcqv
      use m_rbcs
      use m_rbcs0
      use m_shiftgx
      use m_shiftgy
      use m_shiftsx
      use m_shiftsy
      use m_vbcqcg
      use m_vbcs

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: steps, s_steps

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface steps

        module procedure s_steps

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic max
      intrinsic min

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_steps(fpgpvvar,fpexbvar,fpexbopt,fpgwmopt,fpadvopt,  &
     &                   fpcphopt,fphaiopt,fpqcgopt,fpaslopt,fptrkopt,  &
     &                   fptubopt,fmois,dtb,gtinc,atinc,ni,nj,nk,       &
     &                   nqw,nnw,nqi,nni,nqa,qvbr,rst,ptp,ptpp,qv,qvp,  &
     &                   qwtr,qwtrp,nwtr,nwtrp,qice,qicep,nice,nicep,   &
     &                   qcwtr,qcwtrp,qcice,qcicep,qasl,qaslp,qt,qtp,   &
     &                   tke,tkep,ptcpx,ptcpy,qvcpx,qvcpy,qwcpx,qwcpy,  &
     &                   nwcpx,nwcpy,qicpx,qicpy,nicpx,nicpy,           &
     &                   qcwcpx,qcwcpy,qcicpx,qcicpy,qacpx,qacpy,       &
     &                   qtcpx,qtcpy,tkecpx,tkecpy,ptpgpv,ptptd,        &
     &                   qvgpv,qvtd,qwgpv,qwtd,qigpv,qitd,qagpv,qatd,   &
     &                   ptfrc,qvfrc,qwtrf,nwtrf,qicef,nicef,           &
     &                   qcwtrf,qcicef,qaslf,qtf,tkef,                  &
     &                   ptpf,qvf,dtdrst)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: fpgpvvar
                       ! Formal parameter of unique index of gpvvar

      integer, intent(in) :: fpexbvar
                       ! Formal parameter of unique index of exbvar

      integer, intent(in) :: fpexbopt
                       ! Formal parameter of unique index of exbopt

      integer, intent(in) :: fpgwmopt
                       ! Formal parameter of unique index of gwmopt

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

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: gtinc
                       ! Lapse of forecast time
                       ! from GPV data reading

      real, intent(in) :: atinc
                       ! Lapse of forecast time
                       ! from aerosol data reading

      real, intent(in) :: qvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state water vapor mixing ratio

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at present

      real, intent(in) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at present

      real, intent(in) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(in) :: qwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at present

      real, intent(in) :: qwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at past

      real, intent(in) :: nwtr(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at present

      real, intent(in) :: nwtrp(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at past

      real, intent(in) :: qice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at present

      real, intent(in) :: qicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at past

      real, intent(in) :: nice(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at present

      real, intent(in) :: nicep(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at past

      real, intent(in) :: qcwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at present

      real, intent(in) :: qcwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at past

      real, intent(in) :: qcice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at present

      real, intent(in) :: qcicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at past

      real, intent(in) :: qasl(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at present

      real, intent(in) :: qaslp(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at past

      real, intent(in) :: qt(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at present

      real, intent(in) :: qtp(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at past

      real, intent(in) :: tke(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at present

      real, intent(in) :: tkep(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at past

      real, intent(in) :: ptcpx(1:nj,1:nk,1:2)
                       ! Phase speed of potential temperature
                       ! on west and east boundary

      real, intent(in) :: ptcpy(1:ni,1:nk,1:2)
                       ! Phase speed of potential temperature
                       ! on south and north boundary

      real, intent(in) :: qvcpx(1:nj,1:nk,1:2)
                       ! Phase speed of water vapor mixing ratio
                       ! on west and east boundary

      real, intent(in) :: qvcpy(1:ni,1:nk,1:2)
                       ! Phase speed of water vapor mixing ratio
                       ! on south and north boundary

      real, intent(in) :: qwcpx(1:nj,1:nk,1:2,1:nqw)
                       ! Phase speed of water hydrometeor
                       ! on west and east boundary

      real, intent(in) :: qwcpy(1:ni,1:nk,1:2,1:nqw)
                       ! Phase speed of water hydrometeor
                       ! on south and north boundary

      real, intent(in) :: nwcpx(1:nj,1:nk,1:2,1:nnw)
                       ! Phase speed of water concentrations
                       ! on west and east boundary

      real, intent(in) :: nwcpy(1:ni,1:nk,1:2,1:nnw)
                       ! Phase speed of water concentrations
                       ! on south and north boundary

      real, intent(in) :: qicpx(1:nj,1:nk,1:2,1:nqi)
                       ! Phase speed of ice hydrometeor
                       ! on west and east boundary

      real, intent(in) :: qicpy(1:ni,1:nk,1:2,1:nqi)
                       ! Phase speed of ice hydrometeor
                       ! on south and north boundary

      real, intent(in) :: nicpx(1:nj,1:nk,1:2,1:nni)
                       ! Phase speed of ice concentrations
                       ! on west and east boundary

      real, intent(in) :: nicpy(1:ni,1:nk,1:2,1:nni)
                       ! Phase speed of ice concentrations
                       ! on south and north boundary

      real, intent(in) :: qcwcpx(1:nj,1:nk,1:2,1:nqw)
                       ! Phase speed of charging distribution for water
                       ! on west and east boundary

      real, intent(in) :: qcwcpy(1:ni,1:nk,1:2,1:nqw)
                       ! Phase speed of charging distribution for water
                       ! on south and north boundary

      real, intent(in) :: qcicpx(1:nj,1:nk,1:2,1:nqi)
                       ! Phase speed of charging distribution for ice
                       ! on west and east boundary

      real, intent(in) :: qcicpy(1:ni,1:nk,1:2,1:nqi)
                       ! Phase speed of charging distribution for ice
                       ! on south and north boundary

      real, intent(in) :: qacpx(1:nj,1:nk,1:2,1:nqa(0))
                       ! Phase speed of aerosol mixing ratio
                       ! on west and east boundary

      real, intent(in) :: qacpy(1:ni,1:nk,1:2,1:nqa(0))
                       ! Phase speed of aerosol mixing ratio
                       ! on south and north boundary

      real, intent(in) :: qtcpx(1:nj,1:nk,1:2)
                       ! Phase speed of tracer mixing ratio
                       ! on west and east boundary

      real, intent(in) :: qtcpy(1:ni,1:nk,1:2)
                       ! Phase speed of tracer mixing ratio
                       ! on south and north boundary

      real, intent(in) :: tkecpx(1:nj,1:nk,1:2)
                       ! Phase speed of turbulent kinetic energy
                       ! on west and east boundary

      real, intent(in) :: tkecpy(1:ni,1:nk,1:2)
                       ! Phase speed of turbulent kinetic energy
                       ! on south and north boundary

      real, intent(in) :: ptpgpv(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation of GPV data
                       ! at marked time

      real, intent(in) :: ptptd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! potential temperature perturbation of GPV data

      real, intent(in) :: qvgpv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio of GPV data
                       ! at marked time

      real, intent(in) :: qvtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! water vapor mixing ratio of GPV data

      real, intent(in) :: qwgpv(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor of GPV data at marked time

      real, intent(in) :: qwtd(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Time tendency of water hydrometeor of GPV data

      real, intent(in) :: qigpv(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor of GPV data at marked time

      real, intent(in) :: qitd(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Time tendency of ice hydrometeor of GPV data

      real, intent(in) :: qagpv(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Mixing ratio of aerosol data at marked time

      real, intent(in) :: qatd(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Time tendency of mixing ratio of aerosol data

      real, intent(in) :: ptfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in potential temperature equation

      real, intent(in) :: qvfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in water vapor mixing ratio equation

! Input and output variables

      real, intent(inout) :: qwtrf(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at future

      real, intent(inout) :: nwtrf(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at future

      real, intent(inout) :: qicef(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at future

      real, intent(inout) :: nicef(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at future

      real, intent(inout) :: qcwtrf(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at future

      real, intent(inout) :: qcicef(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at future

      real, intent(inout) :: qaslf(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at future

      real, intent(inout) :: qtf(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at future

      real, intent(inout) :: tkef(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at future

! Output variables

      real, intent(out) :: ptpf(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at future

      real, intent(out) :: qvf(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at future

! Internal shared variables

      character(len=108) gpvvar
                       ! Control flag of input GPV data variables

      character(len=108) exbvar
                       ! Control flag of
                       ! extrenal boundary forced variables

      integer exbopt   ! Option for external boundary forcing
      integer gwmopt   ! Option for gravity wave mode integration
      integer advopt   ! Option for advection scheme
      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes
      integer qcgopt   ! Option for charging distribution
      integer aslopt   ! Option for aerosol processes
      integer trkopt   ! Option for mixing ratio tracking
      integer tubopt   ! Option for turbulent mixing

      integer n        ! Array index in 4th direction

      integer ib       ! Exchanging variables number
      integer nb       ! Number of exchanging variables

      real dtb2        ! 2.0 x dtb

      real, intent(inout) :: dtdrst(0:ni+1,0:nj+1,1:nk)
                       ! 2.0 x dtb / rst

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer n_sub    ! Substitute for n

!-----7--------------------------------------------------------------7--

! Initialize the character variables.

      call inichar(gpvvar)
      call inichar(exbvar)

! -----

! Get the required namelist variables.

      call getcname(fpgpvvar,gpvvar)
      call getcname(fpexbvar,exbvar)
      call getiname(fpexbopt,exbopt)
      call getiname(fpgwmopt,gwmopt)
      call getiname(fpadvopt,advopt)
      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)
      call getiname(fpqcgopt,qcgopt)
      call getiname(fpaslopt,aslopt)
      call getiname(fptrkopt,trkopt)
      call getiname(fptubopt,tubopt)

! -----

! Set common used variable.

      dtb2=2.e0*dtb

! -----

!! Solve the scalar variables to the next time step.

!$omp parallel default(shared) private(k,n_sub)

! Set common used variable.

      if(fmois(1:5).eq.'moist'.or.gwmopt.eq.0                           &
     &  .or.aslopt.ge.1.or.trkopt.ge.1.or.tubopt.ge.2                   &
     &  .or.(abs(cphopt).ge.1.and.abs(cphopt).lt.20)) then

        if(advopt.le.3) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-2
              dtdrst(i,j,k)=dtb2/rst(i,j,k)
            end do
            end do

!$omp end do

          end do

        else

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-2
              dtdrst(i,j,k)=dtb/rst(i,j,k)
            end do
            end do

!$omp end do

          end do

        end if

      end if

! -----

! Solve the potential temperature perturbation to the next time step.

      if(gwmopt.eq.0) then

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            ptpf(i,j,k)=ptpp(i,j,k)+ptfrc(i,j,k)*dtdrst(i,j,k)
          end do
          end do

!$omp end do

        end do

      end if

! -----

! Solve the hydrometeor to the next time step.

      if(fmois(1:5).eq.'moist') then

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            qvf(i,j,k)=max(0.e0,qvp(i,j,k)+qvfrc(i,j,k)*dtdrst(i,j,k))
          end do
          end do

!$omp end do

        end do

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-2
                qwtrf(i,j,k,1)                                          &
     &            =max(0.e0,qwtrp(i,j,k,1)+qwtrf(i,j,k,1)*dtdrst(i,j,k))

                qwtrf(i,j,k,2)                                          &
     &            =max(0.e0,qwtrp(i,j,k,2)+qwtrf(i,j,k,2)*dtdrst(i,j,k))

              end do
              end do

!$omp end do

            end do

          end if

          if(abs(cphopt).eq.4) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-2
                nwtrf(i,j,k,1)                                          &
     &            =nwtrp(i,j,k,1)+nwtrf(i,j,k,1)*dtdrst(i,j,k)

                nwtrf(i,j,k,2)                                          &
     &            =nwtrp(i,j,k,2)+nwtrf(i,j,k,2)*dtdrst(i,j,k)

              end do
              end do

!$omp end do

            end do

          end if

          if(abs(cphopt).ge.2) then

            if(haiopt.eq.0) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  qicef(i,j,k,1)=max(0.e0,                              &
     &              qicep(i,j,k,1)+qicef(i,j,k,1)*dtdrst(i,j,k))

                  qicef(i,j,k,2)=max(0.e0,                              &
     &              qicep(i,j,k,2)+qicef(i,j,k,2)*dtdrst(i,j,k))

                  qicef(i,j,k,3)=max(0.e0,                              &
     &              qicep(i,j,k,3)+qicef(i,j,k,3)*dtdrst(i,j,k))

                end do
                end do

!$omp end do

              end do

            else

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  qicef(i,j,k,1)=max(0.e0,                              &
     &              qicep(i,j,k,1)+qicef(i,j,k,1)*dtdrst(i,j,k))

                  qicef(i,j,k,2)=max(0.e0,                              &
     &              qicep(i,j,k,2)+qicef(i,j,k,2)*dtdrst(i,j,k))

                  qicef(i,j,k,3)=max(0.e0,                              &
     &              qicep(i,j,k,3)+qicef(i,j,k,3)*dtdrst(i,j,k))

                  qicef(i,j,k,4)=max(0.e0,                              &
     &              qicep(i,j,k,4)+qicef(i,j,k,4)*dtdrst(i,j,k))

                end do
                end do

!$omp end do

              end do

            end if

          end if

          if(abs(cphopt).eq.2) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-2
                nicef(i,j,k,1)                                          &
     &            =nicep(i,j,k,1)+nicef(i,j,k,1)*dtdrst(i,j,k)
              end do
              end do

!$omp end do

            end do

          else if(abs(cphopt).ge.3) then

            if(haiopt.eq.0) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  nicef(i,j,k,1)                                        &
     &              =nicep(i,j,k,1)+nicef(i,j,k,1)*dtdrst(i,j,k)

                  nicef(i,j,k,2)                                        &
     &              =nicep(i,j,k,2)+nicef(i,j,k,2)*dtdrst(i,j,k)

                  nicef(i,j,k,3)                                        &
     &              =nicep(i,j,k,3)+nicef(i,j,k,3)*dtdrst(i,j,k)

                end do
                end do

!$omp end do

              end do

            else

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  nicef(i,j,k,1)                                        &
     &              =nicep(i,j,k,1)+nicef(i,j,k,1)*dtdrst(i,j,k)

                  nicef(i,j,k,2)                                        &
     &              =nicep(i,j,k,2)+nicef(i,j,k,2)*dtdrst(i,j,k)

                  nicef(i,j,k,3)                                        &
     &              =nicep(i,j,k,3)+nicef(i,j,k,3)*dtdrst(i,j,k)

                  nicef(i,j,k,4)                                        &
     &              =nicep(i,j,k,4)+nicef(i,j,k,4)*dtdrst(i,j,k)

                end do
                end do

!$omp end do

              end do

            end if

          end if

          if(cphopt.lt.0) then

            if(qcgopt.eq.2) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  qcwtrf(i,j,k,1)                                       &
     &              =qcwtrp(i,j,k,1)+qcwtrf(i,j,k,1)*dtdrst(i,j,k)

                  qcwtrf(i,j,k,2)                                       &
     &              =qcwtrp(i,j,k,2)+qcwtrf(i,j,k,2)*dtdrst(i,j,k)

                end do
                end do

!$omp end do

              end do

            end if

            if(haiopt.eq.0) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  qcicef(i,j,k,1)                                       &
     &              =qcicep(i,j,k,1)+qcicef(i,j,k,1)*dtdrst(i,j,k)

                  qcicef(i,j,k,2)                                       &
     &              =qcicep(i,j,k,2)+qcicef(i,j,k,2)*dtdrst(i,j,k)

                  qcicef(i,j,k,3)                                       &
     &              =qcicep(i,j,k,3)+qcicef(i,j,k,3)*dtdrst(i,j,k)

                end do
                end do

!$omp end do

              end do

            else

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  qcicef(i,j,k,1)                                       &
     &              =qcicep(i,j,k,1)+qcicef(i,j,k,1)*dtdrst(i,j,k)

                  qcicef(i,j,k,2)                                       &
     &              =qcicep(i,j,k,2)+qcicef(i,j,k,2)*dtdrst(i,j,k)

                  qcicef(i,j,k,3)                                       &
     &              =qcicep(i,j,k,3)+qcicef(i,j,k,3)*dtdrst(i,j,k)

                  qcicef(i,j,k,4)                                       &
     &              =qcicep(i,j,k,4)+qcicef(i,j,k,4)*dtdrst(i,j,k)

                end do
                end do

!$omp end do

              end do

            end if

          end if

        else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

          if(abs(cphopt).ge.11) then

            do n_sub=1,nqw

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  qwtrf(i,j,k,n_sub)                                    &
     &              =qwtrp(i,j,k,n_sub)+qwtrf(i,j,k,n_sub)*dtdrst(i,j,k)
                end do
                end do

!$omp end do

              end do

            end do

            do n_sub=1,nnw

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  nwtrf(i,j,k,n_sub)                                    &
     &              =nwtrp(i,j,k,n_sub)+nwtrf(i,j,k,n_sub)*dtdrst(i,j,k)
                end do
                end do

!$omp end do

              end do

            end do

          end if

          if(abs(cphopt).eq.12) then

            do n_sub=1,nqi

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  qicef(i,j,k,n_sub)                                    &
     &              =qicep(i,j,k,n_sub)+qicef(i,j,k,n_sub)*dtdrst(i,j,k)
                end do
                end do

!$omp end do

              end do

            end do

            do n_sub=1,nni

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  nicef(i,j,k,n_sub)                                    &
     &              =nicep(i,j,k,n_sub)+nicef(i,j,k,n_sub)*dtdrst(i,j,k)
                end do
                end do

!$omp end do

              end do

            end do

          end if

        end if

      end if

! -----

! Solve the aerosol to the next time step.

      if(aslopt.ge.1) then

        do n_sub=1,nqa(0)

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-2
              qaslf(i,j,k,n_sub)=max(0.e0,                              &
     &          qaslp(i,j,k,n_sub)+qaslf(i,j,k,n_sub)*dtdrst(i,j,k))
            end do
            end do

!$omp end do

          end do

        end do

      end if

! -----

! Solve the tracer to the next time step.

      if(trkopt.ge.1) then

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            qtf(i,j,k)=max(0.e0,qtp(i,j,k)+qtf(i,j,k)*dtdrst(i,j,k))
          end do
          end do

!$omp end do

        end do

      end if

! -----

! Solve the turbulent kinetic energy to the next time step.

      if(tubopt.ge.2) then

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            tkef(i,j,k)=max(0.e0,tkep(i,j,k)+tkef(i,j,k)*dtdrst(i,j,k))
          end do
          end do

!$omp end do

        end do

      end if

! -----

!$omp end parallel

!! -----

! Set the lateral boundary conditions.

      if(gwmopt.eq.0) then

        if(exbopt.ge.1.and.exbvar(5:5).ne.'x') then

          call exbcpt(idexbvar,idwbc,idebc,idadvopt,idexnews,1,dtb,     &
     &                gtinc,ni,nj,nk,ptp,ptpp,ptcpx,ptcpy,              &
     &                ptpgpv,ptptd,ptpf)

        else

          call rbcpt(idlbcvar,idwbc,idebc,idsbc,idnbc,idnggopt,idlspopt,&
     &               idvspopt,idadvopt,idlbnews,1,dtb,gtinc,ni,nj,nk,   &
     &               ptp,ptpp,ptcpx,ptcpy,ptpgpv,ptptd,ptpf)

          call lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,ptpf)

        end if

      end if

      if(fmois(1:5).eq.'moist') then

        if(exbopt.ge.1.and.exbvar(6:6).ne.'x') then

          call exbcq(idexbvar,idwbc,idebc,idadvopt,idexnews,6,1,dtb,    &
     &               gtinc,ni,nj,nk,qv,qvp,qvcpx,qvcpy,qvgpv,qvtd,qvf)

        else

          call rbcqv(idgpvvar,idlbcvar,idwbc,idebc,idsbc,idnbc,idnggopt,&
     &               idlspopt,idvspopt,idadvopt,idlbnews,1,dtb,gtinc,   &
     &               ni,nj,nk,qvbr,qv,qvp,qvcpx,qvcpy,qvgpv,qvtd,qvf)

          call lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,qvf)

        end if

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            if(exbopt.ge.1                                              &
     &        .and.gpvvar(3:3).eq.'o'.and.exbvar(7:7).ne.'x') then

              call s_exbcq(idexbvar,idwbc,idebc,idadvopt,idexnews,      &
     &                     7,1,dtb,gtinc,ni,nj,nk,qwtr(0,0,1,1),        &
     &                     qwtrp(0,0,1,1),qwcpx(1,1,1,1),qwcpy(1,1,1,1),&
     &                     qwgpv(0,0,1,1),qwtd(0,0,1,1),qwtrf(0,0,1,1))

            else

              call s_rbcq(idgpvvar,idlbcvar,idwbc,idebc,idsbc,idnbc,    &
     &                    idnggopt,idlspopt,idvspopt,idadvopt,idlbnews, &
     &                    3,7,1,dtb,gtinc,ni,nj,nk,qwtr(0,0,1,1),       &
     &                    qwtrp(0,0,1,1),qwcpx(1,1,1,1),qwcpy(1,1,1,1), &
     &                    qwgpv(0,0,1,1),qwtd(0,0,1,1),qwtrf(0,0,1,1))

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                    qwtrf(0,0,1,1))

            end if

            if(exbopt.ge.1                                              &
     &        .and.gpvvar(4:4).eq.'o'.and.exbvar(7:7).ne.'x') then

              call s_exbcq(idexbvar,idwbc,idebc,idadvopt,idexnews,      &
     &                     7,1,dtb,gtinc,ni,nj,nk,qwtr(0,0,1,2),        &
     &                     qwtrp(0,0,1,2),qwcpx(1,1,1,2),qwcpy(1,1,1,2),&
     &                     qwgpv(0,0,1,2),qwtd(0,0,1,2),qwtrf(0,0,1,2))

            else

              call s_rbcq(idgpvvar,idlbcvar,idwbc,idebc,idsbc,idnbc,    &
     &                    idnggopt,idlspopt,idvspopt,idadvopt,idlbnews, &
     &                    4,7,1,dtb,gtinc,ni,nj,nk,qwtr(0,0,1,2),       &
     &                    qwtrp(0,0,1,2),qwcpx(1,1,1,2),qwcpy(1,1,1,2), &
     &                    qwgpv(0,0,1,2),qwtd(0,0,1,2),qwtrf(0,0,1,2))

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                    qwtrf(0,0,1,2))

            end if

          end if

          if(abs(cphopt).eq.4) then

            call s_rbcs(idlbcvar,idwbc,idebc,idsbc,idnbc,               &
     &                  idadvopt,idlbnews,7,dtb,ni,nj,nk,nwtr(0,0,1,1), &
     &                  nwtrp(0,0,1,1),nwcpx(1,1,1,1),nwcpy(1,1,1,1),   &
     &                  nwtrf(0,0,1,1))

            call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,      &
     &                  nwtrf(0,0,1,1))

            call s_rbcs(idlbcvar,idwbc,idebc,idsbc,idnbc,               &
     &                  idadvopt,idlbnews,7,dtb,ni,nj,nk,nwtr(0,0,1,2), &
     &                  nwtrp(0,0,1,2),nwcpx(1,1,1,2),nwcpy(1,1,1,2),   &
     &                  nwtrf(0,0,1,2))

            call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,      &
     &                  nwtrf(0,0,1,2))

          end if

          if(abs(cphopt).ge.2) then

            if(exbopt.ge.1                                              &
     &        .and.gpvvar(5:5).eq.'o'.and.exbvar(7:7).ne.'x') then

              call s_exbcq(idexbvar,idwbc,idebc,idadvopt,idexnews,      &
     &                     7,1,dtb,gtinc,ni,nj,nk,qice(0,0,1,1),        &
     &                     qicep(0,0,1,1),qicpx(1,1,1,1),qicpy(1,1,1,1),&
     &                     qigpv(0,0,1,1),qitd(0,0,1,1),qicef(0,0,1,1))

            else

              call s_rbcq(idgpvvar,idlbcvar,idwbc,idebc,idsbc,idnbc,    &
     &                    idnggopt,idlspopt,idvspopt,idadvopt,idlbnews, &
     &                    5,7,1,dtb,gtinc,ni,nj,nk,qice(0,0,1,1),       &
     &                    qicep(0,0,1,1),qicpx(1,1,1,1),qicpy(1,1,1,1), &
     &                    qigpv(0,0,1,1),qitd(0,0,1,1),qicef(0,0,1,1))

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                    qicef(0,0,1,1))

            end if

            if(exbopt.ge.1                                              &
     &        .and.gpvvar(6:6).eq.'o'.and.exbvar(7:7).ne.'x') then

              call s_exbcq(idexbvar,idwbc,idebc,idadvopt,idexnews,      &
     &                     7,1,dtb,gtinc,ni,nj,nk,qice(0,0,1,2),        &
     &                     qicep(0,0,1,2),qicpx(1,1,1,2),qicpy(1,1,1,2),&
     &                     qigpv(0,0,1,2),qitd(0,0,1,2),qicef(0,0,1,2))

            else

              call s_rbcq(idgpvvar,idlbcvar,idwbc,idebc,idsbc,idnbc,    &
     &                    idnggopt,idlspopt,idvspopt,idadvopt,idlbnews, &
     &                    6,7,1,dtb,gtinc,ni,nj,nk,qice(0,0,1,2),       &
     &                    qicep(0,0,1,2),qicpx(1,1,1,2),qicpy(1,1,1,2), &
     &                    qigpv(0,0,1,2),qitd(0,0,1,2),qicef(0,0,1,2))

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                    qicef(0,0,1,2))

            end if

            if(exbopt.ge.1                                              &
     &        .and.gpvvar(7:7).eq.'o'.and.exbvar(7:7).ne.'x') then

              call s_exbcq(idexbvar,idwbc,idebc,idadvopt,idexnews,      &
     &                     7,1,dtb,gtinc,ni,nj,nk,qice(0,0,1,3),        &
     &                     qicep(0,0,1,3),qicpx(1,1,1,3),qicpy(1,1,1,3),&
     &                     qigpv(0,0,1,3),qitd(0,0,1,3),qicef(0,0,1,3))

            else

              call s_rbcq(idgpvvar,idlbcvar,idwbc,idebc,idsbc,idnbc,    &
     &                    idnggopt,idlspopt,idvspopt,idadvopt,idlbnews, &
     &                    7,7,1,dtb,gtinc,ni,nj,nk,qice(0,0,1,3),       &
     &                    qicep(0,0,1,3),qicpx(1,1,1,3),qicpy(1,1,1,3), &
     &                    qigpv(0,0,1,3),qitd(0,0,1,3),qicef(0,0,1,3))

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                    qicef(0,0,1,3))

            end if

            if(haiopt.eq.1) then

              if(exbopt.ge.1                                            &
     &          .and.gpvvar(8:8).eq.'o'.and.exbvar(7:7).ne.'x') then

                call s_exbcq(idexbvar,idwbc,idebc,idadvopt,idexnews,    &
     &                     7,1,dtb,gtinc,ni,nj,nk,qice(0,0,1,4),        &
     &                     qicep(0,0,1,4),qicpx(1,1,1,4),qicpy(1,1,1,4),&
     &                     qigpv(0,0,1,4),qitd(0,0,1,4),qicef(0,0,1,4))

              else

                call s_rbcq(idgpvvar,idlbcvar,idwbc,idebc,idsbc,idnbc,  &
     &                     idnggopt,idlspopt,idvspopt,idadvopt,idlbnews,&
     &                     8,7,1,dtb,gtinc,ni,nj,nk,qice(0,0,1,4),      &
     &                     qicep(0,0,1,4),qicpx(1,1,1,4),qicpy(1,1,1,4),&
     &                     qigpv(0,0,1,4),qitd(0,0,1,4),qicef(0,0,1,4))

                call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,  &
     &                     qicef(0,0,1,4))

              end if

            end if

          end if

          if(abs(cphopt).eq.2) then

            call s_rbcs(idlbcvar,idwbc,idebc,idsbc,idnbc,               &
     &                  idadvopt,idlbnews,7,dtb,ni,nj,nk,nice(0,0,1,1), &
     &                  nicep(0,0,1,1),nicpx(1,1,1,1),nicpy(1,1,1,1),   &
     &                  nicef(0,0,1,1))

            call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,      &
     &                  nicef(0,0,1,1))

          else if(abs(cphopt).ge.3) then

            call s_rbcs(idlbcvar,idwbc,idebc,idsbc,idnbc,               &
     &                  idadvopt,idlbnews,7,dtb,ni,nj,nk,nice(0,0,1,1), &
     &                  nicep(0,0,1,1),nicpx(1,1,1,1),nicpy(1,1,1,1),   &
     &                  nicef(0,0,1,1))

            call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,      &
     &                  nicef(0,0,1,1))

            call s_rbcs(idlbcvar,idwbc,idebc,idsbc,idnbc,               &
     &                  idadvopt,idlbnews,7,dtb,ni,nj,nk,nice(0,0,1,2), &
     &                  nicep(0,0,1,2),nicpx(1,1,1,2),nicpy(1,1,1,2),   &
     &                  nicef(0,0,1,2))

            call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,      &
     &                  nicef(0,0,1,2))

            call s_rbcs(idlbcvar,idwbc,idebc,idsbc,idnbc,               &
     &                  idadvopt,idlbnews,7,dtb,ni,nj,nk,nice(0,0,1,3), &
     &                  nicep(0,0,1,3),nicpx(1,1,1,3),nicpy(1,1,1,3),   &
     &                  nicef(0,0,1,3))

            call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,      &
     &                  nicef(0,0,1,3))

            if(haiopt.eq.1) then

              call s_rbcs(idlbcvar,idwbc,idebc,idsbc,idnbc,             &
     &                   idadvopt,idlbnews,7,dtb,ni,nj,nk,nice(0,0,1,4),&
     &                   nicep(0,0,1,4),nicpx(1,1,1,4),nicpy(1,1,1,4),  &
     &                   nicef(0,0,1,4))

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                   nicef(0,0,1,4))

            end if

          end if

          if(cphopt.lt.0) then

            if(qcgopt.eq.2) then

              call s_rbcs(idlbcvar,idwbc,idebc,idsbc,idnbc,             &
     &                  idadvopt,idlbnews,7,dtb,ni,nj,nk,qcwtr(0,0,1,1),&
     &                  qcwtrp(0,0,1,1),qcwcpx(1,1,1,1),qcwcpy(1,1,1,1),&
     &                  qcwtrf(0,0,1,1))

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                  qcwtrf(0,0,1,1))

              call s_rbcs(idlbcvar,idwbc,idebc,idsbc,idnbc,             &
     &                  idadvopt,idlbnews,7,dtb,ni,nj,nk,qcwtr(0,0,1,2),&
     &                  qcwtrp(0,0,1,2),qcwcpx(1,1,1,2),qcwcpy(1,1,1,2),&
     &                  qcwtrf(0,0,1,2))

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                  qcwtrf(0,0,1,2))

            end if

            call s_rbcs(idlbcvar,idwbc,idebc,idsbc,idnbc,               &
     &                  idadvopt,idlbnews,7,dtb,ni,nj,nk,qcice(0,0,1,1),&
     &                  qcicep(0,0,1,1),qcicpx(1,1,1,1),qcicpy(1,1,1,1),&
     &                  qcicef(0,0,1,1))

            call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,      &
     &                  qcicef(0,0,1,1))

            call s_rbcs(idlbcvar,idwbc,idebc,idsbc,idnbc,               &
     &                  idadvopt,idlbnews,7,dtb,ni,nj,nk,qcice(0,0,1,2),&
     &                  qcicep(0,0,1,2),qcicpx(1,1,1,2),qcicpy(1,1,1,2),&
     &                  qcicef(0,0,1,2))

            call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,      &
     &                  qcicef(0,0,1,2))

            call s_rbcs(idlbcvar,idwbc,idebc,idsbc,idnbc,               &
     &                  idadvopt,idlbnews,7,dtb,ni,nj,nk,qcice(0,0,1,3),&
     &                  qcicep(0,0,1,3),qcicpx(1,1,1,3),qcicpy(1,1,1,3),&
     &                  qcicef(0,0,1,3))

            call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,      &
     &                  qcicef(0,0,1,3))

            if(haiopt.eq.1) then

              call s_rbcs(idlbcvar,idwbc,idebc,idsbc,idnbc,             &
     &                  idadvopt,idlbnews,7,dtb,ni,nj,nk,qcice(0,0,1,4),&
     &                  qcicep(0,0,1,4),qcicpx(1,1,1,4),qcicpy(1,1,1,4),&
     &                  qcicef(0,0,1,4))

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                  qcicef(0,0,1,4))

            end if

          end if

        else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

          if(abs(cphopt).ge.11) then

            do n=1,nqw

              call s_rbcs(idlbcvar,idwbc,idebc,idsbc,idnbc,idadvopt,    &
     &                    idlbnews,7,dtb,ni,nj,nk,qwtr(0,0,1,n),        &
     &                    qwtrp(0,0,1,n),qwcpx(1,1,1,n),qwcpy(1,1,1,n), &
     &                    qwtrf(0,0,1,n))

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                    qwtrf(0,0,1,n))

            end do

            do n=1,nnw

              call s_rbcs(idlbcvar,idwbc,idebc,idsbc,idnbc,idadvopt,    &
     &                    idlbnews,7,dtb,ni,nj,nk,nwtr(0,0,1,n),        &
     &                    nwtrp(0,0,1,n),nwcpx(1,1,1,n),nwcpy(1,1,1,n), &
     &                    nwtrf(0,0,1,n))

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                    nwtrf(0,0,1,n))

            end do

          end if

          if(abs(cphopt).eq.12) then

            do n=1,nqi

              call s_rbcs(idlbcvar,idwbc,idebc,idsbc,idnbc,idadvopt,    &
     &                    idlbnews,7,dtb,ni,nj,nk,qice(0,0,1,n),        &
     &                    qicep(0,0,1,n),qicpx(1,1,1,n),qicpy(1,1,1,n), &
     &                    qicef(0,0,1,n))

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                    qicef(0,0,1,n))

            end do

            do n=1,nni

              call s_rbcs(idlbcvar,idwbc,idebc,idsbc,idnbc,idadvopt,    &
     &                    idlbnews,7,dtb,ni,nj,nk,nice(0,0,1,n),        &
     &                    nicep(0,0,1,n),nicpx(1,1,1,n),nicpy(1,1,1,n), &
     &                    nicef(0,0,1,n))

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                    nicef(0,0,1,n))

            end do

          end if

        end if

      end if

      if(aslopt.ge.1) then

        do n=1,nqa(0)

          if(exbopt.ge.1.and.exbvar(8:8).ne.'x') then

            call s_exbcq(idexbvar,idwbc,idebc,idadvopt,idexnews,        &
     &                   8,1,dtb,atinc,ni,nj,nk,qasl(0,0,1,n),          &
     &                   qaslp(0,0,1,n),qacpx(1,1,1,n),qacpy(1,1,1,n),  &
     &                   qagpv(0,0,1,n),qatd(0,0,1,n),qaslf(0,0,1,n))

          else

            call s_rbcq(idgpvvar,idlbcvar,idwbc,idebc,idsbc,idnbc,      &
     &                  idnggopt,idlspopt,idvspopt,idadvopt,idlbnews,   &
     &                  9,8,1,dtb,atinc,ni,nj,nk,qasl(0,0,1,n),         &
     &                  qaslp(0,0,1,n),qacpx(1,1,1,n),qacpy(1,1,1,n),   &
     &                  qagpv(0,0,1,n),qatd(0,0,1,n),qaslf(0,0,1,n))

            call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,      &
     &                  qaslf(0,0,1,n))

          end if

        end do

      end if

      if(trkopt.ge.1) then

        call rbcs0(idlbcvar,idwbc,idebc,idsbc,idnbc,idadvopt,idlbnews,  &
     &             9,dtb,ni,nj,nk,qt,qtp,qtcpx,qtcpy,qtf)

        call lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,qtf)

      end if

      if(tubopt.ge.2) then

        call rbcs0(idlbcvar,idwbc,idebc,idsbc,idnbc,idadvopt,idlbnews,  &
     &             10,dtb,ni,nj,nk,tke,tkep,tkecpx,tkecpy,tkef)

        call lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,tkef)

      end if

! -----

!!!! Exchange the value horizontally.

! Count the number of variables.

      nb=0

      if(gwmopt.eq.0) then
        nb=nb+1
      end if

      if(fmois(1:5).eq.'moist') then

        nb=nb+1

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then
            nb=nb+nqw
          end if

          if(abs(cphopt).ge.2) then

            if(abs(cphopt).eq.2) then
              nb=nb+nqi+1
            else
              nb=nb+nqi+nni
            end if

          end if

          if(abs(cphopt).eq.4) then
            nb=nb+nnw
          end if

          if(cphopt.lt.0) then

            if(qcgopt.eq.2) then
              nb=nb+nqw
            end if

            nb=nb+nqi

          end if

        else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

          if(abs(cphopt).ge.11) then
            nb=nb+nqw+nnw
          end if

          if(abs(cphopt).eq.12) then
            nb=nb+nqi+nni
          end if

        end if

      end if

      if(aslopt.ge.1) then
        nb=nb+nqa(0)
      end if

      if(trkopt.ge.1) then
        nb=nb+1
      end if

      if(tubopt.ge.2) then
        nb=nb+1
      end if

! -----

!!! Exchange the value horizontally between sub domain.

!! in x direction.

! Put the exchanging buffer.

      ib=0

      if(gwmopt.eq.0) then

        ib=ib+1

        call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,ptpf,         &
     &                  ib,nb,sbuf)

      end if

      if(fmois(1:5).eq.'moist') then

        ib=ib+1

        call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,qvf,          &
     &                  ib,nb,sbuf)

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            ib=ib+1

            call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qwtrf(0,0,1,1),ib,nb,sbuf)

            ib=ib+1

            call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qwtrf(0,0,1,2),ib,nb,sbuf)

          end if

          if(abs(cphopt).eq.4) then

            ib=ib+1

            call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      nwtrf(0,0,1,1),ib,nb,sbuf)

            ib=ib+1

            call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      nwtrf(0,0,1,2),ib,nb,sbuf)

          end if

          if(abs(cphopt).ge.2) then

            ib=ib+1

            call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qicef(0,0,1,1),ib,nb,sbuf)

            ib=ib+1

            call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qicef(0,0,1,2),ib,nb,sbuf)

            ib=ib+1

            call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qicef(0,0,1,3),ib,nb,sbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qicef(0,0,1,4),ib,nb,sbuf)

            end if

          end if

          if(abs(cphopt).eq.2) then

            ib=ib+1

            call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      nicef(0,0,1,1),ib,nb,sbuf)

          else if(abs(cphopt).ge.3) then

            ib=ib+1

            call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      nicef(0,0,1,1),ib,nb,sbuf)

            ib=ib+1

            call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      nicef(0,0,1,2),ib,nb,sbuf)

            ib=ib+1

            call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      nicef(0,0,1,3),ib,nb,sbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nicef(0,0,1,4),ib,nb,sbuf)

            end if

          end if

          if(cphopt.lt.0) then

            if(qcgopt.eq.2) then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qcwtrf(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qcwtrf(0,0,1,2),ib,nb,sbuf)

            end if

            ib=ib+1

            call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qcicef(0,0,1,1),ib,nb,sbuf)

            ib=ib+1

            call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qcicef(0,0,1,2),ib,nb,sbuf)

            ib=ib+1

            call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qcicef(0,0,1,3),ib,nb,sbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qcicef(0,0,1,4),ib,nb,sbuf)

            end if

          end if

        else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

          if(abs(cphopt).ge.11) then

            do n=1,nqw

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qwtrf(0,0,1,n),ib,nb,sbuf)

            end do

            do n=1,nnw

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nwtrf(0,0,1,n),ib,nb,sbuf)

            end do

          end if

          if(abs(cphopt).eq.12) then

            do n=1,nqi

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qicef(0,0,1,n),ib,nb,sbuf)

            end do

            do n=1,nni

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nicef(0,0,1,n),ib,nb,sbuf)

            end do

          end if

        end if

      end if

      if(aslopt.ge.1) then

        do n=1,nqa(0)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,            &
     &                    qaslf(0,0,1,n),ib,nb,sbuf)

        end do

      end if

      if(trkopt.ge.1) then

        ib=ib+1

        call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,qtf,          &
     &                  ib,nb,sbuf)

      end if

      if(tubopt.ge.2) then

        ib=ib+1

        call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,tkef,         &
     &                  ib,nb,sbuf)

      end if

! -----

! Call the exchanger.

      if(nb.ne.0) then

        call s_shiftsx(idwbc,idebc,'all',nj,nk,nb,sbuf,rbuf)

      end if

! -----

! Get the exchanging buffer.

      ib=0

      if(gwmopt.eq.0) then

        ib=ib+1

        call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,ptpf,         &
     &                  ib,nb,rbuf)

      end if

      if(fmois(1:5).eq.'moist') then

        ib=ib+1

        call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,qvf,          &
     &                  ib,nb,rbuf)

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            ib=ib+1

            call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qwtrf(0,0,1,1),ib,nb,rbuf)

            ib=ib+1

            call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qwtrf(0,0,1,2),ib,nb,rbuf)

          end if

          if(abs(cphopt).eq.4) then

            ib=ib+1

            call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      nwtrf(0,0,1,1),ib,nb,rbuf)

            ib=ib+1

            call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      nwtrf(0,0,1,2),ib,nb,rbuf)

          end if

          if(abs(cphopt).ge.2) then

            ib=ib+1

            call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qicef(0,0,1,1),ib,nb,rbuf)

            ib=ib+1

            call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qicef(0,0,1,2),ib,nb,rbuf)

            ib=ib+1

            call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qicef(0,0,1,3),ib,nb,rbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qicef(0,0,1,4),ib,nb,rbuf)

            end if

          end if

          if(abs(cphopt).eq.2) then

            ib=ib+1

            call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      nicef(0,0,1,1),ib,nb,rbuf)

          else if(abs(cphopt).ge.3) then

            ib=ib+1

            call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      nicef(0,0,1,1),ib,nb,rbuf)

            ib=ib+1

            call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      nicef(0,0,1,2),ib,nb,rbuf)

            ib=ib+1

            call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      nicef(0,0,1,3),ib,nb,rbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nicef(0,0,1,4),ib,nb,rbuf)

            end if

          end if

          if(cphopt.lt.0) then

            if(qcgopt.eq.2) then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qcwtrf(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qcwtrf(0,0,1,2),ib,nb,rbuf)

            end if

            ib=ib+1

            call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qcicef(0,0,1,1),ib,nb,rbuf)

            ib=ib+1

            call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qcicef(0,0,1,2),ib,nb,rbuf)

            ib=ib+1

            call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qcicef(0,0,1,3),ib,nb,rbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qcicef(0,0,1,4),ib,nb,rbuf)

            end if

          end if

        else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

          if(abs(cphopt).ge.11) then

            do n=1,nqw

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qwtrf(0,0,1,n),ib,nb,rbuf)

            end do

            do n=1,nnw

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nwtrf(0,0,1,n),ib,nb,rbuf)

            end do

          end if

          if(abs(cphopt).eq.12) then

            do n=1,nqi

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qicef(0,0,1,n),ib,nb,rbuf)

            end do

            do n=1,nni

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nicef(0,0,1,n),ib,nb,rbuf)

            end do

          end if

        end if

      end if

      if(aslopt.ge.1) then

        do n=1,nqa(0)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,            &
     &                    qaslf(0,0,1,n),ib,nb,rbuf)

        end do

      end if

      if(trkopt.ge.1) then

        ib=ib+1

        call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,qtf,          &
     &                  ib,nb,rbuf)

      end if

      if(tubopt.ge.2) then

        ib=ib+1

        call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,tkef,         &
     &                  ib,nb,rbuf)

      end if

! -----

!! -----

!! In y direction.

! Put the exchanging buffer.

      ib=0

      if(gwmopt.eq.0) then

        ib=ib+1

        call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,ptpf,         &
     &                  ib,nb,sbuf)

      end if

      if(fmois(1:5).eq.'moist') then

        ib=ib+1

        call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,qvf,          &
     &                  ib,nb,sbuf)

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            ib=ib+1

            call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      qwtrf(0,0,1,1),ib,nb,sbuf)

            ib=ib+1

            call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      qwtrf(0,0,1,2),ib,nb,sbuf)

          end if

          if(abs(cphopt).eq.4) then

            ib=ib+1

            call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      nwtrf(0,0,1,1),ib,nb,sbuf)

            ib=ib+1

            call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      nwtrf(0,0,1,2),ib,nb,sbuf)

          end if

          if(abs(cphopt).ge.2) then

            ib=ib+1

            call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      qicef(0,0,1,1),ib,nb,sbuf)

            ib=ib+1

            call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      qicef(0,0,1,2),ib,nb,sbuf)

            ib=ib+1

            call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      qicef(0,0,1,3),ib,nb,sbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qicef(0,0,1,4),ib,nb,sbuf)

            end if

          end if

          if(abs(cphopt).eq.2) then

            ib=ib+1

            call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      nicef(0,0,1,1),ib,nb,sbuf)

          else if(abs(cphopt).ge.3) then

            ib=ib+1

            call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      nicef(0,0,1,1),ib,nb,sbuf)

            ib=ib+1

            call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      nicef(0,0,1,2),ib,nb,sbuf)

            ib=ib+1

            call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      nicef(0,0,1,3),ib,nb,sbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        nicef(0,0,1,4),ib,nb,sbuf)

            end if

          end if

          if(cphopt.lt.0) then

            if(qcgopt.eq.2) then

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qcwtrf(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qcwtrf(0,0,1,2),ib,nb,sbuf)

            end if

            ib=ib+1

            call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      qcicef(0,0,1,1),ib,nb,sbuf)

            ib=ib+1

            call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      qcicef(0,0,1,2),ib,nb,sbuf)

            ib=ib+1

            call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      qcicef(0,0,1,3),ib,nb,sbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qcicef(0,0,1,4),ib,nb,sbuf)

            end if

          end if

        else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

          if(abs(cphopt).ge.11) then

            do n=1,nqw

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qwtrf(0,0,1,n),ib,nb,sbuf)

            end do

            do n=1,nnw

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        nwtrf(0,0,1,n),ib,nb,sbuf)

            end do

          end if

          if(abs(cphopt).eq.12) then

            do n=1,nqi

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qicef(0,0,1,n),ib,nb,sbuf)

            end do

            do n=1,nni

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        nicef(0,0,1,n),ib,nb,sbuf)

            end do

          end if

        end if

      end if

      if(aslopt.ge.1) then

        do n=1,nqa(0)

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,            &
     &                    qaslf(0,0,1,n),ib,nb,sbuf)

        end do

      end if

      if(trkopt.ge.1) then

        ib=ib+1

        call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,qtf,          &
     &                  ib,nb,sbuf)

      end if

      if(tubopt.ge.2) then

        ib=ib+1

        call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,tkef,         &
     &                  ib,nb,sbuf)

      end if

! -----

! Call the exchanger.

      if(nb.ne.0) then

        call s_shiftsy(idsbc,idnbc,'all',ni,nk,nb,sbuf,rbuf)

      end if

! -----

! Get the exchanging buffer.

      ib=0

      if(gwmopt.eq.0) then

        ib=ib+1

        call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,ptpf,         &
     &                  ib,nb,rbuf)

      end if

      if(fmois(1:5).eq.'moist') then

        ib=ib+1

        call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,qvf,          &
     &                  ib,nb,rbuf)

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            ib=ib+1

            call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      qwtrf(0,0,1,1),ib,nb,rbuf)

            ib=ib+1

            call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      qwtrf(0,0,1,2),ib,nb,rbuf)

          end if

          if(abs(cphopt).eq.4) then

            ib=ib+1

            call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      nwtrf(0,0,1,1),ib,nb,rbuf)

            ib=ib+1

            call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      nwtrf(0,0,1,2),ib,nb,rbuf)

          end if

          if(abs(cphopt).ge.2) then

            ib=ib+1

            call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      qicef(0,0,1,1),ib,nb,rbuf)

            ib=ib+1

            call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      qicef(0,0,1,2),ib,nb,rbuf)

            ib=ib+1

            call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      qicef(0,0,1,3),ib,nb,rbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qicef(0,0,1,4),ib,nb,rbuf)

            end if

          end if

          if(abs(cphopt).eq.2) then

            ib=ib+1

            call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      nicef(0,0,1,1),ib,nb,rbuf)

          else if(abs(cphopt).ge.3) then

            ib=ib+1

            call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      nicef(0,0,1,1),ib,nb,rbuf)

            ib=ib+1

            call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      nicef(0,0,1,2),ib,nb,rbuf)

            ib=ib+1

            call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      nicef(0,0,1,3),ib,nb,rbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        nicef(0,0,1,4),ib,nb,rbuf)

            end if

          end if

          if(cphopt.lt.0) then

            if(qcgopt.eq.2) then

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qcwtrf(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qcwtrf(0,0,1,2),ib,nb,rbuf)

            end if

            ib=ib+1

            call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      qcicef(0,0,1,1),ib,nb,rbuf)

            ib=ib+1

            call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      qcicef(0,0,1,2),ib,nb,rbuf)

            ib=ib+1

            call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      qcicef(0,0,1,3),ib,nb,rbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qcicef(0,0,1,4),ib,nb,rbuf)

            end if

          end if

        else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

          if(abs(cphopt).ge.11) then

            do n=1,nqw

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qwtrf(0,0,1,n),ib,nb,rbuf)

            end do

            do n=1,nnw

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        nwtrf(0,0,1,n),ib,nb,rbuf)

            end do

          end if

          if(abs(cphopt).eq.12) then

            do n=1,nqi

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qicef(0,0,1,n),ib,nb,rbuf)

            end do

            do n=1,nni

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        nicef(0,0,1,n),ib,nb,rbuf)

            end do

          end if

        end if

      end if

      if(aslopt.ge.1) then

        do n=1,nqa(0)

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,            &
     &                    qaslf(0,0,1,n),ib,nb,rbuf)

        end do

      end if

      if(trkopt.ge.1) then

        ib=ib+1

        call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,qtf,          &
     &                  ib,nb,rbuf)

      end if

      if(tubopt.ge.2) then

        ib=ib+1

        call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,tkef,         &
     &                  ib,nb,rbuf)

      end if

! -----

!! -----

!!! -----

!!! Exchange the value horizontally between group domain.

!! in x direction.

! Put the exchanging buffer.

      ib=0

      if(gwmopt.eq.0) then

        ib=ib+1

        call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,ptpf,         &
     &                  ib,nb,sbuf)

      end if

      if(fmois(1:5).eq.'moist') then

        ib=ib+1

        call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,qvf,          &
     &                  ib,nb,sbuf)

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qwtrf(0,0,1,1),ib,nb,sbuf)

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qwtrf(0,0,1,2),ib,nb,sbuf)

          end if

          if(abs(cphopt).eq.4) then

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      nwtrf(0,0,1,1),ib,nb,sbuf)

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      nwtrf(0,0,1,2),ib,nb,sbuf)

          end if

          if(abs(cphopt).ge.2) then

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qicef(0,0,1,1),ib,nb,sbuf)

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qicef(0,0,1,2),ib,nb,sbuf)

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qicef(0,0,1,3),ib,nb,sbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qicef(0,0,1,4),ib,nb,sbuf)

            end if

          end if

          if(abs(cphopt).eq.2) then

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      nicef(0,0,1,1),ib,nb,sbuf)

          else if(abs(cphopt).ge.3) then

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      nicef(0,0,1,1),ib,nb,sbuf)

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      nicef(0,0,1,2),ib,nb,sbuf)

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      nicef(0,0,1,3),ib,nb,sbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nicef(0,0,1,4),ib,nb,sbuf)

            end if

          end if

          if(cphopt.lt.0) then

            if(qcgopt.eq.2) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qcwtrf(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qcwtrf(0,0,1,2),ib,nb,sbuf)

            end if

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qcicef(0,0,1,1),ib,nb,sbuf)

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qcicef(0,0,1,2),ib,nb,sbuf)

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qcicef(0,0,1,3),ib,nb,sbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qcicef(0,0,1,4),ib,nb,sbuf)

            end if

          end if

        else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

          if(abs(cphopt).ge.11) then

            do n=1,nqw

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qwtrf(0,0,1,n),ib,nb,sbuf)

            end do

            do n=1,nnw

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nwtrf(0,0,1,n),ib,nb,sbuf)

            end do

          end if

          if(abs(cphopt).eq.12) then

            do n=1,nqi

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qicef(0,0,1,n),ib,nb,sbuf)

            end do

            do n=1,nni

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nicef(0,0,1,n),ib,nb,sbuf)

            end do

          end if

        end if

      end if

      if(aslopt.ge.1) then

        do n=1,nqa(0)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,            &
     &                    qaslf(0,0,1,n),ib,nb,sbuf)

        end do

      end if

      if(trkopt.ge.1) then

        ib=ib+1

        call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,qtf,          &
     &                  ib,nb,sbuf)

      end if

      if(tubopt.ge.2) then

        ib=ib+1

        call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,tkef,         &
     &                  ib,nb,sbuf)

      end if

! -----

! Call the exchanger.

      if(nb.ne.0) then

        call s_shiftgx(idwbc,idebc,'all',nj,nk,nb,sbuf,rbuf)

      end if

! -----

! Get the exchanging buffer.

      ib=0

      if(gwmopt.eq.0) then

        ib=ib+1

        call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,ptpf,         &
     &                  ib,nb,rbuf)

      end if

      if(fmois(1:5).eq.'moist') then

        ib=ib+1

        call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,qvf,          &
     &                  ib,nb,rbuf)

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qwtrf(0,0,1,1),ib,nb,rbuf)

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qwtrf(0,0,1,2),ib,nb,rbuf)

          end if

          if(abs(cphopt).eq.4) then

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      nwtrf(0,0,1,1),ib,nb,rbuf)

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      nwtrf(0,0,1,2),ib,nb,rbuf)

          end if

          if(abs(cphopt).ge.2) then

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qicef(0,0,1,1),ib,nb,rbuf)

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qicef(0,0,1,2),ib,nb,rbuf)

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qicef(0,0,1,3),ib,nb,rbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qicef(0,0,1,4),ib,nb,rbuf)

            end if

          end if

          if(abs(cphopt).eq.2) then

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      nicef(0,0,1,1),ib,nb,rbuf)

          else if(abs(cphopt).ge.3) then

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      nicef(0,0,1,1),ib,nb,rbuf)

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      nicef(0,0,1,2),ib,nb,rbuf)

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      nicef(0,0,1,3),ib,nb,rbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nicef(0,0,1,4),ib,nb,rbuf)

            end if

          end if

          if(cphopt.lt.0) then

            if(qcgopt.eq.2) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qcwtrf(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qcwtrf(0,0,1,2),ib,nb,rbuf)

            end if

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qcicef(0,0,1,1),ib,nb,rbuf)

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qcicef(0,0,1,2),ib,nb,rbuf)

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qcicef(0,0,1,3),ib,nb,rbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qcicef(0,0,1,4),ib,nb,rbuf)

            end if

          end if

        else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

          if(abs(cphopt).ge.11) then

            do n=1,nqw

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qwtrf(0,0,1,n),ib,nb,rbuf)

            end do

            do n=1,nnw

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nwtrf(0,0,1,n),ib,nb,rbuf)

            end do

          end if

          if(abs(cphopt).eq.12) then

            do n=1,nqi

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qicef(0,0,1,n),ib,nb,rbuf)

            end do

            do n=1,nni

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nicef(0,0,1,n),ib,nb,rbuf)

            end do

          end if

        end if

      end if

      if(aslopt.ge.1) then

        do n=1,nqa(0)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,            &
     &                    qaslf(0,0,1,n),ib,nb,rbuf)

        end do

      end if

      if(trkopt.ge.1) then

        ib=ib+1

        call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,qtf,          &
     &                  ib,nb,rbuf)

      end if

      if(tubopt.ge.2) then

        ib=ib+1

        call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,tkef,         &
     &                  ib,nb,rbuf)

      end if

! -----

!! -----

!! In y direction.

! Put the exchanging buffer.

      ib=0

      if(gwmopt.eq.0) then

        ib=ib+1

        call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,ptpf,         &
     &                  ib,nb,sbuf)

      end if

      if(fmois(1:5).eq.'moist') then

        ib=ib+1

        call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,qvf,          &
     &                  ib,nb,sbuf)

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            ib=ib+1

            call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      qwtrf(0,0,1,1),ib,nb,sbuf)

            ib=ib+1

            call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      qwtrf(0,0,1,2),ib,nb,sbuf)

          end if

          if(abs(cphopt).eq.4) then

            ib=ib+1

            call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      nwtrf(0,0,1,1),ib,nb,sbuf)

            ib=ib+1

            call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      nwtrf(0,0,1,2),ib,nb,sbuf)

          end if

          if(abs(cphopt).ge.2) then

            ib=ib+1

            call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      qicef(0,0,1,1),ib,nb,sbuf)

            ib=ib+1

            call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      qicef(0,0,1,2),ib,nb,sbuf)

            ib=ib+1

            call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      qicef(0,0,1,3),ib,nb,sbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qicef(0,0,1,4),ib,nb,sbuf)

            end if

          end if

          if(abs(cphopt).eq.2) then

            ib=ib+1

            call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      nicef(0,0,1,1),ib,nb,sbuf)

          else if(abs(cphopt).ge.3) then

            ib=ib+1

            call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      nicef(0,0,1,1),ib,nb,sbuf)

            ib=ib+1

            call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      nicef(0,0,1,2),ib,nb,sbuf)

            ib=ib+1

            call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      nicef(0,0,1,3),ib,nb,sbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        nicef(0,0,1,4),ib,nb,sbuf)

            end if

          end if

          if(cphopt.lt.0) then

            if(qcgopt.eq.2) then

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qcwtrf(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qcwtrf(0,0,1,2),ib,nb,sbuf)

            end if

            ib=ib+1

            call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      qcicef(0,0,1,1),ib,nb,sbuf)

            ib=ib+1

            call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      qcicef(0,0,1,2),ib,nb,sbuf)

            ib=ib+1

            call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      qcicef(0,0,1,3),ib,nb,sbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qcicef(0,0,1,4),ib,nb,sbuf)

            end if

          end if

        else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

          if(abs(cphopt).ge.11) then

            do n=1,nqw

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qwtrf(0,0,1,n),ib,nb,sbuf)

            end do

            do n=1,nnw

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        nwtrf(0,0,1,n),ib,nb,sbuf)

            end do

          end if

          if(abs(cphopt).eq.12) then

            do n=1,nqi

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qicef(0,0,1,n),ib,nb,sbuf)

            end do

            do n=1,nni

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        nicef(0,0,1,n),ib,nb,sbuf)

            end do

          end if

        end if

      end if

      if(aslopt.ge.1) then

        do n=1,nqa(0)

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,            &
     &                    qaslf(0,0,1,n),ib,nb,sbuf)

        end do

      end if

      if(trkopt.ge.1) then

        ib=ib+1

        call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,qtf,          &
     &                  ib,nb,sbuf)

      end if

      if(tubopt.ge.2) then

        ib=ib+1

        call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,tkef,         &
     &                  ib,nb,sbuf)

      end if

! -----

! Call the exchanger.

      if(nb.ne.0) then

        call s_shiftgy(idsbc,idnbc,'all',ni,nk,nb,sbuf,rbuf)

      end if

! -----

! Get the exchanging buffer.

      ib=0

      if(gwmopt.eq.0) then

        ib=ib+1

        call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,ptpf,         &
     &                  ib,nb,rbuf)

      end if

      if(fmois(1:5).eq.'moist') then

        ib=ib+1

        call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,qvf,          &
     &                  ib,nb,rbuf)

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            ib=ib+1

            call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      qwtrf(0,0,1,1),ib,nb,rbuf)

            ib=ib+1

            call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      qwtrf(0,0,1,2),ib,nb,rbuf)

          end if

          if(abs(cphopt).eq.4) then

            ib=ib+1

            call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      nwtrf(0,0,1,1),ib,nb,rbuf)

            ib=ib+1

            call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      nwtrf(0,0,1,2),ib,nb,rbuf)

          end if

          if(abs(cphopt).ge.2) then

            ib=ib+1

            call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      qicef(0,0,1,1),ib,nb,rbuf)

            ib=ib+1

            call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      qicef(0,0,1,2),ib,nb,rbuf)

            ib=ib+1

            call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      qicef(0,0,1,3),ib,nb,rbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qicef(0,0,1,4),ib,nb,rbuf)

            end if

          end if

          if(abs(cphopt).eq.2) then

            ib=ib+1

            call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      nicef(0,0,1,1),ib,nb,rbuf)

          else if(abs(cphopt).ge.3) then

            ib=ib+1

            call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      nicef(0,0,1,1),ib,nb,rbuf)

            ib=ib+1

            call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      nicef(0,0,1,2),ib,nb,rbuf)

            ib=ib+1

            call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      nicef(0,0,1,3),ib,nb,rbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        nicef(0,0,1,4),ib,nb,rbuf)

            end if

          end if

          if(cphopt.lt.0) then

            if(qcgopt.eq.2) then

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qcwtrf(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qcwtrf(0,0,1,2),ib,nb,rbuf)

            end if

            ib=ib+1

            call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      qcicef(0,0,1,1),ib,nb,rbuf)

            ib=ib+1

            call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      qcicef(0,0,1,2),ib,nb,rbuf)

            ib=ib+1

            call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      qcicef(0,0,1,3),ib,nb,rbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qcicef(0,0,1,4),ib,nb,rbuf)

            end if

          end if

        else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

          if(abs(cphopt).ge.11) then

            do n=1,nqw

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qwtrf(0,0,1,n),ib,nb,rbuf)

            end do

            do n=1,nnw

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        nwtrf(0,0,1,n),ib,nb,rbuf)

            end do

          end if

          if(abs(cphopt).eq.12) then

            do n=1,nqi

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qicef(0,0,1,n),ib,nb,rbuf)

            end do

            do n=1,nni

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        nicef(0,0,1,n),ib,nb,rbuf)

            end do

          end if

        end if

      end if

      if(aslopt.ge.1) then

        do n=1,nqa(0)

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,            &
     &                    qaslf(0,0,1,n),ib,nb,rbuf)

        end do

      end if

      if(trkopt.ge.1) then

        ib=ib+1

        call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,qtf,          &
     &                  ib,nb,rbuf)

      end if

      if(tubopt.ge.2) then

        ib=ib+1

        call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,tkef,         &
     &                  ib,nb,rbuf)

      end if

! -----

!! -----

!! in x direction again.

! Put the exchanging buffer.

      ib=0

      if(gwmopt.eq.0) then

        ib=ib+1

        call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,ptpf,         &
     &                  ib,nb,sbuf)

      end if

      if(fmois(1:5).eq.'moist') then

        ib=ib+1

        call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,qvf,          &
     &                  ib,nb,sbuf)

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qwtrf(0,0,1,1),ib,nb,sbuf)

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qwtrf(0,0,1,2),ib,nb,sbuf)

          end if

          if(abs(cphopt).eq.4) then

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      nwtrf(0,0,1,1),ib,nb,sbuf)

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      nwtrf(0,0,1,2),ib,nb,sbuf)

          end if

          if(abs(cphopt).ge.2) then

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qicef(0,0,1,1),ib,nb,sbuf)

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qicef(0,0,1,2),ib,nb,sbuf)

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qicef(0,0,1,3),ib,nb,sbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qicef(0,0,1,4),ib,nb,sbuf)

            end if

          end if

          if(abs(cphopt).eq.2) then

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      nicef(0,0,1,1),ib,nb,sbuf)

          else if(abs(cphopt).ge.3) then

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      nicef(0,0,1,1),ib,nb,sbuf)

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      nicef(0,0,1,2),ib,nb,sbuf)

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      nicef(0,0,1,3),ib,nb,sbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nicef(0,0,1,4),ib,nb,sbuf)

            end if

          end if

          if(cphopt.lt.0) then

            if(qcgopt.eq.2) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qcwtrf(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qcwtrf(0,0,1,2),ib,nb,sbuf)

            end if

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qcicef(0,0,1,1),ib,nb,sbuf)

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qcicef(0,0,1,2),ib,nb,sbuf)

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qcicef(0,0,1,3),ib,nb,sbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qcicef(0,0,1,4),ib,nb,sbuf)

            end if

          end if

        else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

          if(abs(cphopt).ge.11) then

            do n=1,nqw

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qwtrf(0,0,1,n),ib,nb,sbuf)

            end do

            do n=1,nnw

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nwtrf(0,0,1,n),ib,nb,sbuf)

            end do

          end if

          if(abs(cphopt).eq.12) then

            do n=1,nqi

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qicef(0,0,1,n),ib,nb,sbuf)

            end do

            do n=1,nni

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nicef(0,0,1,n),ib,nb,sbuf)

            end do

          end if

        end if

      end if

      if(aslopt.ge.1) then

        do n=1,nqa(0)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,            &
     &                    qaslf(0,0,1,n),ib,nb,sbuf)

        end do

      end if

      if(trkopt.ge.1) then

        ib=ib+1

        call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,qtf,          &
     &                  ib,nb,sbuf)

      end if

      if(tubopt.ge.2) then

        ib=ib+1

        call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,tkef,         &
     &                  ib,nb,sbuf)

      end if

! -----

! Call the exchanger.

      if(nb.ne.0) then

        call s_shiftgx(idwbc,idebc,'all',nj,nk,nb,sbuf,rbuf)

      end if

! -----

! Get the exchanging buffer.

      ib=0

      if(gwmopt.eq.0) then

        ib=ib+1

        call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,ptpf,         &
     &                  ib,nb,rbuf)

      end if

      if(fmois(1:5).eq.'moist') then

        ib=ib+1

        call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,qvf,          &
     &                  ib,nb,rbuf)

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qwtrf(0,0,1,1),ib,nb,rbuf)

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qwtrf(0,0,1,2),ib,nb,rbuf)

          end if

          if(abs(cphopt).eq.4) then

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      nwtrf(0,0,1,1),ib,nb,rbuf)

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      nwtrf(0,0,1,2),ib,nb,rbuf)

          end if

          if(abs(cphopt).ge.2) then

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qicef(0,0,1,1),ib,nb,rbuf)

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qicef(0,0,1,2),ib,nb,rbuf)

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qicef(0,0,1,3),ib,nb,rbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qicef(0,0,1,4),ib,nb,rbuf)

            end if

          end if

          if(abs(cphopt).eq.2) then

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      nicef(0,0,1,1),ib,nb,rbuf)

          else if(abs(cphopt).ge.3) then

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      nicef(0,0,1,1),ib,nb,rbuf)

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      nicef(0,0,1,2),ib,nb,rbuf)

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      nicef(0,0,1,3),ib,nb,rbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nicef(0,0,1,4),ib,nb,rbuf)

            end if

          end if

          if(cphopt.lt.0) then

            if(qcgopt.eq.2) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qcwtrf(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qcwtrf(0,0,1,2),ib,nb,rbuf)

            end if

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qcicef(0,0,1,1),ib,nb,rbuf)

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qcicef(0,0,1,2),ib,nb,rbuf)

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qcicef(0,0,1,3),ib,nb,rbuf)

            if(haiopt.eq.1) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qcicef(0,0,1,4),ib,nb,rbuf)

            end if

          end if

        else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

          if(abs(cphopt).ge.11) then

            do n=1,nqw

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qwtrf(0,0,1,n),ib,nb,rbuf)

            end do

            do n=1,nnw

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nwtrf(0,0,1,n),ib,nb,rbuf)

            end do

          end if

          if(abs(cphopt).eq.12) then

            do n=1,nqi

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qicef(0,0,1,n),ib,nb,rbuf)

            end do

            do n=1,nni

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nicef(0,0,1,n),ib,nb,rbuf)

            end do

          end if

        end if

      end if

      if(aslopt.ge.1) then

        do n=1,nqa(0)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,            &
     &                    qaslf(0,0,1,n),ib,nb,rbuf)

        end do

      end if

      if(trkopt.ge.1) then

        ib=ib+1

        call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,qtf,          &
     &                  ib,nb,rbuf)

      end if

      if(tubopt.ge.2) then

        ib=ib+1

        call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,tkef,         &
     &                  ib,nb,rbuf)

      end if

! -----

!! -----

!!! -----

!!!! -----

! Set the periodic boundary conditions.

      if(gwmopt.eq.0) then

        call bcycle(idwbc,idebc,idsbc,idnbc,                            &
     &              2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,ptpf)

      end if

      if(fmois(1:5).eq.'moist') then

        call bcycle(idwbc,idebc,idsbc,idnbc,                            &
     &              2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,qvf)

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            call s_bcycle(idwbc,idebc,idsbc,idnbc,                      &
     &                    2,1,ni-2,ni-1,2,1,nj-2,nj-1,                  &
     &                    ni,nj,nk,qwtrf(0,0,1,1))

            call s_bcycle(idwbc,idebc,idsbc,idnbc,                      &
     &                    2,1,ni-2,ni-1,2,1,nj-2,nj-1,                  &
     &                    ni,nj,nk,qwtrf(0,0,1,2))

          end if

          if(abs(cphopt).eq.4) then

            call s_bcycle(idwbc,idebc,idsbc,idnbc,                      &
     &                    2,1,ni-2,ni-1,2,1,nj-2,nj-1,                  &
     &                    ni,nj,nk,nwtrf(0,0,1,1))

            call s_bcycle(idwbc,idebc,idsbc,idnbc,                      &
     &                    2,1,ni-2,ni-1,2,1,nj-2,nj-1,                  &
     &                    ni,nj,nk,nwtrf(0,0,1,2))

          end if

          if(abs(cphopt).ge.2) then

            call s_bcycle(idwbc,idebc,idsbc,idnbc,                      &
     &                    2,1,ni-2,ni-1,2,1,nj-2,nj-1,                  &
     &                    ni,nj,nk,qicef(0,0,1,1))

            call s_bcycle(idwbc,idebc,idsbc,idnbc,                      &
     &                    2,1,ni-2,ni-1,2,1,nj-2,nj-1,                  &
     &                    ni,nj,nk,qicef(0,0,1,2))

            call s_bcycle(idwbc,idebc,idsbc,idnbc,                      &
     &                    2,1,ni-2,ni-1,2,1,nj-2,nj-1,                  &
     &                    ni,nj,nk,qicef(0,0,1,3))

            if(haiopt.eq.1) then

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,qicef(0,0,1,4))

            end if

          end if

          if(abs(cphopt).eq.2) then

            call s_bcycle(idwbc,idebc,idsbc,idnbc,                      &
     &                    2,1,ni-2,ni-1,2,1,nj-2,nj-1,                  &
     &                    ni,nj,nk,nicef(0,0,1,1))

          else if(abs(cphopt).ge.3) then

            call s_bcycle(idwbc,idebc,idsbc,idnbc,                      &
     &                    2,1,ni-2,ni-1,2,1,nj-2,nj-1,                  &
     &                    ni,nj,nk,nicef(0,0,1,1))

            call s_bcycle(idwbc,idebc,idsbc,idnbc,                      &
     &                    2,1,ni-2,ni-1,2,1,nj-2,nj-1,                  &
     &                    ni,nj,nk,nicef(0,0,1,2))

            call s_bcycle(idwbc,idebc,idsbc,idnbc,                      &
     &                    2,1,ni-2,ni-1,2,1,nj-2,nj-1,                  &
     &                    ni,nj,nk,nicef(0,0,1,3))

            if(haiopt.eq.1) then

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,nicef(0,0,1,4))

            end if

          end if

          if(cphopt.lt.0) then

            if(qcgopt.eq.2) then

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,qcwtrf(0,0,1,1))

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,qcwtrf(0,0,1,2))

            end if

            call s_bcycle(idwbc,idebc,idsbc,idnbc,                      &
     &                    2,1,ni-2,ni-1,2,1,nj-2,nj-1,                  &
     &                    ni,nj,nk,qcicef(0,0,1,1))

            call s_bcycle(idwbc,idebc,idsbc,idnbc,                      &
     &                    2,1,ni-2,ni-1,2,1,nj-2,nj-1,                  &
     &                    ni,nj,nk,qcicef(0,0,1,2))

            call s_bcycle(idwbc,idebc,idsbc,idnbc,                      &
     &                    2,1,ni-2,ni-1,2,1,nj-2,nj-1,                  &
     &                    ni,nj,nk,qcicef(0,0,1,3))

            if(haiopt.eq.1) then

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,qcicef(0,0,1,4))

            end if

          end if

        else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

          if(abs(cphopt).ge.11) then

            do n=1,nqw

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,qwtrf(0,0,1,n))

            end do

            do n=1,nnw

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,nwtrf(0,0,1,n))

            end do

          end if

          if(abs(cphopt).eq.12) then

            do n=1,nqi

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,qicef(0,0,1,n))

            end do

            do n=1,nni

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,nicef(0,0,1,n))

            end do

          end if

        end if

      end if

      if(aslopt.ge.1) then

        do n=1,nqa(0)

          call s_bcycle(idwbc,idebc,idsbc,idnbc,                        &
     &              2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,qaslf(0,0,1,n))

        end do

      end if

      if(trkopt.ge.1) then

        call bcycle(idwbc,idebc,idsbc,idnbc,                            &
     &              2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,qtf)

      end if

      if(tubopt.ge.2) then

        call bcycle(idwbc,idebc,idsbc,idnbc,                            &
     &              2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,tkef)

      end if

! -----

! Set the boundary conditions at the four corners.

      if(gwmopt.eq.0) then

        call bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,             &
     &               ni,nj,nk,ptpf)

      end if

      if(fmois(1:5).eq.'moist') then

        call bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,             &
     &               ni,nj,nk,qvf)

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,       &
     &                     ni,nj,nk,qwtrf(0,0,1,1))

            call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,       &
     &                     ni,nj,nk,qwtrf(0,0,1,2))

          end if

          if(abs(cphopt).eq.4) then

            call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,       &
     &                     ni,nj,nk,nwtrf(0,0,1,1))

            call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,       &
     &                     ni,nj,nk,nwtrf(0,0,1,2))

          end if

          if(abs(cphopt).ge.2) then

            call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,       &
     &                     ni,nj,nk,qicef(0,0,1,1))

            call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,       &
     &                     ni,nj,nk,qicef(0,0,1,2))

            call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,       &
     &                     ni,nj,nk,qicef(0,0,1,3))

            if(haiopt.eq.1) then

              call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,     &
     &                       ni,nj,nk,qicef(0,0,1,4))

            end if

          end if

          if(abs(cphopt).eq.2) then

            call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,       &
     &                     ni,nj,nk,nicef(0,0,1,1))

          else if(abs(cphopt).ge.3) then

            call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,       &
     &                     ni,nj,nk,nicef(0,0,1,1))

            call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,       &
     &                     ni,nj,nk,nicef(0,0,1,2))

            call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,       &
     &                     ni,nj,nk,nicef(0,0,1,3))

            if(haiopt.eq.1) then

              call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,     &
     &                       ni,nj,nk,nicef(0,0,1,4))

            end if

          end if

          if(cphopt.lt.0) then

            if(qcgopt.eq.2) then

              call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,     &
     &                       ni,nj,nk,qcwtrf(0,0,1,1))

              call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,     &
     &                       ni,nj,nk,qcwtrf(0,0,1,2))

            end if

            call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,       &
     &                     ni,nj,nk,qcicef(0,0,1,1))

            call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,       &
     &                     ni,nj,nk,qcicef(0,0,1,2))

            call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,       &
     &                     ni,nj,nk,qcicef(0,0,1,3))

            if(haiopt.eq.1) then

              call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,     &
     &                       ni,nj,nk,qcicef(0,0,1,4))

            end if

          end if

        else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

          if(abs(cphopt).ge.11) then

            do n=1,nqw

              call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,     &
     &                       ni,nj,nk,qwtrf(0,0,1,n))

            end do

            do n=1,nnw

              call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,     &
     &                       ni,nj,nk,nwtrf(0,0,1,n))

            end do

          end if

          if(abs(cphopt).eq.12) then

            do n=1,nqi

              call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,     &
     &                       ni,nj,nk,qicef(0,0,1,n))

            end do

            do n=1,nni

              call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,     &
     &                       ni,nj,nk,nicef(0,0,1,n))

            end do

          end if

        end if

      end if

      if(aslopt.ge.1) then

        do n=1,nqa(0)

          call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,         &
     &                   ni,nj,nk,qaslf(0,0,1,n))

        end do

      end if

      if(trkopt.ge.1) then

        call bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,             &
     &               ni,nj,nk,qtf)

      end if

      if(tubopt.ge.2) then

        call bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,             &
     &               ni,nj,nk,tkef)

      end if

! -----

! Set the bottom and the top boundary conditions.

      if(gwmopt.eq.0) then

        call vbcs(ni,nj,nk,ptpf)

      end if

      if(fmois(1:5).eq.'moist') then

        call vbcs(ni,nj,nk,qvf)

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            call s_vbcs(ni,nj,nk,qwtrf(0,0,1,1))
            call s_vbcs(ni,nj,nk,qwtrf(0,0,1,2))

          end if

          if(abs(cphopt).eq.4) then

            call s_vbcs(ni,nj,nk,nwtrf(0,0,1,1))
            call s_vbcs(ni,nj,nk,nwtrf(0,0,1,2))

          end if

          if(abs(cphopt).ge.2) then

            call s_vbcs(ni,nj,nk,qicef(0,0,1,1))
            call s_vbcs(ni,nj,nk,qicef(0,0,1,2))
            call s_vbcs(ni,nj,nk,qicef(0,0,1,3))

            if(haiopt.eq.1) then

              call s_vbcs(ni,nj,nk,qicef(0,0,1,4))

            end if

          end if

          if(abs(cphopt).eq.2) then

            call s_vbcs(ni,nj,nk,nicef(0,0,1,1))

          else if(abs(cphopt).ge.3) then

            call s_vbcs(ni,nj,nk,nicef(0,0,1,1))
            call s_vbcs(ni,nj,nk,nicef(0,0,1,2))
            call s_vbcs(ni,nj,nk,nicef(0,0,1,3))

            if(haiopt.eq.1) then

              call s_vbcs(ni,nj,nk,nicef(0,0,1,4))

            end if

          end if

          if(cphopt.lt.0) then

            if(qcgopt.eq.2) then

              call s_vbcqcg(ni,nj,nk,qcwtrf(0,0,1,1))
              call s_vbcqcg(ni,nj,nk,qcwtrf(0,0,1,2))

            end if

            call s_vbcqcg(ni,nj,nk,qcicef(0,0,1,1))
            call s_vbcqcg(ni,nj,nk,qcicef(0,0,1,2))
            call s_vbcqcg(ni,nj,nk,qcicef(0,0,1,3))

            if(haiopt.eq.1) then

              call s_vbcqcg(ni,nj,nk,qcicef(0,0,1,4))

            end if

          end if

        else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

          if(abs(cphopt).ge.11) then

            do n=1,nqw

              call s_vbcs(ni,nj,nk,qwtrf(0,0,1,n))

            end do

            do n=1,nnw

              call s_vbcs(ni,nj,nk,nwtrf(0,0,1,n))

            end do

          end if

          if(abs(cphopt).eq.12) then

            do n=1,nqi

              call s_vbcs(ni,nj,nk,qicef(0,0,1,n))

            end do

            do n=1,nni

              call s_vbcs(ni,nj,nk,nicef(0,0,1,n))

            end do

          end if

        end if

      end if

      if(aslopt.ge.1) then

        do n=1,nqa(0)

          call s_vbcs(ni,nj,nk,qaslf(0,0,1,n))

        end do

      end if

      if(trkopt.ge.1) then

        call vbcs(ni,nj,nk,qtf)

      end if

      if(tubopt.ge.2) then

        call vbcs(ni,nj,nk,tkef)

      end if

! -----

      end subroutine s_steps

!-----7--------------------------------------------------------------7--

      end module m_steps
