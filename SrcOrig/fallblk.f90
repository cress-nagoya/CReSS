!***********************************************************************
      module m_fallblk
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2000/10/27, 2000/11/17, 2001/04/15,
!                   2001/05/29, 2001/06/29, 2001/12/11, 2002/01/15,
!                   2002/04/02, 2002/12/02, 2003/03/28, 2003/04/30,
!                   2003/05/19, 2003/10/31, 2003/11/28, 2003/12/12,
!                   2004/03/05, 2004/03/22, 2004/04/01, 2004/04/15,
!                   2004/05/07, 2004/05/31, 2004/06/10, 2004/08/01,
!                   2004/08/20, 2004/09/01, 2004/09/10, 2004/09/25,
!                   2004/10/12, 2004/12/17, 2005/01/31, 2006/01/10,
!                   2006/02/13, 2006/04/03, 2006/05/12, 2006/07/21,
!                   2006/09/30, 2007/05/14, 2007/10/19, 2007/11/26,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2009/11/05,
!                   2011/01/14, 2011/06/01, 2011/09/22, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform the fall out.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkfall
      use m_comindx
      use m_commath
      use m_getiname
      use m_getrname
      use m_temparam
      use m_upwnp
      use m_upwqcg
      use m_upwqp

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: fallblk, s_fallblk

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface fallblk

        module procedure s_fallblk

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic int
      intrinsic min
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_fallblk(fpcphopt,fphaiopt,fpqcgopt,fpdz,             &
     &                     dtb,ni,nj,nk,jcb,rbr,rst,ucq,urq,            &
     &                     uiq,usq,ugq,uhq,ucn,urn,uin,usn,ugn,uhn,     &
     &                     qcf,qrf,qif,qsf,qgf,qhf,nccf,ncrf,ncif,ncsf, &
     &                     ncgf,nchf,qccf,qrcf,qicf,qscf,qgcf,qhcf,     &
     &                     prc,prr,pri,prs,prg,prh,tmp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fphaiopt
                       ! Formal parameter of unique index of haiopt

      integer, intent(in) :: fpqcgopt
                       ! Formal parameter of unique index of qcgopt

      integer, intent(in) :: fpdz
                       ! Formal parameter of unique index of dz

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: ucq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of cloud water

      real, intent(in) :: urq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of rain water

      real, intent(in) :: uiq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of cloud ice

      real, intent(in) :: usq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of snow

      real, intent(in) :: ugq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of graupel

      real, intent(in) :: uhq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of hail

      real, intent(in) :: ucn(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of cloud water concentrations

      real, intent(in) :: urn(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of rain water concentrations

      real, intent(in) :: uin(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of cloud ice concentrations

      real, intent(in) :: usn(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of snow concentrations

      real, intent(in) :: ugn(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of graupel concentrations

      real, intent(in) :: uhn(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of hail concentrations

! Input and output variables

      real, intent(inout) :: qcf(0:ni+1,0:nj+1,1:nk)
                       ! Cloud water mixnig ratio at future

      real, intent(inout) :: qrf(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixnig ratio at future

      real, intent(inout) :: qif(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixnig ratio at future

      real, intent(inout) :: qsf(0:ni+1,0:nj+1,1:nk)
                       ! Snow mixnig ratio at future

      real, intent(inout) :: qgf(0:ni+1,0:nj+1,1:nk)
                       ! Graupel mixnig ratio at future

      real, intent(inout) :: qhf(0:ni+1,0:nj+1,1:nk)
                       ! Hail mixnig ratio at future

      real, intent(inout) :: nccf(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud water at future

      real, intent(inout) :: ncrf(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of rain water at future

      real, intent(inout) :: ncif(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud ice at future

      real, intent(inout) :: ncsf(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of snow at future

      real, intent(inout) :: ncgf(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of graupel at future

      real, intent(inout) :: nchf(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of hail at future

      real, intent(inout) :: qccf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of cloud water at future

      real, intent(inout) :: qrcf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of rain water at future

      real, intent(inout) :: qicf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of cloud ice at future

      real, intent(inout) :: qscf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of snow at future

      real, intent(inout) :: qgcf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of graupel at future

      real, intent(inout) :: qhcf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of hail at future

      real, intent(inout) :: prc(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation for cloud water

      real, intent(inout) :: prr(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation for rain

      real, intent(inout) :: pri(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation for cloud ice

      real, intent(inout) :: prs(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation for snow

      real, intent(inout) :: prg(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation for graupel

      real, intent(inout) :: prh(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation for hail

! Internal shared variables

      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes
      integer qcgopt   ! Option for charging distribution

      integer npstp    ! Number of steps for fall out time integration

      integer ipstp    ! Index of steps for fall out time integration

      real dz          ! Grid distance in z direction

      real dtp         ! Time steps interval of fall out integration

      real dtpc        ! Time steps interval of fall out integration
                       ! for cloud water

      real dtpr        ! Time steps interval of fall out integration
                       ! for rain water

      real dtpi        ! Time steps interval of fall out integration
                       ! for cloud ice

      real dtps        ! Time steps interval of fall out integration
                       ! for snow

      real dtpg        ! Time steps interval of fall out integration
                       ! for graupel

      real dtph        ! Time steps interval of fall out integration
                       ! for hail

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real dzjcb       ! dz x jcb

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)
      call getiname(fpqcgopt,qcgopt)
      call getrname(fpdz,dz)

! -----

!! Perform the fall out.

! Initialize the processed variables.

      dtpc=dtb
      dtpr=dtb
      dtpi=dtb
      dtps=dtb
      dtpg=dtb
      dtph=dtb

! -----

! Get the minimum time interval.

!$omp parallel default(shared) private(k)

      if(flqcqi_opt.eq.0) then

        if(haiopt.eq.0) then

          do k=1,nk-1

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,dzjcb) reduction(min: dtpr,dtps,dtpg)

            do j=1,nj-1
            do i=1,ni-1
              dzjcb=dz*jcb(i,j,k)

              dtpr=min(dzjcb/(urq(i,j,k)+eps),dtpr)
              dtps=min(dzjcb/(usq(i,j,k)+eps),dtps)
              dtpg=min(dzjcb/(ugq(i,j,k)+eps),dtpg)

            end do
            end do

!$omp end do

          end do

        else

          do k=1,nk-1

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,dzjcb) reduction(min: dtpr,dtps,dtpg,dtph)

            do j=1,nj-1
            do i=1,ni-1
              dzjcb=dz*jcb(i,j,k)

              dtpr=min(dzjcb/(urq(i,j,k)+eps),dtpr)
              dtps=min(dzjcb/(usq(i,j,k)+eps),dtps)
              dtpg=min(dzjcb/(ugq(i,j,k)+eps),dtpg)
              dtph=min(dzjcb/(uhq(i,j,k)+eps),dtph)

            end do
            end do

!$omp end do

          end do

        end if

      else

        if(haiopt.eq.0) then

          do k=1,nk-1

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,dzjcb) reduction(min: dtpc,dtpr,dtpi,dtps,dtpg)

            do j=1,nj-1
            do i=1,ni-1
              dzjcb=dz*jcb(i,j,k)

              dtpc=min(dzjcb/(ucq(i,j,k)+eps),dtpc)
              dtpr=min(dzjcb/(urq(i,j,k)+eps),dtpr)
              dtpi=min(dzjcb/(uiq(i,j,k)+eps),dtpi)
              dtps=min(dzjcb/(usq(i,j,k)+eps),dtps)
              dtpg=min(dzjcb/(ugq(i,j,k)+eps),dtpg)

            end do
            end do

!$omp end do

          end do

        else

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j,dzjcb)                           &
!$omp&   reduction(min: dtpc,dtpr,dtpi,dtps,dtpg,dtph)

            do j=1,nj-1
            do i=1,ni-1
              dzjcb=dz*jcb(i,j,k)

              dtpc=min(dzjcb/(ucq(i,j,k)+eps),dtpc)
              dtpr=min(dzjcb/(urq(i,j,k)+eps),dtpr)
              dtpi=min(dzjcb/(uiq(i,j,k)+eps),dtpi)
              dtps=min(dzjcb/(usq(i,j,k)+eps),dtps)
              dtpg=min(dzjcb/(ugq(i,j,k)+eps),dtpg)
              dtph=min(dzjcb/(uhq(i,j,k)+eps),dtph)

            end do
            end do

!$omp end do

          end do

        end if

      end if

!$omp end parallel

! -----

! Get the time interval and number of steps.

      if(flqcqi_opt.eq.0) then

        if(haiopt.eq.0) then
          dtp=min(dtpr,dtps,dtpg)
        else
          dtp=min(dtpr,dtps,dtpg,dtph)
        end if

      else

        if(haiopt.eq.0) then
          dtp=min(dtpc,dtpr,dtpi,dtps,dtpg)
        else
          dtp=min(dtpc,dtpr,dtpi,dtps,dtpg,dtph)
        end if

      end if

      call chkfall(dtp)

      if(dtp.lt.dtb) then
        npstp=int((dtb+.001e0)/dtp)+1
        dtp=dtb/real(npstp)

      else
        npstp=1
        dtp=dtb

      end if

! -----

! Calculate the sedimentation and precipitation.

      if(flqcqi_opt.ne.0) then

        do ipstp=1,npstp

          call upwqp(idadvopt,iddziv,dtp,ni,nj,nk,rbr,rst,ucq,qcf,prc,  &
     &               tmp1)

          call upwqp(idadvopt,iddziv,dtp,ni,nj,nk,rbr,rst,uiq,qif,pri,  &
     &               tmp1)

          if(abs(cphopt).eq.4) then

            call upwnp(iddziv,dtp,ni,nj,nk,rbr,rst,ucn,nccf,tmp1)

          end if

          call upwnp(iddziv,dtp,ni,nj,nk,rbr,rst,uin,ncif,tmp1)

          if(cphopt.lt.0) then

            if(qcgopt.eq.2) then

              call upwqcg(iddziv,dtp,ni,nj,nk,rbr,rst,ucq,qccf,tmp1)

            end if

            call upwqcg(iddziv,dtp,ni,nj,nk,rbr,rst,uiq,qicf,tmp1)

          end if

        end do

      end if

      do ipstp=1,npstp

        call upwqp(idadvopt,iddziv,dtp,ni,nj,nk,rbr,rst,urq,qrf,prr,    &
     &             tmp1)

        call upwqp(idadvopt,iddziv,dtp,ni,nj,nk,rbr,rst,usq,qsf,prs,    &
     &             tmp1)

        call upwqp(idadvopt,iddziv,dtp,ni,nj,nk,rbr,rst,ugq,qgf,prg,    &
     &             tmp1)

        if(haiopt.eq.1) then

          call upwqp(idadvopt,iddziv,dtp,ni,nj,nk,rbr,rst,uhq,qhf,prh,  &
     &               tmp1)

        end if

        if(abs(cphopt).eq.4) then

          call upwnp(iddziv,dtp,ni,nj,nk,rbr,rst,urn,ncrf,tmp1)

        end if

        if(abs(cphopt).ge.3) then

          call upwnp(iddziv,dtp,ni,nj,nk,rbr,rst,usn,ncsf,tmp1)
          call upwnp(iddziv,dtp,ni,nj,nk,rbr,rst,ugn,ncgf,tmp1)

          if(haiopt.eq.1) then

            call upwnp(iddziv,dtp,ni,nj,nk,rbr,rst,uhn,nchf,tmp1)

          end if

        end if

        if(cphopt.lt.0) then

          if(qcgopt.eq.2) then

            call upwqcg(iddziv,dtp,ni,nj,nk,rbr,rst,urq,qrcf,tmp1)

          end if

          call upwqcg(iddziv,dtp,ni,nj,nk,rbr,rst,usq,qscf,tmp1)
          call upwqcg(iddziv,dtp,ni,nj,nk,rbr,rst,ugq,qgcf,tmp1)

          if(haiopt.eq.1) then

            call upwqcg(iddziv,dtp,ni,nj,nk,rbr,rst,uhq,qhcf,tmp1)

          end if

        end if

      end do

! -----

!! -----

      end subroutine s_fallblk

!-----7--------------------------------------------------------------7--

      end module m_fallblk
