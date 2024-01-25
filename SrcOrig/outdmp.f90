!***********************************************************************
      module m_outdmp
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/04/18
!     Modification: 2000/06/01, 2000/07/05, 2000/12/18, 2001/01/15,
!                   2001/03/13, 2001/04/15, 2001/05/29, 2001/06/29,
!                   2001/12/11, 2002/01/07, 2002/04/02, 2002/06/18,
!                   2002/08/15, 2002/12/02, 2003/03/28, 2003/04/30,
!                   2003/05/19, 2003/09/01, 2003/11/28, 2003/12/12,
!                   2004/02/01, 2004/03/05, 2004/05/31, 2004/06/10,
!                   2004/08/01, 2004/08/20, 2004/09/25, 2006/01/10,
!                   2006/02/13, 2006/07/21, 2006/09/11, 2006/09/21,
!                   2006/11/06, 2007/01/05, 2007/01/20, 2007/01/31,
!                   2007/05/14, 2007/05/21, 2007/07/30, 2007/10/19,
!                   2007/11/26, 2008/01/11, 2008/04/17, 2008/05/02,
!                   2008/07/01, 2008/08/25, 2008/10/10, 2009/01/30,
!                   2009/02/27, 2011/06/01, 2011/08/18, 2011/09/22,
!                   2011/11/10, 2013/01/28, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the inferior procedures to read in the variables to the
!     dumped file.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comblk
      use m_comcapt
      use m_comdmp
      use m_comindx
      use m_commath
      use m_getcname
      use m_getiname
      use m_getrname
      use m_inichar
      use m_outdmp2d
      use m_outdmp3d
      use m_rotuvm2s
      use m_setcst2d
      use m_setcst3d
      use m_setproj
      use m_temparam
      use m_totalq
      use m_totals
      use m_var8u8s
      use m_var8v8s
      use m_var8w8s
      use m_vint31s

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: outdmp, s_outdmp

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outdmp

        module procedure s_outdmp

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic mod
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_outdmp(fpdmpvar,fpcphopt,fphaiopt,fpqcgopt,          &
     &                    fpaslopt,fptrkopt,fptubopt,fpdmplev,fpdz,     &
     &                    fmois,ni,nj,nk,nqw,nnw,nqi,nni,nqa,zph,zsth,  &
     &                    lon,ubr,vbr,pbr,ptbr,qvbr,u,v,w,pp,ptp,qv,    &
     &                    qwtr,nwtr,qice,nice,qcwtr,qcice,qasl,qt,tke,  &
     &                    maxvl,prwtr,price,z1d,zph8s,var,              &
     &                    tmp1,tmp2,tmp3,tmp4)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: fpdmpvar
                       ! Formal parameter of unique index of dmpvar

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

      integer, intent(in) :: fpdmplev
                       ! Formal parameter of unique index of dmplev

      integer, intent(in) :: fpdz
                       ! Formal parameter of unique index of dz

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

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(in) :: zsth(1:nk)
                       ! 1 dimensional stretched z coordinates

      real, intent(in) :: lon(0:ni+1,0:nj+1)
                       ! Longitude

      real, intent(in) :: ubr(0:ni+1,0:nj+1,1:nk)
                       ! Base state x components of velocity

      real, intent(in) :: vbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state y components of velocity

      real, intent(in) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: qvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state water vapor mixing ratio

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

      real, intent(in) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation

      real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

      real, intent(in) :: qwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor

      real, intent(in) :: nwtr(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations

      real, intent(in) :: qice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor

      real, intent(in) :: nice(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations

      real, intent(in) :: qcwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water

      real, intent(in) :: qcice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice

      real, intent(in) :: qasl(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio

      real, intent(in) :: qt(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio

      real, intent(in) :: tke(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy

! Input and output variables

      real, intent(inout) :: maxvl(0:ni+1,0:nj+1,1:nk)
                       ! Maximum instantaneous wind velocity

      real, intent(inout) :: prwtr(0:ni+1,0:nj+1,1:2,1:nqw)
                       ! Precipitation and accumulation for water

      real, intent(inout) :: price(0:ni+1,0:nj+1,1:2,1:nqi)
                       ! Precipitation and accumulation for ice

! Internal shared variables

      character(len=108) dmpvar
                       ! Control flag of dumped variables

      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes
      integer qcgopt   ! Option for charging distribution
      integer aslopt   ! Option for aerosol processes
      integer trkopt   ! Option for mixing ratio tracking
      integer tubopt   ! Option for turbulent mixing

      integer dmplev   ! Option for z coordinates of dumped variables

      integer nqw_sub  ! Substitute for nqw
      integer nnw_sub  ! Substitute for nnw
      integer nqi_sub  ! Substitute for nqi
      integer nni_sub  ! Substitute for nni

      integer istr     ! Start index of types of aerosol array
      integer iend     ! End index of types of aerosol array

      real dz          ! Grid distance in z direction

      real cpj(1:7)    ! Map projection parameters

      real, intent(inout) :: z1d(1:nk)
                       ! Constant height at scalar points

      real, intent(inout) :: zph8s(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates at scalar points

      real, intent(inout) :: var(0:ni+1,0:nj+1,1:nk)
                       ! Temporary output variable

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp4(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variable

      integer k        ! Array index in z drection

!-----7--------------------------------------------------------------7--

!!!! Control the inferior procedures to read in the variables to the
!!!! dumped file.

      if(fdmp(1:3).eq.'act'.or.fmon(1:3).eq.'act') then

! Initialize the character variable.

        call inichar(dmpvar)

! -----

! Get the required namelist variables.

        call getcname(fpdmpvar,dmpvar)
        call getiname(fpcphopt,cphopt)
        call getiname(fphaiopt,haiopt)
        call getiname(fpqcgopt,qcgopt)
        call getiname(fpaslopt,aslopt)
        call getiname(fptrkopt,trkopt)
        call getiname(fptubopt,tubopt)
        call getiname(fpdmplev,dmplev)
        call getrname(fpdz,dz)

! -----

! Set the substituted variables.

        nqw_sub=nqw
        nnw_sub=nnw
        nqi_sub=nqi
        nni_sub=nni

! -----

! Calculate the constant height.

!$omp parallel default(shared)

        if(fdmp(1:3).eq.'act') then

          if(mod(dmplev,10).eq.2) then

!$omp do schedule(runtime) private(k)

            do k=2,nk-2
              z1d(k)=dz*(real(k)-1.5e0)
            end do

!$omp end do

          else if(mod(dmplev,10).eq.3) then

!$omp do schedule(runtime) private(k)

            do k=2,nk-2
              z1d(k)=.5e0*(zsth(k)+zsth(k+1))
            end do

!$omp end do

          end if

        end if

!$omp end parallel

! -----

! Get the z physical coordinates at scalar points.

        if(fdmp(1:3).eq.'act') then

          if(mod(dmplev,10).eq.2.or.mod(dmplev,10).eq.3) then

            call var8w8s(ni,nj,nk,zph,zph8s)

          end if

        end if

! -----

!!! Read in the 3 dimensional data to the dumped file.

        if(fdmp(1:3).eq.'act') then

!! Read in the data to the dumped file to terrain following plane.

          if(mod(dmplev,10).eq.1) then

! Dump the x and y components of velocity.

            if(dmplev.eq.1) then

              if(dmpvar(1:1).eq.'o') then

                call outdmp3d('ubar  ',4,capt(28),ncpt(28),'oxx',       &
     &                        ni,nj,nk,ubr)

                call outdmp3d('u     ',1,capt(1),ncpt(1),'oxx',         &
     &                        ni,nj,nk,u)

              else if(dmpvar(1:1).eq.'-') then

                call outdmp3d('u     ',1,capt(1),ncpt(1),'oxx',         &
     &                        ni,nj,nk,u)

              end if

              if(dmpvar(2:2).eq.'o') then

                call outdmp3d('vbar  ',4,capt(29),ncpt(29),'xox',       &
     &                        ni,nj,nk,vbr)

                call outdmp3d('v     ',1,capt(2),ncpt(2),'xox',         &
     &                        ni,nj,nk,v)

              else if(dmpvar(2:2).eq.'-') then

                call outdmp3d('v     ',1,capt(2),ncpt(2),'xox',         &
     &                        ni,nj,nk,v)

              end if

! -----

! Dump the zonal and meridional velocity.

            else

              call setproj(idmpopt,idnspol,idtlat1,idtlat2,cpj)

              if(dmpvar(1:1).eq.'o'.or.dmpvar(2:2).eq.'o') then

                call var8u8s(ni,nj,nk,ubr,tmp1)
                call var8v8s(ni,nj,nk,vbr,tmp2)

                call var8u8s(ni,nj,nk,u,tmp3)
                call var8v8s(ni,nj,nk,v,tmp4)

                call rotuvm2s(idmpopt,idnspol,idtlon,                   &
     &                        2,ni-2,2,nj-2,2,nk-2,cpj,                 &
     &                        0,ni+1,0,nj+1,1,nk,lon,tmp1,tmp2)

                call rotuvm2s(idmpopt,idnspol,idtlon,                   &
     &                        2,ni-2,2,nj-2,2,nk-2,cpj,                 &
     &                        0,ni+1,0,nj+1,1,nk,lon,tmp3,tmp4)

                if(dmpvar(1:1).eq.'o') then

                  call outdmp3d('ubar  ',4,capt(32),ncpt(32),'xxx',     &
     &                          ni,nj,nk,tmp1)

                  call outdmp3d('u     ',1,capt(30),ncpt(30),'xxx',     &
     &                          ni,nj,nk,tmp3)

                else if(dmpvar(1:1).eq.'-') then

                  call outdmp3d('u     ',1,capt(30),ncpt(30),'xxx',     &
     &                          ni,nj,nk,tmp3)

                end if

                if(dmpvar(2:2).eq.'o') then

                  call outdmp3d('vbar  ',4,capt(33),ncpt(33),'xxx',     &
     &                          ni,nj,nk,tmp2)

                  call outdmp3d('v     ',1,capt(31),ncpt(31),'xxx',     &
     &                          ni,nj,nk,tmp4)

                else if(dmpvar(2:2).eq.'-') then

                  call outdmp3d('v     ',1,capt(31),ncpt(31),'xxx',     &
     &                          ni,nj,nk,tmp4)

                end if

              else if(dmpvar(1:1).ne.'x'.and.dmpvar(2:2).ne.'x') then

                call var8u8s(ni,nj,nk,u,tmp1)
                call var8v8s(ni,nj,nk,v,tmp2)

                call rotuvm2s(idmpopt,idnspol,idtlon,                   &
     &                        2,ni-2,2,nj-2,2,nk-2,cpj,                 &
     &                        0,ni+1,0,nj+1,1,nk,lon,tmp1,tmp2)

                if(dmpvar(1:1).eq.'-') then

                  call outdmp3d('u     ',1,capt(30),ncpt(30),'xxx',     &
     &                          ni,nj,nk,tmp1)

                end if

                if(dmpvar(2:2).eq.'-') then

                  call outdmp3d('v     ',1,capt(31),ncpt(31),'xxx',     &
     &                          ni,nj,nk,tmp2)

                end if

              end if

            end if

! -----

! Dump the z components of velocity.

            if(dmpvar(3:3).eq.'o') then

              call outdmp3d('w     ',1,capt(3),ncpt(3),'xxo',           &
     &                      ni,nj,nk,w)

            end if

! -----

! Dump the base state pressure and its perturbation.

            if(dmpvar(4:4).eq.'o') then

              call outdmp3d('pbar  ',4,capt(4),ncpt(4),'xxx',           &
     &                      ni,nj,nk,pbr)

              call outdmp3d('pp    ',2,capt(5),ncpt(5),'xxx',           &
     &                      ni,nj,nk,pp)

            else if(dmpvar(4:4).eq.'-') then

              call totals(ni,nj,nk,pbr,pp,var)

              call outdmp3d('p     ',1,capt(26),ncpt(26),'xxx',         &
     &                      ni,nj,nk,var)

            end if

! -----

! Dump the base state potential temperature and its perturbation.

            if(dmpvar(5:5).eq.'o') then

              call outdmp3d('ptbar ',5,capt(6),ncpt(6),'xxx',           &
     &                      ni,nj,nk,ptbr)

              call outdmp3d('ptp   ',3,capt(7),ncpt(7),'xxx',           &
     &                      ni,nj,nk,ptp)

            else if(dmpvar(5:5).eq.'-') then

              call totals(ni,nj,nk,ptbr,ptp,var)

              call outdmp3d('pt    ',2,capt(27),ncpt(27),'xxx',         &
     &                      ni,nj,nk,var)

            end if

! -----

! Dump the water vapor mixing ratio.

            if(fmois(1:5).eq.'moist') then

              if(dmpvar(6:6).eq.'o') then

                call outdmp3d('qvbar ',5,capt(8),ncpt(8),'xxx',         &
     &                        ni,nj,nk,qvbr)

                call outdmp3d('qv    ',2,capt(9),ncpt(9),'xxx',         &
     &                        ni,nj,nk,qv)

              else if(dmpvar(6:6).eq.'-') then

                call outdmp3d('qv    ',2,capt(9),ncpt(9),'xxx',         &
     &                        ni,nj,nk,qv)

              end if

            end if

! -----

! Dump the water hydrometeor mixing ratio and the ice hydrometeor mixing
! ratio and concentrations.

            if(fmois(1:5).eq.'moist') then

              if(abs(cphopt).lt.10) then

                if(abs(cphopt).ge.1.and.dmpvar(7:7).eq.'o') then

                  call s_outdmp3d('qc    ',2,capt(10),ncpt(10),'xxx',   &
     &                            ni,nj,nk,qwtr(0,0,1,1))

                  call s_outdmp3d('qr    ',2,capt(11),ncpt(11),'xxx',   &
     &                            ni,nj,nk,qwtr(0,0,1,2))

                end if

                if(abs(cphopt).ge.2.and.dmpvar(7:7).eq.'o') then

                  call s_outdmp3d('qi    ',2,capt(14),ncpt(14),'xxx',   &
     &                            ni,nj,nk,qice(0,0,1,1))

                  call s_outdmp3d('qs    ',2,capt(15),ncpt(15),'xxx',   &
     &                            ni,nj,nk,qice(0,0,1,2))

                  call s_outdmp3d('qg    ',2,capt(16),ncpt(16),'xxx',   &
     &                            ni,nj,nk,qice(0,0,1,3))

                  if(haiopt.eq.1) then

                    call s_outdmp3d('qh    ',2,capt(99),ncpt(99),       &
     &                              'xxx',ni,nj,nk,qice(0,0,1,4))

                  end if

                end if

                if(abs(cphopt).eq.4.and.dmpvar(8:8).eq.'o') then

                  call s_outdmp3d('ncc   ',3,capt(74),ncpt(74),'xxx',   &
     &                            ni,nj,nk,nwtr(0,0,1,1))

                  call s_outdmp3d('ncr   ',3,capt(75),ncpt(75),'xxx',   &
     &                            ni,nj,nk,nwtr(0,0,1,2))

                end if

                if(abs(cphopt).ge.3.and.dmpvar(8:8).eq.'o') then

                  call s_outdmp3d('nci   ',3,capt(21),ncpt(21),'xxx',   &
     &                            ni,nj,nk,nice(0,0,1,1))

                  call s_outdmp3d('ncs   ',3,capt(22),ncpt(22),'xxx',   &
     &                            ni,nj,nk,nice(0,0,1,2))

                  call s_outdmp3d('ncg   ',3,capt(23),ncpt(23),'xxx',   &
     &                            ni,nj,nk,nice(0,0,1,3))

                  if(haiopt.eq.1) then

                    call s_outdmp3d('nch   ',3,capt(102),ncpt(102),     &
     &                              'xxx',ni,nj,nk,nice(0,0,1,4))

                  end if

                end if

                if(abs(cphopt).ge.2.and.dmpvar(7:7).eq.'o') then

                  if(flout_opt.eq.1) then

                    if(flqcqi_opt.ne.0) then

                      call outdmp3d('ucq   ',3,capt(82),ncpt(82),       &
     &                              'xxx',ni,nj,nk,ucq)

                    end if

                    call outdmp3d('urq   ',3,capt(83),ncpt(83),'xxx',   &
     &                            ni,nj,nk,urq)

                    if(flqcqi_opt.ne.0) then

                      call outdmp3d('uiq   ',3,capt(84),ncpt(84),       &
     &                              'xxx',ni,nj,nk,uiq)

                    end if

                    call outdmp3d('usq   ',3,capt(85),ncpt(85),'xxx',   &
     &                            ni,nj,nk,usq)

                    call outdmp3d('ugq   ',3,capt(86),ncpt(86),'xxx',   &
     &                            ni,nj,nk,ugq)

                    if(haiopt.eq.1) then

                      call outdmp3d('uhq   ',3,capt(104),ncpt(104),     &
     &                              'xxx',ni,nj,nk,uhq)

                    end if

                  end if

                end if

                if(cphopt.lt.0.and.dmpvar(9:9).eq.'o') then

                  if(qcgopt.eq.2) then

                    call s_outdmp3d('qcc   ',3,capt(43),ncpt(43),       &
     &                              'xxx',ni,nj,nk,qcwtr(0,0,1,1))

                    call s_outdmp3d('qrc   ',3,capt(44),ncpt(44),       &
     &                              'xxx',ni,nj,nk,qcwtr(0,0,1,2))

                  end if

                  call s_outdmp3d('qic   ',3,capt(45),ncpt(45),'xxx',   &
     &                            ni,nj,nk,qcice(0,0,1,1))

                  call s_outdmp3d('qsc   ',3,capt(46),ncpt(46),'xxx',   &
     &                            ni,nj,nk,qcice(0,0,1,2))

                  call s_outdmp3d('qgc   ',3,capt(47),ncpt(47),'xxx',   &
     &                            ni,nj,nk,qcice(0,0,1,3))

                  if(haiopt.eq.1) then

                    call s_outdmp3d('qhc   ',3,capt(103),ncpt(103),     &
     &                              'xxx',ni,nj,nk,qcice(0,0,1,4))

                  end if

                end if

              else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

                if(abs(cphopt).ge.11) then

                  if(dmpvar(7:7).eq.'o') then

                    call totalq(1,nqw,ni,nj,nk,nqw_sub,qwtr,var)

                    call outdmp3d('qwtr  ',4,capt(35),ncpt(35),'xxx',   &
     &                            ni,nj,nk,var)

                  end if

                  if(dmpvar(8:8).eq.'o') then

                    call totalq(1,nnw,ni,nj,nk,nnw_sub,nwtr,var)

                    call outdmp3d('nwtr  ',4,capt(36),ncpt(36),'xxx',   &
     &                            ni,nj,nk,var)

                  end if

                end if

                if(abs(cphopt).eq.12) then

                  if(dmpvar(7:7).eq.'o') then

                    call totalq(1,nqi,ni,nj,nk,nqi_sub,qice,var)

                    call outdmp3d('qice  ',4,capt(39),ncpt(39),'xxx',   &
     &                            ni,nj,nk,var)

                  end if

                  if(dmpvar(8:8).eq.'o') then

                    call totalq(1,nni,ni,nj,nk,nni_sub,nice,var)

                    call outdmp3d('nice  ',4,capt(40),ncpt(40),'xxx',   &
     &                            ni,nj,nk,var)

                  end if

                end if

              end if

            end if

! -----

! Dump the aerosol mixing ratio.

            if(aslopt.ge.1.and.dmpvar(10:10).eq.'o') then

              istr=1
              iend=nqa(1)-nqw-nqi

              call s_totalq(istr,iend,ni,nj,nk,nqa(0),qasl,var)

              call outdmp3d('drydu ',5,capt(87),ncpt(87),'xxx',         &
     &                      ni,nj,nk,var)

              istr=iend+1
              iend=iend+nqw+nqi

              call s_totalq(istr,iend,ni,nj,nk,nqa(0),qasl,var)

              call outdmp3d('wetdu ',5,capt(91),ncpt(91),'xxx',         &
     &                      ni,nj,nk,var)

              istr=iend+1
              iend=iend+nqa(2)-nqw-nqi

              call s_totalq(istr,iend,ni,nj,nk,nqa(0),qasl,var)

              call outdmp3d('dryca ',5,capt(88),ncpt(88),'xxx',         &
     &                      ni,nj,nk,var)

              istr=iend+1
              iend=iend+nqw+nqi

              call s_totalq(istr,iend,ni,nj,nk,nqa(0),qasl,var)

              call outdmp3d('wetca ',5,capt(92),ncpt(92),'xxx',         &
     &                      ni,nj,nk,var)

              istr=iend+1
              iend=iend+nqa(3)-nqw-nqi

              call s_totalq(istr,iend,ni,nj,nk,nqa(0),qasl,var)

              call outdmp3d('drysu ',5,capt(89),ncpt(89),'xxx',         &
     &                      ni,nj,nk,var)

              istr=iend+1
              iend=iend+nqw+nqi

              call s_totalq(istr,iend,ni,nj,nk,nqa(0),qasl,var)

              call outdmp3d('wetsu ',5,capt(93),ncpt(93),'xxx',         &
     &                      ni,nj,nk,var)

              istr=iend+1
              iend=iend+nqa(4)-nqw-nqi

              call s_totalq(istr,iend,ni,nj,nk,nqa(0),qasl,var)

              call outdmp3d('drysa ',5,capt(90),ncpt(90),'xxx',         &
     &                      ni,nj,nk,var)

              istr=iend+1
              iend=iend+nqw+nqi

              call s_totalq(istr,iend,ni,nj,nk,nqa(0),qasl,var)

              call outdmp3d('wetsa ',5,capt(94),ncpt(94),'xxx',         &
     &                      ni,nj,nk,var)

            end if

! -----

! Dump the tracer mixing ratio.

            if(trkopt.ge.1.and.dmpvar(11:11).eq.'o') then

              call outdmp3d('qt    ',2,capt(34),ncpt(34),'xxx',         &
     &                      ni,nj,nk,qt)

            end if

! -----

! Dump the turbulent kinetic energy and the maximum instantaneous wind
! velocity.

            if(tubopt.ge.2) then

              if(dmpvar(12:12).eq.'+'.or.dmpvar(12:12).eq.'o') then

                call outdmp3d('tke   ',3,capt(25),ncpt(25),'xxx',       &
     &                        ni,nj,nk,tke)

              end if

              if(dmpvar(12:12).eq.'+'.or.dmpvar(12:12).eq.'-') then

                call outdmp3d('maxvl ',5,capt(73),ncpt(73),'xxx',       &
     &                        ni,nj,nk,maxvl)

                call setcst3d(0,ni+1,0,nj+1,1,nk,lim35n,maxvl)

              end if

            end if

! -----

!! -----

!! Read in the data to the dumped file to the constant height.

          else if(mod(dmplev,10).eq.2.or.mod(dmplev,10).eq.3) then

! Dump the x and y components of velocity.

            if(dmplev.lt.10) then

              if(dmpvar(1:1).eq.'o') then

                call var8u8s(ni,nj,nk,ubr,tmp1)

                call vint31s(ni,nj,nk,z1d,zph8s,tmp1,var)

                call outdmp3d('ubar  ',4,capt(28),ncpt(28),'xxx',       &
     &                        ni,nj,nk,var)

                call var8u8s(ni,nj,nk,u,tmp1)

                call vint31s(ni,nj,nk,z1d,zph8s,tmp1,var)

                call outdmp3d('u     ',1,capt(1),ncpt(1),'xxx',         &
     &                        ni,nj,nk,var)

              else if(dmpvar(1:1).eq.'-') then

                call var8u8s(ni,nj,nk,u,tmp1)

                call vint31s(ni,nj,nk,z1d,zph8s,tmp1,var)

                call outdmp3d('u     ',1,capt(1),ncpt(1),'xxx',         &
     &                        ni,nj,nk,var)

              end if

              if(dmpvar(2:2).eq.'o') then

                call var8v8s(ni,nj,nk,vbr,tmp1)

                call vint31s(ni,nj,nk,z1d,zph8s,tmp1,var)

                call outdmp3d('vbar  ',4,capt(29),ncpt(29),'xxx',       &
     &                        ni,nj,nk,var)

                call var8v8s(ni,nj,nk,v,tmp1)

                call vint31s(ni,nj,nk,z1d,zph8s,tmp1,var)

                call outdmp3d('v     ',1,capt(2),ncpt(2),'xxx',         &
     &                        ni,nj,nk,var)

              else if(dmpvar(2:2).eq.'-') then

                call var8v8s(ni,nj,nk,v,tmp1)

                call vint31s(ni,nj,nk,z1d,zph8s,tmp1,var)

                call outdmp3d('v     ',1,capt(2),ncpt(2),'xxx',         &
     &                        ni,nj,nk,var)

              end if

! -----

! Dump the zonal and meridional velocity.

            else

              call setproj(idmpopt,idnspol,idtlat1,idtlat2,cpj)

              if(dmpvar(1:1).eq.'o'.or.dmpvar(2:2).eq.'o') then

                call var8u8s(ni,nj,nk,ubr,tmp1)
                call var8v8s(ni,nj,nk,vbr,tmp2)

                call var8u8s(ni,nj,nk,u,tmp3)
                call var8v8s(ni,nj,nk,v,tmp4)

                call rotuvm2s(idmpopt,idnspol,idtlon,                   &
     &                        2,ni-2,2,nj-2,2,nk-2,cpj,                 &
     &                        0,ni+1,0,nj+1,1,nk,lon,tmp1,tmp2)

                call rotuvm2s(idmpopt,idnspol,idtlon,                   &
     &                        2,ni-2,2,nj-2,2,nk-2,cpj,                 &
     &                        0,ni+1,0,nj+1,1,nk,lon,tmp3,tmp4)

                if(dmpvar(1:1).eq.'o') then

                  call vint31s(ni,nj,nk,z1d,zph8s,tmp1,var)

                  call outdmp3d('ubar  ',4,capt(32),ncpt(32),'xxx',     &
     &                          ni,nj,nk,var)

                  call vint31s(ni,nj,nk,z1d,zph8s,tmp3,var)

                  call outdmp3d('u     ',1,capt(30),ncpt(30),'xxx',     &
     &                          ni,nj,nk,var)

                else if(dmpvar(1:1).eq.'-') then

                  call vint31s(ni,nj,nk,z1d,zph8s,tmp3,var)

                  call outdmp3d('u     ',1,capt(30),ncpt(30),'xxx',     &
     &                          ni,nj,nk,var)

                end if

                if(dmpvar(2:2).eq.'o') then

                  call vint31s(ni,nj,nk,z1d,zph8s,tmp2,var)

                  call outdmp3d('vbar  ',4,capt(33),ncpt(33),'xxx',     &
     &                          ni,nj,nk,var)

                  call vint31s(ni,nj,nk,z1d,zph8s,tmp4,var)

                  call outdmp3d('v     ',1,capt(31),ncpt(31),'xxx',     &
     &                          ni,nj,nk,var)

                else if(dmpvar(2:2).eq.'-') then

                  call vint31s(ni,nj,nk,z1d,zph8s,tmp4,var)

                  call outdmp3d('v     ',1,capt(31),ncpt(31),'xxx',     &
     &                          ni,nj,nk,var)

                end if

              else if(dmpvar(1:1).ne.'x'.and.dmpvar(2:2).ne.'x') then

                call var8u8s(ni,nj,nk,u,tmp1)
                call var8v8s(ni,nj,nk,v,tmp2)

                call rotuvm2s(idmpopt,idnspol,idtlon,                   &
     &                        2,ni-2,2,nj-2,2,nk-2,cpj,                 &
     &                        0,ni+1,0,nj+1,1,nk,lon,tmp1,tmp2)

                if(dmpvar(1:1).eq.'-') then

                  call vint31s(ni,nj,nk,z1d,zph8s,tmp1,var)

                  call outdmp3d('u     ',1,capt(30),ncpt(30),'xxx',     &
     &                          ni,nj,nk,var)

                end if

                if(dmpvar(2:2).eq.'-') then

                  call vint31s(ni,nj,nk,z1d,zph8s,tmp2,var)

                  call outdmp3d('v     ',1,capt(31),ncpt(31),'xxx',     &
     &                          ni,nj,nk,var)

                end if

              end if

            end if

! -----

! Dump the z components of velocity.

            if(dmpvar(3:3).eq.'o') then

              call var8w8s(ni,nj,nk,w,tmp1)

              call vint31s(ni,nj,nk,z1d,zph8s,tmp1,var)

              call outdmp3d('w     ',1,capt(3),ncpt(3),'xxx',           &
     &                      ni,nj,nk,var)

            end if

! -----

! Dump the base state pressure and its perturbation.

            if(dmpvar(4:4).eq.'o') then

              call vint31s(ni,nj,nk,z1d,zph8s,pbr,var)

              call outdmp3d('pbar  ',4,capt(4),ncpt(4),'xxx',           &
     &                      ni,nj,nk,var)

              call vint31s(ni,nj,nk,z1d,zph8s,pp,var)

              call outdmp3d('pp    ',2,capt(5),ncpt(5),'xxx',           &
     &                      ni,nj,nk,var)

            else if(dmpvar(4:4).eq.'-') then

              call totals(ni,nj,nk,pbr,pp,tmp1)

              call vint31s(ni,nj,nk,z1d,zph8s,tmp1,var)

              call outdmp3d('p     ',1,capt(26),ncpt(26),'xxx',         &
     &                      ni,nj,nk,var)

            end if

! -----

! Dump the base state potential temperature and its perturbation.

            if(dmpvar(5:5).eq.'o') then

              call vint31s(ni,nj,nk,z1d,zph8s,ptbr,var)

              call outdmp3d('ptbar ',5,capt(6),ncpt(6),'xxx',           &
     &                      ni,nj,nk,var)

              call vint31s(ni,nj,nk,z1d,zph8s,ptp,var)

              call outdmp3d('ptp   ',3,capt(7),ncpt(7),'xxx',           &
     &                      ni,nj,nk,var)

            else if(dmpvar(5:5).eq.'-') then

              call totals(ni,nj,nk,ptbr,ptp,tmp1)

              call vint31s(ni,nj,nk,z1d,zph8s,tmp1,var)

              call outdmp3d('pt    ',2,capt(27),ncpt(27),'xxx',         &
     &                      ni,nj,nk,var)

            end if

! -----

! Dump the water vapor mixing ratio.

            if(fmois(1:5).eq.'moist') then

              if(dmpvar(6:6).eq.'o') then

                call vint31s(ni,nj,nk,z1d,zph8s,qvbr,var)

                call outdmp3d('qvbar ',5,capt(8),ncpt(8),'xxx',         &
     &                        ni,nj,nk,var)

                call vint31s(ni,nj,nk,z1d,zph8s,qv,var)

                call outdmp3d('qv    ',2,capt(9),ncpt(9),'xxx',         &
     &                        ni,nj,nk,var)

              else if(dmpvar(6:6).eq.'-') then

                call vint31s(ni,nj,nk,z1d,zph8s,qv,var)

                call outdmp3d('qv    ',2,capt(9),ncpt(9),'xxx',         &
     &                        ni,nj,nk,var)

              end if

            end if

! -----

! Dump the water hydrometeor mixing ratio and the ice hydrometeor mixing
! ratio and concentrations.

            if(fmois(1:5).eq.'moist') then

              if(abs(cphopt).lt.10) then

                if(abs(cphopt).ge.1.and.dmpvar(7:7).eq.'o') then

                  call s_vint31s(ni,nj,nk,z1d,zph8s,qwtr(0,0,1,1),var)

                  call outdmp3d('qc    ',2,capt(10),ncpt(10),'xxx',     &
     &                          ni,nj,nk,var)

                  call s_vint31s(ni,nj,nk,z1d,zph8s,qwtr(0,0,1,2),var)

                  call outdmp3d('qr    ',2,capt(11),ncpt(11),'xxx',     &
     &                          ni,nj,nk,var)

                end if

                if(abs(cphopt).ge.2.and.dmpvar(7:7).eq.'o') then

                  call s_vint31s(ni,nj,nk,z1d,zph8s,qice(0,0,1,1),var)

                  call outdmp3d('qi    ',2,capt(14),ncpt(14),'xxx',     &
     &                          ni,nj,nk,var)

                  call s_vint31s(ni,nj,nk,z1d,zph8s,qice(0,0,1,2),var)

                  call outdmp3d('qs    ',2,capt(15),ncpt(15),'xxx',     &
     &                          ni,nj,nk,var)

                  call s_vint31s(ni,nj,nk,z1d,zph8s,qice(0,0,1,3),var)

                  call outdmp3d('qg    ',2,capt(16),ncpt(16),'xxx',     &
     &                          ni,nj,nk,var)

                  if(haiopt.eq.1) then

                    call s_vint31s(ni,nj,nk,z1d,zph8s,qice(0,0,1,4),var)

                    call outdmp3d('qh    ',2,capt(99),ncpt(99),'xxx',   &
     &                            ni,nj,nk,var)

                  end if

                end if

                if(abs(cphopt).eq.4.and.dmpvar(8:8).eq.'o') then

                  call s_vint31s(ni,nj,nk,z1d,zph8s,nwtr(0,0,1,1),var)

                  call outdmp3d('ncc   ',3,capt(74),ncpt(74),'xxx',     &
     &                          ni,nj,nk,var)

                  call s_vint31s(ni,nj,nk,z1d,zph8s,nwtr(0,0,1,2),var)

                  call outdmp3d('ncr   ',3,capt(75),ncpt(75),'xxx',     &
     &                          ni,nj,nk,var)

                end if

                if(abs(cphopt).ge.3.and.dmpvar(8:8).eq.'o') then

                  call s_vint31s(ni,nj,nk,z1d,zph8s,nice(0,0,1,1),var)

                  call outdmp3d('nci   ',3,capt(21),ncpt(21),'xxx',     &
     &                          ni,nj,nk,var)

                  call s_vint31s(ni,nj,nk,z1d,zph8s,nice(0,0,1,2),var)

                  call outdmp3d('ncs   ',3,capt(22),ncpt(22),'xxx',     &
     &                          ni,nj,nk,var)

                  call s_vint31s(ni,nj,nk,z1d,zph8s,nice(0,0,1,3),var)

                  call outdmp3d('ncg   ',3,capt(23),ncpt(23),'xxx',     &
     &                          ni,nj,nk,var)

                  if(haiopt.eq.1) then

                    call s_vint31s(ni,nj,nk,z1d,zph8s,nice(0,0,1,4),var)

                    call outdmp3d('nch   ',3,capt(102),ncpt(102),       &
     &                            'xxx',ni,nj,nk,var)

                  end if

                end if

                if(abs(cphopt).ge.2.and.dmpvar(7:7).eq.'o') then

                  if(flout_opt.eq.1) then

                    if(flqcqi_opt.ne.0) then

                      call vint31s(ni,nj,nk,z1d,zph8s,ucq,var)

                      call outdmp3d('ucq   ',3,capt(82),ncpt(82),       &
     &                              'xxx',ni,nj,nk,var)

                    end if

                    call vint31s(ni,nj,nk,z1d,zph8s,urq,var)

                    call outdmp3d('urq   ',3,capt(83),ncpt(83),'xxx',   &
     &                            ni,nj,nk,var)

                    if(flqcqi_opt.ne.0) then

                      call vint31s(ni,nj,nk,z1d,zph8s,uiq,var)

                      call outdmp3d('uiq   ',3,capt(84),ncpt(84),       &
     &                              'xxx',ni,nj,nk,var)

                    end if

                    call vint31s(ni,nj,nk,z1d,zph8s,usq,var)

                    call outdmp3d('usq   ',3,capt(85),ncpt(85),'xxx',   &
     &                            ni,nj,nk,var)

                    call vint31s(ni,nj,nk,z1d,zph8s,ugq,var)

                    call outdmp3d('ugq   ',3,capt(86),ncpt(86),'xxx',   &
     &                            ni,nj,nk,var)

                    if(haiopt.eq.1) then

                      call vint31s(ni,nj,nk,z1d,zph8s,uhq,var)

                      call outdmp3d('uhq   ',3,capt(104),ncpt(104),     &
     &                              'xxx',ni,nj,nk,var)

                    end if

                  end if

                end if

                if(cphopt.lt.0.and.dmpvar(9:9).eq.'o') then

                  if(qcgopt.eq.2) then

                   call s_vint31s(ni,nj,nk,z1d,zph8s,qcwtr(0,0,1,1),var)

                   call outdmp3d('qcc   ',3,capt(43),ncpt(43),'xxx',    &
     &                           ni,nj,nk,var)

                   call s_vint31s(ni,nj,nk,z1d,zph8s,qcwtr(0,0,1,2),var)

                   call outdmp3d('qrc   ',3,capt(44),ncpt(44),'xxx',    &
     &                           ni,nj,nk,var)

                  end if

                  call s_vint31s(ni,nj,nk,z1d,zph8s,qcice(0,0,1,1),var)

                  call outdmp3d('qic   ',3,capt(45),ncpt(45),'xxx',     &
     &                          ni,nj,nk,var)

                  call s_vint31s(ni,nj,nk,z1d,zph8s,qcice(0,0,1,2),var)

                  call outdmp3d('qsc   ',3,capt(46),ncpt(46),'xxx',     &
     &                          ni,nj,nk,var)

                  call s_vint31s(ni,nj,nk,z1d,zph8s,qcice(0,0,1,3),var)

                  call outdmp3d('qgc   ',3,capt(47),ncpt(47),'xxx',     &
     &                          ni,nj,nk,var)

                  if(haiopt.eq.1) then

                   call s_vint31s(ni,nj,nk,z1d,zph8s,qcice(0,0,1,4),var)

                   call outdmp3d('qhc   ',3,capt(103),ncpt(103),        &
     &                           'xxx',ni,nj,nk,var)

                  end if

                end if

              else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

                if(abs(cphopt).ge.11) then

                  if(dmpvar(7:7).eq.'o') then

                    call totalq(1,nqw,ni,nj,nk,nqw_sub,qwtr,tmp1)

                    call vint31s(ni,nj,nk,z1d,zph8s,tmp1,var)

                    call outdmp3d('qwtr  ',4,capt(35),ncpt(35),'xxx',   &
     &                            ni,nj,nk,var)

                  end if

                  if(dmpvar(8:8).eq.'o') then

                    call totalq(1,nnw,ni,nj,nk,nnw_sub,nwtr,tmp1)

                    call vint31s(ni,nj,nk,z1d,zph8s,tmp1,var)

                    call outdmp3d('nwtr  ',4,capt(36),ncpt(36),'xxx',   &
     &                            ni,nj,nk,var)

                  end if

                end if

                if(abs(cphopt).eq.12) then

                  if(dmpvar(7:7).eq.'o') then

                    call totalq(1,nqi,ni,nj,nk,nqi_sub,qice,tmp1)

                    call vint31s(ni,nj,nk,z1d,zph8s,tmp1,var)

                    call outdmp3d('qice  ',4,capt(39),ncpt(39),'xxx',   &
     &                            ni,nj,nk,var)

                  end if

                  if(dmpvar(8:8).eq.'o') then

                    call totalq(1,nni,ni,nj,nk,nni_sub,nice,tmp1)

                    call vint31s(ni,nj,nk,z1d,zph8s,tmp1,var)

                    call outdmp3d('nice  ',4,capt(40),ncpt(40),'xxx',   &
     &                            ni,nj,nk,var)

                  end if

                end if

              end if

            end if

! -----

! Dump the aerosol mixing ratio.

            if(aslopt.ge.1.and.dmpvar(10:10).eq.'o') then

              istr=1
              iend=nqa(1)-nqw-nqi

              call s_totalq(istr,iend,ni,nj,nk,nqa(0),qasl,tmp1)

              call vint31s(ni,nj,nk,z1d,zph8s,tmp1,var)

              call outdmp3d('drydu ',5,capt(87),ncpt(87),'xxx',         &
     &                      ni,nj,nk,var)

              istr=iend+1
              iend=iend+nqw+nqi

              call s_totalq(istr,iend,ni,nj,nk,nqa(0),qasl,tmp1)

              call vint31s(ni,nj,nk,z1d,zph8s,tmp1,var)

              call outdmp3d('wetdu ',5,capt(91),ncpt(91),'xxx',         &
     &                      ni,nj,nk,var)

              istr=iend+1
              iend=iend+nqa(2)-nqw-nqi

              call s_totalq(istr,iend,ni,nj,nk,nqa(0),qasl,tmp1)

              call vint31s(ni,nj,nk,z1d,zph8s,tmp1,var)

              call outdmp3d('dryca ',5,capt(88),ncpt(88),'xxx',         &
     &                      ni,nj,nk,var)

              istr=iend+1
              iend=iend+nqw+nqi

              call s_totalq(istr,iend,ni,nj,nk,nqa(0),qasl,tmp1)

              call vint31s(ni,nj,nk,z1d,zph8s,tmp1,var)

              call outdmp3d('wetca ',5,capt(92),ncpt(92),'xxx',         &
     &                      ni,nj,nk,var)

              istr=iend+1
              iend=iend+nqa(3)-nqw-nqi

              call s_totalq(istr,iend,ni,nj,nk,nqa(0),qasl,tmp1)

              call vint31s(ni,nj,nk,z1d,zph8s,tmp1,var)

              call outdmp3d('drysu ',5,capt(89),ncpt(89),'xxx',         &
     &                      ni,nj,nk,var)

              istr=iend+1
              iend=iend+nqw+nqi

              call s_totalq(istr,iend,ni,nj,nk,nqa(0),qasl,tmp1)

              call vint31s(ni,nj,nk,z1d,zph8s,tmp1,var)

              call outdmp3d('wetsu ',5,capt(93),ncpt(93),'xxx',         &
     &                      ni,nj,nk,var)

              istr=iend+1
              iend=iend+nqa(4)-nqw-nqi

              call s_totalq(istr,iend,ni,nj,nk,nqa(0),qasl,tmp1)

              call vint31s(ni,nj,nk,z1d,zph8s,tmp1,var)

              call outdmp3d('drysa ',5,capt(90),ncpt(90),'xxx',         &
     &                      ni,nj,nk,var)

              istr=iend+1
              iend=iend+nqw+nqi

              call s_totalq(istr,iend,ni,nj,nk,nqa(0),qasl,tmp1)

              call vint31s(ni,nj,nk,z1d,zph8s,tmp1,var)

              call outdmp3d('wetsa ',5,capt(94),ncpt(94),'xxx',         &
     &                      ni,nj,nk,var)

            end if

! -----

! Dump the tracer mixing ratio.

            if(trkopt.ge.1.and.dmpvar(11:11).eq.'o') then

              call vint31s(ni,nj,nk,z1d,zph8s,qt,var)

              call outdmp3d('qt    ',2,capt(34),ncpt(34),'xxx',         &
     &                      ni,nj,nk,var)

            end if

! -----

! Dump the turbulent kinetic energy and the maximum instantaneous wind
! velocity.

            if(tubopt.ge.2) then

              if(dmpvar(12:12).eq.'+'.or.dmpvar(12:12).eq.'o') then

                call vint31s(ni,nj,nk,z1d,zph8s,tke,var)

                call outdmp3d('tke   ',3,capt(25),ncpt(25),'xxx',       &
     &                        ni,nj,nk,var)

              end if

              if(dmpvar(12:12).eq.'+'.or.dmpvar(12:12).eq.'-') then

                call vint31s(ni,nj,nk,z1d,zph8s,maxvl,var)

                call outdmp3d('maxvl ',5,capt(73),ncpt(73),'xxx',       &
     &                        ni,nj,nk,var)

                call setcst3d(0,ni+1,0,nj+1,1,nk,lim35n,maxvl)

              end if

            end if

! -----

          end if

!! -----

        end if

!!! -----

! Dump the precipitation variables.

        if(fmon(1:3).eq.'act') then

          if(fmois(1:5).eq.'moist') then

            if(abs(cphopt).lt.10) then

              if(abs(cphopt).ge.1.and.dmpvar(14:14).eq.'o') then

                if(flqcqi_opt.ne.0) then

                  call s_outdmp2d('pcr   ',3,capt(95),ncpt(95),'xx',    &
     &                            ni,nj,prwtr(0,0,1,1))

                  call s_outdmp2d('pca   ',3,capt(96),ncpt(96),'xx',    &
     &                            ni,nj,prwtr(0,0,2,1))

                  call s_setcst2d(0,ni+1,0,nj+1,0.e0,prwtr(0,0,2,1))

                end if

                call s_outdmp2d('prr   ',3,capt(12),ncpt(12),'xx',      &
     &                          ni,nj,prwtr(0,0,1,2))

                call s_outdmp2d('pra   ',3,capt(13),ncpt(13),'xx',      &
     &                          ni,nj,prwtr(0,0,2,2))

                call s_setcst2d(0,ni+1,0,nj+1,0.e0,prwtr(0,0,2,2))

              end if

              if(abs(cphopt).ge.2.and.dmpvar(14:14).eq.'o') then

                if(flqcqi_opt.ne.0) then

                  call s_outdmp2d('pir   ',3,capt(97),ncpt(97),'xx',    &
     &                            ni,nj,price(0,0,1,1))

                  call s_outdmp2d('pia   ',3,capt(98),ncpt(98),'xx',    &
     &                            ni,nj,price(0,0,2,1))

                  call s_setcst2d(0,ni+1,0,nj+1,0.e0,price(0,0,2,1))

                end if

                call s_outdmp2d('psr   ',3,capt(17),ncpt(17),'xx',      &
     &                          ni,nj,price(0,0,1,2))

                call s_outdmp2d('psa   ',3,capt(18),ncpt(18),'xx',      &
     &                          ni,nj,price(0,0,2,2))

                call s_outdmp2d('pgr   ',3,capt(19),ncpt(19),'xx',      &
     &                          ni,nj,price(0,0,1,3))

                call s_outdmp2d('pga   ',3,capt(20),ncpt(20),'xx',      &
     &                          ni,nj,price(0,0,2,3))

                call s_setcst2d(0,ni+1,0,nj+1,0.e0,price(0,0,2,2))
                call s_setcst2d(0,ni+1,0,nj+1,0.e0,price(0,0,2,3))

                if(haiopt.eq.1) then

                  call s_outdmp2d('phr   ',3,capt(100),ncpt(100),       &
     &                            'xx',ni,nj,price(0,0,1,4))

                  call s_outdmp2d('pha   ',3,capt(101),ncpt(101),       &
     &                            'xx',ni,nj,price(0,0,2,4))

                  call s_setcst2d(0,ni+1,0,nj+1,0.e0,price(0,0,2,4))

                end if

              end if

            else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

              if(abs(cphopt).ge.11.and.dmpvar(14:14).eq.'o') then

                call s_outdmp2d('pwr   ',3,capt(37),ncpt(37),'xx',      &
     &                          ni,nj,prwtr(0,0,1,1))

                call s_outdmp2d('pwa   ',3,capt(38),ncpt(38),'xx',      &
     &                          ni,nj,prwtr(0,0,2,1))

                call s_setcst2d(0,ni+1,0,nj+1,0.e0,prwtr(0,0,2,1))

              end if

              if(abs(cphopt).eq.12.and.dmpvar(14:14).eq.'o') then

                call s_outdmp2d('pir   ',3,capt(41),ncpt(41),'xx',      &
     &                          ni,nj,price(0,0,1,1))

                call s_outdmp2d('pia   ',3,capt(42),ncpt(42),'xx',      &
     &                          ni,nj,price(0,0,2,1))

                call s_setcst2d(0,ni+1,0,nj+1,0.e0,price(0,0,2,1))

              end if

            end if

          end if

        end if

! -----

! Dump the z physical coordinates on terrain following coordinates.

        if(fdmp(1:3).eq.'act') then

          if(mod(dmplev,10).eq.1.and.dmpvar(15:15).eq.'o') then

            call outdmp3d('zph   ',3,capt(24),ncpt(24),'xxo',           &
     &                    ni,nj,nk,zph)

          end if

        end if

! -----

      end if

!!!! -----

      end subroutine s_outdmp

!-----7--------------------------------------------------------------7--

      end module m_outdmp
