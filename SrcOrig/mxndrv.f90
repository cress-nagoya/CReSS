!***********************************************************************
      module m_mxndrv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/03/23
!     Modification: 2000/04/18, 2000/06/01, 2000/07/05, 2001/01/15,
!                   2001/02/13, 2001/04/15, 2001/05/29, 2001/06/29,
!                   2002/04/02, 2002/07/15, 2002/12/02, 2003/03/28,
!                   2003/04/30, 2003/05/19, 2003/11/28, 2003/12/12,
!                   2004/04/15, 2004/05/31, 2004/06/10, 2004/09/01,
!                   2005/02/10, 2006/01/10, 2006/02/13, 2006/07/21,
!                   2006/09/11, 2006/09/21, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/07/30, 2007/11/26, 2008/05/02,
!                   2008/06/09, 2008/07/25, 2008/08/25, 2009/01/30,
!                   2009/02/27, 2011/08/18, 2011/09/22, 2011/11/10,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read in the maximum and minimum value of optional prognostic
!     variable to the standard i/o when the current forecast time
!     reaches marked time.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkstd
      use m_comindx
      use m_comkind
      use m_commpi
      use m_getcname
      use m_getiname
      use m_getrname
      use m_inichar
      use m_outmxn
      use m_outstd12
      use m_totalq

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: mxndrv, s_mxndrv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface mxndrv

        module procedure s_mxndrv

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic int
      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_mxndrv(fpmxnvar,fpcphopt,fphaiopt,fpqcgopt,          &
     &                    fpaslopt,fptrkopt,fptubopt,fpetime,fpmxnitv,  &
     &                    fmois,ctime,ni,nj,nk,nqw,nnw,nqi,nni,nqa,     &
     &                    u,v,w,pp,ptp,qv,qwtr,nwtr,qice,nice,          &
     &                    qcwtr,qcice,qasl,qt,tke,var)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: fpmxnvar
                       ! Formal parameter of unique index of mxnvar

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

      integer, intent(in) :: fpetime
                       ! Formal parameter of unique index of etime

      integer, intent(in) :: fpmxnitv
                       ! Formal parameter of unique index of mxnitv

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

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

! Internal shared variables

      character(len=108) mxnvar
                       ! Control flag of maximum and mininum output

      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes
      integer qcgopt   ! Option for charging distribution
      integer aslopt   ! Option for aerosol processes
      integer trkopt   ! Option for mixing ratio tracking
      integer tubopt   ! Option for turbulent mixing

      integer outcnt   ! Counter of output variables

      integer ni_sub   ! Substitute for ni
      integer nj_sub   ! Substitute for nj
      integer nk_sub   ! Substitute for nk

      integer nqw_sub  ! Substitute for nqw
      integer nnw_sub  ! Substitute for nnw
      integer nqi_sub  ! Substitute for nqi
      integer nni_sub  ! Substitute for nni

      integer istr     ! Start index of types of aerosol array
      integer iend     ! End index of types of aerosol array

      real etime       ! Forecast stop time

      real mxnitv      ! Time interval
                       ! of max/min output to standard i/o

      real, intent(inout) :: var(0:ni+1,0:nj+1,1:nk)
                       ! Temporary output variable

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getrname(fpetime,etime)
      call getrname(fpmxnitv,mxnitv)

! -----

!!!!! Read in the maximum and minimum value of optional prognostic
!!!!! variable to the standard i/o when the current forecast time
!!!!! reaches marked time.

      if(mod(ctime,1000_i8*int(mxnitv+.1e0,i8)).eq.0_i8                 &
     &  .or.ctime.eq.1000_i8*int(etime+.1e0,i8)) then

! Initialize the character variable.

        call inichar(mxnvar)

! -----

! Get the required namelist variables.

        call getcname(fpmxnvar,mxnvar)
        call getiname(fpcphopt,cphopt)
        call getiname(fphaiopt,haiopt)
        call getiname(fpqcgopt,qcgopt)
        call getiname(fpaslopt,aslopt)
        call getiname(fptrkopt,trkopt)
        call getiname(fptubopt,tubopt)

! -----

! Read in the message to standard i/o.

        if(mype.eq.root) then

          call outstd12(1,'    ',4,'       ',                           &
     &                  ctime,0,0,1,0.e0,0,0,1,0.e0)

        end if

        call chkstd(root)

! -----

! Initialize the outcnt.

        outcnt=0

! -----

! Set the substituted variables.

        ni_sub=ni
        nj_sub=nj
        nk_sub=nk

        nqw_sub=nqw
        nnw_sub=nnw
        nqi_sub=nqi
        nni_sub=nni

! -----

!!!! Read in the maximum and minimum value of optional prognostic
!!!! variable to the standard i/o.

! For x components of velocity.

        if(mxnvar(1:1).eq.'o') then

          call outmxn('u   ',1,'[m/s]  ',ctime,ni,nj,nk,                &
     &                1,ni_sub,1,nj_sub-1,1,nk_sub-1,outcnt,u)

        end if

! -----

! For y components of velocity.

        if(mxnvar(2:2).eq.'o') then

          call outmxn('v   ',1,'[m/s]  ',ctime,ni,nj,nk,                &
     &                1,ni_sub-1,1,nj_sub,1,nk_sub-1,outcnt,v)

        end if

! -----

! For z components of velocity.

        if(mxnvar(3:3).eq.'o') then

          call outmxn('w   ',1,'[m/s]  ',ctime,ni,nj,nk,                &
     &                1,ni_sub-1,1,nj_sub-1,1,nk_sub,outcnt,w)

        end if

! -----

! For pressure perturbation.

        if(mxnvar(4:4).eq.'o') then

          call outmxn('pp  ',2,'[Pa]   ',ctime,ni,nj,nk,                &
     &                1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,outcnt,pp)

        end if

! -----

! For potential temperature perturbation.

        if(mxnvar(5:5).eq.'o') then

          call outmxn('ptp ',3,'[K]    ',ctime,ni,nj,nk,                &
     &                1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,outcnt,ptp)

        end if

! -----

!!! For water vapor mixing ratio and water and ice hydrometeor.

        if(fmois(1:5).eq.'moist') then

! For water vapor mixing ratio.

          if(mxnvar(6:6).eq.'o') then

            call outmxn('qv  ',2,'[kg/kg]',ctime,ni,nj,nk,              &
     &                  1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,outcnt,qv)

          end if

! -----

!! For the bulk categories.

          if(abs(cphopt).lt.10) then

! For water mixing ratio.

            if(abs(cphopt).ge.1.and.mxnvar(7:7).eq.'o') then

              call s_outmxn('qc  ',2,'[kg/kg]',ctime,ni,nj,nk,          &
     &                      1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,           &
     &                      outcnt,qwtr(0,0,1,1))

              call s_outmxn('qr  ',2,'[kg/kg]',ctime,ni,nj,nk,          &
     &                      1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,           &
     &                      outcnt,qwtr(0,0,1,2))

            end if

! -----

! For ice mixing ratio.

            if(abs(cphopt).ge.2.and.mxnvar(7:7).eq.'o') then

              call s_outmxn('qi  ',2,'[kg/kg]',ctime,ni,nj,nk,          &
     &                      1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,           &
     &                      outcnt,qice(0,0,1,1))

              call s_outmxn('qs  ',2,'[kg/kg]',ctime,ni,nj,nk,          &
     &                      1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,           &
     &                      outcnt,qice(0,0,1,2))

              call s_outmxn('qg  ',2,'[kg/kg]',ctime,ni,nj,nk,          &
     &                      1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,           &
     &                      outcnt,qice(0,0,1,3))

              if(haiopt.eq.1) then

                call s_outmxn('qh  ',2,'[kg/kg]',ctime,ni,nj,nk,        &
     &                        1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,         &
     &                        outcnt,qice(0,0,1,4))

              end if

            end if

! -----

! For water concentrations.

            if(abs(cphopt).eq.4.and.mxnvar(8:8).eq.'o') then

              call s_outmxn('ncc ',3,'[1/kg] ',ctime,ni,nj,nk,          &
     &                      1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,           &
     &                      outcnt,nwtr(0,0,1,1))

              call s_outmxn('ncr ',3,'[1/kg] ',ctime,ni,nj,nk,          &
     &                      1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,           &
     &                      outcnt,nwtr(0,0,1,2))

            end if

! -----

! For ice concentrations.

            if(abs(cphopt).ge.3.and.mxnvar(8:8).eq.'o') then

              call s_outmxn('nci ',3,'[1/kg] ',ctime,ni,nj,nk,          &
     &                      1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,           &
     &                      outcnt,nice(0,0,1,1))

              call s_outmxn('ncs ',3,'[1/kg] ',ctime,ni,nj,nk,          &
     &                      1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,           &
     &                      outcnt,nice(0,0,1,2))

              call s_outmxn('ncg ',3,'[1/kg] ',ctime,ni,nj,nk,          &
     &                      1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,           &
     &                      outcnt,nice(0,0,1,3))

              if(haiopt.eq.1) then

                call s_outmxn('nch ',3,'[1/kg] ',ctime,ni,nj,nk,        &
     &                        1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,         &
     &                        outcnt,nice(0,0,1,4))

              end if

            end if

! -----

! For charging distributions of water hydrometeor.

            if(cphopt.lt.0.and.mxnvar(9:9).eq.'o') then

              if(qcgopt.eq.2) then

                call s_outmxn('qcc ',3,'[fC/kg]',ctime,ni,nj,nk,        &
     &                        1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,         &
     &                        outcnt,qcwtr(0,0,1,1))

                call s_outmxn('qrc ',3,'[fC/kg]',ctime,ni,nj,nk,        &
     &                        1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,         &
     &                        outcnt,qcwtr(0,0,1,2))

              end if

            end if

! -----

! For charging distributions of ice hydrometeor.

            if(cphopt.lt.0.and.mxnvar(9:9).eq.'o') then

              call s_outmxn('qic ',3,'[fC/kg]',ctime,ni,nj,nk,          &
     &                      1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,           &
     &                      outcnt,qcice(0,0,1,1))

              call s_outmxn('qsc ',3,'[fC/kg]',ctime,ni,nj,nk,          &
     &                      1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,           &
     &                      outcnt,qcice(0,0,1,2))

              call s_outmxn('qgc ',3,'[fC/kg]',ctime,ni,nj,nk,          &
     &                      1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,           &
     &                      outcnt,qcice(0,0,1,3))

              if(haiopt.eq.1) then

                call s_outmxn('qhc ',3,'[fC/kg]',ctime,ni,nj,nk,        &
     &                        1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,         &
     &                        outcnt,qcice(0,0,1,4))

              end if

            end if

! -----

!! -----

!! For the bin categories.

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

! For water mixing ratio.

            if(abs(cphopt).ge.11.and.mxnvar(7:7).eq.'o') then

              call totalq(1,nqw,ni,nj,nk,nqw_sub,qwtr,var)

              call outmxn('qwtr',4,'[kg/kg]',ctime,ni,nj,nk,            &
     &                    1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,outcnt,var)

            end if

! -----

! For water concentrations.

            if(abs(cphopt).ge.11.and.mxnvar(8:8).eq.'o') then

              call totalq(1,nnw,ni,nj,nk,nnw_sub,nwtr,var)

              call outmxn('nwtr',4,'[1/kg] ',ctime,ni,nj,nk,            &
     &                    1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,outcnt,var)

            end if

! -----

! For ice mixing ratio.

            if(abs(cphopt).eq.12.and.mxnvar(7:7).eq.'o') then

              call totalq(1,nqi,ni,nj,nk,nqi_sub,qice,var)

              call outmxn('qice',4,'[kg/kg]',ctime,ni,nj,nk,            &
     &                    1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,outcnt,var)

            end if

! -----

! For ice concentrations.

            if(abs(cphopt).eq.12.and.mxnvar(8:8).eq.'o') then

              call totalq(1,nni,ni,nj,nk,nni_sub,nice,var)

              call outmxn('nice',4,'[1/kg] ',ctime,ni,nj,nk,            &
     &                    1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,outcnt,var)

            end if

! -----

          end if

!! -----

        end if

!!! -----

! For aerosol.

        if(aslopt.ge.1.and.mxnvar(10:10).eq.'o') then

          istr=1
          iend=nqa(1)

          call s_totalq(istr,iend,ni,nj,nk,nqa(0),qasl,var)

          call outmxn('du  ',2,'[kg/kg]',ctime,ni,nj,nk,                &
     &                1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,outcnt,var)

          istr=iend+1
          iend=iend+nqa(2)

          call s_totalq(istr,iend,ni,nj,nk,nqa(0),qasl,var)

          call outmxn('ca  ',2,'[kg/kg]',ctime,ni,nj,nk,                &
     &                1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,outcnt,var)

          istr=iend+1
          iend=iend+nqa(3)

          call s_totalq(istr,iend,ni,nj,nk,nqa(0),qasl,var)

          call outmxn('su  ',2,'[kg/kg]',ctime,ni,nj,nk,                &
     &                1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,outcnt,var)

          istr=iend+1
          iend=iend+nqa(4)

          call s_totalq(istr,iend,ni,nj,nk,nqa(0),qasl,var)

          call outmxn('sa  ',2,'[kg/kg]',ctime,ni,nj,nk,                &
     &                1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,outcnt,var)

        end if

! -----

! For tracer.

        if(trkopt.ge.1.and.mxnvar(11:11).eq.'o') then

          call outmxn('qt  ',2,'[kg/kg]',ctime,ni,nj,nk,                &
     &                1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,outcnt,qt)

        end if

! -----

! For turbulent kinetic energy.

        if(tubopt.ge.2.and.mxnvar(12:12).eq.'o') then

          call outmxn('tke ',3,'[J/kg] ',ctime,ni,nj,nk,                &
     &                1,ni_sub-1,1,nj_sub-1,1,nk_sub-1,outcnt,tke)

        end if

! -----

!!!! -----

! Read in the message to standard i/o.

        if(outcnt.eq.0) then

          if(mype.eq.root) then

            call outstd12(3,'    ',4,'       ',                         &
     &                    ctime,0,0,1,0.e0,0,0,1,0.e0)

          end if

        end if

        call chkstd(root)

! -----

      end if

!!!!! -----

      end subroutine s_mxndrv

!-----7--------------------------------------------------------------7--

      end module m_mxndrv
