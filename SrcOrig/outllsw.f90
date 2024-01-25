!***********************************************************************
      module m_outllsw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/04/11
!     Modification: 2007/07/30, 2007/08/24, 2007/09/28, 2008/01/11,
!                   2008/04/17, 2008/05/02, 2008/08/25, 2008/10/10,
!                   2009/01/05, 2009/02/27, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read in the latitude and the longitude at south-west corner to
!     standard i/o.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_comkind
      use m_commpi
      use m_currpe
      use m_getcname
      use m_getiname
      use m_getllsw
      use m_inichar
      use m_outctl
      use m_outstd04
      use m_outstd05
      use m_outstd17

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: outllsw, s_outllsw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outllsw

        module procedure s_outllsw

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
      subroutine s_outllsw(fpfltyp_uni,fpdmpmon,fpuniopt_uni,fproc,     &
     &                     nstp0,nstp1,ni,nj,nk,tmp1,tmp2,tmp3,tmp4)
!***********************************************************************

! Input variables

      character(len=3), intent(in) :: fproc
                       ! Control flag of processing type

      integer, intent(in) :: fpfltyp_uni
                       ! Formal parameter of unique index of fltyp_uni

      integer, intent(in) :: fpdmpmon
                       ! Formal parameter of unique index of dmpmon

      integer, intent(in) :: fpuniopt_uni
                       ! Formal parameter of unique index of uniopt_uni

      integer(kind=i8), intent(in) :: nstp0
                       ! Start index of main do loop

      integer(kind=i8), intent(in) :: nstp1
                       ! End index of main do loop

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

! Internal shared variables

      character(len=108) fltyp_uni
                       ! Control flag of processed file type

      integer dmpmon   ! Option for monitor variables output

      integer uniopt_uni
                       ! Option for uniting process

      integer fout     ! Control flag of processing output routines

      real latsw       ! Latitude at south-west corner
      real lonsw       ! Longitude at south-west corner

      real, intent(inout) :: tmp1(1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp4(1:(nj-3)*njgrp*njsub)
                       ! Temporary array

!-----7--------------------------------------------------------------7--

! Initialize character variable.

      call inichar(fltyp_uni)

! -----

! Get the required namelist variables.

      call getcname(fpfltyp_uni,fltyp_uni)
      call getiname(fpdmpmon,dmpmon)
      call getiname(fpuniopt_uni,uniopt_uni)

! -----

! Set the control flag, fout.

      if(fltyp_uni(1:3).eq.'all') then

        if((fproc(1:3).eq.'dmp'.and.dmpmon.eq.0).or.                    &
     &     (fproc(1:3).eq.'mon'.and.dmpmon.eq.1)) then

          fout=1

        else

          fout=0

        end if

      else

        if(fltyp_uni(1:3).eq.fproc(1:3)) then

          fout=1

        else

          fout=0

        end if

      end if

      if(uniopt_uni.lt.-10) then

        fout=-fout

      end if

! -----

! Read in the message to standard i/o.

      if((fltyp_uni(1:3).eq.'all'.and.fproc(1:3).eq.'dmp').or.          &
     &   (fltyp_uni(1:3).eq.'dmp'.and.fproc(1:3).eq.'dmp').or.          &
     &   (fltyp_uni(1:3).eq.'mon'.and.fproc(1:3).eq.'mon')) then

        call outstd04(0,0_i8)

      end if

! -----

! Create the GrADS control file.

      if((mod(uniopt_uni,10).ge.-4.and.uniopt_uni.lt.0).or.             &
     &  (ngrp.eq.1.and.mod(uniopt_uni,10).le.-5)) then

        call outctl(idexprim,idcrsdir,idncexp,idnccrs,idmpopt,idnspol,  &
     &              idsthopt,iddmplev,iduniopt_uni,iddx,iddy,iddz,      &
     &              idtlat1,idtlat2,idtlon,idflitv_uni,                 &
     &              fproc,nstp0,nstp1,ni,nj,nk,                         &
     &              tmp1,tmp2,tmp3,tmp4)

      end if

! -----

! Calculate the latitude and the longitude at south-west corner.

      if(fout.eq.1) then

        if(abs(mod(uniopt_uni,10)).eq.1                                 &
     &    .or.abs(mod(uniopt_uni,10)).eq.3) then

          mygrp=0

          call currpe('unite   ',5,'ijgrp')

        else if(abs(mod(uniopt_uni,10)).eq.2                            &
     &    .or.abs(mod(uniopt_uni,10)).eq.4) then

          myred=0

          call currpe('unite   ',5,'ijred')

        else if(ngrp.eq.1.and.abs(mod(uniopt_uni,10)).ge.5) then

          mygrp=0

          call currpe('unite   ',5,'ijgrp')

        end if

        call getllsw(idmpopt,idnspol,idtlon,iddx,iddy,ni,nj,latsw,lonsw)

      end if

! -----

! Read in the message to standard i/o.

      if(fout.ne.0) then

        if(fout.eq.1) then

          call outstd17(latsw,lonsw)

        end if

        call outstd05(0)

      end if

! -----

      end subroutine s_outllsw

!-----7--------------------------------------------------------------7--

      end module m_outllsw
