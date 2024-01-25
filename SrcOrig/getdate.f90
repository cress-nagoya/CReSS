!***********************************************************************
      module m_getdate
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/05/10
!     Modification: 1999/06/07, 1999/06/21, 1999/11/01, 2000/01/05,
!                   2000/01/17, 2001/04/15, 2001/05/29, 2001/11/20,
!                   2002/06/18, 2003/04/30, 2003/05/19, 2003/11/05,
!                   2005/02/10, 2006/09/21, 2007/01/20, 2007/07/30,
!                   2008/05/02, 2008/08/25, 2008/10/10, 2009/01/05,
!                   2009/02/27, 2009/11/13, 2011/11/10

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the current forecast date from the lapse of forecast
!     time with second format, ssssssss.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_comdays
      use m_comkind
      use m_commpi
      use m_cpondpe
      use m_destroy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getdate, s_getdate

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getdate

        module procedure s_getdate

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_getdate(idate,ctime,cdate)
!***********************************************************************

! Input variables

      character(len=12), intent(in) :: idate
                       ! Forecast start date
                       ! with Gregorian calendar, yyyymmddhhmm

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time
                       ! with second format, ssssssss

! Output variable

      character(len=12), intent(out) :: cdate
                       ! Current forecast date
                       ! with Gregorian calendar, yyyymmddhhmm

! Internal shared variables

      integer stat     ! Runtime status

      integer i        ! Index of do loops

      integer iyr      ! Year of forecast start date
      integer imo      ! Month of forecast start date
      integer idy      ! Day of forecast start date
      integer ihr      ! Hour of forecast start date
      integer imn      ! Minite of forecast start date

      integer cyr      ! Year of current forecast date
      integer cmo      ! Month of current forecast date
      integer cdy      ! Day of current forecast date
      integer chr      ! Hour of current forecast date
      integer cmn      ! Minite of current forecast date

      integer cela     ! Temporary variable
      integer crem     ! Temporary variable

!-----7--------------------------------------------------------------7--

! Read out the integer variables from the forecast start date with
! Gregorian calendar, yyyymmddhhmm.

      read(idate(1:12),'(i4.4,4i2.2)') iyr,imo,idy,ihr,imn

! -----

!! Calculate the current forecast date from the lapse of forecast time.

      if(mod(ctime,60000_i8).eq.0_i8.and.ctime.le.99996000000_i8) then

! Get the current forecast date.

        cmn=ctime/60000_i8+imn

        chr=cmn/60+ihr
        cmn=mod(cmn,60)

        cdy=chr/24+idy
        chr=mod(chr,24)

        if(mod(iyr,400).eq.0                                            &
     &    .or.(mod(iyr,4).eq.0.and.mod(iyr,100).ne.0)) then

          crem=remitc(imo)

        else

          crem=rem(imo)

        end if

        if((cdy-crem).le.0) then

          cyr=iyr

          if(mod(cyr,400).eq.0                                          &
     &      .or.(mod(cyr,4).eq.0.and.mod(cyr,100).ne.0)) then

            if(imo.gt.1) then
              cela=elaitc(imo-1)
            else
              cela=0
            end if

          else

            if(imo.gt.1) then
              cela=ela(imo-1)
            else
              cela=0
            end if

          end if

          cmo=1
          cdy=cdy+cela

        else if((cdy-crem).gt.0) then

          cyr=iyr+1
          cmo=1
          cdy=cdy-crem

        end if

        imo=13

        if(mod(cyr,400).eq.0                                            &
     &    .or.(mod(cyr,4).eq.0.and.mod(cyr,100).ne.0)) then

          do i=1,12

            if(cdy.le.elaitc(i)) then

              if(i.lt.imo) then
                imo=i
              end if

            end if

          end do

          if(imo.gt.1) then

            cmo=cmo+imo-1
            cdy=cdy-elaitc(imo-1)

          end if

        else

          do i=1,12

            if(cdy.le.ela(i)) then

              if(i.lt.imo) then
                imo=i
              end if

            end if

          end do

          if(imo.gt.1) then

            cmo=cmo+imo-1
            cdy=cdy-ela(imo-1)

          end if

        end if

! -----

! Read in the current forecast date with Gregorian calendar,
! yyyymmddhhmm to the integer variables.

        write(cdate(1:12),'(i4.4,4i2.2)') cyr,cmo,cdy,chr,cmn

! -----

      end if

!! -----

! If error occured, call the procedure destroy.

      if(mod(ctime,60000_i8).eq.0_i8.and.ctime.le.99996000000_i8) then
        stat=0
      else
        stat=1
      end if

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('getdate ',7,'cont',7,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('getdate ',7,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

      end subroutine s_getdate

!-----7--------------------------------------------------------------7--

      end module m_getdate
