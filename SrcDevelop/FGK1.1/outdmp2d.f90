!***********************************************************************
      module m_outdmp2d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/03/25, 1999/04/06, 1999/08/23, 1999/09/30,
!                   1999/10/12, 1999/11/19, 2000/01/17, 2000/04/18,
!                   2000/07/05, 2001/01/15, 2001/03/13, 2001/04/15,
!                   2001/06/29, 2002/04/02, 2002/06/18, 2003/03/28,
!                   2003/04/30, 2003/05/19, 2003/06/27, 2004/06/10,
!                   2004/09/25, 2005/02/10, 2006/09/21, 2006/12/04,
!                   2007/01/05, 2007/01/20, 2008/01/11, 2008/04/17,
!                   2008/05/02, 2008/07/01, 2008/08/25, 2008/10/10,
!                   2009/02/27, 2011/09/22

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read in the 2 dimensional variables to the dumped file.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_comdmp
      use m_comindx
      use m_commpi
      use m_cpondpe
      use m_destroy
      use m_getiname
      use m_outcap
      use m_outstd09

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: outdmp2d, s_outdmp2d

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outdmp2d

        module procedure s_outdmp2d

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
      subroutine s_outdmp2d(vname,ncvn,vcap,ncvc,xo,ni,nj,var2d)
!***********************************************************************

! Input variables

      character(len=6), intent(in) :: vname
                       ! Optional variable name

      character(len=60), intent(in) :: vcap
                       ! Caption for dumped variable

      character(len=2), intent(in) :: xo
                       ! Control flag of variable arrangement

      integer, intent(in) :: ncvn
                       ! Number of character of vname

      integer, intent(in) :: ncvc
                       ! Number of character of vcap

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      real, intent(in) :: var2d(0:ni+1,0:nj+1)
                       ! 2 dimensional optional variable

! Internal shared variables

      integer dmpfmt   ! Option for dumped file format
      integer dmpmon   ! Option for monitor variables output

      integer stat     ! Runtime status

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      integer cntdmp   ! Counter of dumped variables

!-----7--------------------------------------------------------------7--

!!! Read in the 2 dimensional variables to the dumped file.

      if(fmon(1:3).eq.'act') then

! Get the required namelist variables.

        call getiname(iddmpfmt,dmpfmt)
        call getiname(iddmpmon,dmpmon)

! -----

! Initialize the runtime status.

        stat=1

! -----

!! Read in the 2 dimensional monitor variables to the 3 dimensional
!! dumped file.

        if(dmpmon.eq.0) then

! Increse the cnt3d.

          cnt3d=cnt3d+1

! -----

! Read in the variable at the u points to the dumped file.

          if(xo(1:2).eq.'ox') then

            if(dmpfmt.eq.1) then

              write(io3d,*,iostat=stat,err=100)                         &
     &             ((.5e0*(var2d(i,j)+var2d(i+1,j)),i=2,ni-2),j=2,nj-2)

            else if(dmpfmt.eq.2) then

              rec3d=rec3d+1

              write(io3d,rec=rec3d,iostat=stat,err=100)                 &
     &             ((.5e0*(var2d(i,j)+var2d(i+1,j)),i=2,ni-2),j=2,nj-2)

            end if

! -----

! Read in the variable at the v points to the dumped file.

          else if(xo(1:2).eq.'xo') then

            if(dmpfmt.eq.1) then

              write(io3d,*,iostat=stat,err=100)                         &
     &             ((.5e0*(var2d(i,j)+var2d(i,j+1)),i=2,ni-2),j=2,nj-2)

            else if(dmpfmt.eq.2) then

              rec3d=rec3d+1

              write(io3d,rec=rec3d,iostat=stat,err=100)                 &
     &             ((.5e0*(var2d(i,j)+var2d(i,j+1)),i=2,ni-2),j=2,nj-2)

            end if

! -----

! Read in the variable at the scalar points to the dumped file.

          else if(xo(1:2).eq.'xx') then

            if(dmpfmt.eq.1) then

              write(io3d,*,iostat=stat,err=100)                         &
     &             ((var2d(i,j),i=2,ni-2),j=2,nj-2)

            else if(dmpfmt.eq.2) then

              rec3d=rec3d+1

              write(io3d,rec=rec3d,iostat=stat,err=100)                 &
     &             ((var2d(i,j),i=2,ni-2),j=2,nj-2)

            end if

          end if

! -----

! If error occured, call the procedure destroy.

  100     call chkerr(stat)

          if(stat.lt.0) then

            if(mype.eq.-stat-1) then

              call destroy('outdmp2d',8,'cont',4,'              ',14,   &
     &                     io3d,stat)

            end if

            call cpondpe

            call destroy('outdmp2d',8,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

! -----

! Read in the messages to the dumped data checking file.

          call outcap('dmp',vname,vcap,ncvc,io3c,3)

! -----

! Read in the messages to the standard i/o.

          if(mype.eq.root) then

            call outstd09(vname,ncvn,vcap,ncvc,2,cnt3d)

          end if

          call cpondpe

! -----

!! -----

!! Read in the 2 dimensional monitor variables to the 2 dimensional
!! dumped file.

        else

! Increse the cnt2d.

          cnt2d=cnt2d+1

! -----

! Read in the variable at the u points to the dumped file.

          if(xo(1:2).eq.'ox') then

            if(dmpfmt.eq.1) then

              write(io2d,*,iostat=stat,err=110)                         &
     &             ((.5e0*(var2d(i,j)+var2d(i+1,j)),i=2,ni-2),j=2,nj-2)

            else if(dmpfmt.eq.2) then

              rec2d=rec2d+1

              write(io2d,rec=rec2d,iostat=stat,err=110)                 &
     &             ((.5e0*(var2d(i,j)+var2d(i+1,j)),i=2,ni-2),j=2,nj-2)

            end if

! -----

! Read in the variable at the v points to the dumped file.

          else if(xo(1:2).eq.'xo') then

            if(dmpfmt.eq.1) then

              write(io2d,*,iostat=stat,err=110)                         &
     &             ((.5e0*(var2d(i,j)+var2d(i,j+1)),i=2,ni-2),j=2,nj-2)

            else if(dmpfmt.eq.2) then

              rec2d=rec2d+1

              write(io2d,rec=rec2d,iostat=stat,err=110)                 &
     &             ((.5e0*(var2d(i,j)+var2d(i,j+1)),i=2,ni-2),j=2,nj-2)

            end if

! -----

! Read in the variable at the scalar points to the dumped file.

          else if(xo(1:2).eq.'xx') then

            if(dmpfmt.eq.1) then

              write(io2d,*,iostat=stat,err=110)                         &
     &             ((var2d(i,j),i=2,ni-2),j=2,nj-2)

            else if(dmpfmt.eq.2) then

              rec2d=rec2d+1

              write(io2d,rec=rec2d,iostat=stat,err=110)                 &
     &             ((var2d(i,j),i=2,ni-2),j=2,nj-2)

            end if

          end if

! -----

! If error occured, call the procedure destroy.

  110     call chkerr(stat)

          if(stat.lt.0) then

            if(mype.eq.-stat-1) then

              call destroy('outdmp2d',8,'cont',4,'              ',14,   &
     &                     io2d,stat)

            end if

            call cpondpe

            call destroy('outdmp2d',8,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

! -----

! Read in the messages to the dumped data checking file.

          call outcap('dmp',vname,vcap,ncvc,io2c,3)

! -----

! Read in the messages to the standard i/o.

          if(mype.eq.root) then

            cntdmp=cnt3d+cnt2d

            call outstd09(vname,ncvn,vcap,ncvc,2,cntdmp)

          end if

          call cpondpe

! -----

        end if

!! -----

      end if

!!! -----

      end subroutine s_outdmp2d

!-----7--------------------------------------------------------------7--

      end module m_outdmp2d
