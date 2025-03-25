!***********************************************************************
      module m_outdmp3d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/03/25, 1999/04/06, 1999/08/23, 1999/09/30,
!                   1999/10/12, 1999/11/01, 1999/11/19, 2000/01/17,
!                   2000/04/18, 2000/07/05, 2001/01/15, 2001/03/13,
!                   2001/04/15, 2001/06/29, 2002/04/02, 2002/06/18,
!                   2003/03/28, 2003/04/30, 2003/05/19, 2003/06/27,
!                   2004/06/10, 2004/09/25, 2005/02/10, 2006/09/21,
!                   2007/01/05, 2007/01/20, 2008/01/11, 2008/04/17,
!                   2008/05/02, 2008/08/25, 2008/10/10, 2009/02/27,
!                   2011/09/22

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read in the 3 dimensional variables to the dumped file.

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

      public :: outdmp3d, s_outdmp3d

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outdmp3d

        module procedure s_outdmp3d

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
      subroutine s_outdmp3d(vname,ncvn,vcap,ncvc,xo,ni,nj,nk,var3d)
!***********************************************************************

! Input variables

      character(len=6), intent(in) :: vname
                       ! Optional variable name

      character(len=60), intent(in) :: vcap
                       ! Caption for dumped variable

      character(len=3), intent(in) :: xo
                       ! Control flag of variable arrangement

      integer, intent(in) :: ncvn
                       ! Number of character of vname

      integer, intent(in) :: ncvc
                       ! Number of character of vcap

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: var3d(0:ni+1,0:nj+1,1:nk)
                       ! 3 dimensional optional variable

! Internal shared variables

      integer dmpfmt   ! Option for dumped file format
      integer dmpmon   ! Option for monitor variables output

      integer stat     ! Runtime status

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer cntdmp   ! Counter of dumped variables

!-----7--------------------------------------------------------------7--

!! Read in the 3 dimensional variables to the dumped file.

      if(fdmp(1:3).eq.'act') then

! Get the required namelist variables.

        call getiname(iddmpfmt,dmpfmt)
        call getiname(iddmpmon,dmpmon)

! -----

! Initialize the runtime status.

        stat=1

! -----

! Increse the cnt3d.

        cnt3d=cnt3d+1

! -----

! Read in the variable at the u points to the dumped file.

        if(xo(1:3).eq.'oxx') then

          if(dmpfmt.eq.1) then

            do k=2,nk-2

              write(io3d,*,iostat=stat,err=100)                         &
     &          ((.5e0*(var3d(i,j,k)+var3d(i+1,j,k)),i=2,ni-2),j=2,nj-2)

            end do

          else if(dmpfmt.eq.2) then

            do k=2,nk-2
              rec3d=rec3d+1

              write(io3d,rec=rec3d,iostat=stat,err=100)                 &
     &          ((.5e0*(var3d(i,j,k)+var3d(i+1,j,k)),i=2,ni-2),j=2,nj-2)

            end do

          end if

! -----

! Read in the variable at the v points to the dumped file.

        else if(xo(1:3).eq.'xox') then

          if(dmpfmt.eq.1) then

            do k=2,nk-2

              write(io3d,*,iostat=stat,err=100)                         &
     &          ((.5e0*(var3d(i,j,k)+var3d(i,j+1,k)),i=2,ni-2),j=2,nj-2)

            end do

          else if(dmpfmt.eq.2) then

            do k=2,nk-2
              rec3d=rec3d+1

              write(io3d,rec=rec3d,iostat=stat,err=100)                 &
     &          ((.5e0*(var3d(i,j,k)+var3d(i,j+1,k)),i=2,ni-2),j=2,nj-2)

            end do

          end if

! -----

! Read in the variable at the w points to the dumped file.

        else if(xo(1:3).eq.'xxo') then

          if(dmpfmt.eq.1) then

            do k=2,nk-2

              write(io3d,*,iostat=stat,err=100)                         &
     &          ((.5e0*(var3d(i,j,k)+var3d(i,j,k+1)),i=2,ni-2),j=2,nj-2)

            end do

          else if(dmpfmt.eq.2) then

            do k=2,nk-2
              rec3d=rec3d+1

              write(io3d,rec=rec3d,iostat=stat,err=100)                 &
     &          ((.5e0*(var3d(i,j,k)+var3d(i,j,k+1)),i=2,ni-2),j=2,nj-2)

            end do

          end if

! -----

! Read in the variable at the scalar points to the dumped file.

        else if(xo(1:3).eq.'xxx') then

          if(dmpfmt.eq.1) then

            do k=2,nk-2

              write(io3d,*,iostat=stat,err=100)                         &
     &             ((var3d(i,j,k),i=2,ni-2),j=2,nj-2)

            end do

          else if(dmpfmt.eq.2) then

            do k=2,nk-2
              rec3d=rec3d+1

              write(io3d,rec=rec3d,iostat=stat,err=100)                 &
     &             ((var3d(i,j,k),i=2,ni-2),j=2,nj-2)

            end do

          end if

        end if

! -----

! If error occured, call the procedure destroy.

  100   call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('outdmp3d',8,'cont',4,'              ',14,     &
     &                   io3d,stat)

          end if

          call cpondpe

          call destroy('outdmp3d',8,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

! -----

! Read in the messages to the dumped data checking file.

        call outcap('dmp',vname,vcap,ncvc,io3c,nk)

! -----

! Read in the messages to the standard i/o.

        if(mype.eq.root) then

          if(dmpmon.eq.0) then

            call outstd09(vname,ncvn,vcap,ncvc,3,cnt3d)

          else

            cntdmp=cnt3d+cnt2d

            call outstd09(vname,ncvn,vcap,ncvc,3,cntdmp)

          end if

        end if

        call cpondpe

! -----

      end if

!! -----

      end subroutine s_outdmp3d

!-----7--------------------------------------------------------------7--

      end module m_outdmp3d
