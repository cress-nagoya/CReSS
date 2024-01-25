!***********************************************************************
      module m_abortslv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/10/06
!     Modification: 2004/05/31, 2004/09/01, 2005/02/10, 2007/01/20,
!                   2007/07/30, 2008/05/02, 2008/07/01, 2008/07/25,
!                   2008/08/25, 2008/10/10, 2009/01/05, 2009/02/27,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     check the file to abort the solver.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_comfile
      use m_comkind
      use m_commpi
      use m_cpondpe
      use m_destroy
      use m_getcname
      use m_getiname
      use m_getrname
      use m_inichar

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: abortslv, s_abortslv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface abortslv

        module procedure s_abortslv

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic int
      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_abortslv(fpcrsdir,fpnccrs,fpdmpmon,fpresopt,fpmxnopt,&
     &                      fpdmpitv,fpmonitv,fpresitv,fpmxnitv,ctime)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpcrsdir
                       ! Formal parameter of unique index of crsdir

      integer, intent(in) :: fpnccrs
                       ! Formal parameter of unique index of nccrs

      integer, intent(in) :: fpdmpmon
                       ! Formal parameter of unique index of dmpmon

      integer, intent(in) :: fpresopt
                       ! Formal parameter of unique index of resopt

      integer, intent(in) :: fpmxnopt
                       ! Formal parameter of unique index of mxnopt

      integer, intent(in) :: fpdmpitv
                       ! Formal parameter of unique index of dmpitv

      integer, intent(in) :: fpmonitv
                       ! Formal parameter of unique index of monitv

      integer, intent(in) :: fpresitv
                       ! Formal parameter of unique index of resitv

      integer, intent(in) :: fpmxnitv
                       ! Formal parameter of unique index of mxnitv

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

! Internal shared variables

      logical finq     ! Descriptor of file existance

      character(len=108) crsdir
                       ! User specified directory for CReSS files

      integer nccrs    ! Number of character of crsdir

      integer dmpmon   ! Option for monitor variables output
      integer resopt   ! Option for restart output
      integer mxnopt   ! Option for maxmum and minimum output

      integer fmot     ! Descriptor to put into motion

      integer stat     ! Runtime status

      real dmpitv      ! Time interval of dumped files
      real monitv      ! Time interval of dumped files
                       ! for monitor variables

      real resitv      ! Time interval of restart files

      real mxnitv      ! Time interval
                       ! of maximum and minimum output to standard i/o

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(crsdir)

! -----

! Get the required namelist variables.

      call getcname(fpcrsdir,crsdir)
      call getiname(fpnccrs,nccrs)
      call getiname(fpdmpmon,dmpmon)
      call getiname(fpresopt,resopt)
      call getiname(fpmxnopt,mxnopt)
      call getrname(fpdmpitv,dmpitv)
      call getrname(fpmonitv,monitv)
      call getrname(fpresitv,resitv)
      call getrname(fpmxnitv,mxnitv)

! -----

! Distinguish to put into motion.

      fmot=0

      if(mod(ctime,1000_i8*int(dmpitv+.1e0,i8)).eq.0_i8) then
        fmot=fmot+1
      end if

      if(dmpmon.eq.1                                                    &
     &  .and.mod(ctime,1000_i8*int(monitv+.1e0,i8)).eq.0_i8) then

        fmot=fmot+1

      end if

      if(resopt.eq.1                                                    &
     &  .and.mod(ctime,1000_i8*int(resitv+.1e0,i8)).eq.0_i8) then

        fmot=fmot+1

      end if

      if(mxnopt.eq.1                                                    &
     &  .and.mod(ctime,1000_i8*int(mxnitv+.1e0,i8)).eq.0_i8) then

        fmot=fmot+1

      end if

! -----

!! Check the file to abort the solver.

      if(fmot.ge.1) then

! Inquire the file, kill.solver.

        stat=0

        if(mype.eq.root) then

          inquire(file=crsdir(1:nccrs)//klfl(1:11),exist=finq)

          if(finq.eqv..true.) then
            stat=stat+1
          end if

          inquire(file=curdir(1:2)//klfl(1:11),exist=finq)

          if(finq.eqv..true.) then
            stat=stat+1
          end if

        end if

! -----

! If user create the file to abort the solver, call the procedure
! destroy.

        call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('abortslv',8,'cont',9,'              ',14,101, &
     &                   stat)

          end if

          call cpondpe

          call destroy('abortslv',8,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

! -----

      end if

!! -----

      end subroutine s_abortslv

!-----7--------------------------------------------------------------7--

      end module m_abortslv
