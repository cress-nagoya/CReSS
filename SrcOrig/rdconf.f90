!***********************************************************************
      module m_rdconf
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/05/31
!     Modification: 2004/06/10, 2004/08/01, 2004/08/20, 2004/09/10,
!                   2005/01/14, 2005/02/10, 2005/08/05, 2006/04/03,
!                   2007/01/20, 2008/05/02, 2008/08/25, 2008/10/10,
!                   2008/12/11, 2009/02/27, 2009/03/31, 2011/08/18,
!                   2011/09/22

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read out the namelist variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_defname
      use m_destroy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: rdconf, s_rdconf

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rdconf

        module procedure s_rdconf

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
      subroutine s_rdconf(stat)
!***********************************************************************

! Output variable

      integer, intent(out) :: stat
                       ! Runtime status

!-----7--------------------------------------------------------------7--

! Initialize the runtime status.

      stat=0

! -----

! Read out the namelist variables.

      read(5,nml=sysdep,iostat=stat,end=100,err=100)
      read(5,nml=runame,iostat=stat,end=100,err=100)
      read(5,nml=drname,iostat=stat,end=100,err=100)
      read(5,nml=dimset,iostat=stat,end=100,err=100)
      read(5,nml=project,iostat=stat,end=100,err=100)
      read(5,nml=gridset,iostat=stat,end=100,err=100)
      read(5,nml=gridsth,iostat=stat,end=100,err=100)
      read(5,nml=terrain,iostat=stat,end=100,err=100)
      read(5,nml=flength,iostat=stat,end=100,err=100)
      read(5,nml=boundry,iostat=stat,end=100,err=100)
      read(5,nml=gpvpram,iostat=stat,end=100,err=100)
      read(5,nml=rdrpram,iostat=stat,end=100,err=100)
      read(5,nml=sfcphys,iostat=stat,end=100,err=100)
      read(5,nml=initype,iostat=stat,end=100,err=100)
      read(5,nml=gridmove,iostat=stat,end=100,err=100)
      read(5,nml=ptinicon,iostat=stat,end=100,err=100)
      read(5,nml=integrat,iostat=stat,end=100,err=100)
      read(5,nml=smoother,iostat=stat,end=100,err=100)
      read(5,nml=mapfcter,iostat=stat,end=100,err=100)
      read(5,nml=coriolis,iostat=stat,end=100,err=100)
      read(5,nml=earthcrv,iostat=stat,end=100,err=100)
      read(5,nml=buoyancy,iostat=stat,end=100,err=100)
      read(5,nml=diabatic,iostat=stat,end=100,err=100)
      read(5,nml=ddamping,iostat=stat,end=100,err=100)
      read(5,nml=cloudphy,iostat=stat,end=100,err=100)
      read(5,nml=asolproc,iostat=stat,end=100,err=100)
      read(5,nml=mixtrace,iostat=stat,end=100,err=100)
      read(5,nml=turbulen,iostat=stat,end=100,err=100)
      read(5,nml=outfomat,iostat=stat,end=100,err=100)

      read(5,nml=project_gpv,iostat=stat,end=100,err=100)
      read(5,nml=gridset_gpv,iostat=stat,end=100,err=100)
      read(5,nml=datconf_gpv,iostat=stat,end=100,err=100)

      read(5,nml=project_asl,iostat=stat,end=100,err=100)
      read(5,nml=gridset_asl,iostat=stat,end=100,err=100)
      read(5,nml=datconf_asl,iostat=stat,end=100,err=100)

      read(5,nml=project_rdr,iostat=stat,end=100,err=100)
      read(5,nml=gridset_rdr,iostat=stat,end=100,err=100)
      read(5,nml=datconf_rdr,iostat=stat,end=100,err=100)

      read(5,nml=project_trn,iostat=stat,end=100,err=100)
      read(5,nml=gridset_trn,iostat=stat,end=100,err=100)
      read(5,nml=datconf_trn,iostat=stat,end=100,err=100)

      read(5,nml=project_lnd,iostat=stat,end=100,err=100)
      read(5,nml=gridset_lnd,iostat=stat,end=100,err=100)
      read(5,nml=datconf_lnd,iostat=stat,end=100,err=100)

      read(5,nml=project_sst,iostat=stat,end=100,err=100)
      read(5,nml=gridset_sst,iostat=stat,end=100,err=100)

      read(5,nml=project_ice,iostat=stat,end=100,err=100)
      read(5,nml=gridset_ice,iostat=stat,end=100,err=100)

      read(5,nml=uniconf_uni,iostat=stat,end=100,err=100)

      read(5,nml=rstconf_rst,iostat=stat,end=100,err=100)

! -----

! If error occured, call the procedure destroy.

  100 if(stat.ne.0) then

        call destroy('rdconf  ',6,'cont',3,'              ',14,5,stat)

        return

      end if

! -----

      end subroutine s_rdconf

!-----7--------------------------------------------------------------7--

      end module m_rdconf
