!***********************************************************************
      module m_commpi
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/25
!     Modification: 1999/04/06, 1999/09/16, 1999/10/12, 2000/01/17,
!                   2000/02/07, 2000/03/23, 2001/02/13, 2001/06/06,
!                   2002/04/02, 2002/07/31, 2002/08/15, 2003/05/19,
!                   2004/08/20, 2005/08/05, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/04/11, 2008/05/02, 2008/08/25,
!                   2009/01/05, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the parameters of parallelizing.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      public

! Exceptional access control

!     none

!-----7--------------------------------------------------------------7--

! Module variables

      character(len=6) fpara
                       ! Control flag of parallelizing

      integer mype     ! My processor element number

      integer npe      ! Number of entried processor elements

      integer root     ! Root processor element number

      integer tag      ! Message tag

      integer mysub    ! My sub domain number in group domain

      integer nsub     ! Number of sub domain in group domain

      integer nisub    ! Number of sub domain in group domain
                       ! in x direction

      integer njsub    ! Number of sub domain in group domain
                       ! in y direction

      integer isub     ! Index of sub domain in group domain
                       ! in x direction

      integer jsub     ! Index of sub domain in group domain
                       ! in y direction

      integer mygrp    ! My group domain number in entire domain

      integer ngrp     ! Number of group domain in entire domain

      integer nigrp    ! Number of group domain in entire domain
                       ! in x direction

      integer njgrp    ! Number of group domain in entire domain
                       ! in y direction

      integer igrp     ! Index of group domain in entire domain
                       ! in x direction

      integer jgrp     ! Index of group domain in entire domain
                       ! in y direction

      integer mysrl    ! My serial number of group domain
                       ! in entire domain

      integer nsrl     ! Maximum serial number of group domain
                       ! in entire domain

      integer myred    ! My group domain number
                       ! in reductional entire domain

      integer nred     ! Number of group domain
                       ! in reductional entire domain

      integer nired    ! Number of group domain
                       ! in reductional entire domain in x direction

      integer njred    ! Number of group domain
                       ! in reductional entire domain in y direction

      integer ired     ! Index of group domain
                       ! in reductional entire domain in x direction

      integer jred     ! Index of group domain
                       ! in reductional entire domain in y direction

      integer iwred    ! Index of group domain at west boundary
                       ! in reductional entire domain

      integer jsred    ! Index of group domain at south boundary
                       ! in reductional entire domain

      integer ebw      ! Descriptor of closing to west boundary
                       ! in group domain location

      integer ebe      ! Descriptor of closing to east boundary
                       ! in group domain location

      integer ebs      ! Descriptor of closing to south boundary
                       ! in group domain location

      integer ebn      ! Descriptor of closing to north boundary
                       ! in group domain location

      integer ebsw     ! Descriptor of closing to south-west boundary
                       ! in group domain location

      integer ebse     ! Descriptor of closing to south-east boundary
                       ! in group domain location

      integer ebnw     ! Descriptor of closing to north-west boundary
                       ! in group domain location

      integer ebne     ! Descriptor of closing to north-east boundary
                       ! in group domain location

      integer dstw_sub ! West sending destination between sub domain
      integer dste_sub ! East sending destination between sub domain
      integer dsts_sub ! South sending destination between sub domain
      integer dstn_sub ! North sending destination between sub domain

      integer srcw_sub ! West receiving source between sub domain
      integer srce_sub ! East receiving source between sub domain
      integer srcs_sub ! South receiving source between sub domain
      integer srcn_sub ! North receiving source between sub domain

      integer dstw_grp ! West sending destination between group domain
      integer dste_grp ! East sending destination between group domain
      integer dsts_grp ! South sending destination between group domain
      integer dstn_grp ! North sending destination between group domain

      integer srcw_grp ! West receiving source between group domain
      integer srce_grp ! East receiving source between group domain
      integer srcs_grp ! South receiving source between group domain
      integer srcn_grp ! North receiving source between group domain

      integer dstw_sub_bnd
                       ! West sending destination
                       ! between sub domain for boundary processing

      integer dste_sub_bnd
                       ! East sending destination
                       ! between sub domain for boundary processing

      integer dsts_sub_bnd
                       ! South sending destination
                       ! between sub domain for boundary processing

      integer dstn_sub_bnd
                       ! North sending destination
                       ! between sub domain for boundary processing

      integer srcw_sub_bnd
                       ! West receiving source
                       ! between sub domain for boundary processing

      integer srce_sub_bnd
                       ! East receiving source
                       ! between sub domain for boundary processing

      integer srcs_sub_bnd
                       ! South receiving source
                       ! between sub domain for boundary processing

      integer srcn_sub_bnd
                       ! North receiving source
                       ! between sub domain for boundary processing

      integer dstw_grp_bnd
                       ! West sending destination
                       ! between group domain for boundary processing

      integer dste_grp_bnd
                       ! East sending destination
                       ! between group domain for boundary processing

      integer dsts_grp_bnd
                       ! South sending destination
                       ! between group domain for boundary processing

      integer dstn_grp_bnd
                       ! North sending destination
                       ! between group domain for boundary processing

      integer srcw_grp_bnd
                       ! West receiving source
                       ! between group domain for boundary processing

      integer srce_grp_bnd
                       ! East receiving source
                       ! between group domain for boundary processing

      integer srcs_grp_bnd
                       ! South receiving source
                       ! between group domain for boundary processing

      integer srcn_grp_bnd
                       ! North receiving source
                       ! between group domain for boundary processing

      integer mysub_rst
                       ! My sub domain number for rstruct

      integer nsub_rst ! Number of sub domain for rstruct

      integer nisub_rst
                       ! Number of sub domain in x direction for rstruct

      integer njsub_rst
                       ! Number of sub domain in y direction for rstruct

      integer isub_rst ! Index of sub domain in x direction for rstruct
      integer jsub_rst ! Index of sub domain in y direction for rstruct

      integer mpi_comm_cress
                       ! Communication world for all processor elements

      integer mpi_comm_cress_sub
                       ! Communication world
                       ! for processor elements in each groped domain

! Module procedure

!     none

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

!     none

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

!     none

!-----7--------------------------------------------------------------7--

      end module m_commpi
