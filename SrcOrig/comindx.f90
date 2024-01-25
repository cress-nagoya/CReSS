!***********************************************************************
      module m_comindx
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/12/18
!     Modification: 2001/01/15, 2001/03/13, 2001/04/15, 2001/06/06,
!                   2001/08/07, 2001/10/18, 2001/11/14, 2001/11/20,
!                   2002/02/05, 2002/04/02, 2002/06/18, 2002/07/03,
!                   2002/08/15, 2002/08/27, 2002/09/02, 2002/09/09,
!                   2002/10/15, 2002/10/31, 2002/11/11, 2003/03/28,
!                   2003/05/19, 2003/07/15, 2003/08/08, 2003/09/01,
!                   2003/10/10, 2003/11/05, 2003/12/12, 2004/01/09,
!                   2004/03/05, 2004/04/15, 2004/05/07, 2004/05/31,
!                   2004/08/01, 2004/08/20, 2004/09/01, 2005/01/14,
!                   2005/02/10, 2005/08/05, 2005/12/13, 2006/04/03,
!                   2006/07/21, 2006/09/30, 2006/11/27, 2006/12/04,
!                   2007/01/05, 2007/01/20, 2007/03/23, 2007/04/24,
!                   2007/05/21, 2007/07/30, 2008/01/11, 2008/04/17,
!                   2008/05/02, 2008/07/01, 2008/08/25, 2008/10/10,
!                   2008/12/11, 2009/01/05, 2009/01/30, 2009/02/27,
!                   2011/05/16, 2011/08/09, 2011/08/18, 2011/09/22,
!                   2011/11/10

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the unique indices for each namelist variable in the
!     namelist table.

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

      integer, parameter :: idwlngth=2
                       ! Unique index of wlngth in namelist table

      integer, parameter :: idsavmem=170
                       ! Unique index of savmem in namelist table

      integer, parameter :: idexprim=1
                       ! Unique index of exprim in namelist table

      integer, parameter :: idncexp=1
                       ! Unique index of ncexp in namelist table

      integer, parameter :: idcrsdir=19
                       ! Unique index of crsdir in namelist table

      integer, parameter :: iddatdir=20
                       ! Unique index of datdir in namelist table

      integer, parameter :: idnccrs=181
                       ! Unique index of nccrs in namelist table

      integer, parameter :: idncdat=182
                       ! Unique index of ncdat in namelist table

      integer, parameter :: idnumpe=209
                       ! Unique index of numpe in namelist table

      integer, parameter :: idxgroup=210
                       ! Unique index of xgroup in namelist table

      integer, parameter :: idygroup=211
                       ! Unique index of ygroup in namelist table

      integer, parameter :: idxsub=212
                       ! Unique index of xsub in namelist table

      integer, parameter :: idysub=213
                       ! Unique index of ysub in namelist table

      integer, parameter :: idxdim=183
                       ! Unique index of xdim in namelist table

      integer, parameter :: idydim=184
                       ! Unique index of ydim in namelist table

      integer, parameter :: idzdim=185
                       ! Unique index of zdim in namelist table

      integer, parameter :: idmpopt=6
                       ! Unique index of mpopt in namelist table

      integer, parameter :: idnspol=7
                       ! Unique index of nspol in namelist table

      integer, parameter :: idtlat1=21
                       ! Unique index of tlat1 in namelist table

      integer, parameter :: idtlat2=22
                       ! Unique index of tlat2 in namelist table

      integer, parameter :: idtlon=23
                       ! Unique index of tlon in namelist table

      integer, parameter :: iddisr=783
                       ! Unique index of disr in namelist table

      integer, parameter :: iddx=5
                       ! Unique index of dx in namelist table

      integer, parameter :: iddy=6
                       ! Unique index of dy in namelist table

      integer, parameter :: iddz=7
                       ! Unique index of dz in namelist table

      integer, parameter :: iddxiv=8
                       ! Unique index of dxiv in namelist table

      integer, parameter :: iddyiv=9
                       ! Unique index of dyiv in namelist table

      integer, parameter :: iddziv=10
                       ! Unique index of dziv in namelist table

      integer, parameter :: idulat=17
                       ! Unique index of ulat in namelist table

      integer, parameter :: idulon=18
                       ! Unique index of ulon in namelist table

      integer, parameter :: idriu=19
                       ! Unique index of riu in namelist table

      integer, parameter :: idrju=20
                       ! Unique index of rju in namelist table

      integer, parameter :: idzsfc=12
                       ! Unique index of zsfc in namelist table

      integer, parameter :: idzflat=11
                       ! Unique index of zflat in namelist table

      integer, parameter :: idsthopt=5
                       ! Unique index of sthopt in namelist table

      integer, parameter :: iddzmin=13
                       ! Unique index of dzmin in namelist table

      integer, parameter :: idlayer1=14
                       ! Unique index of layer1 in namelist table

      integer, parameter :: idlayer2=15
                       ! Unique index of layer2 in namelist table

      integer, parameter :: idtrnopt=17
                       ! Unique index of trnopt in namelist table

      integer, parameter :: idmnthgh=24
                       ! Unique index of mnthgh in namelist table

      integer, parameter :: idmntwx=26
                       ! Unique index of mntwx in namelist table

      integer, parameter :: idmntwy=27
                       ! Unique index of mntwy in namelist table

      integer, parameter :: idmntcx=28
                       ! Unique index of mntcx in namelist table

      integer, parameter :: idmntcy=29
                       ! Unique index of mntcy in namelist table

      integer, parameter :: ididate=2
                       ! Unique index of idate in namelist table

      integer, parameter :: idstime=1
                       ! Unique index of stime in namelist table

      integer, parameter :: idetime=2
                       ! Unique index of etime in namelist table

      integer, parameter :: idwbc=8
                       ! Unique index of wbc in namelist table

      integer, parameter :: idebc=9
                       ! Unique index of ebc in namelist table

      integer, parameter :: idsbc=10
                       ! Unique index of sbc in namelist table

      integer, parameter :: idnbc=11
                       ! Unique index of nbc in namelist table

      integer, parameter :: idbbc=12
                       ! Unique index of bbc in namelist table

      integer, parameter :: idtbc=13
                       ! Unique index of tbc in namelist table

      integer, parameter :: idlbcvar=17
                       ! Unique index of lbcvar in namelist table

      integer, parameter :: idlbnews=447
                       ! Unique index of lbnews in namelist table

      integer, parameter :: idlbnorm=463
                       ! Unique index of lbnorm in namelist table

      integer, parameter :: idgwave=770
                       ! Unique index of gwave in namelist table

      integer, parameter :: idgpvvar=5
                       ! Unique index of gpvvar in namelist table

      integer, parameter :: idgpvitv=54
                       ! Unique index of gpvitv in namelist table

      integer, parameter :: idaslitv=784
                       ! Unique index of aslitv in namelist table

      integer, parameter :: idgsmopt=171
                       ! Unique index of gsmopt in namelist table

      integer, parameter :: idgsmcnt=172
                       ! Unique index of gsmcnt in namelist table

      integer, parameter :: idgsmcoe=448
                       ! Unique index of gsmcoe in namelist table

      integer, parameter :: idnggopt=23
                       ! Unique index of nggopt in namelist table

      integer, parameter :: idnggvar=6
                       ! Unique index of nggvar in namelist table

      integer, parameter :: idnggcoe=66
                       ! Unique index of nggcoe in namelist table

      integer, parameter :: idnggdlt=67
                       ! Unique index of nggdlt in namelist table

      integer, parameter :: idnggstr=38
                       ! Unique index of nggstr in namelist table

      integer, parameter :: idnggend=68
                       ! Unique index of nggend in namelist table

      integer, parameter :: idnggc20=59
                       ! Unique index of nggc20 in namelist table

      integer, parameter :: idexbopt=30
                       ! Unique index of exbopt in namelist table

      integer, parameter :: idexbvar=11
                       ! Unique index of exbvar in namelist table

      integer, parameter :: idexnews=446
                       ! Unique index of exnews in namelist table

      integer, parameter :: idexnorm=462
                       ! Unique index of exnorm in namelist table

      integer, parameter :: idexbwid=59
                       ! Unique index of exbwid in namelist table

      integer, parameter :: idlspopt=58
                       ! Unique index of lspopt in namelist table

      integer, parameter :: idlspvar=7
                       ! Unique index of lspvar in namelist table

      integer, parameter :: idlspsmt=445
                       ! Unique index of lspsmt in namelist table

      integer, parameter :: idlsnews=58
                       ! Unique index of lsnews in namelist table

      integer, parameter :: idlsnorm=460
                       ! Unique index of lsnorm in namelist table

      integer, parameter :: idwdnews=32
                       ! Unique index of wdnews in namelist table

      integer, parameter :: idwdnorm=178
                       ! Unique index of wdnorm in namelist table

      integer, parameter :: idvspopt=31
                       ! Unique index of vspopt in namelist table

      integer, parameter :: idvspvar=10
                       ! Unique index of vspvar in namelist table

      integer, parameter :: idvspgpv=56
                       ! Unique index of vspgpv in namelist table

      integer, parameter :: idvspbar=112
                       ! Unique index of vspbar in namelist table

      integer, parameter :: idbotgpv=57
                       ! Unique index of botgpv in namelist table

      integer, parameter :: idbotbar=114
                       ! Unique index of botbar in namelist table

      integer, parameter :: idrdrvar=13
                       ! Unique index of rdrvar in namelist table

      integer, parameter :: idrdritv=416
                       ! Unique index of rdritv in namelist table

      integer, parameter :: idngropt=162
                       ! Unique index of ngropt in namelist table

      integer, parameter :: idngrvar=14
                       ! Unique index of ngrvar in namelist table

      integer, parameter :: idngrcoe=417
                       ! Unique index of ngrcoe in namelist table

      integer, parameter :: idngrdlt=418
                       ! Unique index of ngrdlt in namelist table

      integer, parameter :: idngrstr=419
                       ! Unique index of ngrstr in namelist table

      integer, parameter :: idngrend=420
                       ! Unique index of ngrend in namelist table

      integer, parameter :: idngrc20=780
                       ! Unique index of ngrc20 in namelist table

      integer, parameter :: idngraff=421
                       ! Unique index of ngraff in namelist table

      integer, parameter :: idsfcdat=9
                       ! Unique index of sfcdat in namelist table

      integer, parameter :: idsfcopt=46
                       ! Unique index of sfcopt in namelist table

      integer, parameter :: idlevpbl=47
                       ! Unique index of levpbl in namelist table

      integer, parameter :: idlevund=60
                       ! Unique index of levund in namelist table

      integer, parameter :: iddtgrd=61
                       ! Unique index of dtgrd in namelist table

      integer, parameter :: iddzgrd=71
                       ! Unique index of dzgrd in namelist table

      integer, parameter :: iddzsea=461
                       ! Unique index of dzsea in namelist table

      integer, parameter :: idtgdeep=78
                       ! Unique index of tgdeep in namelist table

      integer, parameter :: idprvres=12
                       ! Unique index of prvres in namelist table

      integer, parameter :: idncprv=61
                       ! Unique index of ncprv in namelist table

      integer, parameter :: idlnduse=179
                       ! Unique index of lnduse in namelist table

      integer, parameter :: idgralbe=73
                       ! Unique index of gralbe in namelist table

      integer, parameter :: idgrbeta=72
                       ! Unique index of grbeta in namelist table

      integer, parameter :: idgrz0m=74
                       ! Unique index of grz0m in namelist table

      integer, parameter :: idgrz0h=465
                       ! Unique index of grz0h in namelist table

      integer, parameter :: idgrcap=766
                       ! Unique index of grcap in namelist table

      integer, parameter :: idgrnuu=767
                       ! Unique index of grnuu in namelist table

      integer, parameter :: idsstcst=75
                       ! Unique index of sstcst in namelist table

      integer, parameter :: idsstitv=796
                       ! Unique index of sstitv in namelist table

      integer, parameter :: iddstopt=169
                       ! Unique index of dstopt in namelist table

      integer, parameter :: idiniopt=14
                       ! Unique index of iniopt in namelist table

      integer, parameter :: idsnddim=186
                       ! Unique index of snddim in namelist table

      integer, parameter :: idsndtyp=4
                       ! Unique index of sndtyp in namelist table

      integer, parameter :: idzsnd0=36
                       ! Unique index of zsnd0 in namelist table

      integer, parameter :: idpsnd0=37
                       ! Unique index of psnd0 in namelist table

      integer, parameter :: idmasopt=165
                       ! Unique index of masopt in namelist table

      integer, parameter :: idmaseps=449
                       ! Unique index of maseps in namelist table

      integer, parameter :: idalpha1=65
                       ! Unique index of alpha1 in namelist table

      integer, parameter :: idalpha2=16
                       ! Unique index of alpha2 in namelist table

      integer, parameter :: idmovopt=16
                       ! Unique index of movopt in namelist table

      integer, parameter :: idumove=55
                       ! Unique index of umove in namelist table

      integer, parameter :: idvmove=60
                       ! Unique index of vmove in namelist table

      integer, parameter :: idpt0opt=18
                       ! Unique index of pt0opt in namelist table

      integer, parameter :: idpt0num=34
                       ! Unique index of pt0num in namelist table

      integer, parameter :: idptp0=768
                       ! Unique index of ptp0 in namelist table

      integer, parameter :: idpt0rx=30
                       ! Unique index of pt0rx in namelist table

      integer, parameter :: idpt0ry=31
                       ! Unique index of pt0ry in namelist table

      integer, parameter :: idpt0rz=32
                       ! Unique index of pt0rz in namelist table

      integer, parameter :: idpt0cx=33
                       ! Unique index of pt0cx in namelist table

      integer, parameter :: idpt0cy=34
                       ! Unique index of pt0cy in namelist table

      integer, parameter :: idpt0cz=35
                       ! Unique index of pt0cz in namelist table

      integer, parameter :: idpt0ds=69
                       ! Unique index of pt0ds in namelist table

      integer, parameter :: iddtbig=3
                       ! Unique index of dtbig in namelist table

      integer, parameter :: iddtsml=4
                       ! Unique index of dtsml in namelist table

      integer, parameter :: idgwmopt=48
                       ! Unique index of gwmopt in namelist table

      integer, parameter :: idimpopt=29
                       ! Unique index of impopt in namelist table

      integer, parameter :: idadvopt=38
                       ! Unique index of advopt in namelist table

      integer, parameter :: idgsdeps=450
                       ! Unique index of gsdeps in namelist table

      integer, parameter :: idweicoe=64
                       ! Unique index of weicoe in namelist table

      integer, parameter :: idfilcoe=39
                       ! Unique index of filcoe in namelist table

      integer, parameter :: iddtvcul=769
                       ! Unique index of dtvcul in namelist table

      integer, parameter :: idsmtopt=37
                       ! Unique index of smtopt in namelist table

      integer, parameter :: idsmhcoe=62
                       ! Unique index of smhcoe in namelist table

      integer, parameter :: idsmvcoe=63
                       ! Unique index of smvcoe in namelist table

      integer, parameter :: idsmqcoe=62
                       ! Unique index of smqcoe in namelist table

      integer, parameter :: idnlhcoe=443
                       ! Unique index of nlhcoe in namelist table

      integer, parameter :: idnlvcoe=444
                       ! Unique index of nlvcoe in namelist table

      integer, parameter :: idnlqcoe=443
                       ! Unique index of nlqcoe in namelist table

      integer, parameter :: idmfcopt=19
                       ! Unique index of mfcopt in namelist table

      integer, parameter :: idcoropt=20
                       ! Unique index of coropt in namelist table

      integer, parameter :: idcrvopt=166
                       ! Unique index of crvopt in namelist table

      integer, parameter :: idbuyopt=21
                       ! Unique index of buyopt in namelist table

      integer, parameter :: iddiaopt=43
                       ! Unique index of diaopt in namelist table

      integer, parameter :: iddivopt=24
                       ! Unique index of divopt in namelist table

      integer, parameter :: idcphopt=22
                       ! Unique index of cphopt in namelist table

      integer, parameter :: idhaiopt=208
                       ! Unique index of haiopt in namelist table

      integer, parameter :: idthresq=464
                       ! Unique index of thresq in namelist table

      integer, parameter :: iddtcmph=774
                       ! Unique index of dtcmph in namelist table

      integer, parameter :: idncbinw=36
                       ! Unique index of ncbinw in namelist table

      integer, parameter :: idbbinw=772
                       ! Unique index of bbinw in namelist table

      integer, parameter :: idsbinw=773
                       ! Unique index of sbinw in namelist table

      integer, parameter :: idqcgopt=3
                       ! Unique index of qcgopt in namelist table

      integer, parameter :: ideledlt=771
                       ! Unique index of eledlt in namelist table

      integer, parameter :: idaslopt=201
                       ! Unique index of aslopt in namelist table

      integer, parameter :: idtrkopt=173
                       ! Unique index of trkopt in namelist table

      integer, parameter :: idqt0opt=174
                       ! Unique index of qt0opt in namelist table

      integer, parameter :: idqt0num=175
                       ! Unique index of qt0num in namelist table

      integer, parameter :: idqt0=451
                       ! Unique index of qt0 in namelist table

      integer, parameter :: idqt0rx=452
                       ! Unique index of qt0rx in namelist table

      integer, parameter :: idqt0ry=453
                       ! Unique index of qt0ry in namelist table

      integer, parameter :: idqt0rz=454
                       ! Unique index of qt0rz in namelist table

      integer, parameter :: idqt0cx=455
                       ! Unique index of qt0cx in namelist table

      integer, parameter :: idqt0cy=456
                       ! Unique index of qt0cy in namelist table

      integer, parameter :: idqt0cz=457
                       ! Unique index of qt0cz in namelist table

      integer, parameter :: idqt0ds=458
                       ! Unique index of qt0ds in namelist table

      integer, parameter :: idqtdt=459
                       ! Unique index of qtdt in namelist table

      integer, parameter :: idqtstr=781
                       ! Unique index of qtstr in namelist table

      integer, parameter :: idqtend=782
                       ! Unique index of qtend in namelist table

      integer, parameter :: idtubopt=33
                       ! Unique index of tubopt in namelist table

      integer, parameter :: idisoopt=44
                       ! Unique index of isoopt in namelist table

      integer, parameter :: iddmpfmt=25
                       ! Unique index of dmpfmt in namelist table

      integer, parameter :: iddmplev=45
                       ! Unique index of dmplev in namelist table

      integer, parameter :: iddmpmon=180
                       ! Unique index of dmpmon in namelist table

      integer, parameter :: iddmpvar=15
                       ! Unique index of dmpvar in namelist table

      integer, parameter :: iddmpitv=40
                       ! Unique index of dmpitv in namelist table

      integer, parameter :: idmonitv=76
                       ! Unique index of monitv in namelist table

      integer, parameter :: idresopt=176
                       ! Unique index of resopt in namelist table

      integer, parameter :: idresitv=41
                       ! Unique index of resitv in namelist table

      integer, parameter :: idmxnopt=177
                       ! Unique index of mxnopt in namelist table

      integer, parameter :: idmxnvar=18
                       ! Unique index of mxnvar in namelist table

      integer, parameter :: idmxnitv=42
                       ! Unique index of mxnitv in namelist table

      integer, parameter :: idmpopt_gpv=26
                       ! Unique index of mpopt_gpv in namelist table

      integer, parameter :: idnspol_gpv=27
                       ! Unique index of nspol_gpv in namelist table

      integer, parameter :: idtlat1_gpv=51
                       ! Unique index of tlat1_gpv in namelist table

      integer, parameter :: idtlat2_gpv=52
                       ! Unique index of tlat2_gpv in namelist table

      integer, parameter :: idtlon_gpv=53
                       ! Unique index of tlon_gpv in namelist table

      integer, parameter :: idxdim_gpv=187
                       ! Unique index of xdim_gpv in namelist table

      integer, parameter :: idydim_gpv=188
                       ! Unique index of ydim_gpv in namelist table

      integer, parameter :: idzdim_gpv=189
                       ! Unique index of zdim_gpv in namelist table

      integer, parameter :: iddx_gpv=43
                       ! Unique index of dx_gpv in namelist table

      integer, parameter :: iddy_gpv=44
                       ! Unique index of dy_gpv in namelist table

      integer, parameter :: iddxiv_gpv=45
                       ! Unique index of dxiv_gpv in namelist table

      integer, parameter :: iddyiv_gpv=46
                       ! Unique index of dyiv_gpv in namelist table

      integer, parameter :: idulat_gpv=47
                       ! Unique index of ulat_gpv in namelist table

      integer, parameter :: idulon_gpv=48
                       ! Unique index of ulon_gpv in namelist table

      integer, parameter :: idriu_gpv=49
                       ! Unique index of riu_gpv in namelist table

      integer, parameter :: idrju_gpv=50
                       ! Unique index of rju_gpv in namelist table

      integer, parameter :: idintopt_gpv=28
                       ! Unique index of intopt_gpv in namelist table

      integer, parameter :: idrotopt_gpv=15
                       ! Unique index of rotopt_gpv in namelist table

      integer, parameter :: iddatype_gpv=3
                       ! Unique index of datype_gpv in namelist table

      integer, parameter :: idrefsfc_gpv=35
                       ! Unique index of refsfc_gpv in namelist table

      integer, parameter :: idetrvar_gpv=16
                       ! Unique index of etrvar_gpv in namelist table

      integer, parameter :: idmpopt_asl=202
                       ! Unique index of mpopt_asl in namelist table

      integer, parameter :: idnspol_asl=203
                       ! Unique index of nspol_asl in namelist table

      integer, parameter :: idtlat1_asl=785
                       ! Unique index of tlat1_asl in namelist table

      integer, parameter :: idtlat2_asl=786
                       ! Unique index of tlat2_asl in namelist table

      integer, parameter :: idtlon_asl=787
                       ! Unique index of tlon_asl in namelist table

      integer, parameter :: idxdim_asl=204
                       ! Unique index of xdim_asl in namelist table

      integer, parameter :: idydim_asl=205
                       ! Unique index of ydim_asl in namelist table

      integer, parameter :: idzdim_asl=206
                       ! Unique index of zdim_asl in namelist table

      integer, parameter :: iddx_asl=788
                       ! Unique index of dx_asl in namelist table

      integer, parameter :: iddy_asl=789
                       ! Unique index of dy_asl in namelist table

      integer, parameter :: iddxiv_asl=790
                       ! Unique index of dxiv_asl in namelist table

      integer, parameter :: iddyiv_asl=791
                       ! Unique index of dyiv_asl in namelist table

      integer, parameter :: idulat_asl=792
                       ! Unique index of ulat_asl in namelist table

      integer, parameter :: idulon_asl=793
                       ! Unique index of ulon_asl in namelist table

      integer, parameter :: idriu_asl=794
                       ! Unique index of riu_asl in namelist table

      integer, parameter :: idrju_asl=795
                       ! Unique index of rju_asl in namelist table

      integer, parameter :: idintopt_asl=207
                       ! Unique index of intopt_asl in namelist table

      integer, parameter :: idmpopt_rdr=163
                       ! Unique index of mpopt_rdr in namelist table

      integer, parameter :: idnspol_rdr=164
                       ! Unique index of nspol_rdr in namelist table

      integer, parameter :: idtlat1_rdr=422
                       ! Unique index of tlat1_rdr in namelist table

      integer, parameter :: idtlat2_rdr=423
                       ! Unique index of tlat2_rdr in namelist table

      integer, parameter :: idtlon_rdr=424
                       ! Unique index of tlon_rdr in namelist table

      integer, parameter :: idxdim_rdr=190
                       ! Unique index of xdim_rdr in namelist table

      integer, parameter :: idydim_rdr=191
                       ! Unique index of ydim_rdr in namelist table

      integer, parameter :: idzdim_rdr=192
                       ! Unique index of zdim_rdr in namelist table

      integer, parameter :: iddx_rdr=425
                       ! Unique index of dx_rdr in namelist table

      integer, parameter :: iddy_rdr=426
                       ! Unique index of dy_rdr in namelist table

      integer, parameter :: iddxiv_rdr=427
                       ! Unique index of dxiv_rdr in namelist table

      integer, parameter :: iddyiv_rdr=428
                       ! Unique index of dyiv_rdr in namelist table

      integer, parameter :: idulat_rdr=429
                       ! Unique index of ulat_rdr in namelist table

      integer, parameter :: idulon_rdr=430
                       ! Unique index of ulon_rdr in namelist table

      integer, parameter :: idriu_rdr=113
                       ! Unique index of riu_rdr in namelist table

      integer, parameter :: idrju_rdr=115
                       ! Unique index of rju_rdr in namelist table

      integer, parameter :: idrotopt_rdr=4
                       ! Unique index of rotopt_rdr in namelist table

      integer, parameter :: iddatype_rdr=8
                       ! Unique index of datype_rdr in namelist table

      integer, parameter :: idrdrcoe_rdr=70
                       ! Unique index of rdrcoe_rdr in namelist table

      integer, parameter :: idrdrexp_rdr=431
                       ! Unique index of rdrexp_rdr in namelist table

      integer, parameter :: idmpopt_trn=49
                       ! Unique index of mpopt_trn in namelist table

      integer, parameter :: idnspol_trn=50
                       ! Unique index of nspol_trn in namelist table

      integer, parameter :: idtlat1_trn=79
                       ! Unique index of tlat1_trn in namelist table

      integer, parameter :: idtlat2_trn=80
                       ! Unique index of tlat2_trn in namelist table

      integer, parameter :: idtlon_trn=81
                       ! Unique index of tlon_trn in namelist table

      integer, parameter :: idxdim_trn=193
                       ! Unique index of xdim_trn in namelist table

      integer, parameter :: idydim_trn=194
                       ! Unique index of ydim_trn in namelist table

      integer, parameter :: iddx_trn=82
                       ! Unique index of dx_trn in namelist table

      integer, parameter :: iddy_trn=83
                       ! Unique index of dy_trn in namelist table

      integer, parameter :: iddxiv_trn=84
                       ! Unique index of dxiv_trn in namelist table

      integer, parameter :: iddyiv_trn=85
                       ! Unique index of dyiv_trn in namelist table

      integer, parameter :: idulat_trn=86
                       ! Unique index of ulat_trn in namelist table

      integer, parameter :: idulon_trn=87
                       ! Unique index of ulon_trn in namelist table

      integer, parameter :: idriu_trn=88
                       ! Unique index of riu_trn in namelist table

      integer, parameter :: idrju_trn=89
                       ! Unique index of rju_trn in namelist table

      integer, parameter :: idintopt_trn=51
                       ! Unique index of intopt_trn in namelist table

      integer, parameter :: idmpopt_lnd=52
                       ! Unique index of mpopt_lnd in namelist table

      integer, parameter :: idnspol_lnd=53
                       ! Unique index of nspol_lnd in namelist table

      integer, parameter :: idtlat1_lnd=90
                       ! Unique index of tlat1_lnd in namelist table

      integer, parameter :: idtlat2_lnd=91
                       ! Unique index of tlat2_lnd in namelist table

      integer, parameter :: idtlon_lnd=92
                       ! Unique index of tlon_lnd in namelist table

      integer, parameter :: idxdim_lnd=195
                       ! Unique index of xdim_lnd in namelist table

      integer, parameter :: idydim_lnd=196
                       ! Unique index of ydim_lnd in namelist table

      integer, parameter :: iddx_lnd=93
                       ! Unique index of dx_lnd in namelist table

      integer, parameter :: iddy_lnd=94
                       ! Unique index of dy_lnd in namelist table

      integer, parameter :: iddxiv_lnd=95
                       ! Unique index of dxiv_lnd in namelist table

      integer, parameter :: iddyiv_lnd=96
                       ! Unique index of dyiv_lnd in namelist table

      integer, parameter :: idulat_lnd=97
                       ! Unique index of ulat_lnd in namelist table

      integer, parameter :: idulon_lnd=98
                       ! Unique index of ulon_lnd in namelist table

      integer, parameter :: idriu_lnd=99
                       ! Unique index of riu_lnd in namelist table

      integer, parameter :: idrju_lnd=100
                       ! Unique index of rju_lnd in namelist table

      integer, parameter :: idintopt_lnd=54
                       ! Unique index of intopt_lnd in namelist table

      integer, parameter :: idnumctg_lnd=55
                       ! Unique index of numctg_lnd in namelist table

      integer, parameter :: idlnduse_lnd=62
                       ! Unique index of lnduse_lnd in namelist table

      integer, parameter :: idalbe_lnd=216
                       ! Unique index of albe_lnd in namelist table

      integer, parameter :: idbeta_lnd=116
                       ! Unique index of beta_lnd in namelist table

      integer, parameter :: idz0m_lnd=316
                       ! Unique index of z0m_lnd in namelist table

      integer, parameter :: idz0h_lnd=466
                       ! Unique index of z0h_lnd in namelist table

      integer, parameter :: idcap_lnd=566
                       ! Unique index of cap_lnd in namelist table

      integer, parameter :: idnuu_lnd=666
                       ! Unique index of nuu_lnd in namelist table

      integer, parameter :: idmpopt_sst=56
                       ! Unique index of mpopt_sst in namelist table

      integer, parameter :: idnspol_sst=57
                       ! Unique index of nspol_sst in namelist table

      integer, parameter :: idtlat1_sst=101
                       ! Unique index of tlat1_sst in namelist table

      integer, parameter :: idtlat2_sst=102
                       ! Unique index of tlat2_sst in namelist table

      integer, parameter :: idtlon_sst=103
                       ! Unique index of tlon_sst in namelist table

      integer, parameter :: idxdim_sst=197
                       ! Unique index of xdim_sst in namelist table

      integer, parameter :: idydim_sst=198
                       ! Unique index of ydim_sst in namelist table

      integer, parameter :: iddx_sst=104
                       ! Unique index of dx_sst in namelist table

      integer, parameter :: iddy_sst=105
                       ! Unique index of dy_sst in namelist table

      integer, parameter :: iddxiv_sst=106
                       ! Unique index of dxiv_sst in namelist table

      integer, parameter :: iddyiv_sst=107
                       ! Unique index of dyiv_sst in namelist table

      integer, parameter :: idulat_sst=108
                       ! Unique index of ulat_sst in namelist table

      integer, parameter :: idulon_sst=109
                       ! Unique index of ulon_sst in namelist table

      integer, parameter :: idriu_sst=110
                       ! Unique index of riu_sst in namelist table

      integer, parameter :: idrju_sst=111
                       ! Unique index of rju_sst in namelist table

      integer, parameter :: idmpopt_ice=167
                       ! Unique index of mpopt_ice in namelist table

      integer, parameter :: idnspol_ice=168
                       ! Unique index of nspol_ice in namelist table

      integer, parameter :: idtlat1_ice=432
                       ! Unique index of tlat1_ice in namelist table

      integer, parameter :: idtlat2_ice=433
                       ! Unique index of tlat2_ice in namelist table

      integer, parameter :: idtlon_ice=434
                       ! Unique index of tlon_ice in namelist table

      integer, parameter :: idxdim_ice=199
                       ! Unique index of xdim_ice in namelist table

      integer, parameter :: idydim_ice=200
                       ! Unique index of ydim_ice in namelist table

      integer, parameter :: iddx_ice=435
                       ! Unique index of dx_ice in namelist table

      integer, parameter :: iddy_ice=436
                       ! Unique index of dy_ice in namelist table

      integer, parameter :: iddxiv_ice=437
                       ! Unique index of dxiv_ice in namelist table

      integer, parameter :: iddyiv_ice=438
                       ! Unique index of dyiv_ice in namelist table

      integer, parameter :: idulat_ice=439
                       ! Unique index of ulat_ice in namelist table

      integer, parameter :: idulon_ice=440
                       ! Unique index of ulon_ice in namelist table

      integer, parameter :: idriu_ice=441
                       ! Unique index of riu_ice in namelist table

      integer, parameter :: idrju_ice=442
                       ! Unique index of rju_ice in namelist table

      integer, parameter :: idfltyp_uni=21
                       ! Unique index of fltyp_uni in namelist table

      integer, parameter :: idflitv_uni=797
                       ! Unique index of flitv_uni in namelist table

      integer, parameter :: idrmopt_uni=216
                       ! Unique index of rmopt_uni in namelist table

      integer, parameter :: iduniopt_uni=217
                       ! Unique index of uniopt_uni in namelist table

      integer, parameter :: idugroup_uni=218
                       ! Unique index of ugroup_uni in namelist table

      integer, parameter :: idxsub_rst=214
                       ! Unique index of xsub_rst in namelist table

      integer, parameter :: idysub_rst=215
                       ! Unique index of ysub_rst in namelist table

      integer, parameter :: idflitv_rst=798
                       ! Unique index of flitv_rst in namelist table

      integer, parameter :: idrmopt_rst=219
                       ! Unique index of rmopt_rst in namelist table

      integer, parameter :: idiwest=39
                       ! Unique index of iwest in namelist table

      integer, parameter :: idieast=40
                       ! Unique index of ieast in namelist table

      integer, parameter :: idjsouth=41
                       ! Unique index of jsouth in namelist table

      integer, parameter :: idjnorth=42
                       ! Unique index of jnorth in namelist table

      integer, parameter :: idzeropt=0
                       ! Unique index of zeropt in namelist table

      integer, parameter :: idoneopt=-1
                       ! Unique index of oneopt in namelist table

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

      end module m_comindx
