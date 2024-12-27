!***********************************************************************
      module m_setref
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/11/05
!     Modification: 2004/10/12, 2004/12/17, 2006/03/06, 2006/08/08,
!                   2006/09/30, 2007/01/20, 2007/04/11, 2007/05/14,
!                   2007/07/30, 2007/11/26, 2008/03/12, 2008/05/02,
!                   2008/07/01, 2008/08/25, 2008/10/10, 2009/01/30,
!                   2009/02/27, 2009/03/12, 2009/08/20, 2009/11/05,
!                   2011/06/01, 2011/08/18, 2011/09/22, 2011/11/10

!     Author      : Satoki Tsujino
!     Modification: 2024/12/25 

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the referenced table.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comcapt
      use m_comdays
      use m_commath
      use m_comtable

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: setref, s_setref

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface setref

        module procedure s_setref

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_setref
!***********************************************************************

!-----7--------------------------------------------------------------7--

! Set the table for the number of elapse of days form the start of year.

      ela(0)=0
      ela(1)=31
      ela(2)=59
      ela(3)=90
      ela(4)=120
      ela(5)=151
      ela(6)=181
      ela(7)=212
      ela(8)=243
      ela(9)=273
      ela(10)=304
      ela(11)=334
      ela(12)=365

! -----

! Set the table for the number of elapse of days form the start of
! intercalary year.

      elaitc(0)=0
      elaitc(1)=31
      elaitc(2)=60
      elaitc(3)=91
      elaitc(4)=121
      elaitc(5)=152
      elaitc(6)=182
      elaitc(7)=213
      elaitc(8)=244
      elaitc(9)=274
      elaitc(10)=305
      elaitc(11)=335
      elaitc(12)=366

! -----

! Set the table for the number of remainder of days to the end of year.

      rem(1)=365
      rem(2)=334
      rem(3)=306
      rem(4)=275
      rem(5)=245
      rem(6)=214
      rem(7)=184
      rem(8)=153
      rem(9)=122
      rem(10)=92
      rem(11)=61
      rem(12)=31

! -----

! Set the table for the number of remainder of days to the end of
! intercalary year.

      remitc(1)=366
      remitc(2)=335
      remitc(3)=306
      remitc(4)=275
      remitc(5)=245
      remitc(6)=214
      remitc(7)=184
      remitc(8)=153
      remitc(9)=122
      remitc(10)=92
      remitc(11)=61
      remitc(12)=31

! -----

! Set the standard referenced cloud cover.

      rcdl(0)=0.e0
      rcdl(1)=0.e0
      rcdl(2)=0.e0
      rcdl(3)=0.e0
      rcdl(4)=0.e0
      rcdl(5)=0.e0
      rcdl(6)=0.e0
      rcdl(7)=0.e0
      rcdl(8)=0.e0
      rcdl(9)=0.e0
      rcdl(10)=0.e0
      rcdl(11)=0.e0
      rcdl(12)=0.e0
      rcdl(13)=0.e0
      rcdl(14)=0.e0
      rcdl(15)=0.e0
      rcdl(16)=0.e0
      rcdl(17)=0.e0
      rcdl(18)=0.e0
      rcdl(19)=0.e0
      rcdl(20)=0.e0
      rcdl(21)=0.e0
      rcdl(22)=0.e0
      rcdl(23)=0.e0
      rcdl(24)=0.e0
      rcdl(25)=0.e0
      rcdl(26)=0.e0
      rcdl(27)=0.e0
      rcdl(28)=0.e0
      rcdl(29)=0.e0
      rcdl(30)=0.e0
      rcdl(31)=0.e0
      rcdl(32)=0.e0
      rcdl(33)=0.e0
      rcdl(34)=0.e0
      rcdl(35)=0.e0
      rcdl(36)=0.e0
      rcdl(37)=0.e0
      rcdl(38)=0.e0
      rcdl(39)=0.e0
      rcdl(40)=0.e0
      rcdl(41)=0.e0
      rcdl(42)=0.e0
      rcdl(43)=0.e0
      rcdl(44)=0.e0
      rcdl(45)=0.e0
      rcdl(46)=0.e0
      rcdl(47)=0.e0
      rcdl(48)=0.e0
      rcdl(49)=0.e0
      rcdl(50)=0.e0
      rcdl(51)=0.e0
      rcdl(52)=0.e0
      rcdl(53)=0.e0
      rcdl(54)=0.e0
      rcdl(55)=0.e0
      rcdl(56)=0.e0
      rcdl(57)=0.e0
      rcdl(58)=0.e0
      rcdl(59)=0.e0
      rcdl(60)=0.e0
      rcdl(61)=0.e0
      rcdl(62)=0.e0
      rcdl(63)=0.e0
      rcdl(64)=0.e0
      rcdl(65)=0.e0
      rcdl(66)=.014e0
      rcdl(67)=.028e0
      rcdl(68)=.042e0
      rcdl(69)=.056e0
      rcdl(70)=.07e0
      rcdl(71)=.078e0
      rcdl(72)=.086e0
      rcdl(73)=.094e0
      rcdl(74)=.102e0
      rcdl(75)=.11e0
      rcdl(76)=.126e0
      rcdl(77)=.142e0
      rcdl(78)=.158e0
      rcdl(79)=.174e0
      rcdl(80)=.19e0
      rcdl(81)=.232e0
      rcdl(82)=.274e0
      rcdl(83)=.316e0
      rcdl(84)=.358e0
      rcdl(85)=.4e0
      rcdl(86)=.49e0
      rcdl(87)=.58e0
      rcdl(88)=.67e0
      rcdl(89)=.76e0
      rcdl(90)=.85e0
      rcdl(91)=.88e0
      rcdl(92)=.91e0
      rcdl(93)=.94e0
      rcdl(94)=.97e0
      rcdl(95)=1.e0
      rcdl(96)=1.e0
      rcdl(97)=1.e0
      rcdl(98)=1.e0
      rcdl(99)=1.e0
      rcdl(100)=1.e0
      rcdl(101)=1.e0

      rcdm(0)=0.e0
      rcdm(1)=0.e0
      rcdm(2)=0.e0
      rcdm(3)=0.e0
      rcdm(4)=0.e0
      rcdm(5)=0.e0
      rcdm(6)=0.e0
      rcdm(7)=0.e0
      rcdm(8)=0.e0
      rcdm(9)=0.e0
      rcdm(10)=0.e0
      rcdm(11)=0.e0
      rcdm(12)=0.e0
      rcdm(13)=0.e0
      rcdm(14)=0.e0
      rcdm(15)=0.e0
      rcdm(16)=0.e0
      rcdm(17)=0.e0
      rcdm(18)=0.e0
      rcdm(19)=0.e0
      rcdm(20)=0.e0
      rcdm(21)=0.e0
      rcdm(22)=0.e0
      rcdm(23)=0.e0
      rcdm(24)=0.e0
      rcdm(25)=0.e0
      rcdm(26)=0.e0
      rcdm(27)=0.e0
      rcdm(28)=0.e0
      rcdm(29)=0.e0
      rcdm(30)=0.e0
      rcdm(31)=0.e0
      rcdm(32)=0.e0
      rcdm(33)=0.e0
      rcdm(34)=0.e0
      rcdm(35)=0.e0
      rcdm(36)=0.e0
      rcdm(37)=0.e0
      rcdm(38)=0.e0
      rcdm(39)=0.e0
      rcdm(40)=0.e0
      rcdm(41)=0.e0
      rcdm(42)=0.e0
      rcdm(43)=0.e0
      rcdm(44)=0.e0
      rcdm(45)=0.e0
      rcdm(46)=0.e0
      rcdm(47)=0.e0
      rcdm(48)=0.e0
      rcdm(49)=0.e0
      rcdm(50)=0.e0
      rcdm(51)=.01e0
      rcdm(52)=.02e0
      rcdm(53)=.03e0
      rcdm(54)=.04e0
      rcdm(55)=.05e0
      rcdm(56)=.064e0
      rcdm(57)=.078e0
      rcdm(58)=.092e0
      rcdm(59)=.106e0
      rcdm(60)=.12e0
      rcdm(61)=.156e0
      rcdm(62)=.192e0
      rcdm(63)=.228e0
      rcdm(64)=.264e0
      rcdm(65)=.3e0
      rcdm(66)=.32e0
      rcdm(67)=.34e0
      rcdm(68)=.36e0
      rcdm(69)=.38e0
      rcdm(70)=.4e0
      rcdm(71)=.42e0
      rcdm(72)=.44e0
      rcdm(73)=.46e0
      rcdm(74)=.48e0
      rcdm(75)=.5e0
      rcdm(76)=.54e0
      rcdm(77)=.58e0
      rcdm(78)=.62e0
      rcdm(79)=.66e0
      rcdm(80)=.7e0
      rcdm(81)=.75e0
      rcdm(82)=.8e0
      rcdm(83)=.85e0
      rcdm(84)=.9e0
      rcdm(85)=.95e0
      rcdm(86)=.96e0
      rcdm(87)=.97e0
      rcdm(88)=.98e0
      rcdm(89)=.99e0
      rcdm(90)=1.e0
      rcdm(91)=1.e0
      rcdm(92)=1.e0
      rcdm(93)=1.e0
      rcdm(94)=1.e0
      rcdm(95)=1.e0
      rcdm(96)=1.e0
      rcdm(97)=1.e0
      rcdm(98)=1.e0
      rcdm(99)=1.e0
      rcdm(100)=1.e0
      rcdm(101)=1.e0

      rcdh(0)=0.e0
      rcdh(1)=0.e0
      rcdh(2)=0.e0
      rcdh(3)=0.e0
      rcdh(4)=0.e0
      rcdh(5)=0.e0
      rcdh(6)=0.e0
      rcdh(7)=0.e0
      rcdh(8)=0.e0
      rcdh(9)=0.e0
      rcdh(10)=0.e0
      rcdh(11)=0.e0
      rcdh(12)=0.e0
      rcdh(13)=0.e0
      rcdh(14)=0.e0
      rcdh(15)=0.e0
      rcdh(16)=0.e0
      rcdh(17)=0.e0
      rcdh(18)=0.e0
      rcdh(19)=0.e0
      rcdh(20)=0.e0
      rcdh(21)=0.e0
      rcdh(22)=0.e0
      rcdh(23)=0.e0
      rcdh(24)=0.e0
      rcdh(25)=0.e0
      rcdh(26)=0.e0
      rcdh(27)=0.e0
      rcdh(28)=0.e0
      rcdh(29)=0.e0
      rcdh(30)=0.e0
      rcdh(31)=0.e0
      rcdh(32)=0.e0
      rcdh(33)=0.e0
      rcdh(34)=0.e0
      rcdh(35)=0.e0
      rcdh(36)=0.e0
      rcdh(37)=0.e0
      rcdh(38)=0.e0
      rcdh(39)=0.e0
      rcdh(40)=0.e0
      rcdh(41)=.01e0
      rcdh(42)=.02e0
      rcdh(43)=.03e0
      rcdh(44)=.04e0
      rcdh(45)=.05e0
      rcdh(46)=.07e0
      rcdh(47)=.09e0
      rcdh(48)=.11e0
      rcdh(49)=.13e0
      rcdh(50)=.15e0
      rcdh(51)=.18e0
      rcdh(52)=.21e0
      rcdh(53)=.24e0
      rcdh(54)=.27e0
      rcdh(55)=.3e0
      rcdh(56)=.35e0
      rcdh(57)=.4e0
      rcdh(58)=.45e0
      rcdh(59)=.5e0
      rcdh(60)=.55e0
      rcdh(61)=.59e0
      rcdh(62)=.63e0
      rcdh(63)=.67e0
      rcdh(64)=.71e0
      rcdh(65)=.75e0
      rcdh(66)=.774e0
      rcdh(67)=.798e0
      rcdh(68)=.822e0
      rcdh(69)=.846e0
      rcdh(70)=.87e0
      rcdh(71)=.886e0
      rcdh(72)=.902e0
      rcdh(73)=.918e0
      rcdh(74)=.934e0
      rcdh(75)=.95e0
      rcdh(76)=.96e0
      rcdh(77)=.97e0
      rcdh(78)=.98e0
      rcdh(79)=.99e0
      rcdh(80)=1.e0
      rcdh(81)=1.e0
      rcdh(82)=1.e0
      rcdh(83)=1.e0
      rcdh(84)=1.e0
      rcdh(85)=1.e0
      rcdh(86)=1.e0
      rcdh(87)=1.e0
      rcdh(88)=1.e0
      rcdh(89)=1.e0
      rcdh(90)=1.e0
      rcdh(91)=1.e0
      rcdh(92)=1.e0
      rcdh(93)=1.e0
      rcdh(94)=1.e0
      rcdh(95)=1.e0
      rcdh(96)=1.e0
      rcdh(97)=1.e0
      rcdh(98)=1.e0
      rcdh(99)=1.e0
      rcdh(100)=1.e0
      rcdh(101)=1.e0

!ORIG cdrat(1)=.273872e0
!ORIG cdrat(2)=.294858e0
!ORIG cdrat(3)=.431270e0

      cdrat(1)=.4375e0
      cdrat(2)=.3750e0
      cdrat(3)=.1875e0

! -----

! Set the Koenig temperature dependent parameters.

      ckoe(0)=0.e0
      ckoe(1)=.7939e-10
      ckoe(2)=.7841e-9
      ckoe(3)=.3369e-8
      ckoe(4)=.4336e-8
      ckoe(5)=.5285e-8
      ckoe(6)=.3728e-8
      ckoe(7)=.1852e-8
      ckoe(8)=.2991e-9
      ckoe(9)=.4248e-9
      ckoe(10)=.7434e-9
      ckoe(11)=.1812e-8
      ckoe(12)=.4394e-8
      ckoe(13)=.9145e-8
      ckoe(14)=.1725e-7
      ckoe(15)=.3348e-7
      ckoe(16)=.1725e-7
      ckoe(17)=.9175e-8
      ckoe(18)=.4412e-8
      ckoe(19)=.2252e-8
      ckoe(20)=.9115e-9
      ckoe(21)=.4876e-9
      ckoe(22)=.3473e-9
      ckoe(23)=.4758e-9
      ckoe(24)=.6306e-9
      ckoe(25)=.8573e-9
      ckoe(26)=.7868e-9
      ckoe(27)=.7192e-9
      ckoe(28)=.6513e-9
      ckoe(29)=.5956e-9
      ckoe(30)=.5333e-9
      ckoe(31)=.4834e-9

      ckoe(32)=.4335e-9
      ckoe(33)=.3836e-9
      ckoe(34)=.3337e-9
      ckoe(35)=.2838e-9
      ckoe(36)=.2339e-9
      ckoe(37)=.1840e-9
      ckoe(38)=.1341e-9
      ckoe(39)=.8420e-10
      ckoe(40)=.3430e-10

      pkoe(0)=0.e0
      pkoe(1)=.4006e0
      pkoe(2)=.4831e0
      pkoe(3)=.5320e0
      pkoe(4)=.5307e0
      pkoe(5)=.5319e0
      pkoe(6)=.5249e0
      pkoe(7)=.4888e0
      pkoe(8)=.3894e0
      pkoe(9)=.4047e0
      pkoe(10)=.4318e0
      pkoe(11)=.4771e0
      pkoe(12)=.5183e0
      pkoe(13)=.5463e0
      pkoe(14)=.5651e0
      pkoe(15)=.5813e0
      pkoe(16)=.5655e0
      pkoe(17)=.5478e0
      pkoe(18)=.5203e0
      pkoe(19)=.4906e0
      pkoe(20)=.4447e0
      pkoe(21)=.4126e0
      pkoe(22)=.3960e0
      pkoe(23)=.4149e0
      pkoe(24)=.4320e0
      pkoe(25)=.4506e0
      pkoe(26)=.4483e0
      pkoe(27)=.4460e0
      pkoe(28)=.4433e0
      pkoe(29)=.4413e0
      pkoe(30)=.4382e0
      pkoe(31)=.4361e0

      pkoe(32)=.4340e0
      pkoe(33)=.4319e0
      pkoe(34)=.4298e0
      pkoe(35)=.4277e0
      pkoe(36)=.4256e0
      pkoe(37)=.4235e0
      pkoe(38)=.4214e0
      pkoe(39)=.4193e0
      pkoe(40)=.4172e0

      ckoe(0)=ckoe(0)*exp(3.e0*pkoe(0)*ln10)
      ckoe(1)=ckoe(1)*exp(3.e0*pkoe(1)*ln10)
      ckoe(2)=ckoe(2)*exp(3.e0*pkoe(2)*ln10)
      ckoe(3)=ckoe(3)*exp(3.e0*pkoe(3)*ln10)
      ckoe(4)=ckoe(4)*exp(3.e0*pkoe(4)*ln10)
      ckoe(5)=ckoe(5)*exp(3.e0*pkoe(5)*ln10)
      ckoe(6)=ckoe(6)*exp(3.e0*pkoe(6)*ln10)
      ckoe(7)=ckoe(7)*exp(3.e0*pkoe(7)*ln10)
      ckoe(8)=ckoe(8)*exp(3.e0*pkoe(8)*ln10)
      ckoe(9)=ckoe(9)*exp(3.e0*pkoe(9)*ln10)
      ckoe(10)=ckoe(10)*exp(3.e0*pkoe(10)*ln10)
      ckoe(11)=ckoe(11)*exp(3.e0*pkoe(11)*ln10)
      ckoe(12)=ckoe(12)*exp(3.e0*pkoe(12)*ln10)
      ckoe(13)=ckoe(13)*exp(3.e0*pkoe(13)*ln10)
      ckoe(14)=ckoe(14)*exp(3.e0*pkoe(14)*ln10)
      ckoe(15)=ckoe(15)*exp(3.e0*pkoe(15)*ln10)
      ckoe(16)=ckoe(16)*exp(3.e0*pkoe(16)*ln10)
      ckoe(17)=ckoe(17)*exp(3.e0*pkoe(17)*ln10)
      ckoe(18)=ckoe(18)*exp(3.e0*pkoe(18)*ln10)
      ckoe(19)=ckoe(19)*exp(3.e0*pkoe(19)*ln10)
      ckoe(20)=ckoe(20)*exp(3.e0*pkoe(20)*ln10)
      ckoe(21)=ckoe(21)*exp(3.e0*pkoe(21)*ln10)
      ckoe(22)=ckoe(22)*exp(3.e0*pkoe(22)*ln10)
      ckoe(23)=ckoe(23)*exp(3.e0*pkoe(23)*ln10)
      ckoe(24)=ckoe(24)*exp(3.e0*pkoe(24)*ln10)
      ckoe(25)=ckoe(25)*exp(3.e0*pkoe(25)*ln10)
      ckoe(26)=ckoe(26)*exp(3.e0*pkoe(26)*ln10)
      ckoe(27)=ckoe(27)*exp(3.e0*pkoe(27)*ln10)
      ckoe(28)=ckoe(28)*exp(3.e0*pkoe(28)*ln10)
      ckoe(29)=ckoe(29)*exp(3.e0*pkoe(29)*ln10)
      ckoe(30)=ckoe(30)*exp(3.e0*pkoe(30)*ln10)
      ckoe(31)=ckoe(31)*exp(3.e0*pkoe(31)*ln10)
      ckoe(32)=ckoe(32)*exp(3.e0*pkoe(32)*ln10)
      ckoe(33)=ckoe(33)*exp(3.e0*pkoe(33)*ln10)
      ckoe(34)=ckoe(34)*exp(3.e0*pkoe(34)*ln10)
      ckoe(35)=ckoe(35)*exp(3.e0*pkoe(35)*ln10)
      ckoe(36)=ckoe(36)*exp(3.e0*pkoe(36)*ln10)
      ckoe(37)=ckoe(37)*exp(3.e0*pkoe(37)*ln10)
      ckoe(38)=ckoe(38)*exp(3.e0*pkoe(38)*ln10)
      ckoe(39)=ckoe(39)*exp(3.e0*pkoe(39)*ln10)
      ckoe(40)=ckoe(40)*exp(3.e0*pkoe(40)*ln10)

! -----

! Set the standard referenced charging rate.

      refch(0,0)=0.e0
      refch(1,0)=0.e0
      refch(2,0)=0.e0
      refch(3,0)=0.e0
      refch(4,0)=0.e0
      refch(5,0)=0.e0
      refch(6,0)=0.e0
      refch(7,0)=0.e0
      refch(8,0)=0.e0
      refch(9,0)=0.e0
      refch(10,0)=0.e0
      refch(11,0)=0.e0

      refch(0,1)=0.e0
      refch(1,1)=.33e0
      refch(2,1)=3.3e0
      refch(3,1)=3.33e0
      refch(4,1)=3.33e0
      refch(5,1)=3.33e0
      refch(6,1)=3.33e0
      refch(7,1)=3.33e0
      refch(8,1)=3.33e0
      refch(9,1)=3.33e0
      refch(10,1)=3.33e0
      refch(11,1)=.33e0

      refch(0,2)=0.e0
      refch(1,2)=3.33e0
      refch(2,2)=13.3e0
      refch(3,2)=13.3e0
      refch(4,2)=13.3e0
      refch(5,2)=13.3e0
      refch(6,2)=13.3e0
      refch(7,2)=13.3e0
      refch(8,2)=13.3e0
      refch(9,2)=13.3e0
      refch(10,2)=10.e0
      refch(11,2)=3.33e0

      refch(0,3)=0.e0
      refch(1,3)=15.e0
      refch(2,3)=20.e0
      refch(3,3)=20.e0
      refch(4,3)=20.e0
      refch(5,3)=20.e0
      refch(6,3)=20.e0
      refch(7,3)=20.e0
      refch(8,3)=20.e0
      refch(9,3)=13.3e0
      refch(10,3)=3.33e0
      refch(11,3)=3.33e0

      refch(0,4)=0.e0
      refch(1,4)=20.e0
      refch(2,4)=26.7e0
      refch(3,4)=26.7e0
      refch(4,4)=23.4e0
      refch(5,4)=22.e0
      refch(6,4)=21.e0
      refch(7,4)=20.e0
      refch(8,4)=20.e0
      refch(9,4)=3.33e0
      refch(10,4)=-.8e0
      refch(11,4)=-.7e0

      refch(0,5)=0.e0
      refch(1,5)=30.e0
      refch(2,5)=33.4e0
      refch(3,5)=30.e0
      refch(4,5)=26.7e0
      refch(5,5)=23.4e0
      refch(6,5)=20.e0
      refch(7,5)=16.7e0
      refch(8,5)=5.84e0
      refch(9,5)=-1.11e0
      refch(10,5)=-1.65e0
      refch(11,5)=-1.4e0

      refch(0,6)=0.e0
      refch(2,6)=33.4e0
      refch(1,6)=33.4e0
      refch(3,6)=33.4e0
      refch(4,6)=30.e0
      refch(5,6)=33.4e0
      refch(6,6)=16.7e0
      refch(7,6)=5.84e0
      refch(8,6)=-1.11e0
      refch(9,6)=-2.22e0
      refch(10,6)=-2.5e0
      refch(11,6)=-2.2e0

      refch(0,7)=0.e0
      refch(1,7)=33.4e0
      refch(2,7)=33.4e0
      refch(3,7)=33.4e0
      refch(4,7)=30.e0
      refch(5,7)=20.e0
      refch(6,7)=8.35e0
      refch(7,7)=1.66e0
      refch(8,7)=-2.22e0
      refch(9,7)=-3.33e0
      refch(10,7)=-3.33e0
      refch(11,7)=-3.e0

      refch(0,8)=0.e0
      refch(1,8)=33.4e0
      refch(2,8)=33.4e0
      refch(3,8)=30.e0
      refch(4,8)=26.7e0
      refch(5,8)=16.7e0
      refch(6,8)=3.33e0
      refch(7,8)=-1.66e0
      refch(8,8)=-3.33e0
      refch(9,8)=-6.67e0
      refch(10,8)=-5.93e0
      refch(11,8)=-3.33e0

      refch(0,9)=0.e0
      refch(1,9)=33.4e0
      refch(2,9)=33.4e0
      refch(3,9)=26.7e0
      refch(4,9)=16.7e0
      refch(5,9)=8.35e0
      refch(6,9)=-3.33e0
      refch(7,9)=-3.33e0
      refch(8,9)=-7.83e0
      refch(9,9)=-10.e0
      refch(10,9)=-8.53e0
      refch(11,9)=-5.84e0

      refch(0,10)=0.e0
      refch(1,10)=33.4e0
      refch(2,10)=33.4e0
      refch(3,10)=16.7e0
      refch(4,10)=8.35e0
      refch(5,10)=-.33e0
      refch(6,10)=-5.84e0
      refch(7,10)=-10.e0
      refch(8,10)=-12.3e0
      refch(9,10)=-13.4e0
      refch(10,10)=-11.2e0
      refch(11,10)=-8.35e0

      refch(0,11)=0.e0
      refch(1,11)=33.4e0
      refch(2,11)=26.7e0
      refch(3,11)=8.35e0
      refch(4,11)=-8.35e0
      refch(5,11)=-13.3e0
      refch(6,11)=-16.7e0
      refch(7,11)=-16.7e0
      refch(8,11)=-16.7e0
      refch(9,11)=-16.7e0
      refch(10,11)=-14.1e0
      refch(11,11)=-10.e0

      refch(0,12)=0.e0
      refch(1,12)=30.e0
      refch(2,12)=16.7e0
      refch(3,12)=13.3e0
      refch(4,12)=3.33e0
      refch(5,12)=-5.84e0
      refch(6,12)=-8.35e0
      refch(7,12)=-13.3e0
      refch(8,12)=-16.7e0
      refch(9,12)=-16.7e0
      refch(10,12)=-16.7e0
      refch(11,12)=-16.7e0

      refch(0,13)=0.e0
      refch(1,13)=20.e0
      refch(2,13)=16.7e0
      refch(3,13)=16.7e0
      refch(4,13)=16.7e0
      refch(5,13)=16.7e0
      refch(6,13)=8.35e0
      refch(7,13)=5.84e0
      refch(8,13)=-.33e0
      refch(9,13)=-5.84e0
      refch(10,13)=-8.35e0
      refch(11,13)=-8.35e0

      refch(0,14)=0.e0
      refch(1,14)=16.7e0
      refch(2,14)=16.7e0
      refch(3,14)=16.7e0
      refch(4,14)=16.7e0
      refch(5,14)=16.7e0
      refch(6,14)=16.7e0
      refch(7,14)=13.3e0
      refch(8,14)=8.35e0
      refch(9,14)=3.33e0
      refch(10,14)=3.33e0
      refch(11,14)=3.33e0

      refch(0,15)=0.e0
      refch(1,15)=16.7e0
      refch(2,15)=16.7e0
      refch(3,15)=16.7e0
      refch(4,15)=16.7e0
      refch(5,15)=16.7e0
      refch(6,15)=16.7e0
      refch(7,15)=16.7e0
      refch(8,15)=16.7e0
      refch(9,15)=13.3e0
      refch(10,15)=8.35e0
      refch(11,15)=8.35e0

      refch(0,16)=0.e0
      refch(1,16)=16.7e0
      refch(2,16)=16.7e0
      refch(3,16)=16.7e0
      refch(4,16)=16.7e0
      refch(5,16)=16.7e0
      refch(6,16)=16.7e0
      refch(7,16)=16.7e0
      refch(8,16)=16.7e0
      refch(9,16)=16.7e0
      refch(10,16)=16.7e0
      refch(11,16)=13.3e0

      refch(0,17)=0.e0
      refch(1,17)=16.7e0
      refch(2,17)=16.7e0
      refch(3,17)=16.7e0
      refch(4,17)=16.7e0
      refch(5,17)=16.7e0
      refch(6,17)=16.7e0
      refch(7,17)=16.7e0
      refch(8,17)=16.7e0
      refch(9,17)=16.7e0
      refch(10,17)=16.7e0
      refch(11,17)=16.7e0

! -----

! Set the cloud condensation nuclei.

      nccn(1)=6.3e1
      nccn(2)=8.23e1
      nccn(3)=9.09e1
      nccn(4)=1.2e2
      nccn(5)=1.58e2

! -----

! Set the standard coalescence efficiency between water bins.

      rrdbw(1)=6.e-4
      rrdbw(2)=8.e-4
      rrdbw(3)=10.e-4
      rrdbw(4)=15.e-4
      rrdbw(5)=20.e-4
      rrdbw(6)=25.e-4
      rrdbw(7)=30.e-4
      rrdbw(8)=40.e-4
      rrdbw(9)=50.e-4
      rrdbw(10)=60.e-4
      rrdbw(11)=70.e-4
      rrdbw(12)=100.e-4
      rrdbw(13)=150.e-4
      rrdbw(14)=200.e-4
      rrdbw(15)=300.e-4

      rrcbw(1)=0.e0
      rrcbw(2)=.05e0
      rrcbw(3)=.1e0
      rrcbw(4)=.15e0
      rrcbw(5)=.2e0
      rrcbw(6)=.25e0
      rrcbw(7)=.3e0
      rrcbw(8)=.35e0
      rrcbw(9)=.4e0
      rrcbw(10)=.45e0
      rrcbw(11)=.5e0
      rrcbw(12)=.55e0
      rrcbw(13)=.6e0
      rrcbw(14)=.65e0
      rrcbw(15)=.7e0
      rrcbw(16)=.75e0
      rrcbw(17)=.8e0
      rrcbw(18)=.85e0
      rrcbw(19)=.9e0
      rrcbw(20)=.95e0
      rrcbw(21)=1.e0

      rewbw(1,1)=.001e0
      rewbw(2,1)=.001e0
      rewbw(3,1)=.001e0
      rewbw(4,1)=.001e0
      rewbw(5,1)=.001e0
      rewbw(6,1)=.001e0
      rewbw(7,1)=.001e0
      rewbw(8,1)=.001e0
      rewbw(9,1)=.001e0
      rewbw(10,1)=.001e0
      rewbw(11,1)=.001e0
      rewbw(12,1)=.001e0
      rewbw(13,1)=.001e0
      rewbw(14,1)=.001e0
      rewbw(15,1)=.001e0

      rewbw(1,2)=.003e0
      rewbw(2,2)=.003e0
      rewbw(3,2)=.003e0
      rewbw(4,2)=.004e0
      rewbw(5,2)=.005e0
      rewbw(6,2)=.005e0
      rewbw(7,2)=.005e0
      rewbw(8,2)=.01e0
      rewbw(9,2)=.1e0
      rewbw(10,2)=.05e0
      rewbw(11,2)=.2e0
      rewbw(12,2)=.5e0
      rewbw(13,2)=.77e0
      rewbw(14,2)=.87e0
      rewbw(15,2)=.97e0

      rewbw(1,3)=.007e0
      rewbw(2,3)=.007e0
      rewbw(3,3)=.007e0
      rewbw(4,3)=.008e0
      rewbw(5,3)=.009e0
      rewbw(6,3)=.01e0
      rewbw(7,3)=.01e0
      rewbw(8,3)=.07e0
      rewbw(9,3)=.4e0
      rewbw(10,3)=.43e0
      rewbw(11,3)=.58e0
      rewbw(12,3)=.79e0
      rewbw(13,3)=.93e0
      rewbw(14,3)=.96e0
      rewbw(15,3)=1.e0

      rewbw(1,4)=.009e0
      rewbw(2,4)=.009e0
      rewbw(3,4)=.009e0
      rewbw(4,4)=.012e0
      rewbw(5,4)=.015e0
      rewbw(6,4)=.01e0
      rewbw(7,4)=.02e0
      rewbw(8,4)=.28e0
      rewbw(9,4)=.6e0
      rewbw(10,4)=.64e0
      rewbw(11,4)=.75e0
      rewbw(12,4)=.91e0
      rewbw(13,4)=.97e0
      rewbw(14,4)=.98e0
      rewbw(15,4)=1.e0

      rewbw(1,5)=.014e0
      rewbw(2,5)=.014e0
      rewbw(3,5)=.014e0
      rewbw(4,5)=.015e0
      rewbw(5,5)=.016e0
      rewbw(6,5)=.03e0
      rewbw(7,5)=.06e0
      rewbw(8,5)=.5e0
      rewbw(9,5)=.7e0
      rewbw(10,5)=.77e0
      rewbw(11,5)=.84e0
      rewbw(12,5)=.95e0
      rewbw(13,5)=.97e0
      rewbw(14,5)=1.e0
      rewbw(15,5)=1.e0

      rewbw(1,6)=.017e0
      rewbw(2,6)=.017e0
      rewbw(3,6)=.017e0
      rewbw(4,6)=.02e0
      rewbw(5,6)=.022e0
      rewbw(6,6)=.06e0
      rewbw(7,6)=.1e0
      rewbw(8,6)=.62e0
      rewbw(9,6)=.78e0
      rewbw(10,6)=.84e0
      rewbw(11,6)=.88e0
      rewbw(12,6)=.95e0
      rewbw(13,6)=1.e0
      rewbw(14,6)=1.e0
      rewbw(15,6)=1.e0

      rewbw(1,7)=.03e0
      rewbw(2,7)=.03e0
      rewbw(3,7)=.024e0
      rewbw(4,7)=.022e0
      rewbw(5,7)=.032e0
      rewbw(6,7)=.062e0
      rewbw(7,7)=.2e0
      rewbw(8,7)=.68e0
      rewbw(9,7)=.83e0
      rewbw(10,7)=.87e0
      rewbw(11,7)=.9e0
      rewbw(12,7)=.95e0
      rewbw(13,7)=1.e0
      rewbw(14,7)=1.e0
      rewbw(15,7)=1.e0

      rewbw(1,8)=.025e0
      rewbw(2,8)=.025e0
      rewbw(3,8)=.025e0
      rewbw(4,8)=.036e0
      rewbw(5,8)=.043e0
      rewbw(6,8)=.13e0
      rewbw(7,8)=.27e0
      rewbw(8,8)=.74e0
      rewbw(9,8)=.86e0
      rewbw(10,8)=.89e0
      rewbw(11,8)=.92e0
      rewbw(12,8)=1.e0
      rewbw(13,8)=1.e0
      rewbw(14,8)=1.e0
      rewbw(15,8)=1.e0

      rewbw(1,9)=.027e0
      rewbw(2,9)=.027e0
      rewbw(3,9)=.027e0
      rewbw(4,9)=.04e0
      rewbw(5,9)=.052e0
      rewbw(6,9)=.2e0
      rewbw(7,9)=.4e0
      rewbw(8,9)=.78e0
      rewbw(9,9)=.88e0
      rewbw(10,9)=.9e0
      rewbw(11,9)=.94e0
      rewbw(12,9)=1.e0
      rewbw(13,9)=1.e0
      rewbw(14,9)=1.e0
      rewbw(15,9)=1.e0

      rewbw(1,10)=.03e0
      rewbw(2,10)=.03e0
      rewbw(3,10)=.03e0
      rewbw(4,10)=.047e0
      rewbw(5,10)=.064e0
      rewbw(6,10)=.25e0
      rewbw(7,10)=.5e0
      rewbw(8,10)=.8e0
      rewbw(9,10)=.9e0
      rewbw(10,10)=.91e0
      rewbw(11,10)=.95e0
      rewbw(12,10)=1.e0
      rewbw(13,10)=1.e0
      rewbw(14,10)=1.e0
      rewbw(15,10)=1.e0

      rewbw(1,11)=.04e0
      rewbw(2,11)=.04e0
      rewbw(3,11)=.033e0
      rewbw(4,11)=.037e0
      rewbw(5,11)=.068e0
      rewbw(6,11)=.24e0
      rewbw(7,11)=.55e0
      rewbw(8,11)=.8e0
      rewbw(9,11)=.9e0
      rewbw(10,11)=.91e0
      rewbw(11,11)=.95e0
      rewbw(12,11)=1.e0
      rewbw(13,11)=1.e0
      rewbw(14,11)=1.e0
      rewbw(15,11)=1.e0

      rewbw(1,12)=.035e0
      rewbw(2,12)=.035e0
      rewbw(3,12)=.035e0
      rewbw(4,12)=.055e0
      rewbw(5,12)=.079e0
      rewbw(6,12)=.29e0
      rewbw(7,12)=.58e0
      rewbw(8,12)=.8e0
      rewbw(9,12)=.9e0
      rewbw(10,12)=.91e0
      rewbw(11,12)=.95e0
      rewbw(12,12)=1.e0
      rewbw(13,12)=1.e0
      rewbw(14,12)=1.e0
      rewbw(15,12)=1.e0

      rewbw(1,13)=.037e0
      rewbw(2,13)=.037e0
      rewbw(3,13)=.037e0
      rewbw(4,13)=.062e0
      rewbw(5,13)=.082e0
      rewbw(6,13)=.29e0
      rewbw(7,13)=.59e0
      rewbw(8,13)=.78e0
      rewbw(9,13)=.9e0
      rewbw(10,13)=.91e0
      rewbw(11,13)=.95e0
      rewbw(12,13)=1.e0
      rewbw(13,13)=1.e0
      rewbw(14,13)=1.e0
      rewbw(15,13)=1.e0

      rewbw(1,14)=.037e0
      rewbw(2,14)=.037e0
      rewbw(3,14)=.037e0
      rewbw(4,14)=.06e0
      rewbw(5,14)=.08e0
      rewbw(6,14)=.29e0
      rewbw(7,14)=.58e0
      rewbw(8,14)=.77e0
      rewbw(9,14)=.89e0
      rewbw(10,14)=.91e0
      rewbw(11,14)=.95e0
      rewbw(12,14)=1.e0
      rewbw(13,14)=1.e0
      rewbw(14,14)=1.e0
      rewbw(15,14)=1.e0

      rewbw(1,15)=.037e0
      rewbw(2,15)=.037e0
      rewbw(3,15)=.037e0
      rewbw(4,15)=.041e0
      rewbw(5,15)=.075e0
      rewbw(6,15)=.25e0
      rewbw(7,15)=.54e0
      rewbw(8,15)=.76e0
      rewbw(9,15)=.88e0
      rewbw(10,15)=.92e0
      rewbw(11,15)=.95e0
      rewbw(12,15)=1.e0
      rewbw(13,15)=1.e0
      rewbw(14,15)=1.e0
      rewbw(15,15)=1.e0

      rewbw(1,16)=.037e0
      rewbw(2,16)=.037e0
      rewbw(3,16)=.037e0
      rewbw(4,16)=.052e0
      rewbw(5,16)=.067e0
      rewbw(6,16)=.25e0
      rewbw(7,16)=.51e0
      rewbw(8,16)=.77e0
      rewbw(9,16)=.88e0
      rewbw(10,16)=.93e0
      rewbw(11,16)=.97e0
      rewbw(12,16)=1.e0
      rewbw(13,16)=1.e0
      rewbw(14,16)=1.e0
      rewbw(15,16)=1.e0

      rewbw(1,17)=.037e0
      rewbw(2,17)=.037e0
      rewbw(3,17)=.037e0
      rewbw(4,17)=.047e0
      rewbw(5,17)=.057e0
      rewbw(6,17)=.25e0
      rewbw(7,17)=.49e0
      rewbw(8,17)=.77e0
      rewbw(9,17)=.89e0
      rewbw(10,17)=.95e0
      rewbw(11,17)=1.e0
      rewbw(12,17)=1.e0
      rewbw(13,17)=1.e0
      rewbw(14,17)=1.e0
      rewbw(15,17)=1.e0

      rewbw(1,18)=.036e0
      rewbw(2,18)=.036e0
      rewbw(3,18)=.036e0
      rewbw(4,18)=.042e0
      rewbw(5,18)=.048e0
      rewbw(6,18)=.23e0
      rewbw(7,18)=.47e0
      rewbw(8,18)=.78e0
      rewbw(9,18)=.92e0
      rewbw(10,18)=1.e0
      rewbw(11,18)=1.02e0
      rewbw(12,18)=1.02e0
      rewbw(13,18)=1.02e0
      rewbw(14,18)=1.02e0
      rewbw(15,18)=1.02e0

      rewbw(1,19)=.04e0
      rewbw(2,19)=.04e0
      rewbw(3,19)=.035e0
      rewbw(4,19)=.033e0
      rewbw(5,19)=.04e0
      rewbw(6,19)=.112e0
      rewbw(7,19)=.45e0
      rewbw(8,19)=.79e0
      rewbw(9,19)=1.01e0
      rewbw(10,19)=1.03e0
      rewbw(11,19)=1.04e0
      rewbw(12,19)=1.04e0
      rewbw(13,19)=1.04e0
      rewbw(14,19)=1.04e0
      rewbw(15,19)=1.04e0

      rewbw(1,20)=.033e0
      rewbw(2,20)=.033e0
      rewbw(3,20)=.033e0
      rewbw(4,20)=.033e0
      rewbw(5,20)=.033e0
      rewbw(6,20)=.119e0
      rewbw(7,20)=.47e0
      rewbw(8,20)=.95e0
      rewbw(9,20)=1.3e0
      rewbw(10,20)=1.7e0
      rewbw(11,20)=2.3e0
      rewbw(12,20)=2.3e0
      rewbw(13,20)=2.3e0
      rewbw(14,20)=2.3e0
      rewbw(15,20)=2.3e0

      rewbw(1,21)=.027e0
      rewbw(2,21)=.027e0
      rewbw(3,21)=.027e0
      rewbw(4,21)=.027e0
      rewbw(5,21)=.027e0
      rewbw(6,21)=.125e0
      rewbw(7,21)=.52e0
      rewbw(8,21)=1.4e0
      rewbw(9,21)=2.3e0
      rewbw(10,21)=3.e0
      rewbw(11,21)=4.e0
      rewbw(12,21)=4.e0
      rewbw(13,21)=4.e0
      rewbw(14,21)=4.e0
      rewbw(15,21)=4.e0

! -----

! Set the referenced caption for dumped variables.

      ncpt(1)=30

      write(capt(1)(1:60),'(a60)')                                      &
     &    'x components of velocity [m/s]                              '

      ncpt(2)=30

      write(capt(2)(1:60),'(a60)')                                      &
     &    'y components of velocity [m/s]                              '

      ncpt(3)=30

      write(capt(3)(1:60),'(a60)')                                      &
     &    'z components of velocity [m/s]                              '

      ncpt(4)=24

      write(capt(4)(1:60),'(a60)')                                      &
     &    'base state pressure [Pa]                                    '

      ncpt(5)=26

      write(capt(5)(1:60),'(a60)')                                      &
     &    'pressure perturbation [Pa]                                  '

      ncpt(6)=36

      write(capt(6)(1:60),'(a60)')                                      &
     &    'base state potential temperature [K]                        '

      ncpt(7)=38

      write(capt(7)(1:60),'(a60)')                                      &
     &    'potential temperature perturbation [K]                      '

      ncpt(8)=43

      write(capt(8)(1:60),'(a60)')                                      &
     &    'base state water vapor mixing ratio [kg/kg]                 '

      ncpt(9)=32

      write(capt(9)(1:60),'(a60)')                                      &
     &    'water vapor mixing ratio [kg/kg]                            '

      ncpt(10)=32

      write(capt(10)(1:60),'(a60)')                                     &
     &    'cloud water mixing ratio [kg/kg]                            '

      ncpt(11)=31

      write(capt(11)(1:60),'(a60)')                                     &
     &    'rain water mixing ratio [kg/kg]                             '

      ncpt(12)=20

      write(capt(12)(1:60),'(a60)')                                     &
     &    'rain fall rate [m/s]                                        '

      ncpt(13)=25

      write(capt(13)(1:60),'(a60)')                                     &
     &    'accumulated rain fall [m]                                   '

      ncpt(14)=30

      write(capt(14)(1:60),'(a60)')                                     &
     &    'cloud ice mixing ratio [kg/kg]                              '

      ncpt(15)=25

      write(capt(15)(1:60),'(a60)')                                     &
     &    'snow mixing ratio [kg/kg]                                   '

      ncpt(16)=28

      write(capt(16)(1:60),'(a60)')                                     &
     &    'graupel mixing ratio [kg/kg]                                '

      ncpt(17)=20

      write(capt(17)(1:60),'(a60)')                                     &
     &    'snow fall rate [m/s]                                        '

      ncpt(18)=25

      write(capt(18)(1:60),'(a60)')                                     &
     &    'accumulated snow fall [m]                                   '

      ncpt(19)=23

      write(capt(19)(1:60),'(a60)')                                     &
     &    'graupel fall rate [m/s]                                     '

      ncpt(20)=28

      write(capt(20)(1:60),'(a60)')                                     &
     &    'accumulated graupel fall [m]                                '

      ncpt(21)=31

      write(capt(21)(1:60),'(a60)')                                     &
     &    'cloud ice concentrations [1/kg]                             '

      ncpt(22)=26

      write(capt(22)(1:60),'(a60)')                                     &
     &    'snow concentrations [1/kg]                                  '

      ncpt(23)=29

      write(capt(23)(1:60),'(a60)')                                     &
     &    'graupel concentrations [1/kg]                               '

      ncpt(24)=26

      write(capt(24)(1:60),'(a60)')                                     &
     &    'z physical coordinates [m]                                  '

      ncpt(25)=31

      write(capt(25)(1:60),'(a60)')                                     &
     &    'turbulent kinetic energy [J/kg]                             '

      ncpt(26)=13

      write(capt(26)(1:60),'(a60)')                                     &
     &    'pressure [Pa]                                               '

      ncpt(27)=25

      write(capt(27)(1:60),'(a60)')                                     &
     &    'potential temperature [K]                                   '

      ncpt(28)=41

      write(capt(28)(1:60),'(a60)')                                     &
     &    'base state x components of velocity [m/s]                   '

      ncpt(29)=41

      write(capt(29)(1:60),'(a60)')                                     &
     &    'base state y components of velocity [m/s]                   '

      ncpt(30)=20

      write(capt(30)(1:60),'(a60)')                                     &
     &    'zonal velocity [m/s]                                        '

      ncpt(31)=25

      write(capt(31)(1:60),'(a60)')                                     &
     &    'meridional velocity [m/s]                                   '

      ncpt(32)=31

      write(capt(32)(1:60),'(a60)')                                     &
     &    'base state zonal velocity [m/s]                             '

      ncpt(33)=36

      write(capt(33)(1:60),'(a60)')                                     &
     &    'base state meridional velocity [m/s]                        '

      ncpt(34)=19

      write(capt(34)(1:60),'(a60)')                                     &
     &    'tracer mixing ratio                                         '

      ncpt(35)=32

      write(capt(35)(1:60),'(a60)')                                     &
     &    'total water mixing ratio [kg/kg]                            '

      ncpt(36)=33

      write(capt(36)(1:60),'(a60)')                                     &
     &    'total water concentrations [1/kg]                           '

      ncpt(37)=27

      write(capt(37)(1:60),'(a60)')                                     &
     &    'total water fall rate [m/s]                                 '

      ncpt(38)=32

      write(capt(38)(1:60),'(a60)')                                     &
     &    'accumulated total water fall [m]                            '

      ncpt(39)=30

      write(capt(39)(1:60),'(a60)')                                     &
     &    'total ice mixing ratio [kg/kg]                              '

      ncpt(40)=31

      write(capt(40)(1:60),'(a60)')                                     &
     &    'total ice concentrations [1/kg]                             '

      ncpt(41)=25

      write(capt(41)(1:60),'(a60)')                                     &
     &    'total ice fall rate [m/s]                                   '

      ncpt(42)=30

      write(capt(42)(1:60),'(a60)')                                     &
     &    'accumulated total ice fall [m]                              '

      ncpt(43)=41

      write(capt(43)(1:60),'(a60)')                                     &
     &    'cloud water charging distribution [fC/kg]                   '

      ncpt(44)=40

      write(capt(44)(1:60),'(a60)')                                     &
     &    'rain water charging distribution [fC/kg]                    '

      ncpt(45)=39

      write(capt(45)(1:60),'(a60)')                                     &
     &    'cloud ice charging distribution [fC/kg]                     '

      ncpt(46)=34

      write(capt(46)(1:60),'(a60)')                                     &
     &    'snow charging distribution [fC/kg]                          '

      ncpt(47)=37

      write(capt(47)(1:60),'(a60)')                                     &
     &    'graupel charging distribution [fC/kg]                       '

      ncpt(48)=52

      write(capt(48)(1:60),'(a60)')                                     &
     &    'x components of velocity at an altitude of 10m [m/s]        '

      ncpt(49)=52

      write(capt(49)(1:60),'(a60)')                                     &
     &    'y components of velocity at an altitude of 10m [m/s]        '

      ncpt(50)=36

      write(capt(50)(1:60),'(a60)')                                     &
     &    'pressure at an altitude of 1.5m [Pa]                        '

      ncpt(51)=48

      write(capt(51)(1:60),'(a60)')                                     &
     &    'potential temperature at an altitude of 1.5m [K]            '

      ncpt(52)=55

      write(capt(52)(1:60),'(a60)')                                     &
     &    'water vapor mixing ratio at an altitude of 1.5m [kg/kg]     '

      ncpt(53)=42

      write(capt(53)(1:60),'(a60)')                                     &
     &    'zonal velocity at an altitude of 10m [m/s]                  '

      ncpt(54)=47

      write(capt(54)(1:60),'(a60)')                                     &
     &    'meridional velocity at an altitude of 10m [m/s]             '

      ncpt(55)=36

      write(capt(55)(1:60),'(a60)')                                     &
     &    'soil and sea surface temperature [K]                        '

      ncpt(56)=34

      write(capt(56)(1:60),'(a60)')                                     &
     &    'sensible heat over surface [W/m^2]                          '

      ncpt(57)=32

      write(capt(57)(1:60),'(a60)')                                     &
     &    'latent heat over surface [W/m^2]                            '

      ncpt(58)=41

      write(capt(58)(1:60),'(a60)')                                     &
     &    'net downward short wave radiation [W/m^2]                   '

      ncpt(59)=36

      write(capt(59)(1:60),'(a60)')                                     &
     &    'downward long wave radiation [W/m^2]                        '

      ncpt(60)=34

      write(capt(60)(1:60),'(a60)')                                     &
     &    'upward long wave radiation [W/m^2]                          '

      ncpt(61)=20

      write(capt(61)(1:60),'(a60)')                                     &
     &    'averaged cloud cover                                        '

      ncpt(62)=58

      write(capt(62)(1:60),'(a60)')                                     &
     &    'surface momentum flux for x components of velocity [N/m^2]  '

      ncpt(63)=58

      write(capt(63)(1:60),'(a60)')                                     &
     &    'surface momentum flux for y components of velocity [N/m^2]  '

      ncpt(64)=34

      write(capt(64)(1:60),'(a60)')                                     &
     &    'surface heat flux [(kg K)/(m^2 s)]                          '

      ncpt(65)=34

      write(capt(65)(1:60),'(a60)')                                     &
     &    'surface moisture flux [kg/(m^2 s)]                          '

      ncpt(66)=18

      write(capt(66)(1:60),'(a60)')                                     &
     &    'terrain height [m]                                          '

      ncpt(67)=17

      write(capt(67)(1:60),'(a60)')                                     &
     &    'latitude [degree]                                           '

      ncpt(68)=18

      write(capt(68)(1:60),'(a60)')                                     &
     &    'longitude [degree]                                          '

      ncpt(69)=16

      write(capt(69)(1:60),'(a60)')                                     &
     &    'map scale factor                                            '

      ncpt(70)=18

      write(capt(70)(1:60),'(a60)')                                     &
     &    'Coriolis parameter                                          '

      ncpt(71)=18

      write(capt(71)(1:60),'(a60)')                                     &
     &    'Coriolis parameter                                          '

      ncpt(72)=24

      write(capt(72)(1:60),'(a60)')                                     &
     &    'real land use categories                                    '

      ncpt(73)=41

      write(capt(73)(1:60),'(a60)')                                     &
     &    'maximum instantaneous wind velocity [m/s]                   '

      ncpt(74)=33

      write(capt(74)(1:60),'(a60)')                                     &
     &    'cloud water concentrations [1/kg]                           '

      ncpt(75)=32

      write(capt(75)(1:60),'(a60)')                                     &
     &    'rain water concentrations [1/kg]                            '

      ncpt(76)=30

      write(capt(76)(1:60),'(a60)')                                     &
     &    'number of water super droplets                              '

      ncpt(77)=28

      write(capt(77)(1:60),'(a60)')                                     &
     &    'number of ice super droplets                                '

      ncpt(78)=30

      write(capt(78)(1:60),'(a60)')                                     &
     &    'global solar radiation [W/m^2]                              '

      ncpt(79)=26

      write(capt(79)(1:60),'(a60)')                                     &
     &    'cloud cover in lower layer                                  '

      ncpt(80)=27

      write(capt(80)(1:60),'(a60)')                                     &
     &    'cloud cover in middle layer                                 '

      ncpt(81)=26

      write(capt(81)(1:60),'(a60)')                                     &
     &    'cloud cover in upper layer                                  '

      ncpt(82)=38

      write(capt(82)(1:60),'(a60)')                                     &
     &    'terminal velocity of cloud water [m/s]                      '

      ncpt(83)=37

      write(capt(83)(1:60),'(a60)')                                     &
     &    'terminal velocity of rain water [m/s]                       '

      ncpt(84)=36

      write(capt(84)(1:60),'(a60)')                                     &
     &    'terminal velocity of cloud ice [m/s]                        '

      ncpt(85)=31

      write(capt(85)(1:60),'(a60)')                                     &
     &    'terminal velocity of snow [m/s]                             '

      ncpt(86)=34

      write(capt(86)(1:60),'(a60)')                                     &
     &    'terminal velocity of graupel [m/s]                          '

      ncpt(87)=39

      write(capt(87)(1:60),'(a60)')                                     &
     &    'total mixing ratio of soil dust [kg/kg]                     '

      ncpt(88)=36

      write(capt(88)(1:60),'(a60)')                                     &
     &    'total mixing ratio of carbon [kg/kg]                        '

      ncpt(89)=37

      write(capt(89)(1:60),'(a60)')                                     &
     &    'total mixing ratio of sulfate [kg/kg]                       '

      ncpt(90)=38

      write(capt(90)(1:60),'(a60)')                                     &
     &    'total mixing ratio of sea salt [kg/kg]                      '

      ncpt(91)=43

      write(capt(91)(1:60),'(a60)')                                     &
     &    'total mixing ratio of wet soil dust [kg/kg]                 '

      ncpt(92)=40

      write(capt(92)(1:60),'(a60)')                                     &
     &    'total mixing ratio of wet carbon [kg/kg]                    '

      ncpt(93)=41

      write(capt(93)(1:60),'(a60)')                                     &
     &    'total mixing ratio of wet sulfate [kg/kg]                   '

      ncpt(94)=42

      write(capt(94)(1:60),'(a60)')                                     &
     &    'total mixing ratio of wet sea salt [kg/kg]                  '

      ncpt(95)=27

      write(capt(95)(1:60),'(a60)')                                     &
     &    'cloud water fall rate [m/s]                                 '

      ncpt(96)=32

      write(capt(96)(1:60),'(a60)')                                     &
     &    'accumulated cloud water fall [m]                            '

      ncpt(97)=25

      write(capt(97)(1:60),'(a60)')                                     &
     &    'cloud ice fall rate [m/s]                                   '

      ncpt(98)=30

      write(capt(98)(1:60),'(a60)')                                     &
     &    'accumulated cloud ice fall [m]                              '

      ncpt(99)=25

      write(capt(99)(1:60),'(a60)')                                     &
     &    'hail mixing ratio [kg/kg]                                   '

      ncpt(100)=20

      write(capt(100)(1:60),'(a60)')                                    &
     &    'hail fall rate [m/s]                                        '

      ncpt(101)=25

      write(capt(101)(1:60),'(a60)')                                    &
     &    'accumulated hail fall [m]                                   '

      ncpt(102)=26

      write(capt(102)(1:60),'(a60)')                                    &
     &    'hail concentrations [1/kg]                                  '

      ncpt(103)=34

      write(capt(103)(1:60),'(a60)')                                    &
     &    'hail charging distribution [fC/kg]                          '

      ncpt(104)=31

      write(capt(104)(1:60),'(a60)')                                    &
     &    'terminal velocity of hail [m/s]                             '

      ncpt(151)=51

      write(capt(151)(1:60),'(a60)')                                    &
     &    'total heating rate by cloud physics processes [K/s]         '

      ncpt(152)=57

      write(capt(152)(1:60),'(a60)')                                    &
     &    'accumulated heating rate by cloud physics processes [K/s]   '

      ncpt(153)=50

      write(capt(153)(1:60),'(a60)')                                    &
     &    'x component of frictional force for turbulence [N]          '

      ncpt(154)=50

      write(capt(154)(1:60),'(a60)')                                    &
     &    'y component of frictional force for turbulence [N]          '

      ncpt(155)=38

      write(capt(155)(1:60),'(a60)')                                    &
     &    'heating for turbulence [K/s]                                '

      ncpt(156)=37

      write(capt(156)(1:60),'(a60)')                                    &
     &    'moistening for turbulence [kg/(s m3)]                       '

      ncpt(157)=45

      write(capt(157)(1:60),'(a60)')                                    &
     &    'u component of change by Asselin Filter [m/s]               '

      ncpt(158)=45

      write(capt(158)(1:60),'(a60)')                                    &
     &    'v component of change by Asselin Filter [m/s]               '

      ncpt(159)=45

      write(capt(159)(1:60),'(a60)')                                    &
     &    'w component of change by Asselin Filter [m/s]               '

      ncpt(160)=41

      write(capt(160)(1:60),'(a60)')                                    &
     &    'pressure of change by Asselin Filter [Pa]                   '

      ncpt(161)=53

      write(capt(161)(1:60),'(a60)')                                    &
     &    'potential temperature of change by Asselin Filter [K]       '

      ncpt(162)=60

      write(capt(162)(1:60),'(a60)')                                    &
     &    'water vapor mixing ratio of change by Asselin Filter [kg/kg]'

      ncpt(163)=38

      write(capt(163)(1:60),'(a60)')                                    &
     &    'x component of numerical diffusion [N]                      '

      ncpt(164)=38

      write(capt(164)(1:60),'(a60)')                                    &
     &    'y component of numerical diffusion [N]                      '

      ncpt(165)=38

      write(capt(165)(1:60),'(a60)')                                    &
     &    'z component of numerical diffusion [N]                      '

      ncpt(166)=48

      write(capt(166)(1:60),'(a60)')                                    &
     &    'numerical diffusion of pressure [(kg Pa)/(m3 s)]            '

      ncpt(167)=60

      write(capt(167)(1:60),'(a60)')                                    &
     &    'numerical diffusion of potential temperature [K/s]          '

      ncpt(168)=59

      write(capt(168)(1:60),'(a60)')                                    &
     &    'numerical diffusion of water vapor mixing ratio [kg/(m3 s)] '

! -----

      end subroutine s_setref

!-----7--------------------------------------------------------------7--

      end module m_setref
