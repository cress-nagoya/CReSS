#!/bin/csh -f

################################################################################
#                                                                              #
#     C shell script, compile.csh                                              #
#                                                                              #
#     This script generates the executable files.                              #
#                                                                              #
#     Author      : Sakakibara Atsushi                                         #
#     Date        : 1999/01/20                                                 #
#     Modification: 1999/01/25, 1999/03/25, 1999/04/06, 1999/05/10,            #
#                   1999/05/20, 1999/06/07, 1999/06/21, 1999/06/28,            #
#                   1999/07/05, 1999/08/03, 1999/08/09, 1999/08/18,            #
#                   1999/09/06, 1999/09/30, 1999/10/07, 1999/10/12,            #
#                   1999/11/01, 1999/11/19, 1999/11/24, 2000/02/02,            #
#                   2000/02/07, 2000/03/17, 2000/03/23, 2000/04/18,            #
#                   2000/06/01, 2000/07/05, 2001/02/13, 2001/03/13,            #
#                   2001/04/15, 2001/06/29, 2002/04/09, 2002/07/03,            #
#                   2002/08/27, 2002/09/09, 2003/05/19, 2003/11/05,            #
#                   2004/03/05, 2004/04/01, 2004/04/15, 2004/05/31,            #
#                   2004/06/10, 2004/07/01, 2004/08/01, 2004/08/20,            #
#                   2004/09/25, 2004/12/17, 2005/01/07, 2005/02/10,            #
#                   2005/04/04, 2005/08/05, 2005/10/05, 2006/01/10,            #
#                   2006/02/13, 2006/04/03, 2006/05/12, 2006/06/21,            #
#                   2006/09/30, 2006/12/04, 2007/01/05, 2007/01/20,            #
#                   2007/04/11, 2007/05/07, 2007/05/14, 2007/06/27,            #
#                   2007/07/30, 2007/11/26, 2008/01/11, 2008/04/17,            #
#                   2008/05/02, 2008/06/09, 2008/07/01, 2008/08/19,            #
#                   2008/10/10, 2010/02/01, 2010/03/01, 2010/03/12,            #
#                   2010/05/17, 2010/12/02, 2010/12/06, 2010/12/08,            #
#                   2010/12/09, 2011/01/14, 2011/03/29, 2011/08/18,            #
#                   2011/09/22, 2011/11/10, 2011/12/17, 2012/01/25,            #
#                   2013/02/13, 2013/10/08                                     #
#                                                                              #
#     Author      : Hasegawa Koichi                                            #
#     Modification: 2014/06/10                                                 #
#                                                                              #
#     Author      : Masaya Kato                                                #
#     Modification: 2016/04/25                                                 #
#                                                                              #
#     Author      : Satoki Tsujino                                             #
#     Modification: 2017/06/14                                                 #
#                                                                              #
################################################################################

umask 027

unalias -a

set srcdir = Src

set prflow = none

set target = none
set conffl = none

set maxtarget = 0
set ncompara  = 1  ## set ncompara>=1

if($#argv == 1 || $#argv == 2) then

  goto start_command-search

  return_command-search:

  set curdir = `$pwd_`

  $ls_ ${curdir}/${srcdir}/Make_* >& /dev/null

  if($status == 0) then

    if($1 == 'all') then

      foreach makefl (${curdir}/${srcdir}/Make_*)

        set target = `$basename_ $makefl | $cut_ -c 6-`

        if(!($target == 'clean' || $target == 'object' ||\
             $target == 'module' || $target == 'include' ||\
             $target == 'radcln' || $target == 'crscln' ||\
             $target == 'spncln')) then

          @ maxtarget = $maxtarget + 1

        endif

      end

      if($maxtarget > 0) then

        set prflow = allcompile

      endif

    else

      foreach makefl (${curdir}/${srcdir}/Make_*)

        set target = `$basename_ $makefl | $cut_ -c 6-`

        if($target == $1) then

          if(!($target == 'clean' || $target == 'object' ||\
               $target == 'module' || $target == 'include'||\
               $target == 'radcln' || $target == 'crscln' ||\
               $target == 'spncln')) then

            @ maxtarget = $maxtarget + 1

            set prflow = onecompile

          else

            if($target == 'clean'||$target == 'radcln'||\
               $target == 'crscln'||$target == 'spncln') then

              set prflow = clean

            endif

          endif

          break

        endif

      end

    endif

  else

    if($1 == 'clean'||$1 == 'radcln'||$1 == 'crscln'||$1 == 'spncln') then

      set prflow = clean

    else

      echo "error: files, Make_* in directory, $srcdir not found."\
           "need to configure to construct directory, ${srcdir}."

      exit

    endif

  endif

echo "test" $prflow $target

  if($#argv == 2) then

    set conffl = $2

    if($conffl != 'compile.conf') then

      echo 'usage: compile.csh target_program_name [compile.conf]'
      echo 'usage: compile.csh clean'

      exit

    endif

    if($prflow == 'clean') then

      echo 'usage: compile.csh target_program_name [compile.conf]'
      echo 'usage: compile.csh clean'

      exit

    endif

  endif

  if($prflow == 'allcompile' || $prflow == 'onecompile') then

    set cnttarget = 0

    start_next:

    @ cnttarget = $cnttarget + 1

    if($prflow == 'allcompile') then

      set curtarget = 0

      foreach makefl (${curdir}/${srcdir}/Make_*)

        set target = `$basename_ $makefl | $cut_ -c 6-`

        if(!($target == 'clean' || $target == 'object' ||\
             $target == 'module' || $target == 'include'||\
             $target == 'radcln' || $target == 'crscln' ||\
             $target == 'spncln')) then

          @ curtarget = $curtarget + 1

          if($curtarget == $cnttarget) then

            break

          endif

        endif

      end

    else

      set target = $1

    endif

    goto start_compiler-search

    return_compiler-search:

    goto start_compile

    return_compile:

    if($cnttarget < $maxtarget) then

      goto start_next

    endif

    echo 'compilation completed.'

    exit

  else if($prflow == 'clean') then

    set target = $1

    goto start_cleaning

    return_cleaning:

    echo 'cleaning up completed.'

    exit

  else

    echo 'usage: compile.csh target_program_name [compile.conf]'
    echo 'usage: compile.csh clean'

    exit

  endif

else

  echo 'usage: compile.csh target_program_name [compile.conf]'
  echo 'usage: compile.csh clean'

  exit

endif

################################################################################
      start_command-search:
################################################################################

set cmdlst = (ls rm cp pwd basename cut diff make touch date ln)

set cmddir = (/dir/ls /dir/rm /dir/cp /dir/pwd /dir/basename /dir/cut\
              /dir/diff /dir/make /dir/touch /dir/date /dir/ln)

set cmdnum = 1
set pthnum = 1
set errnum = $#path

foreach getcmd ($cmdlst)

  while(!(-x ${path[${pthnum}]}/$getcmd))

    if($pthnum == $errnum) then

      echo "error: command, $getcmd not found."\
           'expecting command search path not correct.'

      exit

    endif

    @ pthnum = $pthnum + 1

  end

  set cmddir[${cmdnum}] = ${path[${pthnum}]}

  @ pthnum = 1

  @ cmdnum = $cmdnum + 1

end

set ls_       = ${cmddir[1]}/$cmdlst[1]
set rm_       = ${cmddir[2]}/$cmdlst[2]
set cp_       = ${cmddir[3]}/$cmdlst[3]
set pwd_      = ${cmddir[4]}/$cmdlst[4]
set basename_ = ${cmddir[5]}/$cmdlst[5]
set cut_      = ${cmddir[6]}/$cmdlst[6]
set diff_     = ${cmddir[7]}/$cmdlst[7]
set make_     = ${cmddir[8]}/$cmdlst[8]
set touch_    = ${cmddir[9]}/$cmdlst[9]
set date_     = ${cmddir[10]}/$cmdlst[10]
set ln_       = ${cmddir[11]}/$cmdlst[11]

goto return_command-search

################################################################################
      start_compiler-search:
################################################################################

set curdir = `$pwd_`

set comcmd = default

if($target == 'solver' || $target == 'radlib' || $target == 'spnlib') then

#  set fname = mpif90
 set fname = sxmpif90
  set cname = cc

else

#  set fname = f90
 set fname = sxf90
  set cname = cc

endif

if($#argv == 2 && -e ${curdir}/$conffl) then

  set comcmd = specified

else

  set pthnum = 1
  set errnum = $#path

  while(!(-x ${path[${pthnum}]}/$fname))

    if($pthnum == $errnum) then

      @ pthnum = $pthnum + 1

      break

    endif

    @ pthnum = $pthnum + 1

  end

  if($pthnum <= $errnum) then

    set fname = ${path[${pthnum}]}/$fname

  endif

  @ pthnum = 1

  while(!(-x ${path[${pthnum}]}/$cname))

    if($pthnum == $errnum) then

      @ pthnum = $pthnum + 1

      break

    endif

    @ pthnum = $pthnum + 1

  end

  if($pthnum <= $errnum) then

    set cname = ${path[${pthnum}]}/$cname

  endif

endif

goto return_compiler-search

################################################################################
      start_compile:
################################################################################

set curdir = `$pwd_`

set devdir = SrcDevelop

$ls_ ${curdir}/$srcdir >& /dev/null

if($status != 0)  then

  echo "error: directory, $srcdir not found."\
       "need to configure to construct directory, ${srcdir}."

  exit

endif

$ls_ ${curdir}/BuiltVersion=* >& /dev/null

if($status == 0) then

  set targetdir = `$ls_ ${curdir}/BuiltVersion=*`

  set targetdir = `echo $targetdir | $cut_ -f 1 -d ' '`

  set targetdir = `$basename_ $targetdir | $cut_ -c 14-`

  if($targetdir == 'ORIG') then

    set targetdir = 0rig

  endif

  $ls_ ${curdir}/${devdir}/${targetdir}/* >& /dev/null

  if($status == 0) then

    foreach devcnts (${curdir}/${devdir}/${targetdir}/*)

      $ls_ ${devcnts}/* >& /dev/null

      if($status == 0) then

        foreach devcntfl (${devcnts}/*)

          $diff_ $devcntfl ${curdir}/$srcdir >& /dev/null

          set swpstat = $status

          if($swpstat == 1) then

            $cp_ -p $devcntfl ${curdir}/$srcdir

            echo "cp -p $devcntfl ${curdir}/${srcdir}"

          else if($swpstat >= 2) then

            echo 'error: unexpected directory or file found.'\
                 "need to reconfigure to construct directory, ${srcdir}."

            exit

          endif

        end

      else

        if(-f $devcnts) then

          $diff_ $devcnts ${curdir}/$srcdir >& /dev/null

          set swpstat = $status

          if($swpstat == 1) then

            $cp_ -p $devcnts ${curdir}/$srcdir

            echo "cp -p $devcnts ${curdir}/${srcdir}"

          else if($swpstat >= 2) then

            echo 'error: unexpected file found.'\
                 "need to reconfigure to construct directory, ${srcdir}."

            exit

          endif

        endif

      endif

    end

  else

    echo "error: no content found in directory, ${devdir}/${targetdir}."\
         "need to check directory, ${devdir}/${targetdir}."

    exit

  endif

else

  $ls_ ${curdir}/${devdir}/* >& /dev/null

  if($status == 0) then

    foreach devcnts (${curdir}/${devdir}/*)

      $ls_ ${devcnts}/* >& /dev/null

      if($status == 0) then

        foreach devcntfl (${devcnts}/*)

          $diff_ $devcntfl ${curdir}/$srcdir >& /dev/null

          set swpstat = $status

          if($swpstat == 1) then

            $cp_ -p $devcntfl ${curdir}/$srcdir

            echo "cp -p $devcntfl ${curdir}/${srcdir}"

          else if($swpstat >= 2) then

            echo 'error: unexpected directory or file found.'\
                 "need to reconfigure to construct directory, ${srcdir}."

            exit

          endif

        end

      else

        if(-f $devcnts) then

          $diff_ $devcnts ${curdir}/$srcdir >& /dev/null

          set swpstat = $status

          if($swpstat == 1) then

            $cp_ -p $devcnts ${curdir}/$srcdir

            echo "cp -p $devcnts ${curdir}/${srcdir}"

          else if($swpstat >= 2) then

            echo 'error: unexpected file found.'\
                 "need to reconfigure to construct directory, ${srcdir}."

            exit

          endif

        endif

      endif

    end

  else

    echo "error: no content found in directory, ${devdir}."\
         "need to check directory, ${devdir}."

    exit

  endif

endif

cd ${curdir}/$srcdir

echo "cd ${curdir}/$srcdir"

if($comcmd == 'default') then

  set mkcmd = "FC=$fname CC=$cname TOUCH=$touch_ DATE=$date_"

else

  set mkcmd = "TOUCH=$touch_ DATE=$date_"

endif

set mkmcr = "TARGET=$target CURDIR=$curdir"

#$make_ $target $mkcmd $mkmcr
$make_ $target $mkcmd $mkmcr -j $ncompara 

if(-e ${curdir}/${target}.exe) then

  $rm_ -f ${curdir}/${target}.exe

  echo "rm -f ${curdir}/${target}.exe"

endif

if(-e ${curdir}/${srcdir}/${target}.exe) then

  $ln_ -fs ${curdir}/${srcdir}/${target}.exe ${curdir}/${target}.exe

  echo "ln -fs ${curdir}/${srcdir}/${target}.exe ${curdir}/${target}.exe"

endif

cd $curdir

echo "cd $curdir"

goto return_compile

################################################################################
      start_cleaning:
################################################################################

set curdir = `$pwd_`

$ls_ ${curdir}/*.exe >& /dev/null

if($status == 0) then

  if( $target != 'radcln' && $target != 'spncln' ) then
    $rm_ -f ${curdir}/*.exe

    echo "rm -f ${curdir}/*.exe"
  endif

else

  echo 'information: no executable file(s) to remove.'

endif

$ls_ ${curdir}/$srcdir >& /dev/null

if($status == 0) then

  cd ${curdir}/$srcdir

  echo "cd ${curdir}/$srcdir"

  set mkcmd = "RM=$rm_"
  set mkmcr = "TARGET=$target CURDIR=$curdir"
  if( $target == 'crscln' ) then
    set mkmcr = "TARGET=clean CURDIR=$curdir"
  endif

  $make_ $target $mkcmd $mkmcr

  cd $curdir

  echo "cd $curdir"

else

  echo "information: directory, $srcdir not found. nothing to remove."

endif

goto return_cleaning
