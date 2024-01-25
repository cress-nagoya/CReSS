#!/bin/csh -f

################################################################################
#                                                                              #
#     C shell script, configure.csh                                            #
#                                                                              #
#     This script build Src directory.                                         #
#                                                                              #
#     Author      : Sakakibara Atsushi                                         #
#     Date        : 2013/10/08                                                 #
#     Modification:                                                            #
#                                                                              #
################################################################################

umask 027

unalias -a

set srcdir = Src

set orgdir = SrcOrig
set devdir = SrcDevelop

set extdir = Extra

set target = none

if($#argv == 0 || $#argv == 1) then

  goto start_command-search

  return_command-search:

  if($#argv == 0) then

    set target = orig

    goto start_cleaning

    return_cleaning_orig:

    goto start_building

    return_building_orig:

    echo "configuration completed."

    exit

  else

    if($1 == 'clean') then

      set target = clean

    else

      set curdir = `$pwd_`

      $ls_ ${curdir}/${devdir}/* >& /dev/null

      if($status == 0) then

        foreach devcntdir (${curdir}/${devdir}/*)

          set devcntdir = `$basename_ ${devcntdir}`

          if($devcntdir == '0rig') then

             set target = orig

          else

             set target =\
                   `echo $devcntdir | $tr_ -d '[0-9].' | $tr_ '[A-Z]' '[a-z]'`

          endif

          if($target == $1) then

            break

          else

            set target = none

          endif

        end 

      endif

    endif

    if($target == 'none') then

      echo 'usage: configure.csh'
      echo 'usage: configure.csh target_name'
      echo 'usage: configure.csh clean'

      exit

    else if($target == 'clean') then

      goto start_cleaning

      return_cleaning_clean:

      echo "configuration completed."

      exit

    else

      goto start_cleaning

      return_cleaning_other:

      goto start_building

      return_building_other:

      echo "configuration completed."

      exit

    endif

  endif

else

  echo 'usage: configure.csh'
  echo 'usage: configure.csh target_name'
  echo 'usage: configure.csh clean'

  exit

endif

################################################################################
      start_command-search:
################################################################################

set cmdlst = (ls rm cp pwd basename cut tr mkdir touch)

set cmddir = (/dir/ls /dir/rm /dir/cp /dir/pwd /dir/basename /dir/cut\
              /dir/tr /dir/mkdir /dir/touch)

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
set tr_       = ${cmddir[7]}/$cmdlst[7]
set mkdir_    = ${cmddir[8]}/$cmdlst[8]
set touch_    = ${cmddir[9]}/$cmdlst[9]

goto return_command-search

################################################################################
      start_building:
################################################################################

set curdir = `$pwd_`

$mkdir_ ${curdir}/$srcdir

echo "mkdir ${curdir}/${srcdir}"

$ls_ ${curdir}/${orgdir}/* >& /dev/null

if($status == 0) then

  foreach dircheck (${curdir}/${orgdir}/*)

    if(-d $dircheck) then

      echo "error: unexpected directory found in directory, ${orgdir}."\
           "need to check directory, ${orgdir}."

      exit

    endif

  end

  $cp_ -p ${curdir}/${orgdir}/* ${curdir}/$srcdir

  echo "cp -p ${curdir}/${orgdir}/* ${curdir}/${srcdir}"

else

  echo 'error: original source codes not found.'\
       "need to check directory, ${orgdir}."

  exit

endif

$ls_ ${curdir}/${devdir}/* >& /dev/null

if($status == 0) then

  foreach devcntdir (`$ls_ -r ${curdir}/${devdir}`)

    set devcntdir = `$basename_ ${devcntdir}`

    if($devcntdir == '0rig') then

       set targetdir = orig

    else

       set targetdir =\
                `echo $devcntdir | $tr_ -d '[0-9].' | $tr_ '[A-Z]' '[a-z]'`

    endif

    if($targetdir == $target) then

      set targetdir = $devcntdir

      break

    else

      set targetdir = none

    endif

  end 

  $ls_ ${curdir}/${devdir}/${targetdir}/* >& /dev/null

  if($status == 0) then

    foreach devcnts (${curdir}/${devdir}/${targetdir}/*)

      $ls_ ${devcnts}/* >& /dev/null

      if($status == 0) then

        foreach dircheck (${devcnts}/*)

          if(-d $dircheck) then

            echo 'error: unexpected directory found in directory,'\
                 "${devdir}/${targetdir}. need to check directory,"\
                 "${devdir}/${targetdir}."

            exit

          endif

        end

        $cp_ -p ${devcnts}/* ${curdir}/$srcdir

        echo "cp -p ${devcnts}/* ${curdir}/${srcdir}"

      else

        if(-f $devcnts) then

          $cp_ -p $devcnts ${curdir}/$srcdir

          echo "cp -p $devcnts ${curdir}/${srcdir}"

        endif

      endif

    end

  else

    echo "error: no content found in directory, ${devdir}/${targetdir}."\
         "need to check directory, ${devdir}/${targetdir}."

    exit

  endif

else

  echo "error: no sub directory found in directory, ${devdir}."\
       "need to check directory, ${devdir}."

  exit

endif

$ls_ ${curdir}/${extdir}/${targetdir}/* >& /dev/null

if($status == 0) then

  foreach extcnts (${curdir}/${extdir}/${targetdir}/*)

    $cp_ -rp $extcnts $curdir

    echo "cp -rp $extcnts ${curdir}"

  end

else

  echo 'information: can not search any directories and files'\
       "in directory, ${extdir}/${targetdir}. do nothing."

endif

if($target == 'orig') then

  $touch_ ${curdir}/BuiltVersion=ORIG

  echo "touch ${curdir}/BuiltVersion=ORIG"

else

  $touch_ ${curdir}/BuiltVersion=$targetdir

  echo "touch ${curdir}/BuiltVersion=${targetdir}"

endif

if($target == 'none') then

  exit

else if($target == 'orig') then

  goto return_building_orig

else

  goto return_building_other

endif

################################################################################
      start_cleaning:
################################################################################

set curdir = `$pwd_`

$ls_ ${curdir}/$srcdir >& /dev/null

if($status == 0) then

  $rm_ -rf ${curdir}/$srcdir

  echo "rm -rf ${curdir}/${srcdir}"

else

  echo "information: directory, $srcdir not found. nothing to remove."

endif

$ls_ ${curdir}/BuiltVersion=* >& /dev/null

if($status == 0) then

  set targetdir = `$ls_ ${curdir}/BuiltVersion=*`

  set targetdir = `echo $targetdir | $cut_ -f 1 -d ' '`

  set targetdir = `$basename_ $targetdir | $cut_ -c 14-`

  if($targetdir == 'ORIG') then

    set targetdir = 0rig

  endif

  $ls_ ${curdir}/${extdir}/${targetdir}/* >& /dev/null

  if($status == 0) then

    foreach extcnts (${curdir}/${extdir}/${targetdir}/*)

      $ls_ ${curdir}/`$basename_ ${extcnts}` >& /dev/null

      if($status == 0) then

        $rm_ -rf ${curdir}/`$basename_ ${extcnts}`

        echo "rm -rf ${curdir}/`$basename_ ${extcnts}`"

      else

        echo 'information: directory or file,'\
             "`$basename_ ${extcnts}` not found. nothing to remove."

      endif

    end

  else

    echo 'information: can not search any directories and files'\
         "in directory, ${extdir}/${targetdir}. nothing to remove."

  endif

  $rm_ -f ${curdir}/BuiltVersion=*

  echo "rm -f ${curdir}/BuiltVersion=*"

else

  echo 'information:'\
       'can not search any directories and files to remove. do nothing.'

endif

if($target == 'none') then

  exit

else if($target == 'clean') then

  goto return_cleaning_clean

else if($target == 'orig') then

  goto return_cleaning_orig

else

  goto return_cleaning_other

endif
