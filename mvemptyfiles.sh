#! /bin/csh -f
foreach x (*.f90)

if (! -z $x ) then
echo $x
else
echo " file is empty \c"
mv $x ./00empty_file/
endif

end
