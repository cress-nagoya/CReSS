################################################################################
#                                                                              #
#     Read me file, readme_first.txt                                           #
#                                                                              #
#     This file explain how to compile and run the CReSS and its preprocessor  #
#     programs.                                                                #
#                                                                              #
#     Author      : Sakakibara Atsushi                                         #
#     Date        : 2000/05/08                                                 #
#     Modification: 2000/12/18, 2001/01/15, 2001/03/13, 2001/10/18,            #
#                   2002/06/17, 2002/07/03, 2002/09/09, 2002/10/31,            #
#                   2003/05/19, 2003/07/15, 2004/03/05, 2004/05/31,            #
#                   2004/09/10, 2004/12/17, 2005/01/14, 2006/04/03,            #
#                   2006/12/04, 2007/05/07, 2007/09/28, 2008/10/10,            #
#                   2008/12/11, 2009/01/30, 2009/02/04, 2009/03/31,            #
#                   2009/06/17, 2010/02/01, 2010/05/17, 2010/09/22,            #
#                   2011/01/19, 2011/03/29, 2011/06/01, 2011/07/15,            #
#                   2011/08/09, 2011/08/18, 2011/09/22, 2011/11/10,            #
#                   2011/12/17, 2012/01/25, 2012/06/19, 2013/01/27,            #
#                   2013/01/28, 2013/02/05, 2013/02/13, 2013/03/27,            #
#                   2013/10/08                                                 #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
#                                                                              #
# 1. Directory structure                                                       #
#                                                                              #
#------------------------------------------------------------------------------#

    The CReSS3.4.3m_20131008 directory contains the following directories and
files.

      Doc             : directory; Contains the documentations.
      Extra           : directory; Contains the user extra data directory.
      Form            : directory; Contains the formats of simple simulations.
      SrcDevelop      : directory; Contains the modified source codes.
      SrcOrig         : directory; Contains the original source codes.
      compile.conf    : file     ; Set the compile options.
      compile.csh     : file     ; Generates the executable files.
      configure.csh   : file     ; Construct Src directory.

#------------------------------------------------------------------------------#
#                                                                              #
# 2. How to compile the CReSS and its preprocessor and postprocessor programs  #
#                                                                              #
#------------------------------------------------------------------------------#

    First, You should construct Src directory. you should input the following
command in the CReSS3.4.3m_20131008 directory.

      % configure.csh

    Next, you should modify the file, compile.conf to adjust your machine's
compiler. After modifying, you should input the following commands in the
CReSS3.4.3m_20131008 directory.

      % compile.csh solver [compile.conf]
      % compile.csh gridata [compile.conf]
      % compile.csh asldata [compile.conf]
      % compile.csh radata [compile.conf]
      % compile.csh terrain [compile.conf]
      % compile.csh surface [compile.conf]
      % compile.csh check [compile.conf]
      % compile.csh unite [compile.conf]
      % compile.csh rstruct [compile.conf]

    The first arguments are solver, gridata, asldata, radata, terrain, surface,
check, unite or rstruct. After compiling, generate the objective and the
executable files in the Src directory and the simbolic link of the executable
files in the CReSS3.4.3m_20131008 directory.

    If you want to clean up the executable and the objective files, input the
following command in the CReSS3.4.3m_20131008 directory.

      % compile.csh clean

    And if you want to clean up the Src and/or extra directories, input the
following command in the CReSS3.4.3m_20131008 directory.

      % configure.csh clean

#------------------------------------------------------------------------------#
#                                                                              #
# 3. How to run the CReSS and its preprocessor and postprocessor programs      #
#                                                                              #
#------------------------------------------------------------------------------#

    The user configuration file is read in running. For example you will run the
solver, you should input the following command.

      % solver.exe < user_configuration_file > standard_output &

    You can move the executable files, solver.exe, gridata.exe, asldata.exe,
radata.exe, terrain.exe, surface.exe, check.exe, unite.exe and rstruct.exe and
run them in any directories. The required files are the executable files,
the user configuration file and the input files for the executable programs.
The input files for the executable programs are:

      runname.sounding.txt

      runname.gpv.check.txt
      runname.gpvXXXXXXXX.grpYYYY-subZZZZ.bin

      runname.asl.check.txt
      runname.aslXXXXXXXX.grpYYYY-subZZZZ.bin

      runname.rdr.check.txt
      runname.rdrXXXXXXXX.grpYYYY-subZZZZ.bin

      runname.terrain.check.txt
      runname.terrain.grpYYYY-subZZZZ.bin
      runname.terrain.damp.grpYYYY-subZZZZ.bin

      runname.surface.check.txt
      runname.surface.grpYYYY-subZZZZ.bin

      runname.sst.check.txt
      runname.sstXXXXXXXX.grpYYYY-subZZZZ.bin

      runname.resXXXXXXXX.check.txt
      runname.resXXXXXXXX.grpYYYY-subZZZZ.bin

    And these files have to be in your specified directory.

#------------------------------------------------------------------------------#
#                                                                              #
# 4. File format                                                               #
#                                                                              #
#------------------------------------------------------------------------------#

    Sounding file format: You can use following formats. The references are in
the Form directory.

      1st column: height [m] or pressure [Pa]
      2nd column: temperature [K] or potential temperature [K]
      3rd column: x components of velocity [m/s]
      4th column: y components of velocity [m/s]
      5th column: water vapor mixing ratio [kg/kg] or relative humidity [%]

    Dumped file format : The dumped files are output by text or direct accsess
binary format. You can dump optional variables. The velocity variables u, v and
w and zph are interpolated to the scalar points. And the output index ranges for
each processor element are following.

      text format;

        do k=2,nk-2

          write(ionum,*,err=errnum)                                       &
       &       ((variable(i,j,k),i=2,ni-2),j=2,nj-2)

        end do

      binary format;

        do k=2,nk-2
          recdmp=recdmp+1

          write(ionum,rec=recdmp,err=errnum)                              &
       &       ((variable(i,j,k),i=2,ni-2),j=2,nj-2)

        end do

      To display these files, you have to unite to one file for each damped
forecast time.

    External file format: If you can read in the external data to the model,
input the following command and modify the module rdgpv, rdaero, rdradar,
rdheight, rdland, rdsst and/or rdice in the SrcDevelop/UserModAdd directory.

      % grep '#####\ ' *.f90
      rdaero.f90:! ##### You will have to modify the following lines. #####
      rdaero.f90:! ##### You will have to modify the following lines. #####
      rdgpv.f90:! ##### You will have to modify the following lines. #####
      rdgpv.f90:! ##### You will have to modify the following lines. #####
      rdheight.f90:! ##### You will have to modify the following lines. #####
      rdheight.f90:! ##### You will have to modify the following lines. #####
      rdice.f90:! ##### You will have to modify the following lines. #####
      rdice.f90:! ##### You will have to modify the following lines. #####
      rdland.f90:! ##### You will have to modify the following lines. #####
      rdland.f90:! ##### You will have to modify the following lines. #####
      rdradar.f90:! ##### You will have to modify the following lines. #####
      rdradar.f90:! ##### You will have to modify the following lines. #####
      rdsst.f90:! ##### You will have to modify the following lines. #####
      rdsst.f90:! ##### You will have to modify the following lines. #####

#------------------------------------------------------------------------------#
#                                                                              #
# 5. Sample simulation                                                         #
#                                                                              #
#------------------------------------------------------------------------------#

    Now we explain the cats eye simulation whose formats are in the Form
directory.

    First, you should construct the Src directory in the CReSS3.4.3m_20131008
directory.

      % pwd
      /....../CReSS3.4.3m_20131008

      % ls
      Doc            Form           SrcOrig        compile.csh
      Extra          SrcDevelop     compile.conf   configure.csh

      % configure.csh

      % ls
      Doc            Form           SrcDevelop     compile.conf   configure.csh
      Extra          Src            SrcOrig        compile.csh

    Second, you should copy the required files in the CReSS3.4.3m_20131008
directory.

      % cp Form/form_cats-eye.sounding.txt cats-eye.sounding.txt
      % cp Form/form_cats-eye.user.conf user.conf

      % ls
      Doc                    SrcDevelop             compile.csh
      Extra                  SrcOrig                configure.csh
      Form                   cats-eye.sounding.txt  user.conf
      Src                    compile.conf

    Third, you should compile the source codes and generate the executable
file, solver.exe. So you should input the following command.

      % compile.csh solver
                  :
                  :

    After compiling, generate the file, solver.exe in the CReSS3.4.3m_20131008
directory.

      % ls
      Doc                    SrcDevelop             compile.csh
      Extra                  SrcOrig                configure.csh
      Form                   cats-eye.sounding.txt  solver.exe
      Src                    compile.conf           user.conf

    Next, run that executable file, solver.exe.

      % solver.exe < user.conf > solver.log &

    When runnig is over, generate the following files in your specified
directory:

      cats-eye.geography.check.txt            : geography cheking files
      cats-eye.geography.grpYYYY-subZZZZ.bin  : dumped geography files
      cats-eye.dmp.check.txt                  : history cheking files
      cats-eye.dmpXXXXXXXX.grpYYYY-subZZZZ.bin: dumped history files
      cats-eye.resXXXXXXXX.check.txt          : restart checking file
      cats-eye.resXXXXXXXX.grpYYYY-subZZZZ.bin: restart files

    And the dumped geography and history files should be used for the graphic
applications.

####################################################
####                                            ####
####  More information is in the User's Guide.  ####
####                                            ####
####################################################
