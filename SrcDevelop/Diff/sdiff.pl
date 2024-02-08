#!/usr/bin/perl

$BASEDIR1 = 'RAD1.2.3';
$BASEDIR2 = 'RAD1.2.1';

#===================================================================

open(FILE, "list.txt") or die "$!" ;
@file = <FILE> ;
close(FILE);

foreach $line (@file) {
  
  chomp $line;
  
  open(FILE2, ">./diff__$line.txt") or die "$!" ;
  
  print FILE2 "\n ***** $BASEDIR1 <===> $BASEDIR2 *****\n";
  print FILE2 "\n  [ target program : $line ]\n";

  $samepg = -1;

  if( -e "./Src_$BASEDIR1/$line" ){
  
     if( -e "./Src_$BASEDIR2/$line" ){
  
        system("diff ./Src_$BASEDIR1/$line ./Src_$BASEDIR2/$line > tmp");
  
        if( -z './tmp' ){
           print FILE2 "\n  diff size=0 => same program\n\n";
           $samepg = 1;
        }
        else{
           print FILE2 "\n  diff size>0 => different program\n\n";
##         print FILE2 "  sdiff -b -B -E -W -w 200 $BASEDIR1/$line $BASEDIR2/$line\n";
           print FILE2 "  sdiff -b -B -W -w 200 $BASEDIR1/$line $BASEDIR2/$line\n";
##         system("sdiff -b -B -E -W -w 200 Src_$BASEDIR1/$line Src_$BASEDIR2/$line >> ./diff__$line.txt");
           system("sdiff -b -B -W -w 200 Src_$BASEDIR1/$line Src_$BASEDIR2/$line >> ./diff__$line.txt");
        }
        system("rm tmp");
  
     }
     else{
        print FILE2 "\n  $BASEDIR2/$line is non-existent.\n";
     }

  }
  elsif( -e "./Src_$BASEDIR2/$line" ){
  
     print FILE2 "\n  $BASEDIR1/$line is non-existent.\n";
  
  }
  else{
  
     print FILE2 "\n  both of $line are non-existent.\n";
  
  }
  
  print FILE2 "\n";
  close(FILE2);

  if( $samepg==1 ){
    system("rm ./diff__$line.txt");
  }
  
}

