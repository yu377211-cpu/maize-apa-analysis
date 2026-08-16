#!/usr/bin/perl -w

require ("/share/pub/xingsl/shilai/project/APA/PAC_3seq_map/PAC_identification_scripts/funclib.pl");
use strict;
use Getopt::Long;


##########################################################################################
#  parameters & usage
##########################################################################################
my %opts=();
GetOptions(\%opts,"h","in=s","poly=s","ml=i","mp=i","mg=i","mm=i","mtail:i","mper:s","mr:i","odir:s","suf:s","reg:i","deep:s","oraw:s","debug:s","bar:s","review:s");
my $USAGE=<< "_USAGE_";
Require: funclib.pl
Usage: 
  0) 测试不同参数 -debug T 开启
  -debug T -oraw T

  1) 对于PATseq推荐的参数 
  MAP_findTailAT.pl -in E:/sys/code/testdata/arab1.fastq -poly T -ml 25 -mp 6 -mg 5 -mm 2 -mr 2 -mper 0.75 -mtail 8 -deep T -reg 1 -odir "f:/" -suf "" 
  MAP_findTailAT.pl -in E:/sys/code/testdata/arab1.fastq -poly T -ml 25 -mp 8 -mg 8 -mm 2 -odir "f:/" -suf "" -reg 1
  MAP_findTailAT.pl -in E:/sys/code/testdata/arab1.fastq -poly "A|T" -ml 25 -mp 8 -mg 8 -mm 2 -odir "f:/" -suf "AT" -reg 1

  不输出原始序列
  MAP_findTailAT.pl -in E:/sys/code/testdata/arab1.fastq -poly "A&T" -ml 25 -mp 8 -mg 8 -mm 2 -odir "f:/" -suf "A&T" -reg 1 -oraw F

  debug 输出各种情况
  MAP_findTailAT.pl -in E:/sys/code/testdata/arab1.fastq -poly "A&T" -ml 25 -mp 8 -mg 8 -mm 2 -odir "f:/" -suf "A&T" -reg 1 -oraw F -debug T

-h=help
-in=input fa or fq file
-poly=A/T/A&T/A|T
  A|T 若同时有A和T，当tail_lenA>=lenT时，保留A的结果；一条序列只能含有A或含有T，不能同时有
  A&T 同时找A和T, 二者不相干，一条序列可能同时有A和同时有T
-ml=min length after trim
-mg=margin from the start (poly=T) or to the end (poly=A) (=5)
-mm=mismatch between TxxTTT (=2)
-mr=minT in reg (=3)
-mp=min length of succesive poly (=8)
-mtail=min length of trimmed tail (=8)
-mper=min percent of A/T in trimmed tail (=0.75)
-deep=T(default)/F 若T则在reg不匹配时，再深度查找，比如这种类型 TTTTTCTTTTCTCTTTTTTTT
-odir=output path (default is the same as input)
-suf=suffix, default: xx.suf.T/A.fq
-reg=1(loose, default)/2
  1='^.{0,mg}?(T{mr,}[^T]{0,mm})*T{mp,}([^T]{0,mm}T{mr,})*'
  2='^.{0,mg}?(T{mp,}([^T]{0,mm}T{mr,})*)'
-oraw=T(default)/F; output raw reads for trimmed A/T file (.A.raw.fq)
-bar=Last bast position of three barcode 
-review=T(default)/F若T则在reg/deep不匹配时，再深度查找，比如AACCCCTTTTTTTTTTTTTTTTT#seven
_USAGE_


#############################################################################
#  invoke the program                               
#############################################################################
die $USAGE if $opts{'h'}||!$opts{'in'}||!$opts{'poly'}||!$opts{'ml'}||!$opts{'mp'}||!$opts{'mg'}||!$opts{'mm'};

my $poly=$opts{'poly'};
my $in=$opts{'in'};
my $ml=$opts{'ml'};
my $mp=$opts{'mp'};
my $mg=$opts{'mg'};
my $mm=$opts{'mm'};
my $mr=$opts{'mr'};
my $mtail=$opts{'mtail'};
my $mper=$opts{'mper'};
my $deep=$opts{'deep'};
my $debug=$opts{'debug'};
my $suf=$opts{'suf'};
my $odir=$opts{'odir'};
my $reg=$opts{'reg'};
my $oraw=$opts{'oraw'};
my $bar=$opts{'bar'};
my $review=$opts{'review'};
$reg=1 if !$reg;
die "reg=1/2" if $reg!=1 and $reg!=2;


if (!$mtail) {
  $mtail=8;
}
if (!$mper) {
  $mper=0.75;
}
if (!$mr) {
  $mr=3;
}
if (!$deep | $deep eq 'T') {
  $deep=1;
} else {
  $deep=0;
}
if (!$review | $review eq 'T') {
  $review=1;
} else {
  $review=0;
}
if (!$debug | $debug ne 'T') {
  $debug=0;
} else {
  $debug=1;
}
if(!$bar){
 $bar=0;	
}
$oraw=1 if !$oraw or $oraw ne 'F';
$oraw=0 if $oraw eq 'F';

print "poly=$poly\ninfile=$in\noutput dir=$odir\nmin length after trim(ml)=$ml\nmin continue tail length(mp)=$mp\nmax margin to the end(mg)=$mg\nmax mismatches between xTTT(mm)=$mm\n";
print "min reg T(mr)=$mr\nmin tail length(mtail)=$mtail\nmin T/A percent(mper)=$mper\n";
print "deep=$deep\n";
print "DEBUG MODE, output to .test.file too \n" if $debug;

die "error poly(=A T A|T A&T)" if $poly ne 'A' and $poly ne 'T' and $poly ne 'A|T' and $poly ne 'A&T';
my ($findT,$findA)=(0,0);
$findA=1 if ($poly=~ m/A/);
$findT=1 if ($poly=~ m/T/);
#print "$findT\n";

#if ($deep & $findA) {
#  die "TODO: deep和findA的程序还没实现...";
#}

die "only support fastq file" if seqFormat($in) ne 'fq';
open my $fh1, "<$in" or die "cannot read input file ($in)\n";
#open IN, "<$in" or die "cannot read input file ($in)\n";
$suf=".$suf" if $suf;

my ($cntNotailT,$cntShortT,$cntTotal,$cntFinalT,$cntNotailA,$cntShortA,$cntFinalA,$cntMissA,$cntMissT,$cntBadT,$cntBadA)=(0,0,0,0,0,0,0,0,0,0,0);

my $inpath=getDir($in,0); # x/y/
my $infname=getFileName($in,1); #文件前缀
if (!$odir) {
  $odir=$inpath;
}

if ($findA) {
  open OA,">$odir$infname$suf.A.fq";
  if($debug){
  	open TESTA,">$odir$infname$suf.testA.fq";
  	}
  if($oraw){
  	open RAWA,">$odir$infname$suf.A.raw.fq" ;
  	}
}


if ($findT) {
  open OT,">$odir$infname$suf.T.fq";#有效的结果，剔除了polyT尾巴的#seven
  if($debug){
  	open TESTT,">$odir$infname$suf.T.DEBUG.fq";
  	open BEDT,">$odir$infname$suf.badT.fq";
  	open NOT,">$odir$infname$suf.noTail.fq";
  	open SHORTT,">$odir$infname$suf.shortT.fq";
  	}
  if($oraw){
  	open RAWT,">$odir$infname$suf.T.raw.fq";
  	}
}

#print "###############$findT\t$debug\t$oraw\n";

my $regT="^.{0,$mg}?(T{$mr,}[^T]{0,$mm})*T{$mp,}([^T]{0,$mm}T{$mr,})*";
if ($reg==2) {
  $regT="^.{0,$mg}?T{$mp,}([^T]{0,$mm}T{$mr,})*";
}
print "regT=$regT\n";
my $regA="(A{$mr,}[^A]{0,$mm})*A{$mp,}([^A]{0,$mm}A{$mr,})*.{0,$mg}?\$";
if ($reg==2) {
  $regA="A{$mp,}([^A]{0,$mm}A{$mr,})*.{0,$mg}?\$";
}
print "regA=$regA\n";

if ($deep) {
  print "deep reg=^.{0,$mg}?(T{1,}[^T]{0,2})*T{8,}\n";
}
if($review){
   print "review_reg=T{$mp,}([^T]{0,$mm}T{$mr,}){0,}\n";	
}

my ($md1,$md2,$md3)=(0,0,0);


my @fqkeys = qw/id seq pl qual/;
while(!eof $fh1) {
	my (%e1);
	@e1{@fqkeys} = readfq($fh1);#hash同时存放reads的id、seq、pl、qual
	$cntTotal++;
	my ($seqT,$seqA,$qualT,$qualA)=('','','','');
	my ($haveT,$haveA,$shortT,$notailT,$shortA,$notailA,$badTailT,$badTailA)=(0,0,0,0,0,0,0,0);
	if ($findT) {
	    $seqT=$e1{'seq'};
 		if ($seqT=~s/$regT//) {#把匹配的模式替换了
		  my $length=length($seqT);
		  #my $Ts=substr($e1{'seq'},0,length($e1{'seq'})-$length); #截取下来的T片段#WXH
		  my $Ts=substr($e1{'seq'},$bar,length($e1{'seq'})-$length-$bar);#应该剔除掉barcode，adapter，Seven 
		  #print "$Ts\n";#Seven test
		  #my $Tcnt = $Ts =~ tr/T/T/; #计算T数
		  my $Tcnt = $Ts=~ tr/T/T/;#计算T数，Seven
		  ##print "Yes: $Ts.$Tcnt\n";
		  if ((length($e1{'seq'})-$length-$bar)<$mtail or $Tcnt/length($Ts)<$mper) { #tail过短 or T%<mper
			$badTailT=1;
			print TESTT $Ts."***BAD***".$seqT."\n"  if $debug;
		  } elsif ($length<$ml) { #序列过短
			$shortT=1;
			print TESTT $Ts."***SHR***".$seqT."\n"  if $debug;
		  } else {
		  	++$md1;
			$haveT=1;
		    $qualT=substr($e1{'qual'},length($e1{'qual'})-$length,$length);#满足条件呢seven
		  }
		} 

       if (!$haveT & $deep) { #深度查找,判断固定的正则式，找到再截取到开头，再判断mper.. (TTTTTCTTTTCTCTTTTTTTT)
          $seqT=$e1{'seq'};
		      my $deepreg="^.{0,$mg}?(T{1,}[^T]{0,2})*T{$mp,}";
		      if ($seqT=~s/$deepreg//){
			   		 my $length=length($seqT);
			 			 #my $Ts=substr($e1{'seq'},0,length($e1{'seq'})-$length); #截取下来的T片段 WXH
			 			 my $Ts=substr($e1{'seq'},$bar,length($e1{'seq'})-$length-$bar);#应该剔除掉barcode，adapter，Seven
			 			   #print "$Ts\n";
			  	  my $Tcnt = $Ts =~ tr/T/T/; #计算T数
		      		##print "Deep: $Ts.$Tcnt\n";
			      if ((length($e1{'seq'})-$length-$bar)>=$mtail and $Tcnt/length($Ts)>=$mper and $length>=$ml) { 
				      ++$md2;
				      $haveT=1;
				      $qualT=substr($e1{'qual'},length($e1{'qual'})-$length,$length);
				      print TESTT $Ts."***DEP***".$seqT."\n"  if $debug;
				      $badTailT=0;
				      $shortT=0;
			       }	  
		     }
	     }
	     
	     if(!$haveT & $review){#在度深度查找，只要匹配到TTTTTTTT的,判断mper、ml
	     	   $seqT=$e1{'seq'};
	     	   my $deepmod="T{$mp,}([^T]{0,$mm}T{$mr,}){0,}";
	     	   if($seqT=~ m/$deepmod/){
	     	   	
	     	     #print "1$seqT\n";
	     	     my $pos_start=$-[0];#返回匹配的首位置，从0开始计数
	     	     my $match_length =length$&;#返回第一个匹配的长度 TTTTTT，或TTTTTTXXTTTXXTTTT
             my $pos_end = $pos_start+$match_length;#返回首次匹配的末index，从0开始计数
  
	     	   	 $seqT =substr($seqT,$pos_end,length($seqT)-$pos_end);
	     	   	 #print "2$seqT\n";
	     	   	 my $length=length($seqT);
	     	   	 
	     	   	 my $Ts=substr($e1{'seq'},$bar,$pos_end-$bar);#polyT的序列
	     	   	 #print "3$Ts\n";
	     	   	 my $Tcnt = $Ts =~ tr/T/T/; #计算T数
	     	   	 #print "3$Tcnt\n";
	     	   	 if(length($Ts)>=$mtail and $Tcnt/length($Ts)>=$mper and $length>=$ml){
	     	   	 	  #print "####yes\n";
	     	   	 	  ++$md3;
	     	   	 	  $haveT=1;
	     	   	 	  $qualT=substr($e1{'qual'},$pos_end,length($e1{'qual'})-$pos_end);
                                  if($debug){
                                    print TESTT $Ts."***DEPSEVEN***".$seqT."\n"  if $review;
                                  }
	     	   	 	  
				        $badTailT=0;
				        $shortT=0;
				        #print "$seqT\n";
	     	   	 	}
	     	   	}
	     	
	     	}
	     

	   if (!$haveT) {
		  $notailT=1;	
		  if($debug){
		   # print TESTT "****NOT***$e1{'seq'}\n";#seven
		  }	  		  
	   } elsif ($debug) {
		  # print TESTT "****YES***$e1{'seq'}\n";
	   }
	}#end find poly(T)

	if ($findA) {
	    $seqA=$e1{'seq'};
		if ($seqA=~s/$regA//) {
	      my $length=length($seqA);
		  print TESTA "$seqA*******".substr($e1{'seq'},$length,length($e1{'seq'})-$length)."\n" if $debug;
		  my $As=substr($e1{'seq'},$length,length($e1{'seq'})-$length); #截取下来的A片段
		  my $Acnt = $As =~ tr/A/A/; #计算A数
		  if (length($e1{'seq'})-$length<$mtail or $Acnt/length($As)<$mper) { #tail过短
			$badTailA=1;
		  } elsif ($length<$ml) {
			 $shortA=1;
		  } else {
			$haveA=1;
			$qualA=substr($e1{'qual'},0,$length);
		  }
		} else {
		  $notailA=1;
		  print TESTA "*******$e1{'seq'}\n" if $debug;
		}
	} #~findA


      #输出reads结果
      if ($poly eq 'A|T') {
		  		if ($haveA and $haveT) { #若同时有A和T，则lenA>=lenT时，保留A的结果；
			 				if (length($seqT)<length($seqA)) {
									print OT "$e1{'id'}\n$seqT\n$e1{'pl'}\n$qualT\n";
									$cntFinalT++;
									$cntMissA++;
									if ($oraw) {
				  				print RAWT "$e1{'id'}\n$e1{'seq'}\n$e1{'pl'}\n$e1{'qual'}\n";
									}
			 				} else {
								print OA "$e1{'id'}\n$seqA\n$e1{'pl'}\n$qualA\n";
								if ($oraw) {
				  				print RAWA "$e1{'id'}\n$e1{'seq'}\n$e1{'pl'}\n$e1{'qual'}\n";
								}
								$cntFinalA++;
								$cntMissT++;
						 }
		  		}elsif ($shortA or $shortT) { #有A或T太短，则都不输出
			  				$cntShortA++;
			 					$cntShortT++;
		 			}elsif ($badTailT or $badTailA) { #有A或T太短，则都不输出
			 					 $cntBadA++;
			  				 $cntBadT++;
		      }elsif ($haveA and !$haveT) { #有一方没有，则只输出相应方
			  				 print OA "$e1{'id'}\n$seqA\n$e1{'pl'}\n$qualA\n";
								if ($oraw) {
				  					print RAWA "$e1{'id'}\n$e1{'seq'}\n$e1{'pl'}\n$e1{'qual'}\n";
								}
			   				$cntFinalA++;
			          $cntNotailT++;
		     }elsif ($haveT and !$haveA) { #有一方没有，则只输出相应方
			   				print OT "$e1{'id'}\n$seqT\n$e1{'pl'}\n$qualT\n";
				        if ($oraw) {
				  				print RAWT "$e1{'id'}\n$e1{'seq'}\n$e1{'pl'}\n$e1{'qual'}\n";
								}
			         $cntFinalT++;
			         $cntNotailA++;
		     } elsif (!$haveT and !$haveA) {
			  		 $cntNotailA++;
			   			$cntNotailT++;
		     }
      } else {
         if ($haveA) {
						print OA "$e1{'id'}\n$seqA\n$e1{'pl'}\n$qualA\n";
			     if ($oraw) {
							print RAWA "$e1{'id'}\n$e1{'seq'}\n$e1{'pl'}\n$e1{'qual'}\n";
						}
					 $cntFinalA++;
         } elsif ($shortA) {
			      $cntShortA++;
         }elsif ($badTailA) {
			      $cntBadA++;
         }elsif ($notailA) {
			      $cntNotailA++;
         }#for poly A tail #seven
         
         if ($haveT) {
						print OT "$e1{'id'}\n$seqT\n$e1{'pl'}\n$qualT\n";
			     if ($oraw) {
				    #print RAWT "$e1{'id'}\n$e1{'seq'}\n$e1{'pl'}\n$e1{'qual'}\n";
				    print RAWT "$e1{'seq'}\n";
			     }
			    $cntFinalT++;
         } elsif ($shortT) {
			      $cntShortT++;
			      #print SHORTT "$e1{'id'}\n$e1{'seq'}\n$e1{'pl'}\n$e1{'qual'}\n" if $debug; #seven
			      print SHORTT "$e1{'seq'}\n" if $debug; #seven
         }elsif ($badTailT) {
		       	$cntBadT++;
             #print BEDT "$e1{'id'}\n$e1{'seq'}\n$e1{'pl'}\n$e1{'qual'}\n" if $debug; #seven
             print BEDT "$e1{'seq'}\n" if $debug; #seven
         }	elsif ($notailT) {
           # print NOT "$e1{'id'}\n$e1{'seq'}\n$e1{'pl'}\n$e1{'qual'}\n" if $debug; #seven
           print NOT "$e1{'seq'}\n" if $debug; #seven
		       	$cntNotailT++;
         }  
	    }

}#end while 
close($fh1);

if($findA){
   close(OA);
   if($debug){
   	close(TESTA);
   	}	
   	if($oraw){
   	close(RAWA);	
   	}
}


if($findT){
   close(OT);
   if($debug){
   	close(TESTT);
   	close(NOT);
   	close(BEDT);
   	close(SHORTT);
   	}	
   	if($oraw){
   	close(RAWT);	
   	}
}



print "\ntotal\t$cntTotal\n";
if ($findA) {
  print "[polyA]\nfinal\t$cntFinalA\nnotail\t$cntNotailA\ntooshort\t$cntShortA\nbadTail\t$cntBadA\nmissby(A|T)\t$cntMissA\n";
  if ($cntTotal!=$cntFinalA+$cntNotailA+$cntShortA+$cntMissA+$cntBadA) {
	print "cntTotal!=cntFinalA+cntNotailA+cntShortA+cntMissA+cntBadA\n";
  }
}
if ($findT) {
  print "[polyT]\nfinal\t$cntFinalT\nnotail\t$cntNotailT\ntooshort\t$cntShortT\nbadTail\t$cntBadT\nmissby(A|T)\t$cntMissT\n";
  print "region\t$md1\ndeep\t$md2\nreview\t$md3\n";
  if ($cntTotal!=$cntFinalT+$cntNotailT+$cntShortT+$cntMissT+$cntBadT) {
	print "cntTotal!=cntFinalT+cntNotailT+cntShortT+cntMissT+cntBadT\n";
  }
}


###############################################
sub readfq {
        my $fh = pop @_;
        my @entry = ();
        for (qw/id seq pl qual/) {
                my $line = readline($fh); chomp $line; 
                push @entry, $line; warn "line empty" if !$line;
        }
        return @entry;
}

