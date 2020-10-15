#!/usr/bin/perl
#This script expands the windows and reports a region with signal abs(m.value)>=0.6 and p-value=0.001.
#Nima Rafati: 20130603:20130622:20130721:20130909
#my $interMapFile=$ARGV[0];
my $cntr=0;
my $outCut=2;
my $pValCut=0.001;
my $mValCut=0.6;
my $distanceCut=1000;
my $prvStart=-1;
my $openLocus=1;
my $newStart=1;
my $outlierCntr=0;
my $chr="";
my $openningLine="";
use Getopt::Long;
&GetOptions('pval_col=i' =>\$PVal_col, 'mval_col=i' =>\$MVal_col, 'f=s' =>\$interMapFile, 'pval_cut=f' => \$pValCut, 'mval_cut=f' => \$mValCut, 'distance_cut=i' =>\$distanceCut);
my $contact = "nimarafati\@gmail.com";
my $usage = "detect-regions-by-abs-M.Val-V5-with-dynamic-M.Val-Flag.pl 
-pval_col column of p-values
-mval_col column of m-values
-pval_cut p-value cut-off (default 0.001)
-mval_cut m-value cut off (default 0.6 which the absolute value will be used)
-distance_cut the distance between signals/windows/SNPs (default 1000 bp)
-f input file which should be in bed format
20130603
20130622
20130721
20130909
20140528
20140618 adding distance & openning site to flags
20140623 correct to report single windows 
20140809 do not print empty lines after loopping back to extend the locus
$constact\n";

if($PVal_col eq "" || $MVal_col eq "" || $interMapFile eq "" || $pValCut >=1 || $distanceCut == 0)
{
	print $usage;
	exit
}
#print "pval_col=$PVal_col\tmval_col=$MVal_col\tpValCut=$pValCut\tmValCut=$mValCut\n";<STDIN>;
my $startMValCut=$mValCut;
my $outOpenFile="$interMapFile-Openning.txt";
open(outFile,">$outOpenFile");
open($inFile,$interMapFile) || die "No input file\nUsage:perl *.pl intersected-file.txt\n";
while (<$inFile>)
{
	$cntr++;
	chomp($_);
	$_=~ s/\s+/\t/g;
	$chrArr[$cntr]=$_;
	my @lineArr=split("\t",$_);
	$chr=$lineArr[0];
	$start=$lineArr[1];
	$end=$lineArr[2];
	$pValue=$lineArr[$PVal_col-1];
	$mValue=$lineArr[$MVal_col-1];
	$nowSign=GiveSign($mValue);


	if($prvStart == 0 ||$prvStart == -1 || $prvChr eq "" || $chr ne $prvChr)
	{
		$distance=1;
	}
	elsif($chr eq $prvChr)
	{
		$distance=$start-$prvStart;
	}
#if($chr eq "scaffold699"){
#print "Got you\nopenLocus?$openLocus\t$outlierCntr\tmVal:$mValue\tpVal:$pValue\tstart:$start\tend:$end\tdistance:$distance\n$_\n";<STDIN>;
#print "size:",scalar(@outArr),"\n";foreach my $sd (@outArr){print $sd;<STDIN>;}$check=1}else{$check=0;}

#if($check==1){print "Chr=$chr\tStart=$start\tdsitance=$distance\tPval: $pValue\tMval:$mValue\t$distanceCut\n$_";<STDIN>;}
	##First window
	if($prvChr eq "" && abs($mValue)>=$mValCut && $pValue<=$pValCut && $openLocus==1)#Open the locus if the first window/line of the file is significant
	{
		$openningLine=$_;
		$openLocus=0;
		$forward=1;
		$openSign=$nowSign;
		$locusPos=$cntr;
		$openSign=$nowSign;
		$openMVal=$pValue;
		$newLine="$_";
		push(@outArr,$newLine);
#		$locCntrOpenning=$locCntr+1;
#		print outFile "$chr\t$start\t$end\tLocus-$locCntrOpenning\t$mValue\t$pValue 1st\n"; #Save the opening sitei
		print outFile "$chr\t$start\t$end\t$mValue\t$pValue\n"; #Save the opening site
	}
	if($chr ne "" && $chr eq $prvChr && $prvStart != -1 && $distance<=$distanceCut && $distance >0)#Avoid gaps on the same chromosomes
	{
		if (abs($mValue)>=$mValCut && $pValue<=$pValCut && $openLocus==1)  #Open the locus: filter for absoulte M.Value larger equal to $mValCut
		{
				@ouArr=();
				$openningLine=$_;
				$openLocus=0;
				$forward=1;
				$openSign=$nowSign;
				$locusPos=$cntr;
				$openSign=$nowSign;
				$openMVal=$pValue;
				$newLine="$_";
				push(@outArr,$newLine);
#				$locCntrOpenning=$locCntr+2;
#				print outFile "$chr\t$start\t$end\tLocus-$locCntrOpenning\t$mValue\t$pValue 2nd\n";#Save the opening site
				print outFile "$chr\t$start\t$end\t$mValue\t$pValue\n";#Save the opening site
#				$locCntr++;
		}
		elsif($openLocus==0 && $forward==1 && $outlierCntr <= $outCut && $start-$prvStart<=$distanceCut) #Go forward after opening the locus with second p-value cut off if the number of outliers is less than equal to 2
		{
#                if($check==1){print "HÄR==222222==$distance=>$newLine\n",$_;<STDIN>;foreach my $sd (@outArr){print $sd;<STDIN>;}}
#				print "Go Forward $mValue\t",$start-$prvStart,"\n@outArr\n";<STDIN>;
				$newLine="$_";
				if (abs($mValue) >= $mValCut/2 && $nowSign == $openSign && $outlierCntr <= $outCut && $pValue<=$pValCut*50) #Check the sign of the window after opening site
				{
					push(@outArr,$newLine);
					$outliercntr =0; 
					if(abs($mValue)>$mValCut)#Reset M-Value cutoff if reaches a larger M-Value thatn default M-Value cutoff
					{
						$mValCut=abs($mValue);
					}
#					if($check==1){print "cond1-$outlierCntr:\n$_";<STDIN>;}
				}
#				elsif(abs($mValue) <= $mValCut/2 && $outlierCntr <= $outCut) #count the number of outlier
				else
				{
					push(@outArr,$newLine);
					$outlierCntr++;
#					if($check==1){print "cond2-$outlierCntr:\n$_";<STDIN>;}
				}
		}
		if($outlierCntr==$outCut+1) #Go Backward since the number of outliers reached more than 2 after removing outliers from the end of array
		{
			$size=scalar(@outArr);
			$outlierCntr=0;
			$mValCut=$startMValCut;
			if ($size!=0)# remove the last three windows because they are outliers
			{
				@tstArray=@outArr;
				pop(@outArr);
				pop(@outArr);
				pop(@outArr);
#				$locCntr++;
				if(scalar(@outArr)==0)#if @outArr became empty after removing outliers (which in principle it shouldn't) we add opening site to array
				{
					push(@outArr,$openningLine);
				}		

			}
			for (my $i=$locusPos-1;$i>=$newStart;$i--)#NOW start looping back in @chrArr and update @outArr
			{
				$newLine=$chrArr[$i];
				@newLineArr=split("\t",$newLine);
				$newChr=$newLineArr[0];
				$newPValue=$newLineArr[$PVal_col-1];
				$newMValue=$newLineArr[$MVal_col-1];
				$newSign=GiveSign($newMValue);	
				if($prvChr eq $newChr){
				if (abs($newMValue) >= $mValCut/2 && $newSign == $openSign  && $outlierCntr<=$outCut && $newPValue<=$pValCut*50)#Check the sign of the flanking window
				{
					unshift(@outArr,$newLine);
					$outlierCntr=0;
				}
#				elsif(abs($newMValue) <= $mValCut/2 && $outlierCntr<=$outCut)#count the number of outlier
				else #Add element evenif it is deviating the criteria. it will be removed after reaching 3 outliers
				{
					unshift(@outArr,$newLine);
					$outlierCntr++;
				}}
				if($outlierCntr == $outCut+1) #remove the first three windows because they are outliers
				{
#					print "REMOVE======>$outArr[0]\n";
					shift(@outArr);
					shift(@outArr);
					shift(@outArr);
           	        if(scalar(@outArr)==0)#if @outArr became empty after removing outliers (which in principle it shouldn't) we add opening site to array
					{
						push(@outArr,$openningLine);
					}
					$outlierCntr=0;
					$newStart=$cntr;
					$openSign="";
					last;
				}
			}
			$size=scalar(@outArr);
			if($size>=1)#print the loci (@outArr)
			{
				$locCntr++;
				foreach (@outArr)
				{
					if($_ ne "")
					{
						print "$_\tLocus-$locCntr\n";#<STDIN>;
					}
				}
				$mValCut=$startMValCut;
			}
			@outArr=();
			$openLocus=1;
			$openningLine="";
			#Check the current window by filtering for absoulte M.Value larger equal to $mValCut in order to open a new locus
			#It's less likely since current line has been called as outlier
		  	if (abs($mValue)>=$mValCut && $pValue<=$pValCut && $openLocus==1)  
		    {
#				print "Detected first\t$_";<STDIN>;
	        	$openningLine=$_;
				$openLocus=0;
				$forward=1;
				$openSign=$nowSign;
				$locusPos=$cntr;
				$newLine="$_";
				push(@outArr,$newLine);
#				$locCntrOpenning=$locCntr+2;
#				print outFile "$chr\t$start\t$end\tLocus-$locCntrOpenning\t$mValue\t$pValue 3rd\n";#Save the opening site
				print outFile "$chr\t$start\t$end\t$mValue\t$pValue\n";#Save the opening site
				$locCntr++;
			}
		}
	}
	elsif($chr ne "" && $chr eq $prvChr && $start-$prvStart > $distanceCut )#If there is a gap in chromosome then print the previous identified locus
	{
#		print "GAP==>$_";<STDIN>;
		$size=scalar(@outArr);
		$outlierCntr=0;
		$mValCut=$startMValCut;
		if ($size!=0) #Go Backward since we reached a gap
		{
			$locCntr++;
			for(my $n=$size-1;$n>=$size-4;$n--)
			{
#				$outArr[$n]=~ s/\s+/\t/g;
				@thisArr=split("\t",$outArr[$n]);
				$thisMValue=$thisArr[$MVal_col-1];
				$thisSign=GiveSign($thisMValue);
#				if($check==1){print "+++++++LOOPBACK-opensign=$openSign/thissign=$thisSign:$outArr[$n]";<STDIN>;}
				if ($thisSign!=$openSign)
				{
					if($check==1){foreach $sd (@outArr){print "##",$sd;<STDIN>;}}
					pop(@outArr);
				}
				else
				{
					last;
				}
			}
#			print "$locusPos----$newStart\t$cntr";<STDIN>;
			for (my $i=$locusPos-1;$i>=$newStart;$i--)#Loop back in @outArr
			{
#				if($check==1){print "$i==$locusPos==$outlierCntr=$distance=>$chrArr[$i]==",scalar(@outArr),"\n$chrArr[$i]";<STDIN>;}
				$newLine=$chrArr[$i];
				@newLineArr=split("\t",$newLine);
				$newChr=$newLineArr[0];
				$newPValue=$newLineArr[$PVal_col-1];
				$newMValue=$newLineArr[$MVal_col-1];
				$newSign=GiveSign($newMValue);
				if($newChr eq $prvChr)
				{	
					if (abs($newMValue)>=$mValCut/2 && $newSign == $openSign && $outlierCntr<=$outCut && $newPValue<=$pValCut*50)
					{
						if($newChr eq $prvChr)
						{
							unshift(@outArr,$newLine);#if($check==1){print "A $outlierCntr $outCut\n",$newLine;<STDIN>;}
							$outlierCntr=0;
						}
					}
#					elsif(abs($newMValue) <= $mValCut/2 && $outlierCntr<=$outCut)
					else
					{
						if($newChr eq $prvChr)
						{
							unshift(@outArr,$newLine);#if($check==1){print "B $outlierCntr $outCut\n",$newLine;<STDIN>;}
							$outlierCntr++;
						}
					}
					if($newChr eq $prvChr && $outlierCntr==$outCut+1) #remove the first three windows because they are outliers
					{
#						if($check==1){ print "REMOVE======>$outArr[0]\n";}
						shift(@outArr);
						shift(@outArr);
						shift(@outArr);
						if(scalar(@outArr)==0)
						{
							push(@outArr,$openningLine);
						}
						$outlierCntr=0;
						$newStart=$cntr;
						$openSign="";
					}
				}
			}
			$size=scalar(@outArr);
			if($size>=1)
			{
				$locCntr++;
				foreach (@outArr)
				{
					print "$_\tLocus-$locCntr\n";#if($check==1){<STDIN>;}
				}
				@outArr=();
			}
			$openLocus=1;
			$openningLine="";
		}
	   	if (abs($mValue)>=$mValCut && $pValue<=$pValCut && $openLocus==1)  #Call new loci at current window in new chromosome
		{
		    $openLocus=0;
			$forward=1;
			$openSign=$nowSign;
			$locusPos=$cntr;
			$newLine="$_";
			$openningLine=$_;
			push(@outArr,$newLine);
#			$locCntrOpenning=$locCntr;
#			print outFile "$chr\t$start\t$end\tLocus-$locCntrOpenning\t$mValue\t$pValue 4th\n";#Save the opening site
			print outFile "$chr\t$start\t$end\t$mValue\t$pValue\n";#Save the opening site
		}
		$newStart=$cntr;
	}
	elsif ($prvChr eq "" ||  $chr ne $prvChr) ##On the same chromosome open new locus after printing the previous one
	{
                        $size=scalar(@outArr);
#print "CNTR====>$_\n$size @ line $cntr";<STDIN>;
                        if($size>=1)
                        {
				$locCntr++;
                                foreach (@outArr)
                                {
                                        print "$_\tLocus-$locCntr\n";#<STDIN>;
					$printedLine=$_;
                                }
                        }
                        $openLocus=1;
                        @outArr=();

        	if (abs($mValue)>=$mValCut && $pValue<=$pValCut && $openLocus==1)  #filter for absoulte M.Value larger equal to $mValCut
	        {
	#	print "First=",scalar(@outArr),"\n$_";<STDIN>;
#if($check==1){print "First=",scalar(@outArr),"\n$_";<STDIN>;}
                        $openLocus=0;
                        $forward=1;
                        $openSign=$nowSign;
                        $locusPos=$cntr;
                        $newLine="$_";
                        push(@outArr,$newLine);
#			$locCntrOpenning=$locCntr;
#			print outFile "$chr\t$start\t$end\tLocus-$locCntrOpenning\t$mValue\t$pValue 5th\n";#Save the opening site
			print outFile "$chr\t$start\t$end\t$mValue\t$pValue\n";#Save the opening site
        	}
	        elsif($openLocus==0 && $forward==1 && $outlierCntr<=$outCut) #Go forward after opening the locus with second p-value cut off if the number of outliers is less than equal to 2
        	{       
                    $newLine="$_";
	#		print "Second=",scalar(@outArr),"\n$_";<STDIN>;
                   	if (abs($mValue)>=$mValCut && $newSign == $openSign && $outlierCntr<=$outCut) #Check the sign of the flanking window
			{
				push(@outArr, $newLine);
		       	        $outlierCntr=0;
            		}
	        }
	}
	##Reset the parameters:
	$prvSign=$nowSign;
	$prvStart=$start;
	$prvChr=$chr;
}
close $inFile;

##Check the last @outArr for printing
$size=scalar(@outArr);
$outlierCntr=0;
if ($size!=0)
{
	$locCntr++;
#	print @outArr,"\n";<STDIN>;
	for(my $n=$size-1;$n>=$size-4;$n--)
	{
#		$outArr[$n]=~ s/\s+/\t/g;
		@thisArr=split("\t",$outArr[$n]);
		$thisMValue=$thisArr[$MVal_col-1];
		$thisSign=GiveSign($thisMValue);
		if ($thisSign!=$openSign)
		{
			pop(@outArr);
		}
		else
		{
			last;
		}
	}
	for (my $i=$locusPos-1;$i>=$newStart;$i--)
	{
		$newLine=$chrArr[$i];
		@newLineArr=split("\t",$newLine);
		$newChr=$newLineArr[0];
#		print "looping back chr: $newChr\n$chr\t$prvChr";<STDIN>;
		$newPValue=$newLineArr[$PVal_col-1];
		$newMValue=$newLineArr[$MVal_col-1];
		$newSign=GiveSign($newMValue);	
		if($prvChr eq $newChr){
		if (abs($newMValue)>=$mValCut/2 && $newSign == $openSign  && $outlierCntr<=$outCut)
		{
			unshift(@outArr,$newLine);#print $newLine;<STDIN>;
			$outlierCntr=0;
		}
#		elsif(abs($newMValue) <= $mValCut/2 && $outlierCntr<=$outCut)
		else
		{
			unshift(@outArr,$newLine);#print $newLine;<STDIN>;
			$outlierCntr++;
		}}
		if($outlierCntr==$outCut+1)
		{
			shift(@outArr);
			shift(@outArr);
			shift(@outArr);
			$outlierCntr=0;
			$newStart=$cntr;
			$openSign="";
                        if(scalar(@outArr)==0)
                        {
	                        push(@outArr,$openningLine);
                        }
		}
		$size=scalar(@outArr);
		}
		if($size>=1)
		{
			foreach (@outArr)
			{
				print "$_\tLocus-$locCntr\n";#<STDIN>;
			}
			@outArr=();
		}
}
	
sub GiveSign()
{
		##check the sign of M.Value		
		$nowMValue=$_[0];
		if ($nowMValue<0)
		{
			$sign=-1;
		}
		else
		{
			$sign=1;
		}
		return $sign;
		$sign="";
}
