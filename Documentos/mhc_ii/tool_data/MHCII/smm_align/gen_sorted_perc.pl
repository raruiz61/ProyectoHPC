open(IN, "HLA_list.txt") or die "dd";
while($line=<IN>){
    if($line=~/HLA_(\S+)-(\S+)/){
	$allele=$1."_".$2;
    }
    elsif($line=~/(\S+)/){
	$allele=$1;
    }
	$scan="/home/pwang/programs/smm_align/netMHCII-1.0a/smm_scan_".$allele.".txt";
    $line=~/(\S+)/;
    $oo="perc_".$1.".txt";
    open(SW, "$scan") or die "$scan";
    @arr=();
    $temp=0;
    $cutoff=1000000;
    while($ll=<SW> and $temp<$cutoff){
	if($ll=~/(\S+)/){
	    push(@arr, $1);
	    $temp++;
	}
    }
    @sorted=sort {$a <=> $b} @arr;
    open(OUT, ">$oo") or die "$oo";
    $step=$cutoff/10000;
    for($i=0; $i<=$#sorted;$i=$i+$step){
	print OUT "$sorted[$i]\n";
    }
}

