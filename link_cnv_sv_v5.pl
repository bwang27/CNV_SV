#!/usr/bin/perl -w

##usage: perl link_cnv_sv_v5.pl -s sv -c cnv -o out.txt
#
# $header="sample\ttype\tmarker\tcolor\tchromosome\tstart\tend\tsv1.dis1\tsv1.dis2\tsv1\ttype1\tsv2.dis1\tsv2.dis2\tsv2\ttype2\n";

use Getopt::Std; 

getopt("sco",\%args);

$sv=$args{s}?$args{s}:"sv";
$cnv=$args{c}?$args{c}:"cnv";
$out=$args{o}?$args{o}:"out.txt";

open SV,"$sv"||die $!;
$all=join'',<SV>;
@sv=split(/\n/,$all);
close SV;

open CNV,"$cnv"||die $!;
open CS,">$out" || die $!;;
my @cnvsv;
my @final_sv;
while(<CNV>){
	chomp;
	print CS "$_";
	my %sv_type;
	@aa=split(/\t/,$_);
	@cnv_sj=split(/\_/,$aa[0]);
	$cnv_sjid=$cnv_sj[0];
	my($cnv_chr,$cnv_start,$cnv_end)=split(/\./,$aa[2]);
	foreach $line(@sv){
		#my %sv_type;
		@bb=split(/\t/,$line);
		$sv_sjid=$bb[0];
		my($svbk1_chr,$sv_start,$sv1_strand,$svbk2_chr,$sv_end,$sv2_strand,$seq)=split(/\./,$bb[2]);
		if($cnv_sjid eq $sv_sjid){
			if(($cnv_chr eq $svbk1_chr)&&($cnv_chr ne $svbk2_chr)){
				$dis_svstart_cnvstart=abs($cnv_start-$sv_start);
				$dis_svstart_cnvend=abs($cnv_end-$sv_start);
				push@itssv,$dis_svstart_cnvstart,$dis_svstart_cnvend,$bb[2];
				$onesv=join(":",@itssv);
				$B="B-type";
				$sv_type{$onesv}=$B;
				push @final_sv,$onesv,"B-type";
				#print CS "\t$onesv";
				undef @itssv;
				#print CS "\tB-type";
			}
			if(($cnv_chr eq $svbk2_chr)&&($cnv_chr ne $svbk1_chr)){
				$dis_svend_cnvstart=abs($cnv_start-$sv_end);
				$dis_svend_cnvend=abs($cnv_end-$sv_end);
				push@itssv,$dis_svend_cnvstart,$dis_svend_cnvend,$bb[2];
				$onesv=join(":",@itssv);
				$B="B-type";
				$sv_type{$onesv}=$B;
				push @final_sv,$onesv,"B-type";
				#print CS "\t$onesv";
				undef @itssv;
				#print CS "\tB-type";
			}
			if(($cnv_chr eq $svbk1_chr)&&($cnv_chr eq $svbk2_chr)){
				$dis_svstart_cnvstart=abs($cnv_start-$sv_start);
				$dis_svend_cnvend=abs($cnv_end-$sv_end);
				$dis_svstart_cnvend=abs($cnv_end-$sv_start);
				$dis_svend_cnvstart=abs($cnv_start-$sv_end);
				@myids=($dis_svstart_cnvstart,$dis_svend_cnvend,$dis_svstart_cnvend,$dis_svend_cnvstart);
				@st_nd=sort { $a <=> $b} @myids;
				push@itssv,$st_nd[0],$st_nd[1],$bb[2];
				$onesv=join(":",@itssv);
				$A="A-type";
				$sv_type{$onesv}=$A;
				push @final_sv,$onesv,"A-type";
				#print CS "\t$onesv";
                                undef @itssv;
				#print CS "\tA-type";
			}
	}
}
	my %rehas;
	my @alldis;
	if(@final_sv<3){
		foreach $ele(keys%sv_type){	
			#$new=join("\t",$ele,$sv_type{$ele});
			@sort_pos=split(/\:/,$ele);
			#@sort_pos=grep {defined && /\S/ } @sort_pos;
			if($sort_pos[0]<$sort_pos[1]){
				print CS "\t$sort_pos[0]\t$sort_pos[1]\t$sort_pos[2]\t$sv_type{$ele}";
			}	
			else{
				print CS "\t$sort_pos[1]\t$sort_pos[0]\t$sort_pos[2]\t$sv_type{$ele}";#$la"; #\t$sort_pos[2]$type";}
			}	
		}
	}
    	elsif(@final_sv>2){
		my $sv1;
		my $sv2;
		foreach $dis(keys%sv_type){
			$con=join("\t",$dis,$sv_type{$dis});
			@ss=split(/\:/,$dis);
			@first_two=splice(@ss,0,2);
			@sorted_ftwo=sort {$a <=> $b} @first_two;
			$min=$sorted_ftwo[0];
			$rehas{$min}=$con;
			push @alldis,$min;
		}
		my @sorted_ascending=sort {$a <=> $b} @alldis;
		$sv1=$rehas{$sorted_ascending[0]}; #=$sv1;
		$sv2=$rehas{$sorted_ascending[1]};#=$sv2;
		$sv1=~s/\:/\t/g;
		$sv2=~s/\:/\t/g;
		if($sv1 eq $sv2){
			@sort_sv1pos=split(/\t/,$sv1);
			if($sort_sv1pos[0]<$sort_sv1pos[1]){
				print CS "\t$sort_sv1pos[0]\t$sort_sv1pos[1]\t$sort_sv1pos[2]\t$sort_sv1pos[3]";}
			else{
				print CS "\t$sort_sv1pos[1]\t$sort_sv1pos[0]\t$sort_sv1pos[2]\t$sort_sv1pos[3]";}
		}
		else{
			@sort_sv1pos=split(/\t/,$sv1);
			@sort_sv2pos=split(/\t/,$sv2);
			if($sort_sv1pos[0]<$sort_sv1pos[1]){
                                print CS "\t$sort_sv1pos[0]\t$sort_sv1pos[1]\t$sort_sv1pos[2]\t$sort_sv1pos[3]";}
			else{
				print CS "\t$sort_sv1pos[1]\t$sort_sv1pos[0]\t$sort_sv1pos[2]\t$sort_sv1pos[3]";}
			if($sort_sv2pos[0]<$sort_sv2pos[1]){
                                print CS "\t$sort_sv2pos[0]\t$sort_sv2pos[1]\t$sort_sv2pos[2]\t$sort_sv1pos[3]";}
			else{
				print CS "\t$sort_sv2pos[1]\t$sort_sv2pos[0]\t$sort_sv2pos[2]\t$sort_sv1pos[3]";}
		}
	}
	#my %seen;
	#my@uniq_arr=grep {!$seen{$_}++} @itssv;
	#$mysv=join(";",@uniq_arr);
	#$mysv=join(":",@itssv);
	print CS "\n";
	#undef %seen;
	undef @itssv;
	undef @alldis;
	undef %rehas;
	undef %sv_type;
	undef @final_sv;
}

close CNV;
close CS;

