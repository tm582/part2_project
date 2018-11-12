#!/usr/bin/perl

if (!$ARGV[0]){
	die "Error: Please provide a GTF file"
}

open(FILE,$ARGV[0]);
while(<FILE>){
chomp;

if (/gene_id \"([^"]+)\"/){
	$gene=$1;
	if (/gene_source \"([^"]+)\"/){
		$source=$1;
		$sources{$gene}=$source;
	}
	if (/gene_name \"([^"]+)\"/){
		$name=$1;
		$names{$gene}=$name;
	}
	if (/gene_biotype \"([^"]+)\"/){
		$biotype=$1;
		$biotypes{$gene}=$biotype;
	}
}


}

foreach $thing (keys(%names)){
	print "$thing\t$names{$thing}\t$sources{$thing}\t$biotypes{$thing}\n";
}

