#!/usr/bin/perl -w

# parseoproteinas.pl
# Carina Paola Cornejo Paramo & Saraí Reyes

use strict;

# global variables
my $aa= ();
my $n= ();
my $j=0;
my $i= ();
my $a= ();
my $b= ();
my $segmento= ();
my @secuencia_1letra= ();
my @secuencia_3letras= ();
my %codigo    = ( 
	'ALA' => 'A',
	'ARG' => 'R',
	'ASN' => 'N',
	'ASP' => 'D',
	'CYS' => 'C',
	'GLN' => 'Q',
	'GLU' => 'E',
	'GLY' => 'G',
	'HIS' => 'H',
	'ILE' => 'I',
	'LEU' => 'L',
	'LYS' => 'K',
	'MET' => 'M',
	'PHE' => 'F',
	'PRO' => 'P',
	'SER' => 'S',
	'THR' => 'T',
	'TRP' => 'W',
	'TYR' => 'Y',
	'VAL' => 'V');
my $infile = $ARGV[0] || die "# usage: $0 <promoters file>\n";
open(SEQ, $infile) || die "# cannot open input $infile : $!\n";
while(<SEQ>){
	#print"entra a while\n";
	if(/^ATOM\s+\d{1,4}\s+CA\s+(\w{3})/)
	{
		$aa=$1;
		push (@secuencia_3letras,$aa);
		push (@secuencia_1letra, $codigo{$aa})
	}
	
}
close(SEQ);
print"Secuencia con codigo de 3 letras\n";
foreach $n(@secuencia_3letras){
	print "$n ";
}
print"\nSecuencia con codigo de 1 letra\n";
foreach $i(@secuencia_1letra){
	$j=lc($i);
	print "$j";
}
print"\n";

$a=scalar(@secuencia_3letras);
$b=scalar(@secuencia_1letra);
print"$a\n";
print"$b\n";

my $infile2 = $ARGV[1] || die "# usage: $0 <promoters file>\n";
open(SEQ2, $infile2) || die "# cannot open input $infile : $!\n";
while(<SEQ2>){
	chomp;
	print"while1";
	}
while(<SEQ2>){
	print"while2";
	if(/^\w+/){print"$1";}
}
	
#print "$segmento";



















