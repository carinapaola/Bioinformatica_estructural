#!/usr/bin/perl
# prog1.1 
# Carina Paola Cornejo Paramo, Claudia Sarai Reyes Avila
use strict;
=pod
=head1 DESCRIPCION
Este programa implementa el predictor de promotores de Kanhere y Bansal, que se basa en el hecho de que hay un cambio en la energia libre de la secuencia en la region del promotor ya que disminuye la estabilidad en esta zona pues es el lugar en el que deben abrirse las cadenas de DNA para iniciar la transcripcion. http://www.biomedcentral.com/1471-2105/6/1. El script requiere un parametro, el nombre del archivo que contiene las secuencias en las que se buscaran los promotores. Devuelve un archivo con las coordenadas de los promotores que fueron predichos por el programa.
=head2 VARIABLES
$j, $i, $k	: son contadores usados en las iteraciones (for, foreach).
$cutoff1	:variable que almacena el valor de corte establecido para la diferencia ($Dif) entre los valores de E1 y E2.
$cutoff2	:almacena el valor de corte para E1.
%hseqs		:hash que almancena los identificadores y sus respectivas secuencias extraidos del archivo que se proporciona como parametro.
$T	: temperatura en grados Celcius.
$window	:variable que almacena una subcadena de longitud $windowL, a partir de las cuales se calcula el delta G.
$dG_nt	:almacena el delta G de un nucleotido, calculado a partir del delta G de una ventana ($window).
@dG_seq	:arreglo que almacena temporalmente los valores de delta G para todos los nucleotidos de una secuencia.
$sumE1	:almacena la sumatoria de los valores de delta G para un segmento de 50 nucleotidos de las secuancias, con lo que posteriormente se calcula E1. 
$sumE2	:almacena la sumatoria de los valores de delta G para un segmento de 100 nucleotidos de las secuencias, a partir de este valor de calcula E2.
$E1	:almacena el promedio de delta G de un segmento de 50 nucleotidos para la posicion n, el segmento tiene coordenadas n,(n+49).	 $E2	:almacena el promedio de delta G de un segmento de 100 nucleotidos para la posicion n, el segmento tiene coordenadas n+99,(n+199). $Dif	:almacena la diferencia entre los valores de E1 y E2 para cada nucleotido.
%n_values	:hash cuyas llaves son los identificadores de las secuencias originales y los valores son arreglos que almacenan las posiciones n que pasaron los valores de corte cutoff1 y cutoff2.
@significant_n	:arreglo que almacena temporalmente las posiciones n de una secuencia que pasaron los valores de corte cutoff1 y cutoff2.
$windowL	:longitud de las ventanas a partir de las cuales se calculan los valores de delta G.	
%dG_allSequences	:hash que almacena los valores de deltaG para todas las posiciones de las secuencias.
%NNparams	:hash que almacena los valores de entropia y entalpia calculados experimentalmente para cada par de dinucleotidos asi como los valores por costos de iniciacion y para secuencias palindromicas.	
=cut
# variables globales
my @N		=();
my $j		=0;
my $i		= 0;
my $k		=0;
my $cutoff1	=3.4;
my $cutoff2	= -15.99; 
my $coordenadas=0;
my %hseqs	= ();
my $T           = 37; 
my $window	=();	
my $dG_nt  	= 0;
my @dG_seq	= ();
my $N		=0;
my $sumE1	=0;
my $sumE2	=0;
my $E1		=0;
my $E2		=0;
my $Dif		=0;
my @E1_seq	=();
my @E2_seq	=();
my @Dif_seq	=();
my %n_values= ();
my %all_E1	=();
my %all_E2	=();
my %all_Dif	=();
my @final_nvalues=();
my %hfinal_nvalues=();
my @significant_n	=();
my $windowL     	= 15; 
my %dG_allSequences	=	();
my %NNparams    = ( 
	# SantaLucia J (1998) PNAS 95(4): 1460-1465.
	# [NaCl] 1M, 37C & pH=7 
	# H(enthalpy): kcal/mol	, S(entropy): cal/k�mol
	# stacking dinucleotides
	'AA/TT' , {'H',-7.9, 'S',-22.2},
	'AT/TA' , {'H',-7.2, 'S',-20.4},
	'TA/AT' , {'H',-7.2, 'S',-21.3},
	'CA/GT' , {'H',-8.5, 'S',-22.7},
	'GT/CA' , {'H',-8.4, 'S',-22.4},
	'CT/GA' , {'H',-7.8, 'S',-21.0},
	'GA/CT' , {'H',-8.2, 'S',-22.2},
	'CG/GC' , {'H',-10.6,'S',-27.2},
	'GC/CG' , {'H',-9.8, 'S',-24.4},
	'GG/CC' , {'H',-8.0, 'S',-19.9},
	# costos de iniciacion
	'G'     , {'H', 0.1, 'S',-2.8 },
	'A'     , {'H', 2.3, 'S',4.1  },
	# correccion por simetria (palindromos)
	'sym'   , {'H',   0, 'S',-1.4 } );
my $infile = $ARGV[0] || die "# usage: $0 <promoters file>\n";
print "# parameters: Temperature=$T Window=$windowL\n\n";  
open(SEQ, $infile) || die "# cannot open input $infile : $!\n";
open(MIFICH,">salida.txt")|| die "lo siento, no puedo crear fichero.txt\n"; 
#$seq almacena temporalmente la secuencia de nucleotidos leida a partir del archivo que se da como parametro al programa.
#$name almacena temporalmente el identificador de la secuencia.
#Leyendo el archivo por lineas y extrayendo identificadores y secuencias.
while(<SEQ>){
	if(/^(b\d{4}) \\ ([ATGC]+)/){
		my ($name,$seq) = ($1,$2); 
		printf("sequence %s (%d nts)\n",$name,length($seq));
		# Guardando las secuencias y sus identificadores en el hash $hseqs 
		$hseqs{$name}=$seq;}}
close(SEQ);
#calculando ventanas para cada secuencia guardada
foreach $j(keys %hseqs){
	for ($i=0; ($i+14)<length($hseqs{$j});$i++){
		$window=substr($hseqs{$j},$i,15);
		$dG_nt=duplex_deltaG($window,37);#Guardando deltaG de la ventana, que es el dG del nucleotido devuelto por la subrutina.
		push(@dG_seq, $dG_nt);}#guardando en el arreglo @dG_seq el delta G de los nucleotidos para una secuencia.
	$dG_allSequences{$j}=[@dG_seq];#almacenando los delta G de cada nucleotido para todas las secuencias.
	undef(@dG_seq);}
#calculando valores de E1, E2 y diferencia para cada posición i
foreach $j (keys %dG_allSequences){#recorriendo las secuencias
	for($i=0; ($i+199)<=scalar(@{$dG_allSequences{$j}});$i++){
		#calculando E1
		for($k=$i; $k<=($i+49); $k++){
			$sumE1+=$dG_allSequences{$j}[$k];}
		$E1=$sumE1/50;
		push(@E1_seq, $E1);#Guardanto E1
		#calculando E2
		for($k=$i+99;$k<=($i+199);$k++){
			$sumE2+=$dG_allSequences{$j}[$k];}
		$E2=$sumE2/100;
		push(@E2_seq,$E2);#guardando E2
		$Dif=$E1-$E2;#calculando diferencia entre E1 y E2
		if($Dif>$cutoff1&&$E1>$cutoff2){#evaluando si los valores de E1 y diferencia de una posicion sobrepasan los umbrales
			push(@significant_n, $i);}#almacenando las posiciones de una secuancia que pasan los umbrales
		push(@Dif_seq,$Dif);
		#reiniciando variables
		$sumE1=0;
		$sumE2=0;
		$E1=0;
		$E2=0;
		$Dif=0;}
	#almacenando en el hash %n_values todas las posiciones que pasaron el umbral para cada secuencia.
	$n_values{$j}=[@significant_n];
	$all_E1{$j}=[@E1_seq];
	$all_E2{$j}=[@E2_seq];
	$all_Dif{$j}=[@Dif_seq];
	undef(@Dif_seq);
	undef(@E1_seq);
	undef(@E2_seq);
	undef(@significant_n);}
foreach $j(keys %n_values){#identificando las posiciones significativas de cada segmento de 25 nucleotidos de valores que pasaron los umbrales
	if(defined($n_values{$j})){
		for($i=0;$i<(scalar(@{$n_values{$j}})-1);$i++){
			if($n_values{$j}[$i+1]<($n_values{$j}[$i]+25)){
				$N=$n_values{$j}[$i+1];
			}else{
				$N=$n_values{$j}[$i];
				push(@N,$N);
				$N=$n_values{$j}[$i+1];}}
		push(@N,$N);
		$hfinal_nvalues{$j}=[@N];
		undef(@N);}}
#guardando resultados en el archivo de salida para posteriormente graficar.
print MIFICH "VALORES DE E1\n";
foreach $i (keys %all_E1){
	print MIFICH "$i\n";
	print MIFICH "@{$all_E1{$i}}\n";}
print MIFICH "VALORES DE E2\n";
foreach $i (keys %all_E2){
	print MIFICH "$i\n";
	print MIFICH "@{$all_E2{$i}}";
	print MIFICH "\n";}
print MIFICH "VALORES DE DIFERENCIA\n";
foreach $i (keys %all_Dif){
	print MIFICH "$i\n";
	print MIFICH "@{$all_Dif{$i}}\n";}
print MIFICH "VALORES DE N QUE PASAN LOS UMBRALES\n";
foreach $i (keys %n_values){
	print MIFICH "$i\n";
	print MIFICH "@{$n_values{$i}}\n";}
print MIFICH "POSICIONES SIGNIFICATIVAS DE CADA SEGMENTO\n";
foreach $i (keys %hfinal_nvalues){
	print MIFICH "$i\t";
	print MIFICH "@{$hfinal_nvalues{$i}}\n";}
print MIFICH "VALORES QUE CAEN EN LA REGION PROMOTORA\n";
foreach $i (keys %hfinal_nvalues){
	for($j=0;$j<scalar(@{$hfinal_nvalues{$i}});$j++){
		if (@{$hfinal_nvalues{$i}}[$j]>200&&@{$hfinal_nvalues{$i}}[$j]<400){
			print MIFICH "$i\t";
			print MIFICH "$hfinal_nvalues{$i}[$j]\n";}}}
close(MIFICH);
# calcula la energia libre para un NN del duplex de DNA , dG(t) = (1000*dH - t*dS) / 1000
# parametros: 1) cadena con una secuencia de DNA; 2) Temperatura en grados Celsius
# regresa; 1) escalar con un valor de energia libre.
# usa el hash global %NNparams
sub duplex_deltaG{
   	my ($seq,$tCelsius) = @_; 
	my ($DNAstep,$nt,$dG,$total_dG) = ('','',0,0);
	my @sequence = split(//,uc($seq));
	my $tK = 273.15 + $tCelsius;
	my $seq_temporal=0;
	sub complement{ $_[0]=~ tr/ATGC/TACG/; return  $_[0]; }
	# add dG for overlapping dinculeotides
	for(my $n=0;$n<$#sequence;$n++){
			$seq_temporal=$sequence[$n].$sequence[$n+1];
			$seq_temporal=complement($seq_temporal);
			$DNAstep = $sequence[$n].$sequence[$n+1].'/'.$seq_temporal;
			if(!defined($NNparams{$DNAstep})){
				$DNAstep = reverse($DNAstep);}
			$dG = ((1000*$NNparams{$DNAstep}{'H'})-($tK*$NNparams{$DNAstep}{'S'}))/1000;
			$total_dG += $dG; }
	# add correction for helix initiation
	$nt = $sequence[0]; # first pair
	if(!defined($NNparams{$nt})){ $nt = complement($nt); } 
	$total_dG += ((1000*$NNparams{$nt}{'H'})-($tK*$NNparams{$nt}{'S'}))/1000; 
	$nt = $sequence[$#sequence]; # last pair
	if(!defined($NNparams{$nt})){ $nt = complement($nt); }
	$total_dG += ((1000*$NNparams{$nt}{'H'})-($tK*$NNparams{$nt}{'S'}))/1000;		
	#symetry correction
	while ($inicio<7&&$final>7){
		$temporal=$sequence[$final];
		if($sequence[$inicio] ne complement($temporal)){
			$sim=0;}
		$inicio++;
		$final--;
		}
	if($sim==1){$total_dG+=((1000*$NNparams{'sym'}{'H'})-($tK*$NNparams{'sym'}{'S'}))/1000;}
	return $total_dG;}
