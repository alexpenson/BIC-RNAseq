#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 
use Cwd;
use Cwd 'abs_path';
my ($map, $pre, $config, $help, $species, $cufflinks, $dexseq, $htseq, $chimerascan, $samplekey, $comparisons, $deseq, $star_fusion, $mapsplice, $defuse, $fusioncatcher, $detectFusions, $allfusions, $tophat, $star, $pass1, $lncrna, $lincrna_BROAD, $output, $strand, $r1adaptor, $r2adaptor, $transcript, $no_replicates);


my $bpipe_path = "";
my $software_config_file = "$Bin/software_config.txt";
my $annotation_config = "$Bin/annotation_config.groovy";
$output = "results";


GetOptions ('map=s' => \$map,
	        'pre=s' => \$pre,
            'species=s' => \$species,
            'strand=s' => \$strand,
	        'samplekey=s' => \$samplekey,
	        'comparisons=s' => \$comparisons,
	        'h|help' => \$help,
	        'star' => \$star,
	        'pass1' => \$pass1,
	        'tophat|tophat2' => \$tophat,
	        'cufflinks' => \$cufflinks,
	        'dexseq' => \$dexseq,
	        'htseq' => \$htseq,
	        'deseq' => \$deseq,
	        'chimerascan' => \$chimerascan,
	        'star_fusion' => \$star_fusion,
	        'mapsplice' => \$mapsplice,
	        'defuse' => \$defuse,
	        'fusioncatcher' => \$fusioncatcher,
	        'allfusions' => \$allfusions,
	        'transcript' => \$transcript,
            'lncrna' => \$lncrna,
 	        'output|out|o=s' => \$output,
	        'r1adaptor=s' => \$r1adaptor,
	        'r2adaptor=s' => \$r2adaptor,
            'lincrna_BROAD' => \$lincrna_BROAD,
            'no_replicates' => \$no_replicates) or exit(1);


if(!$map || !$pre || !$species || !$strand || $help){
    print <<HELP;

    USAGE: rnaseq_pipeline.pl -map MAP -species SPECIES -strand STRAND -config CONFIG -pre PRE -samplekey SAMPLEKEY -comparisons COMPARISONS -scheduler SCHEDULER
	* MAP: file listing sample mapping information for processing (REQUIRED)
	* SPECIES: only human(hg19), mouse (mm9/mm10) and human-mouse hybrid (hybrid), fly(dm3), zebra fish(zv9) currently supported (REQUIRED)
	* STRAND: library strand; valid options are none, forward, reverse (REQUIRED)
	* PRE: output prefix (REQUIRED)
	* SAMPLEKEY: tab-delimited file listing sampleName in column A and condition in column B (if -deseq, REQUIRED)
	* COMPARISONS: tab-delimited file listing the conditions to compare in columns A/B (if -deseq, REQUIRED)
	* R1ADAPTOR/R2ADAPTOR: if provided, will trim adaptor sequences; NOTE: if provided for only one end, will also assign it to the other end
	* ALIGNERS SUPPORTED: star (-star), defaults to 2pass method unless -pass1 specified; tophat2 (-tophat); if no aligner specifed, will default to STAR
	* ANALYSES SUPPORTED: cufflinks (-cufflinks); htseq (-htseq); dexseq (-dexseq); deseq (-deseq; must specify samplekey and comparisons)
	* FUSION DETECTION: supported fusion callers chimerascan (-chimerascan), star_fusion (-star_fusion), mapsplice (-mapsplice), defuse (-defuse), fusioncatcher (-fusioncatcher); -allfusions will run all supported fusion detection programs
	* TRANSCRIPT ANALYSIS: enable transcript analysis using express and kallisto (-transcript)
	* OUTPUT: output results directory (default: results)
	* OPTIONS: lncRNA analysis (-lncrna) runs all analyses based on lncRNA GTF (hg19 only); 
HELP
exit;
}


if($samplekey || $comparisons)
{
    $deseq = 1;
    $htseq = 1;
}

die "ERROR: Species must be hg19, mm9, mm10, hybrid, dm3, zv9\n" if($species !~ /hg19|mm9|mm10|hybrid|dm3|zv9/);
die "ERROR: -samplekey and -comparisons file are needed for -deseq\n" if ($deseq && (!$samplekey || !$comparisons));
die "ERROR: only human(hg19) support -lncrna\n" if($lncrna && $species !~/hg19/);
die "ERROR: only human(hg19) support fusion detection\n" if(($allfusions || $chimerascan || $star_fusion || $mapsplice || $defuse || $fusioncatcher) && $species !~/hg19/) ;


if($pre =~ /^\d+/){
    $pre = "s_$pre";
}

if($r1adaptor && !$r2adaptor){
	$r2adaptor = $r1adaptor;
}
elsif($r2adaptor && !$r1adaptor){
	$r1adaptor = $r2adaptor;
}

my $htseq_stranded = '';
my $picard_strand_specificity = '';
if($strand =~ /none/i){
    $htseq_stranded = 'no';
    $picard_strand_specificity = 'NONE';
}
elsif($strand =~ /forward/i){
    $htseq_stranded = 'yes';
    $picard_strand_specificity = 'FIRST_READ_TRANSCRIPTION_STRAND';
}
elsif($strand =~ /reverse/i){
    $htseq_stranded = 'reverse';
    $picard_strand_specificity = 'SECOND_READ_TRANSCRIPTION_STRAND';
}
else
{
    die "ERROR: unrecognized stand specified $strand\n";
}




### prepare file/directory
die "[ERROR]: Mapping file does not exist: $map\n" if(!-e $map);
die "[ERROR]: Samplekey file does not exist: $samplekey\n" if($samplekey && !-e $samplekey);
die "[ERROR]: Comparisons file file does not exist: $comparisons\n" if($comparisons && !-e $comparisons);
$map = abs_path($map);
$samplekey = abs_path($samplekey);
$comparisons = abs_path($comparisons);
`mkdir -p $output`;
die "[ERROR]: Could not create output directory: $output\n" if(!-d $output);
$output = abs_path($output);
`cp $Bin/rnaseq_pipeline.bpipe $output`;
`cp $Bin/bpipe.config $output`;
`sed -i s/PROJECT_NAME_TOBE_REPLACED/$pre/g $output/bpipe.config`;




### get bpipe path from software config file
open(SC, "$software_config_file") or die "Can't open software config file $software_config_file $!";
while(<SC>){
    chomp;
    my @conf = split(/\s+/, $_);
    if($conf[0] =~ /BPIPE/i){
        if(!-e "$conf[1]/bin/bpipe"){
            die "[ERROR]: Can not find bpipe executable IN $conf[1]/bin/ $!";
        }
        $bpipe_path = $conf[1];
    }
}
close SC;
die "[ERROR]: Could not find bpipe path in software config file: $software_config_file\n" if(!$bpipe_path);




### check if directories in mapping file exist
my %mapping_samples = ();
open(MA, "$map") or die "Can't open mapping file $map $!";
while(<MA>){
    chomp;
    my @data = split(/\s+/, $_);
    $mapping_samples{$data[1]} = 1;
    if(!-d $data[3]){
        die "$data[3] in mapping file does not exist\n";
    }
}
close MA;




### check samplekey file and comparisons files
if($deseq){
    my %sample_comparisons = ();
    open(SC, "$comparisons") or die "Can't open comparisons file $comparisons $!";
    while(<SC>){
        chomp;
        my @data = split(/\s+/, $_);
        $sample_comparisons{$data[0]} = 1;
        $sample_comparisons{$data[1]} = 1;
    }
    close SC;
    my %samplekey_samples = ();
    my %samplekey_conditions = ();
    open(SK, "$samplekey") or die "Can't open key file $samplekey $!";
    while(<SK>){
        chomp;
        my @data = split(/\s+/, $_);
        $samplekey_samples{$data[0]} = 1;
        $samplekey_conditions{$data[1]} = 1;
        #print "data[0]: $data[0]\tdata[1]: $data[1]\n";
        if(!$mapping_samples{$data[0]} || !$sample_comparisons{$data[1]}){
            die "either sample $data[0] cannot be found in $map and/or condition $data[1] cannot be found in $comparisons $!";
        }
    }
    close SK;
    foreach my $ms (keys %mapping_samples){
        if(!$samplekey_samples{$ms}){
            die "sample $ms is in mapping file $map, but not in sample key file $samplekey $!";
        }
    }
    foreach my $sc (keys %sample_comparisons){
        if(!$samplekey_conditions{$sc}){
            die "condition $sc is in sample comparisons file $comparisons, but not in sample key file $samplekey $!";
        }
    }
}



### call RNA bpipe pipeline
my $extra_para = "";
if($r1adaptor)
{
    $extra_para .=  " -p flag_trim_read=true -p r1adaptor=$r1adaptor -p r2adaptor=$r2adaptor";
}
if($tophat)
{
    $extra_para .= " -p flag_aligner=tophat";
}
if($pass1)
{
    $extra_para .= " -p flag_star_2p=false";
}
if($cufflinks)
{
    $extra_para .= " -p flag_cufflinks=true";
}
if($htseq)
{
    $extra_para .= " -p flag_htseq=true";
}
if($dexseq)
{
    $extra_para .= " -p flag_dexseq=true";
}
if($deseq)
{
    $extra_para .= " -p flag_deseq=true -p COMPARISON_FILE=$comparisons -p KEY_FILE=$samplekey ";
}
if($no_replicates)
{
    $extra_para .= " -p flag_no_replicates=true";
}
if($allfusions || $chimerascan || $star_fusion || $mapsplice || $defuse || $fusioncatcher)
{
    $extra_para .= " -p flag_detectfusion=true";
}
if($allfusions)
{
    $extra_para .= " -p flag_fusion_chimerascan=true -p flag_fusion_star=true -p flag_fusion_mapsplice=true -p flag_fusion_defuse=true -p flag_fusion_fusioncatcher=true";
}
if($chimerascan)
{
    $extra_para .= " -p flag_fusion_chimerascan=true";
}
if($star_fusion)
{
    $extra_para .= " -p flag_fusion_star=true";
}
if($mapsplice)
{
    $extra_para .= " -p flag_fusion_mapsplice=true";
}
if($defuse)
{  
    $extra_para .= " -p flag_fusion_defuse=true";
}
if($fusioncatcher)
{
    $extra_para .= " -p flag_fusion_fusioncatcher=true";
}
if($transcript)
{
    $extra_para .= " -p flag_transcript=true";
}
if($lncrna)
{
    $extra_para .= " -p flag_lncrna=true";
}

chdir $output;

`BPIPE_BACKGROUND=1 $bpipe_path/bin/bpipe run -p Bin=$Bin -p species=$species -p pre=$pre -p htseq_stranded=$htseq_stranded -p picard_strand_specificity=$picard_strand_specificity -p MAPPING_FILE=$map -p SOFTWARE_CONFIG_FILE=$software_config_file -p ANNOTATION_CONFIG_FILE=$annotation_config $extra_para rnaseq_pipeline.bpipe`;

print "Start running RNA-seq bpipe workflow\n";

























