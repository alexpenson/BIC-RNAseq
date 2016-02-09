#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 
use Cwd;
use Cwd 'abs_path';
my ($map, $pre, $config, $help, $species, $cufflinks, $dexseq, $htseq, $chimerascan, $samplekey, $comparisons, $deseq, $star_fusion, $mapsplice, $defuse, $fusioncatcher, $detectFusions, $allfusions, $tophat, $star, $pass1, $lncrna, $lincrna_BROAD, $output, $strand, $r1adaptor, $r2adaptor, $transcript, $no_replicates, $max_bpipe_job, $rsem, $kallisto, $express, $standard, $rsync_path);


$max_bpipe_job = 1000;
my $uID = `/usr/bin/id -u -n`;
chomp $uID;
my $software_config_file = "$Bin/software_config.txt";
my $annotation_config = "$Bin/annotation_config.groovy";
$rsync_path = "/ifs/solres/$uID/";
$output = "results";
my $bpipe_path = "";



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
            'kallisto' => \$kallisto,
            'rsem' => \$rsem,
            'express' => \$express,
            'standard' => \$standard,
            'lncrna' => \$lncrna,
 	        'output|out|o=s' => \$output,
	        'r1adaptor=s' => \$r1adaptor,
	        'r2adaptor=s' => \$r2adaptor,
            'lincrna_BROAD' => \$lincrna_BROAD,
            'no_replicates' => \$no_replicates,
            'max_bpipe_job=i' => \$max_bpipe_job,
            'rsync=s' => \$rsync_path) or exit(1);


if(!$map || !$pre || !$species || !$strand || $help){
    print <<HELP;

    USAGE: rnaseq_pipeline.pl -map MAP -species SPECIES -strand STRAND -pre PRE -samplekey SAMPLEKEY -comparisons COMPARISONS
	* MAP: file listing sample mapping information for processing (REQUIRED)
	* SPECIES: only human(hg19), mouse (mm9/mm10) and human-mouse hybrid (hybrid), fly(dm3), zebra fish(zv9) currently supported (REQUIRED)
	* STRAND: library strand; valid options are none, forward, reverse (REQUIRED)
	* PRE: output prefix (REQUIRED)
	* SAMPLEKEY: tab-delimited file listing sampleName in column A and condition in column B (if -deseq, REQUIRED)
	* COMPARISONS: tab-delimited file listing the conditions to compare in columns A/B (if -deseq, REQUIRED)
	* R1ADAPTOR/R2ADAPTOR: if provided, will trim adaptor sequences; NOTE: if provided for only one end, will also assign it to the other end
	* ALIGNERS SUPPORTED: star (-star), defaults to 2pass method unless (-pass1) is specified; tophat2 (-tophat); if no aligner specifed, will default to STAR
	* ANALYSES SUPPORTED: cufflinks (-cufflinks); htseq (-htseq); dexseq (-dexseq); deseq (-deseq; must specify samplekey and comparisons)
	* FUSION DETECTION: supported fusion callers chimerascan (-chimerascan), star_fusion (-star_fusion), mapsplice (-mapsplice), defuse (-defuse), fusioncatcher (-fusioncatcher); -allfusions will run all supported fusion detection programs
	* TRANSCRIPT ANALYSIS:  transcript analysis using express (-express), kallisto (-kallisto) and rsem (-rsem), or all (-transcript)
	* STANDARD: standard analysis (-standard): star alignment, htseq gene count, rsem and kallisto transcript counts, counts normalization, and clustering
    * OUTPUT: output results directory (default: results)
	* OPTIONS: lncRNA analysis (-lncrna) runs all analyses based on lncRNA GTF (hg19 only); 
	* max_bpipe_job: use -max_bpipe_job <int> to change the maximum number of cluster jobs that the current bpipe job is allowed to run at the same time. Defualt to 1,000
    * -rsync: path to rsync after the pipeline finishes. Default is /ifs/solres/USER_ID/. Use -rsync NULL to disable automatic rsync.
HELP
exit;
}


if($standard)
{                                                                                                                                                            
    $star = 1;                                                                                                                                                      
    $htseq = 1;                                                                                                                                                           
    $kallisto = 1;                                                                                                                                                        
    $rsem = 1;                                                                                                                                                            
}

if($samplekey || $comparisons)
{
    $htseq = 1;
    $deseq = 1;
}

if(($star || $tophat) && !$htseq && !$dexseq){                                                                                                                                         
    $htseq = 1;                                                                                                                                                           
}

if($transcript){                                                                                                                                                          
    $kallisto = 1;                                                                                                                                                        
    $rsem = 1;                                                                                                                                                            
    ###$express = 1;                                                                                                                                                      
}    

if($allfusions)
{
    $chimerascan = 1;
    $star_fusion = 1;
    $mapsplice = 1;
    $defuse = 1;
    $fusioncatcher = 1;
}

if(($deseq || $dexseq || $htseq || $cufflinks) && (!$star && !$tophat) ){
    $star = 1;
}

die "ERROR: -max_bpipe_job must >= 1\n" if($max_bpipe_job < 1);
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
die "[ERROR]: -rsync path does not exist: $rsync_path. Or use -rsync NULL to disable automatic rsync\n" if($rsync_path !~ /^null$/i && !-d $rsync_path);
$map = abs_path($map);
$samplekey = abs_path($samplekey);
$comparisons = abs_path($comparisons);
if($rsync_path !~ /^null$/i)
{
    $rsync_path = abs_path($rsync_path);
}
`mkdir -p $output`;
die "[ERROR]: Could not create output directory: $output\n" if(!-d $output);
$output = abs_path($output);
`cp $Bin/rnaseq_pipeline.bpipe $output`;
`cp $Bin/bpipe.config $output`;
`sed -i s/PROJECT_NAME_TOBE_REPLACED/$pre/g $output/bpipe.config`;
`sed -i s/USERID_TOBE_REPLACED/$uID/g $output/bpipe.config`;



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
`ln -s -f $bpipe_path/bin/bpipe $output/`;



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
if($star)
{
    $extra_para .= " -p flag_star=true";
}
if($tophat)
{
    $extra_para .= " -p flag_tophat=true";
}
if($pass1)
{
    $extra_para .= " -p flag_star_1p=true";
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
if($chimerascan || $star_fusion || $mapsplice || $defuse || $fusioncatcher)
{
    $extra_para .= " -p flag_detectfusion=true";
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
if($kallisto || $rsem || $express)
{
    $extra_para .= " -p flag_transcript=true";
}
if($kallisto)
{
    $extra_para .= " -p flag_kallisto=true";
}
if($rsem)
{
    $extra_para .= " -p flag_rsem=true";
}
if($express)
{
    $extra_para .= " -p flag_express=true";
}
if($lncrna)
{
    $extra_para .= " -p flag_lncrna=true";
}
if($rsync_path !~ /^null$/i)
{
    $extra_para .= " -p rsync_path=$rsync_path ";
}

chdir $output;

`BPIPE_BACKGROUND=1 MAX_JAVA_MEM=2g $bpipe_path/bin/bpipe run -n $max_bpipe_job -p Bin=$Bin -p species=$species -p pre=$pre -p htseq_stranded=$htseq_stranded -p picard_strand_specificity=$picard_strand_specificity -p MAPPING_FILE=$map -p SOFTWARE_CONFIG_FILE=$software_config_file -p ANNOTATION_CONFIG_FILE=$annotation_config $extra_para rnaseq_pipeline.bpipe`;

print "\n";
print "-------------------------------------------------------------------------\n";
print "Start running RNA-seq bpipe workflow\n";
print "To stop the pipeline, go to the result directory and do \"./bpipe stop\"\n";
print "To check the log, go to the result directory and do \"./bpipe log\"\n";
print "-------------------------------------------------------------------------\n";
print "\n";
























