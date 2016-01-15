environments {
 hg19 {
    REF_SEQ = "/ifs/depot/assemblies/H.sapiens/hg19/hg19.fasta"
    DEXSEQ_GTF = "BIN_DIRECTORY_TOBE_SUBSTITUTE/data/gencode.v18.annotation_dexseq.gtf"
    CHIMERASCAN_INDEX = "/ifs/depot/assemblies/H.sapiens/hg19/index/chimerascan/0.4.5a"
    BOWTIE_INDEX = "/ifs/depot/assemblies/H.sapiens/hg19/index/bowtie/1.0.0/hg19_bowtie"
    BOWTIE2_INDEX = "/ifs/depot/assemblies/H.sapiens/hg19/index/bowtie/2.2.4/hg19_bowtie2"
    chrSplits = "/ifs/depot/assemblies/H.sapiens/hg19/chromosomes"
    RIBOSOMAL_INTERVALS = "BIN_DIRECTORY_TOBE_SUBSTITUTE/data/ribosomal_hg19.interval_file"
    REF_FLAT = "BIN_DIRECTORY_TOBE_SUBSTITUTE/data/refFlat__hg19.txt.gz"
    TRANS_INDEX = "/ifs/depot/assemblies/H.sapiens/hg19/index/bowtie/2.2.4/transcriptome/gencode/v18/gencode.v18.annotation"
    TRANS_INDEX_DEDUP = "/ifs/depot/assemblies/H.sapiens/hg19/index/bowtie/2.2.4/transcriptome/gencode/v18/deduplicated/gencode.v18.annotation.dedup"
    TRANS_FASTA_DEDUP = "/ifs/depot/assemblies/H.sapiens/hg19/index/bowtie/2.2.4/transcriptome/gencode/v18/deduplicated/gencode.v18.annotation.dedup.fa"
    geneNameConversion = "BIN_DIRECTORY_TOBE_SUBSTITUTE/data/gencode18IDToGeneName.txt"
    KALLISTO_INDEX = "/ifs/depot/assemblies/H.sapiens/hg19/index/kallisto/v0.42.1/gencode/v18/gencode.v18.annotation.gtf.fasta.idx"
    GTF = "BIN_DIRECTORY_TOBE_SUBSTITUTE/data/gencode.v18.annotation.gtf"
    GTF_LNCRNA = "BIN_DIRECTORY_TOBE_SUBSTITUTE/data/lncipedia.gtf"
	starDB = "/ifs/depot/assemblies/H.sapiens/hg19/index/star/2.4.1d/gencode/v18/overhang74"
	starDB_LNCRNA = "/ifs/depot/assemblies/H.sapiens/hg19/index/star/2.3.0e_r291/LNCipedia"
	starDB_adaptor = "/ifs/depot/assemblies/H.sapiens/hg19/index/star/2.4.1d/gencode/v18/overhang49"
  }

 mm10 {
    REF_SEQ = "/ifs/depot/assemblies/M.musculus/mm10/mm10.fasta"
    GTF = "BIN_DIRECTORY_TOBE_SUBSTITUTE/data/Mus_musculus.GRCm38.80_canonical_chromosomes.gtf"
    DEXSEQ_GTF = "BIN_DIRECTORY_TOBE_SUBSTITUTE/data/Mus_musculus.GRCm38.80_canonical_chromosomes.dexseq.gtf"
    geneNameConversion = "BIN_DIRECTORY_TOBE_SUBSTITUTE/data/mm10Ensembl80IDToGeneName.txt"
    BOWTIE_INDEX = "/ifs/depot/assemblies/M.musculus/mm10/index/bowtie/1.1.1/mm10_bowtie"
    BOWTIE2_INDEX = "/ifs/depot/assemblies/M.musculus/mm10/index/bowtie/2.2.4/mm10_bowtie2"
    chrSplits = "/ifs/depot/assemblies/M.musculus/mm10/chromosomes"
    TRANS_INDEX = "/ifs/depot/assemblies/M.musculus/mm10/index/bowtie/2.2.4/transcriptome/ensembl/v80/Mus_musculus.GRCm38.80_canonical_chromosomes"
    TRANS_INDEX_DEDUP = ""
    TRANS_FASTA_DEDUP = ""
    RIBOSOMAL_INTERVALS = "BIN_DIRECTORY_TOBE_SUBSTITUTE/data/ribosomal_mm10.interval_file"
    REF_FLAT = "BIN_DIRECTORY_TOBE_SUBSTITUTE/data/refFlat__mm10.txt.gz"
    KALLISTO_INDEX = ""
	starDB = "/ifs/depot/assemblies/M.musculus/mm10/index/star/2.4.1d/ensembl/v80/overhang74"
	starDB_adaptor= "/ifs/depot/assemblies/M.musculus/mm10/index/star/2.4.1d/ensembl/v80/overhang49"
  }

 mm9 {
    REF_SEQ = "/ifs/depot/assemblies/M.musculus/mm9/mm9.fasta"
    GTF = "BIN_DIRECTORY_TOBE_SUBSTITUTE/data/Mus_musculus.NCBIM37.67_ENSEMBL.gtf"
    DEXSEQ_GTF = "BIN_DIRECTORY_TOBE_SUBSTITUTE/data/Mus_musculus.NCBIM37.67_ENSEMBL.dexseq.gtf"
    geneNameConversion = "BIN_DIRECTORY_TOBE_SUBSTITUTE/data/mm9Ensembl67IDToGeneName.txt"
    BOWTIE_INDEX = "/ifs/depot/assemblies/M.musculus/mm9/index/bowtie/1.0.0/mm9_bowtie"
    BOWTIE2_INDEX = "/ifs/depot/assemblies/M.musculus/mm9/index/bowtie/2.1.0/mm9_bowtie2"
    chrSplits = "/ifs/depot/assemblies/M.musculus/mm9/chromosomes"
    TRANS_INDEX = "/ifs/depot/assemblies/M.musculus/mm9/index/bowtie/2.1.0/transcriptome/ensembl/vTBD/ensembl"
    TRANS_INDEX_DEDUP = ""
    TRANS_FASTA_DEDUP = ""
    RIBOSOMAL_INTERVALS = "BIN_DIRECTORY_TOBE_SUBSTITUTE/data/ribosomal_MM9_assemblies.interval_file"
    REF_FLAT = "BIN_DIRECTORY_TOBE_SUBSTITUTE/data/refFlat__mm9.txt.gz"
    KALLISTO_INDEX = ""
	starDB = "/ifs/depot/assemblies/M.musculus/mm9/index/star/2.4.1d/ensembl/v67/overhang74"
	starDB_adaptor = "/ifs/depot/assemblies/M.musculus/mm9/index/star/2.4.1d/ensembl/v67/overhang49"
  }

 hybrid {
    REF_SEQ = "/ifs/depot/assemblies/hybrids/H.sapiens_M.musculus/hg19_mm10/hg19_mm10.fasta"
    GTF = "BIN_DIRECTORY_TOBE_SUBSTITUTE/data/gencode.v18.annotation.gtf"
    DEXSEQ_GTF = "BIN_DIRECTORY_TOBE_SUBSTITUTE/data/gencode.v18.annotation_dexseq.gtf"
    geneNameConversion = "BIN_DIRECTORY_TOBE_SUBSTITUTE/data/gencode18IDToGeneName.txt"
    RIBOSOMAL_INTERVALS = "BIN_DIRECTORY_TOBE_SUBSTITUTE/data/ribosomal_hg19.interval_file"
    REF_FLAT = "BIN_DIRECTORY_TOBE_SUBSTITUTE/data/refFlat__hg19.txt.gz"
	starDB = "/ifs/depot/assemblies/hybrids/H.sapiens_M.musculus/hg19_mm10/index/star/2.4.1d/gencode/v18/overhang74"
	starDB_adaptor = "/ifs/depot/assemblies/hybrids/H.sapiens_M.musculus/hg19_mm10/index/star/2.4.1d/gencode/v18/overhang49"
  }

 zv9 {
    species = "zv9"
    REF_SEQ = ""
    GTF = "BIN_DIRECTORY_TOBE_SUBSTITUTE/data/zv9.gtf"
    starDB = ""
    chrSplits = ""
    geneNameConversion = "Bin/data/zv9EnsemblIDtoGeneName.txt"
  }

 dm3 {
    REF_SEQ = "/ifs/depot/assemblies/D.melanogaster/dm3/dm3.fasta"
    GTF = "BIN_DIRECTORY_TOBE_SUBSTITUTE/data/dm3.flybase_more150bp_CollapseGenes_20140925.gtf"
    DEXSEQ_GTF = ""
    geneNameConversion = ""
    BOWTIE_INDEX = "/ifs/depot/assemblies/D.melanogaster/dm3/index/bowtie/1.1.1/dm3_bowtie"
    BOWTIE2_INDEX = "/ifs/depot/assemblies/D.melanogaster/dm3/index/bowtie/2.2.4/dm3_bowtie2" 
    chrSplits = "/ifs/depot/assemblies/D.melanogaster/dm3/chromosomes"
    TRANS_INDEX = ""
    TRANS_INDEX_DEDUP = ""
    TRANS_FASTA_DEDUP = ""
    RIBOSOMAL_INTERVALS = ""
    REF_FLAT = ""
    KALLISTO_INDEX = ""
    starDB = "/ifs/depot/assemblies/D.melanogaster/dm3/index/star/2.4.1d/flybase/custom20140925/overhang74" 
    starDB_adaptor = "/ifs/depot/assemblies/D.melanogaster/dm3/index/star/2.4.1d/flybase/custom20140925/overhang49" 
  }
 }
























