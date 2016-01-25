# BIC RNAseq pipeline

This is a fork of the BIC RNAseq pipeline that has been coverted to use the bpipe ([github.com/ssadedin/bpipe](https://github.com/ssadedin/bpipe)) pipeline framework.

## Data directory

Do to file size limitations the data directory of this pipeline is not in the repository. To get the pipeline to work you need to link in the data directory on LUNA

```bash
ln -s /ifs/work/socci/Pipelines/CBE/rnaseq_bpipe/data $RNASEQ_BPIPE/data
```
	
Where `$RNASEQ_BPIPE` is the root of the repository.

