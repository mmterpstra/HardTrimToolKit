# HardTrimToolKit
Small toolkit to hardtrim sam reads based on (markIlluminaAdapters) trimming tag.

hardtrimadapters should trim on specified tag (XT if following the gatk best practices) set XO for reverting the read and remove the XT tag

### Pipeline integration notes:

Should work with the following dataflow.

> qname sorted bam | MarkIlluminaAdapters/MarkTrimming |  bwa mem | hardtrimadapters | FixMateInfo/MergeBamAlignment 

# Debug

```bash
(ml Pysam && ml git/2.38.1-GCCcore-12.2.0-nodocs  && export PYTHONPATH="$PWD/tests/install/lib/python3.12/site-packages:""$PYTHONPATH" && pip3 install -f -e . --prefix  tests/install/ ) 

#or ugly:
(ml Pysam && ml git/2.45.1-GCCcore-13.3.0  && export PYTHONPATH="$PWD/tests/install/lib/python3.12/site-packages:""$PYTHONPATH" && pip3 install -f -e . --prefix  tests/install/ && (ml SAMtools && samtools sort -n  path/to/aligned.bam ) | tests/install/bin/hardtrimadapters --input /dev/stdin --ubam-out True --ubam-in True --output /dev/stdout | (ml SAMtools && samtools view -h)  | tail ;)
#only running tests
(ml Pysam/0.22.1-GCC-13.3.0&& export PYTHONPATH="$PYTHONPATH:$PWD/src/" && pytest -o log_cli=true -o log_cli_level=DEBUG -vv) 
```