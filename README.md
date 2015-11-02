![VCFs are Great](http://i.imgur.com/E8gfNu0.gif "Aww Yeah Accuracy Time!")


Install some Libs
=================

* use THIS python for below commands, so we get argparse + all these libs:
* /opt/common/CentOS_6-dev/python/python-2.7.10/bin/python
* download/install cmo package 
* git clone https://github.com/mskcc/cmo.git
* python setup.py install --user
* download/install magic package
* https://pypi.python.org/packages/source/f/filemagic/filemagic-1.6.tar.gz
* tar -xvf and python setup.py install --user
* download/install pyvcf
* https://pypi.python.org/packages/source/P/PyVCF/PyVCF-0.6.7.tar.gz
* tar -xvf and python setup.py install --user


Expected formatting
===================

MAF input requirements
----------------------

* Requires standard TCGA MAF columns (https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification):
  * Hugo_Symbol, Entrez_Gene_ID, Center, NCBI_Build, Chromosome, Start_Position, End_Position,
    Strand, Variant_Classification, Variant_Type, Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2,
    dbSNP_RS, dbSNP_Val_Status, Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, Match_Norm_Seq_Allele1,
    Match_Norm_Seq_Allele2, Tumor_Validation_Allele1, Tumor_Validation_Allele2, Match_Norm_Validation_Allele1,
    Match_Norm_Validation_Allele2, Verification_Status, Validation_Status, Mutation_Status, Sequencing_Phase,
    Sequence_Source, Validation_Method, Score, BAM_File, Sequencer

* Requires MSKCC specific columns:
  * n_ref_count, n_alt_count, n_depth, t_alt_count, t_ref_count, t_depth

* If subsetting the bedfile chromosome labelling must match the variant file style chr1 vs. 1, MT vs. M etc.
* Pre-normalized datasets must be left aligned and processed using VT (http://genome.sph.umich.edu/wiki/Vt).
