#!/usr/bin/env python

import argparse
import cmo
import datetime
import fileinput
import glob
import logging
import magic
import os
import re
import shutil
import subprocess
import sys
import vcf
import atexit
from collections import defaultdict



MAF2VCF_LOCATION = '/opt/common/CentOS_6/vcf2maf/v1.5.4/maf2vcf.pl'
#MAF2VCF_LOCATION = '/opt/common/CentOS_6/vcf2maf/v1.6.2/maf2vcf.pl'
VT_LOCATION = '/home/charris/code/VCF_accuracy_evaluator/vt/vt'
TABIX_LOCATION = '/opt/common/CentOS_6/samtools/samtools-1.2/htslib-1.2.1/tabix'
BGZIP_LOCATION = '/opt/common/CentOS_6/samtools/samtools-1.2/htslib-1.2.1/bgzip'
SORTBED_LOCATION = '/opt/common/CentOS_6/bedtools/bedtools-2.22.0/bin/sortBed'
BEDTOOLS_LOCATION = '/opt/common/CentOS_6/bedtools/bedtools-2.22.0/bin/bedtools'
BCFTOOLS_LOCATION = '/opt/common/CentOS_6/bcftools/bcftools-1.2/bin/bcftools' # python vcf parser doesn't like the vcf output from becftools norm



DIRS_TO_CLEANUP = []
@atexit.register
def cleanup():
    for file in DIRS_TO_CLEANUP:
        if os.path.exists(file) and os.path.isdir(file):
            logger.debug("CLEANUP: Delete %s" % file)
            shutil.rmtree(file);
        else:
            logger.info("Wanted to cleanup %s, but doesn't exist or isn't a directory" % file)



def cleanup_files_later(file):
    global DIRS_TO_CLEANUP
    if file not in DIRS_TO_CLEANUP:
        DIRS_TO_CLEANUP.append(file)


def check_for_programs():

    for f in [MAF2VCF_LOCATION]:

        if not os.path.exists(f):
            logger.critical('Unable to find required tool')
            logger.critical('File name not found: %s'%(f))
            logger.critical('Check constants at the top of compareVariants.py exist')
            sys.exit(1)

#####################################################
#####################################################

#Want to repalce this with the cmo python wrapper
#Uses an older version maf2vcf and doesn't like hg19
 
def convert_maf_to_vcf(maf, vcf_dir, reference):

    reference = reference.encode('ascii','ignore')
    #cmd = ['perl', MAF2VCF_LOCATION, '--input-maf', maf, '--output-dir', vcf_dir, '--ref-fasta', ref_fasta]

    cmd = ['perl', MAF2VCF_LOCATION, '--input-maf', maf, '--output-dir', vcf_dir]
    logger.debug('Running:' + ' '.join(cmd))
    subprocess.call(cmd)


def compare_samples(truth, test):

    truth = [each.split('/')[1] for each in truth]
    test = [each.split('/')[1] for each in test]

    not_in_truth = set(test).difference(set(truth))
    if len(not_in_truth) > 0:
        logger.critical('T/N Samples in test data not found in truth set')
        logger.critical(list(not_in_truth))
        sys.exit(1)
    else:
        return test


def subset_data(vcf, bedfile):
    
    outfile = vcf.replace('.vcf','.subset.vcf')
    cmd = [BEDTOOLS_LOCATION, 'intersect', '-header', '-a', vcf, '-b', bedfile]
    logger.debug('Running' + ' '.join(cmd))
    try:
        rv = subprocess.check_call(cmd, stdout=open(outfile,'w'))
    except subprocess.CalledProcessError, e:
        logger.critical("Non-zero exit code from bedtools intersect! Bailing out.")
        sys.exit(1)


def read_vcf(vcf_file):

    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    vcf_dict = {}
    for record in vcf_reader:
        chrom = record.CHROM
        pos = record.POS
        key = str(chrom) + ":" + str(pos)
        if key not in vcf_dict:
            vcf_dict[key]=record
        else:
            if isinstance(vcf_dict[key], list):
                vcf_dict[key].append(record)
            else:
                list = [vcf_dict[key], record]
                vcf_dict[key]=list
    samples = {}
    for key, record in vcf_dict.items():
        for call in record.samples:
            samples[call.sample]=1
        break
    return (vcf_dict, samples.keys())



def main(ref_vcf, test_vcf, file_type, reference, bedfile, normalize, prefix, plot, log_level):

    global logger
    logger = configure_logging(log_level)
    
    details = open('%s_details.out'%(prefix), 'w')

    check_for_programs()
    if file_type=="MAF":
        convert_maf_to_vcf(ref_vcf, 'truth_vcfs', reference) #would prefer to call on cmo wrappper
        convert_maf_to_vcf(test_vcf, 'test_vcfs', reference) #would prefer to call on cmo wrappper
    else:
        os.mkdir("truth_vcfs")
        os.mkdir("test_vcfs")
        shutil.copy(test_vcf, "test_vcfs")
        shutil.copy(truth_vcf, "truth_vcfs")
    cleanup_files_later(os.path.abspath("truth_vcfs"))
    cleanup_files_later(os.path.abspath("test_vcfs"))
    test = glob.glob('test_vcfs/*.vcf')
    truth = glob.glob('truth_vcfs/*.vcf')
    vcfs = compare_samples(truth, test)

    
    #subsetting
    if bedfile != None:
        for v in vcfs:
            subset_data('truth_vcfs/%s'%(v), bedfile)
            subset_data('test_vcfs/%s'%(v), bedfile)

        test = glob.glob('test_vcfs/*subset.vcf')
        truth = glob.glob('truth_vcfs/*subset.vcf')
        vcfs = compare_samples(truth, test)


    #normalization
    if normalize:
        for v in vcfs:
            test = cmo.util.normalize_vcf('truth_vcfs/%s'%(v), reference)
            truth = cmo.util.normalize_vcf('test_vcfs/%s'%(v), reference)
        test = glob.glob('test_vcfs/*.normalized.vcf.gz')
        truth = glob.glob('truth_vcfs/*.normalized.vcf.gz')
        vcfs = compare_samples(truth, test)
    

    #reading in vcfs
    sample_statistics = defaultdict(dict)
    truth_totals = defaultdict(dict)
    test_totals = defaultdict(dict)
    
    for v in vcfs:
        logger.info("Reading in VCF")
        (truth_chrom_pos_dict, truth_samples) = read_vcf('truth_vcfs/%s'%(v))
        (test_chrom_pos_dict,  test_samples) = read_vcf('test_vcfs/%s'%(v))
        if set(truth_samples) != set(test_samples):
            logger.critical("Samples not matched in truth and test vcf!")
            logger.critical("Truth Samples: %s" % truth_samples)
            logger.critical("Test Samples: %s" % test_samples)
            sys.exit(1)

        logger.info("Doing comparison")
        my_date = datetime.datetime.now().strftime('%Y-%m-%d')
        my_time = datetime.datetime.now().strftime('%H:%M:%S')


        #Loops through each test call to get get the totals for SNPs and INDELs called
        for site, test_record in test_chrom_pos_dict.items():
            if test_record.is_snp:
                site_type="SNP"
            elif test_record.is_indel:
                site_type="INDEL"
            else:
                logger.critical("Site is neither SNP nor INDEL?")
                logger.critical(test_record)
                logger.critical("Bailing out")
                sys.exit(1)

            for sample in test_samples:
                if test_record.genotype(sample).gt_type !=0: #if not reference genotype
                    
                    #add something to the test sample totals:
                    if site_type not in test_totals[sample]:
                        test_totals[sample][site_type]=1
                    else:
                        test_totals[sample][site_type]+=1

            #missed in truth
            if site not in truth_chrom_pos_dict:
                key = "Novel_" + site_type
                for sample in test_samples:
                    if test_record.genotype(sample).gt_type == 0: # this was ref/ref for this sample and is not a missed call:
                        pass #do nothing
                    else:
                        line = '%s\t%s\t%s\t%s\n'%(site_type, site, key, sample)
                        details.write(line)
                        if key not in sample_statistics[sample]:
                            sample_statistics[sample][key]=1
                        else:
                            sample_statistics[sample][key]+=1

        #Loops through each truth call
        for site, truth_record in truth_chrom_pos_dict.items():
            if truth_record.is_snp:
                site_type="SNP"
            elif truth_record.is_indel:
                site_type="INDEL"
            else:
                logger.critical("Site is neither SNP nor INDEL?")
                logger.critical(truth_record)
                logger.critical("Bailing out")
                sys.exit(1)

            for sample in truth_samples:
                if truth_record.genotype(sample).gt_type !=0: #if not reference genotype
                    #add something to the truth sample totals:
                    if site_type not in truth_totals[sample]:
                        truth_totals[sample][site_type]=1
                    else:
                        truth_totals[sample][site_type]+=1

            #missed in truth
            if site not in test_chrom_pos_dict:
                key = "Missed_" + site_type
                for sample in truth_samples:
                    if truth_record.genotype(sample).gt_type == 0: # this was ref/ref for this sample and is not a missed call:
                        pass #do nothing
                    else:
                        line = '%s\t%s\t%s\t%s\n'%(site_type, site, key, sample)
                        details.write(line)

                        if key not in sample_statistics[sample]:
                            sample_statistics[sample][key]=1
                        else:
                            sample_statistics[sample][key]+=1
                        

            else:
                #FIXME this will break if there are multiple records per site, i.e. a list - cough up blood
                test_record = truth_chrom_pos_dict[site]
                if isinstance(test_record, list):
                    logger.critical("Test file has multiple lines for one position!")
                    logger.critical("Change code to handle this!")
                    sys.exit(1)
                for sample in truth_samples:
                    truth_genotype = re.split("[/|]", truth_record.genotype(sample).gt_bases)
                    test_genotype = re.split("[/|]", test_record.genotype(sample).gt_bases)
                    if set(truth_genotype) == set(test_genotype):
                        key = "Correct_%s_Genotype" % site_type
                        line = '%s\t%s\t%s\t%s\n'%(site_type, site, key, sample)
                        details.write(line)
                        if key not in sample_statistics[sample]:
                            sample_statistics[sample][key]=1
                        else:
                            sample_statistics[sample][key]+=1
                    elif test_record.genotype(sample).gt_type !=0:
                        key = "Incorrect_%s_Genotype" % site_type
                        line = '%s\t%s\t%s\t%s\n'%(site_type, site, key, sample)
                        details.write(line)
                        if key not in sample_statistics[sample]:
                            sample_statistics[sample][key]=1
                        else:
                            sample_statistics[sample][key]+=1

    details.close()
    keys = ["Missed_SNP", "Novel_SNP", "Correct_SNP_Genotype", "Incorrect_SNP_Genotype", "Missed_INDEL", "Novel_INDEL", "Correct_INDEL_Genotype", "Incorrect_INDEL_Genotype"]
    ofh = open('%s.out'%(prefix), "w")
    logger.info("Writing stats to file...")

    ofh.write("\t".join(["sample"] + keys + ['Total_Truth_SNPs', 'Total_Truth_INDELs', 'Total_Test_SNPs', 'Total_Test_INDELs', 'Date', 'Time', 'Project_Prefix']) + '\n')
        
    norm_dict = defaultdict(dict)
    tumor_dict = defaultdict(dict)
    for sample, stats in sample_statistics.items():
        line = [sample]
        for key in keys:
            if key not in stats:
                line.append("0")
            else:
                line.append(str(stats[key]))
        for type in ["SNP", "INDEL"]:
            if type in truth_totals[sample]:
                line.append(str(truth_totals[sample][type]))
            else:
                line.append("0")
        for type in ["SNP", "INDEL"]:
            if type in test_totals[sample]:
                line.append(str(test_totals[sample][type]))
            else:
                line.append("0")

        ofh.write('%s\t%s\t%s\t%s\n'%('\t'.join(line), my_date, my_time, prefix))
    ofh.close()
        
    cmd = ['%s/vPlot.R'%(os.path.dirname(os.path.realpath(__file__))), '%s.out'%(prefix), prefix]
    subprocess.call(cmd)



def configure_logging(log_level):

    logger = logging.getLogger("MAF/VCF Accuracy Eval")
    logger.setLevel(log_level)
    ch = logging.StreamHandler()
    ch.setLevel(log_level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Takes in reference and test files and evaluates the test file.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--truth-file', action='store', dest='ref_file', default=None, required=True, help='"Truthful" reference file.')
    parser.add_argument('--test-file', action='store', dest='test_file', default=None, required=True, help='File to be tested.')
    mutex_group = parser.add_mutually_exclusive_group(required=True)
    mutex_group.add_argument("--vcf", dest="file_type", action="store_const", const="VCF", help="Input files are VCF format")
    mutex_group.add_argument("--maf", dest="file_type", action="store_const", const="MAF", help="Input files are MAF format")
    parser.add_argument('--reference', choices=cmo.util.genomes.keys(), required=True)
    parser.add_argument('--bedfile', action='store', dest='bedfile', default=None, help='Optional bedfile to limit the regions of comparison')
    parser.add_argument('--normalize', action='store_true', dest='normalize', default=False, help="Normalize variants with VT?" )
    parser.add_argument('-d', '--debug', action='store_const', const=logging.DEBUG, dest='log_level', default=logging.INFO, help='Turn on debug output')
    parser.add_argument('-p', '--prefix', action='store', dest='prefix', default='comparison_output', help='Prefix for output file and output column')
    parser.add_argument('-P', '--plot', action='store_true', dest='plot', default=True, help='Whether to plot percentages')
    args=parser.parse_args()
    ref_fasta = cmo.util.genomes[args.reference]['fasta']
    main(args.ref_file, args.test_file, args.file_type, ref_fasta, args.bedfile, args.normalize, args.prefix, args.plot, args.log_level)

