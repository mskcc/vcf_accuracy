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

from collections import defaultdict


#1. Check for programs
#2. Convert MAF to VCF
#3. Normalize - Optional
#4. Subset regions - Optional
#For each Tumor Normal pair:
#    4.1 read in vcf
#    4.2 do stats
#    4.3 store results
#5. write results to file
#6. delete temp files


MAF2VCF_LOCATION = '/opt/common/CentOS_6/vcf2maf/v1.5.4/maf2vcf.pl'
VT_LOCATION = '/home/charris/code/VCF_accuracy_evaluator/vt/vt'
TABIX_LOCATION = '/opt/common/CentOS_6/samtools/samtools-1.2/htslib-1.2.1/tabix'
BGZIP_LOCATION = '/opt/common/CentOS_6/samtools/samtools-1.2/htslib-1.2.1/bgzip'
SORTBED_LOCATION = '/opt/common/CentOS_6/bedtools/bedtools-2.22.0/bin/sortBed'
BEDTOOLS_LOCATION = '/opt/common/CentOS_6/bedtools/bedtools-2.22.0/bin/bedtools'
BCFTOOLS_LOCATION = '/opt/common/CentOS_6/bcftools/bcftools-1.2/bin/bcftools' # python vcf parser doesn't like the vcf output from becftools norm


#Not currently used to delete files
FILES_TO_CLEANUP = []
def cleanup():
    for file in FILES_TO_CLEANUP:
        logger.debug("CLEANUP: Delete %s" % file)
        os.unlink(file);


def cleanup_files_later(files):
    global FILES_TO_CLEANUP
    for file in files:
        if file not in FILES_TO_CLEANUP:
            FILES_TO_CLEANUP.append(file)


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
        cleanup_files_later([outfile])
    except subprocess.CalledProcessError, e:
        logger.critical("Non-zero exit code from bedtools intersect! Bailing out.")
        sys.exit(1)


def sort_vcf(vcf):

    outfile = vcf.replace('.vcf', '.sorted.vcf')
    cmd = [SORTBED_LOCATION, '-i', vcf, '-header']
    logger.debug('sortBed command: %s'%(' '.join(cmd)))
    try:
        rv = subprocess.check_call(cmd, stdout=open(outfile,'w'))
        return outfile
    except subprocess.CalledProcessError, e:
        logger.critical("Non-zero exit code from sortBed! Bailing out.")
        sys.exit(1)

    
def bgzip(vcf):

    if re.search('.gz', vcf):
        return vcf
    outfile = '%s.gz'%(vcf)
    cmd = [BGZIP_LOCATION, '-c', vcf]
    logger.debug('BGZIP COMMAND: %s'%(' '.join(cmd)))
    subprocess.call(cmd, stdout=open(outfile, 'w'))
    return outfile


def tabix_file(vcf_file):

    ''' index a vcf file with tabix for random access'''
    with magic.Magic(flags=magic.MAGIC_MIME_TYPE) as m:
        if(m.id_filename(vcf_file).find('gz') == -1):
            logger.critical('VCF File needs to be bgzipped for tabix random access. tabix-0.26/bgzip should be compiled for use')
            sys.exit(1)
    cmd = [TABIX_LOCATION, '-p' , 'vcf', vcf_file]
    logger.debug('Tabix command: %s'%(' '.join(cmd)))
    try:
        rv = subprocess.check_call(cmd)
        cleanup_files_later([vcf_file + '.tbi'])
    except subprocess.CalledProcessError, e:
        logger.critical('Non-zero exit code from Tabix! Bailing out.')
        sys.exit(1)


def normalize_vcf(vcf_file, ref_fasta):
    sorted_vcf = sort_vcf(vcf_file)
    zipped_file = bgzip(sorted_vcf)
    tabix_file(zipped_file)
    cleanup_files_later([sorted_vcf, zipped_file])
    output_vcf = zipped_file.replace('.vcf', '.normalized.vcf')
    cmd = [VT_LOCATION, 'normalize', '-r', ref_fasta, zipped_file, '-o', output_vcf, '-q']
    logger.debug('VT Command: %s'%(' '.join(cmd)))
    #cmd = [BCFTOOLS_LOCATION, 'norm', '-m', '-', '-O', 'b', '-o', output_vcf, zipped_file] #Python vcf parser doesn't like bcftools norm output
    #logger.info('bcftools norm Command: %s'%(' '.join(cmd)))
    try:
        rv = subprocess.check_call(cmd)
        cleanup_files_later([output_vcf])
        return output_vcf
    except subprocess.CalledProcessError, e:
        logger.critical("Non-zero exit code from normalization! Bailing out.")
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



def main(ref_maf, test_maf, reference, outfile, bedfile, normalize, log_level):

    global logger
    logger = configure_logging(log_level)


    check_for_programs()

    convert_maf_to_vcf(ref_maf, 'truth_vcfs', reference) #would prefer to call on cmo wrappper
    convert_maf_to_vcf(test_maf, 'test_vcfs', reference) #would prefer to call on cmo wrappper

    test = glob.glob('test_vcfs/*.vcf')
    truth = glob.glob('truth_vcfs/*.vcf')
    vcfs = compare_samples(truth, test)

    if bedfile != None:
        for v in vcfs:
            subset_data('truth_vcfs/%s'%(v), bedfile)
            subset_data('test_vcfs/%s'%(v), bedfile)

        test = glob.glob('test_vcfs/*subset.vcf')
        truth = glob.glob('truth_vcfs/*subset.vcf')
        vcfs = compare_samples(truth, test)


    if normalize in ('T','t','TRUE','true','True'):
        for v in vcfs:
            test = normalize_vcf('truth_vcfs/%s'%(v), reference)
            truth = normalize_vcf('test_vcfs/%s'%(v), reference)
        test = glob.glob('test_vcfs/*.normalized.vcf.gz')
        truth = glob.glob('truth_vcfs/*.normalized.vcf.gz')
        vcfs = compare_samples(truth, test)

    sample_statistics = defaultdict(dict)
    totals = defaultdict(dict)
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
                if truth_record.genotype(sample).gt_type !=0:
                    #add something to the sample totals:
                    if site_type not in totals[sample]:
                        totals[sample][site_type]=1
                    else:
                        totals[sample][site_type]+=1
            if site not in test_chrom_pos_dict:
                for sample in truth_samples:
                    if truth_record.genotype(sample).gt_type == 0: # this was ref/ref for this sample and is not a missed call:
                        pass #do nothing
                    else:
                        sample_statistics[sample]["Missed " + site_type]+=1
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
                        key = "Correct %s Genotype" % site_type
                        if key not in sample_statistics[sample]:
                            sample_statistics[sample][key]=1
                        else:
                            sample_statistics[sample][key]+=1
                    elif test_record.genotype(sample).gt_type !=0:
                        key = "Incorrect %s Genotype" % site_type
                        if key not in sample_statistics[sample]:
                            sample_statistics[sample][key]=1
                        else:
                            sample_statistics[sample][key]+=1

    keys = ["Missed SNP", "Correct SNP Genotype", "Incorrect SNP Genotype", "Missed INDEL", "Correct INDEL Genotype", "Incorrect INDEL Genotype"]
    ofh = open(outfile, "w")
    logger.info("Writing stats to file...")

    ofh.write("\t".join(["sample"] + keys + ['Total SNPS', 'Total INDELs', 'Date', 'Time']) + '\n')
    for sample, stats in sample_statistics.items():
        line = [sample]
        for key in keys:
            if key not in stats:
                line.append("0")
            else:
                line.append(str(stats[key]))
        for type in ["SNP", "INDEL"]:
            if type in totals[sample]:
                line.append(str(totals[sample][type]))
            else:
                    line.append("0")
        ofh.write('%s\t%s\t%s\n'%('\t'.join(line), my_date, my_time))
    ofh.close()
    
    shutil.rmtree('test_vcfs')
    shutil.rmtree('truth_vcfs')
    



def configure_logging(log_level):

    logger = logging.getLogger("MAF/VCF Accuracy Eval")
    if log_level=="DEBUG":
        log_level=logging.DEBUG
    logger.setLevel(log_level)
    ch = logging.StreamHandler()
    ch.setLevel(log_level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Takes in reference and test MAFs and evaluates the test MAF.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--truth-maf', action='store', dest='ref_maf', default=None, required=True, help='"Truthful" reference MAF.')
    parser.add_argument('--test-maf', action='store', dest='test_maf', default=None, required=True, help='MAF to be tested.')
    parser.add_argument('--reference', choices=cmo.util.genomes.keys(), required=True)
    parser.add_argument('--outfile', action='store', dest='outfile', default='comparison_output.txt', help='Comparison output file.')
    parser.add_argument('--bedfile', action='store', dest='bedfile', default=None, help='Optional bedfile to limit the regions of comparison.')
    parser.add_argument('--normalize', action='store', dest='normalize', default='F')
    parser.add_argument('--log-level', action='store', dest='log_level', default=logging.INFO, help='INFO for basic, DEBUG for (very) detailed debug output.')
    args=parser.parse_args()
    ref_fasta = cmo.util.genomes[args.reference]['fasta']
    main(args.ref_maf, args.test_maf, ref_fasta, args.outfile, args.bedfile, args.normalize, args.log_level)

