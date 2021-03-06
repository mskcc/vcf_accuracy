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

# Use cmo.util to locate the default maf2vcf script and bedtools binary
MAF2VCF_LOCATION = cmo.util.programs['vcf2maf']['default'] + "maf2vcf.pl"
BEDTOOLS_LOCATION = cmo.util.programs['bedtools']['default'] + "bedtools"

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
    for f in [MAF2VCF_LOCATION, BEDTOOLS_LOCATION]:
        if not os.path.exists(f):
            logger.critical('Unable to find required tool: %s'%(f))
            sys.exit(1)

#####################################################
#####################################################

def convert_maf_to_vcf(maf, vcf_dir, reference):

    reference = reference.encode('ascii','ignore')
    cmd = ['perl', MAF2VCF_LOCATION, '--input-maf', maf, '--output-dir', vcf_dir, '--ref-fasta', reference, '--per-tn-vcfs']
    logger.debug('Running:' + ' '.join(cmd))
    subprocess.call(cmd)

###################################################
###################################################

def common_samples(truth, test):

    truth = [each.split('/')[-1] for each in truth]
    test = [each.split('/')[-1] for each in test]

    # If there's a few samples of difference, then limit our analysis to just the subset
    not_in_truth = set(test).difference(set(truth))
    if len(not_in_truth) > 0:
        test_new = [x for x in test if x not in not_in_truth]
        return test_new
    else:
        return test

###################################################
###################################################

def subset_data(vcf, bedfile):

    outfile = vcf.replace('.vcf','.subset.vcf')
    cmd = [BEDTOOLS_LOCATION, 'intersect', '-header', '-a', vcf, '-b', bedfile]
    logger.debug('Running' + ' '.join(cmd))
    try:
        rv = subprocess.check_call(cmd, stdout=open(outfile,'w'))
    except subprocess.CalledProcessError, e:
        logger.critical("Non-zero exit code from bedtools intersect! Bailing out.")
        sys.exit(1)


###################################################
###################################################

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
                rec_list = [vcf_dict[key], record]
                vcf_dict[key]=rec_list
    samples = {}
    for key, record in vcf_dict.items():
        try:
            for call in record.samples:
                samples[call.sample]=1
            break
        except:
            continue
    if not len(samples.keys())>0:
        logger.critical("Unable to find sample names from vcf using pyvcf parser- is vcf formatted correctly?: %s" % vcf_file)
        sys.exit(1)

    return (vcf_dict, samples.keys())


###################################################
###################################################

def get_file_type():

    pass

###################################################
###################################################


def main(ref_file, test_file, first, second, file_type, outdir, reference, bedfile, normalize, prefix, plot, log_level):

    global logger
    logger = configure_logging(log_level)

    # Create the output folder if it doesn't already exist
    if not os.path.exists( outdir ):
        os.mkdir( outdir )
    details = open('%s/%s_details.out'%(outdir,prefix), 'w')

    check_for_programs()
    [test, truth] = [[], []]
    if file_type=="MAF":
        logger.info("Converting MAFs to VCFs")
        convert_maf_to_vcf(ref_file, 'truth_vcfs', reference)
        convert_maf_to_vcf(test_file, 'test_vcfs', reference)
    else:
        for dirname in [ 'truth_vcfs', 'test_vcfs' ]:
            if not os.path.exists( dirname ):
                os.mkdir( dirname )
        shutil.copy(test_file, "test_vcfs/sample.vcf")
        shutil.copy(ref_file, "truth_vcfs/sample.vcf")
    cleanup_files_later(os.path.abspath("truth_vcfs"))
    cleanup_files_later(os.path.abspath("test_vcfs"))
    test = glob.glob('test_vcfs/*.vcf')
    truth = glob.glob('truth_vcfs/*.vcf')
    vcfs = common_samples(truth, test)

    if bedfile is not None:
        logger.info("Subsetting samples to bedfile")
        for v in vcfs:
            subset_data('test_vcfs/%s'%(v), bedfile)
            subset_data('truth_vcfs/%s'%(v), bedfile)

        test = glob.glob('test_vcfs/*subset.vcf')
        truth = glob.glob('truth_vcfs/*subset.vcf')
        vcfs = common_samples(truth, test)

    if normalize:
        logger.info("Normalizing indels using bcftools/vt")
        for v in vcfs:
            test = cmo.util.normalize_vcf('truth_vcfs/%s'%(v), reference)
            truth = cmo.util.normalize_vcf('test_vcfs/%s'%(v), reference)
        test = glob.glob('test_vcfs/*.normalized.vcf.gz')
        truth = glob.glob('truth_vcfs/*.normalized.vcf.gz')
        vcfs = common_samples(truth, test)

    logger.info("Reading in VCFs")
    sample_statistics = defaultdict(dict)
    truth_totals = defaultdict(dict)
    test_totals = defaultdict(dict)

    for v in vcfs:
        (truth_chrom_pos_dict, truth_samples) = read_vcf('truth_vcfs/%s'%(v))
        (test_chrom_pos_dict, test_samples) = read_vcf('test_vcfs/%s'%(v))

        # Run a sanity check that should only apply to VCFs generated by maf2vcf
        if file_type=="MAF" and set(truth_samples) != set(test_samples):
            logger.critical("Samples not matched in truth and test vcf!")
            logger.critical("Truth Samples: %s" % truth_samples)
            logger.critical("Test Samples: %s" % test_samples)
            sys.exit(1)

        logger.info("Doing comparison")
        my_date = datetime.datetime.now().strftime('%Y-%m-%d')
        my_time = datetime.datetime.now().strftime('%H:%M:%S')

        logger.info("Looping through each test call to get the totals SNPs and INDELs called")
        for site, test_record in test_chrom_pos_dict.items():
            if isinstance(test_record, list):
                logger.critical("Skipping site with multiple variant calls")
                logger.critical(site)
                details.write("#Multiple calls: %s" % site)
                continue
            elif test_record.is_snp:
                site_type="SNP"
            elif test_record.is_indel:
                site_type="INDEL"
            else:
                logger.critical("Site is neither SNP nor INDEL?")
                logger.critical(test_record)
                logger.critical("Bailing out")
                sys.exit(1)

            if len(test_record.get_hom_alts()) > 0:
                site_depth = test_record.get_hom_alts()[0]['DP']
                site_cnt = test_record.get_hom_alts()[0]['AD'][1]
            else:
                site_depth = test_record.get_hets()[0]['DP']
                site_cnt = test_record.get_hets()[0]['AD'][1]

            # Only operate on the first sample listed, usually the tumor sample for a TN-pair
            sample = test_record.samples[0].sample
            if test_record.samples[0].gt_type !=0: #if not reference genotype
                if site_type not in test_totals[sample]:
                    test_totals[sample][site_type]=1
                else:
                    test_totals[sample][site_type]+=1

            # If this test call was missed in the truth set
            if site not in truth_chrom_pos_dict:
                key = "Novel_" + site_type
                if test_record.samples[0].gt_type == 0: # this was ref/ref for this sample and is not a missed call
                    pass #do nothing
                else:
                    line = '%s\t%s\t%s\t%s\t%s\t%s\n'%(site_type, site, key, sample, site_depth, site_cnt)
                    details.write(line)
                    if key not in sample_statistics[sample]:
                        sample_statistics[sample][key]=1
                    else:
                        sample_statistics[sample][key]+=1

        logger.info("Loops through each truth call")
        for site, truth_record in truth_chrom_pos_dict.items():
            if isinstance(truth_record, list):
                logger.critical("Skipping site with multiple variant calls")
                logger.critical(site)
                details.write("#Multiple calls: %s" % site)
                continue
            elif truth_record.is_snp:
                site_type="SNP"
            elif truth_record.is_indel:
                site_type="INDEL"
            else:
                logger.critical("Site is neither SNP nor INDEL?")
                logger.critical(truth_record)
                logger.critical("Bailing out")
                sys.exit(1)

            if len(truth_record.get_hom_alts()) > 0:
                site_depth = truth_record.get_hom_alts()[0]['DP']
                site_cnt = truth_record.get_hom_alts()[0]['AD'][1]
            else:
                site_depth = truth_record.get_hets()[0]['DP']
                site_cnt = truth_record.get_hets()[0]['AD'][1]

            # Only operate on the first sample listed, usually the tumor sample for a TN-pair
            sample = truth_record.samples[0].sample
            if truth_record.samples[0].gt_type !=0:
                if site_type not in truth_totals[sample]:
                    truth_totals[sample][site_type]=1
                else:
                    truth_totals[sample][site_type]+=1

            #missed in truth
            if site not in test_chrom_pos_dict:
                key = "Missed_" + site_type
                if truth_record.samples[0].gt_type == 0: # this was ref/ref for this sample and is not a missed call
                    pass #do nothing
                else:
                    line = '%s\t%s\t%s\t%s\t%s\t%s\n'%(site_type, site, key, sample, site_depth, site_cnt)
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
                sample = truth_record.samples[0].sample
                if truth_record.genotype(sample).gt_bases is not None and test_record.genotype(sample).gt_bases is not None:
                    truth_genotype = re.split('[/|]', truth_record.genotype(sample).gt_bases)
                    test_genotype = re.split('[/|]', test_record.genotype(sample).gt_bases)
                    if set(truth_genotype) == set(test_genotype):
                        key = "Correct_%s_Genotype" % site_type
                        line = '%s\t%s\t%s\t%s\t%s\t%s\n'%(site_type, site, key, sample, site_depth, site_cnt)
                        details.write(line)
                        if key not in sample_statistics[sample]:
                            sample_statistics[sample][key]=1
                        else:
                            sample_statistics[sample][key]+=1
                    elif test_record.genotype(sample).gt_type !=0:
                        key = "Incorrect_%s_Genotype" % site_type
                        line = '%s\t%s\t%s\t%s\t%s\t%s\n'%(site_type, site, key, sample, site_depth, site_cnt)
                        details.write(line)
                        if key not in sample_statistics[sample]:
                            sample_statistics[sample][key]=1
                        else:
                            sample_statistics[sample][key]+=1

    details.close()
    keys = ["Missed_SNP", "Novel_SNP", "Correct_SNP_Genotype", "Incorrect_SNP_Genotype", "Missed_INDEL", "Novel_INDEL", "Correct_INDEL_Genotype", "Incorrect_INDEL_Genotype"]
    ofh = open('%s/%s.out'%(outdir,prefix), "w")
    logger.info("Writing stats to file...")

    ofh.write("\t".join(["sample"] + keys + ['Total_Truth_SNPs', 'Total_Truth_INDELs', 'Total_Test_SNPs', 'Total_Test_INDELs', 'Union_SNPs', 'Union_INDELs', 'Date', 'Time', 'Project_Prefix']) + '\n')
    #print "\t".join(["sample"] + keys + ['Total_Truth_SNPs', 'Total_Truth_INDELs', 'Total_Test_SNPs', 'Total_Test_INDELs', 'Date', 'Time', 'Project_Prefix']) + '\n'
    norm_dict = defaultdict(dict)
    tumor_dict = defaultdict(dict)

    for sample, stats in sample_statistics.items():
        snps_union = 0
        indels_union = 0

        line = [sample]
        for key in keys:

            if re.search('SNP', key) and key in stats:
                snps_union = snps_union + stats[key]
                #print key, stats[key], snps_union
            if re.search('INDEL', key) and key in stats:
                indels_union = indels_union + stats[key]
                #print key, stats[key], indels_union
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
        #print sample, snps_union, indels_union

        ofh.write('%s\t%s\t%s\t%s\t%s\t%s\n'%('\t'.join(line), snps_union, indels_union, my_date, my_time, prefix))
    ofh.close()

    # Generate plots if asked for it
    if plot:
        cmd = ['%s/vPlot.R'%(os.path.dirname(os.path.realpath(__file__))), '%s/%s.out'%(outdir,prefix), '%s/%s_details.out'%(outdir,prefix), '%s/%s'%(outdir,prefix)]
        subprocess.call(cmd)

def configure_logging(log_level):

    logger = logging.getLogger("MAF/VCF Parity Evaluator")
    logger.setLevel(log_level)
    ch = logging.StreamHandler()
    ch.setLevel(log_level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger

# If this module is run as a script in CLI, the flow of control starts here
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Takes in a truth and test VCF/MAF to generate concordance metrics/plots', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    vcf_or_maf = parser.add_mutually_exclusive_group(required=True)
    vcf_or_maf.add_argument('--vcf', action='store_const', dest='file_type', const='VCF', help='Input files are VCF format, and only the first sample is compared')
    vcf_or_maf.add_argument('--maf', action='store_const', dest='file_type', const='MAF', help='Input files are MAF format')
    parser.add_argument('--first-file', action='store', dest='file1', required=True, help='"Truthful" reference VCF or MAF file')
    parser.add_argument('--second-file', action='store', dest='file2', required=True, help='VCF or MAF file to be tested for truthiness')
    parser.add_argument('--first-prefix', action='store', dest='truth', default='First', help='Short prefixed label describing the first file')
    parser.add_argument('--second-prefix', action='store', dest='test', default='Second', help='Short prefixed label describing the second file')
    parser.add_argument('--bedfile', action='store', dest='bedfile', help='Optional bedfile to limit the regions of comparison')
    parser.add_argument('--output-dir', action='store', dest='dir', default='parity', help='Writable folder to dump metrics/plots')
    parser.add_argument('--reference', choices=cmo.util.genomes.keys(), default='GRCh37', help='Reference build that these variants use')
    parser.add_argument('--prefix', action='store', dest='prefix', default='comparison_output', help='Prefix for output file and output column')
    parser.add_argument('--debug', action='store_const', const=logging.DEBUG, dest='log_level', default=logging.INFO, help='Set the logging level for debugging')
    parser.add_argument('--normalize', action='store_true', dest='normalize', default=False, help='Normalize calls with bcftools or vt')
    parser.add_argument('--plot', action='store_true', dest='plot', default=True, help='Whether to generate plots')
    args = parser.parse_args()
    ref_fasta = cmo.util.genomes[args.reference]['fasta']
    main(args.file1, args.file2, args.truth, args.test, args.file_type, args.dir, ref_fasta, args.bedfile, args.normalize, args.prefix, args.plot, args.log_level)
