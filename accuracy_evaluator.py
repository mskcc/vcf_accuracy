#!/opt/common/CentOS_6-dev/python/python-2.7.10/bin/python

import sys, os, re, argparse,pprint, logging, magic, gzip
import vcf, cmo, subprocess, atexit
from collections import defaultdict

VT_LOCATION = "/home/charris/code/VCF_accuracy_evaluator/vt/vt"
TABIX_LOCATION = "/opt/common/CentOS_6/samtools/samtools-1.2/htslib-1.2.1/tabix"
BGZIP_LOCATION = "/opt/common/CentOS_6/samtools/samtools-1.2/htslib-1.2.1/bgzip"
SORTBED_LOCATION = "/opt/common/CentOS_6/bedtools/bedtools-2.22.0/bin/sortBed"
logger  = None
FILES_TO_CLEANUP = []
def configure_logging(log_level="INFO"):
    global logger
    logger = logging.getLogger("VCF_Comparator")
    if log_level=="DEBUG":
        log_level=logging.DEBUG
    logger.setLevel(log_level)
    ch = logging.StreamHandler()
    ch.setLevel(log_level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(ch)
    return logger

def tabix_file(vcf_file):
    """ index a vcf file with tabix for random access"""
    with magic.Magic(flags=magic.MAGIC_MIME_TYPE) as m:
        if(m.id_filename(vcf_file).find("gzip") == -1):
            logger.critical("VCF File needs to be bgzipped for tabix random access. tabix-0.26/bgzip should be compiled for use")
            sys.exit(1)
    cmd = [TABIX_LOCATION, "-p" , "vcf", vcf_file]
    logger.debug("Tabix command: %s" % " ".join(cmd))
    try:
        rv = subprocess.check_call(cmd)
        cleanup_files_later([vcf_file + ".tbi"])
    except subprocess.CalledProcessError, e:
        logger.critical("Non-zero exit code from Tabix! Bailing out.")
        sys.exit(1)

def bgzip(file):
    if re.search(".gz", file):
        return file
    cmd = [BGZIP_LOCATION, "-c", file, ">", file+".gz"]
    logger.debug("BGZIP COMMAND: %s" % " ".join(cmd))
    subprocess.call(" ".join(cmd), shell=True)
    return file + ".gz"

def sort_vcf(file):
    outfile = file.replace(".vcf", ".sorted.vcf")
    cmd = [SORTBED_LOCATION, "-i", file, "-header", ">", outfile]
    logger.debug("sortBed command: %s" % " ".join(cmd))
    try:
        rv = subprocess.check_call(" ".join(cmd), shell=True)
        return outfile
    except subprocess.CalledProcessError, e:
        logger.critical("Non-zero exit code from sortBed! Bailing out.")
        sys.exit(1)

@atexit.register
def cleanup():
    for file in FILES_TO_CLEANUP:
        logger.debug("CLEANUP: Delete %s" % file)
        os.unlink(file);


def cleanup_files_later(files):
    global FILES_TO_CLEANUP
    for file in files:
        if file not in FILES_TO_CLEANUP:
            FILES_TO_CLEANUP.append(file)

def normalize_vcf(vcf_file, ref_fasta):
    sorted_vcf = sort_vcf(vcf_file)
    zipped_file = bgzip(sorted_vcf)
    tabix_file(zipped_file)
    cleanup_files_later([sorted_vcf, zipped_file])
    output_vcf = zipped_file.replace(".vcf", ".normalized.vcf")
    cmd = [ VT_LOCATION, "normalize", "-r", ref_fasta, zipped_file, "-o", output_vcf]
    logger.debug("VT Command: %s" % " ".join(cmd))
    try:
        rv = subprocess.check_call(cmd)
        cleanup_files_later([output_vcf])
        return output_vcf
    except subprocess.CalledProcessError, e:
        logger.critical("Non-zero exit code from vt normalize! Bailing out.")
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

        

def main(truth_vcf, test_vcf, ref_fasta, output_file):
    logger.info("Running normalization...")
    normalized_truth_vcf=normalize_vcf(truth_vcf, ref_fasta)
    normalized_test_vcf=normalize_vcf(test_vcf, ref_fasta)
    logger.info("Reading in VCF")
    (truth_chrom_pos_dict, truth_samples) = read_vcf(normalized_truth_vcf)
    (test_chrom_pos_dict,  test_samples) = read_vcf(normalized_test_vcf)
    if set(truth_samples)!=set(test_samples):
        logger.critical("Samples not matched in truth and test vcf!")
        logger.critical("Truth Samples: %s" % truth_samples)
        logger.critical("Test Samples: %s" % test_samples)
        sys.exit(1)
    logger.info("Doing comparison")
    sample_statistics = defaultdict(dict)
    totals = defaultdict(dict)
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
                if truth_record.genotype(sample).gt_type ==0: # this was ref/ref for this sample and is not a missed call:
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
    ofh = open(output_file, "w")
    logger.info("Writing stats to file...")
    ofh.write("\t".join(["sample"] + keys + ["Total SNPS", "Total INDELs"]) + "\n")
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
        ofh.write("\t".join(line) + "\n")




if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Compare two files and spit out stats", add_help=True)
    parser.add_argument("--truth-vcf", action="store", required=True)
    parser.add_argument("--test-vcf", action="store", required=True)
    parser.add_argument("--reference", choices=cmo.util.genomes.keys(), required=True)
    parser.add_argument("--output", help="stats file for output", required=True)
    parser.add_argument("-v", action='store_true', default=False)
    args = parser.parse_args()
    if(args.v):
        configure_logging(logging.DEBUG)
    else:
        configure_logging(logging.INFO)
    test_vcf = args.test_vcf
    truth_vcf =args.truth_vcf
    ref_fasta = cmo.util.genomes[args.reference]['fasta']
    main(truth_vcf, test_vcf, ref_fasta, args.output)


