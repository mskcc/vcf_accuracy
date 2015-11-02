#!/usr/bin/env python

# Convert DMP pipeline output to standard MAF output

####################################################################
####################################################################

import sys
import os
import csv
import pandas as pd
import argparse

def usage():
    print("%prog dmp_file out_prefix")
    sys.exit(-1)

def convert_dmp_to_maf(dmp_file, out_prefix):
    
#     dmp_columns = ['Sample', 'NormalUsed', 'Chrom', 'Start', 'Ref', 'Alt', 'VariantClass', 'Gene', 'Exon', 'Call_Confidence', 'Comments', 
#                    'TranscriptID', 'cDNAchange', 'AAchange', 'dbSNP_ID', 'Cosmic_ID', '1000G_MAF', 'FailureReason', 'CallMethod', 'COSMIC_site',
#                    'N_TotalDepth', 'N_RefCount', 'N_AltCount', 'N_AltFreq', 'T_TotalDepth', 'T_RefCount', 'T_AltCount', 'T_AltFreq', 'T_Ref+', 'T_Ref-', 'T_Alt+', 'T_Alt-', 
#                    'All_N_Aggregate_AlleleDepth', 'All_N_Median_AlleleFreq', 'T_freq/All_N_Freq', 'Occurence_in_Normals', 's_DL_scco_015_N']  
        
    maf_columns = ['Hugo_Symbol', 'Entrez_Gene_ID', 'Center', 'NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS','Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele2', 't_depth', 't_ref_count', 't_alt_count', 'n_depth', 'n_ref_count', 'n_alt_count']
    
#     Read and format DMP input
    dmp_reader = pd.read_csv(dmp_file, delimiter="\t")
    dmp_df = pd.DataFrame(dmp_reader)
    dmp_df = dmp_df.fillna('')
    dmp_df['Sample'] = dmp_df['Sample'].str.replace('-', '_')
    dmp_df['NormalUsed'] = dmp_df['NormalUsed'].str.replace('-', '_')
    dmp_df['NormalUsed'] = dmp_df['NormalUsed'].str.replace('FROZENPOOLEDNORMAL', 'Normal_Pooled_FROZEN_1')

#     Create output MAF
    out_maf = out_prefix + ".maf"
    out_file = open(out_maf, 'wb')
    out_writer = csv.writer(out_file, delimiter="\t")
    out_writer.writerow(maf_columns)
    
#     Convert using vcf2maf logic
    for idx, row in dmp_df.iterrows():
           
        ref = row['Ref']
        alt = row['Alt']
        ref_length = len(ref)
        alt_length = len(alt)
        start_position = row['Start']
#         Handle SNPs, DNPs, TNPs, ONPs
        if ref_length == alt_length:
            end_position = start_position + alt_length - 1
            if alt_length == 1:
                variant_type = "SNP"
            elif alt_length == 2:
                variant_type = "DNP"
            elif alt_length == 3:
                variant_type = "TNP"
            elif alt_length > 3:
                variant_type = "ONP"
#         Handle indels
        elif ref_length != alt_length:
            if ref_length < alt_length: # insertion
                ref = "-"
                alt = alt[ref_length:]
                end_position = start_position + 1
                variant_type = "INS"
            elif ref_length > alt_length: # deletion
                alt = "-"
                ref = ref[alt_length:]
                start_position += 1
                end_position = start_position + ref_length - 1
                variant_type = "DEL"
        
        hugo_symbol = row['Gene']
        entrez_gene_id = "0"
        center = "cmo.mskcc.org"
        ncbi_build = "GRCh37"
        chromosome = row['Chrom']
        strand = "+"
        variant_classification = row['VariantClass']
        reference_allele = ref
        tumor_seq_allele1 = ref
        tumor_seq_allele2 = alt
        if row['dbSNP_ID'] == "":
            dbsnp_rs = "novel"
        else:
            dbsnp_rs = row['dbSNP_ID'] 
        if row['Sample'].startswith('s_'): # prefix 's_' as per CMO convention
            tumor_sample_barcode = row['Sample']
        else:
            tumor_sample_barcode = 's_' + row['Sample']
        if row['NormalUsed'].startswith('s_'):
            matched_norm_sample_barcode = row['NormalUsed']
        else:
            matched_norm_sample_barcode = 's_' + row['NormalUsed'] # prefix 's_' as per CMO convention
        matched_norm_seq_allele1 = ref
        matched_norm_seq_allele2 = ref
        t_depth = row['T_TotalDepth']
        t_ref_count = row['T_RefCount']
        t_alt_count = row['T_AltCount']
        n_depth = row['N_TotalDepth']
        n_ref_count = row['N_RefCount']
        n_alt_count = row['N_AltCount']
        
        maf_variant = [hugo_symbol, entrez_gene_id, center, ncbi_build, chromosome, start_position, end_position, strand, variant_classification, variant_type, reference_allele, tumor_seq_allele1, tumor_seq_allele2, dbsnp_rs, tumor_sample_barcode, matched_norm_sample_barcode, matched_norm_seq_allele1, matched_norm_seq_allele2, t_depth, t_ref_count, t_alt_count, n_depth, n_ref_count, n_alt_count]
        out_writer.writerow(maf_variant)
            

argp = argparse.ArgumentParser( prog = "dmp2maf.py",
    description = "DMP to MAF conversion" )
argp.add_argument( "--dmp-file",        dest = "dmp",     type = file,
    help = "Input DMP file" )
argp.add_argument( "-o",        dest = "out",     type = str,
    help = "Output prefix" )


def _main( ):
    args = argp.parse_args( )
    convert_dmp_to_maf( args.dmp, args.out )
    
if __name__ == "__main__":
    _main( )