#!/usr/bin/env python

import sys
from tempfile import mkdtemp
from shutil import rmtree
from subprocess import call
from os import getcwd, environ, makedirs
from os.path import isfile, exists
from getpass import getuser

# A class full of tests, for nose to find and run, with fixtures for setup and cleanup
class TestsForNose:

    # Nose runs this only once before any of the tests are started
    @classmethod
    def setup_class( self ):
        tool_dir = getcwd()
        tool = tool_dir + '/compare_variant_sets.py'

        # If /scratch is found in this system, use /scratch/$USER as the folder for temp files
        if exists( '/scratch' ):
            tmp_dir = '/scratch/' + getuser()
            if not exists( tmp_dir ):
                makedirs( tmp_dir )
            environ['TMPDIR'] = tmp_dir

        # Build commands to run the tool on datasets included in the repo
        self.outdir = mkdtemp()
        cmd_maf = [ sys.executable, tool,
            '--maf', '--normalize', '--reference', 'GRCh37',
            '--first-file', tool_dir + '/data/truth.maf',
            '--second-file', tool_dir + '/data/test.maf',
            '--first-prefix', 'truth_set',
            '--second-prefix', 'test_set',
            '--output-dir', self.outdir + '/maf_parity' ]
        cmd_vcf = [ sys.executable, tool,
            '--vcf', '--normalize', '--reference', 'GRCh37',
            '--first-file', tool_dir + '/data/truth.vcf',
            '--second-file', tool_dir + '/data/test.vcf',
            '--first-prefix', 'truth_set',
            '--second-prefix', 'test_set',
            '--output-dir', self.outdir + '/vcf_parity' ]

        # Run those commands as shell commands
        try:
            self.cmd_maf_ret = call( cmd_maf, shell = False )
        except OSError as e:
            print >>sys.stderr, "Execution failed:", e
        try:
            self.cmd_vcf_ret = call( cmd_vcf, shell = False )
        except OSError as e:
            print >>sys.stderr, "Execution failed:", e

    # Nose runs this only once after all the tests have completed running
    @classmethod
    def teardown_class( self ):
        rmtree( self.outdir )

    # Check return codes from running the tool
    def test_retcodes( self ):
        assert self.cmd_maf_ret == 0
        assert self.cmd_vcf_ret == 0

    # Check that expected output files exist
    def test_output_files( self ):
        assert isfile( self.outdir + '/maf_parity/comparison_output.out' )
        assert isfile( self.outdir + '/maf_parity/comparison_output_details.out' )
        assert isfile( self.outdir + '/vcf_parity/comparison_output.out' )
        assert isfile( self.outdir + '/vcf_parity/comparison_output_details.out' )
