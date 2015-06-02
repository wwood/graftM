#!/usr/bin/env python

import subprocess
import os
import json
import shutil
import tempfile
import logging
import subprocess32

import graftm.getaxnseq
from Bio import SeqIO
from graftm.hmmer import Hmmer
from graftm.tree_cleaner import TreeCleaner

class Create:
    
    def __init__(self): 
        self.h=Hmmer(None, None)
        self.the_trash=[]
        
    def get_hmm_and_alignment(self, alignment, hmm, base):
        '''Return a HMM file and alignment of sequences to that HMM
        
        Returns
        -------
        * HMM file
        * Alignment of sequences to that HMM
        '''
        if not hmm:
            logging.debug("Building HMM from alignment")
            hmm=self.buildHmm(alignment, base)
            self.the_trash += [hmm]
            
        output_alignment = self.alignSequences(hmm, alignment, base)
        return hmm, output_alignment
    
    def buildHmm(self, alignment, base): 
        counter=0
        if os.path.isfile(base + ".hmm"):
            counter=0
            while os.path.isfile(base + ".hmm"):
                base=base+'_%s' % (str(counter))
                counter+=1
            hmm = base + ".hmm"
        else:   
            hmm = base + ".hmm" # Set a name for a hmm
        cmd = "hmmbuild %s %s >/dev/null" % (hmm, alignment) # Build the command to build the hmm
        logging.debug("Calling command: %s" % (cmd))
        subprocess.check_call(cmd, shell=True) # Call the command
        return hmm
    
    def pipeType(self, hmm):
        logging.debug("Setting pipeline type")
        type=[x.split() for x in open(hmm).readlines() if x.startswith('ALPH') or x.startswith('LENG')]
        for item in type:
            if item[0]=='ALPH':
                if item[1]=='amino':
                    ptype='aa'
                elif item[1]=='DNA' or 'RNA':
                    ptype='na'
                else:
                    raise Exception("Unfamiliar HMM type: %s" % (item[1]))
            elif item[0]=='LENG':
                leng=item[1]
            else:
                raise Exception("Programming Error: Misread HMM file")
        logging.debug("Set pipeline type as: %s " % ptype)
        logging.debug("Found alignment length as: %s" % leng)
        return ptype, leng
    
    def checkAlnLength(self, alignment):
        return len(list(SeqIO.parse(open(alignment, 'r'), 'fasta'))[0].seq)
        
    def alignSequences(self, hmm, sequences, base): 
        stockholm_alignment = base +".aln.sto" # Set an output path for the alignment
        fasta_alignment = base+".insertions.aln.fa" # Set an output path for the alignment
        corrected_fasta_alignment = base+".aln.fa" # Set an output path for the alignment
        cmd = "hmmalign --trim -o %s %s %s" % (stockholm_alignment, hmm, sequences) # Build the command to align the sequences
        logging.debug("Calling command %s" % (cmd))
        subprocess.check_call(cmd, shell=True) # Call the command
        cmd = "seqmagick convert --squeeze %s %s" % (stockholm_alignment, fasta_alignment)
        logging.debug("Calling command %s" % (cmd))
        subprocess.check_call(cmd, shell=True) # Call the command
        logging.debug("Correcting alignment")
        self.h.alignment_correcter([fasta_alignment], corrected_fasta_alignment)
        self.the_trash += [stockholm_alignment, corrected_fasta_alignment, fasta_alignment]
        logging.debug("Wrote corrected alignments to: %s" % (corrected_fasta_alignment))
        return corrected_fasta_alignment
    
    def buildTree(self, alignment, base, ptype): 
        log_file = base + ".tre.log"
        tre_file = base + ".tre"
        if ptype == 'na': # If it's a nucleotide sequence
            cmd = "FastTreeMP -quiet -gtr -nt -log %s -out %s %s" % (log_file, tre_file, alignment)
            logging.debug("Calling command: %s" % (cmd))
            subprocess.check_call(cmd, shell=True) # Call the command
        else: # Or if its an amino acid sequence
            cmd = "FastTreeMP -quiet -log %s -out %s %s" % (log_file, tre_file, alignment)
            logging.debug("Calling command: %s" % (cmd))
            subprocess.check_call(cmd, shell=True) # Call the command
            
        self.the_trash += [log_file, tre_file]
        return log_file, tre_file

    def callTaxitCreate(self, base, aln_file, tre, tre_log, tax, seq, prefix, no_reroot):
        if prefix:
            refpkg = prefix + ".refpkg"
        else:
            refpkg = base + ".refpkg"
            
        cmd = "taxit create -f %s -P %s -t %s -s %s -c -l  %s -T %s -i %s"\
            % (aln_file, refpkg, tre, tre_log, base, tax, seq)
            
        if no_reroot:
            cmd += ' --no-reroot'
            logging.debug("Calling command assuming pre-rerooting: %s" % cmd)
            subprocess32.check_call(cmd, shell=True)
        else:
            logging.debug("Calling command: %s" % cmd)
            logging.info("Attempting to run taxit create with rerooting capabilities")
            try:
                subprocess32.check_call(cmd, shell=True, timeout=20)
            except (subprocess32.TimeoutExpired, subprocess32.CalledProcessError):
                logging.error('''taxit create failed to run in a small amount of time suggesting that
rerooting was unsuccessful. Unfortunately this tree will need to be rerooted 
manually yourself using a tree editor such as ARB or FigTree.
Once you have a rerooted newick format tree, rerun graftm create
specifying the new tree with --rerooted_tree. The tree file to be rerooted is \'%s\'

''' % tre)
                exit(2)
        return refpkg
    
    def compile(self, base, refpkg, hmm, contents, prefix): 
        if prefix:
            gpkg = prefix + ".gpkg"
        else:
            gpkg = base + ".gpkg"
        if os.path.isdir(gpkg): 
            raise Exception("Detected gpkg with name %s" % (gpkg))
        os.mkdir(gpkg)
        os.rename(refpkg, os.path.join(gpkg, refpkg))
        shutil.copyfile(hmm, os.path.join(gpkg, os.path.basename(hmm)))
        json.dump(contents, open(os.path.join(gpkg, 'CONTENTS.json'), 'w'))

    def cleanup(self, the_trashcan):
        for every_piece_of_junk in the_trashcan:
            if os.path.isdir(every_piece_of_junk):
                shutil.rmtree(every_piece_of_junk)
            else:
                os.remove(every_piece_of_junk)
    
    def generate_tree_log_file(self, tree, alignment, output_tree_file_path,
                               output_log_file_path):
        '''Generate the FastTree log file given a tree and the alignment that
        made that tree
        
        Returns
        -------
        Nothing. The log file as parameter is written as the log file.
        '''
        cmd = "FastTree -quiet -nome -mllen -intree '%s' -log %s -out %s %s" %\
                                   (tree, output_log_file_path, 
                                    output_tree_file_path, alignment)
        logging.debug("Running log creation command %s" % cmd)
        subprocess.check_call(['bash','-c',cmd])
        
    def main(self, hmm, alignment, taxonomy, rerooted_tree, tree_log, prefix): 
        base=os.path.basename(alignment).split('.')[0]
            
        logging.info("Building gpkg for %s" % base)
        
        # align sequences to HMM (and potentially build hmm from alignment)
        hmm, output_alignment = self.get_hmm_and_alignment(alignment, hmm, base)
        
        # Create tree unless one was provided
        if not rerooted_tree:
            logging.debug("No tree provided")
            logging.info("Building tree")
            ptype,_ = self.pipeType(hmm)
            log_file, tre_file = self.buildTree(output_alignment, base, ptype)
            no_reroot = False
        else:
            logging.debug("Found pre-rerooted tree file %s" % rerooted_tree)
            tre_file=rerooted_tree
            no_reroot = True
            if tree_log:
                # User specified a log file, go with that
                logging.debug("Using user-specified log file %s" % tree_log)
                log_file = tree_log
            else:
                logging.info("Generating log file")
                log_file_tempfile = tempfile.NamedTemporaryFile(suffix='.tree_log', prefix='graftm') 
                log_file = log_file_tempfile.name
                input_tree_file = tre_file
                tre_file1_tempfile = tempfile.NamedTemporaryFile(suffix='.tree', prefix='graftm')
                tre_file1 = 'clean.tre'#tre_file1_tempfile.name
                # Make the newick file simple (ie. un-arb it) for fasttree
                TreeCleaner().clean_newick_file(input_tree_file, tre_file1)
                tre_file_tempfile = tempfile.NamedTemporaryFile(suffix='.tree', prefix='graftm')
                tre_file = tre_file_tempfile.name
                self.generate_tree_log_file(tre_file1, alignment,
                                            tre_file, log_file)
            
        # Create tax and seqinfo .csv files
        logging.info("Building seqinfo and taxonomy file")
        seq, tax = graftm.getaxnseq.main(base, taxonomy)
        self.the_trash += [seq, tax]
        
        # Create the reference package
        logging.info("Creating reference package")
        refpkg = self.callTaxitCreate(base, output_alignment, tre_file, 
                                      log_file, tax, seq, prefix, no_reroot)

        # Compile the gpkg
        logging.info("Compiling gpkg")
        contents = {"aln_hmm": hmm,
                    "search_hmm": [hmm],
                    "rfpkg": refpkg,
                    "TC":False}
        self.compile(base, refpkg, hmm, contents, prefix)

        logging.info("Cleaning up")
        self.cleanup(self.the_trash)
        
        logging.info("Finished\n")
