import subprocess
import pandas as pd

class Nextclade:
    '''Runs nextclade and stores the information of the run'''
    
    def __init__(self) -> None:
        self.OUTDIR = 'nextclade_runs/'
        self.REF_DIR = 'reference/'
        self.INPUT = 'input.fasta'
        self.META_TSV = self.OUTDIR + 'meta.tsv'
        self.QUERY = self.OUTDIR + 'input.aligned.fasta'
        self.REF = self.REF_DIR + 'reference.fasta'
        
    def run_nextclade(self) -> None: 
        command = ['nextclade', '--in-order', '--input-fasta', self.INPUT, '--input-dataset', self.REF_DIR,
                   '--output-dir', self.OUTDIR, '--output-tsv', self.META_TSV]

        subprocess.call(command)
        
    def make_meta_df(self) -> None:
        return pd.read_csv(self.META_TSV, sep='\t')
        
        
        
