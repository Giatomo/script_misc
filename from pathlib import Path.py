from pathlib import Path
import subprocess as sp
import pandas as pd
from joblib import Parallel, delayed
from parfor import parfor


MODELS_DIR = Path("/home/thomas/Bureau/Model/")
PROTEOME_DIR = Path("/run/media/thomas/DATA/BIOINFORMATICS/new/proteomes/2021-04-06_12-52-25/files")
PROTEOME_INFO_FILE = Path("/run/media/thomas/DATA/BIOINFORMATICS/new/proteomes/assembly_summary.txt")
OUT_DIR = Path("/run/media/thomas/DATA/BIOINFORMATICS/macsyfind_results")
df = pd.read_csv(PROTEOME_INFO_FILE,
                sep='\t',
                header = None,
                index_col = False,
                names = ["RefSeq",
                         "BioProject",
                         "BioSamble",
                         "Genome",
                         "Genome_type",
                         "ID1",
                         "ID2",
                         "Species",
                         "Strain",
                         "Isolate",
                         "Assembly_type",
                         "Assembly_status",
                         "Assembly_release",
                         "Assembly_completion",
                         "Release_date",
                         "Assembly",
                         "Submitter",
                         "GenBank",
                         "Assembly_RefSeq_comparison",
                         "FTP",
                         "Assembly_from"
                ])

def macsyfinder(proteome):
    refseq = "_".join(proteome.name.split("_")[0:2])
    folder_name = df.loc[df['RefSeq'] == refseq].Species.values[0].replace(" ", "_").replace("/", ".")
    (OUT_DIR/folder_name).mkdir(parents= True, exist_ok = True)
    cmd = f"macsyfinder --db-type ordered_replicon --sequence-db {proteome} --models-dir {MODELS_DIR} -o {OUT_DIR/folder_name} --models Zaps all --no-cut-ga --mute"
    sp.run(cmd, shell = True)
    return

results = Parallel(n_jobs=-1, backend="multiprocessing")(
             map(delayed(macsyfinder), PROTEOME_DIR.glob("*.faa")))