from nd2reader import ND2Reader
from pprint import pprint
from pathlib import Path
import subprocess as sp

bf_convert_path = Path("/home/thomas/BioinformaticsTools/bftools/bfconvert")

raw_microscopy_path = Path("/mnt/Thomas/")
subdirs = (x for x in raw_microscopy_path.iterdir() if x.is_dir())

OUT_DIR = Path("/run/media/thomas/DATA/MICROSCOPY/")

for d in subdirs:
    print(d.name)
    out_dir = OUT_DIR/d.name
    out_dir.mkdir(parents = True, exist_ok = True)
    for nd2_file in d.glob("*.nd2"):
        print(nd2_file)
        out_file = out_dir/f"{nd2_file.name}_serie_%s.tiff"
        cmd = f"{bf_convert_path} {nd2_file} {out_file}"
        sp.run(cmd, shell = True)
