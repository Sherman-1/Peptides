from bin.transmembrane import transmembrane 
from tqdm import tqdm
from glob import glob 
import random

path = "input/OPM/7k15.pdb"

try:
    transmembrane(path, None, 15, 0, 20, 70, 1, None, "input/iORFs.csv", verbose = True)
except Exception as e:
    print(e)