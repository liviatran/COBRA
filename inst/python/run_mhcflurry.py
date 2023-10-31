import os, sys, io
from silence_tensorflow import silence_tensorflow
silence_tensorflow()
sys.stdout = open(os.devnull, 'w')
from mhcflurry import Class1PresentationPredictor
#ignore Model.state_updates Tensorflow UserWarning deprecation messages
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
import argparse

def run_mhcflurry(ciwd_file, pep_file):

    predictor = Class1PresentationPredictor.load()

    with open(ciwd_file, mode = "r") as file:
        ciwd_alleles = file.read().split(" ")

    with open(pep_file, mode = "r") as file:
        peps = file.read().split()

    alleles = dict((allele, [allele]) for allele in ciwd_alleles)

    x = predictor.predict(
        peptides=peps,
        alleles = alleles,
        verbose = 0)

    c_num = os.path.basename(os.path.splitext(ciwd_file)[0])
    p_file = os.path.basename(os.path.splitext(pep_file)[0])

    out_name = os.path.dirname(pep_file) + "/" + c_num + "_" + p_file + "_flur.csv"

    x.to_csv(out_name, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--allele")
    parser.add_argument("-p", "--peptide")

    args = parser.parse_args()

    ciwd_alleles = args.allele
    pep_file = args.peptide

    run_mhcflurry(ciwd_alleles, pep_file)
