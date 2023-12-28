import argparse
import os
import re
import io
import subprocess
import numpy as np
import re
from joblib import Parallel, delayed
import multiprocessing


parser = argparse.ArgumentParser(description="Produce baypass format from sorted bams with various filtering specifications")

parser.add_argument("-r", "--ref_freq", type = float,
		    help = "Float between 0.5 and 1.0 to define upper bounder filter of the global reference allele frequency.")

parser.add_argument("-q", "--base_quality", type = int, default = 40,
		    help = "Phred score filter for base quality (defaults to 40)")

parser.add_argument("-Q", "--alignment_quality", type = int, default = 40,
		    help = "Phred score filter for base alignment quality (defaults to 40)")

parser.add_argument("-d", "--minor_count", type = int, 
	            help = "Positive integer for the minimum number of times the minor allele must be seen.")

parser.add_argument("-L", "--min_pop", type = int, 
	            help = "Positive integer for the minimum number of alleles counted in every population (sites with any failing pop are excluded).")

parser.add_argument("-M", "--max_pop", type = int,
	            help = "Positive integer for the maximum number of alleles counted in every population (sites with any failing pop are excluded).")

parser.add_argument('-l', '--positions_file', type=str, 
		    help='Optional file of reference position to pass to samtools mpileup.')

parser.add_argument('-p', '--prefix', type=str, 
		    help='prefix for naming output files (meta data and baypass)')

parser.add_argument("-n", "--num_cores", type = int, metavar = "int cores", default=1, 
		    help = "Number of cores to use. If multiple cores are not available, do not use flag, or set int cores to zero")

parser.add_argument("reference", type=str, help="The indexed reference genome each bam was aligned to")

parser.add_argument("bam_list", type=str, help="file containing path to bam file on each line")

args = parser.parse_args()

nucs = np.array(["A", "T", "G", "C", "N", "a", "t", "g", "c", "n"])

suffix = f"r{args.ref_freq}_d{args.minor_count}_L{args.min_pop}_M{args.max_pop}_q{args.base_quality}_Q{args.alignment_quality}.txt"
vep_str = args.prefix + "_vep_baypass_" + suffix

vep = open(vep_str, "w")

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

bam_count = file_len(args.bam_list)

def parse_seq(seq_str, qual_str, ref):
    
    seq_str = re.sub('\^.', '', seq_str.replace(".", ref).replace(",", ref).replace("$", "").upper())
    
    broken_seq = re.split('[+-]\d+', seq_str)
    
    seq_array =  np.array(list(broken_seq[0] + ''.join(
    [
    broken_seq[i][int(re.findall('\d+', seq_str)[i-1]):]
	    for i in list(range(1, len(broken_seq)))
    ]
    )
    ))

    seq_filtered = ''.join(seq_array)
    nucs = ["A", "T", "G", "C"]
    return [np.char.count(seq_filtered, nuc).item() for nuc in nucs]

def make_baypass(locus):
	bi_site = list(set(np.where(locus > 0)[1]))
	if len(bi_site) == 2:
		bi_pair = [" ".join((i[bi_site]).astype(str)) for i in locus]
		return " ".join(bi_pair)
	else:
		raise ValueError("This site is not bi-alleleic. Check your code!")


def parse_site(pileup):

	#get data parsed
	mp = pileup.strip().split("\t")
	chrom, pos, ref = mp[0:3] #site data
	pop_bam = mp[3:]
	
	#set up indices in groups of 3
	idx = list(range(0,len(pop_bam)-1, 3))

	#process read strings
	dp_array =  np.array([int(pop_bam[i]) for i in idx])

	dp_test = np.min(dp_array) >= args.min_pop and np.max(dp_array) <= args.max_pop

	if dp_test:
		read_list = [pop_bam[i+1] for i in idx]
		read_str = ' '.join(read_list).replace("*", "")
		ref_freq = (read_str.count(",") + read_str.count(".")) / np.sum(dp_array)
		freq_test = (1 - args.ref_freq) <= ref_freq <= args.ref_freq
	
		if freq_test:
			n_pop = list(range(len(read_list)))
			phred_list = [pop_bam[i+2] for i in idx]
			locus = np.array(Parallel(n_jobs = args.num_cores)(delayed(parse_seq)(read_list[i], phred_list[i], ref) for i in n_pop))
			locus_sums = locus.sum(axis = 1)
			locus_max_sums = locus_sums.max()
			locus_min_sums = locus_sums.min()
			dp_test2 = locus_min_sums >= args.min_pop and locus_max_sums <= args.max_pop

			if dp_test2:
			#locus = np.array([parse_seq(read_list[i], phred_list[i], ref) for i in n_pop])
			#print(locus)
				al_counts = np.sum(locus, axis = 0)
				bi_test = len(np.where(al_counts > 0)[0]) == 2
				if bi_test:
					non_zero = np.where(al_counts > 0)
					count_test = np.min(al_counts[non_zero]) >= args.minor_count
					if count_test:	
						l_baypass = make_baypass(locus)
						#print(len(l_baypass.split())
						if len(l_baypass.split()) == int(2*bam_count):
							alleles = nucs[list(set(np.where(locus > 0)[1]))]
							alt = "".join(alleles[alleles != ref].astype(str))
							print(chrom, pos, pos, ref + "/" + alt, ref_freq, l_baypass, file = vep)
							
							#l_baypass = make_baypass(locus)
							#print(l_baypass, file = baypass)
							#print("dp", np.min(dp_array), np.max(dp_array))
							#print("qual", locus_min_sums, locus_max_sums)
							#print("site", l_baypass)
							#print("site", len(l_baypass.split()))
	
							if alt == ref:
								raise ValueError("Invariant site! Abort!")



#construct mpileup command with and without positions argument
if args.positions_file:
	cmd = f"samtools mpileup -f {args.reference} -l {args.position_file} -b {args.bam_list} -q {args.base_quality} -Q {args.alignment_quality}".split()
else:
    cmd = f"samtools mpileup -f {args.reference} -b {args.bam_list} -q {args.base_quality} -Q {args.alignment_quality}".split()

proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
for line in io.TextIOWrapper(proc.stdout, encoding="utf-8"):
	parse_site(line)
	

