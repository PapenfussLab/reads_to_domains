import sys, os
from optparse import OptionParser
from subprocess import check_call
from collections import Counter
from mungo.fasta import FastaReader
from mungo.sequence import sixFrameTranslation

#UPROC stuff
UPROC_LOC = "/home/users/allstaff/tonkin-hill.g/get_domains/uproc-1.1.2/bin/uproc-dna"
MODEL_LOC = "/home/users/allstaff/tonkin-hill.g/get_domains/uproc-1.1.2/model/"
DB_LOC = "/home/users/allstaff/tonkin-hill.g/get_domains/uproc_analysis/uproc_domain_merged_to_2_without_segFilter"

#HMMER stuff
HMMER_SEARCH = "/home/users/allstaff/tonkin-hill.g/rask_based_block_finder2.0/third-party/hmmer-3.1b1/src/hmmsearch"
HMMER_HMM_DB = "/home/users/allstaff/tonkin-hill.g/get_domains/annotating_reads_to_domains/prot_version_ofDNA_HMMs.hmm"


def translate_6_frame(in_fasta, out_fasta):
    with open(out_fasta, 'w') as outfile:
        for h,s in FastaReader(in_fasta):
            for frame, seq in sixFrameTranslation(s).items():
                outfile.write(">"+h+" _frame_"+str(frame)+"\n")
                outfile.write(seq+"\n")
    return out_fasta

def run_uproc(outFile, read1,read2=None):
    if read2:
        uproc_cmd = (UPROC_LOC
            + " -p"
            + " -o " + outFile
            + " -P 3"
            + " -s"
            + " " + DB_LOC
            + " " + MODEL_LOC
            + " " + read1
            + " " + read2)
    else:
        uproc_cmd = (UPROC_LOC
            + " -p"
            + " -o " + outFile
            + " -P 3"
            + " -s"
            + " " + DB_LOC
            + " " + MODEL_LOC
            + " " + read1)

    print uproc_cmd
    check_call(uproc_cmd, shell=True)

    return outFile

def process_uproc_results(listFile, outFile, read1, read2=None):
    #First read in listfile from uproc that contains the read names we want to keep
    reads = set()
    with open(listFile, 'r') as listfile:
        for line in listfile:
            reads.add(line.split(",")[1])

    #now run through both fastq files outputing the found reads into fasta
    with open(outFile,'w') as outfile:
        with open(read1,'r') as fastqfile:
            for line in fastqfile:
                if line[0]=='@':
                    name = line[1:].strip().split()[0]
                    seq = fastqfile.next()
                    seq = seq.strip()
                    if name in reads:
                        outfile.write(">"+line[1:])
                        outfile.write(seq+"\n")
                        
        if read2:               
            with open(read2,'r') as fastqfile:
                for line in fastqfile:
                    if line[0]=='@':
                        name = line[1:].strip().split()[0]
                        seq = fastqfile.next()
                        seq = seq.strip()
                        if name in reads:
                            outfile.write(">"+line[1:])
                            outfile.write(seq+"\n")
    return outFile


def allocate_w_hmmer(infasta, outFile, evalue):
    hmm_cmd = (HMMER_SEARCH
        + " --tblout " + outFile
        + " -E " + str(evalue)
        + " " + HMMER_HMM_DB
        + " " + infasta
        + " > /dev/null")

    print hmm_cmd
    check_call(hmm_cmd, shell=True)

    return outFile

def process_hmmer_results(hmmoutput, outFile):
    #first get the maximum scoreing alignment for each read
    reads = {}

    with open(hmmoutput, 'r') as hmmfile:
        for line in hmmfile:
            if line[0]=="#":
                continue
            else:
                tokens = line.split()
                name = tokens[0]
                score = float(tokens[5])
                if name in reads:
                    if reads[name][1]<score:
                        reads[name] = (tokens[2],score)
                else:
                    reads[name] = (tokens[2],score)

    #now count reads for each domain query
    domain_hits = Counter()
    for read in reads:
        domain_hits[reads[read][0]] += 1

    #now sort and write to file
    domain_hit_list = domain_hits.items()
    domain_hit_list = sorted(domain_hit_list, key=lambda x: x[1]
        , reverse=True)

    with open(outFile,'w') as outfile:
        for domain in domain_hit_list:
            outfile.write(domain[0] + "," + str(domain[1]) + "\n")

    return outFile



def allocate_reads(read1, read2, outdir, dna, evalue):
    #first create output folder if it doesnt exist
    try:
        os.mkdir(outdir)
    except OSError, e:
        if e.errno != 17: #ignores error if folders has already been created.
            raise
        pass

    #now extract prefix from read name
    prefix = os.path.splitext(os.path.basename(read1))[0]

    uprocOut = run_uproc(outdir+prefix+"_uprocList.csv", read1, read2)

    uproc_reads = process_uproc_results(uprocOut, outdir+prefix+"_UprocReads.fa"
        , read1, read2)

    # uproc_reads = outdir+prefix+"_UprocReads.fa"

    if dna:
        hmm_out = allocate_w_hmmer(uproc_reads, outdir+prefix+"_nhmmOut_DNA.txt"
            , evalue)
    else:
        uproc_reads_6frame = translate_6_frame(uproc_reads
            , outdir+prefix+"_Uproc_6frame.fa")
        hmm_out = allocate_w_hmmer(uproc_reads_6frame
            , outdir+prefix+"_nhmmOut.txt", evalue)
    

    process_hmmer_results(hmm_out, outdir+prefix+"_DomainCount.csv")

    return



def main():
    parser = OptionParser()

    parser.add_option("-r", "--read1", dest="read1",
        help="the first fastq file")

    parser.add_option("-R", "--read2", dest="read2", default=None,
        help="the second fastq file")

    parser.add_option("", "--dna", dest="dna", action="store_true"
        , default=False
        , help="specifies that the HMMER models are in DNA form")

    parser.add_option("-E", "--evalue", dest="evalue", type=float, default=0.01
        , help="the E value to be fed to HMMER")

    parser.add_option("-o","--outdir", dest="outdir"
        , help="the output directory for files")

    (options, args) = parser.parse_args()

    allocate_reads(options.read1, options.read2, options.outdir, options.dna
        , options.evalue)


if __name__ == '__main__':
    main()
