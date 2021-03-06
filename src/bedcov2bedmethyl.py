import os
import argparse


alphabet = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}


def _alphabet(letter):
    if letter in alphabet.keys():
        return alphabet[letter]
    return 'N'


def complement_seq(dnaseq):
    rdnaseq = dnaseq[::-1]
    comseq = ''
    try:
        comseq = ''.join([_alphabet(x) for x in rdnaseq])
    except Exception:
        print('something wrong in the dna sequence.')
    return comseq


class DNAReference:
    def __init__(self, reffile):
        self._contignames = []
        self._contigs = {}  # contigname 2 contigseq
        with open(reffile, 'r') as rf:
            contigname = ''
            contigseq = ''
            for line in rf:
                if line.startswith('>'):
                    if contigname != '' and contigseq != '':
                        self._contigs[contigname] = contigseq
                        self._contignames.append(contigname)
                    contigname = line.strip()[1:].split(' ')[0]
                    contigseq = ''
                else:
                    # turn to upper case
                    contigseq += line.strip().upper()
            self._contigs[contigname] = contigseq
            self._contignames.append(contigname)

    def getcontigs(self):
        return self._contigs

    def getcontignames(self):
        return self._contignames


def convert_motif_seq(ori_seq):
    alphabets = {'A': ['A'], 'T': ['T'], 'C': ['C'], 'G': ['G'],
                 'R': ['A', 'G'], 'M': ['A', 'C'], 'S': ['C', 'G'],
                 'Y': ['C', 'T'], 'K': ['G', 'T'], 'W': ['A', 'T'],
                 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
                 'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
                 'N': ['A', 'C', 'G', 'T']}
    outbases = []
    for bbase in ori_seq:
        outbases.append(alphabets[bbase])

    def recursive_permute(bases_list):
        if len(bases_list) == 1:
            return bases_list[0]
        if len(bases_list) == 2:
            pseqs = []
            for fbase in bases_list[0]:
                for sbase in bases_list[1]:
                    pseqs.append(fbase + sbase)
            return pseqs
        else:
            pseqs = recursive_permute(bases_list[1:])
            pseq_list = [bases_list[0], pseqs]
            return recursive_permute(pseq_list)
    return recursive_permute(outbases)


def _coverage2bedmethyl(covfile, outfile, genomefile, motif="C", mloc=0):
    motifset = set(convert_motif_seq(motif))
    motifbase = motif[mloc]
    motiflen = len(motif)

    contigs = DNAReference(genomefile).getcontigs()
    contignames = set(contigs.keys())

    if outfile is None:
        wfile = covfile + ".bed"
    else:
        wfile = outfile
    wf = open(wfile, "w")
    with open(covfile, "r") as rf:
        for line in rf:
            words = line.strip().split("\t")
            chrom, start, end, rmet, cnt_m, cnt_um = words[0], int(words[1]), int(words[2]), \
                float(words[3]), int(words[4]), int(words[5])
            if chrom not in contignames:
                print("skip line cause <chrom not in genome>: {}".format(line.strip()))
                continue
            if start == end:
                start -= 1

            p_s = start - mloc
            p_e = start + motiflen - mloc
            if contigs[chrom][p_s:p_e] in motifset:
                strand = "+"
            else:
                p_s = start - (motiflen - mloc) + 1
                p_e = start + mloc + 1
                if complement_seq(contigs[chrom][p_s:p_e]) in motifset:
                    strand = "-"
                else:
                    print("skip line cause <not targeted motif in genome>: {}".format(line.strip()))
                    continue

            cnt_cov = cnt_m + cnt_um
            metperc = int(round(cnt_m / float(cnt_cov) * 100, 0))
            wf.write("\t".join([chrom, str(start), str(end), ".", str(cnt_cov), strand,
                                str(start), str(end), "0,0,0", str(cnt_cov), str(metperc)]) + "\n")
    wf.flush()
    wf.close()


def main():
    parser = argparse.ArgumentParser("converting bismark coverage file (from bismark2bedGraph cmd) to "
                                     "bedmethyl format")
    parser.add_argument("--cov", type=str, required=True,
                        help="coverage file generated by bismark2bedGraph command. "
                             "e.g., *zero.cov or *.cov")
    parser.add_argument("--genome", type=str, required=True,
                        help="genome reference, .fasta or .fa")
    parser.add_argument("--out", "-o", type=str, required=False, default=None,
                        help="output file path")
    parser.add_argument("--motif", type=str, required=False, default="C",
                        help="targeted motif, default C")
    parser.add_argument('--mloc_in_motif', type=int, required=False,
                        default=0,
                        help='0-based location of the methylation base in the motif, default 0')

    args = parser.parse_args()
    _coverage2bedmethyl(args.cov, args.out, args.genome, args.motif, args.mloc_in_motif)


if __name__ == '__main__':
    main()
