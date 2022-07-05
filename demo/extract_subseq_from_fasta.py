import sys


def parse_region_str(regionstr):
    """
    :param regionstr: chrom:start-end or chrom:start or chrom.
     0-based, half-open region: [start, end)
    :return: chrom, start, end
    """
    try:
        if regionstr is None:
            return None, None, None
        elif ":" in regionstr:
            regioncse = regionstr.strip().split(":")
            assert len(regioncse) == 2
            chrom, se = regioncse[0], regioncse[1]
            if "-" in se:
                sne = se.split("-")
                assert len(sne) == 2
                s, e = int(sne[0]), int(sne[1])
                return chrom, s, e
            else:
                return chrom, int(se), None
        else:
            return regionstr.strip(), None, None
    except Exception:
        raise ValueError("--region not set right!")


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


def main():
    fafile = sys.argv[1]
    region = sys.argv[2]

    dnacontigs = DNAReference(fafile).getcontigs()
    regiontuple = parse_region_str(region)
    chrom, start, end = regiontuple

    subseq = ""
    if chrom is not None:
        chromseq = dnacontigs[chrom]
        if start is None:
            subseq = chromseq
        elif end is None:
            subseq = chromseq[start:]
        else:
            subseq = chromseq[start:end]

    print(">{}".format("_".join([chrom, str(start/1000)+"k", str(end/1000)+"k"])))
    idx = 0
    while idx < len(subseq):
        print(subseq[idx:(idx+60)])
        idx += 60


if __name__ == '__main__':
    main()
