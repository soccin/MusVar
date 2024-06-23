import hashlib
from attrs import define, field

@define
class Genome:
    name = field()
    length = field()

knownGenomes=dict()
knownGenomes['9269f6897a317281255415cd844cd4b9']=Genome("mm10",66)

@define
class SAMHeader:
    tag = field()
    data = field()

    def __str__(self):
        return "\t".join([self.tag,*self.data])

    @classmethod
    def init_from_string(cls, samHeaderString):
        f=samHeaderString.strip().split("\t")
        return cls(f[0],f[1:])

def read_sam_header(fname):
    header=list()
    with open(fname,"r") as fp:
        for line in fp:
            if(line.startswith("@")):
                header.append(SAMHeader.init_from_string(line))
            else:
                return(header)

def compute_genome_hash(header):
    genomeDat=[h for h in header if h.tag=="@SQ"]
    genomeSigStr=";".join([str(g) for g in genomeDat[:10]])
    genomeHash=hashlib.md5(genomeSigStr.encode()).hexdigest()
    return genomeHash

def get_genome(header):
    genome=knownGenomes.get(compute_genome_hash(header),("Unknown",0))
    return genome

def get_transContigs(header):
    genome=get_genome(header)
    if genome.name=="Unknown":
        raise ValueError("Unkown Genome")
    genomeDat=[h for h in header if h.tag=="@SQ"]
    transContigs=[ct.data[0].split(":")[1] for ct in genomeDat[genome.length:]]
    return transContigs

def filterOut_trans_contigs(header,contigs):
    ctags=["SN:"+ct for ct in contigs]
    newHeader=[]
    for hi in header:
        if hi.tag=="@SQ" and hi.data[0] in ctags:
            pass
        else:
            newHeader.append(hi)
    return newHeader

def print_header(header):
    for hh in header:
        print(hh)

if __name__=="__main__":
    header=read_sam_header("s_2-3.sam")
    transContigs=get_transContigs(header)
    headerNew=filterOut_trans_contigs(header,transContigs)
    contigsForComment=["TG:"+ct for ct in transContigs]
    programComment=SAMHeader(
                    tag="@PG",
                    data=["ID:filterTransGenes.py",*contigsForComment]
                    )
    headerNew.append(programComment)
    print_header(headerNew)



