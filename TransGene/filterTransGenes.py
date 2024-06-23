from SAMTools import *
import sys

if __name__=="__main__":

    if len(sys.argv)!=2:
        print("\n\tusage: filterTransGenes.py INPUT.sam\n")
        sys.exit()

    header=read_sam_header(sys.argv[1])
    transContigs=get_transContigs(header)
    headerNew=filterOut_trans_contigs(header,transContigs)
    contigsForComment=["TG:"+ct for ct in transContigs]
    programComment=SAMHeader(
                    tag="@PG",
                    data=["ID:filterTransGenes.py",*contigsForComment]
                    )
    headerNew.append(programComment)
    print_header(headerNew)



