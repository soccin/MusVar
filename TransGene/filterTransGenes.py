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

    transContigs=set(transContigs)

    for alignGrp in read_sam_mapsByGroup(sys.argv[1]):

        ids=set([ai.id for ai in alignGrp])
        if len(ids)>1:
            raise ValueError()

        contigs=set([ai.contig for ai in alignGrp])
        mateContigs=set([ai.mateContig for ai in alignGrp])

        if not (
                transContigs.intersection(contigs)!=set()
                or transContigs.intersection(mateContigs)!=set()
                ):
            for ai in alignGrp:
                print("\t".join(ai.data))
