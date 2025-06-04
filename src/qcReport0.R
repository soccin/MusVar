read_as_metrics<-function(ff) {
    read_tsv(ff,comment="#",show_col_types=F,progress=F) %>%
        filter(CATEGORY=="PAIR") %>%
        select(SAMPLE,everything()) %>%
        select(-LIBRARY,-READ_GROUP,-PCT_HARDCLIP)
}

read_md_metrics<-function(ff) {
    read_tsv(ff,comment="#",show_col_types=F,progress=F,n_max=1, col_types=cols(.default="c")) %>%
        rename(SAMPLE=LIBRARY)
}

read_hs_metrics<-function(ff) {
    sid=basename(dirname(ff))
    read_tsv(ff,comment="#",show_col_types=F,progress=F,n_max=1,col_types=cols(.default="c")) %>%
        mutate(SAMPLE=sid)
}

quibble2 <- function(x, q = c(0.25, 0.5, 0.75)) {
  tibble("{{ x }}" := quantile(x, q), "{{ x }}_q" := q)
}

get_qc_table<-function(tbl) {
    qvals=c(.05/2,.25,.5,.75,(1-.05/2))
    tbl %>%
        select(all_of(names(which(map_vec(tbl,class)=="numeric")))) %>%
        gather(Metric,V) %>%
        group_by(Metric) %>%
        reframe(quibble2(V,qvals)) %>%
        mutate(V_q=cc("q",V_q)) %>%
        mutate(V_q=factor(V_q,levels=cc("q",qvals))) %>%
        spread(V_q,V)#        filter(q_0<q_1 & q_0<q_0.25 & q_1>q_0.75)
}

require(tidyverse)

#halt(".#INCLUDE")

as.metrics=c("SAMPLE", "PCT_PF_READS_ALIGNED", "MEAN_ALIGNED_READ_LENGTH",
    "PCT_READS_ALIGNED_IN_PAIRS", "PCT_PF_READS_IMPROPER_PAIRS",
    "STRAND_BALANCE", "PCT_CHIMERAS", "PCT_ADAPTER", "PCT_SOFTCLIP")

md.metrics=c("SAMPLE","PERCENT_DUPLICATION")

hs.metrics=c("SAMPLE",
"FOLD_ENRICHMENT",
"MEAN_TARGET_COVERAGE",
"PCT_TARGET_BASES_100X",
"PCT_USABLE_BASES_ON_TARGET",
"ZERO_CVG_TARGETS_PCT")

machines=list(
    NovaSeq6000=c("DIANA","MICHELLE","RUTH"),
    NovaSeqX=c("BONO","FAUCI2")
)

pmanifest=read_csv("projManifest_250603.csv")

manifest=read_csv("sarek_input_somatic.csv") %>%
    mutate(igoId=str_extract(fastq_1,"_IGO_([^/]*)/",group=1)) %>%
    mutate(machine=str_extract(fastq_1,"FASTQ.([^_]+)_",group=1)) %>%
    select(-fastq_1,-fastq_2,-lane,-status) %>%
    mutate(machType=case_when(machine %in% machines$NovaSeq6000 ~ "NovaSeq6000", machine %in% machines$NovaSeqX ~ "NovaSeqX", T ~ "Unknown")) %>%
    mutate(projId=gsub("_\\d+","",igoId)) %>%
    left_join(pmanifest) %>%
    mutate(type=ifelse(is.na(type),"invest",type)) %>%
    mutate(sid=cc(patient,sample)) %>%
    distinct

asm=fs::dir_ls("post/metrics",recur=T,regex=".as.txt$") %>%
    map(read_as_metrics,.progress=T) %>%
    bind_rows %>%
    select(all_of(as.metrics)) %>%
    gather(Metric,Value,-SAMPLE) %>%
    rename(sid=SAMPLE) %>%
    left_join(manifest)

mdm=fs::dir_ls("out/reports/markduplicates",recur=T,regex="metrics") %>%
    map(read_md_metrics,.progress=T) %>%
    bind_rows %>%
    type_convert %>%
    select(all_of(md.metrics)) %>%
    gather(Metric,Value,-SAMPLE) %>%
    rename(sample=SAMPLE) %>%
    left_join(manifest)

hsm=fs::dir_ls("post/metrics",recur=T,regex="_hs_metrics.txt") %>%
    map(read_hs_metrics,.progress=T) %>%
    bind_rows %>%
    type_convert %>%
    select(all_of(hs.metrics)) %>%
    gather(Metric,Value,-SAMPLE) %>%
    rename(sample=SAMPLE) %>%
    left_join(manifest)

x11()

boxWidth=.4
p0=bind_rows(list(asm,mdm,hsm)) %>%
    filter(type!="invest") %>%
    ggplot(aes(machType,Value)) + theme_light(20) + geom_boxplot(outlier.shape=NA,width=boxWidth) + facet_wrap(~Metric,scale="free") + geom_jitter(aes(color=type,shape=batch),width=boxWidth/2,size=3,alpha=.5) + ggsci::scale_color_jama()


mm=bind_rows(list(asm,mdm,hsm))

qcAll=mm %>%
    select(sid,Metric,Value) %>%
    spread(Metric,Value) %>%
    get_qc_table
write_csv(qcCtrl,cc("qcAll",DATE(),".csv"))
qca0=qcAll %>% gather(Q,V,-Metric)

mctrl=mm %>%
    filter(type!="invest")
qcCtrl=mctrl %>%
    filter(type=="CONTROLS" & batch=="A") %>%
    select(sid,Metric,Value) %>%
    spread(Metric,Value) %>%
    get_qc_table
write_csv(qcCtrl,cc("qcControls",DATE(),".csv"))

qca1=qcCtrl %>% gather(Q,V,-Metric)

pqc=mm %>% ggplot(aes(Metric,Value)) + theme_light() + coord_flip() + theme(axis.title.y=element_blank(),axis.text.y=element_blank()) + facet_wrap(~Metric,scale="free",ncol=2) + geom_jitter(alpha=.5) + scale_color_manual(values=c("darkred","grey30","grey50","grey30","darkred"))
pqc1=pqc + geom_hline(aes(yintercept=V,color=Q),data=(qca1))

# pg0=asm %>%
#     filter(grepl("ZT",SAMPLE)) %>%
#     select(all_of(as.metrics)) %>%
#     gather(METRIC,VALUE) %>%
#     ggplot(aes(METRIC,VALUE)) + theme_light() + geom_jitter() + coord_flip() +  theme(axis.title.y=element_blank(),axis.text.y=element_blank()) + geom_hline(aes(yintercept=V),data=qca1)
# pg1=pg0 + ggforce::facet_wrap_paginate(~METRIC,scale="free",ncol=1,nrow=6,page=1)

# mdm=fs::dir_ls("out/reports/markduplicates",recur=T,regex="metrics") %>%
#     map(read_md_metrics) %>%
#     bind_rows %>%
#     type_convert

# qcm=get_qc_table(filter(mdm,!grepl("ZT",SAMPLE)))
# qcm1=qcm %>%
#     gather(Q,V,-METRIC)

# pg0=mdm %>%
#     filter(grepl("ZT",SAMPLE)) %>%
#     select(all_of(qcm$METRIC)) %>%
#     gather(METRIC,VALUE) %>%
#     ggplot(aes(METRIC,VALUE)) + theme_light() + geom_jitter() + coord_flip() +  theme(axis.title.y=element_blank(),axis.text.y=element_blank()) + geom_hline(aes(yintercept=V),data=qcm1)
# pg1=pg0 + ggforce::facet_wrap_paginate(~METRIC,scale="free",ncol=1,nrow=6,page=1)




