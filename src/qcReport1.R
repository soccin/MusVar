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

if(exists(".INCLUDED") && .INCLUDED ) halt(".#INCLUDE")
.INCLUDED=TRUE

manifest=read_csv("sarek_input_somatic.csv") %>%
    mutate(igoId=str_extract(fastq_1,"_IGO_([^/]*)/",group=1)) %>%
    mutate(machine=str_extract(fastq_1,"FASTQ.([^_]+)_",group=1)) %>%
    select(-fastq_1,-fastq_2,-lane) %>%
    mutate(machType=case_when(machine %in% machines$NovaSeq6000 ~ "NovaSeq6000", machine %in% machines$NovaSeqX ~ "NovaSeqX", T ~ "Unknown")) %>%
    mutate(projId=gsub("_\\d+","",igoId)) %>%
    mutate(sid=cc(patient,sample)) %>%
    distinct %>%
    mutate(status=ifelse(status==0,"Normal","Tumor"))

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

if(interactive()) x11()

boxWidth=.4

metrics0=c(
    "FOLD_ENRICHMENT", "MEAN_TARGET_COVERAGE",
    "PCT_CHIMERAS", "PCT_PF_READS_ALIGNED", "PCT_PF_READS_IMPROPER_PAIRS",
    "PCT_READS_ALIGNED_IN_PAIRS", "PCT_SOFTCLIP", "PCT_TARGET_BASES_100X",
    "PERCENT_DUPLICATION", "ZERO_CVG_TARGETS_PCT"
)

projNo=readLines("out/pipeline_info/cmd.sh.log") %>% grep("PWD:",.,value=T) %>% strsplit(.,"/") %>% unlist %>% grep("Proj",.,value=T) %>% gsub("Proj_","",.)

if(len(projNo)==0 || is.null(projNo) || projNo=="") {
    projNo="PrjoNo"
}

qcControl=read_csv("MusVar/assets/QC/qcControls_20250604_.csv") %>% filter(Metric %in% metrics0) %>% gather(Q,V,-Metric)
qcControl=qcControl %>% mutate(QQ=as.numeric(gsub("q_","",Q))) %>% mutate(OO=abs(log(QQ/(1-QQ)))/3+.5)

dq=bind_rows(list(asm,mdm,hsm)) %>% filter(Metric %in% metrics0) %>% distinct(sid,Metric,.keep_all=T)
pg=dq %>% ggplot(aes(status,Value,color=status)) + theme_light(16) + geom_hline(aes(yintercept=V,color=Q),data=qcControl,linewidth=qcControl$OO,alpha=.3) + geom_boxplot(outlier.shape=NA,width=boxWidth) + facet_wrap(~Metric,scale="free") + geom_jitter(width=boxWidth/2,size=3,alpha=.5) + scale_color_manual(values=c("darkblue","darkred","grey30","darkred","grey30","darkred","darkgreen")) + theme(panel.grid.minor = element_blank(),panel.grid.major.x=element_blank())
pdf(file=file.path("post/reports",cc("qcPlot_",projNo,".pdf")),width=1.5*14,height=1.5*8.5); print(pg); dev.off()

qcTbl=dq %>% select(patient,sample,status,Metric,Value) %>% spread(Metric,Value) %>% arrange(patient,status)
openxlsx::write.xlsx(qcTbl,file.path("post/reports",cc("qcTable_",projNo,".xlsx")))

quit()
rlang::abort("NOT IMPLEMENTED")

p0=bind_rows(list(asm,mdm,hsm)) %>%
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




