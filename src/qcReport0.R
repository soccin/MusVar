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

quibble2 <- function(x, q = c(0.25, 0.5, 0.75)) {
  tibble("{{ x }}" := quantile(x, q), "{{ x }}_q" := q)
}

get_qc_table<-function(tbl) {
    qvals=c(0,.25,.5,.75,1)
    tbl %>%
        select(all_of(names(which(map_vec(tbl,class)=="numeric")))) %>%
        gather(METRIC,V) %>%
        group_by(METRIC) %>%
        reframe(quibble2(V,qvals)) %>%
        mutate(V_q=cc("q",V_q)) %>%
        mutate(V_q=factor(V_q,levels=cc("q",qvals))) %>%
        spread(V_q,V) %>%
        filter(q_0<q_1 & q_0<q_0.25 & q_1>q_0.75)
}

require(tidyverse)
asm=fs::dir_ls("post/metrics",recur=T,regex=".as.txt$") %>%
    map(read_as_metrics) %>%
    bind_rows

qca=get_qc_table(filter(asm,!grepl("ZT",SAMPLE)))
as.metrics=qca %>% filter(grepl("PCT|BAL",METRIC)) %>% pull(METRIC)
qca1=qca %>% gather(Q,V,-METRIC) %>% filter(METRIC %in% as.metrics)

pg0=asm %>% filter(grepl("ZT",SAMPLE)) %>% select(all_of(as.metrics)) %>% gather(METRIC,VALUE) %>% ggplot(aes(METRIC,VALUE)) + theme_light() + geom_jitter() + coord_flip() +  theme(axis.title.y=element_blank(),axis.text.y=element_blank()) + geom_hline(aes(yintercept=V),data=qca1)
pg1=pg0 + ggforce::facet_wrap_paginate(~METRIC,scale="free",ncol=1,nrow=6,page=1)

mdm=fs::dir_ls("out/reports/markduplicates",recur=T,regex="metrics") %>%
    map(read_md_metrics) %>%
    bind_rows %>%
    type_convert

qcm=get_qc_table(filter(mdm,!grepl("ZT",SAMPLE)))
qcm1=qcm %>% gather(Q,V,-METRIC)

pg0=mdm %>% filter(grepl("ZT",SAMPLE)) %>% select(all_of(qcm$METRIC)) %>% gather(METRIC,VALUE) %>% ggplot(aes(METRIC,VALUE)) + theme_light() + geom_jitter() + coord_flip() +  theme(axis.title.y=element_blank(),axis.text.y=element_blank()) + geom_hline(aes(yintercept=V),data=qcm1)
pg1=pg0 + ggforce::facet_wrap_paginate(~METRIC,scale="free",ncol=1,nrow=6,page=1)




