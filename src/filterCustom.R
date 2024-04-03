x11 = function (...) grDevices::x11(...,type='cairo')

format_maf<-function(maf) {
    maf %>%
        select(1,5,6,9,11,13,16,HGVSp_Short,t_vaf,t_alt_count,n_vaf,n_depth,CALLERS,FILTERS) %>%
        mutate(Tumor_Sample_Barcode=fct_reorder(Tumor_Sample_Barcode,as.numeric(gsub(".*RR","",Tumor_Sample_Barcode)))) %>%
        arrange(Tumor_Sample_Barcode,desc(t_vaf),Start_Position)
}


require(tidyverse)

callerOrder=c("mutect2", "freebayes", "strelka")

mm=fs::dir_ls("int",rec=T,regex=".rda") %>% map(readRDS) %>% bind_rows %>% type_convert

QUAL=mm$vcf_qual
QUAL=as.numeric(ifelse(QUAL==".",-1,QUAL))
QUAL[QUAL<0]=NA
mm$QUAL=QUAL
ps=min(QUAL[QUAL>0],na.rm=T)

qualCut=15000

mm=mm %>%
    mutate(t_vaf=t_alt_count/t_depth,n_vaf=n_alt_count/n_depth) %>%
    mutate(FILTER=case_when(
                    CALLER=="freebayes" & QUAL>qualCut & t_vaf>5*n_vaf ~ "PASS",
                    T ~ FILTER
                )
    ) %>%
    mutate(CALLER=factor(CALLER,levels=callerOrder)) %>%
    arrange(ETAG,CALLER,Tumor_Sample_Barcode)

m1=mm %>%
    group_by(ETAG,Tumor_Sample_Barcode) %>%
    mutate(
        N.CALL=n(),N.PASS=sum(FILTER=="PASS"),
        CALLERS=paste(CALLER,collapse=";"),
        FILTERS=paste(FILTER,collapse="|")
        ) %>%
    ungroup

T_VAF_CUT=0.05

#
# 1 call list but relax mutect to all multiallele and clusters, strelka get rid of lowEVS, FreeBays Q>15000
#
mafCust=m1 %>%
    filter(QUAL>15000 | is.na(QUAL)) %>%
    filter(t_depth>=20 & t_alt_count>=8 & t_vaf>=T_VAF_CUT) %>%
    filter(!grepl("weak|normal|LowEVS",FILTERS)) %>%
    filter(!is.na(HGVSp_Short) & !grepl("=$",HGVSp_Short)) %>%
    distinct(ETAG,Tumor_Sample_Barcode,.keep_all=T)

maf1=m1 %>%
    filter(N.PASS>=1) %>%
    filter(t_vaf >= T_VAF_CUT & t_alt_count >= 8 & t_depth >= 20 & t_vaf > 5*n_vaf) %>%
    distinct(ETAG,Tumor_Sample_Barcode,.keep_all=T) %>%
    filter(!is.na(HGVSp_Short) & !grepl("=$",HGVSp_Short))

tblCust=mafCust %>% format_maf %>% filter(Hugo_Symbol %in% c("Rb1","Pten"))
tbl1=maf1 %>% format_maf %>% filter(Hugo_Symbol %in% c("Rb1","Pten"))

openxlsx::write.xlsx(
    list(HighSens=tbl1,CUSTOM=tblCust),
    "p15529_MusVar_V1_FilterCustom_2024-03-24.xlsx"
)

