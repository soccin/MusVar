getFastqFiles<-function(fdir,read) {
   fs::dir_ls(fdir,recur=T,regex=paste0("_R",read,"_\\d+.fastq.gz")) %>% sort %>% list
}

require(tidyverse)

argv=commandArgs(trailing=T)

fdir=map(argv[1],read_tsv,col_names=F,show_col_types = FALSE) %>%
    bind_rows %>%
    select(sample=X2,fcid=X3,fdir=X4) %>%
    mutate(sample=gsub("^s_","",sample)) %>%
    mutate(sample=ifelse(
                    grepl("Mouse-Pooled-Normal",sample),
                    gsub("-\\d+$","",sample),
                    sample
                )
    ) %>%
    mutate(fcid=str_extract(fcid,"[^_]+_[^_]+_(.*)",group=T)) %>%
    rowwise %>%
    mutate(fastq_1=getFastqFiles(fdir,1),fastq_2=getFastqFiles(fdir,2)) %>%
    unnest(cols=c(fastq_1,fastq_2)) %>%
    select(-fdir)

if(!all(gsub("_R1_","_R2_",fdir$fastq_1)==fdir$fastq_2)) {
    cat("\n\nFATAL ERROR R1,R2 mismatch\n\n")
    rlang:abort("FATAL ERROR")
}

mfile=gsub("mapping","manifest",argv[1]) %>% gsub(".txt",".csv",.)
manifest=read_csv(mfile,show_col_types = FALSE) %>% mutate(status=ifelse(type=="Tumor",1,0)) %>% mutate(sample=gsub("^s_","",sample))

sarekInput=left_join(fdir,manifest) %>% 
    mutate(lane=cc(fcid,str_extract(fastq_1,"_(L\\d+)_",group=1))) %>%
    select(patient,sample,status,lane,matches("fastq"))

write_csv(sarekInput,"sarek_input_somatic.csv")

# #
# # Get normal pool inputs
# #
# poolMapFile=file.path(ROOT,"Proj_06000_MB_sample_mapping.txt")
# poolMap=read_tsv(poolMapFile,col_names=F) %>%
#     mutate(X2="s_Mouse-Pooled-Normal") %>%
#     select(sample=X2,fcid=X3,fdir=X4) %>%
#     mutate(sample=gsub("^s_","",sample)) %>%
#     mutate(fcid=str_extract(fcid,"[^_]+_[^_]+_(.*)",group=T)) %>%
#     rowwise %>%
#     mutate(fastq_1=getFastqFiles(fdir,1),fastq_2=getFastqFiles(fdir,2)) %>%
#     unnest(cols=c(fastq_1,fastq_2)) %>%
#     select(-fdir) %>%
#     mutate(status=0)

# if(!all(gsub("_R1_","_R2_",poolMap$fastq_1)==poolMap$fastq_2)) {
#     cat("\n\nFATAL ERROR R1,R2 mismatch\n\n")
#     rlang:abort("FATAL ERROR")
# }

# mi=list()
# for(fi in group_split(fdir,sample)) {
#     mi[[len(mi)+1]]=bind_rows(fi,poolMap) %>%
#                         mutate(patient=fi$sample[1]) %>%
#                         select(patient,everything()) %>%
#                         mutate(lane=cc(fcid,str_extract(fastq_1,"_(L\\d+)_",group=1))) %>%
#                         select(patient,sample,status,lane,matches("fastq"))

# }

# fs::dir_create(file.path("splits/somatic"))
# for(ii in seq(mi)) {
#     pid=mi[[ii]]$patient[1]
#     write_csv(mi[[ii]],file.path("splits/somatic",cc("pt",pid,"sarek_input_somatic.csv")))
# }