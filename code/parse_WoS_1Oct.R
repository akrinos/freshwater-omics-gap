pacman::p_load(dplyr,ggplot2,ggridges,stringr,
               ggupset,tidyr,data.table,winch,patchwork,
               boot,coin,trend)

overall_marine=read.csv("data/marine_microbial_papers_20261001.txt",sep="\t") %>%
  dplyr::mutate(aquatic_system="marine")
overall_freshwater=read.csv("data/freshwater_microbial_papers_20261001.txt",sep="\t") %>%
  dplyr::mutate(aquatic_system="freshwater")
overall_numbers = overall_marine %>%
    dplyr::bind_rows(overall_freshwater)
ggplot(overall_numbers %>% dplyr::filter(Publication.Years>=1990)) + 
  geom_point(aes(x=Publication.Years,y=Record.Count,
                                         color=aquatic_system)) + 
  geom_line(aes(x=Publication.Years,y=Record.Count,color=aquatic_system))+
  theme_bw(base_size=12) + ylab("Number of papers") + 
  xlab("Publication Year")
wos_results_all = data.frame()
data_dir="data/WoS_final_search"
for (curr in list.files(data_dir)) {
  if (grepl("bib",curr)|grepl("xls",curr)) {
    next
  }
  print(nrow(read.csv(file.path(data_dir,curr),sep="\t",fileEncoding="latin1")))
  print(nrow(data.frame(fread(file.path(data_dir,curr),sep="\t"))))
  curr_fr=data.frame(fread(file.path(data_dir,curr),sep="\t")) %>%
      dplyr::select(AB,TI,PY,PT,DT,AU) %>%
      dplyr::mutate(AB=tolower(AB))
  wos_results_all=wos_results_all %>%
    dplyr::bind_rows(curr_fr)
  print(curr)
}

wos_results_all_articles = wos_results_all %>%
  dplyr::filter(grepl("Article",DT)) %>%
  dplyr::filter(!grepl("here we review progress",AB))
keyword_file = read.csv("keyword_11May2026.csv")
extract_keywords<-function(broad_category_in,specific_category_in,keyword_file_in=keyword_file){
  return(paste(as.character((keyword_file_in %>% 
                             dplyr::filter((broad_category==broad_category_in)&
                                           (specific_category==specific_category_in)))$keyword),collapse="|"))
}
wos_results_all_articles = wos_results_all_articles %>%
    dplyr::mutate(aquatic_system = case_when(str_detect(tolower(AB), 
                                                        extract_keywords("aquatic_system","extreme")) ~ "extreme",
                                             str_detect(tolower(AB), 
                                                        extract_keywords("aquatic_system","marine"))&
                                             !str_detect(tolower(AB), 
                                                         extract_keywords("aquatic_system","sea_urchin"))~ "marine",
                                             str_detect(tolower(AB), 
                                                        extract_keywords("aquatic_system","freshwater")) ~ "freshwater",
                                             str_detect(tolower(AB), 
                                                        extract_keywords("aquatic_system","brackish")) ~ "brackish"))%>%
  dplyr::mutate(cultivation_method = case_when(str_detect(tolower(AB), 
                                                          extract_keywords("cultivation_method","yes")) ~ "cultivation_method"))%>%
  dplyr::mutate(database_method = case_when(str_detect(tolower(AB), 
                                                       extract_keywords("database_method","yes")) ~ "database_method"))%>%
  dplyr::mutate(sequencing_method = case_when(str_detect(tolower(AB), 
                                                         extract_keywords("sequencing_method","yes")) ~ "sequencing_method"))%>%
  dplyr::mutate(sequencing_type = case_when(str_detect(tolower(AB), 
                                                       extract_keywords("sequencing_type","long-read")) ~ "long-read",
                                            str_detect(tolower(AB), 
                                                       extract_keywords("sequencing_type","metatranscriptomics")) ~ "metatranscriptomics",
                                            str_detect(tolower(AB), 
                                                       extract_keywords("sequencing_type","metaproteomics")) ~ "metaproteomics",
                                            str_detect(tolower(AB), 
                                                       extract_keywords("sequencing_type","metagenomics")) ~ "metagenomics",
                                            str_detect(tolower(AB), 
                                                       extract_keywords("sequencing_type","other meta-omics")) ~ "other meta-omics",
                                            str_detect(tolower(AB),  
                                                       extract_keywords("sequencing_type","[meta]barcoding")) ~ "[meta]barcoding",
                                           str_detect(tolower(AB), 
                                                      extract_keywords("sequencing_type","transcriptomics")) ~ "transcriptomics",
                                           str_detect(tolower(AB), 
                                                      extract_keywords("sequencing_type","proteomics")) ~ "proteomics",
                                           str_detect(tolower(AB), 
                                                      extract_keywords("sequencing_type","other DNA/genomics")) ~ "other DNA/genomics"))%>%
  dplyr::mutate(organism_group = case_when(str_detect(tolower(AB), 
                                                      extract_keywords("organism_group","euk_phytoplankton")) ~ "euk_phytoplankton",
                                           str_detect(tolower(AB), 
                                                      extract_keywords("organism_group","cyanobacteria")) ~ "cyanobacteria",
                                            str_detect(tolower(AB), 
                                                       extract_keywords("organism_group","fungus")) ~ "fungus",
                                            str_detect(tolower(AB), 
                                                       extract_keywords("organism_group","other_protist")) ~ "other_protist",
                                            str_detect(tolower(AB), 
                                                       extract_keywords("organism_group","bacteria"))&
                                            !str_detect(tolower(AB), 
                                                        extract_keywords("organism_group","e_coli"))~ "prokaryote"))%>% 
  dplyr::distinct()


passing_articles = wos_results_all_articles %>% 
  dplyr::filter(!is.na(sequencing_method)&!is.na(aquatic_system)&
                  !is.na(organism_group)) %>%
  dplyr::mutate(sequencing_type=case_when(is.na(sequencing_type)~"other sequencing",
                                          TRUE~sequencing_type))
write.csv(passing_articles %>% dplyr::distinct(),"output/passing_articles_13Oct2026.csv")
organized_by_cat = wos_results_all_articles %>%
    dplyr::filter(!is.na(sequencing_method)&!is.na(aquatic_system)&
                    #!is.na(sequencing_type)&
                    !is.na(organism_group)) %>%
    dplyr::group_by(aquatic_system,sequencing_type,PY) %>%
    dplyr::tally() %>%
    dplyr::mutate(sequencing_type=case_when(is.na(sequencing_type)~"other sequencing",
                                            TRUE~sequencing_type))

all_sequencing = wos_results_all_articles %>% 
  dplyr::filter(!is.na(sequencing_method)&!is.na(aquatic_system)&
                  !is.na(sequencing_type)&
                  !is.na(organism_group)) %>%
  dplyr::group_by(aquatic_system,PY) %>%
  dplyr::tally()

ggplot(all_sequencing %>% dplyr::filter(PY<2026)) + 
  geom_point(aes(x=PY,y=n,color=aquatic_system))+
  geom_line(aes(x=PY,y=n,color=aquatic_system,group=aquatic_system)) + 
  theme_bw(base_size=12) + ylab("number of papers") + xlab("") +
  scale_color_manual(values=c("blue","gray","#64E986","red"),
                     breaks=c("marine","brackish","freshwater","extreme"))

ggplot(organized_by_cat) + geom_point(aes(x=PY,y=n,color=sequencing_type,
                                          shape=aquatic_system))+ 
    geom_line(aes(x=PY,y=n,color=sequencing_type,
                  group=paste0(sequencing_type,aquatic_system)))

ggplot(organized_by_cat) + 
    geom_bar(aes(x=PY,y=n,fill=sequencing_type),stat="identity",
             position="fill")+ 
    facet_wrap(~aquatic_system,nrow=3,scales="free_y") + 
    scale_fill_manual(values=c("#e6194b", "#3cb44b", "#ffe119", "#4363d8", 
                               "#f58231", "#911eb4", "#42d4f4", "#f032e6",
                               "#bfef45", "#fabed4"))

ggplot(organized_by_cat %>% 
         dplyr::inner_join(overall_numbers %>%
                            dplyr::mutate(PY=as.numeric(as.character(Final.Publication.Year))))) + 
  geom_bar(aes(x=PY,y=n/Record.Count,fill=factor(sequencing_type,
                                                 levels=c("[meta]barcoding","other DNA/genomics",
                                                          "transcriptomics","proteomics",
                                                          "metagenomics","metatranscriptomics",
                                                          "metaproteomics","other meta-omics",
                                                          "long-read",
                                                          "other sequencing"))),stat="identity")+ 
  facet_wrap(~aquatic_system,nrow=3)  + 
  scale_fill_manual(values=c(
                        "#E69F00", "#D55E00", "#228833", "#CC79A7",
                        "navy", "#0072B2", "#56B4E9", "#44AA99",
                        "#AA3377", "black"
                    ),
                    breaks=c("[meta]barcoding","other DNA/genomics",
                             "transcriptomics","proteomics",
                             "metagenomics","metatranscriptomics",
                             "metaproteomics","other meta-omics",
                             "long-read",
                             "other sequencing"),
                    name="Sequencing Type") + 
  ylab("Proportion of papers") + xlab("Publication year") + 
  theme_bw(base_size=20)

ggsave("figures/freshwater_marine_prop_papers_20261002.pdf",width=8,height=5,units="in")

ggplot(organized_by_cat %>% 
         dplyr::inner_join(overall_numbers %>%
                             dplyr::mutate(PY=as.numeric(as.character(Final.Publication.Year))))%>%
         tidyr::pivot_wider(id_cols = c("sequencing_type","PY",
                                        "Record.Count"),
                            names_from="aquatic_system",values_from="n")%>%
         dplyr::mutate(freshwater=case_when(is.na(freshwater)~0,
                                   TRUE~freshwater),
                       marine=case_when(is.na(marine)~0,
                                        TRUE~marine))%>%
         tidyr::pivot_longer(cols=c("freshwater","marine"),names_to="aquatic_system",
                             values_to="n")) + 
  geom_boxplot(aes(y=n/Record.Count,x=aquatic_system,group=aquatic_system))+ 
  facet_wrap(~sequencing_type,nrow=3,scales="free_y")  +
  ylab("Proportion of papers") + xlab("Publication year") + 
  theme_bw(base_size=20)

prop_no_barcode = ggplot(organized_by_cat %>% 
         dplyr::inner_join(overall_numbers %>%
                             dplyr::mutate(PY=as.numeric(as.character(Final.Publication.Year)))) %>%
         dplyr::filter(!(sequencing_type%in%c("[meta]barcoding",
                                              "other DNA/genomics",
                                              "other sequencing")))) +
  geom_bar(aes(x=PY,y=n/Record.Count,fill=factor(sequencing_type,
                                    levels=c("[meta]barcoding","other DNA/genomics",
                                             "transcriptomics","proteomics",
                                             "metagenomics","metatranscriptomics",
                                             "metaproteomics","other meta-omics",
                                             "long-read","other sequencing"))),stat="identity")+
  facet_wrap(~aquatic_system,nrow=3)  + 
  scale_fill_manual(values=c(
                        "#E69F00", "#D55E00", "#228833", "#CC79A7",
                        "navy", "#0072B2", "#56B4E9", "#44AA99",
                        "#AA3377", "black"),
                    breaks=c("[meta]barcoding","other DNA/genomics",
                             "transcriptomics","proteomics",
                             "metagenomics","metatranscriptomics",
                             "metaproteomics","other meta-omics",
                             "long-read",
                             "other sequencing"),name="Sequencing Type")+ 
  ylab("Proportion of papers") + xlab("Publication year") + 
  theme_bw(base_size=20)

prop_no_barcode_seq = ggplot(organized_by_cat %>% 
                           dplyr::inner_join(organized_by_cat %>% dplyr::group_by(PY,aquatic_system) %>% 
                                               dplyr::summarize(Record.Count=sum(n))) %>%
                           dplyr::filter(aquatic_system%in%c("marine","freshwater"))%>%
                           dplyr::filter(!(sequencing_type%in%c("[meta]barcoding",
                                                                "other DNA/genomics",
                                                                "other sequencing")))) +
  geom_bar(aes(x=PY,y=n/Record.Count,fill=factor(sequencing_type,
                                                 levels=c("[meta]barcoding","other DNA/genomics",
                                                          "transcriptomics","proteomics",
                                                          "metagenomics","metatranscriptomics",
                                                          "metaproteomics","other meta-omics",
                                                          "long-read","other sequencing"))),stat="identity")+
  facet_wrap(~aquatic_system,nrow=1)  + 
  scale_fill_manual(values=c(
    "#E69F00", "#D55E00", "#228833", "#CC79A7",
    "navy", "#0072B2", "#56B4E9", "#44AA99",
    "#AA3377", "black"),
    breaks=c("[meta]barcoding","other DNA/genomics",
             "transcriptomics","proteomics",
             "metagenomics","metatranscriptomics",
             "metaproteomics","other meta-omics",
             "long-read",
             "other sequencing"),name="Sequencing Type")+ 
  ylab("Proportion of sequencing papers") + xlab("Publication year") + 
  theme_bw(base_size=20)

marine_fresh = ggplot(organized_by_cat %>% 
         dplyr::inner_join(overall_numbers %>%
                             dplyr::mutate(PY=as.numeric(as.character(Final.Publication.Year))))%>%
         tidyr::pivot_wider(id_cols = c("sequencing_type","PY",
                                        "Record.Count"),
                            names_from="aquatic_system",values_from="n")%>%
         dplyr::group_by(sequencing_type) %>%
         dplyr::summarise(freshwater=sum(freshwater,na.rm=T),
                          marine=sum(marine,na.rm=T)) %>%
         dplyr::filter(freshwater+marine > 5) %>%
         dplyr::mutate(freshwater_to_marine = marine/freshwater)) + 
  geom_point(aes(y=freshwater_to_marine,x=sequencing_type,group=sequencing_type))+
  geom_hline(aes(yintercept=1)) +
  ylab("Marine papers relative to freshwater") + xlab("") + 
  theme_bw(base_size=20)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
combo_fig=(prop_no_barcode / marine_fresh_pt / marine_fresh )+
  plot_annotation(tag_levels = "A") +plot_layout(heights = c(2, 1,1)) + 
  plot_layout(guides="collect")
ggsave(plot=combo_fig,filename="figures/marine_relative_to_fresh_20261002.pdf",width=10,height=12,units="in")

marine_fresh_pt = ggplot(organized_by_cat %>% 
                        dplyr::filter(sequencing_type!="[meta]barcoding") %>%
                        dplyr::inner_join(overall_numbers %>%
                                            dplyr::mutate(PY=as.numeric(as.character(Final.Publication.Year))))%>%
                        tidyr::pivot_wider(id_cols = c("sequencing_type","PY",
                                                       "Record.Count"),
                                           names_from="aquatic_system",values_from="n")%>%
                        dplyr::group_by(sequencing_type) %>%
                        dplyr::summarise(freshwater=sum(freshwater,na.rm=T),
                                         marine=sum(marine,na.rm=T)) %>%
                        dplyr::filter(freshwater+marine > 5) %>%
                        dplyr::mutate(freshwater_to_marine = marine/freshwater)) + 
  geom_point(aes(x=marine,y=freshwater,color=sequencing_type))+
  geom_abline(aes(intercept=0,slope=1)) + #scale_y_log10() + scale_x_log10() + 
  ylab("total freshwater papers") + xlab("total marine papers") + 
  theme_bw(base_size=20)+
  scale_color_manual(values=c(
    "#E69F00", "#D55E00", "#228833", "#CC79A7",
    "navy", "#0072B2", "#56B4E9", "#44AA99",
    "#AA3377", "black"),
    breaks=c("[meta]barcoding","other DNA/genomics",
             "transcriptomics","proteomics",
             "metagenomics","metatranscriptomics",
             "metaproteomics","other meta-omics",
             "long-read",
             "other sequencing"),name="Sequencing Type")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(organized_by_cat %>% 
         dplyr::inner_join(overall_numbers %>%
                             dplyr::mutate(PY=as.numeric(as.character(Final.Publication.Year))))%>%
         dplyr::group_by(aquatic_system,sequencing_type) %>%
         dplyr::summarise(n=sum(n))) + 
  geom_point(aes(y=n,x=aquatic_system,group=aquatic_system,color=sequencing_type))+ 
  #facet_wrap(~sequencing_type,nrow=3,scales="free_y")  +
  ylab("Proportion of papers") + xlab("Publication year") + 
  theme_bw(base_size=20)

paired_w_bootstrap <- function(x,y){
  # Calculate paired differences
  diffs <- x - y
  
  # Function to extract Wilcoxon V statistic
  wilcox_stat <- function(d) {
    # wilcox.test with paired=TRUE works on the original paired vectors, 
    # but here we pass centered differences or resampled values.
    # For a one-sample signed rank test on differences:
    res <- wilcox.test(d, mu = 0, exact = FALSE, alternative="greater")
    return(res$statistic) # Returns the V statistic
  }
  
  obs_stat <- wilcox_stat(diffs)
  
  # 2. Center the differences under the null hypothesis (median/mean shift = 0)
  centered_diffs <- diffs - median(diffs)
  
  # 3. Define the statistic function for tsboot (must accept data and indices/series)
  boot_stat <- function(d) {
    wilcox.test(d, mu = 0, exact = FALSE, alternative = "greater")$statistic
  }
  
  # 4. Run Moving Block Bootstrap using boot::tsboot
  # l = block length (setting to blocks of 3 years), R = number of replications
  b_reps <- 10000
  mbb_res <- tsboot(
    tseries = centered_diffs, 
    statistic = boot_stat, 
    R = b_reps, 
    l = 3, 
    sim = "fixed"
  )
  
  # 5. Compute the bootstrap p-value
  # how often bootstrap statistics are greater than or equal to the observed statistic
  boot_stats <- mbb_res$t
  p_val_greater <- mean(boot_stats >= obs_stat)
  
  print(paste("Observed V Statistic:", obs_stat))
  print(paste("One-sided Bootstrap p-value (x > y):", p_val_greater))
  
}

paired_perm_test <- function(x,y){
  # 2. Calculate the observed mean difference
  observed_diff <- mean(y - x)
  n <- length(x)
  
  # 3. Setup the permutation loop
  set.seed(42) # Ensures you get the same result every time
  n_permutations <- 10000
  perm_diffs <- numeric(n_permutations)
  
  # 4. Run the permutations (randomly flipping signs of the differences)
  observed_vector_diffs <- y - x
  
  for(i in 1:n_permutations) {
    # Randomly generate 1 or -1 for each row pair
    random_signs <- sample(c(-1, 1), n, replace = TRUE)
    # Apply signs to create a randomized mean difference
    perm_diffs[i] <- mean(observed_vector_diffs * random_signs)
  }
  
  # 5. Calculate the two-sided p-value
  p_value <- (sum(abs(perm_diffs) >= abs(observed_diff)) + 1) / (n_permutations + 1)
  
  # 6. Print the results
  cat("Observed Mean Shift :", observed_diff, "\n")
  cat("Permutation P-value :", p_value, "\n")
}

## paired Wilcoxen signed rank tests to see whether marine sci is consistently
## outperforming freshwater sci

# meta-omics comparison 
meta_omics = organized_by_cat %>% dplyr::filter(grepl("meta",sequencing_type)&
                                                  !grepl("barcoding",sequencing_type)) %>% 
  pivot_wider(id_cols = c("sequencing_type","PY"),names_from=aquatic_system,values_from=n) %>%
  replace_na(list(freshwater=0,marine=0)) %>%
  dplyr::group_by(PY) %>% dplyr::summarise(marine=sum(marine),freshwater=sum(freshwater))
wilcox.test(meta_omics$marine,meta_omics$freshwater,paired=TRUE,alternative="greater")
paired_w_bootstrap(meta_omics$marine,meta_omics$freshwater)
test_result <- paired_perm_test(meta_omics$marine,meta_omics$freshwater)

## do again, but this time as a proportion of all fw/marine papers
meta_omics = organized_by_cat %>% 
  dplyr::inner_join(overall_numbers %>%
                      dplyr::mutate(PY=as.numeric(as.character(Final.Publication.Year)))) %>% 
  dplyr::filter(grepl("meta",sequencing_type)&
                      !grepl("barcoding",sequencing_type)) %>%  
  dplyr::group_by(aquatic_system,Record.Count,PY) %>%
  dplyr::summarise(n=sum(n)) %>%
  dplyr::mutate(frac=n/Record.Count) %>%
  pivot_wider(id_cols = c("PY"),names_from=aquatic_system,values_from=frac) %>%
  replace_na(list(freshwater=0,marine=0))
wilcox.test(meta_omics$marine,meta_omics$freshwater,paired=TRUE,alternative="greater")
paired_w_bootstrap(meta_omics$marine,meta_omics$freshwater)
paired_perm_test(meta_omics$marine,meta_omics$freshwater)


## do again, but this time as a proportion of all fw/marine papers for txpts
meta_transcriptomics = organized_by_cat %>% 
  dplyr::inner_join(overall_numbers %>%
                      dplyr::mutate(PY=as.numeric(as.character(Final.Publication.Year)))) %>% 
  dplyr::filter(grepl("metatranscrip",sequencing_type)&
                  !grepl("barcoding",sequencing_type)) %>%  
  dplyr::group_by(aquatic_system,Record.Count,PY) %>%
  dplyr::summarise(n=sum(n)) %>%
  dplyr::mutate(frac=n/Record.Count) %>%
  pivot_wider(id_cols = c("PY"),names_from=aquatic_system,values_from=frac) %>%
  replace_na(list(freshwater=0,marine=0))
wilcox.test(meta_transcriptomics$marine,meta_transcriptomics$freshwater,paired=TRUE,alternative="greater")
paired_w_bootstrap(meta_transcriptomics$marine,meta_transcriptomics$freshwater)
paired_perm_test(meta_transcriptomics$marine,meta_transcriptomics$freshwater)

## do again, but this time as a proportion of all fw/marine papers for transcripts (not meta)
transcriptomics = organized_by_cat %>% 
  dplyr::inner_join(overall_numbers %>%
                      dplyr::mutate(PY=as.numeric(as.character(Final.Publication.Year)))) %>% 
  dplyr::filter((sequencing_type=="transcriptomics")&
                  !grepl("barcoding",sequencing_type)) %>%  
  dplyr::group_by(aquatic_system,Record.Count,PY) %>%
  dplyr::summarise(n=sum(n)) %>%
  dplyr::mutate(frac=n/Record.Count) %>%
  pivot_wider(id_cols = c("PY"),names_from=aquatic_system,values_from=frac) %>%
  replace_na(list(freshwater=0,marine=0))
wilcox.test(transcriptomics$marine,transcriptomics$freshwater,paired=TRUE,alternative="greater")
paired_w_bootstrap(transcriptomics$marine,transcriptomics$freshwater)
paired_perm_test(transcriptomics$marine,transcriptomics$freshwater)
transcriptomics %>% dplyr::summarise(mean_fw=mean(freshwater),sd_fw=sd(freshwater),
                                     mean_marine=mean(marine),sd_marine=sd(marine))

## do again, but this time as a proportion of all fw/marine papers for long-read
long_read = organized_by_cat %>% 
  dplyr::inner_join(overall_numbers %>%
                      dplyr::mutate(PY=as.numeric(as.character(Final.Publication.Year)))) %>% 
  dplyr::filter((sequencing_type=="long-read")&
                  !grepl("barcoding",sequencing_type)) %>%  
  dplyr::group_by(aquatic_system,Record.Count,PY) %>%
  dplyr::summarise(n=sum(n)) %>%
  dplyr::mutate(frac=n/Record.Count) %>%
  pivot_wider(id_cols = c("PY"),names_from=aquatic_system,values_from=frac) %>%
  replace_na(list(freshwater=0,marine=0))
wilcox.test(long_read$marine,long_read$freshwater,paired=TRUE,alternative="greater")
paired_w_bootstrap(long_read$marine,long_read$freshwater)
paired_perm_test(long_read$marine,long_read$freshwater)
long_read %>% dplyr::summarise(mean_fw=mean(freshwater),sd_fw=sd(freshwater),
                                     mean_marine=mean(marine),sd_marine=sd(marine))

## do again, but this time as a proportion of sequencing papers
meta_omics = organized_by_cat %>% 
  dplyr::inner_join(organized_by_cat %>% dplyr::group_by(PY,aquatic_system) %>%
                      dplyr::summarise(Record.Count=sum(n))) %>% 
  dplyr::filter(grepl("meta",sequencing_type)&
                  !grepl("barcoding",sequencing_type)) %>%  
  dplyr::group_by(aquatic_system,Record.Count,PY) %>%
  dplyr::summarise(n=sum(n)) %>%
  dplyr::mutate(frac=n/Record.Count) %>%
  pivot_wider(id_cols = c("PY"),names_from=aquatic_system,values_from=frac) %>%
  replace_na(list(freshwater=0,marine=0))
wilcox.test(meta_omics$marine,meta_omics$freshwater,paired=TRUE,alternative="greater")
paired_w_bootstrap(meta_omics$marine,meta_omics$freshwater)
paired_perm_test(meta_omics$marine,meta_omics$freshwater)

## now check if sequencing papers as a proportion of total is different.
all_sequencing = organized_by_cat %>% 
  dplyr::inner_join(overall_numbers %>%
              dplyr::mutate(PY=as.numeric(as.character(Final.Publication.Year)))) %>% 
  dplyr::group_by(PY,Record.Count,aquatic_system) %>% dplyr::summarise(n=sum(n)) %>%  
  dplyr::mutate(frac=n/Record.Count) %>%
  pivot_wider(id_cols = c("PY"),names_from=aquatic_system,values_from=frac) %>%
  replace_na(list(freshwater=0,marine=0))
wilcox.test(all_sequencing$marine,all_sequencing$freshwater,paired=TRUE,alternative="greater")
paired_w_bootstrap(all_sequencing$marine,all_sequencing$freshwater)
paired_perm_test(all_sequencing$marine,all_sequencing$freshwater)

## amplicon papers as a proportion of sequencing papers
barcode = organized_by_cat %>% 
  dplyr::inner_join(organized_by_cat %>% dplyr::group_by(PY,aquatic_system) %>%
                      dplyr::summarise(Record.Count=sum(n))) %>% 
  dplyr::filter(sequencing_type=="[meta]barcoding") %>% 
  dplyr::mutate(frac=n/Record.Count) %>%
  pivot_wider(id_cols = c("sequencing_type","PY"),names_from=aquatic_system,values_from=frac) %>%
  replace_na(list(freshwater=0,marine=0))
wilcox.test(barcode$marine,barcode$freshwater,paired=TRUE,alternative="less")
paired_w_bootstrap(barcode$marine,barcode$freshwater)

## test if sequencing is increasing or decreasing over time
frac_sequencing = organized_by_cat %>% 
  dplyr::inner_join(overall_numbers %>%
                      dplyr::mutate(PY=as.numeric(as.character(Final.Publication.Year)))) %>% 
  dplyr::group_by(PY,Record.Count,aquatic_system) %>% dplyr::summarise(n=sum(n)) %>%  
  dplyr::mutate(frac=n/Record.Count)
mk.test((frac_sequencing %>% dplyr::filter(aquatic_system=="marine"))$frac)
mk.test((frac_sequencing %>% dplyr::filter(aquatic_system=="freshwater"))$frac)

mk.test((frac_sequencing %>% dplyr::filter(aquatic_system=="marine") %>% dplyr::filter(PY>2017))$frac)
mk.test((frac_sequencing %>% dplyr::filter(aquatic_system=="freshwater") %>% dplyr::filter(PY>2013))$frac)
mk.test(frac_sequencing$Record.Count)

ggplot(organized_by_cat %>% 
         dplyr::full_join(overall_numbers %>%
                             dplyr::mutate(PY=as.numeric(as.character(Final.Publication.Year)))) %>%
         dplyr::group_by(PY,aquatic_system,Record.Count) %>% 
         dplyr::filter(aquatic_system%in%c("marine","freshwater")) %>%
         dplyr::summarise(n=sum(n)) %>%
         dplyr::mutate(complement=case_when(!is.na(n) ~ Record.Count-n,
                                            TRUE~Record.Count))) + 
  geom_bar(aes(x=PY,y=complement),stat="identity")+ 
  geom_bar(aes(x=PY,y=n,fill="sequencing"),stat="identity")+ 
  facet_wrap(~aquatic_system,nrow=3)  +
  ylab("Total number of papers") + xlab("Publication year") + 
  theme_bw(base_size=20)
ggsave("figures/freshwater_marine_sequencing_papers_20261002.pdf",width=8,height=5,units="in")


organized_by_cat_w_grp_total = organized_by_cat %>%
    dplyr::left_join(all_sequencing %>% dplyr::rename(c("all_n"="n"))) %>%
    dplyr::mutate(prop_in_group = n/all_n)

ggplot(organized_by_cat_w_grp_total %>% dplyr::filter(PY<2026)) + 
  geom_point(aes(x=PY,y=prop_in_group,color=aquatic_system))+ 
  geom_line(aes(x=PY,y=prop_in_group,color=aquatic_system,group=aquatic_system))+ 
  facet_wrap(~factor(sequencing_type,levels=c("DNA","RNA","barcoding",
                                              "long-read","meta-omics","proteomics",
                                              "other sequencing")),nrow=3,scales="free_y")+ 
  ylab("proportion of papers") + xlab("") + theme_bw(base_size=12)+
  scale_color_manual(values=c("blue","gray","#64E986","red"),
                    breaks=c("marine","brackish","freshwater","extreme"))

ggplot(organized_by_cat_w_grp_total %>% dplyr::filter(PY<2026)) + 
  geom_point(aes(x=PY,y=n,color=aquatic_system))+ 
  geom_line(aes(x=PY,y=n,color=aquatic_system,group=aquatic_system))+ 
  facet_wrap(~factor(sequencing_type,levels=c("DNA","RNA","barcoding",
                                              "long-read","meta-omics","metagenomics","proteomics",
                                              "other sequencing")),nrow=3,scales="free_y")+ 
  ylab("proportion of papers") + xlab("") + theme_bw(base_size=12)+
  scale_color_manual(values=c("blue","gray","#64E986","red"),
                     breaks=c("marine","brackish","freshwater","extreme"))

ggplot(organized_by_cat) + 
  geom_bar(aes(x=PY,y=n,fill=aquatic_system),stat="identity",
           position="fill")+ 
  facet_wrap(~sequencing_type,nrow=3,scales="free_y")+ 
  ylab("number of papers") + xlab("") +
  scale_fill_manual(values=c("blue","gray","#64E986","red"),
                     breaks=c("marine","brackish","freshwater","hypersaline"))

ggplot(organized_by_cat) + 
  geom_bar(aes(x=PY,y=n,fill=aquatic_system),stat="identity")+ 
  facet_wrap(~sequencing_type,nrow=3,scales="free_y")

ggplot(organized_by_cat) + 
  geom_bar(aes(x=PY,y=n,fill=sequencing_type),stat="identity")+ 
  facet_wrap(~aquatic_system,nrow=3,scales="free_y")

write.csv(wos_results_all_articles%>%
            dplyr::filter((aquatic_system=="freshwater")&
                          (sequencing_type=="meta-omics")),"wos_results_11May_fw_mo.csv")

write.csv(wos_results_all_articles%>%
            dplyr::filter((aquatic_system=="marine")&
                            (sequencing_type=="meta-omics")),"wos_results_11May_mar_mo.csv")
   
overall_marine=read.csv("data/marine_papers.txt",sep="\t") %>%
    dplyr::mutate(paper_type="marine")
overall_freshwater=read.csv("data/marine_papers.txt",sep="\t") %>%
    dplyr::mutate(paper_type="freshwater")
