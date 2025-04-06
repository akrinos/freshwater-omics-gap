pacman::p_load(dplyr,ggplot2)

wos_results_fresh = data.frame()
data_dir="data/WoS_freshwater_5Apr2025"
for (curr in list.files(data_dir)) {
  wos_results_fresh=wos_results_fresh %>%
    dplyr::bind_rows(read.csv(file.path(data_dir,curr),sep="\t") %>%
                       dplyr::mutate(CY=as.numeric(as.character(CY)),
                                     NR=as.numeric(as.character(NR)),
                                     VL=as.numeric(as.character(VL)),
                                     PN=as.numeric(as.character(PN)),
                                     SU=as.numeric(as.character(SU)),
                                     BP=as.numeric(as.character(BP)),
                                     EP=as.numeric(as.character(EP))))
}
wos_results_fresh = wos_results_fresh %>%
    dplyr::mutate(Ecosystem="Freshwater")

wos_results_marine = data.frame()
data_dir="data/WoS_marine_6Apr2025"
for (curr in list.files(data_dir)) {
  wos_results_marine=wos_results_marine %>%
    dplyr::bind_rows(read.csv(file.path(data_dir,curr),sep="\t") %>%
                       dplyr::mutate(CY=as.numeric(as.character(CY)),
                                     NR=as.numeric(as.character(NR)),
                                     VL=as.numeric(as.character(VL)),
                                     PN=as.numeric(as.character(PN)),
                                     SU=as.numeric(as.character(SU)),
                                     BP=as.numeric(as.character(BP)),
                                     EP=as.numeric(as.character(EP)),
                                     TC=as.numeric(as.character(TC)),
                                     Z9=as.numeric(as.character(Z9))))
}
wos_results_marine = wos_results_marine %>%
  dplyr::mutate(Ecosystem="Marine")

wos_results = wos_results_marine %>%
  dplyr::mutate(U1=as.numeric(as.character(U1))) %>%
  dplyr::mutate(SN=as.numeric(as.character(SN))) %>%
  dplyr::bind_rows(wos_results_fresh)

abstract_categorization = wos_results %>%
  dplyr::mutate(AB=as.character(AB)) %>% rowwise() %>%
  dplyr::mutate(aquatic_system = case_when(grepl("ocean",tolower(AB))~"marine",
                                           grepl("marine",tolower(AB))~"marine",
                                           grepl("estuar",tolower(AB))~"estuarine",
                                           grepl("lake",tolower(AB))~"lake",
                                              grepl("lakes",tolower(AB))~"lake",
                                              grepl("pond",tolower(AB))~"pond",
                                              grepl("river",tolower(AB))~"river",
                                              grepl("stream",tolower(AB))~"stream",
                                              grepl("freshwater",tolower(AB))~"other freshwater",
                                              TRUE~"other")) %>%
  dplyr::mutate(omics_type = case_when(grepl("16s",tolower(AB))|grepl("18s",tolower(AB))|
                                         grepl("amplicon",tolower(AB))~"metabarcoding",
                                       grepl("metatranscript",tolower(AB))~"metatranscriptomics",
                                       grepl("metagenom",tolower(AB))~"metagenomics",
                                       grepl("genom",tolower(AB))~"genomics",
                                       grepl("transcript",tolower(AB))~"transcriptomics",
                                              TRUE~"other"))
  
ggplot(abstract_categorization %>% 
         dplyr::filter(Ecosystem=="Freshwater") %>%
         dplyr::group_by(PY,omics_type,aquatic_system) %>% 
         dplyr::tally() %>% dplyr::filter(freshwater_system!="other")) + 
    geom_bar(aes(x=PY,y=n,group=paste0(PY,freshwater_system),
                 fill=freshwater_system),stat="identity")+
    facet_wrap(~omics_type)

ggplot(abstract_categorization %>% 
         dplyr::filter(Ecosystem=="Freshwater") %>%
         dplyr::group_by(PY,omics_type,aquatic_system) %>% 
         dplyr::tally() %>% dplyr::filter(freshwater_system!="other")) + 
  geom_bar(aes(x=PY,y=n,group=paste0(PY,freshwater_system),
               fill=freshwater_system),stat="identity")+
  facet_wrap(~omics_type)

ggplot(abstract_categorization %>% 
         dplyr::filter(Ecosystem=="Marine") %>%
         dplyr::filter(PY>1950) %>%
         dplyr::group_by(PY,omics_type,aquatic_system) %>% 
         dplyr::tally() %>% dplyr::filter(aquatic_system!="other")) + 
  geom_bar(aes(x=PY,y=n,group=paste0(PY,aquatic_system),
               fill=aquatic_system),stat="identity")+
  facet_wrap(~omics_type)

ggplot(abstract_categorization %>% 
         dplyr::filter(omics_type=="metabarcoding") %>%
         dplyr::filter(PY>1950) %>%
         dplyr::group_by(PY,omics_type,aquatic_system,Ecosystem) %>% 
         dplyr::tally() %>% dplyr::filter(aquatic_system!="other")) + 
  geom_bar(aes(x=PY,y=n,group=paste0(PY,aquatic_system),
               fill=aquatic_system),stat="identity")+
  facet_wrap(~Ecosystem,nrow=2)
