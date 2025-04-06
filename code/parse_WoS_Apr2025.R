pacman::p_load(dplyr,ggplot2)

wos_results = data.frame()
data_dir="data/WoS_freshwater_5Apr2025"
for (curr in list.files(data_dir)) {
  wos_results=wos_results %>%
    dplyr::bind_rows(read.csv(file.path(data_dir,curr),sep="\t") %>%
                       dplyr::mutate(CY=as.numeric(as.character(CY)),
                                     NR=as.numeric(as.character(NR)),
                                     VL=as.numeric(as.character(VL)),
                                     PN=as.numeric(as.character(PN)),
                                     SU=as.numeric(as.character(SU)),
                                     BP=as.numeric(as.character(BP)),
                                     EP=as.numeric(as.character(EP))))
}

abstract_categorization = wos_results %>%
  dplyr::mutate(freshwater_system = case_when(grepl("lake",tolower(AB))~"lake",
                                              grepl("lakes",tolower(AB))~"lake",
                                              grepl("pond",tolower(AB))~"pond",
                                              grepl("river",tolower(AB))~"river",
                                              grepl("stream",tolower(AB))~"stream",
                                              grepl("freshwater",tolower(AB))~"other freshwater",
                                              TRUE~"other"))
  