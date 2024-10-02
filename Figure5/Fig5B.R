#Script for Fig 5B - top subvariants heatmap and bar plot.

#Need to first make aggregated data of 'subvariants', so the +XYZ calls we assigned to samples.

##import v25 data
alldatareload<-read.xlsx("filtered_data_for_paper_v25.xlsx")
alldatareload$Date.of.collection<-as.Date(alldatareload$Date.of.collection, origin = "1899-12-30") #

###get count of all subvars and subvars per main var###
allsubvar<-table(re_filtered_maindataset$Variant.Details)
allsubvar<-allsubvar[order(allsubvar, decreasing = T)]
write.xlsx(allsubvar, file="AllSubvariants_v25.xlsx")

##then break down by main var
#so group by variant name new and then variant details within that
re_filtered_maindataset$Variant_Name_New_Cut<-gsub("\\+", "", re_filtered_maindataset$Variant_Name_New)

data_for_bp<-as.data.frame(re_filtered_maindataset %>%
   group_by(Variant_Name_New_Cut) %>%
   dplyr::count(Variant.Details))
#organized subvar data
write.xlsx(data_for_bp, "organized_subvariants_v25.xlsx")


###Begin heatmap####
