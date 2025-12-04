path1 = "~/Downloads/data_baseline_change_old.Rdata"
old <-  as.data.frame(longForm(get(load(path1))[,,"Composite methylation scores: blood"],
                               colData = c("interventionId", "subjectId", "visitId", "time", "compliance", "smoking_py",
                                           "basename_blood", "age_at_consent",
                                           "nic_replacement",
                                           "smkstop",
                                           "cessation_date"))) |> 
  filter_smk(age_norm = T, type_norm = T)


path = here("data/data_baseline_change.Rdata")
df_change_blood <- as.data.frame(longForm(get(load(path))[,,"Composite methylation scores: blood"],
                                          colData = c("interventionId", "subjectId", "visitId", "time", "compliance", "smoking_py",
                                                      "basename_blood", "age_at_consent",
                                                      "nic_replacement",
                                                      "smkstop",
                                                      "cessation_date"))) |> 
  filter_smk(age_norm = T, type_norm = T)

intersect <- intersect(old$basename_blood, df_change_blood$basename_blood)


old1 <- old |> dplyr::filter(rowname == 'AgeAccelGrimV2') |> 
  dplyr::select(basename_blood, value)
new1 <- df_change_blood |> dplyr::filter(rowname == 'AgeAccelGrimV2') |> 
  dplyr::select(basename_blood, value)

merge <- old1 |> dplyr::left_join(new1, by = 'basename_blood')

identical(merge$value.x, merge$value.y)
merge |> 
  ggplot(aes(x = value.x,
             y = value.y)) +
  geom_point()
