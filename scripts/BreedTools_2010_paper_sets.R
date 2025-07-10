library(breedTools)
library(scales)
library(tidyverse)

# read in data
reference = read.table("https://github.com/josuechinchilla/RHBreedTools/blob/main/data/2010_reference.txt", header = T, row.names = 1, sep = "\t")

validation = read.table("https://github.com/josuechinchilla/RHBreedTools/blob/main/data/2010_validation.txt", header = T, row.names = 1, sep = "\t")

reference_ids = read.table("https://github.com/josuechinchilla/RHBreedTools/blob/main/data/ref_ids.txt", header = T, sep = "\t") 

ref_ids = lapply(as.list(reference_ids),as.character)


# Function to format percentages with one decimal place
format_percent <- function(x) {
  percent_format(accuracy = 0.1)(x)
}

#calculate allele frequency of reference panel
freq2 = allele_freq_poly(reference, ref_ids)

devtools#calculate breed composition of validation set
prediction2 = as.data.frame(solve_composition_poly(validation,freq2))


# Apply the formatting function to all columns in the dataframe and create "assigned_breed" column as the breed with the highest value

columns_to_select <- names(prediction)[-which(names(prediction) == "R2")]  # Exclude column R2

pred_results <- prediction %>%
  mutate(
    across(everything(), ~format_percent(.)),
    predicted_line = columns_to_select[max.col(prediction[columns_to_select], ties.method = "first")],
    ID = row.names(prediction)
  ) 

pred_results = inner_join(pred_results, raw_data, by = "ID")


pred_results_formatted = pred_results %>%
  mutate(
    line_match = ifelse(predicted_line == original_line,"match", "no_match")
  ) %>%
  select(c(5,6,4,7,3,1,2)) #adjust the col order of the final file, the number of columns depends on the number of reference populations.
  

#Write our breed composition results
write.table(pred_results_formatted, "honeybee_regression_results_2010_paper_sets.txt", col.names = T, sep = "\t", quote = F, row.names = F)



