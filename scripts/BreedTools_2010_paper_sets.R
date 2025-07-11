require("BIGr")
require("scales")
require("tidyverse")

# read in data
reference <- as.data.frame(
  read.table(
    "https://raw.githubusercontent.com/josuechinchilla/RHBreedTools/main/data/2010_reference.txt",
    header = TRUE,
    row.names = 1,
    sep = "\t"
  ) %>%
    dplyr::select(-c(1, 2)) %>%
    t()
)

validation <- as.data.frame(
  read.table(
    "~/Desktop/GitHub/RHBreedTools/data/2010_validation.txt",
    header = TRUE,
    row.names = 1,
    sep = "\t"
  ) %>%
    dplyr::select(-c(1, 2)) %>%
    t()
)

reference_ids <- read.table(
    "https://raw.githubusercontent.com/josuechinchilla/RHBreedTools/main/data/ref_ids.txt",
    header = TRUE,
    sep = "\t"
  ) 


ref_ids = lapply(as.list(reference_ids),as.character)


# Function to format percentages with one decimal place
format_percent <- function(x) {
  percent_format(accuracy = 0.1)(x)
}

#calculate allele frequency of reference panel
freq = BIGr:::allele_freq_poly (reference, ref_ids, ploidy = 2)
prediction = as.data.frame(BIGr:::solve_composition_poly(validation,freq))


# Apply the formatting function to all columns in the dataframe and create "assigned_breed" column as the breed with the highest value

columns_to_select <- names(prediction)[-which(names(prediction) == "R2")]  # Exclude column R2

pred_results <- prediction %>%
  as.data.frame() %>%  # in case it's a matrix
  rownames_to_column(var = "ID") %>%
  mutate(
    across(-ID, ~format_percent(.)),
    predicted_line = columns_to_select[max.col(select(., all_of(columns_to_select)), ties.method = "first")]
  )

#create an output plot
# Clean and reshape
pred_results_long <- pred_results %>%
  pivot_longer(cols = c(non.RHB, RHB), names_to = "category", values_to = "percent") %>%
  mutate(
    percent = as.numeric(str_remove(percent, "%"))
  )

# Plot
ggplot(pred_results_long, aes(x = ID, y = percent, fill = category)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1, scale = 1)) +
  labs(
    x = "Individual ID",
    y = "Ancestry Proportion",
    fill = "Ancestry Category"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

#Write our breed composition results
write.table(pred_results_formatted, "honeybee_regression_results_2010_paper_sets.txt", col.names = T, sep = "\t", quote = F, row.names = F)



