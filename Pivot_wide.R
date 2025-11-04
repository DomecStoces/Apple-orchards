library(tidyr)
library(dplyr)
library(readxl)
library(writexl)
library(stringr)

df <- read_excel("Sady_final.xlsx", sheet = "Motýli")

# Pivot to wide format
wide_df <- df %>%
  mutate(Number = as.numeric(Number)) %>%
  # collapse any accidental duplicates first
  group_by(Lokalita, Trap, Month, Treatment, Region, T1, T2, T3, Moisture, Grass, Herbs, Shrubs, Species) %>%
  summarise(Number = sum(Number), .groups = "drop") %>%
  pivot_wider(
    names_from  = Species,
    values_from = Number,
    values_fill = list(Number = 0),   
    values_fn   = list(Number = sum)
  )
# Check result
head(wide_df)

meta_cols <- c("Lokalita", "Trap", "Month", "Treatment", "Region",
               "T1", "T2", "T3", "Moisture", "Grass", "Herbs", "Shrubs")

# Extract species columns and reorder alphabetically
species_cols <- setdiff(colnames(wide_df), meta_cols)
wide_df_sorted <- wide_df[, c(meta_cols, sort(species_cols))]

# Write to Excel file
write_xlsx(wide_df_sorted, "Sady_final_sorted.xlsx")

# Check
file.exists("Sady_final_sorted.xlsx")
# ----Write traits for each species in long format----
traits<-read_excel("traits.xlsx", sheet = "List1")

traits_clean <- traits %>%
  rename(
    Species        = Name,
    Overwinter     = Overwinter,
    Feeding        = Feeding,
    Eggs_laying    = Eggs_laying,
    Pest           = Pest,
    Migration      = Migration,
    Live_in_orchard= `Live in orchard`
  ) %>%
  mutate(
    Species = str_squish(Species),
    Species = str_replace_all(Species, "\\s+", " ")
  )

sady_clean <- Sady_final %>%
  mutate(
    Species = str_squish(Species),
    Species = str_replace_all(Species, "\\s+", " ")
  )
# Duplicated traits should be 0
dup_traits <- traits_clean %>% count(Species) %>% filter(n > 1)

# Join
Sady_with_traits <- sady_clean %>%
  left_join(traits_clean, by = "Species")

# Which species did not match traits?
unmatched_species <- Sady_with_traits %>%
  filter(is.na(Overwinter) & is.na(Feeding) & is.na(Eggs_laying) &
           is.na(Pest) & is.na(Migration) & is.na(Live_in_orchard)) %>%
  distinct(Species) %>%
  arrange(Species)

# Save as .xlsx
writexl::write_xlsx(
  list(
    Sady_with_traits = Sady_with_traits,
    Traits_duplicates = dup_traits,
    Unmatched_species = unmatched_species
  ),
  path = "Sady_final_with_traits.xlsx"
)

# Add phylogenetic relationship to long format
library(readxl)
library(dplyr)
library(stringr)
library(janitor)
library(writexl)
library(janitor)

# --- 1) Load and standardize column names
sady <- read_excel("Sady_final.xlsx")   |> clean_names()   # -> species, number, lokalita, t1, t2, ...
taxa <- read_excel("Fylogeneze.xlsx")  |> clean_names()    # -> name, tribus_subfamily, family

# --- 2) Make the join key the same in both tables
sady <- sady |> mutate(species = str_squish(species))

taxa <- taxa |>
  rename(
    species = name,
    tribe   = tribus_subfamily,
    family  = family
  ) |>
  mutate(species = str_squish(species))

# --- 3) Fix decimals read with commas (only if those columns exist)
num_cols <- c("t1","t2","t3","moisture","grass","herbs","shrubs")
present   <- intersect(num_cols, names(sady))
if (length(present)) {
  sady <- sady |>
    mutate(across(all_of(present),
                  ~ suppressWarnings(as.numeric(str_replace(as.character(.), ",", ".")))))
}

# --- 4) Join taxonomy into every row of the long table
sady_tax <- left_join(sady, taxa, by = "species")

# --- 5) Check what didn’t match (typos / missing taxonomy)
not_matched <- anti_join(sady, taxa, by = "species") |> distinct(species)
print(not_matched)

# --- 6) Save result (plus an 'Unmatched' sheet for quick review)
write_xlsx(
  list(Sady_with_taxa = sady_tax,
       Unmatched      = not_matched),
  "Sady_final_with_taxa.xlsx"
)

# Join beta diversity components to Sady_final
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library/writexl)

clean_key <- function(x) {
  x <- gsub("\\p{Pd}", "-", x, perl = TRUE)
  x <- trimws(x)
  x <- gsub("\\s+", " ", x)
  toupper(x)
}

beta_lookup <- read_excel("beta_components.xlsx", sheet = 1) %>%
  mutate(
    management = recode(management, "EKO" = "BIO"),
    Component  = recode(Component,
                        "βdiversity" = "beta_diversity",
                        "Turnover"   = "turnover",
                        "Nestedness" = "nestedness"),
    # normalize keys
    Locality   = clean_key(Locality),
    management = clean_key(management)
  ) %>%
  # fix the Locality prefix for BIO sites (EKO-xxx -> BIO-xxx)
  mutate(Locality = ifelse(management == "BIO",
                           sub("^EKO-", "BIO-", Locality),
                           Locality)) %>%
  arrange(Locality, management, Component, ID) %>%
  distinct(Locality, management, Component, .keep_all = TRUE) %>%
  select(Locality, management, Component, Distance) %>%
  pivot_wider(names_from = Component, values_from = Distance) %>%
  distinct(Locality, management, .keep_all = TRUE)

Sady_joined <- Sady_final %>%
  mutate(
    lokalita  = clean_key(lokalita),
    treatment = clean_key(treatment)
  ) %>%
  left_join(beta_lookup, by = c("lokalita" = "Locality", "treatment" = "management"))

# Recheck missing (should now be zero)
Sady_joined %>%
  distinct(lokalita, treatment, beta_diversity, turnover, nestedness) %>%
  filter(is.na(beta_diversity) & is.na(turnover) & is.na(nestedness))

Sady_final2 <- Sady_joined %>%
  group_by(lokalita, treatment) %>%
  mutate(across(c(beta_diversity, nestedness, turnover),
                ~ if (all(is.na(.))) NA_real_ else rep(na.omit(.)[1], dplyr::n()))) %>%
  ungroup()
# 5) Save to Excel
write_xlsx(Sady_final2, "Sady_final_with_beta.xlsx")
