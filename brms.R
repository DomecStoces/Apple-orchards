# adding results of βdiversity, Turnover, Nestedness analyses to brms models
library(dplyr)
library(tidyr)
library(readr)
library(readxl)

# --- 1) Read your beta table ---
beta_raw <- read_excel("beta_components.xlsx", sheet = 1)

# --- 2) Clean up and harmonize ---
beta_clean <- beta_raw %>%
  # normalize column names (remove extra spaces)
  rename_with(~trimws(.x)) %>%
  mutate(
    # fix decimal commas (if present)
    Distance = parse_number(as.character(Distance), locale = locale(decimal_mark = ",")),
    # align management with Sady_final: EKO -> BIO
    management = recode(management, "EKO" = "BIO"),
    # clean component names
    Component = recode(Component,
                       "βdiversity" = "beta_diversity",
                       "Turnover"   = "turnover",
                       "Nestedness" = "nestedness")
  )

# --- 3) Aggregate duplicates and pivot wider ---
beta_wide <- beta_clean %>%
  group_by(Locality, management, Component) %>%
  summarise(Distance = mean(Distance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Component, values_from = Distance)

# --- 4) Join with your main dataset ---
Sady_final2 <- Sady_final %>%
  left_join(beta_wide, by = c("lokalita" = "Locality", "treatment" = "management"))

# --- 5) Check for unmatched rows ---
sum(is.na(Sady_final2$beta_diversity) & is.na(Sady_final2$turnover) & is.na(Sady_final2$nestedness))


# Pseudo-phylogenetic tree
library(ape)
library(brms)
library(readxl)

Sady_final <- read_excel("Sady_final_with_taxa.xlsx", sheet = "Sheet1")


# Build a taxonomy-based tree (several routes; one is as.phylo on a classification):
tr  <- as.phylo(~ Family/Genus/Species, data = tax_df)
tr  <- multi2di(tr)                                
tr  <- compute.brlen(tr, method = "Grafen")        
vcv_mat <- vcv(tr, corr = TRUE)                    

mod1<-brm(
  bf(Beta-diversity ~ Treatment + Month + (1|gr(Species, cov = vcv_mat)) + (1|Region/Trap)),
  data = df, data2 = list(vcv_mat = vcv_mat),
  family = gaussian(), control = list(adapt_delta = 0.95))

