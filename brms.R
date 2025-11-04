library(dplyr)
library(tidyr)
library(readxl)

beta_raw <- read_excel("beta_components.xlsx", sheet = 1)

# Make components consistent, keep the first ID within each trio,
# then pivot wider and de-duplicate keys.
beta_wide_unique <- beta_raw %>%
  mutate(
    management = recode(management, "EKO" = "BIO"),
    Component  = recode(Component,
                        "βdiversity" = "beta_diversity",
                        "Turnover"   = "turnover",
                        "Nestedness" = "nestedness")
  ) %>%
  arrange(Locality, management, Component, ID) %>%
  distinct(Locality, management, Component, .keep_all = TRUE) %>%   # first per trio
  select(Locality, management, Component, Distance) %>%
  pivot_wider(names_from = Component, values_from = Distance) %>%
  distinct(Locality, management, .keep_all = TRUE)                  # one row per key

# Sanity checks
stopifnot(!any(duplicated(beta_wide_unique[c("Locality","management")])))

# Now the join is many-to-one (no warning)
Sady_final2 <- Sady_final %>%
  mutate(
    lokalita  = trimws(lokalita),
    treatment = toupper(trimws(treatment))
  ) %>%
  left_join(beta_wide_unique, by = c("lokalita" = "Locality",
                                     "treatment" = "management"))

# Optional: verify BIO-BRN-01 numbers
Sady_final2 %>%
  filter(lokalita == "BIO-BRN-01", treatment == "BIO") %>%
  distinct(lokalita, treatment, beta_diversity, turnover, nestedness)

# 5) Save
write_xlsx(Sady_final2, "Sady_final_with_beta.xlsx")

# Pseudo-phylogenetic tree
# --- packages ---
library(dplyr)
library(tidyr)
library(ape)      
library(brms)
# ========== 1) Build a clean taxonomy table ==========
df <- Sady_final2

# Derive Species as the full binomial (already in df$species)
# Keep only required taxonomic columns
tax0 <- df %>%
  transmute(
    Species_chr = as.character(species),
    Family_chr  = as.character(family),
    Tribe_chr   = as.character(tribe)
  )

# Clean encodings and whitespace; replace missing/blank with placeholders
clean_text <- function(x) {
  x <- ifelse(is.na(x), "", x)
  x <- trimws(x)
  # normalize encoding to ASCII-ish; keep original if conversion fails
  y <- iconv(x, from = "", to = "UTF-8", sub = "")
  y[is.na(y)] <- x[is.na(y)]
  y <- gsub("\\s+", " ", y)
  y
}

tax1 <- tax0 %>%
  mutate(
    Species_chr = clean_text(Species_chr),
    Family_chr  = clean_text(Family_chr),
    Tribe_chr   = clean_text(Tribe_chr),
    # Fill empty taxa with explicit placeholders (as.phylo dislikes blanks)
    Family_chr  = ifelse(Family_chr == "", "FAMILY_UNKNOWN", Family_chr),
    Tribe_chr   = ifelse(Tribe_chr  == "", "TRIBE_UNKNOWN",  Tribe_chr)
  ) %>%
  # One row per species; if a species maps to multiple families/tribes, keep the most frequent
  group_by(Species_chr, Family_chr, Tribe_chr) %>%
  tally() %>%
  group_by(Species_chr) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Species_chr, Family_chr, Tribe_chr)

# Coerce to factors (as.phylo requires factors)
tax_df <- tax1 %>%
  mutate(
    Species = factor(Species_chr),
    Family  = factor(Family_chr),
    Tribe   = factor(Tribe_chr)
  ) %>%
  select(Family, Tribe, Species)

# Quick sanity checks
stopifnot(!any(is.na(tax_df$Species)), !any(is.na(tax_df$Family)), !any(is.na(tax_df$Tribe)))

# ========== 2) Build the pseudo-phylogeny with fallbacks ==========
build_tree_safe <- function(tax_df) {
  # Try Family/Tribe/Species
  tr <- try(as.phylo(~ Family/Tribe/Species, data = tax_df), silent = TRUE)
  if (inherits(tr, "try-error")) {
    # Fallback 1: Family/Species (if Tribe is too sparse)
    tr <- try(as.phylo(~ Family/Species, data = tax_df), silent = TRUE)
  }
  if (inherits(tr, "try-error")) {
    stop("Failed to build classification tree. Check that Family/Tribe/Species are factors with no NAs.")
  }
  tr <- multi2di(tr)
  tr <- compute.brlen(tr, method = "Grafen")
  tr
}

tr <- build_tree_safe(tax_df)

# ========== 3) Correlation VCV and alignment to data ==========
vcv_mat <- vcv(tr, corr = TRUE)

# Align matrix to the species present in your modeling data
keep <- intersect(rownames(vcv_mat), df$species)
vcv_mat <- vcv_mat[keep, keep, drop = FALSE]

# Prepare modeling frame aligned to VCV order
df_mod <- df %>%
  filter(species %in% rownames(vcv_mat)) %>%
  mutate(
    Species   = factor(species, levels = rownames(vcv_mat)),
    Treatment = toupper(treatment),
    Month     = as.numeric(month),
    Region    = factor(region),
    Trap      = factor(trap)
  )
stopifnot(all(levels(df_mod$Species) == rownames(vcv_mat)))

set.seed(123)
mod1 <- brm(
  formula = bf(beta_diversity ~ treatment + (1|Month) + (1 | gr(Species, cov = vcv_mat)) + (1 | region)),
  data    = df_mod,
  data2   = list(vcv_mat = vcv_mat),
  family  = Beta(),
  control = list(adapt_delta = 0.95),
  chains  = 4, iter = 3000
)

# =============== 5) Post-fit checks ===========================
print(mod1)
pp_check(mod1)
summary(mod1)

