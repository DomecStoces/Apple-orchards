library(ape)
library(brms)
# tax_df has columns: Family, Genus, Species (matching your data$Species)
# Build a taxonomy-based tree (several routes; one is as.phylo on a classification):
tr  <- as.phylo(~ Family/Genus/Species, data = tax_df)
tr  <- multi2di(tr)                                # resolve polytomies deterministically
tr  <- compute.brlen(tr, method = "Grafen")        # positive branch lengths, ultrametric-like
vcv_mat <- vcv(tr, corr = TRUE)                    # correlation matrix for brms

mod1<-brm(
  bf(Beta-diversity ~ Treatment + Month + (1|gr(Species, cov = vcv_mat)) + (1|Region/Trap)),
  data = df, data2 = list(vcv_mat = vcv_mat),
  family = gaussian(), control = list(adapt_delta = 0.95))