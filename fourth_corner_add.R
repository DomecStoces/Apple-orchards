# Fourth corner analysis with fuzzy traits
library(readxl)
library(dplyr)
library(tidyr)
library(ade4)

# --- 1) Read -----------------------------------------------------------------
fourth_corner <- lapply(c("sp", "env", "traits"),
                        function(x) read_excel("fourth_corner.xlsx", sheet = x))
names(fourth_corner) <- c("sp", "env", "traits")

# Base data.frames
fourth_corner$sp     <- as.data.frame(fourth_corner$sp)
fourth_corner$env    <- as.data.frame(fourth_corner$env)
fourth_corner$traits <- as.data.frame(fourth_corner$traits)

# --- 2) Checks ---------------------------------------------------------------
stopifnot(nrow(fourth_corner$env) == nrow(fourth_corner$sp))   # sites must match

# --- 3) ENV (R) --------------------------------------------------------------
# characters -> factors, keep numerics numeric
fourth_corner$env[] <- lapply(fourth_corner$env, function(x) if (is.character(x)) factor(x) else x)

# If Month exists, create Month_scaled then DROP Month
if ("Month" %in% names(fourth_corner$env)) {
  # coerce to numeric if needed
  if (is.factor(fourth_corner$env$Month) || is.character(fourth_corner$env$Month)) {
    suppressWarnings(fourth_corner$env$Month <- as.numeric(as.character(fourth_corner$env$Month)))
  }
  fourth_corner$env$Month_scaled <- as.numeric(scale(fourth_corner$env$Month))
  fourth_corner$env$Month <- NULL                # <-- remove raw Month
}

# Coerce other numerics if they arrived as factor/char
for (v in c("Trap")) {
  if (v %in% names(fourth_corner$env)) {
    if (is.factor(fourth_corner$env[[v]]) || is.character(fourth_corner$env[[v]])) {
      suppressWarnings(fourth_corner$env[[v]] <- as.numeric(as.character(fourth_corner$env[[v]])))
    }
  }
}
if ("Treatment" %in% names(fourth_corner$env)) fourth_corner$env$Treatment <- factor(fourth_corner$env$Treatment)
if ("Region"    %in% names(fourth_corner$env)) fourth_corner$env$Region    <- factor(fourth_corner$env$Region)

# --- 4) SP (L) ---------------------------------------------------------------
fourth_corner$sp[is.na(fourth_corner$sp)] <- 0

# --- 5) TRAITS (Q) with fuzzy handling ---------------------------------------
traits_raw <- fourth_corner$traits

# Align rows of traits to species in sp
sp_id_col <- intersect(names(traits_raw), c("Species","species","Taxon","taxon"))
if (length(sp_id_col) == 1) {
  sp_id_col <- sp_id_col[[1]]
  rownames(traits_raw) <- make.names(traits_raw[[sp_id_col]], unique = TRUE)
  traits_num <- traits_raw %>% select(-all_of(sp_id_col))
} else {
  sp_species <- colnames(fourth_corner$sp)
  if (nrow(traits_raw) != length(sp_species)) {
    stop("No species ID column in 'traits', and row count (traits) != species count (sp). Add a species column to 'traits'.")
  }
  rownames(traits_raw) <- sp_species
  traits_num <- traits_raw
}

# Ensure numeric 0/1 (or proportions)
traits_num[] <- lapply(traits_num, function(x) as.numeric(as.character(x)))
traits_num[is.na(traits_num)] <- 0

# ----- helpers for fuzzy handling --------------------------------------------
# Create missing target columns if needed
ensure_cols <- function(df, cols) {
  missing <- setdiff(cols, colnames(df))
  for (m in missing) df[[m]] <- 0
  df
}

# Split a composite column 'comp' equally into two target columns 'to'
split_composite <- function(df, comp, to) {
  if (comp %in% colnames(df)) {
    df <- ensure_cols(df, to)
    w <- df[[comp]]
    df[[to[1]]] <- df[[to[1]]] + 0.5 * w
    df[[to[2]]] <- df[[to[2]]] + 0.5 * w
    df[[comp]]  <- NULL
  }
  df
}

# Row-normalize a set of columns so rows sum to 1 (if row sum > 0)
row_norm <- function(df, cols) {
  ok <- intersect(cols, colnames(df))
  if (length(ok) > 1) {
    rs <- rowSums(df[, ok, drop = FALSE])
    rs[rs == 0] <- 1
    df[, ok] <- sweep(df[, ok, drop = FALSE], 1, rs, "/")
  }
  df
}
# Standardize column names so we can refer to them safely (slashes, spaces, (…) etc.)
std_names <- tolower(gsub("[^a-z0-9]+", "_", names(traits_num)))
std_names <- gsub("_+", "_", gsub("^_|_$", "", std_names))
names(traits_num) <- make.unique(std_names, sep = "_")
base_names <- sub("_[0-9]+$", "", names(traits_num))
collapse_list <- split(seq_along(base_names), base_names)
# Make sure atomic egg-laying targets exist (created as 0 if missing)
traits_num <- do.call(
  cbind,
  lapply(collapse_list, function(idx) {
    m <- as.matrix(traits_num[, idx, drop = FALSE])
    # row-wise max; if you prefer sums, use rowSums(m)
    v <- apply(m, 1, max, na.rm = TRUE)
    as.numeric(v)
  })
)
traits_num <- as.data.frame(traits_num)
names(traits_num) <- names(collapse_list)

traits_num <- ensure_cols(
  traits_num,
  c("aquatic_plants",
    "bark","twig",
    "detritus","fungi","fruit",
    "moss","lichen",
    "leaf","leaf_coniferous_trees","leaf_deciduous_trees",
    "leaf_grass","leaf_herbs","leaf_shrubs",
    "base_grass","base_herbs","stem_grass","stem_herbs")
)

# Split composite egg-laying columns 50/50 into atomic targets (names are AFTER standardization)
traits_num <- split_composite(traits_num, "bark_or_fruit",      c("bark","fruit"))
traits_num <- split_composite(traits_num, "bark_or_leaves",     c("bark","leaf"))
traits_num <- split_composite(traits_num, "bark_or_twigs",      c("bark","twig"))
traits_num <- split_composite(traits_num, "detritus_or_fruit",  c("detritus","fruit"))
traits_num <- split_composite(traits_num, "detritus_or_fungi",  c("detritus","fungi"))
traits_num <- split_composite(traits_num, "detritus_or_moss",   c("detritus","moss"))
traits_num <- split_composite(traits_num, "moss_lichens",       c("moss","lichen"))
traits_num <- split_composite(traits_num, "grass_herbs",        c("grass","herbs"))
traits_num <- split_composite(traits_num, "base_grass_herbs",   c("base_grass","base_herbs"))
traits_num <- split_composite(traits_num, "leaf_or_stem_grass", c("leaf_grass","stem_grass"))
traits_num <- split_composite(traits_num, "leaf_or_stem_herbs", c("leaf_herbs","stem_herbs"))

# Define the final egg-laying set and fuzzy-normalize (rows sum to 1 across this set)
egg_cols <- intersect(
  c("aquatic_plants","bark","twig","detritus","fungi","fruit",
    "moss","lichen","leaf","leaf_coniferous_trees","leaf_deciduous_trees",
    "leaf_grass","leaf_herbs","leaf_shrubs","base_grass","base_herbs",
    "stem_grass","stem_herbs"),
  colnames(traits_num)
)
traits_num <- row_norm(traits_num, egg_cols)

grep("^lichen", names(traits_num), value = TRUE)  # should return just "lichen"
rs <- rowSums(traits_num[, egg_cols, drop = FALSE])
summary(rs[rs > 0])  # ~1 (allow tiny rounding error)
# ----- (2) Fuzzy normalization within trait groups ----------------------------
# Life-stage: keep as multi-label (0/1) unless you want exclusivity
# life_stage_cols <- c("egg","caterpillar","caterpillar/pupa","pupa","pupa / adult","adult")
# (no normalization by default)

# Save back for downstream
fourth_corner$traits_num <- traits_num

# --- 6) Align species between L and Q ----------------------------------------
L <- fourth_corner$sp
Q <- fourth_corner$traits_num
R <- fourth_corner$env

std <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

sp_names_L <- colnames(L)
sp_names_Q <- rownames(Q)
std_L <- std(sp_names_L)
std_Q <- std(sp_names_Q)
match_idx <- match(std_L, std_Q)
if (all(is.na(match_idx))) {
  # 2) fallback: align by position if dimensions agree
  if (nrow(Q) == ncol(L)) {
    rownames(Q) <- colnames(L)
    sp_names_Q  <- rownames(Q)
    std_Q       <- std(sp_names_Q)
    match_idx   <- seq_len(nrow(Q))     # identity
  } else {
    # 3) give a clear diagnostic and stop
    bad_L <- setdiff(unique(std_L), unique(std_Q))
    stop(
      "No overlapping species between L and Q.\n",
      "Examples only in L (std): ", paste(head(bad_L, 10), collapse=", "),
      "\nCheck that 'traits' species rows correspond to 'sp' columns."
    )
  }
}

# Keep only species that matched
keep_sp_in_L <- which(!is.na(match_idx))
if (length(keep_sp_in_L) < 2L) {
  stop("Too few matching species after name harmonization.")
}
L <- L[, keep_sp_in_L, drop = FALSE]
Q <- Q[match_idx[keep_sp_in_L], , drop = FALSE]  # reorder Q rows to matched order
rownames(Q) <- colnames(L)                        # sync names exactly

# 4) drop zero-abundance species (columns all 0) and mirror in Q
nonzero_species <- which(colSums(L, na.rm = TRUE) > 0)
if (length(nonzero_species) < ncol(L)) {
  L <- L[, nonzero_species, drop = FALSE]
  Q <- Q[colnames(L), , drop = FALSE]
}

# 5) drop empty sites (rows all 0) and mirror in R
nonzero_sites <- which(rowSums(L, na.rm = TRUE) > 0)
if (length(nonzero_sites) < nrow(L)) {
  L <- L[nonzero_sites, , drop = FALSE]
  R <- R[nonzero_sites, , drop = FALSE]
}

# final sanity
stopifnot(nrow(R) == nrow(L))
stopifnot(nrow(Q) == ncol(L))
# --- 7) RLQ -------------------------------------------------------------------
afcL <- dudi.coa(L, scannf = FALSE)
acpR <- dudi.hillsmith(R, row.w = afcL$lw, scannf = FALSE)
acpQ <- dudi.pca(Q, row.w = afcL$cw, scale = TRUE, scannf = FALSE)

rlq_res <- rlq(acpR, afcL, acpQ, scannf = FALSE)

# --- 8) Fourth-corner ---------------------------------------------------------
nrepet <- 999
fc <- fourthcorner(
  R, L, Q,
  modeltype = 6,
  p.adjust.method.G = "none",
  p.adjust.method.D = "none",
  nrepet = nrepet
)

plot(fc, alpha = 0.05, stat = "D2")
plot(fc, x.rlq = rlq_res, alpha = 0.05, stat = "D2", type = "biplot")

fc2 <- fourthcorner2(R, L, Q, modeltype = 6, p.adjust.method.G = "fdr", nrepet = nrepet)
fc2$trRLQ

# --- 9) Save figures ----------------------------------------------------------
tiff("table.tiff", units = "in", width = 8, height = 10, res = 600)
plot(fc, alpha = 0.05, stat = "D2")
dev.off()

pdf("table.pdf", width = 8, height = 10)
plot(fc, alpha = 0.05, stat = "D2")
dev.off()