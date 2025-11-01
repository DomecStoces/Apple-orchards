# Fourth corner analysis
library(readxl)
library(ade4)
library(dplyr)

# Load all sheets into a list
fourth_corner <- lapply(c("sp", "env", "traits"), function(x) read_excel("fourth_corner.xlsx", sheet = x))

# Assign names to the list elements
names(fourth_corner) <- c("sp", "env", "traits")

# Access individual sheets with the following syntax
dim(fourth_corner$sp)     # Dimensions of 'sp' sheet
dim(fourth_corner$env)    # Dimensions of 'env' sheet
dim(fourth_corner$traits) # Dimensions of 'traits' sheet

# Add a scaled squared altitude column to env
fourth_corner$env <- as.data.frame(lapply(fourth_corner$env, function(x) {
  if (is.character(x)) as.factor(x) else x
}))
fourth_corner$env <- fourth_corner$env %>%
  mutate(
    Altitude_scaled = scale(Altitude, center = TRUE, scale = TRUE)[,1],
    Altitude_scaled2 = Altitude_scaled^2
  )
fourth_corner$env <- fourth_corner$env %>%
  select(-Altitude) %>%   
  relocate(Altitude_scaled, Altitude_scaled2, .after = Wind)

fourth_corner$traits <- as.data.frame(fourth_corner$traits)
fourth_corner$traits$Body.size <- as.numeric(gsub(",", ".", fourth_corner$traits$Body.size))

# Convert numeric categorical columns to factors
fourth_corner$env$Exposition2 <- as.numeric(fourth_corner$env$Exposition2)
fourth_corner$env$Altitude <- as.numeric(fourth_corner$env$Altitude)
fourth_corner$env$Temperature <- as.numeric(fourth_corner$env$Temperature)
fourth_corner$env$Precipitation <- as.numeric(fourth_corner$env$Precipitation)
fourth_corner$env$Wind <- as.numeric(fourth_corner$env$Wind)


fourth_corner$sp[is.na(fourth_corner$sp)] <- 0
fourth_corner$traits[is.na(fourth_corner$traits)] <- 0
fourth_corner$traits_num <- fourth_corner$traits %>%
  select(-Species) %>%
  mutate(across(everything(), as.numeric))
afcL.aravo <- dudi.coa(fourth_corner$sp, scannf = FALSE)
acpR.aravo <- dudi.hillsmith(fourth_corner$env, row.w = afcL.aravo$lw,
                             scannf = FALSE)
acpQ.aravo <- dudi.pca(fourth_corner$traits_num, 
                       row.w = afcL.aravo$cw, 
                       scannf = FALSE)
rlq.aravo <- rlq(acpR.aravo, afcL.aravo, acpQ.aravo,
                 scannf = FALSE)

nrepet <- 999
four.comb.aravo <- fourthcorner(fourth_corner$env, fourth_corner$sp,
                                fourth_corner$traits_num, modeltype = 6,
                                p.adjust.method.G = "none",
                                p.adjust.method.D = "none", nrepet = nrepet)


plot(four.comb.aravo, alpha = 0.05, stat = "D2")
plot(four.comb.aravo, x.rlq = rlq.aravo, alpha = 0.05,
     stat = "D2", type = "biplot")

Srlq <- fourthcorner2(fourth_corner$env, fourth_corner$sp,
                      fourth_corner$traits,
                      modeltype = 6, p.adjust.method.G = "fdr", nrepet = nrepet)
Srlq$trRLQ

tiff('table.tiff', units="in", width=8, height=10, res=600)
plot(four.comb.aravo, alpha = 0.05, stat = "D2")
dev.off()

pdf('table.pdf', width=8, height=10)
plot(four.comb.aravo, alpha = 0.05, stat = "D2")
dev.off()