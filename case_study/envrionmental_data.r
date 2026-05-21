# Environmental variables and per-species PCA-climate (PC_climate1, PC_climate2)
# Follows Ramirez-Valiente et al. (2023): 19 WorldClim BIOCLIM + soil pH (0-30 cm)
# + annual and summer moisture indices (PET via Hargreaves-Samani 1985),
# one prcomp per species, centred and scaled.
# Soil pH comes from SoilGrids 2.0 instead of Trabucco & Zomer (2010) for reproducibility.

library(terra)
library(geodata)

this_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
setwd(this_dir)

pop_file <- "pop_id_withhabitat.csv"

indiv <- read.csv(pop_file, stringsAsFactors = FALSE, check.names = FALSE)
indiv$Latitude  <- as.numeric(gsub(",", ".", indiv$Latitude))
indiv$Longitude <- as.numeric(gsub(",", ".", indiv$Longitude))

pops <- unique(indiv[, c("Species", "PopulationCode", "Latitude", "Longitude")])
pops <- pops[!is.na(pops$Latitude) & !is.na(pops$Longitude), ]
rownames(pops) <- NULL

pts <- vect(pops, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")


# Soil pH (0-30 cm): thickness-weighted mean of SoilGrids slices 0-5, 5-15, 15-30 cm.
# SoilGrids encodes pH * 10.
soil_path <- "soil"
dir.create(soil_path, showWarnings = FALSE)

soil_depths <- list(list(depth =  5, thick =  5),
                    list(depth = 15, thick = 10),
                    list(depth = 30, thick = 15))

soil_vals_per_layer <- sapply(soil_depths, function(d) {
  r <- soil_world(var = "phh2o", depth = d$depth, path = soil_path)
  terra::extract(r, pts)[, 2]
})
weights <- sapply(soil_depths, function(d) d$thick)
pops$soil_pH <- as.numeric((soil_vals_per_layer %*% weights) / sum(weights) / 10)


# WorldClim 2.1 (1970-2000). Increase resolution for finer extraction (e.g. 0.5).
wc_res <- 10
bio  <- worldclim_global(var = "bio",  res = wc_res, version = "2.1", path = "wc2")
prec <- worldclim_global(var = "prec", res = wc_res, version = "2.1", path = "wc2")
tmin <- worldclim_global(var = "tmin", res = wc_res, version = "2.1", path = "wc2")
tmax <- worldclim_global(var = "tmax", res = wc_res, version = "2.1", path = "wc2")

bio_vals  <- terra::extract(bio,  pts)[, -1]
colnames(bio_vals) <- paste0("bio", 1:19)
pr_vals   <- terra::extract(prec, pts)[, -1]
tmin_vals <- terra::extract(tmin, pts)[, -1]
tmax_vals <- terra::extract(tmax, pts)[, -1]


# Hargreaves-Samani PET (mm/month). Ra computed via FAO/Allen et al. (1998).
deg2rad <- function(x) x * pi / 180

ra_monthly <- function(lat_deg, month) {
  lat     <- deg2rad(lat_deg)
  doy_mid <- c(15, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349)
  J       <- doy_mid[month]
  dr    <- 1 + 0.033 * cos(2 * pi * J / 365)
  delta <- 0.409 * sin(2 * pi * J / 365 - 1.39)
  ws    <- acos(-tan(lat) * tan(delta))
  Gsc <- 0.0820
  (24 * 60 / pi) * Gsc * dr *
    (ws * sin(lat) * sin(delta) + cos(lat) * cos(delta) * sin(ws))
}

days_in_month <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

compute_pet_hs <- function(tmin_mat, tmax_mat, lat_vec) {
  pet_mat <- matrix(NA_real_, nrow(tmin_mat), 12)
  for (i in seq_len(nrow(tmin_mat))) {
    Tmean  <- (tmin_mat[i, ] + tmax_mat[i, ]) / 2
    Trange <- pmax(tmax_mat[i, ] - tmin_mat[i, ], 0)
    Ra_mm  <- sapply(1:12, function(m) 0.408 * ra_monthly(lat_vec[i], m))
    pet_d  <- 0.0023 * Ra_mm * (Tmean + 17.8) * sqrt(Trange)
    pet_mat[i, ] <- pet_d * days_in_month
  }
  pet_mat
}

pet_vals <- compute_pet_hs(as.matrix(tmin_vals),
                           as.matrix(tmax_vals),
                           pops$Latitude)

diff_vals <- as.matrix(pr_vals) - pet_vals
pops$Im_annual <- rowSums(diff_vals)
pops$Im_summer <- rowSums(diff_vals[, 7:9])

pops <- cbind(pops, bio_vals)
env_cols <- c(paste0("bio", 1:19), "soil_pH", "Im_annual", "Im_summer")


# One PCA per species on all populations of that species.
pops$PC_climate1 <- NA_real_
pops$PC_climate2 <- NA_real_
pca_results <- list()

for (sp in sort(unique(pops$Species))) {
  rows    <- pops$Species == sp
  env_mat <- as.matrix(pops[rows, env_cols])

  v    <- apply(env_mat, 2, function(x) stats::var(x, na.rm = TRUE))
  keep <- v > 0
  if (!all(keep)) {
    warning(sprintf("[%s] dropping zero-variance columns: %s",
                    sp, paste(env_cols[!keep], collapse = ", ")))
    env_mat <- env_mat[, keep, drop = FALSE]
  }

  pca_sp <- prcomp(env_mat, center = TRUE, scale. = TRUE)
  pops$PC_climate1[rows] <- pca_sp$x[, 1]
  pops$PC_climate2[rows] <- pca_sp$x[, 2]
  pca_results[[sp]] <- pca_sp

  imp <- summary(pca_sp)$importance
  k   <- min(5, ncol(pca_sp$x))
  cat(sprintf("\nPCA-climate %s (n = %d):\n", sp, sum(rows)))
  print(round(imp[, seq_len(k)], 3))
  cat(sprintf("  PC1+PC2 cumulative: %.1f%%\n", 100 * imp[3, 2]))
}


# Save and merge back per-individual.
write.csv(pops, "pops_with_climate_soilpH_Im.csv", row.names = FALSE)
saveRDS(pca_results, file = "pca_climate_per_species.rds")

new_cols <- setdiff(names(pops),
                    c("Species", "PopulationCode", "Latitude", "Longitude"))
indiv <- indiv[, setdiff(names(indiv), new_cols), drop = FALSE]

indiv_out <- merge(indiv,
                   pops[, c("PopulationCode", new_cols)],
                   by = "PopulationCode", all.x = TRUE, sort = FALSE)
indiv_out <- indiv_out[match(indiv$IndividualCode, indiv_out$IndividualCode), ]

backup <- sub("\\.csv$", "_backup.csv", pop_file)
if (!file.exists(backup)) file.copy(pop_file, backup)
write.csv(indiv_out, pop_file, row.names = FALSE)

cat("\nWrote", length(new_cols), "columns to", pop_file, ":\n  ",
    paste(new_cols, collapse = ", "), "\n")
