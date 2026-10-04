# SPH5105 Assessment 2. Kriging a case surface for Singapore
# runs top to bottom in a clean session from any working directory: it
# creates data/ and output/ there, downloads the boundary, stations and bus
# stops on the first run, and writes every figure and table to output/
#
# AI declaration: an AI coding assistant (Claude) was used to help draft
# and check this script; every line has been reviewed and can be explained
# by the author
#
# data sources
#   case table: the 30-site table printed in the SPH5105 Assessment 2
#     brief (course material, accessed 2026-10-05), copied from the PDF
#     and embedded below; latitude/longitude are checked against x/y
#   Singapore boundary: GADM 4.1 level 0, via geodata
#     https://gadm.org; downloaded 2026-10-05; free for academic use,
#     redistribution not permitted
#   MRT and LRT station exits: Land Transport Authority, via data.gov.sg
#     https://data.gov.sg/datasets/d_b39d3a0871985372d7e1637193335da5/view
#     (data as of August 2025); downloaded 2026-10-05; Singapore Open Data
#     Licence
#   bus stops: Land Transport Authority, via data.gov.sg
#     https://data.gov.sg/datasets/d_3f172c6feb3f4f92a2f47d93eed2908a/view
#     (data as of October 2025); downloaded 2026-10-05; Singapore Open
#     Data Licence

dir.create("data/raw",     recursive = TRUE, showWarnings = FALSE)
dir.create("data/derived", recursive = TRUE, showWarnings = FALSE)
dir.create("output",       recursive = TRUE, showWarnings = FALSE)

packages <- c("sf", "sp", "gstat", "ggplot2", "patchwork", "ggspatial",
              "ggrepel", "spdep", "MASS", "glmmTMB", "jsonlite", "geodata")

to_install <- setdiff(packages, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)

invisible(lapply(packages, library, character.only = TRUE))

set.seed(2026)   # the permutation envelope and the conditional simulations

# every number quoted in the report is collected here and written out at
# the end, so none of them lives only on the console
key_stats <- data.frame(section = character(0), quantity = character(0),
                        value = numeric(0))
add_key <- function(section, quantity, value) {
  key_stats <<- rbind(key_stats,
                      data.frame(section = section, quantity = quantity,
                                 value = signif(as.numeric(value), 4)))
}


# 1. Data preparation and exploration

# the 30-site table from the brief
d <- read.csv(text = "
id,latitude,longitude,x,y,green_space,cases
1,1.35898725,103.80157451,24467.2538,37895.4449,73.5,0
2,1.35919306,103.74172365,17806.5442,37918.3715,73.5,0
3,1.30511305,103.72893254,16382.7875,31938.5296,53.2,0
4,1.34024942,103.80250224,24570.4732,35823.5097,49.6,3
5,1.34136456,103.76005553,19846.6124,35946.9159,51.0,0
6,1.33784567,103.75741086,19552.2773,35557.8244,58.0,0
7,1.35183838,103.76884059,20824.3262,37105.0315,68.8,0
8,1.34141926,103.73409254,16957.2094,35953.0666,71.9,1
9,1.30963239,103.85268296,30155.0733,32438.0223,34.1,4
10,1.2957689,103.83863871,28592.0836,30905.0612,43.7,4
11,1.28075393,103.81096333,25512.042,29244.7929,41.6,0
12,1.28766264,103.85870801,30825.6299,30008.7263,27.1,5
13,1.3117603,103.82594231,27179.0915,32673.3081,38.5,0
14,1.28081029,103.86964069,32042.3559,29251.0435,26.8,2
15,1.28285181,103.86518384,31546.3417,29476.7767,34.1,5
16,1.28999662,103.86824188,31886.6661,30266.818,35.5,6
17,1.36281129,103.89011933,34321.2549,38318.3381,52.5,3
18,1.31854781,103.88275186,33501.4408,33423.8895,26.3,2
19,1.35702013,103.85987322,30955.2262,37677.9234,39.0,2
20,1.35668903,103.86796024,31855.2188,37641.3232,34.6,3
21,1.34803813,103.91398564,36977.351,36684.8726,32.7,1
22,1.37425128,103.9108013,36622.875,39583.3794,43.9,1
23,1.3733223,103.86022265,30994.093,39480.5349,59.3,2
24,1.38411242,103.84399115,29187.7212,40673.6372,58.7,2
25,1.3591962,103.74305249,17954.4289,37918.7131,74.0,3
26,1.4025087,103.90285385,35738.3296,42707.9168,50.7,2
27,1.28869037,103.77990973,22056.0469,30122.4155,47.3,2
28,1.37747661,103.82379617,26940.2743,39939.8831,76.1,1
29,1.32844527,103.80695934,25066.486,34518.2607,49.6,0
30,1.32206898,103.86472879,31495.6439,33813.21,41.7,2
")
stopifnot(nrow(d) == 30, !anyNA(d))
str(d)

# the brief gives every site twice. Transforming latitude/longitude to
# EPSG:3414 must reproduce x/y before either pair is trusted
ll_m   <- sf::st_as_sf(d, coords = c("longitude", "latitude"), crs = 4326) |>
  sf::st_transform(3414) |>
  sf::st_coordinates()
xy_dev <- sqrt((ll_m[, 1] - d$x)^2 + (ll_m[, 2] - d$y)^2)
stopifnot(max(xy_dev) < 1)
add_key("1 data", "largest lat/long vs x/y disagreement (m)", max(xy_dev))

# the spatial object, on the projected pair: SVY21 / Singapore TM, metres
pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 3414, remove = FALSE)
sf::st_crs(pts)$units_gdal               # "metre"

# the boundary, the rail stations and the bus stops, each downloaded once
sg_raw <- geodata::gadm("SGP", level = 0, path = "data/raw")
if (is.null(sg_raw)) stop("GADM server unreachable; re-run later")
sg <- sf::st_as_sf(sg_raw) |>
  sf::st_transform(3414) |>
  sf::st_make_valid()

datagov <- function(id, file) {
  if (!file.exists(file)) {
    poll <- jsonlite::fromJSON(paste0(
      "https://api-open.data.gov.sg/v1/public/api/datasets/", id,
      "/poll-download"))
    download.file(poll$data$url, file, mode = "wb", quiet = TRUE)
  }
  sf::st_read(file, quiet = TRUE) |>
    sf::st_transform(3414)
}

exits <- datagov("d_b39d3a0871985372d7e1637193335da5",
                 "data/raw/lta_mrt_station_exits.geojson")
bus   <- datagov("d_3f172c6feb3f4f92a2f47d93eed2908a",
                 "data/raw/lta_bus_stops.geojson")
nrow(exits); nrow(bus)

# one point per station, the mean position of its exits; exits and bus
# stops within 400 m measure how reachable each station is
ex_name <- tools::toTitleCase(tolower(exits$STATION_NA))
ex_name <- sub("Lrt Station$", "LRT Station",
               sub("Mrt Station$", "MRT Station", ex_name))
mrt <- aggregate(as.data.frame(sf::st_coordinates(exits)),
                 list(name = ex_name), mean) |>
  sf::st_as_sf(coords = c("X", "Y"), crs = 3414)
mrt$n_exits  <- as.integer(table(ex_name)[mrt$name])
mrt$bus_400m <- lengths(sf::st_is_within_distance(mrt, bus, dist = 400))
mrt <- mrt[lengths(sf::st_intersects(mrt, sg)) > 0, ]
nrow(mrt)

# the scale of the problem, which sets every distance choice below
D <- as.matrix(dist(d[, c("x", "y")]))
diag(D) <- NA
hull <- sf::st_convex_hull(sf::st_union(pts))

scale_tab <- data.frame(
  quantity = c("east-west extent (km)", "north-south extent (km)",
               "largest separation (km)", "median nearest neighbour (km)",
               "smallest separation (m)", "convex hull area (km2)",
               "Singapore land area (GADM, km2)",
               "share of Singapore inside the hull (%)"),
  value = c(diff(range(d$x)) / 1000, diff(range(d$y)) / 1000,
            max(D, na.rm = TRUE) / 1000,
            median(apply(D, 1, min, na.rm = TRUE)) / 1000,
            min(D, na.rm = TRUE),
            as.numeric(sf::st_area(hull)) / 1e6,
            as.numeric(sf::st_area(sg)) / 1e6,
            100 * as.numeric(sf::st_area(hull) / sf::st_area(sg))))
scale_tab$value <- round(scale_tab$value, 2)
scale_tab
write.csv(scale_tab, "output/table1a_spatial_scale.csv", row.names = FALSE)

# assumptions, stated before modelling: distances are planar metres in
# EPSG:3414; each count is a value observed at a point (point support),
# not an area rate; there is no population denominator, so the surface is
# of recorded cases, not of risk per person; all sites share one time
# window; green_space is a percentage within an unstated buffer

map_furniture <- list(
  ggspatial::annotation_scale(location = "bl", width_hint = 0.2),
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    height = unit(0.8, "cm"), width = unit(0.8, "cm"),
    style = ggspatial::north_arrow_minimal()),
  theme_void(base_size = 11),
  theme(plot.title = element_text(face = "bold"),
        plot.caption = element_text(size = 7, colour = "grey45"),
        plot.background = element_rect(fill = "white", colour = NA)))

# zero counts are drawn as crosses; a circle of area zero would vanish
p1 <- ggplot() +
  geom_sf(data = sg, fill = "grey97", colour = "grey60", linewidth = 0.3) +
  geom_sf(data = hull, fill = NA, colour = "grey40", linetype = 2,
          linewidth = 0.3) +
  geom_sf(data = pts[pts$cases == 0, ], aes(shape = "no case"),
          size = 2, colour = "grey25", stroke = 0.7) +
  geom_sf(data = pts[pts$cases > 0, ], aes(size = cases), shape = 21,
          fill = "#B2182B", colour = "white", stroke = 0.4, alpha = 0.9) +
  ggrepel::geom_text_repel(data = pts, aes(label = id, geometry = geometry),
                           stat = "sf_coordinates", size = 2.4,
                           colour = "grey30", min.segment.length = 0,
                           segment.size = 0.2, max.overlaps = 50) +
  scale_size_area(max_size = 7, breaks = 1:6, name = "cases") +
  scale_shape_manual(values = c("no case" = 4), name = NULL) +
  labs(title = "Thirty case locations, Singapore",
       subtitle = "Symbol area proportional to the count; crosses are sites with no case; dashed line is the convex hull of the sites",
       caption = "EPSG:3414 (SVY21 / Singapore TM). Boundary: GADM 4.1.") +
  map_furniture

p1
ggsave("output/fig1_case_locations.png", p1, width = 9, height = 6,
       dpi = 300)

# green space: basic statistics, overall and by whether any case was seen
gs <- d$green_space
gsum <- function(z, group) {
  data.frame(group = group, n = length(z), mean = mean(z), sd = sd(z),
             median = median(z), q1 = unname(quantile(z, 0.25)),
             q3 = unname(quantile(z, 0.75)), min = min(z), max = max(z))
}
green_summary <- rbind(gsum(gs, "all sites"),
                       gsum(gs[d$cases == 0], "sites with no case"),
                       gsum(gs[d$cases > 0],  "sites with one or more"))
green_summary[, -1] <- round(green_summary[, -1], 1)
green_summary
write.csv(green_summary, "output/table1b_green_space_summary.csv",
          row.names = FALSE)

p2 <- ggplot(d, aes(green_space)) +
  geom_histogram(binwidth = 5, boundary = 25, fill = "#2E7D32",
                 colour = "white") +
  geom_vline(xintercept = median(gs), linetype = 2) +
  scale_x_continuous(breaks = seq(20, 80, 10)) +
  labs(x = "green space (%)", y = "sites",
       title = "Green space at the 30 sites",
       subtitle = sprintf("mean %.1f, median %.1f (dashed), range %.1f to %.1f",
                          mean(gs), median(gs), min(gs), max(gs))) +
  theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

p2
ggsave("output/fig2_green_space.png", p2, width = 7, height = 4, dpi = 300)


# 2. Variogram construction and fitting

# global spatial autocorrelation first: is there anything to model?
# k nearest neighbours is asymmetric, so the weights are symmetrised for
# the analytical Moran test
lw_pts <- spdep::listw2U(spdep::nb2listw(spdep::knn2nb(spdep::knearneigh(
  sf::st_coordinates(pts), k = 4)), style = "W"))
moran_cases <- spdep::moran.test(pts$cases, lw_pts)
moran_cases
add_key("2 variogram", "Moran's I of cases, k = 4", moran_cases$estimate[1])
add_key("2 variogram", "Moran's I p-value", moran_cases$p.value)

# bins: the cutoff is about half the largest separation, the usual limit
# beyond which too few pairs span the area; 1.5 km bins leave every bin
# but the first with 28 or more pairs
v <- gstat::variogram(cases ~ 1, pts, width = 1500, cutoff = 10000)
v

# a permutation envelope: the same locations with the counts shuffled.
# The band is pointwise, so with seven bins a single bin outside it is
# weak evidence on its own
env <- replicate(199, {
  q <- pts
  q$cases <- sample(q$cases)
  gstat::variogram(cases ~ 1, q, width = 1500, cutoff = 10000)$gamma
})
v$env_lo <- apply(env, 1, quantile, 0.025)
v$env_hi <- apply(env, 1, quantile, 0.975)
v$outside_envelope <- v$gamma < v$env_lo | v$gamma > v$env_hi
v

# a fit that fails to converge still returns parameters, so every fit
# records its warnings
fit_check <- function(v, model, ...) {
  warn <- character(0)
  f <- withCallingHandlers(
    gstat::fit.variogram(v, model, ...),
    warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  attr(f, "converged") <- length(warn) == 0 && !attr(f, "singular")
  f
}

loocv <- function(model, formula = cases ~ 1, data = pts) {
  cv <- gstat::krige.cv(formula, data, model = model,
                        nfold = nrow(data), verbose = FALSE)
  c(rmse = sqrt(mean(cv$residual^2)), mae = mean(abs(cv$residual)),
    z_mean = mean(cv$zscore), z_sd = sd(cv$zscore))
}

# a linear model has a slope, not a sill or a range, and says so
describe_fit <- function(f, label) {
  nug <- sum(f$psill[f$model == "Nug"])
  lin <- any(f$model == "Lin")
  ps  <- if (lin) NA else sum(f$psill[f$model != "Nug"])
  data.frame(model = label, converged = attr(f, "converged"),
             nugget = nug, partial_sill = ps,
             sill = if (lin) NA else nug + ps,
             range_m = if (lin) NA else max(f$range),
             slope_per_km = if (lin) 1000 * f$psill[f$model == "Lin"] else NA,
             t(loocv(f)))
}

# step one: the automatic fits, with the range free
free_fits <- list(
  "Sph, free range, default weights"   = fit_check(v, gstat::vgm(1.5, "Sph", 5000, 1.5)),
  "Exp, free range, default weights"   = fit_check(v, gstat::vgm(1.5, "Exp", 5000, 1.5)),
  "Gau, free range, default weights"   = fit_check(v, gstat::vgm(1.5, "Gau", 5000, 1.5)),
  "Sph, free range, weighted by pairs" = fit_check(v, gstat::vgm(1.5, "Sph", 5000, 1.5),
                                                   fit.method = 1),
  "Exp, free range, weighted by pairs" = fit_check(v, gstat::vgm(1.5, "Exp", 5000, 1.5),
                                                   fit.method = 1))
free_tab <- do.call(rbind, Map(describe_fit, free_fits, names(free_fits)))
free_tab

# none converges: each range runs past the 22 km span of the data, because
# the variogram keeps rising inside the cutoff. The range is not
# identified by 30 points, so it is profiled: held at each value on a
# grid, the nugget and partial sill become a linear least squares problem,
# which rarely fails although it can still be singular.
# Weighting by pairs keeps the 10-pair first bin from dominating, as it
# does under the default weights of pairs over distance squared
ranges  <- seq(3000, 15000, by = 500)
profile <- do.call(rbind, lapply(c("Sph", "Exp"), function(mo) {
  do.call(rbind, lapply(ranges, function(r) {
    f <- fit_check(v, gstat::vgm(1.5, mo, r, 1.5), fit.ranges = FALSE,
                   fit.method = 1)
    data.frame(model = mo, range_m = r, converged = attr(f, "converged"),
               nugget = f$psill[1], partial_sill = f$psill[2],
               at_boundary = f$psill[2] == 0,
               t(loocv(f)))
  }))
}))

# a fixed-range fit is linear in the sills, so it rarely fails, though it
# can be singular; gstat clamps a negative partial sill to zero without
# warning, so a pure nugget at a short range is a boundary solution,
# flagged above. The range itself is
# chosen by cross validation, then checked against the likelihood below

# the profile's choice: the converged fit with the lowest cross validated
# error. Selecting on the same LOOCV makes its RMSE slightly optimistic
best <- profile[profile$converged, ][which.min(
  profile$rmse[profile$converged]), ]
best

fv <- fit_check(v, gstat::vgm(1.5, best$model, best$range_m, 1.5),
                fit.ranges = FALSE, fit.method = 1)
stopifnot(attr(fv, "converged"))
fv

# the competitors, on the same footing
pe <- profile[profile$model == "Exp" & profile$converged, ]
pe <- pe[which.min(pe$rmse), ]
fv_exp <- fit_check(v, gstat::vgm(1.5, "Exp", pe$range_m, 1.5),
                    fit.ranges = FALSE, fit.method = 1)
fv_lin <- fit_check(v, gstat::vgm(0.0002, "Lin", 0, 1.5), fit.method = 1)
# the no-structure control: a nugget equal to the sample variance, so
# prediction is the leave-one-out mean and the z-scores are honest. It is
# never fitted, so the attributes fit_check would set are set by hand
fv_nug <- gstat::vgm(psill = var(d$cases), model = "Nug", range = 0)
attr(fv_nug, "singular")  <- FALSE
attr(fv_nug, "converged") <- TRUE

# the control that matters most: a plane in x and y with independent
# errors, which is universal kriging with a pure nugget
cv_plane <- gstat::krige.cv(cases ~ x + y, pts,
                            model = gstat::vgm(psill = var(d$cases),
                                               model = "Nug", range = 0),
                            nfold = nrow(pts), verbose = FALSE)
plane_row <- data.frame(
  model = "Plane trend in x and y, no spatial correlation",
  converged = TRUE, nugget = NA, partial_sill = NA, sill = NA,
  range_m = NA, slope_per_km = NA,
  rmse = sqrt(mean(cv_plane$residual^2)),
  mae = mean(abs(cv_plane$residual)),
  z_mean = mean(cv_plane$zscore), z_sd = sd(cv_plane$zscore))

vario_tab <- rbind(
  free_tab,
  describe_fit(fv,     sprintf("%s, range profiled (selected)", best$model)),
  describe_fit(fv_exp, "Exp, range profiled"),
  describe_fit(fv_lin, "Linear, no sill (intrinsic)"),
  describe_fit(fv_nug, "Pure nugget, no spatial structure"),
  plane_row)
vario_tab[, -(1:2)] <- round(vario_tab[, -(1:2)], 3)
vario_tab
write.csv(vario_tab, "output/table2_variogram_models.csv", row.names = FALSE)
write.csv(profile,   "output/table2b_range_profile.csv",   row.names = FALSE)

# an independent check on the range: the Gaussian likelihood, constant
# mean, the selected correlation family plus a nugget. The mean and variance are
# concentrated out and the nugget share optimised at each range
# the correlation function follows the family selected above
sph_cor <- function(h, r) ifelse(h < r, 1 - 1.5 * h / r + 0.5 * (h / r)^3, 0)
exp_cor <- function(h, r) exp(-h / r)
cor_fun <- if (best$model == "Exp") exp_cor else sph_cor
Dm <- as.matrix(dist(d[, c("x", "y")]))
prof_ll <- function(r) {
  S <- cor_fun(Dm, r)
  y <- d$cases
  n <- length(y)
  nll <- function(tau, reml) {
    R  <- (1 - tau) * S + tau * diag(n)
    Ri <- solve(R)
    b  <- sum(Ri %*% y) / sum(Ri)
    q  <- as.numeric(t(y - b) %*% Ri %*% (y - b))
    ld <- as.numeric(determinant(R)$modulus)
    m  <- if (reml) n - 1 else n
    0.5 * (m * log(2 * pi * q / m) + ld + m + if (reml) log(sum(Ri)) else 0)
  }
  ml   <- optimize(nll, c(0.001, 0.999), reml = FALSE)
  reml <- optimize(nll, c(0.001, 0.999), reml = TRUE)
  c(range_m = r, ml_loglik = -ml$objective, ml_nugget_share = ml$minimum,
    reml_loglik = -reml$objective, reml_nugget_share = reml$minimum)
}
ll_prof <- as.data.frame(t(sapply(seq(2000, 30000, by = 500), prof_ll)))
write.csv(round(ll_prof, 4), "output/table2c_likelihood_profile.csv",
          row.names = FALSE)

ml_best <- ll_prof[which.max(ll_prof$ml_loglik), ]
ml_ci   <- range(ll_prof$range_m[ll_prof$ml_loglik >=
                                   max(ll_prof$ml_loglik) - qchisq(0.95, 1) / 2])
add_key("2 variogram", "range at maximum likelihood (m)", ml_best$range_m)
# REML is flat with two near-equal peaks, so report how far the selected
# range sits below its maximum rather than a single REML range
add_key("2 variogram", "range at maximum REML (m)",
        ll_prof$range_m[which.max(ll_prof$reml_loglik)])
add_key("2 variogram", "REML log likelihood at the selected range, below its maximum",
        ll_prof$reml_loglik[ll_prof$range_m == fv$range[2]] -
          max(ll_prof$reml_loglik))
add_key("2 variogram", "ML nugget share at its maximum", ml_best$ml_nugget_share)
add_key("2 variogram", "approximate 95% range interval, lower (m)", ml_ci[1])
add_key("2 variogram",
        if (ml_ci[2] == max(ll_prof$range_m))
          "approximate 95% range interval, upper (m; the grid edge, so unbounded)"
        else "approximate 95% range interval, upper (m)", ml_ci[2])

# the trend behind the runaway free fits: a plane in the coordinates is
# significant, and its residuals have a flat variogram
m_plane <- lm(cases ~ I(x / 1000) + I(y / 1000), data = d)
add_key("2 variogram", "trend in x, p", summary(m_plane)$coefficients[2, 4])
add_key("2 variogram", "trend in y, p", summary(m_plane)$coefficients[3, 4])
v_dt <- gstat::variogram(cases ~ x + y, pts, width = 1500, cutoff = 10000)
write.csv(data.frame(dist_m = round(v_dt$dist), np = v_dt$np,
                     gamma = round(v_dt$gamma, 3)),
          "output/table2d_detrended_variogram.csv", row.names = FALSE)

# the first bin is the poorest fit, driven by close pairs with very
# different counts, such as sites 2 and 25, 148 m apart with 0 and 3
add_key("2 variogram", "first bin semivariance, observed", v$gamma[1])
add_key("2 variogram", "first bin semivariance, fitted",
        gstat::variogramLine(fv, dist_vector = v$dist[1])$gamma)

# selecting the range on the same LOOCV that is reported is optimistic.
# Nested LOOCV repeats the whole selection inside every fold without the
# held out site: both families, the same 0.5 km range grid, converged fits
# only. As in the main analysis, the inner LOOCV holds the training
# variogram fixed. It is slow, so the result is kept in data/derived
nested_file <- "data/derived/nested_loocv_full.rds"
if (file.exists(nested_file)) {
  nested <- readRDS(nested_file)
} else {
  nested <- do.call(rbind, lapply(seq_len(nrow(pts)), function(i) {
    tr <- pts[-i, ]
    vi <- gstat::variogram(cases ~ 1, tr, width = 1500, cutoff = 10000)
    fit_r <- function(mo, r) fit_check(vi, gstat::vgm(1.5, mo, r, 1.5),
                                       fit.ranges = FALSE, fit.method = 1)
    gi <- expand.grid(model = c("Sph", "Exp"), range_m = ranges,
                      stringsAsFactors = FALSE)
    gi$rmse <- mapply(function(mo, r) {
      f <- fit_r(mo, r)
      if (attr(f, "converged")) loocv(f, data = tr)["rmse"] else NA_real_
    }, gi$model, gi$range_m)
    b_i <- gi[which.min(gi$rmse), ]
    k_i <- gstat::krige(cases ~ 1, tr, pts[i, ],
                        model = fit_r(b_i$model, b_i$range_m),
                        debug.level = 0)
    data.frame(site = i, model = b_i$model, range_m = b_i$range_m,
               observed = pts$cases[i], pred = k_i$var1.pred,
               var = k_i$var1.var)
  }))
  saveRDS(nested, nested_file)
}
add_key("2 variogram", "nested LOOCV RMSE (honest)",
        sqrt(mean((nested$observed - nested$pred)^2)))
add_key("2 variogram", "nested LOOCV z-score sd",
        sd((nested$observed - nested$pred) / sqrt(nested$var)))
add_key("2 variogram", "folds re-selecting the selected family",
        sum(nested$model == best$model))
add_key("2 variogram", "folds re-selecting a range within 1 km of the selected",
        sum(abs(nested$range_m - fv$range[2]) <= 1000))

# interpretation of the selected parameters
add_key("2 variogram", "nugget (cases squared)",       fv$psill[1])
add_key("2 variogram", "partial sill (cases squared)", fv$psill[2])
add_key("2 variogram", "sill (cases squared)",         sum(fv$psill))
add_key("2 variogram", "range (m)",                    fv$range[2])
add_key("2 variogram", "nugget share of the sill",     fv$psill[1] / sum(fv$psill))
add_key("2 variogram", "sample variance of cases",     var(d$cases))

# the mean ordinary kriging estimates: a generalised least squares mean that
# down-weights clustered sites, and the level extrapolation drifts towards
C_obs <- sum(fv$psill) -
  matrix(gstat::variogramLine(fv, dist_vector = as.vector(as.matrix(
    dist(d[, c("x", "y")]))))$gamma, nrow(d))
diag(C_obs) <- sum(fv$psill)
Ci <- solve(C_obs)
add_key("2 variogram", "mean estimated by ordinary kriging (GLS)",
        sum(Ci %*% d$cases) / sum(Ci))
add_key("2 variogram", "variance of that estimated mean", 1 / sum(Ci))
add_key("2 variogram", "bins outside the permutation envelope",
        sum(v$outside_envelope))

mod_cols <- c(selected = "#B2182B", exponential = "#2166AC",
              linear = "grey50")
curve_of <- function(f, label) {
  cl <- gstat::variogramLine(f, maxdist = 10500, n = 200)
  cl$model <- label
  cl
}
lab_sel <- sprintf("%s, range %.1f km (selected)", best$model,
                   fv$range[2] / 1000)
lab_exp <- sprintf("Exp, range %.1f km", fv_exp$range[2] / 1000)
curves  <- rbind(curve_of(fv, lab_sel), curve_of(fv_exp, lab_exp),
                 curve_of(fv_lin, "Linear, no sill"))
curves$model <- factor(curves$model,
                       levels = c(lab_sel, lab_exp, "Linear, no sill"))

p3a <- ggplot(v, aes(dist / 1000, gamma)) +
  geom_ribbon(aes(ymin = env_lo, ymax = env_hi), fill = "grey88") +
  geom_hline(yintercept = var(d$cases), linetype = 3, colour = "grey40") +
  geom_line(data = curves, aes(dist / 1000, gamma, colour = model,
                               linewidth = model)) +
  geom_point(data = v_dt, aes(dist / 1000, gamma), shape = 2,
             colour = "grey45", size = 1.8) +
  geom_point(aes(size = np), colour = "grey15") +
  geom_text(aes(label = np), nudge_y = 0.3, size = 2.8) +
  annotate("text", x = 0.2, y = fv$psill[1] - 0.25, hjust = 0, size = 3,
           colour = mod_cols[["selected"]],
           label = sprintf("nugget %.2f", fv$psill[1])) +
  annotate("text", x = fv$range[2] / 1000 + 0.2, y = sum(fv$psill) + 0.55,
           hjust = 0, size = 3, colour = mod_cols[["selected"]],
           label = sprintf("sill %.2f at range %.1f km", sum(fv$psill),
                           fv$range[2] / 1000)) +
  scale_size_area(max_size = 4, guide = "none") +
  scale_colour_manual(values = unname(mod_cols)) +
  scale_linewidth_manual(values = c(1.2, 0.6, 0.6), guide = "none") +
  labs(x = "lag distance (km)", y = "semivariance (cases squared)",
       colour = NULL,
       title = "Empirical and fitted variograms",
       subtitle = paste0("Grey band: 95% permutation envelope; dotted: sample variance\n",
                         "Labels: pairs per bin; triangles: after removing an x + y trend")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        plot.background = element_rect(fill = "white", colour = NA))

p3b <- ggplot(profile, aes(range_m / 1000, rmse, colour = model)) +
  geom_line() +
  geom_point(size = 1) +
  geom_point(data = best, size = 3, shape = 21, fill = "white",
             stroke = 1.2) +
  geom_hline(yintercept = loocv(fv_nug)["rmse"], linetype = 2) +
  annotate("text", x = 9, y = loocv(fv_nug)["rmse"] - 0.02,
           label = "pure nugget control", size = 3) +
  scale_colour_manual(values = c(Exp = mod_cols[["exponential"]],
                                 Sph = mod_cols[["selected"]])) +
  labs(x = "range held fixed (km)", y = "LOOCV RMSE (cases)",
       colour = NULL, title = "Profiling the range") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        plot.background = element_rect(fill = "white", colour = NA))

# each curve relative to its own maximum: ML and REML carry different
# constants, so only their shapes are comparable
p3c <- ggplot(ll_prof, aes(range_m / 1000)) +
  geom_line(aes(y = ml_loglik - max(ml_loglik)),
            colour = mod_cols[["selected"]]) +
  geom_line(aes(y = reml_loglik - max(reml_loglik)), colour = "grey50",
            linetype = 2) +
  geom_hline(yintercept = -qchisq(0.95, 1) / 2, linetype = 3) +
  geom_vline(xintercept = fv$range[2] / 1000, linetype = 3,
             colour = mod_cols[["selected"]]) +
  labs(x = "range (km)", y = "log likelihood below its maximum",
       title = "Likelihood check",
       subtitle = "Solid ML, dashed REML; dotted: approximate 95% bound") +
  theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

p3 <- p3a + p3b + p3c + patchwork::plot_layout(widths = c(1.6, 1, 1))
p3
ggsave("output/fig3_variogram.png", p3, width = 15, height = 5.2, dpi = 300)


# 3. Kriging interpolation

# a 1 km grid over Singapore: every cell that touches land. Keeping only
# cells whose centre is on land drops coastal cells, among them the cells
# of the CBD sites on the GADM shoreline, which are the hottest in the data
cells <- sf::st_make_grid(sg, cellsize = 1000, square = TRUE) |>
  sf::st_as_sf()
keep  <- lengths(sf::st_intersects(cells, sg)) > 0
cells <- cells[keep, ]
cent  <- sf::st_centroid(sf::st_geometry(cells)) |> sf::st_as_sf()
nrow(cells)                              # 1 km2 cells
add_key("3 kriging", "1 km2 grid cells covering Singapore", nrow(cells))

k_ok <- gstat::krige(cases ~ 1, pts, cent, model = fv)

cells$pred <- k_ok$var1.pred
cells$var  <- k_ok$var1.var
cells$se   <- sqrt(k_ok$var1.var)

# where the prediction is supported by data. A cell is supported when its
# centre, where it is predicted, lies inside the sites' convex hull, or
# when it contains a site; the second keeps the cells of the sites on the
# hull's corners, whose centres fall just outside it
cent_in    <- lengths(sf::st_intersects(cent, hull)) > 0
cell_site  <- lengths(sf::st_intersects(cells, pts)) > 0
cells$in_hull <- cent_in | cell_site
add_key("3 kriging", "supported cells, centre inside the hull", sum(cent_in))
add_key("3 kriging", "supported cells added because they contain a site",
        sum(cell_site & !cent_in))
cells$dist_nn_km <- apply(as.matrix(sf::st_distance(cent, pts)), 1,
                          min) / 1000
table(cells$in_hull)

# high risk zones: the top fifth of the supported surface
q_brk <- quantile(cells$pred[cells$in_hull], seq(0, 1, 0.2))
q80   <- unname(q_brk[5])
cells$quintile <- NA_integer_
cells$quintile[cells$in_hull] <- as.integer(cut(
  cells$pred[cells$in_hull], q_brk, include.lowest = TRUE, labels = FALSE))
cells$high <- cells$in_hull & cells$pred >= q80

# the probability that a cell exceeds the threshold, across conditional
# simulations: the kriged map is smoothed, so a tail quantity is taken
# from realisations that keep the variability of the data. These are
# Gaussian simulations of the point-support process, so the result is a
# model-based exceedance probability, conditional on the variogram and
# the threshold, not the probability of a true outbreak
# nmax = 30 is the neighbourhood of sequential simulation: the 30 nearest
# of the data and the nodes already simulated, not simply all 30 sites
sim <- gstat::krige(cases ~ 1, pts, cent, model = fv, nsim = 500,
                    nmax = 30, debug.level = 0)
sim_m <- as.matrix(sf::st_drop_geometry(sim))
cells$p_high <- rowMeans(sim_m >= q80)
add_key("3 kriging", "sd of one conditional simulation (cases)",
        sd(sim_m[, 1]))
add_key("3 kriging", "sd of the kriged surface within the sites (cases)",
        sd(cells$pred[cells$in_hull]))
add_key("3 kriging", "sd of the observed cases", sd(d$cases))
add_key("3 kriging", "exceedance probability, median over non-high cells",
        median(cells$p_high[cells$in_hull & !cells$high]))
add_key("3 kriging", "exceedance probability, max over non-high cells",
        max(cells$p_high[cells$in_hull & !cells$high]))

zones <- cells[cells$high, ] |>
  sf::st_union() |>
  sf::st_cast("POLYGON") |>
  sf::st_as_sf()
zones$zone     <- paste0("Z", seq_len(nrow(zones)))
zones$area_km2 <- round(as.numeric(sf::st_area(zones)) / 1e6, 1)
zones
zone_lab <- sf::st_point_on_surface(zones)

# sanity checks on the surface
sanity <- data.frame(
  check = c("observed cases, min", "observed cases, max",
            "predicted, min (all cells)", "predicted, max (all cells)",
            "predicted, min (within the sites)",
            "predicted, max (within the sites)",
            "cells predicted below 0",
            "high risk threshold, 80th percentile within the sites",
            "mean prediction where extrapolated",
            "sample mean of the observed cases",
            "Spearman: standard error against distance to nearest site",
            "the same, within the sites only",
            "kriging standard error, min", "kriging standard error, max",
            "square root of the nugget", "square root of the sill",
            "high risk cells with exceedance probability above 0.5"),
  value = c(min(d$cases), max(d$cases),
            min(cells$pred), max(cells$pred),
            min(cells$pred[cells$in_hull]), max(cells$pred[cells$in_hull]),
            sum(cells$pred < 0), q80,
            mean(cells$pred[!cells$in_hull]), mean(d$cases),
            cor(cells$se, cells$dist_nn_km, method = "spearman"),
            cor(cells$se[cells$in_hull], cells$dist_nn_km[cells$in_hull],
                method = "spearman"),
            min(cells$se), max(cells$se),
            sqrt(fv$psill[1]), sqrt(sum(fv$psill)),
            mean(cells$p_high[cells$high] > 0.5)))
sanity$value <- round(sanity$value, 3)
sanity
write.csv(sanity, "output/table3_sanity_checks.csv", row.names = FALSE)

cv_final <- gstat::krige.cv(cases ~ 1, pts, model = fv, nfold = nrow(pts),
                            verbose = FALSE)
add_key("3 kriging", "LOOCV RMSE, selected model", sqrt(mean(cv_final$residual^2)))
add_key("3 kriging", "LOOCV RMSE, pure nugget", loocv(fv_nug)["rmse"])
add_key("3 kriging", "LOOCV z-score sd", sd(cv_final$zscore))
add_key("3 kriging", "LOOCV correlation, observed and predicted",
        cor(cv_final$observed, cv_final$var1.pred))
add_key("3 kriging", "LOOCV prediction at the site with 6 cases",
        cv_final$var1.pred[which.max(cv_final$observed)])

p_risk <- ggplot() +
  geom_sf(data = cells, aes(fill = pred), colour = NA) +
  geom_sf(data = sg, fill = NA, colour = "grey45", linewidth = 0.25) +
  geom_sf(data = hull, fill = NA, colour = "grey15", linetype = 2,
          linewidth = 0.4) +
  geom_sf(data = zones, fill = NA, colour = "black", linewidth = 0.7) +
  geom_sf(data = pts, shape = 21, size = 1.6, fill = "white",
          colour = "black", stroke = 0.4) +
  geom_sf_label(data = zone_lab, aes(label = zone), size = 3,
                linewidth = 0.2, label.padding = unit(0.12, "lines")) +
  scale_fill_viridis_c(option = "inferno", direction = -1,
                       name = "predicted\ncases per site") +
  labs(title = "(a) Ordinary kriging: predicted cases on a 1 km grid across Singapore",
       subtitle = "Outside the dashed convex hull the surface is extrapolation and drifts towards the mean") +
  map_furniture

q_lab <- sprintf("%d: %.2f to %.2f", 1:5, q_brk[1:5], q_brk[2:6])

p_quint <- ggplot() +
  geom_sf(data = cells[!cells$in_hull, ], fill = "grey92", colour = NA) +
  geom_sf(data = cells[cells$in_hull, ], aes(fill = factor(quintile)),
          colour = NA) +
  geom_sf(data = sg, fill = NA, colour = "grey45", linewidth = 0.25) +
  geom_sf(data = zones, fill = NA, colour = "black", linewidth = 0.7) +
  geom_sf(data = pts, shape = 21, size = 1.6, fill = "white",
          colour = "black", stroke = 0.4) +
  geom_sf_label(data = zone_lab, aes(label = zone), size = 3,
                linewidth = 0.2, label.padding = unit(0.12, "lines")) +
  scale_fill_brewer(palette = "YlOrRd", name = "risk quintile\n(predicted cases)",
                    labels = q_lab) +
  labs(title = "(b) High risk zones: the top quintile within the sites",
       subtitle = "Grey cells are extrapolation and are not classed") +
  map_furniture

p4 <- (p_risk / p_quint) +
  patchwork::plot_annotation(
    caption = sprintf("%s variogram: nugget %.2f, partial sill %.2f, range %.1f km. EPSG:3414. Boundary: GADM 4.1.",
                      best$model, fv$psill[1], fv$psill[2],
                      fv$range[2] / 1000),
    theme = theme(plot.caption = element_text(size = 8, colour = "grey40"),
                  plot.background = element_rect(fill = "white",
                                                 colour = NA)))
p4
ggsave("output/fig4_kriged_risk.png", p4, width = 10, height = 12,
       dpi = 300)

p_se <- ggplot() +
  geom_sf(data = cells, aes(fill = se), colour = NA) +
  geom_sf(data = sg, fill = NA, colour = "grey45", linewidth = 0.25) +
  geom_sf(data = hull, fill = NA, colour = "white", linetype = 2,
          linewidth = 0.3) +
  geom_sf(data = pts, size = 0.9, colour = "white") +
  scale_fill_viridis_c(option = "mako",
                       name = "kriging\nstandard\nerror (cases)") +
  labs(title = "(a) Kriging standard error",
       subtitle = sprintf("Rising with distance from the sites; the scale spans only %.2f to %.2f cases,\nas the nugget alone contributes %.2f even beside a site",
                          min(cells$se), max(cells$se), sqrt(fv$psill[1]))) +
  map_furniture

p_exc <- ggplot() +
  geom_sf(data = cells[cells$in_hull, ], aes(fill = p_high), colour = NA) +
  geom_sf(data = cells[!cells$in_hull, ], fill = "grey92", colour = NA) +
  geom_sf(data = sg, fill = NA, colour = "grey45", linewidth = 0.25) +
  geom_sf(data = zones, fill = NA, colour = "black", linewidth = 0.7) +
  geom_sf(data = pts, size = 0.9) +
  scale_fill_distiller(palette = "RdPu", direction = 1, limits = c(0, 1),
                       name = "P(above the\nhigh risk\nthreshold)") +
  labs(title = "(b) Probability a cell exceeds the high risk threshold",
       subtitle = sprintf("Model-based: share of 500 Gaussian conditional simulations above %.2f",
                          q80)) +
  map_furniture

p5 <- p_se / p_exc
p5
ggsave("output/fig5_kriging_uncertainty.png", p5, width = 10, height = 11,
       dpi = 300)

p6a <- ggplot(sf::st_drop_geometry(cells),
              aes(dist_nn_km, se,
                  colour = ifelse(in_hull, "within the sites",
                                  "extrapolated"))) +
  geom_hline(yintercept = sqrt(fv$psill[1]), linetype = 3) +
  geom_hline(yintercept = sqrt(sum(fv$psill)), linetype = 2) +
  geom_point(size = 0.8, alpha = 0.7) +
  annotate("text", x = 12, y = sqrt(fv$psill[1]) + 0.02, size = 3,
           label = "square root of the nugget") +
  annotate("text", x = 12, y = sqrt(sum(fv$psill)) + 0.02, size = 3,
           label = "square root of the sill") +
  scale_colour_manual(values = c("within the sites" = "#2166AC",
                                 "extrapolated" = "grey60"), name = NULL) +
  labs(x = "distance to nearest site (km)",
       y = "kriging standard error (cases)",
       title = "Uncertainty rises away from the data") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        plot.background = element_rect(fill = "white", colour = NA))

p6b <- ggplot(as.data.frame(cv_final), aes(observed, var1.pred)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(colour = "#B2182B", alpha = 0.7) +
  coord_equal(xlim = c(-0.5, 6.5), ylim = c(-0.5, 6.5)) +
  labs(x = "observed cases", y = "LOOCV prediction (cases)",
       title = "Leave one out cross validation",
       subtitle = sprintf("RMSE %.2f, z-score sd %.2f",
                          sqrt(mean(cv_final$residual^2)),
                          sd(cv_final$zscore))) +
  theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

p6 <- p6a + p6b
p6
ggsave("output/fig6_sanity_checks.png", p6, width = 11, height = 5,
       dpi = 300)


# 4. Spatial correlation with green space

# rank correlation at the sites
sp_cor <- cor.test(d$cases, d$green_space, method = "spearman",
                   exact = FALSE)
sp_cor
add_key("4 green space", "Spearman rho, cases and green space", sp_cor$estimate)
add_key("4 green space", "Spearman p-value (naive: sites treated as independent)",
        sp_cor$p.value)

# the primary model: Poisson, the natural model for a count
m_pois <- glm(cases ~ green_space, family = poisson, data = d)
dispersion <- sum(residuals(m_pois, type = "pearson")^2) /
  m_pois$df.residual
add_key("4 green space", "Poisson dispersion (Pearson)", dispersion)

# sensitivity to overdispersion: quasi-Poisson rescales the errors, the
# negative binomial estimates theta; the LR test is on the boundary, so
# its p-value is halved
m_quasi <- glm(cases ~ green_space, family = quasipoisson, data = d)
m_nb    <- MASS::glm.nb(cases ~ green_space, data = d)
nb_lr_p <- 0.5 * pchisq(2 * as.numeric(logLik(m_nb) - logLik(m_pois)),
                        df = 1, lower.tail = FALSE)
add_key("4 green space", "negative binomial theta", m_nb$theta)
add_key("4 green space", "negative binomial theta, SE", m_nb$SE.theta)
add_key("4 green space", "LR test NB against Poisson, p", nb_lr_p)

# the confounding check. Green space falls towards the city centre, and so
# do population and footfall, which the data cannot measure. Distance to
# Raffles Place MRT stands in for urban centrality; the reference point is
# fixed before the models are fitted
cbd_ref <- mrt[mrt$name == "Raffles Place MRT Station", ]
stopifnot(nrow(cbd_ref) == 1)
d$dcbd_km <- as.numeric(sf::st_distance(pts, cbd_ref)) / 1000
cent_cor <- cor(d$green_space, d$dcbd_km, method = "spearman")
add_key("4 green space", "Spearman rho, green space and distance to CBD",
        cent_cor)

m_dcbd <- glm(cases ~ dcbd_km, family = poisson, data = d)
m_gd   <- glm(cases ~ green_space + dcbd_km, family = poisson, data = d)
m_xy   <- glm(cases ~ I(x / 1000) + I(y / 1000), family = poisson, data = d)
m_gxy  <- glm(cases ~ green_space + I(x / 1000) + I(y / 1000),
              family = poisson, data = d)

# a spatial Poisson model: an exponential random field soaks up the
# spatial structure, leaving green space to explain what remains
d$pos <- glmmTMB::numFactor(d$x / 1000, d$y / 1000)
d$grp <- factor(1)
m_sp  <- glmmTMB::glmmTMB(cases ~ green_space + exp(pos + 0 | grp),
                          family = poisson, data = d)

# with 30 sites this is exploratory; its fit is checked before use
add_key("4 green space", "GLMM optimiser converged (0 = yes)",
        m_sp$fit$convergence)
add_key("4 green space", "GLMM Hessian positive definite (1 = yes)",
        m_sp$sdr$pdHess)

# without a population offset the exponentiated coefficient is a ratio of
# expected counts, not an incidence rate ratio. Quasi-Poisson uses t-based
# inference on its residual degrees of freedom. Residual Moran p-values
# come from fitted residuals and are read as approximate
green_row <- function(m, label, role) {
  quasi <- inherits(m, "glm") && grepl("quasi", family(m)$family)
  crit  <- if (quasi) qt(0.975, m$df.residual) else qnorm(0.975)
  if (inherits(m, "glmmTMB")) {
    b  <- glmmTMB::fixef(m)$cond["green_space"]
    se <- sqrt(vcov(m)$cond["green_space", "green_space"])
    res_I <- NA; res_p <- NA
  } else {
    has_g <- "green_space" %in% names(coef(m))
    b  <- if (has_g) coef(m)["green_space"] else NA
    se <- if (has_g) sqrt(vcov(m)["green_space", "green_space"]) else NA
    mt <- spdep::moran.test(residuals(m, type = "pearson"), lw_pts)
    res_I <- mt$estimate[1]; res_p <- mt$p.value
  }
  data.frame(model = label, role = role,
             count_ratio_per_10pp = exp(10 * b),
             lo95 = exp(10 * (b - crit * se)),
             hi95 = exp(10 * (b + crit * se)),
             p = if (quasi) 2 * pt(-abs(b / se), m$df.residual) else
               2 * pnorm(-abs(b / se)),
             AIC = tryCatch(AIC(m), error = function(e) NA),
             residual_moran_I = res_I, residual_moran_p = res_p,
             row.names = NULL)
}

green_tab <- rbind(
  green_row(m_pois,  "cases ~ green (Poisson)",        "primary"),
  green_row(m_quasi, "cases ~ green (quasi-Poisson)",  "overdispersion"),
  green_row(m_nb,    "cases ~ green (negative binomial)", "overdispersion"),
  green_row(m_dcbd,  "cases ~ distance to CBD",        "confounder alone"),
  green_row(m_gd,    "cases ~ green + distance to CBD", "confounding"),
  green_row(m_xy,    "cases ~ x + y trend",            "confounder alone"),
  green_row(m_gxy,   "cases ~ green + x + y trend",    "confounding"),
  green_row(m_sp,    "cases ~ green + spatial field (GLMM)", "spatial"))
green_tab[, -(1:2)] <- signif(green_tab[, -(1:2)], 3)
green_tab
write.csv(green_tab, "output/table4_green_space_models.csv",
          row.names = FALSE)

# the residual Moran's I falls after any smooth regional covariate, not
# only after green space, so that fall cannot single green space out

# a green space surface for the overlay. Its variogram climbs from 16 to
# over 200 with no plateau, so bounded models chase an ever longer range;
# green space is a regional gradient, which the linear model describes
v_g  <- gstat::variogram(green_space ~ 1, pts, width = 1500, cutoff = 10000)
v_g
fv_g <- fit_check(v_g, gstat::vgm(0.02, "Lin", 0, 50), fit.method = 1)
stopifnot(attr(fv_g, "converged"))
fv_g
add_key("4 green space", "green space LOOCV RMSE, linear variogram (pp)",
        loocv(fv_g, green_space ~ 1)["rmse"])

cells$green <- gstat::krige(green_space ~ 1, pts, cent,
                            model = fv_g)$var1.pred

# descriptive only: both surfaces are smoothed from the same 30 points
add_key("4 green space", "cell-level Spearman, risk and green (descriptive)",
        cor(cells$pred[cells$in_hull], cells$green[cells$in_hull],
            method = "spearman"))

nd <- data.frame(green_space = seq(25, 78, length.out = 100))
pr <- predict(m_pois, nd, type = "link", se.fit = TRUE)
nd$fit <- exp(pr$fit)
nd$lo  <- exp(pr$fit - 1.96 * pr$se.fit)
nd$hi  <- exp(pr$fit + 1.96 * pr$se.fit)

p7a <- ggplot() +
  geom_sf(data = cells[cells$in_hull, ], aes(fill = green), colour = NA) +
  geom_sf(data = cells[!cells$in_hull, ], fill = "grey92", colour = NA) +
  geom_sf(data = sg, fill = NA, colour = "grey45", linewidth = 0.25) +
  geom_sf(data = zones, fill = NA, colour = "#B2182B", linewidth = 0.8) +
  geom_sf(data = pts[pts$cases == 0, ], shape = 4, size = 1.8,
          colour = "black", stroke = 0.6) +
  geom_sf(data = pts[pts$cases > 0, ], aes(size = cases), shape = 21,
          fill = "white", colour = "black", stroke = 0.3) +
  geom_sf_label(data = zone_lab, aes(label = zone), size = 3,
                linewidth = 0.2, label.padding = unit(0.12, "lines")) +
  scale_fill_gradient(low = "#F7FCF5", high = "#00441B",
                      name = "green\nspace (%)") +
  scale_size_area(max_size = 5, breaks = 1:6, name = "cases") +
  labs(title = "(a) Green space against the high risk zones",
       subtitle = "Kriged green space within the sites; red outlines are the high risk zones; crosses had no case") +
  map_furniture

p7b <- ggplot() +
  geom_ribbon(data = nd, aes(green_space, ymin = lo, ymax = hi),
              fill = "grey85") +
  geom_line(data = nd, aes(green_space, fit), colour = "#2E7D32",
            linewidth = 0.9) +
  geom_point(data = d, aes(green_space, cases), size = 2) +
  labs(x = "green space (%)", y = "cases",
       title = "(b) Poisson regression",
       subtitle = sprintf("Count ratio per 10 points of green space %.2f (%.2f to %.2f)",
                          green_tab$count_ratio_per_10pp[1], green_tab$lo95[1],
                          green_tab$hi95[1])) +
  theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

p7 <- p7a + p7b + patchwork::plot_layout(widths = c(1.7, 1))
p7
ggsave("output/fig7_green_space.png", p7, width = 14, height = 5.5,
       dpi = 300)


# 5. Outbreak profiling and intervention targeting

# outbreak size. Kriging predicts counts on the support of the data: the
# expected count recorded at a site, not cases per km2, and summing it
# over cells would invent an outbreak total. So zones are described by
# their block-kriged mean per site with its standard error, and cases per
# km2 is the recorded cases over the zone's area, beside the site density
# that it depends on
krige_block <- function(poly) {
  gstat::krige(cases ~ 1, as(pts, "Spatial"),
               as(sf::st_as_sf(poly), "Spatial"), model = fv)
}
zb <- krige_block(zones)
hb <- krige_block(sf::st_sf(geometry = hull))

site_zone <- sf::st_join(pts, zones, join = sf::st_intersects)

# cells touching the coast include sea, so a zone's footprint and its land
# area differ; cases per km2 uses land
sg_land <- sf::st_union(sf::st_geometry(sg))
land_km2 <- function(g) as.numeric(sf::st_area(
  sf::st_intersection(sf::st_geometry(g), sg_land))) / 1e6
zone_land <- vapply(seq_len(nrow(zones)), function(i)
  sum(land_km2(zones[i, ])), numeric(1))

zone_tab <- data.frame(
  zone = c(zones$zone, "all sites (hull)"),
  footprint_km2 = c(zones$area_km2,
                    round(as.numeric(sf::st_area(hull)) / 1e6, 1)),
  land_km2 = c(zone_land, sum(land_km2(hull))),
  sites = c(vapply(zones$zone, function(z) sum(site_zone$zone %in% z),
                   integer(1)), nrow(pts)),
  observed_cases = c(vapply(zones$zone, function(z)
    sum(site_zone$cases[site_zone$zone %in% z]), numeric(1)), sum(d$cases)),
  block_kriged_mean_per_site = c(zb$var1.pred, hb$var1.pred),
  block_kriged_se = sqrt(c(zb$var1.var, hb$var1.var)))
zone_tab$recorded_cases_per_land_km2 <- zone_tab$observed_cases /
  zone_tab$land_km2
zone_tab$cases_per_site <- zone_tab$observed_cases / zone_tab$sites
zone_tab$excess_over_mean_site <- zone_tab$cases_per_site / mean(d$cases)
zone_tab$sites_per_land_km2 <- zone_tab$sites / zone_tab$land_km2
zone_tab[, -1] <- round(zone_tab[, -1], 2)
zone_tab
write.csv(zone_tab, "output/table5_outbreak_zones.csv", row.names = FALSE)

# siting. Each high risk cell carries its predicted cases as the weight,
# and a station reaches the cells within its catchment, 2 km being about a
# 25 minute walk or a short feeder bus ride
hcent <- cent[cells$high, ]
w     <- cells$pred[cells$high]
hzone <- sf::st_join(hcent, zones, join = sf::st_intersects)$zone

reach_matrix <- function(catch_m) {
  sf::st_is_within_distance(hcent, mrt, dist = catch_m, sparse = FALSE)
}
cover <- reach_matrix(2000)

# one centre per zone. The two zones are separate foci, so each gets one.
# Stations within a pragmatic 10 per cent of the best reach are treated as
# equivalent, and among those the most accessible is chosen: most bus
# stops within 400 m, then most exits
zone_candidates <- function(z, cover) {
  s <- colSums(cover[hzone == z, , drop = FALSE] * w[hzone == z])
  data.frame(zone = z, station = mrt$name, reach = s,
             share_of_best = s / max(s), bus_400m = mrt$bus_400m,
             exits = mrt$n_exits, idx = seq_len(nrow(mrt)))
}
pick_station <- function(cand, tol = 0.10) {
  ok <- cand[cand$share_of_best >= 1 - tol, ]
  ok[order(-ok$bus_400m, -ok$exits, -ok$reach), ][1, ]
}

# the proposal assumes exactly two zones and two distinct stations, and
# stops rather than continue silently if either fails
stopifnot(nrow(zones) == 2)
cand_all <- lapply(zones$zone, zone_candidates, cover = cover)
picked   <- do.call(rbind, lapply(cand_all, pick_station))
stopifnot(length(unique(picked$idx)) == 2)
picked

# the 10 per cent tolerance is a pragmatic choice, not one calibrated to
# the prediction error, so the choice is shown under other tolerances
tol_tab <- do.call(rbind, lapply(c(0, 0.05, 0.10, 0.15), function(tl) {
  pk <- do.call(rbind, lapply(cand_all, pick_station, tol = tl))
  data.frame(tolerance = tl, Z1_station = pk$station[1],
             Z2_station = pk$station[2],
             share_of_high_risk_surface = round(sum(w[
               cover[, pk$idx[1]] | cover[, pk$idx[2]]]) / sum(w), 3))
}))
tol_tab
write.csv(tol_tab, "output/table6d_tolerance_sensitivity.csv",
          row.names = FALSE)

# the top candidates in each zone, so the choice is visible
top_cand <- do.call(rbind, lapply(cand_all, function(cd) {
  cd <- cd[order(-cd$reach), ]
  cd$rank <- seq_len(nrow(cd))
  own <- picked$idx[picked$zone == cd$zone[1]]       # this zone's choice
  cd <- cd[cd$rank <= 6 | cd$idx == own, ]
  cd$chosen <- cd$idx == own
  cd[, c("zone", "rank", "station", "reach", "share_of_best", "bus_400m",
         "exits", "chosen")]
}))
top_cand$reach <- round(top_cand$reach, 2)
top_cand$share_of_best <- round(top_cand$share_of_best, 3)
top_cand
write.csv(top_cand, "output/table6c_candidate_stations.csv",
          row.names = FALSE)

# the alternatives the proposal is compared with
pure_cov <- vapply(cand_all, function(cd) cd$idx[which.max(cd$reach)],
                   integer(1))
pairs <- utils::combn(nrow(mrt), 2)
score <- apply(pairs, 2, function(ij)
  sum(w[cover[, ij[1]] | cover[, ij[2]]]))
max_cover <- pairs[, which.max(score)]

reach_of <- function(idx, cover) {
  cov <- cover[, idx[1]] | cover[, idx[2]]
  data.frame(stations = paste(mrt$name[idx], collapse = " + "),
             high_risk_km2_reached = sum(cov),
             share_of_high_risk_surface = round(sum(w[cov]) / sum(w), 3),
             zones_reached = paste(sort(unique(hzone[cov])),
                                   collapse = ", "),
             km_apart = round(as.numeric(
               sf::st_distance(mrt[idx[1], ], mrt[idx[2], ])) / 1000, 1))
}
cover1 <- reach_matrix(1000)
strategy_tab <- rbind(
  cbind(strategy = "per zone, accessible within 10% (proposed), 2 km",
        reach_of(picked$idx, cover)),
  cbind(strategy = "per zone, coverage only, 2 km",
        reach_of(pure_cov, cover)),
  cbind(strategy = "maximum total coverage, 2 km",
        reach_of(max_cover, cover)),
  cbind(strategy = "per zone, accessible within 10% (proposed), 1 km",
        reach_of(picked$idx, cover1)))
strategy_tab
write.csv(strategy_tab, "output/table6b_siting_strategies.csv",
          row.names = FALSE)

centres <- mrt[picked$idx, ]
centres$label <- paste("Centre", LETTERS[seq_len(nrow(centres))])
centres$serves_zone <- picked$zone

centre_tab <- data.frame(
  centre = centres$label, station = centres$name,
  serves_zone = centres$serves_zone,
  bus_stops_400m = centres$bus_400m, exits = centres$n_exits,
  x = round(sf::st_coordinates(centres)[, 1]),
  y = round(sf::st_coordinates(centres)[, 2]),
  zone_km2_within_2km = vapply(seq_len(nrow(picked)), function(i)
    sum(cover[hzone == picked$zone[i], picked$idx[i]]), integer(1)),
  zone_area_km2 = zones$area_km2[match(picked$zone, zones$zone)],
  dist_to_zone_km = round(vapply(seq_len(nrow(picked)), function(i)
    as.numeric(sf::st_distance(centres[i, ],
                               zones[zones$zone == picked$zone[i], ])),
    numeric(1)) / 1000, 2))
centre_tab
write.csv(centre_tab, "output/table6_intervention_centres.csv",
          row.names = FALSE)

# labels sit off the zones, joined to their centre by a leader line
cxy <- sf::st_coordinates(centres)
lab_pos <- data.frame(
  x = cxy[, 1], y = cxy[, 2],
  lx = cxy[, 1] + 9000,
  ly = cxy[, 2] + ifelse(cxy[, 2] < mean(cxy[, 2]), -3000, 2500),
  text = paste0(centres$label, ": ",
                sub(" (MRT|LRT) Station$", "", centres$name),
                "\nserves ", centres$serves_zone))

p8 <- ggplot() +
  geom_sf(data = cells[cells$in_hull, ], aes(fill = pred), colour = NA) +
  geom_sf(data = sg, fill = NA, colour = "grey45", linewidth = 0.25) +
  geom_sf(data = zones, fill = NA, colour = "black", linewidth = 0.7) +
  geom_sf(data = mrt, size = 0.4, colour = "grey30", alpha = 0.5) +
  geom_sf(data = sf::st_buffer(centres, 2000), fill = NA,
          colour = "#2166AC", linewidth = 0.7, linetype = 2) +
  geom_sf(data = centres, shape = 23, size = 4, fill = "#2166AC",
          colour = "white") +
  geom_segment(data = lab_pos,
               aes(x = x, y = y, xend = x + 0.6 * (lx - x),
                   yend = y + 0.6 * (ly - y)),
               colour = "#2166AC", linewidth = 0.4) +
  geom_label(data = lab_pos, aes(lx, ly, label = text), size = 2.8,
             linewidth = 0.2) +
  geom_sf_label(data = zone_lab, aes(label = zone), size = 3,
                linewidth = 0.2, label.padding = unit(0.12, "lines"),
                nudge_x = 2200, nudge_y = 1200) +
  scale_fill_viridis_c(option = "inferno", direction = -1,
                       name = "predicted\ncases per site") +
  labs(title = "Two temporary intervention centres, one per high risk zone",
       subtitle = "Black: high risk zones. Diamonds: chosen stations with dashed 2 km catchments. Grey dots: other MRT and LRT stations.",
       caption = "Stations and bus stops: LTA via data.gov.sg, Open Data Licence. Surface shown within the sites' convex hull only.") +
  map_furniture

p8
ggsave("output/fig8_intervention_centres.png", p8, width = 10, height = 6.5,
       dpi = 300)


# 6. Robustness, for the interpretation

# do the high risk zones survive other defensible choices? Overlap is the
# Jaccard index of high risk cells against the main analysis
top20 <- function(pred) {
  qq <- quantile(pred[cells$in_hull], 0.80)
  cells$in_hull & pred >= qq
}
jaccard <- function(a, b) sum(a & b) / sum(a | b)

pred_exp <- gstat::krige(cases ~ 1, pts, cent, model = fv_exp)$var1.pred
pred_lin <- gstat::krige(cases ~ 1, pts, cent, model = fv_lin)$var1.pred

# the Anscombe transform stabilises the variance of counts with means
# near 2; its range is profiled on a wider grid, since it runs long
d$ans   <- 2 * sqrt(d$cases + 3 / 8)
pts$ans <- d$ans
v_s  <- gstat::variogram(ans ~ 1, pts, width = 1500, cutoff = 10000)
ranges_wide <- seq(3000, 25000, by = 500)
pr_s <- do.call(rbind, lapply(ranges_wide, function(r) {
  f <- fit_check(v_s, gstat::vgm(0.4, "Sph", r, 0.4), fit.ranges = FALSE,
                 fit.method = 1)
  data.frame(range_m = r, converged = attr(f, "converged"),
             rmse = loocv(f, ans ~ 1)["rmse"])
}))
bs   <- pr_s[pr_s$converged, ][which.min(pr_s$rmse[pr_s$converged]), ]
fv_s <- fit_check(v_s, gstat::vgm(0.4, "Sph", bs$range_m, 0.4),
                  fit.ranges = FALSE, fit.method = 1)
add_key("6 robustness", "Anscombe-scale range selected (m)", bs$range_m)
add_key("6 robustness", "Anscombe range at the edge of its grid",
        bs$range_m %in% range(ranges_wide))
pred_sqrt <- gstat::krige(ans ~ 1, pts, cent, model = fv_s)$var1.pred

# the range is the least identified parameter, so the zones are re-made
# with the sills refitted at other fixed ranges
pred_at_range <- function(r) {
  f <- fit_check(v, gstat::vgm(1.5, "Sph", r, 1.5), fit.ranges = FALSE,
                 fit.method = 1)
  gstat::krige(cases ~ 1, pts, cent, model = f, debug.level = 0)$var1.pred
}
range_alt <- c(4500, 8000, 10000, 15000)
pred_r <- lapply(range_alt, pred_at_range)

# what overlap chance alone gives: two random sets of the same size
n_hi <- sum(cells$high)
chance <- replicate(2000, {
  a <- sample(sum(cells$in_hull), n_hi)
  b <- sample(sum(cells$in_hull), n_hi)
  length(intersect(a, b)) / length(union(a, b))
})

q75 <- cells$in_hull & cells$pred >= quantile(cells$pred[cells$in_hull], 0.75)
robust <- data.frame(
  variant = c("main: Sph profiled, top fifth",
              "Exp variogram, top fifth",
              "linear variogram, top fifth",
              "Anscombe-transformed counts, top fifth",
              sprintf("Sph, range fixed at %.1f km, top fifth",
                      range_alt / 1000),
              "main surface, top quarter instead",
              "chance: two random sets, mean",
              "chance: two random sets, 97.5th percentile"),
  high_risk_km2 = c(sum(cells$high), sum(top20(pred_exp)),
                    sum(top20(pred_lin)), sum(top20(pred_sqrt)),
                    vapply(pred_r, function(p) sum(top20(p)), numeric(1)),
                    sum(q75), n_hi, n_hi),
  jaccard_with_main = c(1, jaccard(cells$high, top20(pred_exp)),
                        jaccard(cells$high, top20(pred_lin)),
                        jaccard(cells$high, top20(pred_sqrt)),
                        vapply(pred_r, function(p)
                          jaccard(cells$high, top20(p)), numeric(1)),
                        jaccard(cells$high, q75),
                        mean(chance), unname(quantile(chance, 0.975))))
robust$jaccard_with_main <- round(robust$jaccard_with_main, 2)
robust
write.csv(robust, "output/table7_robustness.csv", row.names = FALSE)

# block support: the mean over each 1 km2 cell rather than its centre
k_blk <- gstat::krige(cases ~ 1, pts, cent, model = fv, block = c(1000, 1000))
add_key("6 robustness", "mean kriging variance, point support", mean(cells$var))
add_key("6 robustness", "mean kriging variance, 1 km2 block support",
        mean(k_blk$var1.var))

key_stats
write.csv(key_stats, "output/table8_key_statistics.csv", row.names = FALSE)

# the kriged grid, kept for the report
sf::st_write(cells, "data/derived/kriged_grid.gpkg", delete_dsn = TRUE,
             quiet = TRUE)

sessionInfo()
