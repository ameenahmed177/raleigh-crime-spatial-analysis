# ST 533 Project 3
# Part 1A: CSR Assessment (LARCENY in Raleigh)

library(maps)
library(splancs)
library(sf)

# Load crime data
crime = read.csv("WakeCrime.csv")

# Load township shapefile
townships = st_read("Townships-shp/townships.shp")

# Select Raleigh township
raleigh_poly = townships[townships$NAME == "RALEIGH", ]

# Convert polygon to splancs format
coords_poly = st_coordinates(raleigh_poly)[,1:2]
raleigh_sp = as.points(coords_poly)

# Filter crime data: Raleigh + LARCENY
raleigh = subset(crime, City == "RALEIGH")
larceny = subset(raleigh,
                 Crime_Category == "LARCENY" &
                 !is.na(Longitude) & !is.na(Latitude))

# Convert lon/lat points to township CRS for consistent spatial analysis
larceny_sf = st_as_sf(larceny,
                      coords = c("Longitude", "Latitude"),
                      crs = 4326,
                      remove = FALSE)
larceny_sf = st_transform(larceny_sf, st_crs(townships))

# Extract projected coordinates
coords_larceny = st_coordinates(larceny_sf)[, 1:2]

# Downsample to keep the point-pattern calculations fast enough for iteration,
# while still leaving enough events for an interesting pattern analysis.
max_points = 100
if(nrow(coords_larceny) > max_points){
        set.seed(1)
        coords_larceny = coords_larceny[sample(nrow(coords_larceny), max_points), ]
}

# Plot points
png("project3_larceny_map.png", width=1800, height=1800, res=300)
polymap(raleigh_sp)
pointmap(coords_larceny, add=TRUE)
dev.off()

# Number of points
nrow(coords_larceny)

# Choose spatial lags (adjust if needed)
h = seq(0, 2000, by=500)

# K-function
k_larceny = khat(as.points(coords_larceny), raleigh_sp, h)

# L-function
l_larceny = sqrt(k_larceny / pi) - h

# Plot L-function
png("project3_larceny_L.png", width=1800, height=1800, res=300)
plot(h, l_larceny, type="l",
     ylab="L", xlab="h",
     main="L-function: Larceny (Raleigh)")
abline(h=0, col=2)
dev.off()

# CSR envelope
# Keep the Monte Carlo check lightweight so the script finishes quickly.
nsim = 999
sim = matrix(0, nsim, length(h))

for(i in 1:nsim){
  sim[i,] = sqrt(khat(csr(raleigh_sp, nrow(coords_larceny)),
                       raleigh_sp, h) / pi) - h
}

# Envelope (min/max)
env = apply(sim, 2, range)

# Plot envelope
png("project3_larceny_envelope.png", width=1800, height=1800, res=300)
matplot(h, t(env), type="l", lty=2, col=3,
        ylab="L", xlab="h",
        main="CSR Envelope: Larceny (Raleigh)")
lines(h, l_larceny, lwd=2)
abline(h=0, col=2)
legend("topright",
       legend=c("Observed L", "Envelope", "CSR=0"),
       col=c("black","green","red"),
       lty=c(1,2,1),
       lwd=c(2,1,1))
dev.off()

# Part 1B: CSR Assessment (SIMPLE ASSAULT in Raleigh)

# Filter crime data: Raleigh + SIMPLE ASSAULT
simple_assault = subset(raleigh,
                        Crime_Category == "SIMPLE ASSAULT" &
                        !is.na(Longitude) & !is.na(Latitude))

# Convert lon/lat points to township CRS for consistent spatial analysis
simple_assault_sf = st_as_sf(simple_assault,
                             coords = c("Longitude", "Latitude"),
                             crs = 4326,
                             remove = FALSE)
simple_assault_sf = st_transform(simple_assault_sf, st_crs(townships))

# Extract projected coordinates
coords_simple_assault = st_coordinates(simple_assault_sf)[, 1:2]

# Downsample to keep the point-pattern calculations fast enough for iteration,
# while still leaving enough events for an interesting pattern analysis.
max_points = 100
if(nrow(coords_simple_assault) > max_points){
        set.seed(1)
        coords_simple_assault = coords_simple_assault[sample(nrow(coords_simple_assault), max_points), ]
}

# Plot points
png("project3_simple_assault_map.png", width=1800, height=1800, res=300)
polymap(raleigh_sp)
pointmap(coords_simple_assault, add=TRUE)
dev.off()

# Number of points
nrow(coords_simple_assault)

# Choose spatial lags (adjust if needed)
h = seq(0, 2000, by=500)

# K-function
k_simple_assault = khat(as.points(coords_simple_assault), raleigh_sp, h)

# L-function
l_simple_assault = sqrt(k_simple_assault / pi) - h

# Plot L-function
png("project3_simple_assault_L.png", width=1800, height=1800, res=300)
plot(h, l_simple_assault, type="l",
     ylab="L", xlab="h",
     main="L-function: Simple Assault (Raleigh)")
abline(h=0, col=2)
dev.off()

# CSR envelope
# Keep the Monte Carlo check lightweight so the script finishes quickly.
nsim = 999
sim = matrix(0, nsim, length(h))

for(i in 1:nsim){
  sim[i,] = sqrt(khat(csr(raleigh_sp, nrow(coords_simple_assault)),
                       raleigh_sp, h) / pi) - h
}

# Envelope (min/max)
env = apply(sim, 2, range)

# Plot envelope
png("project3_simple_assault_envelope.png", width=1800, height=1800, res=300)
matplot(h, t(env), type="l", lty=2, col=3,
        ylab="L", xlab="h",
        main="CSR Envelope: Simple Assault (Raleigh)")
lines(h, l_simple_assault, lwd=2)
abline(h=0, col=2)
legend("topright",
       legend=c("Observed L", "Envelope", "CSR=0"),
       col=c("black","green","red"),
       lty=c(1,2,1),
       lwd=c(2,1,1))
dev.off()

# Part 2: Comparing LARCENY and SIMPLE ASSAULT in Raleigh

library(MASS)

# Combine the two point patterns
all_coords = rbind(coords_larceny, coords_simple_assault)

# Labels: 1 = LARCENY, 0 = SIMPLE ASSAULT
n1 = nrow(coords_larceny)
n0 = nrow(coords_simple_assault)
labels_obs = c(rep(1, n1), rep(0, n0))

# Create grid over Raleigh township bounding box
bbox = st_bbox(raleigh_poly)
xg = seq(bbox["xmin"], bbox["xmax"], length.out = 100)
yg = seq(bbox["ymin"], bbox["ymax"], length.out = 100)
grid = expand.grid(x = xg, y = yg)

# Keep only grid points inside Raleigh polygon
grid_sf = st_as_sf(grid, coords = c("x", "y"), crs = st_crs(townships))
inside = st_within(grid_sf, raleigh_poly, sparse = FALSE)[,1]
grid_in = grid[inside, ]

# Bandwidth for kernel density
bandwidth = 1000

# Estimate intensity surface on grid
get_logrisk = function(coords1, coords0, grid_pts, bw){
  kd1 = kde2d(coords1[,1], coords1[,2],
              h = c(bw, bw),
              n = c(length(unique(grid_pts$x)), length(unique(grid_pts$y))),
              lims = c(range(grid_pts$x), range(grid_pts$y)))
  
  kd0 = kde2d(coords0[,1], coords0[,2],
              h = c(bw, bw),
              n = c(length(unique(grid_pts$x)), length(unique(grid_pts$y))),
              lims = c(range(grid_pts$x), range(grid_pts$y)))
  
  # Small constant avoids log(0)
  eps = 1e-8
  logrisk = log(kd1$z + eps) - log(kd0$z + eps) - log(nrow(coords1) / nrow(coords0))
  
  list(x = kd1$x, y = kd1$y, z = logrisk)
}

# Observed log relative risk
obs_lr = get_logrisk(coords_larceny, coords_simple_assault, grid_in, bandwidth)

# Plot contour map of observed log relative risk
png("project3_log_relative_risk_contour.png", width=1800, height=1800, res=300)
filled.contour(obs_lr$x, obs_lr$y, obs_lr$z,
               color.palette = terrain.colors,
               plot.title = title(main = "Log Relative Risk:\nLarceny vs Simple Assault",
                                  xlab = "x", ylab = "y"),
               plot.axes = {
                 axis(1)
                 axis(2)
                 contour(obs_lr$x, obs_lr$y, obs_lr$z, add = TRUE)
               })
dev.off()

# Monte Carlo relabeling
nsim_part2 = 199
sim_array = array(NA, dim = c(length(obs_lr$x), length(obs_lr$y), nsim_part2))

set.seed(1)

for(i in 1:nsim_part2){
  perm_labels = sample(labels_obs)
  
  coords1_sim = all_coords[perm_labels == 1, ]
  coords0_sim = all_coords[perm_labels == 0, ]
  
  sim_lr = get_logrisk(coords1_sim, coords0_sim, grid_in, bandwidth)
  sim_array[,,i] = sim_lr$z
}

# Pointwise 95% Monte Carlo intervals
lower = apply(sim_array, c(1,2), quantile, probs = 0.025, na.rm = TRUE)
upper = apply(sim_array, c(1,2), quantile, probs = 0.975, na.rm = TRUE)

# Classify significant regions
sig_map = matrix(0, nrow = nrow(obs_lr$z), ncol = ncol(obs_lr$z))
sig_map[obs_lr$z < lower] = -1
sig_map[obs_lr$z > upper] = 1

# Plot significance map
png("project3_log_relative_risk_significance.png", width=1800, height=1800, res=300)
image(obs_lr$x, obs_lr$y, sig_map,
      col = c("blue", "white", "red"),
      xlab = "x", ylab = "y",
      main = "Significant Regions for Log Relative Risk")
contour(obs_lr$x, obs_lr$y, sig_map, add = TRUE, drawlabels = FALSE)
legend("topright",
       legend = c("Simple Assault higher", "Not significant", "Larceny higher"),
       fill = c("blue", "white", "red"))
dev.off()