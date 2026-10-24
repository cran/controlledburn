## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6, fig.height = 5,
  out.width = "80%",
  dpi = 120
)
library(ggplot2)
library(controlledburn)

# Shared polygon for all illustrations
poly_x <- c(1.3, 6.7, 5.2, 2.8, 1.3)
poly_y <- c(1.5, 2.3, 5.8, 4.2, 1.5)
nx <- 8; ny <- 7

# Point-in-polygon via winding number
pip <- function(px, py, poly_x, poly_y) {
  n <- length(poly_x) - 1
  winding <- 0
  for (i in seq_len(n)) {
    y0 <- poly_y[i]; y1 <- poly_y[i + 1]
    x0 <- poly_x[i]; x1 <- poly_x[i + 1]
    if (y0 <= py) {
      if (y1 > py) {
        cross_x <- x0 + (py - y0) / (y1 - y0) * (x1 - x0)
        if (px < cross_x) winding <- winding + 1
      }
    } else {
      if (y1 <= py) {
        cross_x <- x0 + (py - y0) / (y1 - y0) * (x1 - x0)
        if (px < cross_x) winding <- winding - 1
      }
    }
  }
  winding != 0
}

# Edge-scanline intersections
edge_crossings <- function(poly_x, poly_y, y_mid) {
  n <- length(poly_x) - 1
  crossings <- numeric(0)
  for (i in seq_len(n)) {
    y0 <- poly_y[i]; y1 <- poly_y[i + 1]
    x0 <- poly_x[i]; x1 <- poly_x[i + 1]
    if ((y0 < y_mid & y1 >= y_mid) | (y1 < y_mid & y0 >= y_mid)) {
      t <- (y_mid - y0) / (y1 - y0)
      crossings <- c(crossings, x0 + t * (x1 - x0))
    }
  }
  sort(crossings)
}

# Base grid theme
theme_grid <- function() {
  theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank())
}

## ----scanline-fill, fig.cap = "Scanline fill: for each row, polygon edges are intersected with the row's horizontal scanline (dashed blue). Red dots mark the intersections. Between each pair of crossings, all cells along the blue segment are 'inside' the polygon."----
scanline_segs <- do.call(rbind, lapply(seq(0.5, ny - 0.5), function(y_mid) {
  xints <- edge_crossings(poly_x, poly_y, y_mid)
  if (length(xints) >= 2) {
    pairs <- matrix(xints, ncol = 2, byrow = TRUE)
    data.frame(y = y_mid, xstart = pairs[, 1], xend = pairs[, 2])
  }
}))

scanline_pts <- do.call(rbind, lapply(seq(0.5, ny - 0.5), function(y_mid) {
  xints <- edge_crossings(poly_x, poly_y, y_mid)
  if (length(xints) >= 2) data.frame(y = y_mid, x = xints)
}))

ggplot() +
  geom_vline(xintercept = 0:nx, colour = "grey85", linewidth = 0.3) +
  geom_hline(yintercept = 0:ny, colour = "grey85", linewidth = 0.3) +
  geom_hline(yintercept = seq(0.5, ny - 0.5), colour = "steelblue",
             linewidth = 0.2, linetype = "dashed") +
  geom_polygon(aes(x = poly_x, y = poly_y),
               fill = NA, colour = "black", linewidth = 0.8) +
  geom_segment(data = scanline_segs,
               aes(x = xstart, xend = xend, y = y, yend = y),
               colour = "steelblue", linewidth = 1.5, alpha = 0.5) +
  geom_point(data = scanline_pts, aes(x, y),
             colour = "firebrick", size = 2.5) +
  coord_equal(xlim = c(0, nx), ylim = c(0, ny)) +
  theme_grid()

## ----center-rule, fig.cap = "Cell-centre rule: each cell is classified by whether its centre point (dot) falls inside the polygon. Red = inside, grey = outside. No partial coverage — every cell is either fully in or fully out."----
centers <- expand.grid(
  cx = seq(0.5, nx - 0.5, by = 1),
  cy = seq(0.5, ny - 0.5, by = 1)
)
centers$inside <- mapply(pip, centers$cx, centers$cy,
                         MoreArgs = list(poly_x = poly_x, poly_y = poly_y))

# Highlight "in" cells
in_rects <- centers[centers$inside, ]
in_rects <- data.frame(
  xmin = in_rects$cx - 0.5, xmax = in_rects$cx + 0.5,
  ymin = in_rects$cy - 0.5, ymax = in_rects$cy + 0.5
)

ggplot() +
  geom_rect(data = in_rects,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "steelblue", alpha = 0.4) +
  geom_vline(xintercept = 0:nx, colour = "grey80", linewidth = 0.3) +
  geom_hline(yintercept = 0:ny, colour = "grey80", linewidth = 0.3) +
  geom_polygon(aes(x = poly_x, y = poly_y),
               fill = NA, colour = "black", linewidth = 0.8) +
  geom_point(data = centers, aes(cx, cy, colour = inside), size = 1.8) +
  scale_colour_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey60"),
                      guide = "none") +
  coord_equal(xlim = c(0, nx), ylim = c(0, ny)) +
  theme_grid()

## ----exact-coverage, fig.cap = "Exact coverage: interior cells (blue, fraction = 1.0) are fully inside. Boundary cells (coral) have an exact coverage fraction — the proportion of the cell area that falls inside the polygon edge. Total coverage equals the polygon's true area."----
# Proper segment-cell intersection via parametric clipping
seg_crosses_cell <- function(x0, y0, x1, y1, cxmin, cymin, cxmax, cymax) {
  if (min(x0, x1) > cxmax || max(x0, x1) < cxmin ||
      min(y0, y1) > cymax || max(y0, y1) < cymin) return(FALSE)
  sdx <- x1 - x0; sdy <- y1 - y0
  if (sdx != 0) {
    tx1 <- (cxmin - x0) / sdx; tx2 <- (cxmax - x0) / sdx
    if (tx1 > tx2) { tmp <- tx1; tx1 <- tx2; tx2 <- tmp }
  } else {
    tx1 <- if (x0 >= cxmin && x0 <= cxmax) -Inf else Inf
    tx2 <- if (x0 >= cxmin && x0 <= cxmax) Inf else -Inf
  }
  if (sdy != 0) {
    ty1 <- (cymin - y0) / sdy; ty2 <- (cymax - y0) / sdy
    if (ty1 > ty2) { tmp <- ty1; ty1 <- ty2; ty2 <- tmp }
  } else {
    ty1 <- if (y0 >= cymin && y0 <= cymax) -Inf else Inf
    ty2 <- if (y0 >= cymin && y0 <= cymax) Inf else -Inf
  }
  max(tx1, ty1, 0) < min(tx2, ty2, 1) - 1e-10
}

cell_class <- function(cx, cy) {
  xmin <- cx - 0.5; xmax <- cx + 0.5
  ymin <- cy - 0.5; ymax <- cy + 0.5
  n <- length(poly_x) - 1
  edge_crosses <- FALSE
  for (i in seq_len(n)) {
    if (seg_crosses_cell(poly_x[i], poly_y[i], poly_x[i+1], poly_y[i+1],
                         xmin, ymin, xmax, ymax)) {
      edge_crosses <- TRUE; break
    }
  }
  if (!edge_crosses) {
    return(if (pip(cx, cy, poly_x, poly_y)) "interior" else "exterior")
  }
  "boundary"
}

centers$class <- mapply(cell_class, centers$cx, centers$cy)

cell_rects <- data.frame(
  xmin = centers$cx - 0.5, xmax = centers$cx + 0.5,
  ymin = centers$cy - 0.5, ymax = centers$cy + 0.5,
  class = centers$class
)

ggplot() +
  geom_rect(data = cell_rects[cell_rects$class != "exterior", ],
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = class),
            alpha = 0.5, colour = "grey70", linewidth = 0.2) +
  geom_vline(xintercept = 0:nx, colour = "grey80", linewidth = 0.3) +
  geom_hline(yintercept = 0:ny, colour = "grey80", linewidth = 0.3) +
  geom_polygon(aes(x = poly_x, y = poly_y),
               fill = NA, colour = "black", linewidth = 0.8) +
  scale_fill_manual(values = c("interior" = "steelblue",
                                "boundary" = "coral"),
                    labels = c("interior" = "interior (fraction = 1.0)",
                               "boundary" = "boundary (0 < fraction < 1)")) +
  coord_equal(xlim = c(0, nx), ylim = c(0, ny)) +
  theme_grid()

## ----side-by-side, fig.width = 9, fig.height = 4.5, fig.cap = "The same polygon, two classification schemes. Left: cell-centre rule — every cell is in or out. Right: exact coverage — boundary cells (coral) carry a fractional weight. Interior cells (blue) are the same in both."----
# Side-by-side
centers$approx <- ifelse(centers$inside, "in", "out")
centers$coverage <- centers$class

d1 <- data.frame(
  xmin = centers$cx - 0.5, xmax = centers$cx + 0.5,
  ymin = centers$cy - 0.5, ymax = centers$cy + 0.5,
  fill = ifelse(centers$approx == "in", "approx: in", "none"),
  panel = "Cell-centre rule (approx)"
)
d2 <- data.frame(
  xmin = centers$cx - 0.5, xmax = centers$cx + 0.5,
  ymin = centers$cy - 0.5, ymax = centers$cy + 0.5,
  fill = ifelse(centers$coverage == "interior", "interior",
                ifelse(centers$coverage == "boundary", "boundary", "none")),
  panel = "Exact coverage (coverage)"
)

dd <- rbind(d1, d2)
dd <- dd[dd$fill != "none", ]

poly_df <- data.frame(x = poly_x, y = poly_y)
grid_v <- data.frame(x = 0:nx)
grid_h <- data.frame(y = 0:ny)

ggplot(dd) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = fill),
            alpha = 0.5, colour = "grey70", linewidth = 0.2) +
  geom_vline(data = grid_v, aes(xintercept = x),
             colour = "grey80", linewidth = 0.3) +
  geom_hline(data = grid_h, aes(yintercept = y),
             colour = "grey80", linewidth = 0.3) +
  geom_polygon(data = poly_df, aes(x, y),
               fill = NA, colour = "black", linewidth = 0.8) +
  scale_fill_manual(
    values = c("approx: in" = "steelblue",
               "interior" = "steelblue",
               "boundary" = "coral"),
    labels = c("approx: in" = "cell centre inside",
               "interior" = "interior (1.0)",
               "boundary" = "boundary (fractional)")
  ) +
  facet_wrap(~panel) +
  coord_equal(xlim = c(0, nx), ylim = c(0, ny)) +
  theme_grid() +
  theme(strip.text = element_text(face = "bold", size = 11))

## ----approx-detail, fig.cap = "Approx mode in detail: for each row, polygon edges are intersected with the scanline at the cell-centre y-coordinate. The winding number accumulates left-to-right: it goes nonzero when entering the polygon (green arrow) and returns to zero when leaving (red arrow). Cells whose centre falls in a nonzero-winding span are classified as 'inside' and emitted as runs."----
# Illustration of winding accumulation for one row
y_show <- 3.5
xints <- edge_crossings(poly_x, poly_y, y_show)

# Build winding annotation
winding_df <- data.frame(
  x = c(0, xints),
  xend = c(xints, nx),
  y = y_show, 
  winding = cumsum(c(0, rep(c(1, -1), length.out = length(xints))))
)

ggplot() +
  geom_vline(xintercept = 0:nx, colour = "grey85", linewidth = 0.3) +
  geom_hline(yintercept = 0:ny, colour = "grey85", linewidth = 0.3) +
  # Highlight the row
  geom_rect(aes(xmin = 0, xmax = nx, ymin = y_show - 0.5, ymax = y_show + 0.5),
            fill = "lightyellow", alpha = 0.5) +
  # Polygon
  geom_polygon(aes(x = poly_x, y = poly_y),
               fill = NA, colour = "black", linewidth = 0.8) +
  # Winding spans
  geom_segment(data = winding_df[winding_df$winding != 0, ],
               aes(x = x, xend = xend, y = y, yend = y),
               colour = "steelblue", linewidth = 2, alpha = 0.4) +
  # Crossing points with arrows
  geom_point(data = data.frame(x = xints, y = y_show),
             aes(x, y), colour = "firebrick", size = 3) +
  # Winding labels
  geom_text(data = winding_df,
            aes(x = (x + xend) / 2, y = y + 0.35,
                label = paste0("w=", winding)),
            size = 3.5, colour = "grey30") +
  # Scanline
  geom_hline(yintercept = y_show, colour = "steelblue",
             linewidth = 0.3, linetype = "dashed") +
  coord_equal(xlim = c(0, nx), ylim = c(0, ny)) +
  theme_grid() +
  labs(title = sprintf("Winding number sweep at y = %.1f", y_show))

## ----coverage-detail, fig.cap = "Coverage mode boundary cell detail. The polygon edge enters at the top side and exits at the left. The shaded region is the portion of the cell inside the polygon — its area divided by the cell area gives the coverage fraction (0.667). This is computed analytically from the entry/exit coordinates, not by sampling."----
# Cell (3.5, 4.5) where edge (5.2,5.8) → (2.8,4.2) cuts diagonally
cell_xmin <- 3; cell_xmax <- 4; cell_ymin <- 4; cell_ymax <- 5

# Edge enters top side at (4.0, 5.0), exits left side at (3.0, 4.333)
entry <- c(4.0, 5.0)
exit  <- c(3.0, 4 + 1/3)

# Covered region: polygon interior is below-right of the edge
covered_poly <- data.frame(
  x = c(cell_xmin, cell_xmax, cell_xmax, entry[1], exit[1]),
  y = c(cell_ymin, cell_ymin, cell_ymax, entry[2], exit[2])
)

frac <- abs(sum(covered_poly$x * c(covered_poly$y[-1], covered_poly$y[1]) -
                c(covered_poly$x[-1], covered_poly$x[1]) * covered_poly$y) / 2)

ggplot() +
  geom_rect(aes(xmin = cell_xmin, xmax = cell_xmax,
                ymin = cell_ymin, ymax = cell_ymax),
            fill = "grey95", colour = "black", linewidth = 0.8) +
  geom_polygon(data = covered_poly, aes(x, y),
               fill = "steelblue", alpha = 0.3) +
  geom_segment(aes(x = entry[1] + 0.3,
                   y = entry[2] + 0.3 * (exit[2] - entry[2]) / (exit[1] - entry[1]),
                   xend = exit[1] - 0.3,
                   yend = exit[2] - 0.3 * (exit[2] - entry[2]) / (exit[1] - entry[1])),
               colour = "black", linewidth = 0.8) +
  geom_point(aes(x = c(entry[1], exit[1]), y = c(entry[2], exit[2])),
             colour = "firebrick", size = 3) +
  annotate("text", x = 3.55, y = 4.3, label = "inside\n(covered)",
           colour = "steelblue", fontface = "bold", size = 4) +
  annotate("text", x = 3.25, y = 4.85, label = "outside",
           colour = "grey50", size = 3.5) +
  annotate("text", x = entry[1] + 0.06, y = entry[2] - 0.06, label = "entry",
           colour = "firebrick", size = 3, hjust = 0) +
  annotate("text", x = exit[1] - 0.06, y = exit[2], label = "exit",
           colour = "firebrick", size = 3, hjust = 1) +
  coord_equal(xlim = c(cell_xmin - 0.3, cell_xmax + 0.3),
              ylim = c(cell_ymin - 0.3, cell_ymax + 0.3)) +
  theme_grid() +
  labs(title = "Boundary cell: coverage fraction from edge traversal",
       subtitle = sprintf("fraction = %.3f (shaded area / cell area)", frac))

## ----resolution-scaling, fig.width = 9, fig.height = 3.5, fig.cap = "Resolution scaling: the same polygon at 8×7, 16×14, and 32×28 resolution. Interior cells (blue) grow as area (O(n²)); boundary cells (coral) grow as perimeter (O(n)). Sparse output stores only the blue runs and coral edges — the white cells cost nothing."----
library(geos)
g <- as_geos_geometry(sprintf("POLYGON ((%s))",
  paste(sprintf("%.1f %.1f", poly_x, poly_y), collapse = ", ")))

burns <- lapply(c(8, 16, 32), function(s) {
  r <- burn(g, extent = c(0, nx, 0, ny), dimension = c(as.integer(s), as.integer(s * ny/nx)))
  dx <- nx / s; dy_g <- ny / (s * ny/nx)
  
  cells <- data.frame(
    xmin = numeric(0), xmax = numeric(0),
    ymin = numeric(0), ymax = numeric(0),
    type = character(0)
  )
  
  if (nrow(r$runs) > 0) {
    for (i in seq_len(nrow(r$runs))) {
      row <- r$runs$row[i]
      for (col in r$runs$col_start[i]:r$runs$col_end[i]) {
        cells <- rbind(cells, data.frame(
          xmin = (col - 1) * dx, xmax = col * dx,
          ymin = ny - row * dy_g, ymax = ny - (row - 1) * dy_g,
          type = "interior"
        ))
      }
    }
  }
  
  if (nrow(r$edges) > 0) {
    for (i in seq_len(nrow(r$edges))) {
      row <- r$edges$row[i]; col <- r$edges$col[i]
      cells <- rbind(cells, data.frame(
        xmin = (col - 1) * dx, xmax = col * dx,
        ymin = ny - row * dy_g, ymax = ny - (row - 1) * dy_g,
        type = "boundary"
      ))
    }
  }
  
  cells$panel <- sprintf("%d × %d", s, as.integer(s * ny/nx))
  cells
})

all_cells <- do.call(rbind, burns)
all_cells$panel <- factor(all_cells$panel,
                          levels = c("8 × 7", "16 × 14", "32 × 28"))

ggplot(all_cells) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = type), colour = "grey80", linewidth = 0.1) +
  geom_polygon(data = data.frame(x = poly_x, y = poly_y),
               aes(x, y), fill = NA, colour = "black", linewidth = 0.6) +
  scale_fill_manual(values = c("interior" = "steelblue", "boundary" = "coral")) +
  facet_wrap(~panel) +
  coord_equal(xlim = c(0, nx), ylim = c(0, ny)) +
  theme_grid() +
  theme(strip.text = element_text(face = "bold"))

## ----zonal-stats-setup--------------------------------------------------------
library(terra)

# A 20×20 raster of "temperature" values: a gradient from 10 to 30
r_temp <- rast(nrows = 20, ncols = 20, xmin = 0, xmax = 10,
               ymin = 0, ymax = 10, crs = "")
values(r_temp) <- seq(10, 30, length.out = 400)

# A polygon that doesn't align with cell boundaries
poly <- geos::as_geos_geometry(
  "POLYGON ((1.3 2.1, 7.8 1.5, 6.2 8.3, 2.5 6.9, 1.3 2.1))")

## ----zonal-stats-burn---------------------------------------------------------
# Burn with exact coverage fractions
b <- burn(poly, extent = c(0, 10, 0, 10), dimension = c(20L, 20L))

# Extract raster values as a 20×20 matrix (row 1 = top = ymax)
vals_temp <- matrix(terra::values(r_temp), nrow = 20, ncol = 20, byrow = TRUE)

run_values <- numeric(0)
run_weights <- numeric(0)
for (i in seq_len(nrow(b$runs))) {
  row <- b$runs$row[i]
  cols <- b$runs$col_start[i]:b$runs$col_end[i]
  run_values <- c(run_values, vals_temp[row, cols])
  run_weights <- c(run_weights, rep(1.0, length(cols)))
}

# Extract raster values at edge (boundary) cells — fraction in (0, 1)
edge_values <- vals_temp[cbind(b$edges$row, b$edges$col)]
edge_weights <- b$edges$fraction

# Area-weighted mean: sum(value * weight) / sum(weight)
all_values <- c(run_values, edge_values)
all_weights <- c(run_weights, edge_weights)
weighted_mean <- sum(all_values * all_weights) / sum(all_weights)
cat(sprintf("Weighted mean temperature: %.4f\n", weighted_mean))

## ----zonal-stats-visualise, fig.width = 6, fig.height = 5.5, fig.cap = "Zonal statistics from sparse output. The raster shows a temperature gradient; the polygon boundary is overlaid. Interior cells (weight 1.0) and boundary cells (weight = coverage fraction) are combined for the area-weighted mean — the same result as exactextract."----
# Visualise: raster + polygon + coverage weights
mat <- materialize_chunk(b)
coverage_df <- expand.grid(
  col = seq_len(20), row = seq_len(20)
)
coverage_df$coverage <- as.vector(mat)
coverage_df$x <- (coverage_df$col - 0.5) * 0.5
coverage_df$y <- 10 - (coverage_df$row - 0.5) * 0.5

ggplot() +
  geom_raster(data = coverage_df[coverage_df$coverage > 0, ],
              aes(x, y, fill = coverage)) +
  scale_fill_gradient(low = "lightyellow", high = "steelblue",
                      name = "weight") +
  geom_polygon(data = data.frame(x = poly_x, y = poly_y),
               aes(x, y), fill = NA, colour = "black", linewidth = 0.8) +
  annotate("text", x = 5, y = 9.5,
           label = sprintf("weighted mean = %.2f", weighted_mean),
           fontface = "bold", size = 4) +
  coord_equal(xlim = c(0, 10), ylim = c(0, 10)) +
  theme_grid()

