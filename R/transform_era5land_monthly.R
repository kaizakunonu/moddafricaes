# Transform a CDS ERA5-Land monthly-means NetCDF (one file, N years x 12
# months x 6 vars) into a single analysis Parquet.

.expected_cols_monthly <- c(
  "timestamp", "year", "month",
  "lat", "lon", "elevation_m",
  "t2m_celsius", "tp_mm",
  "wind_speed_mps", "relative_humidity_pct", "vpd_hpa",
  "ssrd_wm2", "ssrd_mj_m2"
)

.expected_schema_monthly <- function() {
  arrow::schema(
    timestamp             = arrow::timestamp("us", "UTC"),
    year                  = arrow::int32(),
    month                 = arrow::int32(),
    lat                   = arrow::float32(),
    lon                   = arrow::float32(),
    elevation_m           = arrow::float32(),
    t2m_celsius           = arrow::float32(),
    tp_mm                 = arrow::float32(),
    wind_speed_mps        = arrow::float32(),
    relative_humidity_pct = arrow::float32(),
    vpd_hpa               = arrow::float32(),
    ssrd_wm2              = arrow::float32(),
    ssrd_mj_m2            = arrow::float32()
  )
}

.default_metadata_monthly <- function(version) {
  c(
    timestamp_convention = "month_start_utc",
    timestamp_meaning    = "value is the monthly aggregate for the calendar month beginning at this timestamp",
    tp_mm_meaning        = "total precipitation accumulation over the calendar month, in mm",
    ssrd_wm2_meaning     = "monthly mean downward shortwave flux at the surface, in W/m^2",
    ssrd_mj_m2_meaning   = "total downward shortwave energy at the surface over the calendar month, in MJ/m^2",
    elevation_m_meaning  = "Surface elevation (m); DEM averaged onto the ERA5-Land 0.1 deg grid; cells with no DEM coverage are filled by nearest-neighbour from the closest valued cell.",
    source_product       = "CDS reanalysis-era5-land-monthly-means / monthly_averaged_reanalysis",
    pipeline_version     = version
  )
}

#' Transform an ERA5-Land monthly-means NetCDF into a Parquet
#'
#' Reads a single CDS `reanalysis-era5-land-monthly-means` NetCDF, converts
#' temperatures to degC, computes wind speed / RH / VPD, derives monthly
#' precipitation (mm) and shortwave totals from the per-day means (the CDS
#' monthly product stores `ssrd` and `tp` as means of daily totals, **not**
#' cumulative since 00 UTC — multiply by `days_in_month`, do **not**
#' de-cumulate), joins a DEM-derived `elevation_m`, and writes one Parquet.
#'
#' Compared to [transform_era5land()] (hourly): the input is a single file,
#' cell x month volume is trivial so there is no chunking, and the
#' de-accumulation logic does not apply.
#'
#' @param input_nc    Path to the monthly-means NetCDF.
#' @param output_path Output Parquet path. Parent directory is created.
#' @param dem_tif     Path to a DEM raster covering the AOI.
#' @param version     Version tag embedded in the Parquet footer metadata.
#' @param compression Parquet compression codec.
#' @param metadata    Named character vector embedded in the Parquet footer.
#'   Defaults to a generic ERA5-Land monthly-means dictionary.
#'
#' @return The output Parquet path (invisibly).
#' @export
#'
#' @examples
#' \dontrun{
#' transform_era5land_monthly(
#'   input_nc    = "data/.../era5land-kinshasa-monthly-2014-2017.nc",
#'   output_path = "data/.../era5land_monthly_kinshasa_2014-2017.parquet",
#'   dem_tif     = "RDC_Elevation_Complete.tif"
#' )
#' }
transform_era5land_monthly <- function(input_nc,
                                       output_path,
                                       dem_tif,
                                       version     = "v0.1.0",
                                       compression = "zstd",
                                       metadata    = NULL) {
  stopifnot(
    is.character(input_nc),    length(input_nc)    == 1L, file.exists(input_nc),
    is.character(output_path), length(output_path) == 1L,
    is.character(dem_tif),     length(dem_tif)     == 1L, file.exists(dem_tif)
  )
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

  if (.parquet_is_readable(output_path)) {
    message(sprintf("Already done: %s", basename(output_path)))
    return(invisible(output_path))
  }

  r <- terra::rast(input_nc)
  vars_present <- unique(sub("_valid_time=.*", "", names(r)))
  for (v in c("ssrd", "d2m", "t2m", "tp", "u10", "v10")) {
    stopifnot(v %in% vars_present)
  }

  elev_all <- .build_elevation_vector(dem_tif, r[[1]])
  stopifnot(length(elev_all) == terra::ncell(r))
  message(sprintf("Elevation vector: %d cells, range [%.1f, %.1f] m",
                  length(elev_all),
                  min(elev_all, na.rm = TRUE),
                  max(elev_all, na.rm = TRUE)))

  vidx <- list(
    d2m  = .layers_for_var(r, "d2m"),
    t2m  = .layers_for_var(r, "t2m"),
    tp   = .layers_for_var(r, "tp"),
    u10  = .layers_for_var(r, "u10"),
    v10  = .layers_for_var(r, "v10"),
    ssrd = .layers_for_var(r, "ssrd")
  )
  ts_ref <- vidx$d2m$ts
  for (v in c("t2m", "tp", "u10", "v10", "ssrd")) {
    stopifnot(identical(vidx[[v]]$ts, ts_ref))
  }
  timestamps  <- as.POSIXct(ts_ref, origin = "1970-01-01", tz = "UTC")
  n_months    <- length(timestamps)
  coords_full <- terra::xyFromCell(r[[1]], seq_len(terra::ncell(r)))

  message(sprintf("File %s: %d months, %d cells",
                  basename(input_nc), n_months, terra::ncell(r)))
  t0 <- Sys.time()

  m_t2m  <- terra::values(r[[vidx$t2m$idx]])
  m_d2m  <- terra::values(r[[vidx$d2m$idx]])
  m_tp   <- terra::values(r[[vidx$tp$idx]])
  m_u10  <- terra::values(r[[vidx$u10$idx]])
  m_v10  <- terra::values(r[[vidx$v10$idx]])
  m_ssrd <- terra::values(r[[vidx$ssrd$idx]])

  if (mean(m_t2m, na.rm = TRUE) > 200) {
    m_t2m <- m_t2m - 273.15
    m_d2m <- m_d2m - 273.15
  }

  # Sanity check: ssrd in the monthly-means product is the mean of *daily*
  # totals in J/m^2, so the cell-wise mean lands in ~[1e6, 5e7]. Hourly
  # accumulation would be ~24x smaller, monthly accumulation ~30x larger --
  # either tripwires this check.
  ssrd_mean <- mean(m_ssrd, na.rm = TRUE)
  if (!(ssrd_mean > 1e6 && ssrd_mean < 5e7)) {
    stop(sprintf(
      "ssrd mean = %.3g J/m^2 is outside the expected daily-total range [1e6, 5e7]. Wrong CDS product?",
      ssrd_mean))
  }

  valid <- !is.na(m_t2m[, 1]) & !is.na(m_d2m[, 1])
  m_t2m  <- m_t2m[valid, , drop = FALSE]
  m_d2m  <- m_d2m[valid, , drop = FALSE]
  m_tp   <- m_tp[valid, , drop = FALSE]
  m_u10  <- m_u10[valid, , drop = FALSE]
  m_v10  <- m_v10[valid, , drop = FALSE]
  m_ssrd <- m_ssrd[valid, , drop = FALSE]
  coords <- coords_full[valid, , drop = FALSE]
  v_elev <- elev_all[valid]
  n_pix  <- nrow(coords)

  days_vec    <- as.integer(lubridate::days_in_month(timestamps))
  m_tp_mm     <- sweep(m_tp,   2, days_vec * 1000, `*`)
  m_ssrd_wm2  <- m_ssrd / 86400
  m_ssrd_mjm2 <- sweep(m_ssrd, 2, days_vec / 1e6, `*`)

  m_es  <- 6.112 * exp((17.67 * m_t2m) / (m_t2m + 243.5))
  m_ea  <- 6.112 * exp((17.67 * m_d2m) / (m_d2m + 243.5))
  m_rh  <- 100 * m_ea / m_es
  m_vpd <- m_es - m_ea
  m_ws  <- sqrt(m_u10^2 + m_v10^2)

  df <- data.frame(
    timestamp             = rep(timestamps, each = n_pix),
    lon                   = rep(coords[, 1], times = n_months),
    lat                   = rep(coords[, 2], times = n_months),
    elevation_m           = round(rep(v_elev, times = n_months), 1),
    t2m_celsius           = round(as.vector(m_t2m),       3),
    tp_mm                 = round(as.vector(m_tp_mm),     3),
    wind_speed_mps        = round(as.vector(m_ws),        3),
    relative_humidity_pct = round(as.vector(m_rh),        2),
    vpd_hpa               = round(as.vector(m_vpd),       3),
    ssrd_wm2              = round(as.vector(m_ssrd_wm2),  2),
    ssrd_mj_m2            = round(as.vector(m_ssrd_mjm2), 2),
    stringsAsFactors      = FALSE
  )
  df$year  <- lubridate::year(df$timestamp)
  df$month <- lubridate::month(df$timestamp)
  df <- df[, .expected_cols_monthly]

  tbl <- arrow::as_arrow_table(df, schema = .expected_schema_monthly())
  tbl$metadata <- as.list(
    if (is.null(metadata)) .default_metadata_monthly(version) else metadata
  )
  arrow::write_parquet(tbl, output_path, compression = compression)

  if (!.parquet_is_readable(output_path)) {
    stop(sprintf("Output %s failed readability check", output_path))
  }
  message(sprintf("Wrote %s: %s rows, %.1fs, %.2f MB",
                  basename(output_path),
                  format(nrow(df), big.mark = ","),
                  as.numeric(Sys.time() - t0, units = "secs"),
                  file.size(output_path) / 1024^2))
  invisible(output_path)
}
