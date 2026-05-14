# Compute monthly min / max 2 m temperature from the multi-variable hourly
# ERA5-Land NetCDFs that extract_era5land() produces.

.raster_is_readable <- function(path) {
  if (!file.exists(path)) return(FALSE)
  ok <- tryCatch(
    {
      terra::rast(path)
      TRUE
    },
    error = function(e) FALSE
  )
  if (!ok) {
    message(sprintf("Removing unreadable raster: %s", path))
    unlink(path)
  }
  ok
}

#' Monthly min / max 2 m temperature from ERA5-Land hourly NetCDFs
#'
#' For each (year, month) NetCDF produced by [extract_era5land()] — which
#' contains six variables interleaved — pulls out only the `t2m` layers,
#' converts to degC, reduces to daily min and max with [terra::tapp()], then
#' reduces those to a single monthly min and a single monthly max raster.
#' Per-month results are cached as GeoTIFFs; once all twelve months for a
#' year are present, a 12-layer annual stack is written.
#'
#' The input is the multi-variable NetCDF from `extract_era5land()`, so the
#' `t2m` layers are selected via the `t2m_valid_time=<unix>` layer-name
#' convention rather than by raster index.
#'
#' @param input_dir    Directory holding `<year>/<nc_template>` NetCDFs.
#' @param output_dir   Directory for per-month and annual GeoTIFFs. A
#'   `<year>/` subdirectory is created for the monthly caches.
#' @param years        Integer vector of years to process.
#' @param nc_template  `sprintf` template for the input NetCDF basename;
#'   receives `(year, month)`. Defaults to the convention used by
#'   [extract_era5land()].
#' @param month_min_template,month_max_template `sprintf` templates for the
#'   cached monthly min and max GeoTIFFs; each receives `(year, month)`.
#' @param annual_min_template,annual_max_template `sprintf` templates for
#'   the 12-layer annual stacks; each receives `(year)`.
#'
#' @return Character vector of annual stack paths written (invisibly).
#' @export
#'
#' @examples
#' \dontrun{
#' transform_era5land_minmax(
#'   input_dir  = "data/climate/era5land/input/nc",
#'   output_dir = "data/climate/era5land/output/minmax",
#'   years      = 2014:2017
#' )
#' }
transform_era5land_minmax <- function(
  input_dir,
  output_dir,
  years,
  nc_template          = "era5land_%d%02d.nc",
  month_min_template   = "era5land_monthly_min_%d%02d.tif",
  month_max_template   = "era5land_monthly_max_%d%02d.tif",
  annual_min_template  = "era5land_annual_minimum_temperature_%d.tif",
  annual_max_template  = "era5land_annual_maximum_temperature_%d.tif"
) {
  stopifnot(
    is.character(input_dir),  length(input_dir)  == 1L,
    is.character(output_dir), length(output_dir) == 1L,
    is.numeric(years),        length(years)      >= 1L
  )

  out_paths <- character()

  for (year in years) {
    year_out_dir <- file.path(output_dir, year)
    dir.create(year_out_dir, recursive = TRUE, showWarnings = FALSE)

    annual_min_path <- file.path(output_dir,
                                 sprintf(annual_min_template, year))
    annual_max_path <- file.path(output_dir,
                                 sprintf(annual_max_template, year))

    annual_min <- terra::rast()
    annual_max <- terra::rast()

    for (month in 1:12) {
      month_str <- sprintf("%02d", month)
      monthly_min_cache <- file.path(year_out_dir,
                                     sprintf(month_min_template, year, month))
      monthly_max_cache <- file.path(year_out_dir,
                                     sprintf(month_max_template, year, month))

      if (.raster_is_readable(monthly_min_cache) &&
            .raster_is_readable(monthly_max_cache)) {
        message(sprintf("Loading cached results for %d-%s", year, month_str))
        annual_min <- c(annual_min, terra::rast(monthly_min_cache))
        annual_max <- c(annual_max, terra::rast(monthly_max_cache))
        next
      }

      nc_path <- file.path(input_dir, year,
                           sprintf(nc_template, year, month))
      if (!file.exists(nc_path)) {
        message(sprintf("SKIP %d-%s: input NC missing (%s)",
                        year, month_str, nc_path))
        next
      }

      r <- terra::rast(nc_path)
      t2m_layers <- .layers_for_var(r, "t2m")
      if (length(t2m_layers$idx) == 0L) {
        message(sprintf("SKIP %d-%s: no t2m layers in %s",
                        year, month_str, basename(nc_path)))
        next
      }

      t2m <- r[[t2m_layers$idx]]
      terra::time(t2m) <- as.POSIXct(t2m_layers$ts,
                                     origin = "1970-01-01", tz = "UTC")
      t2m <- t2m - 273.15

      daily_min <- terra::tapp(t2m, index = "days", fun = min)
      daily_max <- terra::tapp(t2m, index = "days", fun = max)
      monthly_min <- min(daily_min)
      monthly_max <- max(daily_max)

      month_time <- as.POSIXct(
        sprintf("%d-%s-01 00:00:00", year, month_str), tz = "UTC"
      )
      terra::time(monthly_min) <- month_time
      terra::time(monthly_max) <- month_time
      names(monthly_min) <- format(month_time, "%Y-%m-%d")
      names(monthly_max) <- format(month_time, "%Y-%m-%d")

      terra::writeRaster(monthly_min, monthly_min_cache, overwrite = TRUE)
      terra::writeRaster(monthly_max, monthly_max_cache, overwrite = TRUE)

      annual_min <- c(annual_min, monthly_min)
      annual_max <- c(annual_max, monthly_max)
    }

    if (terra::nlyr(annual_min) == 12L && terra::nlyr(annual_max) == 12L) {
      terra::writeRaster(annual_min, annual_min_path, overwrite = TRUE)
      terra::writeRaster(annual_max, annual_max_path, overwrite = TRUE)
      out_paths <- c(out_paths, annual_min_path, annual_max_path)
      message(sprintf("Year %d: wrote annual min/max stacks", year))
    } else {
      message(sprintf(
        "Year %d: %d/12 months processed; annual stacks not written",
        year, terra::nlyr(annual_min)
      ))
    }
  }

  invisible(out_paths)
}
