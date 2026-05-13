# Extract bounding-box coordinates for a given spatial extent.

#' Get the bounding box of a shapefile
#'
#' Reads a shapefile from `shp_path`, reprojects it to `crs`, and returns
#' its bounding box in Copernicus CDS order (N, W, S, E).
#'
#' `shp_path` may be a plain `.shp` path or a GDAL virtual path such as
#' `"/vsizip/<archive.zip>/<file.shp>"` for shapefiles inside a zip.
#'
#' @param shp_path Path to the shapefile (any source `sf::st_read()` accepts).
#' @param crs      Target CRS as an EPSG code.
#'
#' @return Named numeric vector of length 4 in CDS order: `c(N, W, S, E)`.
#' @export
#'
#' @examples
#' \dontrun{
#' get_bbox("/vsizip/data-acquisition/utilities/province26.zip/provinces26/Province26.shp")
#' }
get_bbox <- function(shp_path, crs = 4326) {
  stopifnot(is.character(shp_path), length(shp_path) == 1L)

  shp <- sf::st_read(shp_path, quiet = TRUE)
  shp <- sf::st_transform(shp, crs)
  bb  <- sf::st_bbox(shp)

  c(N = unname(bb[["ymax"]]),
    W = unname(bb[["xmin"]]),
    S = unname(bb[["ymin"]]),
    E = unname(bb[["xmax"]]))
}

# Run with a sample shapefile when invoked as a script.
if (sys.nframe() == 0L) {
  shp <- "/vsizip/data-acquisition/utilities/province26.zip/provinces26/Province26.shp"
  bbox <- get_bbox(shp)
  cat("Bounding box (N, W, S, E for CDS):\n")
  cat(sprintf("c(%.5f, %.5f, %.5f, %.5f)\n",
              bbox[["N"]], bbox[["W"]], bbox[["S"]], bbox[["E"]]))
}
