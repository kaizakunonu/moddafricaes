# Extract bounding-box coordinates for a given spatial extent.

#' Get the combined bounding box for selected DRC provinces
#'
#' Reads the bundled `province26` shapefile, filters to the requested
#' provinces, and returns the combined bounding box in Copernicus CDS
#' order (N, W, S, E). Optionally prints per-province centroids and
#' bounding boxes for inspection.
#'
#' @param provinces  Character vector of province names; matched against
#'   `name_field` in the shapefile.
#' @param shp_zip    Path to the shapefile zip, relative to the project root.
#' @param shp_path   Path inside the zip to the `.shp` file.
#' @param name_field Attribute column holding province names.
#' @param crs        Target CRS as an EPSG code.
#' @param verbose    If `TRUE`, print centroids and per-province bboxes.
#'
#' @return Named numeric vector of length 4 in CDS order: `c(N, W, S, E)`.
#' @export
#'
#' @examples
#' provinces = c("Kinshasa")
get_province_bbox <- function(provinces  = c("Kinshasa"),
                              shp_zip    = "data-acquisition/utilities/province26.zip",
                              shp_path   = "provinces26/Province26.shp",
                              name_field = "NOM",
                              crs        = 4326,
                              verbose    = FALSE) {
  stopifnot(
    is.character(provinces), length(provinces) >= 1L,
    is.character(shp_zip),   length(shp_zip)   == 1L,
    file.exists(shp_zip)
  )

  shp_src <- paste0("/vsizip/", shp_zip, "/", shp_path)
  shp <- sf::st_read(shp_src, quiet = TRUE)
  shp <- sf::st_transform(shp, crs)

  sel <- shp[shp[[name_field]] %in% provinces, ]
  missing <- setdiff(provinces, sel[[name_field]])
  if (length(missing) > 0L) {
    stop("Province(s) not found in shapefile: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  if (verbose) {
    ctr    <- suppressWarnings(sf::st_centroid(sel))
    coords <- sf::st_coordinates(ctr)
    centroids <- data.frame(
      name = sel[[name_field]],
      lon  = coords[, 1L],
      lat  = coords[, 2L]
    )
    print(centroids, row.names = FALSE)

    cat("\nBounding boxes:\n")
    for (i in seq_len(nrow(sel))) {
      bb <- sf::st_bbox(sel[i, ])
      cat(sprintf("%-20s xmin=%.5f ymin=%.5f xmax=%.5f ymax=%.5f\n",
                  sel[[name_field]][i],
                  bb[["xmin"]], bb[["ymin"]], bb[["xmax"]], bb[["ymax"]]))
    }
  }

  bb_all <- sf::st_bbox(sel)
  c(N = unname(bb_all[["ymax"]]),
    W = unname(bb_all[["xmin"]]),
    S = unname(bb_all[["ymin"]]),
    E = unname(bb_all[["xmax"]]))
}

# Run with defaults when invoked as a script (Rscript or top-level),
# but stay silent when sourced for its function definition.
if (sys.nframe() == 0L) {
  bbox <- get_province_bbox(verbose = TRUE)
  cat("\nCombined bounding box (N, W, S, E for CDS):\n")
  cat(sprintf("c(%.5f, %.5f, %.5f, %.5f)\n",
              bbox[["N"]], bbox[["W"]], bbox[["S"]], bbox[["E"]]))
}
