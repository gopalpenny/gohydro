# Interpolation

#' Inverse distance squared interpolation
#' @param precip_df Must contain precip, x, y columns (elev optional)
#' @param out_xy Coordinates (x, y) at which to interpolate (elev optional)
#' @param elev_gradient Linear elevation gradient (see details)
#' @details
#' For elevation gradient, \code{elev_gradient} should be specified in units
#' of precipitation over elevation (i.e., if precipitation is in mm and
#' elevation in m, \code{elev_gradient} should be [mm/m]). Additionally,
#' \code{precip_df} and \code{out_xy} must contain a column named \code{elev}.
#' @examples
#' # Generate topography with sin functions
#' get_elev <- function(x, y) 500*sin(x / 750) + 500*sin(y / 1000) + 1000
#'
#' # Generate random precipitation observations
#' set.seed(100)
#' N_obs<- 10
#' precip_df <- data.frame(precip = rnorm(N_obs, 750, 100),
#'                         x = runif(N_obs, 0, 5000),
#'                         y = runif(N_obs, 0, 5000))
#' precip_df$elev <- get_elev(precip_df$x, precip_df$y)
#'
#' out_xy <- expand.grid(x = seq(0, 5000, length.out = 9),
#'                       y = seq(0, 5000, length.out = 9))
#' out_xy$elev <- get_elev(out_xy$x, out_xy$y)
#'
#' # Plot elevation
#' library(ggplot2)
#' ggplot() +
#'   geom_raster(data = out_xy, aes(x, y, fill = elev)) + coord_equal()
#'
#' # Interpolate precipitation
#' precip_interp <- interpolate_IDS(precip_df, out_xy, elev_gradient = 0)
#' precip_interp <- interpolate_IDS(precip_df, out_xy, elev_gradient = 0.1)
#'
#' # Plot the data
#' ggplot() +
#'   geom_raster(data = precip_interp, aes(x, y, fill = precip)) +
#'   geom_point(data = precip_df, aes(x, y, fill = precip), shape = 21, size = 3) +
#'   scale_fill_viridis_c() +
#'   coord_equal()
interpolate_IDS <- function(precip_df, out_xy, elev_gradient = 0) {
  # If elev_gradient != 0, convert to common elevation (mean(elev))
  elev_common <- mean(precip_df$elev)
  precip_prep <- precip_df %>%
    dplyr::mutate(precip_common_elev = precip + (elev_common - elev) * elev_gradient) %>%
    tidyr::crossing(out_xy %>% dplyr::rename(x_out = x, y_out = y, elev_out = elev))
  precip_interp <- precip_prep %>%
    dplyr::mutate(dist2 = (x_out - x)^2 + (y_out - y)^2,
                  inv_dist2 = 1/dist2)  %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c("x_out", "y_out", "elev_out")))) %>%
    dplyr::summarize(precip_common_elev_interp = weighted.mean(precip_common_elev, inv_dist2)) %>%
    dplyr::mutate(precip_interp = precip_common_elev_interp + (elev_out - precip_common_elev_interp) * elev_gradient)

  return(precip_interp %>% dplyr::select(x = x_out, y = y_out, precip = precip_interp))
}
