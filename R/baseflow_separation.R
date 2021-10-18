# baseflow_separation.R

#' Digital filter
#'
#' Digital filter for baseflow separation
#' @param y Baseflow timeseries
#' @param alpha Filter parameters (0.9 - 0.95)
digital_filter <- function(y,alpha) {
  # x is the streamflow timeseries as a vector (no timestamps)
  f <- rep(0,length(y))
  for (k in 2:length(y)) {
    f[k] <- alpha * f[k-1] + (1+alpha)/2 * (y[k] - y[k-1])
  }
  return(f)
}


#' Recursive digital filter
#'
#' @param y Baseflow timeseries
#' @param alpha Filter parameters (0.9 - 0.95)
#' @export
#' @details
#' This functions executes the recusive digital filter for baseflow separation as
#' described by Nathan and McMahon (1990). The digital filter operates by
#' removing out high frequency components of streamflow. These high frequencies can
#' be considered quickflow which is calculated as:
#'
#' \code{f[k] = alpha * f[k-1] + (1+alpha)/2 * (y[k] - y[k-1])}
#'
#' The filter is run forwards, then backwards, then forwards again. Nathan and McMahon
#' found that a value of 0.9 closely matched baseflow separation via manual
#' techniques.
#'
#' Nathan, R. J., & McMahon, T. A. (1990). Evaluation of automated techniques for base flow
#' and recession analyses. Water resources research, 26(7), 1465-1473.
#' @examples
#'
#' # Generate random streamflow timeseries (not necessarily realistic)
#' kern <- rep(1/15, 15)
#' flow <- data.frame(t = seq(1, 365 * 4 + length(kern)))
#' flow$sin_annual <- cos(flow$t * 2 * pi / 365 - pi)
#' flow$rain_sin <- sin(flow$t * 2 * pi / 365 - )
#' flow$total_flow <- flow$sin_annual + (flow$sin_annual + 1) * as.numeric(stats::filter(exp(rnorm(nrow(flow))), kern)) + 1
#' flow$baseflow90 <- recursive_digital_filter(flow$total_flow, 0.9)
#' flow$baseflow95 <- recursive_digital_filter(flow$total_flow, 0.95)
#'
#' # Plot the results
#' library(ggplot2)
#' ggplot(flow) +
#'   geom_line(aes(t, total_flow, color = "Total flow")) +
#'   geom_line(aes(t, baseflow90, color = "Baseflow, alpha = 0.9")) +
#'   geom_line(aes(t, baseflow95, color = "Baseflow, alpha = 0.95"))
recursive_digital_filter <- function(x,alpha) {
  # x is the streamflow timeseries as a vector (no timestamps)
  cat("Running baseflow separation based on recursive digital filter (Nathan and McMahon, 1990)\n")
  if (max(is.na(x))) {
    cat("Caution, data has NA values!\n")
  }
  na_idx <- is.na(x)
  y <- x[!na_idx]

  f1 <- digital_filter(y,alpha)
  f1_fix <- pmin(pmax(f1,0),y) # filtered quick response cannot be <0 or >flow
  qb1 <- y - f1_fix
  f2 <- rev(digital_filter(rev(qb1),alpha))
  f2_fix <- pmin(pmax(f2,0),qb1) # filtered quick response cannot be <0 or >flow
  qb2 <- qb1 - f2_fix
  f3 <- digital_filter(qb2,alpha)
  f3_fix <- pmin(pmax(f3,0),qb2) # filtered quick response cannot be <0 or >flow
  qb3 <- qb2 - f3_fix

  # ggplot() + geom_line(aes(1:length(y),y),color="blue") +
  #   geom_line(aes(1:length(y),qb1),color="red",size=0.25) +
  #   geom_line(aes(1:length(y),qb2),color="orange",size=0.25) +
  #   geom_line(aes(1:length(y),qb3),color="green",size=0.5) +
  #   xlim(c(1,365*2))

  qb <- rep(NA,length(x))
  qb[!na_idx] <- qb3
  return(qb)
}
