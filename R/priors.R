#' @keywords internal
HGnedin <- function(V, h, gamma = 0.5) {
  exp(lchoose(V, h) + lgamma(h - gamma) - lgamma(1 - gamma) +
        log(gamma) + lgamma(V + gamma - h) - lgamma(V + gamma))
}

#' @keywords internal
urn_DP <- function(v_minus, alpha_PY) c(v_minus, alpha_PY)

#' @keywords internal
urn_PY <- function(v_minus, alpha_PY, sigma_PY) {
  H <- length(v_minus)
  c(v_minus - sigma_PY, alpha_PY + H * sigma_PY)
}

#' @keywords internal
urn_DM <- function(v_minus, beta_DM, H_DM) {
  H <- length(v_minus)
  c(v_minus + beta_DM, beta_DM * (H_DM - H) * (H_DM > H))
}

#' @keywords internal
urn_GN <- function(v_minus, gamma) {
  H <- length(v_minus)
  n_ <- sum(v_minus) + 1
  c((v_minus + 1) * (n_ - H + gamma), H^2 - H * gamma)
}
