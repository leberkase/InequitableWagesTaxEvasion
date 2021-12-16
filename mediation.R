
censor  <- function(.x) pmax(0, pmin(1, .x))



model.matrix.with.condition <- function(.fit, .condition) {
  mf <- model.frame(.fit)

  .condition %>%
    purrr::iwalk(~`<<-`(mf[[.y]], .x))
                 
  model.matrix(formula(.fit), data = mf)
}

#return n*K matrix of draws
sample_mediation <- function(.theta, .K, .censor, ...) {
  mm <- model.matrix.with.condition(...)

  n <- nrow(mm)
  beta <- .theta[-2]
  mu <- mm %*% beta
  sigma <- exp(.theta[2])

  samples <- rnorm(n * .K, mu, sigma)
  
  if(.censor)
    matrix(censor(samples), byrow = TRUE, nrow = n)
  else
    matrix(samples, byrow = TRUE, nrow = n)

}


sample_response <- function(.mediator, .fit, .condition, .theta, .K, .censor, ...) {

  
  mm_list <- .mediator %>%
    purrr::array_branch(2) %>%
    purrr::map(~append(.condition, list(estimate01 = .x))) %>%
    purrr::map(~model.matrix.with.condition(.condition = .x, .fit = .fit))


  mm <- do.call(rbind, mm_list)
  n <- nrow(.mediator)

  beta <- .theta[-2]
  mu <- mm %*% beta
  sigma <- exp(.theta[2])

  samples <- rnorm(n * .K, mu, sigma)
  
  if(.censor)
    matrix(censor(samples), byrow = TRUE, nrow = n)
  else
    matrix(samples, byrow = TRUE, nrow = n)
  
}



sample_mediation_statistic <- function(.theta_med, .theta_resp, .control, .treatment, .med_fit, .resp_fit, .censor_mediator=FALSE, .censor_response = FALSE, ...) {
  m0 <- sample_mediation(.fit = .med_fit, .condition = .control, .theta = .theta_med, .censor=.censor_mediator, ...)
  m1 <- sample_mediation(.fit = .med_fit, .condition = .treatment, .theta = .theta_med,  .censor=.censor_mediator, ...)

  y00 <- sample_response(.fit = .resp_fit, .condition = .control, .mediator = m0, .theta = .theta_resp,  .censor=.censor_response, ...)
  y10 <- sample_response(.fit = .resp_fit, .condition = .treatment, .mediator = m0, .theta = .theta_resp,  .censor=.censor_response, ...)
  y01 <- sample_response(.fit = .resp_fit, .condition = .control, .mediator = m1, .theta = .theta_resp,   .censor=.censor_response,...)
  y11 <- sample_response(.fit = .resp_fit, .condition = .treatment, .mediator = m1, .theta = .theta_resp,   .censor=.censor_response, ...)

  base_out <- data.frame(
    delta1 = mean(y11 - y10),
    delta0 = mean(y01 - y00),
    zeta1 = mean(y11 - y01),
    zeta0 = mean(y10 - y00)
  )


  within(base_out, {
    delta <- 0.5 * (delta1 + delta0)
    zeta <- 0.5 * (zeta1 + zeta0)
    tau <- 0.5 * (delta1 + zeta1 + delta0 + zeta0)
    nu <- 0.5 * (delta0 + delta1)/tau
  })

}




parametric_mediation <- function(.med_fit, .resp_fit, .J, .alpha = 0.05, .seed = NULL, ...) {
  # get J coefs from MVN from fit for med, resp
  theta_med_sample <- rmvnorm( .J, coef(.med_fit), vcov(.med_fit)) %>% array_branch(1)
  theta_resp_sample <- rmvnorm( .J, coef(.resp_fit), vcov(.resp_fit)) %>% array_branch(1)

  furr_opts <- if (is.null(.seed))
                   furrr::furrr_options()
  else
    furrr::furrr_options(seed = as.integer(.seed))

  stats <- purrr::map2_dfr(theta_med_sample, theta_resp_sample,
                              sample_mediation_statistic,
                              .med_fit = .med_fit,
                              .resp_fit = .resp_fit,
                              ## .options = furr_opts,
                              ...) %>%
    as.data.frame() %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "draw") %>%
    group_by(parameter) %>%
    do(mean = mean(.$draw), hdi = HDInterval::hdi(.$draw, credMass = 1-.alpha)) %>%
    transmute(parameter, lower = hdi[1], mean = mean[1], upper = hdi[2])

}















