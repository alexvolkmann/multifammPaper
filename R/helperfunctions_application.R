#' Plot the Snooker Fit per Observation
#'
#' This function plots the observed two-dimensional snooker trajectories against
#' their fitted trajectories.
#'
#' @param obs Named level of the observation to be plotted.
#' @param model Model to be evaluated.
#' @param model_comp Additional model if two models are to be compared.
plot_snooker_fit <- function(obs, model, model_comp = NULL) {

  # If two models are to be compared
  if (!is.null(model_comp)) {
    dat <- cbind(model$data, fit = model$model$fitted.values,
                 fit_comp = model_comp$model$fitted.values)
    cols <- c("y_vec", "fit", "fit_comp")
    col_vals <- c("firebrick3", "steelblue3", "gray44")
    labs <- c("Fit", "Sensitivity Analysis", "Data")
  } else {
    dat <- cbind(model$data, fit = model$model$fitted.values)
    cols <- c("y_vec", "fit")
    col_vals <- c("steelblue3", "gray44")
    labs <- c("Fit", "Data")
  }

  # Prepare the data for plotting
  dat <- dat %>%
    filter(n_long %in% obs) %>%
    droplevels() %>%
    separate(dim, into = c("location", "axis"), "\\.") %>%
    pivot_longer(cols = all_of(cols), names_to = "color",
                 values_to = "val") %>%
    mutate_at(c("color", "location", "axis"), as.factor) %>%
    pivot_wider(id_cols = c("n_long", "t", "location", "color"),
                names_from = axis,
                values_from = val) %>%
    unite(curve_group, location, color, sep = ".", remove = FALSE) %>%
    mutate(n_long = factor(n_long, levels = obs, labels = names(obs)))

  # Prepare a datset to mark the beginnings of the trajectories
  beg <- dat %>%
    group_by(n_long, curve_group) %>%
    filter(color != "0") %>%
    slice(which.min(t)) %>%
    ungroup()

  # Plot the observation with its fit
  plot_list <- lapply(levels(dat$location), function (loc) {
    p <- ggplot2::ggplot(data = dat %>% filter(location == loc),
                         aes(x = x, y = y, group = curve_group, colour = color,
                             linetype = color)) +
      geom_path() +
      theme_bw(base_size = 8) +
      theme(plot.title = element_text(size=8),
            legend.position = "bottom") +
      ylab("Y") +
      xlab("X") +
      geom_point(data = beg %>% filter(location == loc), aes(x = x, y = y),
                 shape = 8, size = 0.8, colour = "black", alpha = 0.5) +
      facet_grid(location~n_long) +
      scale_color_manual(name = "", values = col_vals, labels = labs) +
      scale_linetype_manual(name = "", values = c(1, 1, 2), labels = labs)
    if(loc != "shoulder") {
      p <- p + guides(color = FALSE, linetype = FALSE)
    }
    p
  })

}


#' Function to Estimate the Covariance Operator
#'
#' This function conducts a fast symmetric covariance estimation on some
#' data.table as proposed by Cederbaum et al. (2018). The parameters in the
#' function definition depend on the implementation of the sparseflmm(). Given
#' the estimated covariance surface, the function also conducts an MFPCA thus
#' giving FPCs.
#'
#' @param data Data.table containig the information needed for the sparseflmm()
#'   function.
#' @param m_mean Order of penalty for basis function (as in sparseFLMM).
#' @param covariate Covariate effects (as in sparseFLMM).
#' @param num_covariates Number of covariates included in the model (as in
#'   sparseFLMM).
#' @param covariate_form Vector of strings for type of covariate (as in
#'   sparseFLMM).
#' @param interaction TRUE if there are interactions between covariates (as in
#'   sparseFLMM). Defaults to \code{FALSE}.
#' @param which_interaction Symmetric matrix specifying the interaction terms
#'   (as in sparseFLMM).
estim_overall_cov <- function(data, m_mean = c(2, 1),
                              covariate = FALSE, num_covariates = 4,
                              covariate_form = rep("by", times = 4),
                              interaction = FALSE,
                              which_interaction = matrix(NA)) {

  model_list <- multifamm:::apply_sparseFLMM(
    fRI_B = FALSE, fRI_C = FALSE, nested = FALSE, bs = "ps", bf_mean = 8,
    bf_covariates = 8, m_mean = m_mean, covariate = covariate,
    num_covariates = num_covariates, covariate_form = covariate_form,
    interaction = interaction, which_interaction = which_interaction,
    bf_covs = 5, m_covs = list(m_mean), var_level = 1, use_famm = FALSE,
    save_model_famm = FALSE, one_dim = NULL, data = data)

  mfpca_info <- multifamm:::prepare_mfpca(model_list = model_list,
                                          fRI_B = FALSE, mfpc_weight = FALSE)

  multifamm:::conduct_mfpca(mfpca_info, mfpc_weight)
}


#' Function to Approximate an Autocorrelation Function Based on a Surface
#'
#' This function takes a data set given by the covariance_surf_plot_helper()
#' function and approximates an autocorrelation function by taking the mean of
#' the off-diagonals.
#'
#' @param cov_dat Data set containing the surface information as given by
#'   covariance_surf_plot_helper().
#' @param lag Maximal lag possible. Default is 100 such as the evaluation points
#'   of the surface.
approx_autocor <- function(cov_dat, lag = 100) {
  dat <- cov_dat %>% filter(row_dim == col_dim)
  dat <- split(dat, dat$row_dim)
  acl <- lapply(dat, function (dim) {
    V <- matrix(dim$value, nrow = 100, ncol = 100, byrow = FALSE)
    data.frame(x = seq(0, 1, length.out = 100),
               acl = sapply(seq_len(lag) - 1, function (l) {
                 mean(V[row(V) == col(V) - l])
               }))
  })
  dat <- bind_rows(acl, .id = "dim") %>%
    mutate(dim = as.factor(dim))
}
