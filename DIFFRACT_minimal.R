# DIFFRACT: minimal WT versus KO example .
# Run from this directory: source("DIFFRACT_minimal.R")
# Dependencies: install.packages(c("igraph", "infotheo", "energy"))
# Reference: Hause et al. (2026): "Local fractal dimension as an indicator of functional reorganization in molecular interaction networks" 
# https://doi.org/10.64898/2026.10.27.762149

diffract <- function(input_dir = "example_data", output_dir = "output") {
  required <- c("igraph", "infotheo", "energy")
  missing <- required[!base::vapply(required, base::requireNamespace,
                                  base::logical(1), quietly = TRUE)]
  if (base::length(missing)) {
    base::stop("Install required packages: ", base::paste(missing, collapse = ", "))
  }
  backbone <- igraph::read_graph(base::file.path(input_dir, "STRINGfull_v12_sample.graphml"),
                                  format = "graphml")
  expression <- utils::read.csv(base::file.path(input_dir, "proteinabundance_WT_KO.csv"),
                                check.names = FALSE, stringsAsFactors = FALSE)
  groups <- base::list(
    WT = base::grep("^WT[0-9]+$", base::names(expression), value = TRUE),
    KO = base::grep("^KO[0-9]+$", base::names(expression), value = TRUE)
  )
  if (!"node" %in% base::names(expression) ||
      base::anyDuplicated(expression$node) ||
      base::anyNA(expression$node) || base::any(!base::nzchar(expression$node))) {
    base::stop("Expression data require unique, nonempty neutral node identifiers.")
  }
  if (base::any(base::lengths(groups) < 3L)) {
    base::stop("At least three WT and three KO samples are required.")
  }
  genes <- igraph::vertex_attr(backbone, "name")
  if (base::is.null(genes) || base::anyDuplicated(genes) ||
      !base::setequal(genes, expression$node) ||
      igraph::is_directed(backbone) || !igraph::is_simple(backbone) ||
      !igraph::is_connected(backbone)) {
    base::stop("The backbone must be simple, undirected, connected and match node exactly.")
  }
  expression <- expression[base::match(genes, expression$node), , drop = FALSE]
  values <- base::as.matrix(expression[, base::unlist(groups), drop = FALSE])
  if (!base::is.numeric(values) || base::any(!base::is.finite(values))) {
    base::stop("Example expression values must be finite numeric values.")
  }
  base::rownames(values) <- genes
  base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  parameters <- base::list(
    methods = c("MI", "Pearson", "dCor", "VarMiss"),
    epsilon = 1e-2, alpha = 2.0, beta = 4,
    gamma = 1.3, n_bins = 100L,
    mi_bins = 2L, varmiss_mi_bins = 3L, discretization = "equalfreq",
    min_valid = 3L, min_unique = 3L, varmiss_min_valid = 2L,
    use_weights = TRUE, invert_weights = FALSE,
    contrast_sd_epsilon = 1e-6, p_cutoff = 0.05, min_component_size = 5L
  )
  parameter_table <- base::data.frame(
    parameter = base::names(parameters),
    value = base::vapply(parameters, function(x) base::paste(x, collapse = ", "),
                         base::character(1))
  )
  utils::write.csv(parameter_table, base::file.path(output_dir, "parameters.csv"),
                   row.names = FALSE)
  utils::capture.output(utils::sessionInfo(),
                        file = base::file.path(output_dir, "sessionInfo.txt"))

  edges <- igraph::as_edgelist(backbone, names = TRUE)
  build_weighted_graph <- function(condition, method) {
    xmat <- values[, groups[[condition]], drop = FALSE]
    variances <- base::apply(xmat, 1L, stats::var, na.rm = TRUE)
    association <- base::rep(NA_real_, base::nrow(edges))
    for (i in base::seq_len(base::nrow(edges))) {
      x <- xmat[edges[i, 1L], ]
      y <- xmat[edges[i, 2L], ]
      if (method == "VarMiss") {
        if (base::sum(!base::is.na(x) & x != 0) < parameters$varmiss_min_valid ||
            base::sum(!base::is.na(y) & y != 0) < parameters$varmiss_min_valid) next
        association[i] <- base::tryCatch({
          xd <- infotheo::discretize(x, disc = parameters$discretization,
                                     nbins = parameters$varmiss_mi_bins)
          yd <- infotheo::discretize(y, disc = parameters$discretization,
                                     nbins = parameters$varmiss_mi_bins)
          infotheo::mutinformation(xd, yd)
        }, error = function(e) NA_real_)
        penalty <- (variances[edges[i, 1L]] + variances[edges[i, 2L]]) / 2
        association[i] <- association[i] / (penalty + parameters$epsilon)^parameters$beta
      } else {
        keep <- base::is.finite(x) & base::is.finite(y) & x != 0 & y != 0
        x <- base::as.numeric(x[keep])
        y <- base::as.numeric(y[keep])
        if (base::length(x) < parameters$min_valid ||
            base::length(base::unique(x)) < parameters$min_unique ||
            base::length(base::unique(y)) < parameters$min_unique) next
        if (!base::is.finite(stats::var(x)) || !base::is.finite(stats::var(y)) ||
            stats::var(x) == 0 || stats::var(y) == 0) next
        association[i] <- base::tryCatch({
          if (method == "MI") {
            xd <- infotheo::discretize(x, disc = parameters$discretization,
                                       nbins = parameters$mi_bins)
            yd <- infotheo::discretize(y, disc = parameters$discretization,
                                       nbins = parameters$mi_bins)
            infotheo::mutinformation(base::as.integer(xd[[1L]]),
                                     base::as.integer(yd[[1L]]))
          } else if (method == "Pearson") {
            base::abs(stats::cor(x, y, method = "pearson", use = "complete.obs"))
          } else if (method == "dCor") {
            energy::dcor(x, y)
          }
        }, error = function(e) NA_real_)
      }
    }
    
    if (base::any(!base::is.finite(association) | association < 0)) {
      base::stop("Unevaluable backbone edges for ", method, "/", condition,
                 ". Supply complete, variable profiles for all backbone proteins.")
    }
    strength <- association^parameters$gamma
    if (base::diff(base::range(strength)) <= .Machine$double.eps) {
      distance <- base::rep(base::ceiling(parameters$n_bins / 2), base::length(strength))
    } else {
      bins <- base::cut(strength,
                        breaks = base::seq(base::min(strength), base::max(strength),
                                           length.out = parameters$n_bins + 1L),
                        include.lowest = TRUE, labels = FALSE)
      distance <- parameters$n_bins - base::as.integer(bins) + 1L
    }
    graph <- igraph::set_edge_attr(backbone, "weight", value = distance)
    graph
  }

  calculate_nodewise_LFDs <- function(graph, use_weights, invert_weights) {
    comp <- igraph::components(graph)
    giant <- igraph::induced_subgraph(graph,
      base::which(comp$membership == base::which.max(comp$csize)))
    weight_vector <- NA
    if (use_weights) {
      weight_vector <- igraph::edge_attr(giant, "weight")
      if (invert_weights) weight_vector <- 1 / weight_vector
    }
    nfd <- r2 <- base::rep(NA_real_, igraph::vcount(giant))
    for (node_id in base::seq_len(igraph::vcount(giant))) {
      spl <- igraph::distances(giant, v = node_id, mode = "all", weights = weight_vector)
      grow <- base::sort(spl[base::is.finite(spl) & spl > 0])
      counts <- base::table(grow)
      radii <- base::as.numeric(base::names(counts))
      cumulative <- base::cumsum(base::as.numeric(counts))
      if (base::length(radii) >= 2L) {
        x <- base::log(radii)
        y <- base::log(cumulative)
        fit <- stats::lm(y ~ x)
        nfd[node_id] <- stats::coef(fit)[2L]
        r2[node_id] <- base::summary(fit)$r.squared
      }
    }
    base::data.frame(node = igraph::vertex_attr(giant, "name"), nfd = nfd, r2 = r2)
  }

  extract_components <- function(contrast, method_dir) {
    component_dir <- base::file.path(method_dir, "components")
    base::dir.create(component_dir, showWarnings = FALSE)
    sig_nodes <- contrast$node[base::which(contrast$p_value < parameters$p_cutoff)]
    sig_graph <- igraph::induced_subgraph(backbone,
      base::which(genes %in% sig_nodes))
    idx <- base::match(igraph::vertex_attr(sig_graph, "name"), contrast$node)
    sig_graph <- igraph::set_vertex_attr(sig_graph, "p_value", value = contrast$p_value[idx])
    sig_graph <- igraph::set_vertex_attr(sig_graph, "deltaLFD_modulated",
                                        value = contrast$delta_lfd[idx])
    igraph::write_graph(sig_graph, base::file.path(component_dir, "significant.graphml"),
                        format = "graphml")
    comp <- igraph::components(sig_graph)
    membership <- base::data.frame(
      node = igraph::vertex_attr(sig_graph, "name"),
      component = base::as.integer(comp$membership),
      component_size = base::as.integer(comp$csize[comp$membership])
    )
    membership$retained <- membership$component_size >= parameters$min_component_size
    retained <- base::which(comp$csize >= parameters$min_component_size)
    membership$exported_component <- base::match(membership$component, retained)
    utils::write.csv(membership, base::file.path(component_dir, "membership.csv"),
                     row.names = FALSE)
    if (igraph::ecount(sig_graph) > 0L) {
      giant <- igraph::induced_subgraph(sig_graph,
        base::which(comp$membership == base::which.max(comp$csize)))
      igraph::write_graph(giant, base::file.path(component_dir, "GIANT.graphml"),
                          format = "graphml")
    }
    for (j in base::seq_along(retained)) {
      graph <- igraph::induced_subgraph(sig_graph,
        base::which(comp$membership == retained[j]))
      component_nodes <- igraph::vertex_attr(graph, "name")
      table <- contrast[base::match(component_nodes, contrast$node), , drop = FALSE]
      table$component <- j
      stem <- base::file.path(component_dir, base::paste0("Component_", j))
      igraph::write_graph(graph, base::paste0(stem, ".graphml"), format = "graphml")
      utils::write.csv(table, base::paste0(stem, ".csv"), row.names = FALSE)
    }
    base::data.frame(significant_nodes = igraph::vcount(sig_graph),
                     significant_edges = igraph::ecount(sig_graph),
                     retained_components = base::length(retained))
  }

  base::message("Calculating structural LFD on ", igraph::vcount(backbone), " proteins.")
  structural <- calculate_nodewise_LFDs(backbone, use_weights = FALSE, invert_weights = FALSE)
  utils::write.csv(structural, base::file.path(output_dir, "LFD_structural.csv"),
                   row.names = FALSE)
  summary_rows <- base::list()
  for (method in parameters$methods) {
    base::message("Processing ", method)
    method_dir <- base::file.path(output_dir, method)
    base::dir.create(method_dir, showWarnings = FALSE)
    lfd <- base::list()
    for (condition in base::names(groups)) {
      graph <- build_weighted_graph(condition, method)
      igraph::write_graph(graph,
        base::file.path(method_dir, base::paste0("weighted_", condition, ".graphml")),
        format = "graphml")
      table <- calculate_nodewise_LFDs(graph, parameters$use_weights,
                                       parameters$invert_weights)
      table$lfd_structural <- structural$nfd[base::match(table$node, structural$node)]
      table$lfd_modulated <- table$nfd - table$lfd_structural
      lfd[[condition]] <- table
      utils::write.csv(table,
        base::file.path(method_dir, base::paste0("LFD_", condition, ".csv")), row.names = FALSE)
    }
    contrast <- base::data.frame(
      node = lfd$WT$node, WT = lfd$WT$lfd_modulated,
      KO = lfd$KO$lfd_modulated[base::match(lfd$WT$node, lfd$KO$node)]
    )
    contrast$delta_lfd <- contrast$KO - contrast$WT
    contrast$z_score <- (contrast$delta_lfd - base::mean(contrast$delta_lfd, na.rm = TRUE)) /
      (stats::sd(contrast$delta_lfd, na.rm = TRUE) + parameters$contrast_sd_epsilon)
    contrast$p_value <- 2 * stats::pnorm(-base::abs(contrast$z_score))
    utils::write.csv(contrast, base::file.path(method_dir, "delta_LFD_KO_vs_WT.csv"),
                     row.names = FALSE)
    status <- extract_components(contrast, method_dir)
    status$method <- method
    summary_rows[[method]] <- status[, c("method", "significant_nodes",
                                         "significant_edges", "retained_components")]
  }
  summary <- base::do.call(base::rbind, summary_rows)
  utils::write.csv(summary, base::file.path(output_dir, "component_summary.csv"),
                   row.names = FALSE)
  if (base::any(summary$retained_components < 1L)) {
    base::stop("Demonstration check failed: each method must yield a component of at least five nodes.")
  }
  base::message("DIFFRACT complete. Results: ", base::normalizePath(output_dir))
  base::invisible(summary)
}

base::local({
  script_path <- NULL
  for (frame in base::sys.frames()) {
    if (!base::is.null(frame$ofile)) script_path <- frame$ofile
  }
  if (base::is.null(script_path)) {
    arguments <- base::commandArgs(trailingOnly = FALSE)
    script_arg <- base::grep("^--file=", arguments, value = TRUE)
    if (base::length(script_arg)) script_path <- base::sub("^--file=", "", script_arg[1L])
  }
  script_dir <- if (base::is.null(script_path)) base::getwd() else
    base::dirname(base::normalizePath(script_path, mustWork = TRUE))
  diffract(base::file.path(script_dir, "example_data"),
               base::file.path(script_dir, "output"))
})
