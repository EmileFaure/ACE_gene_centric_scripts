# Function adapted from the ggtriplotRDA.R from Francois Gillet, Daniel Borcard & Pierre Legendre, 2021
# to implement colors and shape for the sites
#
# A function to draw a triplot (scaling 1 or scaling 2) from an object
# of class "rda" containing RDA result from vegan's rda() function.
#
# This new version can handle RDA results with a single canonical axis 
# or with any meaningful combination of canonical and/or residual axes,
# as in vegan::plot.cca() function. 
# The following packages are required: vegan, tidyverse, ggrepel.
# 
# ARGUMENTS
#
# ##### General parameters
# res.rda          An rda{vegan} object.
# ax1, ax2         Axes to be drawn as abscissa and ordinate. Defaults: 1 and 2.
# site.sc          Can be set to "lc" (linear constraints or model scores, default)
#                  or "wa" (weighted averages, default in vegan).
# scaling          Scaling type: only 1 or 2 are supported. Default: 2.
#
# ##### Items to be plotted
# plot.sites       If TRUE, the sites will be plotted as small circles.
# plot.spe         If TRUE, the species (or other response variables) will be plotted.
# plot.env         If TRUE, arrows for the explanatory variables will be plotted.
# plot.centr       If TRUE, symbols will be plotted at the centroids of factor levels.
# arrows.only      if TRUE, plot arrows for quant. explanatory var. and factor classes
# label.sites      If TRUE, labels are added to the site symbols.
# label.spe        If TRUE, labels are added to the species arrows.
# label.env        If TRUE, labels are added to the environmental variable arrows.
# label.centr      If TRUE, labels are added to the centroids of factor levels.
#
# ##### Multipliers, selection of species to be plotted
# mult.spe         Multiplier for length of the species arrows. Default: 1.
# mult.arrow       Multiplier for length of the environmental arrows. Default: 1.
# select.spe       Vector containing a selection of the species numbers to be drawn in
#                  the biplot, e.g. c(1,2,5,8,12). Draw all species if select.spe=NULL
#                  (default value). The species that are well represented in the RDA plot
#                  can be identified using goodness(RDA.output.object,display="species")
#
# ##### Position of the plot in frame, margins
# mar.percent      Factor to expand plot size to accommodate all items and labels. Positive
#                  values increase the margins around the plot, negative values reduce
#                  them.
# optimum          If TRUE, the longest species and environmental arrows are stretched to
#                  a length equal to the distance to the origin of the site farthest from
#                  the origin in the plot of (ax1, ax2). This is an optimal combined
#                  representation of the three elements. The lengths of the species and
#                  environmental arrows can be further modified using the arguments
#                  mult.spe and mult.arrow.
# move.origin      Move plot origin right-left and up-down. Default: move.origin=c(0,0).
#                  Ex. move.origin=c(-1,0.5) moves origin by 1 unit left and 0.5 unit up.
#
# ##### Variables
# silent           If FALSE, intermediate computation steps are printed. Default: TRUE.
#
# License: GPL-2
# Authors: Francois Gillet, Daniel Borcard & Pierre Legendre, 2021
# Adapted by Lisa-Marie Delpech


ggRDA <-
  function(res.rda,
           ax1 = 1,
           ax2 = 2,
           site.sc = "lc",
           scaling = 2,
           plot.sites = TRUE,
           plot.spe = TRUE,
           plot.env = TRUE,
           plot.centr = TRUE,
           arrows.only = FALSE,
           arrow.spe = TRUE,
           label.sites = TRUE,
           label.spe = TRUE,
           label.env = TRUE,
           label.centr = TRUE,
           mult.spe = 1,
           mult.arrow = 1,
           select.spe = NULL,
           mar.percent = 0.15,
           optimum = TRUE,
           move.origin = c(0, 0),
           silent = TRUE,
           #large.data = FALSE, # Add this option to plot species in background first if there are many
           vec.col = NULL, # Change to site.col
           vec.shape = NULL, # Change to site.shape
           spe.col = NULL,
           env.names = function(x){x}, # A function to rename variables
           cen.names = function(x){x},
           return.data = FALSE) { 
    
    require(tidyverse)
    require(ggrepel)
    
    ### Internal functions ----
    
    stretch <-
      function(sites, mat, ax1, ax2, n, silent = silent) {
        # Compute stretching factor for the species or environmental arrows
        # First, compute the longest distance to centroid for the sites
        D <- sqrt(sites[,ax1]^2 + sites[,ax2]^2)
        target <- max(D)
        # Then, compute the longest distance to centroid for the species or 
        # environmental arrows
        if (inherits(mat, what = "matrix")) {
          D <- sqrt(mat[,ax1]^2 + mat[,ax2]^2)
          longest <- max(D)
        } else {
          tmp2 <- rbind(c(0, 0), mat[c(ax1, ax2)])
          longest <- dist(tmp2)
          # print(tmp2)
        }  # If a single row left in 'mat'
        #
        if (!silent)
          cat("target =",
              target,
              " longest =",
              longest,
              " fact =",
              target / longest,
              "\n")
        fact <- target / longest
      }
    
    larger.plot <-
      function(sit.sc,
               spe.sc,
               BP.sc,
               percent,
               move.origin,
               ax1,
               ax2) {
        # Internal function to expand plot limits 
        # (adapted from code by Pierre Legendre)
        mat <- rbind(sit.sc, spe.sc, BP.sc)
        range.mat <- apply(mat, 2, range)
        rownames(range.mat) <- c("Min", "Max")
        z <- apply(range.mat, 2, function(x)
          x[2] - x[1])
        range.mat[1, ] <- range.mat[1, ] - z * percent
        range.mat[2, ] <- range.mat[2, ] + z * percent
        if (move.origin[1] != 0)
          range.mat[, ax1] <- range.mat[, ax1] - move.origin[1]
        if (move.origin[2] != 0)
          range.mat[, ax2] <- range.mat[, ax2] - move.origin[2]
        range.mat
      }
    ### End internal functions
    
    if (class(res.rda)[1] != "rda" & class(res.rda)[2] != "rda")
      stop("The input file is not a vegan rda output object")
    if (length(res.rda$colsum) == 1)
      stop("Function triplot.rda is not compatible with results that contain no species scores")
    if (scaling != 1 & scaling != 2)
      stop("Function only available for scaling = 1 or 2")
    if (site.sc == "lc") {
      cat("\n-----------------------------------------------------------------------")
      cat("\nSite constraints (lc) selected. To obtain site scores that are weighted")
      cat("\nsums of species scores (default in vegan), argument site.sc must be set")
      cat("\nto wa.")
      cat("\n-----------------------------------------------------------------------\n")
    }
    
    k <- length(res.rda$CCA$eig)        # number of RDA eigenvalues
    n.sp <- length(res.rda$colsum)      # number of species
    # 'vec' will contain the selection of species to be drawn
    if (is.null(select.spe)) {
      vec <- 1:n.sp } else {
        vec <- select.spe
      }
    
    # Scaling 1: the species scores have norms of 1
    # Scaling 1: the site scores are scaled to variances = can.eigenvalues
    # Scaling 2: the species scores have norms of sqrt(can.eigenvalues)
    # Scaling 2: the site scores are scaled to variances of 1
    
    # This version reconstructs and uses the original RDA output of L&L 2012, 
    # Section 11.1.3
    
    Tot.var <- res.rda$tot.chi            # Total variance in response data Y
    eig.val <- c(res.rda$CCA$eig, res.rda$CA$eig) # all eigenvalues
    Lambda <- diag(eig.val)               # Diagonal matrix of eigenvalues
    eig.val.rel <- eig.val / Tot.var      # Relative eigenvalues of Y-hat
    Diag <- diag(sqrt(eig.val.rel))       # Diagonal matrix of sqrt(relative eigenvalues)
    U.sc1 <- cbind(res.rda$CCA$v, res.rda$CA$v) # All species scores, scaling=1
    U.sc2 <- U.sc1 %*% sqrt(Lambda)       # Species scores, scaling=2
    colnames(U.sc2) <- colnames(U.sc1)
    n <- nrow(res.rda$CCA$u)              # Number of observations
    Z.sc2 <- cbind(res.rda$CCA$u, res.rda$CA$u) * sqrt(n - 1)  # "lc" site scores, scaling=2
    Z.sc1 <- Z.sc2 %*% sqrt(Lambda)       # "lc" site scores, scaling=1
    colnames(Z.sc1) <- colnames(Z.sc2)
    F.sc2 <- cbind(res.rda$CCA$wa, res.rda$CA$u) * sqrt(n - 1)  # "wa" site scores, scaling=2
    F.sc1 <- F.sc2 %*% sqrt(Lambda)       # "wa" site scores, scaling=1
    colnames(F.sc1) <- colnames(F.sc2)
    BP.sc2 <- res.rda$CCA$biplot          # Biplot scores, scaling=2 ; cor(Z.sc1, X)
    BP.sc2 <- cbind(BP.sc2, matrix(0, nrow = nrow(BP.sc2), 
                                   ncol = length(eig.val) - k))
    colnames(BP.sc2) <- colnames(F.sc2)
    BP.sc1 <- BP.sc2 %*% Diag             # Biplot scores, scaling=1
    colnames(BP.sc1) <- colnames(BP.sc2)
    
    if (!is.null(res.rda$CCA$centroids)) {
      centroids.sc2 <- res.rda$CCA$centroids * sqrt(n - 1) # Centroids, scaling=2
      centroids.sc2 <- cbind(centroids.sc2, matrix(0, nrow = nrow(centroids.sc2), 
                                                   ncol = length(eig.val) - k))
      colnames(centroids.sc2) <- colnames(F.sc2)
      centroids.sc1 <- centroids.sc2 %*% sqrt(Lambda)      # Centroids, scaling=1
      colnames(centroids.sc1) <- colnames(centroids.sc2)
    }
    centroids.present <- TRUE
    if (is.null(res.rda$CCA$centroids)) {
      centroids.present <- FALSE
      if (plot.centr | label.centr) {
        cat("\nNo factor, hence levels cannot be plotted with symbols;")
        cat("\n'plot.centr' is set to FALSE\n")
        plot.centr  <- FALSE
        label.centr <- FALSE
      }
    }
    
    if (is.null(select.spe)) {vec <- 1:n.sp} else {vec <- select.spe}
    
    if (scaling == 1) {
      if (site.sc == "lc") {
        sit.sc <- Z.sc1
      } else {
        sit.sc <- F.sc1
      }
      spe.sc <- U.sc1[vec, ]
      BP.sc  <- BP.sc1
      if (centroids.present)
        centroids <- centroids.sc1
    } else {
      # For scaling 2
      if (site.sc == "lc") {
        sit.sc <- Z.sc2
      } else {
        sit.sc <- F.sc2
      }
      spe.sc <- U.sc2[vec, ]
      BP.sc  <- BP.sc2
      if (centroids.present)
        centroids <- centroids.sc2
    }
    
    fact.spe <- 1
    fact.env <- 1
    if (centroids.present & (plot.centr | label.centr)) {
      to.plot <- which(!(rownames(BP.sc) %in% rownames(centroids)))
    } else {
      to.plot <- 1:nrow(BP.sc)
    }
    
    if (optimum) {
      if (plot.spe | label.spe){
        fact.spe <-
          stretch(sit.sc, spe.sc, ax1, ax2, n, silent = silent) # This requires computing dist() for a large matrix if there are many species
        # In that case plot.spe is usually set to FALSE and this step is not needed.
        # If it is set to TRUE, it will return Error: vector memory exhausted (limit reached?) if matrix is too large
      }
      if (arrows.only) {
        fact.env <-
          stretch(sit.sc, BP.sc, ax1, ax2, n, silent = silent)
      } else {
        # arrows only==FALSE
        quant.env.present <- FALSE
        if (length(to.plot) > 0) {
          quant.env.present <- TRUE
          fact.env <-
            stretch(sit.sc, BP.sc[to.plot, ], ax1, ax2, n, silent = silent)
        }
      }
    }
    
    if (!silent)
      cat("fact.spe =", fact.spe, "   fact.env =", fact.env, "\n")
    spe.sc <- spe.sc * fact.spe * mult.spe
    BP.sc <- BP.sc * fact.env * mult.arrow
    lim <-
      larger.plot(
        sit.sc[, ],
        spe.sc[, ],
        BP.sc[, ],
        percent = mar.percent,
        move.origin = move.origin,
        ax1 = ax1,
        ax2 = ax2
      )
    if (!silent)
      print(lim)
    
    
    ### Drawing the triplot begins here ----
    
    # Draw the main plot
    mat <- rbind.data.frame(sit.sc, spe.sc, BP.sc)
    sitsc <- as.data.frame(sit.sc) %>% select(all_of(ax1), all_of(ax2))
    # names(sitsc) <- c("ax1", "ax2")
    axn <- names(sitsc)
    
    if(ax1 <= k){
      titre <- "RDA triplot - Scaling"
    } else{
      titre <- "Biplot of residuals of RDA - Scaling"
    }
    
    pp <- ggplot(sitsc, aes(x = .data[[axn[1]]], y = .data[[axn[2]]])) +
      labs(title = paste(titre, scaling, "-", site.sc),
           x = paste0(names(sitsc)[1], " (", 
                      round(100 * eig.val.rel[ax1], 1), "%)"), 
           y = paste0(names(sitsc)[2], " (", 
                      round(100 * eig.val.rel[ax2], 1), "%)")) +
      xlim(lim[1, ax1], lim[2, ax1]) +
      ylim(lim[1, ax2], lim[2, ax2]) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_vline(xintercept = 0, linetype = "dotted")  # +
    #coord_equal()
    
    # Draw the species scores
    spesc <- as.data.frame(spe.sc) %>% select(all_of(ax1), all_of(ax2))
    names(spesc) <- c("axis1", "axis2")
    
    if (plot.spe) {
      if (arrow.spe) {
        pp <- pp + geom_segment(
          data = spesc,
          aes(
            x = 0,
            xend = axis1,
            y = 0,
            yend = axis2
          ),
          colour = "#787878", # "red"
          arrow = arrow(length = unit(0.01, "npc"))
        )
      } else {
        if (is.null(spe.col)){
          pp <- pp + geom_point(
            data = spesc,
            aes(
              x = axis1,
              y = axis2,
            ),
            shape = 4,
            size = 2,
            colour = "#787878",
            alpha = 0.3
          )
        } else {
          pp <- pp + geom_point(
            data = spesc,
            aes(
              x = axis1,
              y = axis2,
              col = spe.col
            ),
            shape = 4,
            size = 2,
          )
        }
      }
      if (label.spe)
        pp <- pp + geom_text_repel(
          data = spesc,
          aes(x = axis1, y = axis2, label = rownames(spesc)),
          size = 3.5,
          fontface = "italic",
          colour = "#787878")
    } else {
      if (label.spe)
        pp <- pp + geom_text(
          data = spesc,
          aes(x = axis1, y = axis2, label = rownames(spesc)),
          size = 3.5,
          fontface = "italic",
          colour = "#787878")
    }
    
    # Draw the site scores ("lc" or "wa")
    if (plot.sites) {
      pp <- pp + geom_point(aes(col = vec.col, shape = vec.shape),
                            size = 3) # let the user chose aesthetics? # alpha = 0.7
      if (label.sites)
        pp <- pp + geom_text_repel(aes(label = rownames(sitsc)), size = 3)
    } else {
      if (label.sites)
        pp <- pp + geom_text(aes(label = rownames(sitsc)), size = 3)
    }
    
    # Draw the explanatory variables (only if at least ax1 is canonical)
    if(ax1 <= k){
      BPsc <- as.data.frame(BP.sc) %>% select(all_of(ax1), all_of(ax2))
      names(BPsc) <- c("axis1", "axis2")
      
      if (!arrows.only) {
        # 1. Quantitative variables
        if (quant.env.present & plot.env) {
          # Print arrows and labels for quantitative variables
          pp <- pp + geom_segment(
            data = BPsc[to.plot, ],
            aes(
              x = 0,
              xend = axis1 * mult.arrow,
              y = 0,
              yend = axis2 * mult.arrow
            ),
            colour = "#B70034", # "blue"
            # size = 1,
            arrow = arrow(length = unit(0.01, "npc"))
          )
          if (label.env)
            # Print labels for the quantitative variables
            pp <- pp + geom_text_repel(
              data = BPsc[to.plot, ],
              aes(x = axis1, y = axis2, label = env.names(rownames(BPsc)[to.plot])),
              colour = "#B70034") # "blue"
        } else {
          if (quant.env.present & !plot.env & label.env)
            # Only print labels for quantitative variables
            pp <- pp + geom_text(
              data = BPsc[to.plot, ],
              aes(x = axis1, y = axis2, label = env.names(rownames(BPsc)[to.plot])),
              colour = "#B70034") # "blue"
        }
        
        # 2. Centroids and labels of factor levels
        censc <- as.data.frame(centroids) %>% select(all_of(ax1), all_of(ax2))
        names(censc) <- c("axis1", "axis2")
        
        if (centroids.present & plot.centr) {
          # Print symbols and labels for factor classes
          pp <- pp + geom_point(
            data = censc,
            aes(axis1, axis2),
            shape = 1,
            size = 2,
            #alpha = 0.7,
            colour = "black" # purple # #B70034
          )
          if (label.centr)
            pp <- pp + geom_text_repel(
              data = censc,
              aes(x = axis1, y = axis2, label = cen.names(rownames(censc))),
              #fontface = "bold",
              colour = "black") # #B70034
        } else {
          if (centroids.present & !plot.centr & label.centr)
            # Only print labels for classes
            pp <- pp + geom_text(
              data = censc,
              aes(x = axis1, y = axis2, label = cen.names(rownames(censc))),
              colour = "blue")
        }
      }
      
      # 3. All env. var.: plot arrows and labels for all var. in 'BP.sc', quant. and factors
      if (arrows.only) {
        pp <- pp + geom_segment(
          data = BPsc,
          aes(
            x = 0,
            xend = axis1 * mult.arrow,
            y = 0,
            yend = axis2 * mult.arrow
          ),
          colour = "blue",
          # size = 1,
          arrow = arrow(length = unit(0.01, "npc"))
        )
        if (label.env)
          # Print labels for the quantitative variables
          pp <- pp + geom_text_repel(
            data = BPsc,
            aes(x = axis1, y = axis2, label = env.names(rownames(BPsc))),
            colour = "blue")
      }
    }
    
    pp <- pp + theme_linedraw() +
      theme(panel.grid = element_blank())
    pp
    
    # output <- list(Sites = sitsc,
    #                Species = spesc,
    #                Biplot = BPsc,
    #                Centroids = censc,
    #                #Axes = c(ax1,ax2),
    #                Method = data.frame(Scaling = scaling,
    #                                    Site_scaling = site.sc,
    #                                    Stretch.spe = fact.spe,
    #                                    Stretch.env = fact.env),
    #                gg = pp
    # )
    # class(output) <- c("scaled rda sc", "ggplot", "list")
    # 
    # if (return.data){
    #   output
    # } else {
    #   output$gg
    #   }
  }

