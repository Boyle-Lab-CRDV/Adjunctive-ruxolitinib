###################################################### DO MIXED MODEL ##################################################
#
# Script to run linear mixed-effect model regression
# Written by  : Damian Oyong
# Date        : 28 March 2023
# Last update : 16 Sep 2025
# Description : This will calculate linear mixed effect model OR linear model (if random effects = 0)
#               To ensure the model fits into the assumptions, residual bootstrapping is added at 5000
#
# UPDATE: - Based on input and expert advice from Stacey Llewelllyn (QIMR Berghofer)
#         - Transform all your response variables into log to adhere to mixed model assumptions
#         - For data with zero value, add HALF of the smallest non-zero value prior to log transformation
#         - Use lmer instead of lme, AND with REML = TRUE, this tells the model to focus more on the fixed effects
#         - Sometimes you still get non-normal residuals in the model, therefore we apply Residual bootstrapping
#           to obtain consistent estimators of bias and SE
#         - Added predicted values based on bootstrapping
#         - Added if random effect is 0 or "singular", we analyse the data using simple linear regression model
#         - Added a function to get individual Fixed effect estimates, e.g. just RUX overtime
#           - IF contrast = TRUE, this will get contrast for the second level fixed2 variable
#         - Added do.mm2 function to analyse just 1 fixed effect
#         - Fix contrast = TRUE formatting when only two fixed1 variables available. The issue was "term" column is generated when only two variables present
#

#---------------------------------------------------- LIBRARY ----------------------------------------------------
library(lme4)
library(ggResidpanel)
library(lmeresampler)
library(dplyr)
library(tidyr)
library(gtools)
library(openxlsx)
library(cowplot)
library(boot)
library(emmeans)
library(tibble)
library(ggplot2)
library(ggpubr)
library(stringr)

#---------------------------------------------------- FUNCTION ----------------------------------------------------
# JUST RUN THIS FUNCTION
# DO NOT EDIT UNLESS KNOW WHAT YOU ARE DOING

mm.transform <- function(df_count, log = "log10") {
  
  zero.cols <- which(sapply(df_count, function(col) any(col == 0, na.rm = TRUE)))
  zero.cols <- names(df_count)[zero.cols]
  normal.cols <- setdiff(names(df_count), zero.cols)
  
  if (log == "log10") {
    if (length(zero.cols) == 0) {
      df_count[normal.cols] <- lapply(df_count[normal.cols], function(x) {
        x_transformed <- log10(x)
        x_transformed[is.na(x)] <- NA
        return(x_transformed)
      })
    } else {
      df_count[zero.cols] <- lapply(df_count[zero.cols], function(x) {
        x_transformed <- log10(x + (min(x[x != 0], na.rm = TRUE) / 2))
        x_transformed[is.na(x)] <- NA
        return(x_transformed)
      })
      df_count[normal.cols] <- lapply(df_count[normal.cols], function(x) {
        x_transformed <- log10(x)
        x_transformed[is.na(x)] <- NA
        return(x_transformed)
      })
    }
  }
  
  if (log == "log2") {
    if (length(zero.cols) == 0) {
      df_count[normal.cols] <- lapply(df_count[normal.cols], function(x) {
        x_transformed <- log2(x)
        x_transformed[is.na(x)] <- NA
        return(x_transformed)
      })
    } else {
      df_count[zero.cols] <- lapply(df_count[zero.cols], function(x) {
        x_transformed <- log2(x + (min(x[x != 0], na.rm = TRUE) / 2))
        x_transformed[is.na(x)] <- NA
        return(x_transformed)
      })
      df_count[normal.cols] <- lapply(df_count[normal.cols], function(x) {
        x_transformed <- log2(x)
        x_transformed[is.na(x)] <- NA
        return(x_transformed)
      })
    }
  }
  
  return(df_count)
}

mm.impute <- function(df_count){
  zero.cols <- df_count %>% 
    select(where(~ any(. == 0))) %>% 
    colnames()
  
  df_count[zero.cols] <- lapply(df_count[zero.cols], function(x) x + (min(x[x != 0])/2))
  
  return(df_count)
}

do.mm <- function(df_count, df_meta, fixed1, fixed2, random1, lm.opt = FALSE,
                  intercept = FALSE, contrast = FALSE, bts = 5000){
  set.seed(123) ## Pseudonumbering to reproduce consistent random numbers
  i = 1
  
  cat(crayon::green("Analysing data using linear mixed-effect model with residual bootstrapping\n"))
  cat(crayon::red("Will take awhile to run the analysis - be patient\n"))
  
  mm_pred <- list()
  
  suppressWarnings(dir.create(paste0("Residual Diagnostic")))
  residual_dir <- paste0("Residual Diagnostic/")
  
  for(i in 1:ncol(df_count)){
    
    cat(i, "----",colnames(df_count)[i],"\n")
    
    temp_comb <- df_meta %>%
      cbind(df_count[, i]) %>%
      setNames(c(colnames(df_meta), "Response")) %>%
      tibble::remove_rownames() %>%
      mutate(term = as.character(1:n()))
    head(temp_comb)
    
    mm.formula <- as.formula(paste0("Response~1+",fixed1,"*",fixed2,"+ (1|",random1,")"))
    
    mm = try(lmer(mm.formula, data = temp_comb, REML = TRUE))
    
    #### IF random intercept = 0, we use simple linear regression
    
    if(summary(mm)$varcor == 0 | isSingular(mm) | lm.opt == TRUE){
      cat(crayon::red("Random effects were ZERO, fitting regular linear model instead\n"))
      
      temp_comb <- temp_comb %>%
        dplyr::filter(., !is.na(!!rlang::sym(fixed1)) & !is.na(!!rlang::sym(fixed2)))
      
      lm.formula <- as.formula(paste0("Response~1+",fixed1,"*",fixed2)) # LM Formula
      lm = try(lm(lm.formula, data = temp_comb)) # Regular LM
      
      ### Residual Diagnostic
      p = suppressMessages({try(ggResidpanel::resid_panel(lm, smoother = TRUE, 
                                                          qqbands = TRUE, type = "pearson")+
                                  cowplot::draw_figure_label(label = paste0(colnames(df_count)[i]),
                                                             size = 12, fontface = "bold",
                                                             position = "top.right"))
      })
      save_plot(filename = paste0(residual_dir,
                                  colnames(df_count)[i],
                                  "_Residual Diagnostic Plot (LM).pdf"), plot = p)
      
      ### Residual Bootstrapping lm
      
      ## merging fitted and residuals from model into dataset
      temp_comb <- temp_comb %>%
        mutate(
          fit = fitted(lm)[match(term, names(fitted(lm)))],
          resid = resid(lm)[match(term, names(resid(lm)))]
        ) %>% 
        filter(., !is.na(Response))
      
      ## residual resampling
      stat <- function(temp_comb, indices){
        data.star  <- temp_comb %>% mutate(Response = fit + resid[indices])
        model.star <- lm(lm.formula, data = data.star)
        output     <- coef(model.star)
        return(output)
      }
      
      lm.boot <- boot(temp_comb, stat, R=bts)
      lm.est <- summary(lm.boot)
      rownames(lm.est) <- names(lm.boot$t0)
      lm.CI  <- confint(lm.boot, type="perc")
      
      boot.stats <- cbind(lm.est, lm.CI) %>%
        mutate(bootMean = original + bootBias) %>%
        rownames_to_column(., var = "term") %>%
        select(., c("term", "bootMean", "bootSE", "2.5 %", "97.5 %")) %>%
        rename(., term = "term", rep.mean = "bootMean", se = "bootSE" , lower = "2.5 %", upper = "97.5 %") %>%
        mutate(t_statistic = rep.mean / se) %>%
        mutate(pval = 2 * (1 - pt(abs(t_statistic), df = bts - 1)))
      
      # Create dataframe of results for Fixed1 and Fixed2
      if(intercept == FALSE){
        fix1.df <- dplyr::filter(boot.stats, grepl(paste0("^",fixed1), term) & 
                                   !grepl(paste0(fixed2), term)) %>% 
          mutate(term = gsub(paste0(fixed1), "", term),
                 fixed = fixed1)
        fix2.df <- dplyr::filter(boot.stats, grepl(paste0(fixed2), term) & 
                                   grepl(paste0(fixed1), term)) %>% 
          mutate(term = gsub(paste0("^",fixed1,"([^:]+):?.*"), "\\1", term),
                 fixed = fixed2)
        
      } else {
        fix1.df <- dplyr::filter(boot.stats, grepl(paste0("^",fixed1), term) & 
                                   !grepl(paste0(fixed2), term)) %>% 
          mutate(term = gsub(paste0(fixed1), "", term),
                 fixed = fixed1)
        fix2.df <- dplyr::filter(boot.stats, grepl(paste0(fixed2), term)) %>% 
          mutate(term = gsub(paste0("^",fixed1,"([^:]+):?.*"), "\\1", term),
                 fixed = fixed2) %>% 
          mutate(term = c(levels(df_meta[[fixed1]])[1], term[-1]))
      }
      
      # Combine both dataframes
      comb_res <- rbind(fix1.df, fix2.df) %>%
        mutate(Response = paste0(colnames(df_count)[i]))
      
      #### Obtaining predicted values
      stat1 <- function(temp_comb, indices){
        data.star  <- temp_comb %>% mutate(Response = fit + resid[indices])
        model.star <- lm(lm.formula, data=data.star)
        output     <- predict(model.star)
        return(output)
      }
      
      predval <- t(try(boot(temp_comb %>% filter(., !is.na(Response)), stat1, R = bts)$t))
      
      mm_pred[colnames(df_count)[i]] <- list(data.frame(predict = rowMeans(predval), se = apply(predval, 1, sd)) %>%
                                               cbind(., temp_comb) %>% 
                                               select(-c(term, Response))
      )
      
      
      if(contrast == TRUE){
        
        reference = levels(temp_comb[[fixed1]])[1] # Reference comparison
        lm.contrast_stat <- function(temp_comb, indices) {
          # Fit model to the bootstrapped sample
          lm.data_star  <- temp_comb %>% mutate(Response = fit + resid[indices])
          lm.model_star <- lm(lm.formula, data = lm.data_star)
          
          em.formula <- as.formula(paste0("~", fixed1,"*",fixed2))
          em_boot <- emmeans::emmeans(lm.model_star, em.formula, data = lm.data_star)
          
          comparisons <- list()
          for (a in 1:length(levels(summary(em_boot)[[fixed1]]))) {
            vec <- rep(0, length(summary(em_boot)[[fixed1]]))
            vec[length(levels(summary(em_boot)[[fixed1]])) + a] <- 1
            comparisons[[a]] <- vec
          }
          
          contrast_compare <- list()
          for (a in 2:length(comparisons)) {
            # Calculate the difference between comparisons[[i]] and comparisons[[1]]
            contrast_compare[[letters[a-1]]] <- comparisons[[a]] - comparisons[[1]]
          }
          
          names <- c()
          for (a in 2:length(levels(em_boot)[[fixed1]])) {
            name <- paste0(levels(em_boot)[[fixed2]][2], "_", levels(em_boot)[[fixed1]][a], "v", reference)
            names <- c(names, name)
          }
          
          names(contrast_compare) <- names
          
          boot_contrasts <- contrast(em_boot, method = contrast_compare)
          est_values <- summary(boot_contrasts)[["estimate"]]
          names(est_values) <- summary(boot_contrasts)$contrast
          return(est_values)
        }
        
        lm_res <- boot(temp_comb, lm.contrast_stat, R = bts)
        # lm_res$t[,1]
        # Calculate the mean, SE, t-stats and CI of the bootstrap results
        boot_contrast <- data.frame(rep.mean = apply(lm_res$t, 2, mean))
        boot_contrast_se <- data.frame(se =
                                         apply(lm_res$t, 2, function(x) sd(x) / sqrt(length(x))))
        boot_contrast_ci <- t(as.data.frame(apply(lm_res$t, 2, function(x) quantile(x, probs = c(0.025, 0.975)))))
        boot_contrast_t <- data.frame(t_statistic =
                                        apply(lm_res$t, 2, function(x) mean(x) / (sd(x) / sqrt(length(x)))))
        
        contrast.results <- cbind(boot_contrast, boot_contrast_se, boot_contrast_t, boot_contrast_ci)
        
        ## Adding contrast names
        names <- c()
        # drop unused levels
        lev1 <- levels(droplevels(temp_comb[[fixed1]]))
        for (a in seq(2, length(lev1))) {
          name <- paste0(levels(temp_comb[[fixed2]])[2], "_", lev1[a], "v", reference)
          names <- c(names, name)
        }
        rownames(contrast.results) <- names
        
        # Compute p-values ####
        ## FROM HERE: https://www.r-bloggers.com/2019/05/bootstraping-follow-up-contrasts-for-within-subject-anovas/
        
        boot_pvalues <- function(x, side = c(0, -1, 1)) {
          side <- side[1]
          x <- as.data.frame(x$t)
          ##LOOPS OVER EACH COLUMN IN DATA FRAME
          ps <- sapply(x, function(.x) {
            s <- na.omit(.x)
            s0 <- 0  ## reference value for comparison
            N <- length(s)  ## number of bootstrapped estimates (after removing NA)
            # Calculates the minimum of the p-values for both tails and 
            # then multiplies by 2 to adjust for the two-sided test. 
            # The calculation involves counting how many bootstrap estimates are greater than or equal to, 
            # and less than or equal to, the reference value (s0), respectively.
            if (side == 0) {
              min((1 + sum(s >= s0)) / (N + 1),
                  (1 + sum(s <= s0)) / (N + 1)) * 2
            } else if (side < 0) {
              (1 + sum(s <= s0)) / (N + 1)
            } else if (side > 0) {
              (1 + sum(s >= s0)) / (N + 1)
            }
          })
          setNames(ps,colnames(x))
        }
        bootp <- boot_pvalues(lm_res, side=c(0))  ## two-sided.
        
        # Holm correction for multiple comparisons
        p_adjusted_holm <- p.adjust(bootp, method = "holm")
        
        contrast.results <- cbind(contrast.results, as.data.frame(bootp) ,as.data.frame(p_adjusted_holm)) %>% 
          rownames_to_column(., var = "term") %>% 
          mutate(Response = paste0(colnames(df_count)[i])) %>% 
          rename(pvalue = "bootp", p_adj_holm = "p_adjusted_holm", lower = "2.5%", upper = "97.5%")
        
        # mm_contrast[colnames(df_count)[i]] <- list(bootresults)
      }
      rm(lm, boot.stats, lm.est, lm.CI)
    } else {
      
      ### Residual Diagnostic if lmer
      p = suppressMessages({try(ggResidpanel::resid_panel(mm, smoother = TRUE, 
                                                          qqbands = TRUE, type = "pearson")+
                                  cowplot::draw_figure_label(label = paste0(colnames(df_count)[i]),
                                                             size = 12, fontface = "bold",
                                                             position = "top.right"))
      })
      save_plot(filename = paste0(residual_dir,
                                  colnames(df_count)[i],
                                  "_Residual Diagnostic Plot (LMER).pdf"), plot = p)
      
      ### Residual Bootstrapping lmer
      mm.boot <- try(resid_bootstrap(mm, .f = fixef,  B = bts))
      boot.stats <- as.data.frame(mm.boot$stats) %>%
        select(-c("observed","bias")) %>%
        mutate(t_statistic = rep.mean/se) %>%
        mutate(pval = 2 * (1 - pt(abs(t_statistic), df = bts - 1)))
      CI <- as.data.frame( confint(mm.boot, type="perc")) ## using Percentile confidence intervals
      boot.stats <- left_join(boot.stats, CI[,c(1,3,4)], by="term")
      
      if(intercept == FALSE){
        fix1.df <- dplyr::filter(boot.stats, grepl(paste0("^",fixed1), term) & 
                                   !grepl(paste0(fixed2), term)) %>% 
          mutate(term = gsub(paste0(fixed1), "", term),
                 fixed = fixed1)
        fix2.df <- dplyr::filter(boot.stats, grepl(paste0(fixed2), term) & 
                                   grepl(paste0(fixed1), term)) %>% 
          mutate(term = gsub(paste0("^",fixed1,"([^:]+):?.*"), "\\1", term),
                 fixed = fixed2)
        
      } else {
        fix1.df <- dplyr::filter(boot.stats, grepl(paste0("^",fixed1), term) & 
                                   !grepl(paste0(fixed2), term)) %>% 
          mutate(term = gsub(paste0(fixed1), "", term),
                 fixed = fixed1)
        fix2.df <- dplyr::filter(boot.stats, grepl(paste0(fixed2), term)) %>% 
          mutate(term = gsub(paste0("^",fixed1,"([^:]+):?.*"), "\\1", term),
                 fixed = fixed2) %>% 
          mutate(term = c(levels(df_meta[[fixed1]])[1], term[-1]))
      }
      
      comb_res <- rbind(fix1.df, fix2.df) %>%
        mutate(Response = paste0(colnames(df_count)[i]))
      
      #### Obtaining predicted values
      set.seed(123)
      mm_pred[colnames(df_count)[i]] <- list(as.data.frame(try(resid_bootstrap(mm, .f = predict,  B = bts)$stats)) %>%
                                               right_join(., temp_comb, by = "term") %>%
                                               tidyr::drop_na(rep.mean) %>%
                                               select(-c(term, Response)) %>%
                                               rename(predict = "rep.mean")
      )
      
      if(contrast == TRUE){
        
        reference = levels(temp_comb[[fixed1]])[1]
        
        func.contrast <- function(mod){
          em.formula <- as.formula(paste0("~", fixed1,"*",fixed2))
          em <- emmeans::emmeans(mod, em.formula) # Get the expected mean
          
          comparisons <- list()
          for (a in 1:length(levels(summary(em)[[fixed1]]))) {
            vec <- rep(0, length(summary(em)[[fixed1]]))
            vec[length(levels(summary(em)[[fixed1]])) + a] <- 1
            comparisons[[a]] <- vec
          }
          
          contrast_compare <- list()
          for (a in 2:length(comparisons)) {
            contrast_compare[[letters[a-1]]] <- comparisons[[a]] - comparisons[[1]]
          }
          
          names <- c()
          for (a in 2:length(levels(em)[[fixed1]])) {
            name <- paste0(levels(em)[[fixed2]][2], "_", levels(em)[[fixed1]][a], "v", reference)
            names <- c(names, name)
          }
          names(contrast_compare) <- names
          
          c_ <- contrast(em, method = contrast_compare)
          
          est_values <- summary(c_)$estimate
          # names(est_values) <- mgsub::mgsub(paste(summary(c_)$Treatment, summary(c_)$contrast, sep = "_"), 
          #                                   c(" ","\\(","\\)","-"), c("","","","v"))
          names(est_values) <- summary(c_)$contrast
          return(est_values)
        }
        
        mm.boot_contrast <- try(resid_bootstrap(mm, .f = func.contrast,  B = bts))
        contrast.results <- as.data.frame(mm.boot_contrast$stats) %>% 
          select(-c("observed","bias")) %>%
          { if (!"term" %in% names(.)) mutate(., term = confint(mm.boot_contrast, type="perc")$term) else . } %>% 
          mutate(t_statistic = rep.mean/se) %>%
          mutate(pvalue = 2 * (1 - pt(abs(t_statistic), df = bts - 1)))
        
        ##bootstrap pvalues and 95% ci
        CI_contrast <- as.data.frame(confint(mm.boot_contrast, type="perc")) ## using Percentile confidence intervals
        contrast.results <- left_join(contrast.results, CI_contrast[,c(1,3,4)], by="term") %>% 
          mutate(p_adj_holm = p.adjust(pvalue, method = "holm"), # Holm correction for multiple comparisons
                 Response = paste0(colnames(df_count)[i]))  
        
        # mm_contrast[colnames(df_count)[i]] <- list(mm.contrast_results)
      }
      rm(mm, boot.stats, CI)
    }
    
    if(contrast == TRUE){
      if(i == 1){
        mm_res = comb_res
        mm_contrast = contrast.results
      }else{
        mm_res = rbind(mm_res,comb_res)
        mm_contrast = rbind(mm_contrast, contrast.results)
      }
      
    } else {
      if(i == 1){
        mm_res = comb_res
      }else{
        mm_res = rbind(mm_res,comb_res)
      }  
    }
  }
  
  
  rm(comb_res, temp_comb)
  
  rownames(mm_res) = NULL
  
  if(contrast == TRUE){
    mm_res <- mm_res %>%
      dplyr::group_by(term, fixed) %>%
      mutate(fdr = p.adjust(pval, method = "fdr")) %>% 
      ungroup()
    mm_contrast <- mm_contrast %>%
      dplyr::group_by(term) %>%
      mutate(fdr = p.adjust(pvalue, method = "fdr")) %>% 
      ungroup()
  } else {
    mm_res <- mm_res %>%
      dplyr::group_by(term, fixed) %>%
      mutate(fdr = p.adjust(pval, method = "fdr")) %>% 
      ungroup()
  }
  
  # return(mm_res)
  if(contrast == TRUE){
    return(list(mm_res = mm_res, mm_pred = mm_pred, mm_contrast = mm_contrast))
  }else{
    return(list(mm_res = mm_res, mm_pred = mm_pred))
  }
}


do.mm2 <- function(df_count, df_meta, fixed1, random1, bts = 5000, lm.opt = FALSE){
  set.seed(123) ## Pseudonumbering to reproduce consistent random numbers
  i = 1
  
  cat(crayon::green("Analysing data using linear mixed-effect model with residual bootstrapping\n"))
  cat(crayon::red("Will take awhile to run the analysis - be patient\n"))
  
  mm_pred <- list()
  
  suppressWarnings(dir.create(paste0("Residual Diagnostic")))
  residual_dir <- paste0("Residual Diagnostic/")
  
  
  for(i in 1:ncol(df_count)){
    
    cat(i, "----",colnames(df_count)[i],"\n")
    
    temp_comb <- df_meta %>%
      cbind(df_count[, i]) %>%
      setNames(c(colnames(df_meta), "Response")) %>%
      tibble::remove_rownames() %>%
      mutate(term = as.character(1:n()))
    head(temp_comb)
    
    mm.formula <- as.formula(paste0("Response~1+",fixed1,"+ (1|",random1,")"))
    
    mm = try(lmer(mm.formula, data = temp_comb, REML = TRUE))
    
    #### IF random intercept = 0, we use simple linear regression
    
    if(summary(mm)$varcor == 0 | isSingular(mm) | lm.opt == TRUE){
      cat(crayon::red("Random effects were ZERO, fitting regular linear model instead\n"))
      
      temp_comb <- temp_comb %>%
        dplyr::filter(., !is.na(!!rlang::sym(fixed1)))
      
      lm.formula <- as.formula(paste0("Response~1+",fixed1)) # LM Formula
      lm = try(lm(lm.formula, data = temp_comb)) # Regular LM
      
      ### Residual Diagnostic
      p = suppressMessages({try(ggResidpanel::resid_panel(lm, smoother = TRUE, qqbands = TRUE, type = "pearson")+
                                  cowplot::draw_figure_label(label = paste0(colnames(df_count)[i]),
                                                             size = 12, fontface = "bold",
                                                             position = "top.right"))
      })
      save_plot(filename = paste0(residual_dir,
                                  colnames(df_count)[i],
                                  "_Residual Diagnostic Plot (LM).pdf"), plot = p)
      
      ### Residual Bootstrapping lm
      
      ## merging fitted and residuals from model into dataset
      temp_comb <- temp_comb %>%
        mutate(
          fit = fitted(lm)[match(term, names(fitted(lm)))],
          resid = resid(lm)[match(term, names(resid(lm)))]
        ) %>% 
        filter(., !is.na(Response))
      
      ## residual resampling
      stat <- function(temp_comb, indices){
        data.star  <- temp_comb %>% mutate(Response = fit + resid[indices])
        model.star <- lm(lm.formula, data = data.star)
        output     <- coef(model.star)
        return(output)
      }
      
      lm.boot <- boot(temp_comb, stat, R=bts)
      lm.est <- summary(lm.boot)
      rownames(lm.est) <- names(lm.boot$t0)
      lm.CI  <- confint(lm.boot, type="perc")
      
      boot.stats <- cbind(lm.est, lm.CI) %>%
        mutate(bootMean = original + bootBias) %>%
        rownames_to_column(., var = "term") %>%
        select(., c("term", "bootMean", "bootSE", "2.5 %", "97.5 %")) %>%
        rename(., term = "term", rep.mean = "bootMean", se = "bootSE" , lower = "2.5 %", upper = "97.5 %") %>%
        mutate(t_statistic = rep.mean / se) %>%
        mutate(pval = 2 * (1 - pt(abs(t_statistic), df = bts - 1)))
      
      # Create dataframe of results for Fixed1
      fix1.df <- dplyr::filter(boot.stats, grepl(paste0("^",fixed1), term)) %>% 
        mutate(term = gsub(paste0(fixed1), "", term),
               fixed = fixed1)
      
      # Comb_res
      comb_res <- fix1.df %>%
        mutate(Response = paste0(colnames(df_count)[i]))
      
      #### Obtaining predicted values
      stat1 <- function(temp_comb, indices){
        data.star  <- temp_comb %>% mutate(Response = fit + resid[indices])
        model.star <- lm(lm.formula, data=data.star)
        output     <- predict(model.star)
        return(output)
      }
      
      predval <- t(try(boot(temp_comb, stat1, R = bts)$t))
      
      mm_pred[colnames(df_count)[i]] <- list(data.frame(predict = rowMeans(predval), se = apply(predval, 1, sd)) %>%
                                               cbind(., temp_comb) %>% 
                                               select(-c(term, Response))
      )
      
    } else {
      
      ### Residual Diagnostic if lmer
      p = suppressMessages({try(ggResidpanel::resid_panel(mm, smoother = TRUE, qqbands = TRUE, type = "pearson")+
                                  cowplot::draw_figure_label(label = paste0(colnames(df_count)[i]),
                                                             size = 12, fontface = "bold",
                                                             position = "top.right"))
      })
      save_plot(filename = paste0(residual_dir,
                                  colnames(df_count)[i],
                                  "_Residual Diagnostic Plot (LMER).pdf"), plot = p)
      
      ### Residual Bootstrapping lmer
      mm.boot <- try(resid_bootstrap(mm, .f = fixef,  B = bts))
      boot.stats <- as.data.frame(mm.boot$stats) %>%
        select(-c("observed","bias")) %>%
        mutate(t_statistic = rep.mean/se) %>%
        mutate(pval = 2 * (1 - pt(abs(t_statistic), df = bts - 1)))
      CI <- as.data.frame( confint(mm.boot, type="perc")) ## using Percentile confidence intervals
      boot.stats <- left_join(boot.stats, CI[,c(1,3,4)], by="term")
      
      fix1.df <- boot.stats[2:(length(levels(temp_comb[[fixed1]]))),] %>%
        mutate(term = levels(df_meta[[fixed1]])[-1]) %>%
        mutate(fixed = fixed1)
      
      comb_res <- fix1.df %>%
        mutate(Response = paste0(colnames(df_count)[i]))
      
      #### Obtaining predicted values
      mm_pred[colnames(df_count)[i]] <- list(as.data.frame(try(resid_bootstrap(mm, .f = predict,  B = bts)$stats)) %>%
                                               right_join(., temp_comb, by = "term") %>%
                                               tidyr::drop_na(rep.mean) %>%
                                               select(-c(term, Response)) %>%
                                               rename(predict = "rep.mean")
      )
      
      rm(mm, boot.stats, CI)
    }
    
    if(i == 1){
      mm_res = comb_res
    }else{
      mm_res = rbind(mm_res,comb_res)
    }
    
  }
  
  
  rm(comb_res, temp_comb)
  
  rownames(mm_res) = NULL
  
  mm_res <- mm_res %>%
    dplyr::group_by(term, fixed) %>%
    mutate(fdr = p.adjust(pval, method = "fdr")) %>% 
    ungroup()
  
  return(list(mm_res = mm_res, mm_pred = mm_pred))
  
}
