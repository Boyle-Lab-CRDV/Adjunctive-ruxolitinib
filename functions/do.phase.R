mm.transform <- function(df_count, log = "log10") {
  zero.cols <- df_count %>% 
    select(where(~ any(. == 0))) %>% 
    colnames()
  normal.cols <- setdiff(colnames(df_count), zero.cols)
  
  if(log == "log10"){
    if(purrr::is_empty(zero.cols) == TRUE){
      df_count[normal.cols] <- lapply(df_count[normal.cols], function(x) log10(x))
    } else {
      df_count[zero.cols] <- lapply(df_count[zero.cols], function(x) log10(x + (min(x[x != 0])/2)))
      df_count[normal.cols] <- lapply(df_count[normal.cols], function(x) log10(x))
    }
  }
  
  if(log == "log2"){
    if(purrr::is_empty(zero.cols) == TRUE){
      df_count[normal.cols] <- lapply(df_count[normal.cols], function(x) log2(x))
    } else {
      df_count[zero.cols] <- lapply(df_count[zero.cols], function(x) log2(x + (min(x[x != 0])/2)))
      df_count[normal.cols] <- lapply(df_count[normal.cols], function(x) log2(x))
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

do.phase <- function(df_count, df_meta, bts = 5000, model = "mm"){
  set.seed(123) ## Pseudonumbering to reproduce consistent random numbers
  i = 1
  
  cat(crayon::green("Analysing data using linear mixed-effect model with residual bootstrapping\n"))
  cat(crayon::red("Will take awhile to run the analysis - be patient\n"))
  
  mm_pred <- list()
  
  suppressWarnings(dir.create(paste0("Residual Diagnostic")))
  residual_dir <- paste0("Residual Diagnostic/")
  
  for(i in 1:ncol(df_count)){
    
    cat(i, "----",colnames(df_count)[i],"\n")
    
    temp_comb <- cbind(df_meta, df_count[,i])
    colnames(temp_comb)[ncol(temp_comb)] <- "Response"
    temp_comb$term <- rownames(temp_comb)
    head(temp_comb)
    
    # Add other metadata
    temp_comb <- temp_comb %>% 
      mutate(TreatmentPhase = ifelse(Treatment == "RUX" & Phase == "REINFECTION", 1, 0),
             Phase2_Day = case_when(Day < 91 ~ 0,
                                    Day >= 91 ~ Day - 91,
                                    TRUE ~ NA))
    
    mm.formula <- as.formula(paste0("Response ~ ",
                                    "Day + Phase + Phase2_Day + Treatment + Treatment*Phase2_Day",
                                    "+ (1|ID)"))
    
    mm = try(lmer(mm.formula, data = temp_comb, REML = TRUE))
    
    
    if(summary(mm)$varcor == 0 | isSingular(mm) | model == "lm"){
      cat(crayon::red("Random effects were ZERO, fitting regular linear model instead\n"))
      
      lm.formula <- as.formula(paste0("Response~",
                                      "Day + Phase + Phase2_Day + Treatment + Treatment*Phase2_Day")) # LM Formula
      lm1 = try(lm(lm.formula, data = temp_comb)) # Regular LM
      
      ### Residual Bootstrapping lm
      temp_comb1 <- temp_comb %>%
        mutate(
          fit = fitted(lm1)[match(term, names(fitted(lm1)))],
          resid = resid(lm1)[match(term, names(resid(lm1)))]
        ) %>% 
        filter(., !is.na(Response))
      
      ## residual resampling
      stat_coef <- function(temp_comb1, indices){
        data.star  <- temp_comb1 %>% mutate(Response = fit + resid[indices])
        model.star <- lm(lm.formula, data = data.star)
        output     <- coef(model.star)
        return(output)
      }
      
      set.seed(123) ## Pseudonumbering to reproduce consistent random numbers
      lm.boot <- boot(temp_comb1, stat_coef, R=bts)
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
      
      if(boot.stats[boot.stats$term == "Phase2_Day:TreatmentRUX",]$pval > 0.05){
        
        lm.formula2 <- as.formula(paste0("Response~",
                                         "Day + Treatment + Phase + Phase2_Day"))
        
        lm2 = try(lm(lm.formula2, data = temp_comb)) # Regular LM
        
        ## merging fitted and residuals from model into dataset
        temp_comb2 <- temp_comb %>%
          mutate(
            fit = fitted(lm2)[match(term, names(fitted(lm2)))],
            resid = resid(lm2)[match(term, names(resid(lm2)))]
          ) %>% 
          filter(., !is.na(Response))
        
        ## residual resampling
        stat_coef2 <- function(temp_comb2, indices){
          data.star  <- temp_comb2 %>% mutate(Response = fit + resid[indices])
          model.star <- lm(lm.formula2, data = data.star)
          output     <- coef(model.star)
          return(output)
        }
        
        p = suppressMessages({try(ggResidpanel::resid_panel(lm2, smoother = TRUE, 
                                                            qqbands = TRUE, type = "pearson")+
                                    cowplot::draw_figure_label(label = paste0(colnames(df_count)[i]),
                                                               size = 12, fontface = "bold",
                                                               position = "top.right"))
        })
        save_plot(filename = paste0(residual_dir,
                                    colnames(df_count)[i],
                                    "_Residual Diagnostic Plot (LM).pdf"), plot = p)
        
        set.seed(123) ## Pseudonumbering to reproduce consistent random numbers
        lm.boot2 <- boot(temp_comb2, stat_coef2, R=bts)
        lm.est2 <- summary(lm.boot2)
        rownames(lm.est2) <- names(lm.boot2$t0)
        lm.CI2  <- confint(lm.boot2, type="perc")
        
        boot.stats2 <- cbind(lm.est2, lm.CI2) %>%
          mutate(bootMean = original + bootBias) %>%
          rownames_to_column(., var = "term") %>%
          select(., c("term", "bootMean", "bootSE", "2.5 %", "97.5 %")) %>%
          rename(., term = "term", rep.mean = "bootMean", se = "bootSE" , lower = "2.5 %", upper = "97.5 %") %>%
          mutate(t_statistic = rep.mean / se) %>%
          mutate(pval = 2 * (1 - pt(abs(t_statistic), df = bts - 1)))
        
        comb_res <- filter(boot.stats2, term == "Phase2_Day") %>%
          mutate(Response = paste0(colnames(df_count)[i]))
        
        #### Obtaining predicted values
        stat_pred2 <- function(temp_comb2, indices){
          data.star  <- temp_comb2 %>% mutate(Response = fit + resid[indices])
          model.star <- lm(lm.formula2, data=data.star)
          output     <- predict(model.star)
          return(output)
        }
        
        predval2 <- t(try(boot(temp_comb2 %>% filter(., !is.na(Response)), stat_pred2, R = bts)$t))
        
        mm_pred[colnames(df_count)[i]] <- list(data.frame(predict = rowMeans(predval2), se = apply(predval2, 1, sd)) %>%
                                                 cbind(., temp_comb2) %>% 
                                                 select(-c(term, Response))
        )
        
        rm(boot.stats2, lm2, lm.boot2, lm.est2, lm.formula2, lm.CI2, temp_comb2, stat_coef2, predval2, stat_pred2)
        
      } else {
        
        
        ### Residual Diagnostic
        p = suppressMessages({try(ggResidpanel::resid_panel(lm1, smoother = TRUE, 
                                                            qqbands = TRUE, type = "pearson")+
                                    cowplot::draw_figure_label(label = paste0(colnames(df_count)[i]),
                                                               size = 12, fontface = "bold",
                                                               position = "top.right"))
        })
        save_plot(filename = paste0(residual_dir,
                                    colnames(df_count)[i],
                                    "_Residual Diagnostic Plot (LM).pdf"), plot = p)
        
        comb_res <- filter(boot.stats, term %in% c("Phase2_Day:TreatmentRUX", "Phase2_Day")) %>%
          mutate(Response = paste0(colnames(df_count)[i]))
        
        #### Obtaining predicted values
        stat_pred1 <- function(temp_comb1, indices){
          data.star  <- temp_comb1 %>% mutate(Response = fit + resid[indices])
          model.star <- lm(lm.formula, data=data.star)
          output     <- predict(model.star)
          return(output)
        }
        
        predval1 <- t(try(boot(temp_comb1 %>% filter(., !is.na(Response)), stat_pred1, R = bts)$t))
        
        mm_pred[colnames(df_count)[i]] <- list(data.frame(predict = rowMeans(predval1), se = apply(predval1, 1, sd)) %>%
                                                 cbind(., temp_comb1) %>% 
                                                 select(-c(term, Response))
        )
        
        rm(boot.stats, lm1, lm.boot, lm.est, lm.formula, lm.CI, temp_comb1, stat_coef, predval1, stat_pred1)
      }
      
    } else {
      
      set.seed(123) ## Pseudonumbering to reproduce consistent random numbers
      mm.boot <- try(resid_bootstrap(mm, .f = fixef,  B = bts))
      boot.stats <- as.data.frame(mm.boot$stats) %>%
        select(-c("observed","bias")) %>%
        mutate(t_statistic = rep.mean/se) %>%
        mutate(pval = 2 * (1 - pt(abs(t_statistic), df = bts - 1))) %>% 
        left_join(., as.data.frame( confint(mm.boot, type="perc")), by = "term") %>% 
        select(., -c("estimate","type","level"))
      
      if(boot.stats[boot.stats$term == "Phase2_Day:TreatmentRUX",]$pval > 0.05){
        
        rm(boot.stats, mm, mm.boot)
        
        mm.formula2 <- as.formula(paste0("Response~",
                                         "Day + Treatment + Phase + Phase2_Day",
                                         "+ (1|ID)"))
        
        mm <- try(lmer(mm.formula2, data = temp_comb, REML = TRUE))
        
        p = suppressMessages({try(ggResidpanel::resid_panel(mm, smoother = TRUE, 
                                                            qqbands = TRUE, type = "pearson")+
                                    cowplot::draw_figure_label(label = paste0(colnames(df_count)[i]),
                                                               size = 12, fontface = "bold",
                                                               position = "top.right"))
        })
        save_plot(filename = paste0(residual_dir,
                                    colnames(df_count)[i],
                                    "_Residual Diagnostic Plot (LMER).pdf"), plot = p)
        
        set.seed(123) ## Pseudonumbering to reproduce consistent random numbers
        mm.boot <- try(resid_bootstrap(mm, .f = fixef,  B = bts))
        boot.stats <- as.data.frame(mm.boot$stats) %>%
          select(-c("observed","bias")) %>%
          mutate(t_statistic = rep.mean/se) %>%
          mutate(pval = 2 * (1 - pt(abs(t_statistic), df = bts - 1))) %>% 
          left_join(., as.data.frame( confint(mm.boot, type="perc")), by = "term") %>% 
          select(., -c("estimate","type","level"))
        
        comb_res <- filter(boot.stats, term == "Phase2_Day") %>%
          mutate(Response = paste0(colnames(df_count)[i]))
        
        rm(mm.boot, boot.stats)
        
      } else{
        p = suppressMessages({try(ggResidpanel::resid_panel(mm, smoother = TRUE, 
                                                            qqbands = TRUE, type = "pearson")+
                                    cowplot::draw_figure_label(label = paste0(colnames(df_count)[i]),
                                                               size = 12, fontface = "bold",
                                                               position = "top.right"))
        })
        save_plot(filename = paste0(residual_dir,
                                    colnames(df_count)[i],
                                    "_Residual Diagnostic Plot (LMER).pdf"), plot = p)
        
        comb_res <- filter(boot.stats, term %in% c("Phase2_Day:TreatmentRUX", "Phase2_Day")) %>%
          mutate(Response = paste0(colnames(df_count)[i]))
      }
      
      ## setting up function for contrasts
      predvals <- function(mod){
        values<-predict(mod, re.form = NA)
        values
      }
      
      mm_pred[colnames(df_count)[i]] <- list(left_join(
        as.data.frame(try(resid_bootstrap(mm, .f = predict,  B = bts)$stats)) %>%
          right_join(., temp_comb, by = "term") %>%
          tidyr::drop_na(rep.mean) %>%
          rename(predict = "rep.mean") %>% 
          select(c(term, predict)),
        as.data.frame(try(resid_bootstrap(mm, .f = predvals,  B = bts)$stats)) %>%
          right_join(., temp_comb, by = "term") %>%
          tidyr::drop_na(observed) %>%
          rename(predval = "rep.mean"),
        by = "term"
      ) %>% 
        select(., -c("term")))
      
    }
    
    
    if(i == 1){
      mm_res = comb_res
    }else{
      mm_res = rbind(mm_res,comb_res)
    }  
    
  }
  
  
  rm(mm, comb_res, temp_comb)
  
  rownames(mm_res) = NULL
  
  mm_res <- mm_res %>%
    dplyr::group_by(term) %>%
    mutate(fdr = p.adjust(pval, method = "fdr")) %>% 
    ungroup()
  
  return(list(mm_res = mm_res, mm_pred = mm_pred))
}
