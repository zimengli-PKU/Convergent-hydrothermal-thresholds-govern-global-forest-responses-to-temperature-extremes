# Comparison: all-data OLS, spatial-resampled OLS,spatial-resamlped LMM 

library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(cowplot)

# 1.parameter settings
lon_col<-"lon"
lat_col<-"lat"
species_col<-"Species"

block_size<-5
records_per_block<-20
bin_width<-100
B<-1000

min_n<-21
min_blocks<-3
alpha<-0.05
significant_cutoff<-0.90

# 2. dataset import 
raw_dat<-read.csv("response.csv")
load("TP60.RData")

to_numeric<-function(x) {as.numeric(as.character(x))}
raw_dat$MAT<-to_numeric(T60$MAT)
raw_dat$MAP<-to_numeric(P60$MAP)
raw_dat$longitude<-to_numeric(raw_dat[[lon_col]])
raw_dat$latitude<-to_numeric(raw_dat[[lat_col]])
raw_dat$species<-trimws(as.character(raw_dat[[species_col]]))

finite_mean<-function(x) {x<-x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_}

finite_q<-function(x, p) {
  x<-x[is.finite(x)]
  if (length(x) >= 2) {unname(quantile(x, p))
  } else {NA_real_} }

run_analysis<-function(response_var, seed=0) {
  set.seed(seed)
  dat<-raw_dat
  dat$response<-to_numeric(dat[[response_var]])
  dat<-dat %>%
    filter(is.finite(response),
      is.finite(MAT),is.finite(MAP),MAP>0,is.finite(longitude),
      is.finite(latitude),!is.na(species),species != "") %>%
    mutate(longitude=((longitude+180) %% 360)-180,
      block=paste(floor((longitude+180)/block_size),
        pmin(floor((latitude+90)/block_size),
          ceiling(180/block_size)-1),sep="_"),
      species=factor(species),
      MAP_bin=ceiling(MAP/bin_width))
  bin_ids<-sort(unique(dat$MAP_bin))
  sample_table<-dat %>%
    group_by(MAP_bin) %>%
    summarise(MAP_mean=mean(MAP),
      N=n(),N_blocks=n_distinct(block),
      MAT_min=min(MAT),MAT_max=max(MAT),
      .groups="drop")
  
  # 3.1 MAP bin
  fit_one<-function(d, model, iteration, bin) {
    
    out<-data.frame(
      MAP_bin=bin,Model=model,Iteration=iteration,
      Intercept=NA_real_, Slope=NA_real_, Slope_SE=NA_real_,
      Slope_p=NA_real_, Slope_CI_low=NA_real_,  Slope_CI_high=NA_real_,
      Threshold=NA_real_,  Threshold_inverse_SE=NA_real_,
      Threshold_delta_low=NA_real_, Threshold_delta_high=NA_real_,
      Singular=FALSE, Fit_OK=FALSE)
    
    if (nrow(d)<min_n||n_distinct(d$block)< min_blocks||n_distinct(d$MAT)<2) {
      return(out)}
    
    mixed<-model=="Subsample_species_LMM"
    d$species<-droplevels(d$species)
    
    if (mixed &&(nlevels(d$species)<2||nlevels(d$species)>=nrow(d)))
      {return(out)}
    ref<-mean(d$MAT)
    d$MAT_c<-d$MAT-ref

    fit<-tryCatch({if (mixed) {lmerTest::lmer(
            response ~ MAT_c+(1 | species),data=d, 
            REML=TRUE,control=lmerControl(
              optimizer="bobyqa", optCtrl=list(maxfun=100000),
              check.conv.singular="ignore"))
        } else {
          lm(response ~ MAT_c, data=d)
        }
      },error=function(e) NULL
    )
    
    cf<-coef(summary(fit))
    beta<-if (mixed) fixef(fit) else coef(fit)
    
    a<-unname(beta["(Intercept)"])
    b<-unname(beta["MAT_c"])
    
    if (!is.finite(a) || !is.finite(b)) {
      return(out)
    }
    
    slope_se<-cf["MAT_c", "Std. Error"]
    slope_p<-cf["MAT_c", "Pr(>|t|)"]
    
    df<-if (mixed) {
      cf["MAT_c", "df"]
    } else {
      df.residual(fit)
    }
    
    crit<-if (is.finite(df) && df>0) {
      qt(1-alpha/2, df)
    } else {NA_real_}
    
    out$Intercept<-a-ref * b
    out$Slope<-b
    out$Slope_SE<-slope_se
    out$Slope_p<-slope_p
    
    if (model == "All_OLS") {
        slope_ci<-tryCatch( confint(fit, parm="MAT_c", level=1-alpha),
        error=function(e) NULL)
      
      if (!is.null(slope_ci) && all(is.finite(slope_ci))) {
        out$Slope_CI_low<-unname(slope_ci[1])
        out$Slope_CI_high<-unname(slope_ci[2])
      }
    } else if (is.finite(crit) && is.finite(slope_se)) {
      out$Slope_CI_low<-b-crit * slope_se
      out$Slope_CI_high<-b+crit * slope_se
    }
    
    out$Fit_OK<-TRUE
    
    if (b != 0) {
      threshold<-ref-a/b
      
      if (is.finite(threshold)) {
        out$Threshold<-threshold
        
        if (model == "All_OLS") {
          inverse_fit<-tryCatch(
            lm(MAT ~ response, data=d),
            error=function(e) NULL
          )
          
          inverse_prediction<-if (!is.null(inverse_fit)) {
            tryCatch(
              predict(
                inverse_fit,
                newdata=data.frame(response=0),
                se.fit=TRUE
              ),
              error=function(e) NULL
            )
          } else {NULL}
          
          if (!is.null(inverse_prediction)) {
            inverse_se<-unname(inverse_prediction$se.fit[1])
            
            if (is.finite(inverse_se)) {
              out$Threshold_inverse_SE<-inverse_se
              out$Threshold_delta_low<-threshold-1.96 * inverse_se
              out$Threshold_delta_high<-threshold+1.96 * inverse_se
            }
          }
        } else {
           V<-as.matrix(vcov(fit))[
            c("(Intercept)", "MAT_c"),
            c("(Intercept)", "MAT_c")
          ]
          
          gradient<-c(-1/b, a/b^2)
          threshold_var<-as.numeric(
            t(gradient) %*% V %*% gradient
          )
          
          if (
            is.finite(threshold_var) &&
            threshold_var >= 0 &&
            is.finite(crit)
          ) {
            threshold_se<-sqrt(threshold_var)
            out$Threshold_delta_low<-threshold-crit * threshold_se
            out$Threshold_delta_high<-threshold+crit * threshold_se
          }
        }
      }
    }
    
    out
  }
  
 #all-data OLS
  baseline_results<-bind_rows(
    lapply(bin_ids, function(i) {
      fit_one(
        dat[dat$MAP_bin == i, ],
        model="All_OLS",
        iteration=0,
        bin=i
      )
    })
  )
  
  # spatial balance subsampling
  block_indices<-split(seq_len(nrow(dat)), dat$block)
  subsample_parts<-vector("list", B)
  
  for (r in seq_len(B)) {
    ids<-unlist(
      lapply(block_indices, function(ix) {
        sample(
          ix,
          size=min(length(ix), records_per_block),
          replace=FALSE
        )
      }),
      use.names=FALSE
    )
    
    dsub<-dat[ids, ]
    
    subsample_parts[[r]]<-bind_rows(
      lapply(bin_ids, function(i) {
        d<-dsub[dsub$MAP_bin == i, ]
        
        bind_rows(
          fit_one(
            d,
            model="Subsample_OLS",
            iteration=r,
            bin=i
          ),
          fit_one(
            d,
            model="Subsample_species_LMM",
            iteration=r,
            bin=i
          )
        )
      })
    )
    
    if (r %% 50 == 0) {
      message(response_var, " subsampling: ", r, "/", B)
    }
  }
  
  subsample_results<-bind_rows(subsample_parts)
  
  ensemble_summary<-subsample_results %>%
    group_by(MAP_bin, Model) %>%
    group_modify(~ {
      d<-.x
      
      mean_intercept<-finite_mean(d$Intercept)
      mean_slope<-finite_mean(d$Slope)
      
      significant_same_direction <-
        d$Fit_OK &
        is.finite(d$Slope_p) &
        is.finite(d$Slope) &
        d$Slope_p < alpha &
        is.finite(mean_slope) &
        sign(d$Slope) == sign(mean_slope)
      significant_fraction <-
        sum(significant_same_direction, na.rm=TRUE)/B
      
      threshold<-if (
        is.finite(mean_intercept) &&
        is.finite(mean_slope) &&
        mean_slope != 0
      ) {
        -mean_intercept/mean_slope
      } else {
        NA_real_
      }
      
      tibble(
        Intercept=mean_intercept,
        Slope=mean_slope,
        Slope_interval_low=finite_q(d$Slope, alpha/2),
        Slope_interval_high=finite_q(d$Slope, 1-alpha/2),
        Slope_p=NA_real_,
        N_successful=sum(d$Fit_OK),
        Success_fraction=mean(d$Fit_OK),
        Significant_same_direction_fraction=significant_fraction,
        Singular_fraction=finite_mean(as.numeric(d$Singular[d$Fit_OK])),
        Supported=significant_fraction >= significant_cutoff,
        Threshold=threshold,
        Threshold_interval_low=finite_q(d$Threshold, alpha/2),
        Threshold_interval_high=finite_q(d$Threshold, 1-alpha/2)
      )
    }) %>%
    ungroup()
  
  baseline_summary<-baseline_results %>%
    transmute(
      MAP_bin,
      Model,
      Intercept,
      Slope,
      Slope_interval_low=Slope_CI_low,
      Slope_interval_high=Slope_CI_high,
      Slope_p,
      N_successful=as.integer(Fit_OK),
      Success_fraction=as.numeric(Fit_OK),
      Significant_same_direction_fraction=NA_real_,
      Singular_fraction=NA_real_,
      Supported =
        Fit_OK &
        is.finite(Slope_p) &
        Slope_p < alpha,
      Threshold,
      Threshold_interval_low=if_else(
        Supported,
        Threshold_delta_low,
        NA_real_
      ),
      Threshold_interval_high=if_else(
        Supported,
        Threshold_delta_high,
        NA_real_
      )
    )
  
  #Summary -------------------------------------------------------
  model_order<-c("All_OLS",
    "Subsample_OLS",
    "Subsample_species_LMM")
  
  comparison_table<-bind_rows(baseline_summary,ensemble_summary) %>%
    left_join(sample_table, by="MAP_bin") %>%
    mutate(Model=factor(Model, levels=model_order),
      Evidence=if_else(Supported,"Significant","Non-significant"),
      In_MAT_range =is.finite(Threshold)&Threshold>= MAT_min&Threshold <= MAT_max)
  %>%arrange(MAP_bin, Model)
  model_colours<-c(All_OLS="#DE776E",Subsample_OLS="#65B948",Subsample_species_LMM="#7298F7")
  
  model_labels<-c("All-data SLR","Subsampled SLR","Subsampled species LMM")
  
  offsets<-c(All_OLS=-0.1,Subsample_OLS=0,Subsample_species_LMM=0.1) * bin_width
  
  plot_data<-comparison_table %>%
    mutate(MAP_display =MAP_mean +unname(offsets[as.character(Model)]))
  
  make_plot<-function(d, value, lower, upper, y_label) {
    d<-d[is.finite(d[[value]]), ]
    bars<-d[is.finite(d[[lower]]) & is.finite(d[[upper]]),]
    
    ggplot(d,aes(x=MAP_display,y=.data[[value]],colour=Model)) +
      geom_errorbar(data=bars,aes(ymin=.data[[lower]],ymax=.data[[upper]]),
        width=12, alpha=0.5,linewidth=0.55) +
      geom_point(aes(shape=Evidence), size=2.8, stroke=0.9) +
      scale_colour_manual(values=model_colours,
        breaks=model_order,labels=model_labels, name=NULL,drop=FALSE) +
      scale_shape_manual(values=c("Non-significant"=1,"Significant"=16),
        breaks=c("Non-significant", "Significant"),
        name="Slope significance") +
      labs(x="MAP (mm)", y=y_label) +
      theme_classic(base_size=12) +
      theme(axis.text=element_text(colour="black"))}
  
  if (response_var == "hotr") {
    slope_label<-expression("Slope of " * Delta * "RWI"[h] * " ~ MAT")
    threshold_label<-expression("Threshold MAT"[h] * " (" * degree * "C)")
  } else {slope_label<-expression("Slope of " * Delta * "RWI"[c] * " ~ MAT")
    threshold_label<-expression("Threshold MAT"[c] * " (" * degree * "C)")
  }
  
  p_slope<-make_plot(plot_data,value="Slope",
    lower="Slope_interval_low",upper="Slope_interval_high",
    y_label=slope_label) +
    geom_hline(yintercept=0,linetype="dashed", colour="grey50")
  
  p_threshold<-make_plot(filter(plot_data, Supported),value="Threshold",
    lower="Threshold_interval_low",   upper="Threshold_interval_high",
    y_label=threshold_label ) +
    theme(legend.position="none")
  
  list(
    p_slope=p_slope,
    p_threshold=p_threshold,
    plot_data=plot_data,
    comparison_table=comparison_table,
    baseline_results=baseline_results,
    subsample_results=subsample_results
  )
}


hot_result<-run_analysis("hotr", seed=0)
cold_result<-run_analysis("coldr", seed=0)
