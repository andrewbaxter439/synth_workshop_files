# ---
# title: "Synth function testing TPS"
# author: "Andy Baxter"
# date: "17/09/2021"
# output: html_document
# ---
# 
# # Using the Synth package to do synthetic control analysis on England's
# Teenage Pregnancy Strategy
#
# This is some R code to repeat the analyses described in the video and
# contained within the [shiny app](https://phd.andybaxter.me/synth-app). These
# steps are described in Abadie et al.,
# [2011](http://www.jstatsoft.org/v42/i13/).
#
# Setting up and importing data:
# 
# # ----setup ----------------------------------------------------------------------------------------------

load("data/synth_dat.rdata")
library(tidyverse)
library(Synth)

theme_set(theme_minimal())

#
# # # Graphing pre-intervention rates across countries
#
# The intervention year is set to 1999. Looking at the pre-intervention years,
# England and Wales are not an 'outlier', so the Synthetic Control method should
# be possible here.
#
# #
# ----graphing_pre------------------------------------------------------------------------------------------------------

u18_birth_rates %>% 
  filter(Year < 1999) %>% 
  mutate(Country = fct_reorder2(Country, Year, rate)) %>% 
  ggplot(aes(Year, rate, colour = Country)) +
  geom_line() 


#
# # # Preparing data
#
# This step is described by Abadie et al., (2011) and uses the `dataprep`
# function from the Synth package.
#
# #
# ----dataprep----------------------------------------------------------------------------------------------------------

dp_u18 <- dataprep(
  foo = u18_birth_rates,
  dependent = "rate",
  # `synth` can take the mean for each country across all years as one
  # predictor. Here I've opted to split years into clustered periods, as this
  # improves model fit by allowing fitted trend to match England's
  # pre-intervention trend shape
  special.predictors = list(
    a = list("rate", yrs = 1990:1993, "mean"),
    b = list("rate", yrs = 1994L, "mean"),
    c = list("rate", yrs = 1995L, "mean"),
    d = list("rate", yrs = 1996:1998, "mean")
  ),
  time.variable = "Year",
  unit.variable = "Code", 
  unit.names.variable = "Country",
  time.optimize.ssr = 1990:1998,
  time.predictors.prior = 1990:1998,
  treatment.identifier = ccodes$Code[ccodes$Country == "England and Wales"],
  controls.identifier = ccodes$Code[ccodes$Country != "England and Wales"],
  time.plot = 1990:2013
)


#
# Aside - the `dataprep` function creates four tables to pass on to the `synth`
# function. These contain the yearly birth rates across the optimising period
# (1990 to 1998) for England, the yearly birth rates for control countries, and
# the means of birth rates across the four groups of predictor years (see above)
# fo England and for control countries. Uncomment the lines below to see the
# data that `synth` will use when constructing a control series:
#
# #
# ----compute_tables----------------------------------------------------------------------------------------------------

# knitr::kable(dp_u18$Z1)
# knitr::kable(dp_u18$Z0)
# knitr::kable(dp_u18$X1)
# knitr::kable(dp_u18$X0)


#
# It also stores the rates for each country (England and Wales, control
# countries) for the whole of the plotting period, which we will use later
# (`dp_u18$Y1` and `dp_u18$Y0`).
#
#
# # # Fitting model
#
# #
# ----------------------------------------------------------------------------------------------------------------------

# The `synth` function fits the prepared date
# The `so` object we create here can be used to further explore the result
so <- synth(dp_u18)


#
# The `synth` function optimises pre-intervention fit (minimising MSPE) by
# assigning a series of weights to donor countries. The outputs above:
#
# - The weights applied to each country are listed in `soution.w`. - The
# pre-intervention MSPE for our model (MSPE/LOSS V) is `r so$loss.v`.
#
# We can use these weights to get predictions:
#
# #
# ----getting_pred------------------------------------------------------------------------------------------------------

# This multiples the matrix of donor values by year by the vector of country
# weights to produce a time series vector of predicted rates

# dp_u18$Y0 = observed rates across control countries
synthC <- dp_u18$Y0 %*% so$solution.w

# Extracting the mspe for use later
mspe <- so$loss.v[1]


md <-
  tibble(
    Year = 1990:2013,
    Treated = dp_u18$Y1[, 1],  # dp_u18$Y1 = observed rates in England
    Synthetic = synthC[, 1]    # Our predicted rates from above
  ) %>%
  pivot_longer(-Year, names_to = "Group", values_to = "Rate")

md

# 
# # Graphing predictions
# 
# Pre-intervention period only:
# 
# # ----pre-only----------------------------------------------------------------------------------------------------------
md %>% 
  filter(Year < 1999) %>% 
  ggplot(aes(Year, Rate, colour = Group)) +
  geom_line() +
  ylim(0, NA) +
  geom_vline(xintercept = 1998.5, linetype = "dashed")

#
# This seems like a good fit! Let's plot the whole period to see if the strategy
# produced a noticable difference after 1999:
#
# #
# ----whole_graph-------------------------------------------------------------------------------------------------------


md %>%
  ggplot(aes(Year, Rate, colour = Group)) +
  geom_line() +
  ylim(0, NA) +
  geom_vline(xintercept = 1998.5, linetype = "dashed") +
  labs(subtitle = paste0("Pre-intervention MSPE = ", signif(mspe, 3)))


# 
# # Looking at weightngs
# 
# `Synth` can output some tables for us to see what's going on under the surface:
# 
# # ----------------------------------------------------------------------------------------------------------------------


st <- synth.tab(so, dp_u18)

st$tab.w %>% 
  arrange(desc(w.weights)) %>% 
  select(Country = unit.names, Weight = w.weights) %>%
  knitr::kable(caption = "Country weights")


# Scotland has been given the majority of the weight (67.2%). Portugal (29.5%),
# the USA (1.6%) and New Zealand (1.2%) have been included and weighted too.
#
# Note that these weights add up to 1. This is one of the advantages of the
# Synthetic Control method - it doesn't extrapolate outside of the limits set by
# the highest and lowest observations in the donor pool in choosing suitable
# comparators. If England and Wales had shown the highest rates across the
# pre-intervention period, the Synthetic Control method assumes that multiplying
# lower rates by a factor greater than 1 (e.g. doubling the observed rates in a
# low-birthrate country) does not give a valid prediction of the relative rate
# changes in England and Wales, and therefore there would be no suitable
# comparators.
#
# # # Placebo tests by country
#
# The next task is to do country placebo tests. To do this we remove England and
# Wales from the donor pool and re-run for all other countries. This takes time
# to compute. Do set it to run we can write a quick function with all the steps
# we've taken so far:
#
# #
# ----generatePlacebos_func---------------------------------------------------------------------------------------------

generatePlacebos <- function(data,
                             predictors = NULL,
                             special.predictors = NULL, 
                             time.optimize.ssr = 1990:1998,
                             time.plot = 1990:2013, 
                             dependent = "rate", 
                             ...) {
  require(dplyr)
  require(Synth)
  require(tidyr)
  
  data <- data.frame(data %>% 
                       filter(Country != "England and Wales"))
  
  ccodes <- data %>% 
    select(Code, Country) %>% 
    arrange(Code) %>% 
    unique()
  
  placebos <- data.frame()
  
  for(i in 1:nrow(ccodes)){
    
    dp <- dataprep(
      data,
      predictors = predictors,
      predictors.op = "mean",
      special.predictors = special.predictors,
      time.predictors.prior = time.plot[1]:1998,
      dependent = dependent,
      unit.variable = "Code",
      unit.names.variable = "Country",
      time.variable = "Year",
      treatment.identifier = ccodes$Code[i],
      controls.identifier = ccodes$Code[-i],
      time.optimize.ssr = time.optimize.ssr,
      time.plot = time.plot
    )
    

    so <- synth(dp, ...)
    synthC <- dp$Y0 %*% so$solution.w
    
    md <-
      tibble(
        Year = 1990:2013,
        Treated = dp$Y1[, 1],
        Synthetic = synthC[, 1]
      ) %>%
      pivot_longer(-Year, names_to = "Group", values_to = "Rate")
    
    mspe <- so$loss.v
    
    # We'll get this function to calculate gaps for us as well
    gaps <- md %>%
      spread(Group, Rate) %>% 
      mutate(Country = ccodes$Country[i],
             Gap = Treated - Synthetic,
             pre_mspe = mspe)
    
    placebos <- bind_rows(placebos, gaps)
  }
  return(placebos)
}

#
# And we'll use this to get our placebo predictions alongside observations for
# each country:
#
# #
# ----generating_placebos-----------------------------------------------------------------------------------------------
pl_u18 <- generatePlacebos(u18_birth_rates,
                           special.predictors = list(
                             a = list("rate", yrs = 1990:1993, "mean"),
                             b = list("rate", yrs = 1994L, "mean"),
                             c = list("rate", yrs = 1995L, "mean"),
                             d = list("rate", yrs = 1996:1998, "mean")
                           ))

#
# Now we have placebo tests for each country, our next step with this is to plot
# the 'gaps' - the yearly differences between the observed and predicted rate
# for each country. For speed, let's define another function:
#
# #
# ----gg_gaps_func------------------------------------------------------------------------------------------------------

gg_gaps <- function(md, pl, dp = NULL, mspe_limit) {

  placebos_filtered <- pl %>% 
    filter(pre_mspe < 5*mspe_limit)
  
  p <- md %>% 
    spread(Group, Rate) %>% 
    mutate(Gap = Treated - Synthetic) %>%  # Gaps for England and Wales
    ggplot(aes(Year, Gap)) +
    geom_segment(x = min(md$Year), xend = 2013, y = 0, yend = 0) +
    geom_line(data = placebos_filtered, aes(group = Country), col = "grey") +
    geom_line(col = "#951272", size = 2) +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    geom_vline(xintercept = 1998.5, linetype = "dotted") +
    ylab("Gap = Treated - Synthetic Control")
  
  # Print which countries were filtered out
  placebos_filtered %>% 
    pull(Country) %>% 
    unique() %>% 
    paste(collapse = ", ") %>% 
    cat("Countries removed:", .)
  
  # Return graph
  p
  
}

#
# Now pass our model and placebos to the function, and use the `mspe` from the
# fitted model to filter which countries will show (with less than 5 times the
# MSPE of England and Wales):
#
# #
# ----gg_gaps-----------------------------------------------------------------------------------------------------------

gg_gaps(md, pl_u18, mspe_limit = mspe)


#
# Finally, we can calculate MSPE ratios for each country and compare (using
# another rather long function!):
#
# #
# ----gg_pre_post_func--------------------------------------------------------------------------------------------------

gg_pre_postMSPE_tab <- function(md, pl) {
  
  pts <- function (size) return(size * 5/14)
  
  df <- md %>% 
    spread(Group, Rate) %>% 
    mutate(Gap = Treated - Synthetic,
           Country = "England and Wales") %>% 
    select(Year, Country, Gap) %>% 
    bind_rows(pl %>% select(Year, Country, Gap)) %>% 
    mutate(period = ifelse(Year<1999, "pre", "post")) %>% 
    group_by(Country, period) %>% 
    summarise(mspe = mean(Gap**2)) %>% 
    spread(period, mspe) %>% 
    mutate(ratio = post/pre)
  
  df_tidy <- df %>% 
    ungroup() %>% 
    mutate(rank = rank(-ratio),
           Country = fct_reorder(factor(Country), ratio),
           ratio_lab = str_pad(str_replace(as.character(round(ratio, 2)), "(\\.\\d)$", "\\10"), 6, side = "left"),
           EngWa = ifelse(Country == "England and Wales", "E", "C")) %>% 
    arrange(rank) %>% 
    select(Country, ratio, ratio_lab, rank, EngWa)
  
  
  points_plot <- ggplot(df_tidy, aes(ratio, Country, colour = EngWa)) +
    geom_point(size = 3, shape = 'diamond') +
    scale_x_continuous("Post/pre-MSPE ratio") +
    ggtitle(" ") +
    theme(axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = "#e0e0e0"),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "#666666"),
          strip.text = element_blank(),
          panel.spacing.y = unit(1, "cm"),
          plot.margin = margin(l = 0, r = 0),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(size = 12),
          legend.position = "none"
    ) +
    scale_colour_manual(values = c("E" = "#003865", "C" = "darkgrey"))
  
  base_plot <- df_tidy %>% 
    ggplot(aes(x = 0, y = Country, colour = EngWa, fontface = ifelse(EngWa == "E", "bold", "plain"))) +
    ylab(NULL) +
    xlab(" ") +
    theme(
      strip.text = element_blank(),
      plot.margin = margin(l = 0, r = 0),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(color = "white"),
      panel.background = element_blank(),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 12, hjust = 0)
    ) +
    coord_cartesian(clip = "off") +
    scale_x_continuous(expand = expansion(add = c(0,0)), limits = c(0,2)) +
    scale_colour_manual(values = c("E" = "#003865", "C" = "black"))
  
  c_plot <- base_plot + 
    geom_text(aes(label = Country), hjust = 0, size = pts(9)) +
    ggtitle("Country") +
    theme(strip.text = element_text(colour = "black"),
          strip.background = element_blank(),
          strip.placement = "outside" )
  
  
  rank_plot <- base_plot + 
    geom_text(aes(label = rank), hjust = 0, size = pts(9)) +
    ggtitle("Rank") +
    theme(strip.text = element_text(colour = "black"),
          strip.background = element_blank(),
          strip.placement = "outside" )
  
  rat_plot <- base_plot + 
    geom_text(aes(x = Inf, label = ratio_lab), hjust = 1, size = pts(9)) +
    ggtitle("Ratio") +
    theme(strip.text = element_text(colour = "black"),
          strip.background = element_blank(),
          strip.placement = "outside",
          plot.title = element_text(size = 12, hjust = 1))
  
  # To try and lay out everything in a tidy manner
  gridExtra::grid.arrange(c_plot, rank_plot, rat_plot, points_plot,
                                 layout_matrix = matrix(c(1,1,1,1,1,1,2, 2,3,3, 4, 4,4,4,4,4,4,4,4,4,4,4), nrow = 1))
  
}


# 
# Plotting:
# 
# # ----gg_pre_post_mspe--------------------------------------------------------------------------------------------------

gg_pre_postMSPE_tab(md, pl_u18)


#
# England and Wales saw a pretty low relative ratio, close to 1, indicating a
# relatively well-fit series of pre-intervention predictions and/or a relatively
# small deviation from null-effect predictions across the exposure period.
