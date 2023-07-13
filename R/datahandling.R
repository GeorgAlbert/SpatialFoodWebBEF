library(tidyverse)
library(gridExtra)

# 1. data and functions ####
# _ 1.1. functions ####

# collapse realwebIDs for cluster bash script
bash_short.nmbr <- function(val_sel){
  
  val_short = c()
  for(i in 1:length(val_sel)){
    
    if(i == 1){
      
      val_short = c(val_short, as.character(val_sel[i]))
      
    } else if(val_sel[i] - val_sel[i-1] != 1){
      
      if(i == 2){
        val_short = c(val_short, ",", as.character(val_sel[i]))
      } else {
        if((val_sel[i-1] - val_sel[i-2] != 1)){
          val_short = c(val_short, ",", as.character(val_sel[i]))
        } else {
          val_short = c(val_short, "-", as.character(val_sel[i-1]), ",", as.character(val_sel[i]))
        }
      }
      
    } else if(i == length(val_sel)){
      
      if(val_sel[i] - val_sel[i-1] != 1){
        val_short = c(val_short, ",", as.character(val_sel[i]))
      } else {
        val_short = c(val_short, "-", as.character(val_sel[i]))
      }
      
    }
  }
  
  return(val_short)
}

# summary function
data_summary <- function(data, varnames, groupnames){
  # modified from http://www.sthda.com/english/wiki/ggplot2-error-bars-quick-start-guide-r-software-and-data-visualization
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      
      median = median(x[[col]], na.rm=TRUE),
      quantile25 = quantile(x[[col]], probs = 0.25, na.rm=TRUE),
      quantile75 = quantile(x[[col]], probs = 0.75, na.rm=TRUE),
      
      n = length(x[[col]][!is.na(x[[col]])]))
  }
  
  data_sum <- c()
  for(i in 1:length(varnames)){
    data_sum <- rbind(data_sum, 
                      ddply(data, groupnames, .fun=summary_func, varnames[i]))
  }
  
  data_sum <- cbind(varname = rep(varnames, each = nrow(data_sum) / length(varnames)), data_sum)
  data_sum <- rename(data_sum, c("quantile25.25%" = "quantile25",
                                 "quantile75.75%" = "quantile75"))
  return(data_sum)
}


# _ 1.2. data ####

path.to.data <- "dat/"
dat <- read.csv(paste0(path.to.data, "dat.csv"))


# _ 1.3. add columns ####

# - food web treatment
dat$webtreatment = NA
dat$webtreatment[dat$nA == 0] = "none"
dat$webtreatment[dat$nA > 0 & dat$nested == FALSE] = "non-nested"
dat$webtreatment[dat$nA > 0 & dat$nested == TRUE] = "nested"
dat$webtreatment = factor(dat$webtreatment, levels = c("none", "non-nested", "nested"))

# - total biomass density
dat$dens = dat$densA + dat$densP_sum


# _ 1.4. check completeness - SLURM ####

# identify which simulations are done and which need to be simulated for longer !
tmax = 50000 # --> default 50000
done = dat[dat$tmax == tmax | dat$survP_ind == 0, ]
continue = dat[dat$tmax < tmax & dat$survP_ind != 0, ]


# get list of unfinished realwebIDs that can be copied from console and pasted to SLURM script specifying array job IDs
all = 1:5100 # --> all realwebIDs that we want to simulate - 1:5100 by default ! 
all[!(all %in% done$realwebID)] |>
  bash_short.nmbr() |>
  paste0(collapse = "")



# _ 1.5. check unviable food webs ####

# - whenever all plants are dead, food web is considered unviable

nrow(dat[dat$survP_sp == 0 & dat$webtreatment == "none" & dat$nP == 1, ])
nrow(dat[dat$survP_sp == 0 & dat$webtreatment == "none" & dat$nP == 16, ])

nrow(dat[dat$survP_sp == 0 & dat$webtreatment == "non-nested" & dat$nP == 1, ])
nrow(dat[dat$survP_sp == 0 & dat$webtreatment == "non-nested" & dat$nP == 16, ])

nrow(dat[dat$survP_sp == 0 & dat$webtreatment == "nested" & dat$nP == 1, ])
nrow(dat[dat$survP_sp == 0 & dat$webtreatment == "nested" & dat$nP == 16, ])




# 2. Figures ####
# _ 2.1. Fig.2 - plant diversity-productivity relationships ####

FIG2 <- data_summary(dat[dat$prodP_sum != 0, ], 
             varname = "prodP_sum", 
             groupnames = c("nP", "N_ikk", "webtreatment")) %>%
ggplot(aes(nP, median, col = webtreatment)) +
  geom_line(position = position_dodge(width = 4)) +
  geom_errorbar(aes(ymin=quantile25, ymax=quantile75), 
                position = position_dodge(width = 4), width=4) +
  geom_point(position = position_dodge(width = 4)) +
  facet_grid(.~factor(N_ikk, levels = seq(1, 0.2, by = -0.2))) +
  scale_color_manual(values = c("forestgreen", "blue4", "lightskyblue")) +
  scale_x_continuous(name = "Plant species richness", 
                     limits = c(-3, 20), breaks = c(1, 16)) +
  scale_y_continuous(name = "Plant productivity", breaks = c(700, 1200)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")


ggsave(paste0(path.to.data, "/FIG2_raw.png"), 
       plot = FIG2, width = 173, height = 70, unit = "mm", dpi = 600)



# _ 2.2. Fig.3 - plant diversity and survival ####

# surviving plant species
p1 <- data_summary(dat[dat$nP == 16, ],
             varname = "survP_sp", 
             groupnames = c("N_ikk", "webtreatment")) %>%
  ggplot(aes(N_ikk, median, col = webtreatment)) +
  geom_line(position = position_dodge(width = 0.1), lty = 3) +
  geom_errorbar(aes(ymin=quantile25, ymax=quantile75), 
                position = position_dodge(width = 0.1), width=.1) +
  geom_point(position = position_dodge(width = 0.1)) +
  scale_color_manual(values = c("forestgreen", "blue4", "lightskyblue")) +
  scale_x_continuous(name = "", 
                     breaks = c(.2, .4, .6, .8, 1), labels = c(1, "", "", "", 0),
                     trans = "reverse") +
  scale_y_continuous(name = "Realized plant species richness", 
                     breaks = c(7, 10, 13, 16), limits = c(6.5, 16)) +
  theme_bw() +
  theme(legend.position = "none")


# surviving plant individuals
p2 <- data_summary(dat[dat$nP == 16, ],
                   varname = "survP_ind", 
                   groupnames = c("N_ikk", "webtreatment")) %>%
  ggplot(aes(N_ikk, median, col = webtreatment)) +
  geom_line(position = position_dodge(width = 0.1), lty = 3) +
  geom_errorbar(aes(ymin=quantile25, ymax=quantile75), 
                position = position_dodge(width = 0.1), width=.1) +
  geom_point(position = position_dodge(width = 0.1)) +
  scale_color_manual(values = c("forestgreen", "blue4", "lightskyblue")) +
  scale_x_continuous(name = "Spatial resource overlap", 
                     breaks = c(.2, .4, .6, .8, 1), labels = c(1, "", "", "", 0),
                     trans = "reverse") +
  scale_y_continuous(name = "Realized plant density", breaks = c(32, 48, 64)) +
  theme_bw() +
  theme(legend.position = "none")


# shannon diversity
p3 <- data_summary(dat[dat$nP == 16, ],
             varname = "shannon_div_prod", 
             groupnames = c("N_ikk", "webtreatment")) %>%
  ggplot(aes(N_ikk, median, col = webtreatment)) +
  geom_line(position = position_dodge(width = 0.1), lty = 3) +
  geom_errorbar(aes(ymin=quantile25, ymax=quantile75), 
                position = position_dodge(width = 0.1), width=.1) +
  geom_point(position = position_dodge(width = 0.1)) +
  scale_color_manual(values = c("forestgreen", "blue4", "lightskyblue")) +
  scale_x_continuous(name = "", 
                     breaks = c(.2, .4, .6, .8, 1), labels = c(1, "", "", "", 0),
                     trans = "reverse") +
  scale_y_continuous(name = "Shannon diversity") +
  theme_bw() +
  theme(legend.position = "none")



FIG3 <- arrangeGrob(p1, p2, p3, ncol = 3)


ggsave(paste0(path.to.data, "/FIG3_raw.png"), 
       plot = FIG3, width = 173, height = 80, unit = "mm", dpi = 600)



# _ 2.3. Fig.S1 productivity split for monocultures and mixtures ####

FIGS1 <- data_summary(dat[dat$prodP_sum != 0, ], 
                      varname = "prodP_sum", 
                      groupnames = c("N_ikk", "webtreatment", "nP")) %>%
  ggplot(aes(N_ikk, median, col = webtreatment)) +
  geom_line(position = position_dodge(width = 0.1), lty = 3) +
  geom_errorbar(aes(ymin=quantile25, ymax=quantile75), 
                position = position_dodge(width = 0.1), width=.1) +
  geom_point(position = position_dodge(width = 0.1)) +
  scale_color_manual(values = c("forestgreen", "blue4", "lightskyblue")) +
  scale_x_continuous(name = "Spatial resource overlap", 
                     breaks = c(.2, .4, .6, .8, 1), labels = c(1, "", "", "", 0),
                     trans = "reverse") +
  scale_y_continuous(name = "Plant productivity", breaks = c(700, 1200), limits = c(700, 1200)) +
  facet_grid(.~factor(nP, labels = c("monoculture", "16-species mixture"))) +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.position = "none")


ggsave(paste0(path.to.data, "/FIGS1_raw.png"), 
       plot = FIGS1, width = 150, height = 90, unit = "mm", dpi = 600)




# _ 2.4. Fig.S2 herbivory ####

FIGS2 <- data_summary(dat[dat$prodP_sum != 0 & dat$webtreatment != "none", ], 
                      varname = "fluxHin", 
                      groupnames = c("nP", "N_ikk", "webtreatment")) %>%
  ggplot(aes(nP, median, col = webtreatment)) +
  geom_line(position = position_dodge(width = 4)) +
  geom_errorbar(aes(ymin=quantile25, ymax=quantile75), 
                position = position_dodge(width = 4), width=4) +
  geom_point(position = position_dodge(width = 4)) +
  facet_grid(.~factor(N_ikk, levels = seq(1, 0.2, by = -0.2))) +
  scale_color_manual(values = c("blue4", "lightskyblue")) +
  scale_x_continuous(name = "Plant species richness", 
                     limits = c(-3, 20), breaks = c(1, 16)) +
  scale_y_continuous(name = "Herbivory") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")


ggsave(paste0(path.to.data, "/FIGS2_raw.png"), 
       plot = FIGS2, width = 173, height = 70, unit = "mm", dpi = 600)



# _ 2.5. Fig.S3 total biomass ####

FIGS3 <- data_summary(dat[dat$prodP_sum != 0, ], 
                      varname = "dens", 
                      groupnames = c("nP", "N_ikk", "webtreatment")) %>%
  ggplot(aes(nP, median, col = webtreatment)) +
  geom_line(position = position_dodge(width = 4)) +
  geom_errorbar(aes(ymin=quantile25, ymax=quantile75), 
                position = position_dodge(width = 4), width=4) +
  geom_point(position = position_dodge(width = 4)) +
  facet_grid(.~factor(N_ikk, levels = seq(1, 0.2, by = -0.2))) +
  scale_color_manual(values = c("forestgreen", "blue4", "lightskyblue")) +
  scale_x_continuous(name = "Plant species richness", 
                     limits = c(-3, 20), breaks = c(1, 16)) +
  scale_y_continuous(name = "Total biomass") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

ggsave(paste0(path.to.data, "/FIGS3_raw.png"), 
       plot = FIGS3, width = 173, height = 70, unit = "mm", dpi = 600)



# _ 2.6. Fig.S4 animal diversity and survival ####

dat_nA <- rbind(
  data_summary(dat[dat$webtreatment != "none", ], 
               varname = "nA", 
               groupnames = c("nP", "N_ikk", "webtreatment")),
  data_summary(dat[dat$webtreatment != "none", ], 
               varname = "survA", 
               groupnames = c("nP", "N_ikk", "webtreatment")))

dat_nA$nAXwebtreatment <- paste(dat_nA$varname, dat_nA$webtreatment)
dat_nA$nAXwebtreatment <- factor(dat_nA$nAXwebtreatment,
                                 levels = c("nA non-nested", "nA nested", "survA non-nested", "survA nested"))


FIGS4 <- ggplot(dat_nA, aes(nP, median, col = nAXwebtreatment)) +
  geom_line(aes(lty = varname), 
            position = position_dodge(width = 2)) +
  geom_errorbar(aes(ymin=quantile25, ymax=quantile75),
                position = position_dodge(width = 2), width=2) +
  geom_point(position = position_dodge(width = 2)) +
  scale_color_manual(values = c("blue4", "lightskyblue", "blue4", "lightskyblue")) +
  scale_x_continuous(name = "Plant species richness", 
                     breaks = c(1, 16), limits = c(0, 17)) +
  scale_y_continuous(name = "Animal species richness", limits = c()) +
  scale_linetype_manual(values = c(2, 1)) +
  
  facet_grid(.~factor(N_ikk, levels = seq(1, 0.2, by = -0.2))) +
  
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

ggsave(paste0(path.to.data, "/FIGS4_raw.png"), 
       plot = FIGS4, width = 173, height = 70, unit = "mm", dpi = 600)

