source('../utils.R')
usePackage("ggplot2")
usePackage("dplyr")

##############################################################################################
################################## Read Results ##############################################
##############################################################################################

resfolder <- './results'
files <- list.files(resfolder)

res1 <- read.csv(paste0(resfolder, '/', files[1]))
res1 <- na.omit(res1)
results <- data.frame(matrix(NA, nrow=nrow(res1)*length(files), ncol=ncol(res1)))
names(results) <- names(res1)

i <- 1
for(file in files){
  res <- read.csv(paste0(resfolder, '/', file))
  res <- res %>% filter(!is.na(stop_time))
  results[i:(i+nrow(res)-1),] <- res
  i <- i + nrow(res)
  print((i-1)/nrow(res))
}
results <- results %>% arrange(sim)

##############################################################################################
####################################### Plots ################################################
##############################################################################################

hr <- results %>%
        group_by(lambda1_maj, lambda1) %>% 
         summarise(HR=mean(HR, na.rm=TRUE))

plot_tte <- results %>% group_by(type, ncolX, lambda1_maj, lambda1, eff.cols, method) %>% 
              summarise(stop1 = mean(stop_time <= 12 & stop_type=="Decrease"), 
                        stop2 = mean(stop_time <= 18 & stop_type=="Decrease"),
                        stop3 = mean(stop_time <= 24 & stop_type=="Decrease"), 
                        var1 = mean(stop_time <= 12 & stop_type=="Decrease")*(1-mean(stop_time <= 12 & stop_type=="Decrease")) / n(), 
                        var2 = mean(stop_time <= 18 & stop_type=="Decrease")*(1-mean(stop_time <= 18 & stop_type=="Decrease")) / n(),
                        var3 = mean(stop_time <= 24 & stop_type=="Decrease")*(1-mean(stop_time <= 24 & stop_type=="Decrease")) / n())
plot_tte <- plot_tte %>% left_join(hr, by=c("lambda1_maj", "lambda1"))
plot_tte$base <- ifelse(plot_tte$lambda1_maj == 0.1, 'Treatment Has No Effect on Majority', 'Treatment Benefits Majority')

plot_tte <- plot_tte %>% ungroup() %>% 
              dplyr::select(type, ncolX, eff.cols, base, HR, method, 
                            stop1, stop2, stop3, 
                            var1, var2, var3)

plot_tte$dashed <- ifelse(plot_tte$type=="known", 'yes', 'no')
plot_tte$method <- factor(plot_tte$method, levels = c('OF', 'Pocock'))
levels(plot_tte$method) <- c("O'Brien-Fleming", "Pocock")
plot_tte$type <- factor(plot_tte$type, levels = c('known', 'est', 'unknown'))

############### Figure S10 #################################

ggplot(plot_tte %>% mutate(size = paste0("Minority Group Size: ",round(100/(2*eff.cols),0), '%')) %>% 
         filter(eff.cols %in% c(1,2),
                ncolX == 5,
                method=="O'Brien-Fleming"),
       aes(x=HR)) + 
  facet_grid(cols=vars(size),rows = vars(base)) + 
  geom_ribbon(aes(ymin=stop3-1.96*sqrt(var3),
                  ymax=stop3+1.96*sqrt(var3), 
                  fill=type, 
                  linetype=dashed), alpha=0.2, 
              outline.type= "full") +
  geom_line(aes(y=stop3, col=type,linetype=dashed)) + 
  geom_point(aes(y=stop3, col=type)) + 
  theme_classic() + 
  theme(legend.title = element_blank()) + 
  scale_fill_manual(values = c("unknown" = "#619CFF",
                               "est" ="#F8766D",
                               "known"="#00BA38"), 
                    labels = c("Oracle","CLASH", "Homogeneous")) + 
  scale_color_manual(values = c("unknown" = "#619CFF",
                                "est" ="#F8766D",
                                "known"="#00BA38"), 
                     labels = c("Oracle", "CLASH", "Homogeneous")) +
  xlab('Hazard Ratio for Minority Group') + 
  ylab('Probability of Stopping Early (at any interim checkpoint)') + 
  scale_linetype(guide = 'none') + 
  theme(legend.position = 'bottom', text = element_text(size = 16))  + 
  scale_x_continuous(breaks = seq(0,1,0.2)) + 
  scale_y_continuous(breaks = seq(0,1,0.1))

############### Figure S11 #################################

plot_improvement <- plot_tte %>% 
  pivot_longer(cols = -c(1:6,13), 
               names_pattern = "([a-z]+)([0-9]+)$", 
               names_to = c("variable", "time")) %>% 
  pivot_wider(names_from = "variable",values_from = "value") %>% 
  filter(type %in% c("unknown", "est")) %>% 
  pivot_wider(names_from = "type", 
              values_from = c("stop", "var"))%>% 
  mutate(improvement = stop_est - stop_unknown, 
         variance = var_est + var_unknown) %>% 
  mutate(time = paste0("Checkpoint #", as.numeric(time)), 
         ncolX = paste0(ncolX, " covariates"),
         ncolX = factor(ncolX, levels = paste0(c(3,5,10), " covariates")))

# Figure S11a
ggplot(plot_improvement %>% filter(base == 'Treatment Benefits Majority',
                                   eff.cols==2, 
                                   method=="O'Brien-Fleming",
                                   ncolX %in% paste0(c(3, 5,10)," covariates")),
       aes(x=HR)) + 
  facet_grid(rows=vars(ncolX), cols = vars(time)) + 
  geom_bar(aes(y = improvement),stat='identity', fill='lightcyan2') + 
  geom_errorbar(aes(ymin=improvement-1.96*sqrt(variance),
                    ymax=improvement+1.96*sqrt(variance)), col='black') +
  theme_classic() + 
  ylab("Increase in Stopping Probability (CLASH vs. Homogeneous)") + 
  xlab('Hazard Ratio for Minority Group') + 
  scale_x_continuous(breaks=seq(0,1,0.2)) + 
  scale_y_continuous(breaks=seq(0,1,0.1), limits =c(-0.1,1)) + 
  theme(text = element_text(size = 18), 
        axis.title.y = element_text(size = 14), 
        axis.text.y = element_text(size = 12))

# Figure S11b
ggplot(plot_improvement %>% filter(base == 'Treatment Benefits Majority',
                                   eff.cols==1, 
                                   method=="O'Brien-Fleming",
                                   ncolX %in% paste0(c(3, 5,10)," covariates")),
       aes(x=HR)) + 
  facet_grid(rows=vars(ncolX), cols = vars(time)) + 
  geom_bar(aes(y = improvement),stat='identity', fill='lightcyan2') + 
  geom_errorbar(aes(ymin=improvement-1.96*sqrt(variance),
                    ymax=improvement+1.96*sqrt(variance)), col='black') +
  theme_classic() + 
  ylab("Increase in Stopping Probability (CLASH vs. Homogeneous)") + 
  xlab('Hazard Ratio for Minority Group') + 
  scale_x_continuous(breaks=seq(0,1,0.2)) + 
  scale_y_continuous(breaks=seq(0,1,0.1), limits =c(-0.1,1)) + 
  theme(text = element_text(size = 18), 
        axis.title.y = element_text(size = 14), 
        axis.text.y = element_text(size = 12))

############### Figure S12 #################################

# Figure S12a
ggplot(plot_improvement %>% filter(base == 'Treatment Has No Effect on Majority',
                                   eff.cols==2, 
                                   method=="O'Brien-Fleming",
                                   ncolX %in% paste0(c(3, 5,10)," covariates")),
       aes(x=HR)) + 
  facet_grid(rows=vars(ncolX), cols = vars(time)) + 
  geom_bar(aes(y = improvement),stat='identity', fill='lightcyan2') + 
  geom_errorbar(aes(ymin=improvement-1.96*sqrt(variance),
                    ymax=improvement+1.96*sqrt(variance)), col='black') +
  theme_classic() + 
  ylab("Increase in Stopping Probability (CLASH vs. Homogeneous)") + 
  xlab('Hazard Ratio for Minority Group') + 
  scale_x_continuous(breaks=seq(0,1,0.2)) + 
  scale_y_continuous(breaks=seq(0,1,0.1), limits =c(-0.1,1)) + 
  theme(text = element_text(size = 18), 
        axis.title.y = element_text(size = 14), 
        axis.text.y = element_text(size = 12))

# Figure S12b
ggplot(plot_improvement %>% filter(base == 'Treatment Has No Effect on Majority',
                                   eff.cols==1, 
                                   method=="O'Brien-Fleming",
                                   ncolX %in% paste0(c(3, 5,10)," covariates")),
       aes(x=HR)) + 
  facet_grid(rows=vars(ncolX), cols = vars(time)) + 
  geom_bar(aes(y = improvement),stat='identity', fill='lightcyan2') + 
  geom_errorbar(aes(ymin=improvement-1.96*sqrt(variance),
                    ymax=improvement+1.96*sqrt(variance)), col='black') +
  theme_classic() + 
  ylab("Increase in Stopping Probability (CLASH vs. Homogeneous)") + 
  xlab('Hazard Ratio for Minority Group') + 
  scale_x_continuous(breaks=seq(0,1,0.2)) + 
  scale_y_continuous(breaks=seq(0,1,0.1), limits =c(-0.1,1)) + 
  theme(text = element_text(size = 18), 
        axis.title.y = element_text(size = 14), 
        axis.text.y = element_text(size = 12))

