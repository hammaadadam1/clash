source('../utils.R')
usePackage("ggplot2")
usePackage("dplyr")

##############################################################################################
################################## Read Results ##############################################
##############################################################################################

resfolder <- './results_new'
files <- list.files(resfolder)

res1 <- read.csv(paste0(resfolder, '/', files[1]))
res1 <- res1 %>% filter(!is.na(stop_time))
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

plot_norm <- results %>% group_by(type, ncolX, effect, effect_maj, eff.cols, method) %>% 
              summarise(stop25 = mean(stop_time <= 500 & stop_type=="Increase"), 
                        stop50 = mean(stop_time <= 1000 & stop_type=="Increase"),
                        stop75 = mean(stop_time <= 1500 & stop_type=="Increase"), 
                        var25 = mean(stop_time <= 500 & stop_type=="Increase")*(1-mean(stop_time <= 500 & stop_type=="Increase")) / n(), 
                        var50 = mean(stop_time <= 1000 & stop_type=="Increase")*(1-mean(stop_time <= 1000 & stop_type=="Increase")) / n(),
                        var75 = mean(stop_time <= 1500 & stop_type=="Increase")*(1-mean(stop_time <= 1500 & stop_type=="Increase")) / n())

plot_norm$base <- ifelse(plot_norm$effect_maj ==0, 'Treatment Has No Effect on Majority', 'Treatment Benefits Majority')
plot_norm$dashed <- ifelse(plot_norm$type=="known", 'yes', 'no')
plot_norm$method <- factor(plot_norm$method, levels = c('OF', 'mixsprt', 'maxsprt', 'bayesian'))
levels(plot_norm$method) <- c("O'Brien-Fleming", "mSPRT", "MaxSPRT", 'Bayesian Estimation')
plot_norm$type <- factor(plot_norm$type, levels = c('known', 'est', 'subtle','unknown'))

############### Figure 2 #################################
ggplot(plot_norm %>% filter(eff.cols==3, 
                           ncolX == 5),
       aes(x=effect)) + 
  facet_grid(cols=vars(method),rows = vars(base)) + 
  geom_ribbon(aes(ymin=stop75-1.96*sqrt(var75),
                  ymax=stop75+1.96*sqrt(var75), 
                  fill=type, 
                  linetype=dashed), alpha=0.2, 
              outline.type= "full") +
  geom_line(aes(y=stop75, col=type, linetype=dashed), linewidth=1) + 
  geom_point(aes(y=stop75, col=type, shape=type), size=4) + 
  theme_classic() + 
  theme(legend.title = element_blank()) + 
  scale_fill_manual(values = c("unknown" = "#56B4E9",
                               "est" ="#D55E00",
                               "subtle" ="#CC79A7",
                               "known"="#009E73"), 
                    labels = c("Oracle","CLASH","SUBTLE", "Homogeneous")) + 
  scale_color_manual(values = c("unknown" = "#56B4E9",
                                "est" ="#D55E00",
                                "subtle" ="#CC79A7",
                                "known"="#009E73"), 
                     labels = c("Oracle", "CLASH","SUBTLE", "Homogeneous")) +
  scale_shape_manual(values = c("unknown"= 16,
                                "est"=15, 
                                "subtle"=17, 
                                "known"=5), 
                     labels = c("Oracle", "CLASH","SUBTLE", "Homogeneous")) +
  xlab('Harmful Treatment Effect on Minority Group') + 
  ylab('Probability of Stopping Early (at any checkpoint)') + 
  scale_linetype(guide = 'none') +
  theme(legend.position = 'bottom', 
        text = element_text(size = 16), 
        strip.text.x = element_text(size = 18), 
        legend.text = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18), 
        legend.key.width= unit(1, 'cm'))  + 
  scale_x_continuous(breaks = seq(0,1,0.2)) + 
  scale_y_continuous(breaks = seq(0,1,0.1))
ggsave('./figures/fig2.pdf',height=7.5, width=14)

############### Figure 3 #################################
plot_improvement_subtle <- plot_norm %>% 
                    pivot_longer(cols = -c(1:6,13:14), 
                                 names_pattern = "([a-z]+)([0-9]+)$", 
                                 names_to = c("variable", "time")) %>% 
            pivot_wider(names_from = "variable",values_from = "value") %>% 
            filter(type %in% c("subtle", "est")) %>% 
            pivot_wider(names_from = "type", 
                        values_from = c("stop", "var")) %>% 
            mutate(improvement = stop_est - stop_subtle, 
                   variance = var_est + var_subtle) %>% 
            mutate(time = paste0("Checkpoint #", as.numeric(time)/25), 
                   ncolX = paste0(ncolX, " covariates"),
                   ncolX = factor(ncolX, levels = paste0(c(3,5,10), " covariates")))

# Figure 3a
ggplot(plot_improvement_subtle %>% filter(effect_maj == -0.1,
                                   eff.cols==3, 
                                   method=='mSPRT',
                                   ncolX %in% paste0(c(3, 5,10)," covariates")),
       aes(x=effect)) + 
  facet_grid(rows=vars(ncolX), cols = vars(time)) + 
  geom_bar(aes(y = improvement),stat='identity', fill='lightcyan2') + 
  geom_errorbar(aes(ymin=improvement-1.96*sqrt(variance),
                    ymax=improvement+1.96*sqrt(variance)), col='black') +
  theme_classic() + 
  ylab("Increase in Stopping Probability (CLASH vs. SUBTLE)") + 
  xlab('Harmful Treatment Effect on Minority Group') + 
  scale_x_continuous(breaks=seq(0,1,0.2)) + 
  scale_y_continuous(breaks=seq(0,1,0.2), limits =c(-0.1,1)) + 
  theme(text = element_text(size = 18), 
        axis.title.y = element_text(size = 14), 
        axis.text.y = element_text(size = 12))
ggsave('./figures/fig3a.pdf',height=5.5, width=8)

# Figure 3b
ggplot(plot_improvement_subtle %>% filter(effect_maj == -0.1,
                                          eff.cols==2, 
                                          method=='mSPRT',
                                   ncolX %in% paste0(c(3, 5,10)," covariates")),
       aes(x=effect)) + 
  facet_grid(rows=vars(ncolX), cols = vars(time)) + 
  geom_bar(aes(y = improvement),stat='identity', fill='lightcyan2') + 
  geom_errorbar(aes(ymin=improvement-1.96*sqrt(variance),
                    ymax=improvement+1.96*sqrt(variance)), col='black') +
  theme_classic() + 
  ylab("Increase in Stopping Probability (CLASH vs. SUBTLE)") + 
  xlab('Harmful Treatment Effect on Minority Group') + 
  scale_x_continuous(breaks=seq(0,1,0.2)) + 
  scale_y_continuous(breaks=seq(0,1,0.2), limits =c(-0.1,1)) + 
  theme(text = element_text(size = 18), 
        axis.title.y = element_text(size = 14), 
        axis.text.y = element_text(size = 12))
ggsave('./figures/fig3b.pdf',height=5.5, width=8)

############### Figure S6 #################################

# Figure S6a
ggplot(plot_norm %>% filter(eff.cols==2, 
                            ncolX == 5),
       aes(x=effect)) + 
  facet_grid(cols=vars(method),rows = vars(base)) + 
  geom_ribbon(aes(ymin=stop75-1.96*sqrt(var75),
                  ymax=stop75+1.96*sqrt(var75), 
                  fill=type, 
                  linetype=dashed), alpha=0.2, 
              outline.type= "full") +
  geom_line(aes(y=stop75, col=type, linetype=dashed), linewidth=1) + 
  geom_point(aes(y=stop75, col=type, shape=type), size=4) + 
  theme_classic() + 
  theme(legend.title = element_blank()) + 
  scale_fill_manual(values = c("unknown" = "#56B4E9",
                               "est" ="#D55E00",
                               "subtle" ="#CC79A7",
                               "known"="#009E73"), 
                    labels = c("Oracle","CLASH","SUBTLE", "Homogeneous")) + 
  scale_color_manual(values = c("unknown" = "#56B4E9",
                                "est" ="#D55E00",
                                "subtle" ="#CC79A7",
                                "known"="#009E73"), 
                     labels = c("Oracle", "CLASH","SUBTLE", "Homogeneous")) +
  scale_shape_manual(values = c("unknown"= 16,
                                "est"=15, 
                                "subtle"=17, 
                                "known"=5), 
                     labels = c("Oracle", "CLASH","SUBTLE", "Homogeneous")) +
  xlab('Harmful Treatment Effect on Minority Group') + 
  ylab('Probability of Stopping Early (at any checkpoint)') + 
  scale_linetype(guide = 'none') + 
  theme(legend.position = 'bottom', 
        text = element_text(size = 16), 
        strip.text.x = element_text(size = 18), 
        legend.text = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18), 
        legend.key.width= unit(1, 'cm'))  + 
  scale_x_continuous(breaks = seq(0,1,0.2)) + 
  scale_y_continuous(breaks = seq(0,1,0.1))
ggsave('./figures/app/figs6a.pdf',height=7.5, width=14)

# Figure S6b
ggplot(plot_norm %>% filter(eff.cols==1, 
                            ncolX == 5),
       aes(x=effect)) + 
  facet_grid(cols=vars(method),rows = vars(base)) + 
  geom_ribbon(aes(ymin=stop75-1.96*sqrt(var75),
                  ymax=stop75+1.96*sqrt(var75), 
                  fill=type, 
                  linetype=dashed), alpha=0.2, 
              outline.type= "full") +
  geom_line(aes(y=stop75, col=type, linetype=dashed), linewidth=1) + 
  geom_point(aes(y=stop75, col=type, shape=type), size=4) + 
  theme_classic() + 
  theme(legend.title = element_blank()) + 
  scale_fill_manual(values = c("unknown" = "#56B4E9",
                               "est" ="#D55E00",
                               "subtle" ="#CC79A7",
                               "known"="#009E73"), 
                    labels = c("Oracle","CLASH","SUBTLE", "Homogeneous")) + 
  scale_color_manual(values = c("unknown" = "#56B4E9",
                                "est" ="#D55E00",
                                "subtle" ="#CC79A7",
                                "known"="#009E73"), 
                     labels = c("Oracle", "CLASH","SUBTLE", "Homogeneous")) +
  scale_shape_manual(values = c("unknown"= 16,
                                "est"=15, 
                                "subtle"=17, 
                                "known"=5), 
                     labels = c("Oracle", "CLASH","SUBTLE", "Homogeneous")) +
  xlab('Harmful Treatment Effect on Minority Group') + 
  ylab('Probability of Stopping Early (at any checkpoint)') + 
  scale_linetype(guide = 'none') + 
  theme(legend.position = 'bottom', 
        text = element_text(size = 16), 
        strip.text.x = element_text(size = 18), 
        legend.text = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18), 
        legend.key.width= unit(1, 'cm'))  + 
  scale_x_continuous(breaks = seq(0,1,0.2)) + 
  scale_y_continuous(breaks = seq(0,1,0.1))
ggsave('./figures/app/figs6b.pdf',height=7.5, width=14)

############### Figure S7 #################################
# Figure S7a
ggplot(plot_improvement_subtle %>% filter(effect_maj == 0,
                                          eff.cols==3, 
                                          method=='mSPRT',
                                          ncolX %in% paste0(c(3, 5,10)," covariates")),
       aes(x=effect)) + 
  facet_grid(rows=vars(ncolX), cols = vars(time)) + 
  geom_bar(aes(y = improvement),stat='identity', fill='lightcyan2') + 
  geom_errorbar(aes(ymin=improvement-1.96*sqrt(variance),
                    ymax=improvement+1.96*sqrt(variance)), col='black') +
  theme_classic() + 
  ylab("Increase in Stopping Probability (CLASH vs. SUBTLE)") + 
  xlab('Harmful Treatment Effect on Minority Group') + 
  scale_x_continuous(breaks=seq(-0.2,1,0.2)) + 
  scale_y_continuous(breaks=seq(-0.2,1,0.2), limits =c(-0.2,1)) + 
  theme(text = element_text(size = 18), 
        axis.title.y = element_text(size = 14), 
        axis.text.y = element_text(size = 12))
ggsave('./figures/app/figs7a.pdf',height=5.5, width=8)

# Figure S7b
ggplot(plot_improvement_subtle %>% filter(effect_maj == 0,
                                          eff.cols==2, 
                                          method=='mSPRT',
                                          ncolX %in% paste0(c(3, 5,10)," covariates")),
       aes(x=effect)) + 
  facet_grid(rows=vars(ncolX), cols = vars(time)) + 
  geom_bar(aes(y = improvement),stat='identity', fill='lightcyan2') + 
  geom_errorbar(aes(ymin=improvement-1.96*sqrt(variance),
                    ymax=improvement+1.96*sqrt(variance)), col='black') +
  theme_classic() + 
  ylab("Increase in Stopping Probability (CLASH vs. SUBTLE)") + 
  xlab('Harmful Treatment Effect on Minority Group') + 
  scale_x_continuous(breaks=seq(-0.2,1,0.2)) + 
  scale_y_continuous(breaks=seq(-0.2,1,0.2), limits =c(-0.2,1)) + 
  theme(text = element_text(size = 18), 
        axis.title.y = element_text(size = 14), 
        axis.text.y = element_text(size = 12))
ggsave('./figures/app/figs7b.pdf',height=5.5, width=8)

############### Figure S8 #################################
plot_improvement_homo <- plot_norm %>% 
  pivot_longer(cols = -c(1:6,13:14), 
               names_pattern = "([a-z]+)([0-9]+)$", 
               names_to = c("variable", "time")) %>% 
  pivot_wider(names_from = "variable",values_from = "value") %>% 
  filter(type %in% c("unknown", "est")) %>% 
  pivot_wider(names_from = "type", 
              values_from = c("stop", "var")) %>% 
  mutate(improvement = stop_est - stop_unknown, 
         variance = var_est + var_unknown) %>% 
  mutate(time = paste0("Checkpoint #", as.numeric(time)/25), 
         ncolX = paste0(ncolX, " covariates"),
         ncolX = factor(ncolX, levels = paste0(c(3,5,10), " covariates")))

# Figure S8a
ggplot(plot_improvement_homo %>% filter(effect_maj == -0.1,
                                          eff.cols==3, 
                                          method=="O'Brien-Fleming",
                                          ncolX %in% paste0(c(3, 5,10)," covariates")),
       aes(x=effect)) + 
  facet_grid(rows=vars(ncolX), cols = vars(time)) + 
  geom_bar(aes(y = improvement),stat='identity', fill='lightcyan2') + 
  geom_errorbar(aes(ymin=improvement-1.96*sqrt(variance),
                    ymax=improvement+1.96*sqrt(variance)), col='black') +
  theme_classic() + 
  ylab("Increase in Stopping Probability (CLASH vs. Homogeneous)") + 
  xlab('Harmful Treatment Effect on Minority Group') + 
  scale_x_continuous(breaks=seq(0,1,0.2)) + 
  scale_y_continuous(breaks=seq(0,1,0.2), limits =c(-0.1,1)) + 
  theme(text = element_text(size = 18), 
        axis.title.y = element_text(size = 14), 
        axis.text.y = element_text(size = 12))
ggsave('./figures/app/figs8a.pdf',height=5.5, width=8)

# Figure S8b
ggplot(plot_improvement_homo %>% filter(effect_maj == -0.1,
                                        eff.cols==2, 
                                        method=="O'Brien-Fleming",
                                        ncolX %in% paste0(c(3, 5,10)," covariates")),
       aes(x=effect)) + 
  facet_grid(rows=vars(ncolX), cols = vars(time)) + 
  geom_bar(aes(y = improvement),stat='identity', fill='lightcyan2') + 
  geom_errorbar(aes(ymin=improvement-1.96*sqrt(variance),
                    ymax=improvement+1.96*sqrt(variance)), col='black') +
  theme_classic() + 
  ylab("Increase in Stopping Probability (CLASH vs. Homogeneous)") + 
  xlab('Harmful Treatment Effect on Minority Group') + 
  scale_x_continuous(breaks=seq(0,1,0.2)) + 
  scale_y_continuous(breaks=seq(0,1,0.2), limits =c(-0.1,1)) + 
  theme(text = element_text(size = 18), 
        axis.title.y = element_text(size = 14), 
        axis.text.y = element_text(size = 12))
ggsave('./figures/app/figs8b.pdf',height=5.5, width=8)

############### Figure S9 #################################

# Figure S9a
ggplot(plot_improvement_homo %>% filter(effect_maj == 0,
                                        eff.cols==3, 
                                        method=="O'Brien-Fleming",
                                        ncolX %in% paste0(c(3, 5,10)," covariates")),
       aes(x=effect)) + 
  facet_grid(rows=vars(ncolX), cols = vars(time)) + 
  geom_bar(aes(y = improvement),stat='identity', fill='lightcyan2') + 
  geom_errorbar(aes(ymin=improvement-1.96*sqrt(variance),
                    ymax=improvement+1.96*sqrt(variance)), col='black') +
  theme_classic() + 
  ylab("Increase in Stopping Probability (CLASH vs. Homogeneous)") + 
  xlab('Harmful Treatment Effect on Minority Group') + 
  scale_x_continuous(breaks=seq(0,1,0.2)) + 
  scale_y_continuous(breaks=seq(0,1,0.2), limits =c(-0.1,1)) + 
  theme(text = element_text(size = 18), 
        axis.title.y = element_text(size = 14), 
        axis.text.y = element_text(size = 12))
ggsave('./figures/app/figs9a.pdf',height=5.5, width=8)

# Figure S9b
ggplot(plot_improvement_homo %>% filter(effect_maj == 0,
                                        eff.cols==2, 
                                        method=="O'Brien-Fleming",
                                        ncolX %in% paste0(c(3, 5,10)," covariates")),
       aes(x=effect)) + 
  facet_grid(rows=vars(ncolX), cols = vars(time)) + 
  geom_bar(aes(y = improvement),stat='identity', fill='lightcyan2') + 
  geom_errorbar(aes(ymin=improvement-1.96*sqrt(variance),
                    ymax=improvement+1.96*sqrt(variance)), col='black') +
  theme_classic() + 
  ylab("Increase in Stopping Probability (CLASH vs. Homogeneous)") + 
  xlab('Harmful Treatment Effect on Minority Group') + 
  scale_x_continuous(breaks=seq(0,1,0.2)) + 
  scale_y_continuous(breaks=seq(0,1,0.2), limits =c(-0.1,1)) + 
  theme(text = element_text(size = 18), 
        axis.title.y = element_text(size = 14), 
        axis.text.y = element_text(size = 12))
ggsave('./figures/app/figs9b.pdf',height=5.5, width=8)


