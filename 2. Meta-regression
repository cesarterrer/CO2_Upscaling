library(devtools)
install_github("cjvanlissa/metaforest")
install_github("wviechtb/metafor")
source("packages.R")
source("colours.R")
data <- read.csv("AGB_effects2019.csv",na.strings=c("",NA))
data <- data[complete.cases(data$id),]
data$obs <- 1:nrow(data)
levels(data$Biome) <- list("Bo" = "Boreal_Forest","Cr"="Cropland", "Gr"="Grassland", "Sh"="Shrubland",
                           "Te"="Temperate_Forest","Tr"="Tropical_Forest")


# Use reciprocal transformations for variables that show ceiling/floor effects
data$CNr_recip <- 1/data$CNr
data$P_recip <- 1/data$P
data$MAP_recip <- 1/data$MAP2
data$MAT_recip <- 1/data$MAT

# Drop unused levels, otherwise rma() goes haywire
data$Fumigation.type <- droplevels(data$Fumigation.type)
contrasts(data$Fumigation.type) <- cbind(
  FaceG_OTC = c(-1, -1,  2),
  Face_G    = c(-1,  1,  0))

contrasts(data$Myc) <- cbind(
  AM_ECM = c(-1, 1))

final_aicc <- rma(yi, vi, data=data, mods= ~ Myc*CNr_recip + Myc*P_recip + Fumigation.type,
                  knha=TRUE, control=list(stepadj=0.5))
anova(final_aicc, btt=5:6) # Fumigation.type

final <- rma.mv(yi, vi, data=data, mods= ~  Myc*CNr_recip + Myc*P_recip + Fumigation.type,
       random = ~ 1 | Site.Name / obs)
final
final_results <- round(coef(summary(final)),4)
write.csv(final_results,"summary_results.csv")

############### C:N ##################
CNrange <- range(data$CNr_recip)
newCN <- seq(CNrange[1],CNrange[2],by=0.001)
am.df <- filter(data, Myc=="AM")
am.mv <- rma.mv(yi, vi, data=am.df, mods= ~ CNr_recip + P_recip + Fumigation.type,
                 random = ~ 1 | Site.Name / obs)
CNr.new <- data.frame(CNr_recip = newCN, P_recip = mean(data$P_recip), Fumigation.type=factor("FACE", levels=c("FACE","G","OTC")))
CNr_mods <- model.matrix(~ CNr_recip + P_recip + Fumigation.type, CNr.new)[,-1]
AMCNpred <- as.data.frame.list(predict(am.mv, newmods = CNr_mods, addx=T, transf=make_pct))
AMCNpred$CNr <- 1/(AMCNpred$X.CNr_recip)

AMCNplot <- ggplot(data=AMCNpred, aes(CNr, pred)) + 
  geom_point(data=am.df,aes(x = CNr,y = make_pct(yi), size=1/(am.mv$tau2+vi), shape=Fumigation.type), fill=my_bold[1], alpha=myalpha) +
  #scale_shape_discrete(solid=F) +
  scale_shape_manual(values=c(21,22,24), name=element_blank()) +
  geom_line (size=0.8, aes(col="FACE")) + 
  scale_colour_manual(values="black", name=element_blank()) +
  geom_ribbon(aes(ymax= ci.ub, ymin=ci.lb),alpha=0.1) + 
  scale_y_continuous(breaks=seq(-25,90,25), limits=c(-25, 90)) +
  xlim(6,23) + ylab(expression(paste("Biomass response to ", CO[2]," enrichment (%)", sep=""))) + xlab ("Soil C:N") +
  geom_hline(yintercept=0, linetype = "longdash") + 
  guides(size=FALSE, shape = guide_legend(order=1, override.aes = list(size=3, fill="transparent")),
         col = guide_legend(order=2, override.aes = list(size=1))) + 
  ylab(expression(paste("Biomass response to ", CO[2]," enrichment (%)", sep=""))) + xlab ("Soil C:N") +
  geom_hline(yintercept=0,linetype = "longdash") + 
  theme_classic(base_size = 12) + 
  theme(legend.position = c(1, 1), legend.justification = c(1, 1),
        legend.spacing.y = unit(0, "cm"))
AMCNplot

ecm.df <- filter(data, Myc=="ECM")
ecm.mv <- rma.mv(yi, vi, data=ecm.df, mods= ~ P_recip + CNr_recip + Fumigation.type,
                 random = ~ 1 | Site.Name / obs)
ECMCNpred <- as.data.frame.list(predict(ecm.mv, newmods = CNr_mods, addx=T, transf=make_pct))
ECMCNpred$CNr <- 1/(AMCNpred$X.CNr_recip)

ECMCNplot <- ggplot(data=ECMCNpred, aes(CNr, pred)) + 
  geom_point(data=ecm.df,aes(x = CNr,y = make_pct(yi), size=1/(ecm.mv$tau2+vi), shape=Fumigation.type), fill=my_bold[2], alpha=myalpha) +
  #scale_shape_discrete(solid=F) +
  scale_shape_manual(values=c(21,22,24), name=element_blank()) +
  geom_line (size=0.8, aes(col="FACE")) + 
  scale_colour_manual(values="black", name=element_blank()) +
  geom_ribbon(aes(ymax= ci.ub, ymin=ci.lb),alpha=0.1) + 
  scale_y_continuous(breaks=seq(-25,90,25), limits=c(-25, 90)) +
  ylab(expression(paste("Biomass response to ", CO[2]," enrichment (%)", sep=""))) + xlab ("Soil C:N") +
  geom_hline(yintercept=0, linetype = "longdash") + xlim(min(ecm.df$CNr),max(ecm.df$CNr)) +
  guides(size=FALSE, shape = guide_legend(order=1, override.aes = list(size=3, fill="transparent")),
         col = guide_legend(order=2, override.aes = list(size=1))) + 
  ylab(expression(paste("Biomass response to ", CO[2]," enrichment (%)", sep=""))) + xlab ("Soil C:N") +
  geom_hline(yintercept=0,linetype = "longdash") + 
  theme_classic(base_size = 12) + 
  theme(legend.position = c(1, 1), legend.justification = c(1, 1),
        legend.spacing.y = unit(0, "cm"))
ECMCNplot

######## P #########
Prange <- range(data$P_recip)
newP <- seq(Prange[1],Prange[2],by=0.001)
P.new <- data.frame(P_recip = newP, CNr_recip= mean(data$CNr_recip), Fumigation.type=factor("FACE", levels=c("FACE","G","OTC")))
P_mods <- model.matrix(~ P_recip + CNr_recip + Fumigation.type, P.new)[,-1]
ECMPpred <- as.data.frame.list(predict(ecm.mv, newmods = P_mods, addx=T, transf=make_pct))
ECMPpred$P <- 1/(ECMPpred$X.P_recip)

ECMPplot <- ggplot(data=ECMPpred, aes(P, pred)) + 
  geom_point(data=ecm.df,aes(x = P,y = make_pct(yi), size=1/(ecm.mv$tau2+vi), shape=Fumigation.type), fill=my_bold[2], alpha=myalpha) +
  #scale_shape_discrete(solid=F) +
  scale_shape_manual(values=c(21,22,24), name=element_blank()) +
  geom_line (size=0.8, aes(col="FACE")) + 
  scale_colour_manual(values="black", name=element_blank()) +
  geom_ribbon(aes(ymax= ci.ub, ymin=ci.lb),alpha=0.1) + 
  scale_y_continuous(breaks=seq(-25,90,25), limits=c(-25, 90)) +
  ylab(expression(paste("Biomass response to ", CO[2]," enrichment (%)", sep=""))) + xlab ("Soil available P (ppm)") +
  geom_hline(yintercept=0, linetype = "longdash") + 
  guides(size=FALSE, shape = guide_legend(order=1, override.aes = list(size=3, fill="transparent")),
         col = guide_legend(order=2, override.aes = list(size=1))) + 
  ylab(expression(paste("Biomass response to ", CO[2]," enrichment (%)", sep=""))) +
  geom_hline(yintercept=0,linetype = "longdash") + 
  theme_classic(base_size = 12) + 
  theme(legend.position = c(1, 1), legend.justification = c(1, 1),
        legend.spacing.y = unit(0, "cm"))
ECMPplot

AMPpred <- as.data.frame.list(predict(am.mv, newmods = P_mods, addx=T, transf=make_pct))
AMPpred$P <- 1/(AMPpred$X.P_recip)

AMPplot <- ggplot(data=AMPpred, aes(P, pred)) + 
  geom_point(data=am.df,aes(x = P,y = make_pct(yi), size=1/(am.mv$tau2+vi), shape=Fumigation.type), fill=my_bold[1], alpha=myalpha) +
  #scale_shape_discrete(solid=F) +
  scale_shape_manual(values=c(21,22,24), name=element_blank()) +
  geom_line (size=0.8, aes(col="FACE")) + 
  scale_colour_manual(values="black", name=element_blank()) +
  geom_ribbon(aes(ymax= ci.ub, ymin=ci.lb),alpha=0.1) + 
  scale_y_continuous(breaks=seq(-25,90,25), limits=c(-25, 90)) +
  ylab(expression(paste("Biomass response to ", CO[2]," enrichment (%)", sep=""))) + xlab ("Soil available P (ppm)") +
  geom_hline(yintercept=0, linetype = "longdash") + 
  guides(size=FALSE, shape = guide_legend(order=1, override.aes = list(size=3, fill="transparent")),
         col = guide_legend(order=2, override.aes = list(size=1))) + 
  ylab(expression(paste("Biomass response to ", CO[2]," enrichment (%)", sep=""))) +
  geom_hline(yintercept=0,linetype = "longdash") + 
  theme_classic(base_size = 12) + 
  theme(legend.position = c(1, 1), legend.justification = c(1, 1),
        legend.spacing.y = unit(0, "cm"))
AMPplot

######## MAP #########
MAPrange <- range(data$MAP_recip)
newMAP <- seq(MAPrange[1],MAPrange[2],by=0.00001)
am.map <- rma.mv(yi, vi, data=am.df, mods= ~ CNr_recip + P_recip + MAP_recip + Fumigation.type,
                random = ~ 1 | Site.Name / obs)
ecm.map <- rma.mv(yi, vi, data=ecm.df, mods= ~ CNr_recip + P_recip + MAP_recip + Fumigation.type,
                 random = ~ 1 | Site.Name / obs)
MAP.new <- data.frame(P_recip = mean(data$P_recip), CNr_recip= mean(data$CNr_recip), MAP_recip= newMAP, Fumigation.type=factor("FACE", levels=c("FACE","G","OTC")))
MAP_mods <- model.matrix(~ P_recip + CNr_recip + MAP_recip + Fumigation.type, MAP.new)[,-1]
AMMAPpred <- as.data.frame.list(predict(am.map, newmods = MAP_mods, addx=T, transf=make_pct))
AMMAPpred$MAP <- 1/(AMMAPpred$X.MAP_recip)
ECMMAPpred <- as.data.frame.list(predict(ecm.map, newmods = MAP_mods, addx=T, transf=make_pct))
ECMMAPpred$MAP <- 1/(ECMMAPpred$X.MAP_recip)

AMMAPplot <- ggplot(data=AMMAPpred, aes(MAP, pred)) + 
  geom_point(data=am.df,aes(x = MAP2,y = make_pct(yi), size=1/(am.mv$tau2+vi), shape=Fumigation.type), fill=my_bold[1], alpha=myalpha) +
  #scale_shape_discrete(solid=F) +
  scale_shape_manual(values=c(21,22,24), name=element_blank()) +
  geom_line (size=0.8, aes(col="FACE")) + 
  scale_colour_manual(values="black", name=element_blank()) +
  geom_ribbon(aes(ymax= ci.ub, ymin=ci.lb),alpha=0.1) + 
  scale_y_continuous(breaks=seq(-25,90,25), limits=c(-25, 90)) +
  ylab(expression(paste("Biomass response to ", CO[2]," enrichment (%)", sep=""))) + xlab ("MAP (mm)") +
  geom_hline(yintercept=0, linetype = "longdash") + 
  guides(size=FALSE, shape = guide_legend(order=1, override.aes = list(size=3, fill="transparent")),
         col = guide_legend(order=2, override.aes = list(size=1))) + 
  ylab(expression(paste("Biomass response to ", CO[2]," enrichment (%)", sep=""))) +
  geom_hline(yintercept=0,linetype = "longdash") + 
  theme_classic(base_size = 12) + 
  theme(legend.position = c(1, 1), legend.justification = c(1, 1),
        legend.spacing.y = unit(0, "cm"))

ECMMAPplot <- ggplot(data=ECMMAPpred, aes(MAP, pred)) + 
  geom_point(data=ecm.df,aes(x = MAP2,y = make_pct(yi), size=1/(ecm.mv$tau2+vi), shape=Fumigation.type), fill=my_bold[2], alpha=myalpha) +
  #scale_shape_discrete(solid=F) +
  scale_shape_manual(values=c(21,22,24), name=element_blank()) +
  geom_line (size=0.8, aes(col="FACE")) + 
  scale_colour_manual(values="black", name=element_blank()) +
  geom_ribbon(aes(ymax= ci.ub, ymin=ci.lb),alpha=0.1) + 
  scale_y_continuous(breaks=seq(-25,90,25), limits=c(-25, 90)) +
  ylab(expression(paste("Biomass response to ", CO[2]," enrichment (%)", sep=""))) + xlab ("MAP (mm)") +
  geom_hline(yintercept=0, linetype = "longdash") + 
  guides(size=FALSE, shape = guide_legend(order=1, override.aes = list(size=3, fill="transparent")),
         col = guide_legend(order=2, override.aes = list(size=1))) + 
  ylab(expression(paste("Biomass response to ", CO[2]," enrichment (%)", sep=""))) +
  geom_hline(yintercept=0,linetype = "longdash") + 
  theme_classic(base_size = 12) + 
  xlim(min(ecm.df$MAP2), max(ecm.df$MAP2)) +
  theme(legend.position = c(1, 1), legend.justification = c(1, 1),
        legend.spacing.y = unit(0, "cm"))


# AM-CNr / ECM-P plot #
prow <- plot_grid( AMCNplot,
                   ECMPplot+ ylab("") + theme(legend.position="none"),
                   align = 'hv',
                   labels = c("a", "b"),
                   hjust = -2,
                   nrow = 1
)
save_plot("figures/AGB_CNP.pdf", prow, ncol=2, nrow=1, useDingbats=F, base_width=3.5)
save_plot("figures/AGB_CNP.png", prow, ncol=2, nrow=1, dpi= 300, base_width=3.5)

# AM-P / ECM-CNr plot #
prow2 <- plot_grid(AMPplot + annotation_custom(xmin = -Inf, xmax = Inf, ymin=80, ymax = 80, 
                                               textGrob("AM", gp=gpar(fontface="bold"))),
                   ECMCNplot+ ylab("") + theme(legend.position="none")+
                     annotation_custom(xmin = -Inf, xmax = Inf, ymin=80, ymax = 80, 
                                       textGrob("ECM", gp=gpar(fontface="bold"))),
                   AMMAPplot + theme(legend.position="none"),
                   ECMMAPplot + ylab("") + theme(legend.position="none"),
                   align = 'hv',
                   labels = "auto",
                   hjust = -2,
                   nrow = 2
)
save_plot("figures/AGB_CNP2.pdf", prow2, ncol=2, nrow=2, useDingbats=F)
save_plot("figures/AGB_CNP2.png", prow2, ncol=2, nrow=2, dpi= 300)
