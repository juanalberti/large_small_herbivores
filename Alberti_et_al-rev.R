# REQUIRED PACKAGES ####

# Data manipulation and visualization
library(dplyr) 
library(tidyr)
library(ggplot2)
library(forcats)
library(cowplot)
library(paletteer)

# Statistical modeling and evaluation
library(glmmTMB)   # Generalized Linear Mixed Models
library(performance)# Model diagnostics (check_model)
library(car)       # Wald chi-square tests (Anova)
library(emmeans)   # Post-hoc pairwise comparisons
library(ggeffects) # Extract marginal effects/predictions
library(vegan)     # Multivariate analysis for community ecology


# 1. PLANT RICHNESS AND COVER OVER TIME ####


# Read richness and cover data
riqueza <- read.csv("rich_cov.csv") |> 
  mutate(
    yr = as.character(yr),
    Plot = as.character(Plot),
    Fenced_Area = ifelse(Herbivory == "Control", 
                         paste0("fence_n_", Plot), 
                         paste0("fence_y_", Plot))
  )


# Stats for richness

# Full model: Gaussian distribution (log link), testing interaction between year and Herbivory.
# Includes random intercept for Plot and an AR1 autocorrelation structure for repeated measures.
modelo0 <- glmmTMB(
  riq ~ yr * Herbivory + (1 | Plot)+ (1|Fenced_Area) + ar1(yr + 0 | uidplot),
  family = gaussian("log"),
  data = riqueza,
  REML = T
)

# Test if removing the AR1 autocorrelation structure improves/worsens model fit
modelo1 <- update(modelo0, . ~ . - ar1(yr + 0 | uidplot))
anova(modelo0, modelo1)

# Test if removing the random intercept (Plot) is justified
modelo2 <- update(modelo0, . ~ . - (1 | Plot))
anova(modelo0, modelo2)

# Test if removing the random intercept (Fenced_Area) is justified
modelo3 <- update(modelo0, . ~ . - (1 | Fenced_Area))
anova(modelo0, modelo3)

# Removing Fenced_Area is fine, let's re-check the others
modelo1 <- update(modelo3, . ~ . - ar1(yr + 0 | uidplot))
anova(modelo3, modelo1)

# Test if removing the random intercept (Plot) is justified
modelo2 <- update(modelo3, . ~ . - (1 | Plot))
anova(modelo3, modelo2)

# Model diagnostics and final evaluation on the selected model (modelo0)
check_model(modelo3) # Visual check of assumptions
Anova(modelo3)       # Type II Wald chi-square tests for fixed effects
emmeans(modelo3, pairwise ~ Herbivory | yr) # Pairwise comparisons by year

# Manually created dataframe containing Tukey post-hoc significance letters for plotting
tukey_riq <- data.frame(
  yr = c(rep(2021:2025, each = 3)),
  Herbivory = rep(c("Control", "Partial", "Exclosure"), 5),
  value = c(14, 13, 11, 
            12.2, 12, 8,
            14.3, 14, 11.5,
            14.3, 15.4, 11.7, 
            16, 13, 9.5),
  txt = c(
    "a","ab","b",
    "a","a","b",
    "a","b*","c",
    "a","a","b",
    "a","b","c"
  ),
  name = "riq"
)


# Stats for cover

# AR1 structure evaluated. 
modelo0 <- glmmTMB(
  cov ~ yr * Herbivory + ar1(yr + 0 | uidplot)+ (1 | Plot) + (1|Fenced_Area), 
  data = riqueza,
  REML = T
)

# Structure selection via Likelihood Ratio Tests
modelo1 <- update(modelo0, . ~ . - ar1(yr + 0 | uidplot))
anova(modelo0, modelo1)
modelo2 <- update(modelo0, . ~ . - (1 | Plot))
anova(modelo0, modelo2)
modelo3 <- update(modelo0, . ~ . - (1 | Fenced_Area))
anova(modelo0, modelo3)

# Removing Fenced_Area is fine, let's re-check the others
modelo1 <- update(modelo3, . ~ . - ar1(yr + 0 | uidplot))
anova(modelo3, modelo1)

# Test if removing the random intercept (Plot) is justified
modelo2 <- update(modelo3, . ~ . - (1 | Plot))
anova(modelo3, modelo2)


# Final evaluation
check_model(modelo2)
Anova(modelo2)
# Post-hoc tests for main effects (removing interaction for global comparisons)
emmeans(modelo2 |> update(. ~ . - yr:Herbivory), pairwise ~ Herbivory)
emmeans(modelo2 |> update(. ~ . - yr:Herbivory), pairwise ~ yr)

tukey_cov <- NULL # No significant interaction letters to plot for cover

# Combine formatting dataframes for plotting richness and cover together
tukey <- tukey_riq %>%
  bind_rows(tukey_cov)

riqueza_l <- riqueza %>%
  pivot_longer(cols = c(riq, cov))

# Plot: Figure 2 (Richness and Cover)
ggplot(
  riqueza_l |>
    ungroup() |>
    mutate(Herbivory = fct_relevel(Herbivory, "Control", "Partial")),
  aes(as.numeric(yr), value, color = Herbivory, group = Herbivory)
) +
  stat_summary(position = position_dodge(width = 0.6), size = 1.15) +
  theme_cowplot() +
  geom_text(
    data = tukey |>
      ungroup() |>
      mutate(Herbivory = fct_relevel(Herbivory, "Control", "Partial")),
    aes(yr, value, label = txt, color = Herbivory),
    inherit.aes = F,
    show.legend = FALSE,
    position = position_dodge(width = 0.6),
    size = 7
  ) +
  scale_color_brewer(
    palette = "Set1",
    breaks = c("Control", "Partial", "Exclosure"),
    labels = c("Large + small herbivores", "Large excluded", "All excluded")
  ) +
  facet_wrap(
    . ~ factor(name, levels = c("riq", "cov")),
    ncol = 1,
    scales = "free",
    labeller = as_labeller(c(cov = "Plant cover (%)", riq = "Plant richness"))
  ) +
  labs(x = "", y = "", color = "") +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "white", color = "black"),
    legend.justification = "center",
    text = element_text(size = 20),
    axis.text = element_text(size = 16))



# 2. PLANT LITTER ####

litter <- read.csv("litter.csv") |> 
  mutate(
    yr = as.character(yr),
    Plot = as.character(Plot),
    Fenced_Area = ifelse(Herbivory == "Control", 
                         paste0("fence_n_", Plot), 
                         paste0("fence_y_", Plot))
  )

# Model selection for log10-transformed litter cover
modelo0 <- glmmTMB(
  log10(Litter) ~ yr * Herbivory + ar1(yr + 0 | uidplot) + (1 | Plot) + (1 | Fenced_Area),
  data = litter,
  REML = T)

# Structure selection via Likelihood Ratio Tests
modelo1 <- update(modelo0, . ~ . - ar1(yr + 0 | uidplot))
anova(modelo0, modelo1)
modelo2 <- update(modelo0, . ~ . - (1 | Plot))
anova(modelo0, modelo2)
modelo3 <- update(modelo0, . ~ . - (1 | Fenced_Area))
anova(modelo0, modelo3)

# Removing Fenced_Area is fine, let's re-check the others
modelo1 <- update(modelo3, . ~ . - ar1(yr + 0 | uidplot))
anova(modelo3, modelo1)

# Test if removing the random intercept (Plot) is justified
modelo2 <- update(modelo3, . ~ . - (1 | Plot))
anova(modelo3, modelo2)


check_model(modelo3)
Anova(modelo3)
emmeans(modelo3, pairwise ~ Herbivory|yr)

# Significance letters dataframe
tukey_lit <- data.frame(
  yr = c(rep(2021:2025, each = 3)),
  Herbivory = rep(c("Control", "Partial", "Exclosure"), 5),
  cover = c(11, 14, 23, 
            14, 17, 29,
            13.5, 21, 46,
            13, 16, 28, 
            15, 20, 55),
  txt = c(
    "a", "ab", "b*",
    "a","ab","b",
    "a","a","b",
    "a","ab","b",
    "a","a","b"
  ),
  name = "Litter"
)


# Plot: Figure S1 (Litter cover)

ggplot(
  litter |> 
    mutate(Herbivory = fct_relevel(Herbivory, "Control", "Partial")),
  aes(as.numeric(yr), Litter, color = Herbivory, group = Herbivory)
) +
  stat_summary(position = position_dodge(width = 0.6), size = 1.15) +
  theme_cowplot() +
  geom_text(
    data = tukey_lit |>
      ungroup() |>
      mutate(Herbivory = fct_relevel(Herbivory, "Control", "Partial")),
    aes(yr, cover, label = txt, color = Herbivory),
    inherit.aes = F,
    show.legend = FALSE,
    position = position_dodge(width = 0.6),
    size = 7
  ) +
  scale_color_brewer(
    palette = "Set1",
    breaks = c("Control", "Partial", "Exclosure"),
    labels = c("Large + small herbivores", "Large excluded", "All excluded")
  ) +
  labs(x = "", y = "Plant litter cover (%)", color = "") +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "white", color = "black"),
    legend.justification = "center",
    text = element_text(size = 20),
    axis.text = element_text(size = 16))



# 3. SPATIAL HETEROGENEITY IN COMPOSITION ####

beta <- read.csv("beta.csv") |> 
  mutate(
    yr = as.character(yr),
    Plot = as.character(Plot),
    Fenced_Area = ifelse(Herbivory == "Control", 
                         paste0("fence_n_", Plot), 
                         paste0("fence_y_", Plot))
  )

# Replace NAs and isolate community matrix
beta[is.na(beta)] <- 0
comu <- beta[, 6:(ncol(beta)-1)]

# Calculate Bray-Curtis distances on fourth-root transformed data (downweighting dominant species)
vd <- vegdist(comu^(1/4), method = "bray")
# Non-metric Multidimensional Scaling (NMDS) for visualization
mds <- metaMDS(vd, trymax = 1000, k = 3)
mds$stress

coordenadas_nmds <- as.data.frame(scores(mds, display = "sites"))

datos_grafico <- beta %>%
  mutate(
    MDS1 = coordenadas_nmds$NMDS1,
    MDS2 = coordenadas_nmds$NMDS2,
    Herbivory = fct_relevel(Herbivory, "Control", "Partial")
  )


# Plot: Figure S5 (NMDS Community Composition)

ggplot(datos_grafico, aes(x = MDS1, y = MDS2, color = Herbivory)) +
  geom_point() +
  stat_ellipse() +
  theme_cowplot() +
  scale_color_brewer(
    palette = "Set1",
    breaks = c("Control", "Partial", "Exclosure"),
    labels = c("Large + small herbivores", "Large excluded", "All excluded")
  ) +
  theme(
    legend.position = "top",
    legend.justification = "center",
    text = element_text(size = 20),
    axis.text = element_text(size = 16)
  ) +
  labs(color = "")


# Statistical testing for Spatial Heterogeneity

# Calculate multivariate dispersion (distance to group centroid) as a measure of beta diversity
beta$dist <- betadisper(
  vegdist(comu^(1 / 4), method = "bray"),
  paste0(beta$yr, beta$Herbivory)
)$distances

# Model selection for distance to centroid
modelo0 <- glmmTMB(
  dist ~ yr * Herbivory + ar1(yr + 0 | uidplot)+(1 | Plot)+(1|Fenced_Area),
  family = gaussian("log"),
  data = beta,
  REML = T)
modelo1 <- update(modelo0, . ~ . - ar1(yr + 0 | uidplot))
anova(modelo0, modelo1)
modelo2 <- update(modelo0, . ~ . - (1 | Plot))
anova(modelo0, modelo2)
modelo3 <- update(modelo0, . ~ . - (1 | Fenced_Area))
anova(modelo0, modelo3)

# Removing Fenced_Area is fine, let's re-check the others
modelo1 <- update(modelo3, . ~ . - ar1(yr + 0 | uidplot))
anova(modelo3, modelo1)

# Test if removing the random intercept (Plot) is justified
modelo2 <- update(modelo3, . ~ . - (1 | Plot))
anova(modelo3, modelo2)

check_model(modelo2)
Anova(modelo2)
emmeans(modelo2, pairwise ~ Herbivory)


# Plot: Figure 3 (Spatial Heterogeneity over time)
ggplot(
  beta |>
    mutate(Herbivory = fct_relevel(Herbivory, "Control", "Partial")),
  aes(as.numeric(yr), dist, color = Herbivory, group = Herbivory)
) +
  stat_summary(position = position_dodge(width = 0.6), size = 1.15) +
  theme_cowplot() +
  scale_color_brewer(
    palette = "Set1",
    breaks = c("Control", "Partial", "Exclosure"),
    labels = c("Large + small herbivores", "Large excluded", "All excluded")
  ) +
  labs(x = "", y = "Spatial heterogeneity in species composition", color = "") +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "white", color = "black"),
    legend.justification = "center",
    text = element_text(size = 20),
    axis.text = element_text(size = 16))



# 4. LIGHT INTERCEPTION (PAR) ####

par <- read.csv("par.csv") |> 
  mutate(
    Plot = as.character(Plot),
    Fenced_Area = ifelse(Herbivory == "Control", 
                         paste0("fence_n_", Plot), 
                         paste0("fence_y_", Plot))
  )

# Single-year evaluation, therefore time and AR1 are excluded
modelo0 <- glmmTMB(
  inter ~ Herbivory + (1 | Plot) + (1|Fenced_Area),
  data = par, REML = T
)
modelo1 <- update(modelo0, . ~ . - (1 | Plot))
anova(modelo0, modelo1)
modelo2 <- update(modelo0, . ~ . - (1 | Fenced_Area))
anova(modelo0, modelo2)

modelo3 <- update(modelo2, . ~ . - (1 | Plot))
anova(modelo2, modelo3)



check_model(modelo2)
Anova(modelo2)
emmeans(modelo2, pairwise ~ Herbivory)


# Plot: Figure S2 (PAR interception)

tukey_par <- data.frame(
  Herbivory = rep(c("Control", "Partial", "Exclosure"), 7),
  inter = c(
    .35, .33, .58),
  txt = c(
    "a","a","b")
)

ggplot(par |> 
         mutate(Herbivory = fct_relevel(Herbivory, "Control", "Partial")),
       aes(Herbivory, 1-inter, color = Herbivory)) +
  stat_summary()+
  geom_text(
    data = tukey_par |>
      mutate(Herbivory = fct_relevel(Herbivory, "Control", "Partial")),
    aes(Herbivory, inter, label = txt, color = Herbivory),
    inherit.aes = F,
    show.legend = FALSE,
    size = 7
  ) +
  theme_cowplot() +
  scale_color_brewer(
    palette = "Set1",
    breaks = c("Control", "Partial", "Exclosure"),
    labels = c("Large + small herbivores", "Large excluded", "All excluded")
  ) +
  labs(x = "", y = "PAR interception", color = "") +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "white", color = "black"),
    legend.justification = "center",
    text = element_text(size = 18),
    axis.text = element_text(size = 14))



# 5. TEMPORAL HETEROGENEITY ####

beta_temp <- read.csv("beta_temp.csv") |> 
  mutate(yr = as.character(yr))

# Calculate interannual compositional turnover per permanent plot
distancias_globales <- vegdist((beta_temp %>% select(-uidplot, -yr))^(1/4), method = "bray")
dispersion_temporal <- betadisper(distancias_globales, group = beta_temp$uidplot)

rango <- tibble(
  uidplot = dispersion_temporal$group.distances |> names(),
  distance = dispersion_temporal$group.distances |> as.numeric()
) |> 
  mutate(Herbivory = substr(uidplot, 1,(nchar(uidplot)-1)),
         Plot = substr(uidplot, nchar(uidplot), nchar(uidplot)),
         Fenced_Area = ifelse(Herbivory == "Control", 
                              paste0("fence_n_", Plot), 
                              paste0("fence_y_", Plot)))

# Model selection. Note 'dispformula' is used to explicitly model variance heteroscedasticity 
modelo0 <- glmmTMB(
  distance ~ Herbivory + (1 | Plot) + (1|Fenced_Area),
  dispformula = ~Herbivory,
  data = rango, REML = T
)
modelo1 <- update(modelo0, . ~ . - (1 | Fenced_Area))
anova(modelo0, modelo1)
modelo2 <- update(modelo1, . ~ . - (1 | Plot))
anova(modelo2, modelo1)

# Testing if modeling dispersion significantly improves fit vs homoscedastic assumption
modelo3 <- update(modelo2, dispformula = ~1)
anova(modelo2, modelo3)

check_model(modelo2)
Anova(modelo2)
emmeans(modelo2, pairwise ~ Herbivory)


# Plot: Figure S3 (Temporal Heterogeneity)

tukey_tmp <- data.frame(
  Herbivory = rep(c("Control", "Partial", "Exclosure"), 7),
  distance = c(
    .155, .185, .225),
  txt = c(
    "a","b","b")
)

ggplot(rango |> 
         mutate(Herbivory = fct_relevel(Herbivory, "Control", "Partial")),
       aes(Herbivory, distance, color = Herbivory)) +
  stat_summary()+
  geom_text(
    data = tukey_tmp |>
      mutate(Herbivory = fct_relevel(Herbivory, "Control", "Partial")),
    aes(Herbivory, distance, label = txt, color = Herbivory),
    inherit.aes = F,
    show.legend = FALSE,
    size = 7
  ) +
  theme_cowplot() +
  scale_color_brewer(
    palette = "Set1",
    breaks = c("Control", "Partial", "Exclosure"),
    labels = c("Large + small herbivores", "Large excluded", "All excluded")
  ) +
  labs(x = "", y = "Temporal heterogeneity in species composition", color = "") +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "white", color = "black"),
    legend.justification = "center",
    text = element_text(size = 18),
    axis.text = element_text(size = 14))



# 6. GRAZING INTENSITY EFFECTS ON KEY TAXA ####

intensity <- read.csv("intensity.csv")


# Plot: Figure S6 (Difference in relative cover)

ggplot(intensity, aes(x = Herbivory, y = cover, color = taxa))+
  geom_rect(xmin= 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf, inherit.aes = F, fill = "gray80")+
  geom_hline(yintercept = 0, linetype =2)+
  stat_summary(position = position_dodge(0.25), size = 1.2)+
  theme_cowplot()+
  scale_x_discrete(limits = c("fa_so", "pa", "ex"), 
                   name = "Comparison",
                   labels= c("Non-territorial vs.\nTerritorial",
                             "Large excluded vs.\nLarge + small herbivores",
                             "All excluded vs.\nLarge + small herbivores"))+
  labs( y = "Difference in relative cover", color = "Species")+
  theme(text = element_text(size = 18), 
        axis.text = element_text(size = 14),
        legend.text = element_text(face = "italic"))



# Plot: Figure S4 (Community composition temporal trends by species)

composition <- read.csv("composition.csv")

ggplot(composition |> 
         mutate(Herbivory = fct_relevel(Herbivory, "Control", "Partial")),
       aes(as.numeric(yr), cover, color = taxa)
) +
  facet_grid(
    . ~ Herbivory,
    labeller = labeller(
      .cols = c("Control" = "Large + small herbivores", "Partial" = "Large excluded", "Exclosure" = "All excluded")
    )
  ) +
  stat_summary() +
  stat_summary(fun = "mean", geom = "line") +
  theme_cowplot() +
  scale_color_paletteer_d(name = "Especie", "ggthemes::colorblind") +
  labs(x = "", y = "Cover (%)", color = "Species") +
  panel_border() +
  background_grid(major = 'y', minor = "none") +
  theme(
    legend.text = element_text(face = "italic"),
    legend.position = "top",
    strip.background = element_rect(fill = "white", color = "black"),
    legend.justification = "center",
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = "lightgray")
  )



# 7. MAGNITUDE OF EFFECT (ENVIRONMENTAL MODULATORS) ####

mods <- read.csv("mods.csv") |> 
  mutate(yr = as.character(yr),
         Plot = as.character(Plot))

# Standardize (Z-transform) predictors for direct effect size comparison
datos_modelos <- mods |> 
  mutate(
    densidad_z = scale(densidad),
    ie_anual_z = scale(1-IE.annual)
  )


## 7.1 Cover LRR: Large excluded ####

mod_cov_par <- glmmTMB(lrr_cov ~ densidad_z * ie_anual_z + (1|Plot), 
                       data = datos_modelos |> filter(Tratamiento == "Partial"), REML = T)
mod_cov_par1 <- update(mod_cov_par, . ~ . - (1 | Plot))
anova(mod_cov_par, mod_cov_par1)

check_model(mod_cov_par)

# Re-fit to ML to test fixed effects structure
mod_cov_par2 <- update(mod_cov_par, REML = F)
mod_cov_par3 <-update(mod_cov_par2, .~. - densidad_z:ie_anual_z)
anova(mod_cov_par2, mod_cov_par3) # Testing interaction

# Final evaluations reverted to REML = T for accurate parameter/variance estimates
Anova(mod_cov_par3 |> update(REML = T))
check_model(mod_cov_par3 |> update(REML = T))
summary(mod_cov_par3 |> update(REML = T))


## 7.2 Cover LRR: All excluded ####

mod_cov_exc <- glmmTMB(lrr_cov ~ densidad_z * ie_anual_z + (1|Plot), 
                       data = filter(datos_modelos, Tratamiento == "Exclosure"))
mod_cov_exc1 <- update(mod_cov_exc, . ~ . - (1 | Plot))
anova(mod_cov_exc, mod_cov_exc1)

mod_cov_exc_2<-update(mod_cov_exc, REML = F)
mod_cov_exc_3<-update(mod_cov_exc_2, . ~ . - densidad_z:ie_anual_z)
anova(mod_cov_exc_2, mod_cov_exc_3)

check_model(mod_cov_exc_3 |> update(REML = T))
Anova(mod_cov_exc_3 |> update(REML = T))
summary(mod_cov_exc_3 |> update(REML = T))


# 7.3 Richness LRR: Large excluded ####

mod_riq_par <- glmmTMB(lrr_riq ~ densidad_z * ie_anual_z + (1|Plot), 
                       data = datos_modelos |> filter(Tratamiento == "Partial"), REML = T)
mod_riq_par1 <- update(mod_riq_par, . ~ . - (1 | Plot))
anova(mod_riq_par, mod_riq_par1)

mod_riq_par2<-update(mod_riq_par1, REML = F)
mod_riq_par3 <-update(mod_riq_par2, .~. - densidad_z:ie_anual_z)
anova(mod_riq_par2, mod_riq_par3)

Anova(mod_riq_par3 |> update(REML = T))
check_model(mod_riq_par3 |> update(REML = T))
summary(mod_riq_par3 |> update(REML = T))

## 7.4 Richness LRR: All excluded ####

mod_riq_exc <- glmmTMB(lrr_riq ~ densidad_z * ie_anual_z + (1|Plot), 
                       data = filter(datos_modelos, Tratamiento == "Exclosure"))
mod_riq_exc1 <- update(mod_riq_exc, . ~ . - (1 | Plot))
anova(mod_riq_exc, mod_riq_exc1)

mod_riq_exc2<-update(mod_riq_exc1, REML = F)
mod_riq_exc3 <-update(mod_riq_exc2, .~. - densidad_z:ie_anual_z)
anova(mod_riq_exc2, mod_riq_exc3)

Anova(mod_riq_exc3 |> update(REML = T))
check_model(mod_riq_exc3 |> update(REML = T))
summary(mod_riq_exc3 |> update(REML = T))


# Plot: Figure 4 (Magnitude of effect predictions)

# Generate marginal predictions holding other variables constant using ggpredict
pred_riq_par <- ggpredict(mod_riq_par3 |> update(REML = T), terms = "densidad_z") %>% 
  as.data.frame() %>% mutate(Tratamiento = "Partial")

pred_riq_tot <- ggpredict(mod_riq_exc3 |> update(REML = T), terms = "densidad_z") %>% 
  as.data.frame() %>% mutate(Tratamiento = "Exclosure")

plot_data_riq <- bind_rows(pred_riq_par, pred_riq_tot) %>%
  mutate(Tratamiento = factor(Tratamiento, levels = c("Partial", "Exclosure")))

pred_riq_par_w <- ggpredict(mod_riq_par3 |> update(REML = T), terms = "ie_anual_z") %>% 
  as.data.frame() %>% mutate(Tratamiento = "Partial")

pred_riq_tot_w <- ggpredict(mod_riq_exc3 |> update(REML = T), terms = "ie_anual_z") %>% 
  as.data.frame() %>% mutate(Tratamiento = "Exclosure")

plot_data_riq_w <- bind_rows(pred_riq_par_w, pred_riq_tot_w) %>%
  mutate(Tratamiento = factor(Tratamiento, levels = c("Partial", "Exclosure")))

pred_cov_par <- ggpredict(mod_cov_par3 |> update(REML = T), terms = "ie_anual_z") %>% 
  as.data.frame() %>% mutate(Tratamiento = "Partial")

pred_cov_tot <- ggpredict(mod_cov_exc_3 |> update(REML = T), terms = "ie_anual_z") %>% 
  as.data.frame() %>% mutate(Tratamiento = "Exclosure")

plot_data_cov <- bind_rows(pred_cov_par, pred_cov_tot) %>%
  mutate(Tratamiento = factor(Tratamiento, levels = c("Partial", "Exclosure")))

pred_cov_par_d <- ggpredict(mod_cov_par3 |> update(REML = T), terms = "densidad_z") %>% 
  as.data.frame() %>% mutate(Tratamiento = "Partial")

pred_cov_tot_d <- ggpredict(mod_cov_exc_3 |> update(REML = T), terms = "densidad_z") %>% 
  as.data.frame() %>% mutate(Tratamiento = "Exclosure")

plot_data_cov_d <- bind_rows(pred_cov_par_d, pred_cov_tot_d) %>%
  mutate(Tratamiento = factor(Tratamiento, levels = c("Partial", "Exclosure")))


# Fixed color palette definition to ensure consistency when subsetting data
colors <- setNames(RColorBrewer::brewer.pal(n = 3, name = "Set1")[2:3], c("Partial", "Exclosure"))

# Panel assembly (subsetting layers to plot only significant/marginal relationships)
panel_A <- ggplot(plot_data_riq, aes(x = x, y = predicted, color = Tratamiento, fill = Tratamiento, linetype = Tratamiento)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  scale_color_manual(values = colors, labels = c("Large excluded", "All excluded")) +
  scale_fill_manual(values = colors, labels = c("Large excluded", "All excluded")) +
  scale_linetype_manual(values = c("Partial" = "solid", "Exclosure" = "solid"), guide = "none") +
  theme_cowplot() +
  labs(
    x = bquote("Guanaco density (Z-score)"),
    y = "Log-Response Ratio (Richness)",
    color = "", fill = ""
  ) +
  theme(legend.position = "top")

panel_B <- ggplot(plot_data_riq_w |> 
                    filter(Tratamiento=="Exclosure"), aes(x = x, y = predicted, color = Tratamiento, fill = Tratamiento, linetype = Tratamiento)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  scale_color_manual(values = colors, labels = c("Large excluded", "All excluded")) +
  scale_fill_manual(values = colors, labels = c("Large excluded", "All excluded")) +
  scale_linetype_manual(values = c("Exclosure" = "dashed"), guide = "none") +
  theme_cowplot() +
  labs(
    x = "Water stress (Z-score)",
    y = "Log-Response Ratio (Richness)",
    color = "Treatment", fill = "Treatment"
  ) +
  theme(legend.position = "none")

panel_C <- ggplot(plot_data_cov_d, aes(x = x, y = predicted, color = Tratamiento, fill = Tratamiento, linetype = Tratamiento)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(data = subset(plot_data_cov_d, Tratamiento == "Partial"), size = 1.2) +
  geom_ribbon(data = subset(plot_data_cov_d, Tratamiento == "Partial"), aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  scale_color_manual(values = colors, labels = c("Large excluded", "All excluded")) +
  scale_fill_manual(values = colors, labels = c("Large excluded", "All excluded")) +
  scale_linetype_manual(values = c("Partial" = "dashed"), guide = "none") +
  theme_cowplot() +
  labs(
    x = "Guanaco density (Z-score)",
    y = "Log-Response Ratio (Plant Cover)"
  ) +
  theme(legend.position = "none")

panel_D <- ggplot(plot_data_cov, aes(x = x, y = predicted, color = Tratamiento, fill = Tratamiento, linetype = Tratamiento)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(data = subset(plot_data_cov, Tratamiento == "Partial"), size = 1.2) +
  geom_ribbon(data = subset(plot_data_cov, Tratamiento == "Partial"), aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  scale_color_manual(values = colors, labels = c("Large excluded", "All excluded")) +
  scale_fill_manual(values = colors, labels = c("Large excluded", "All excluded")) +
  scale_linetype_manual(values = c("Partial" = "solid"), guide = "none") +
  theme_cowplot() +
  labs(
    x = "Water stress (Z-score)",
    y = "Log-Response Ratio (Plant Cover)"
  ) +
  theme(legend.position = "none")

# Combine panels using cowplot's plot_grid 
plot_grid(
  get_plot_component(panel_A + theme(
    legend.justification = "center",
    legend.spacing = unit(1, "cm"),
    legend.text = element_text(margin = margin(r = 10))
  ), "guide-box-top"),
  
  plot_grid(
    panel_A + theme(legend.position = "none", axis.title.x = element_blank()), 
    panel_B + theme(axis.title.y = element_blank(), axis.title.x = element_blank()), 
    panel_C, 
    panel_D + theme(axis.title.y = element_blank()), 
    ncol = 2, 
    labels = "AUTO",
    label_x = 0.17,
    align = "vh",
    axis = "tblr"
  ),
  rel_heights = c(.05, .95), 
  ncol = 1
)


# Plot: Figure S7 (Temporal overlap of plant metrics and environmental drivers)

env_data <- read.csv("env_data.csv")

# Scaling multipliers to fit density and water stress on the secondary y-axis
mult_riq <- 20
mult_cov <- 100

p_riq <- ggplot() +
  stat_summary(
    data = riqueza_l |> filter(name == "riq") |> ungroup() |> mutate(Herbivory = fct_relevel(Herbivory, "Control", "Partial")),
    aes(as.numeric(yr), value, color = Herbivory, group = Herbivory),
    position = position_dodge(width = 0.6), size = 1.05, fun.data = mean_se, geom = "pointrange"
  ) +
  geom_line(data = env_data, aes(x = as.numeric(yr), y = dens_ha * mult_riq), color = "black", linetype = "dashed", linewidth = 1) +
  geom_line(data = env_data, aes(x = as.numeric(yr), y = (1-IE.annual) * mult_riq), color = "gray50", linetype = "dotted", linewidth = 1) +
  geom_text(
    data = tukey |> filter(name == "riq") |> ungroup() |> mutate(Herbivory = fct_relevel(Herbivory, "Control", "Partial")),
    aes(yr, value, label = txt, color = Herbivory),
    position = position_dodge(width = 0.6), size = 6, show.legend = FALSE
  ) +
  scale_color_brewer(
    palette = "Set1", breaks = c("Control", "Partial", "Exclosure"), labels = c("Large + small herbivores", "Large excluded", "All excluded")
  ) +
  scale_y_continuous(
    name = "Plant richness",
    sec.axis = sec_axis(~ . / mult_riq, name = "")
  ) +
  theme_cowplot() +
  labs(x = "", color = "") +
  theme(
    legend.position = "top", legend.justification = "center",
    legend.text = element_text(margin = margin(r = 25)),
    text = element_text(size = 18), axis.text = element_text(size = 14),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(), # Hidden X axis to merge with lower panel
    plot.margin = margin(0, 30, 10, 10)
  )

p_cov <- ggplot() +
  stat_summary(
    data = riqueza_l |> filter(name == "cov") |> ungroup() |> mutate(Herbivory = fct_relevel(Herbivory, "Control", "Partial")),
    aes(as.numeric(yr), value, color = Herbivory, group = Herbivory),
    position = position_dodge(width = 0.6), size = 1.05, fun.data = mean_se, geom = "pointrange"
  ) +
  geom_line(data = env_data, aes(x = as.numeric(yr), y = dens_ha * mult_cov), color = "black", linetype = "dashed", linewidth = 1) +
  geom_line(data = env_data, aes(x = as.numeric(yr), y = (1-IE.annual) * mult_cov), color = "gray50", linetype = "dotted", linewidth = 1) +
  geom_text(
    data = tukey |> filter(name == "cov") |> ungroup() |> mutate(Herbivory = fct_relevel(Herbivory, "Control", "Partial")),
    aes(yr, value, label = txt, color = Herbivory),
    position = position_dodge(width = 0.6), size = 6, show.legend = FALSE
  ) +
  scale_color_brewer(palette = "Set1", guide = "none") +
  scale_y_continuous(
    name = "Plant cover (%)",
    sec.axis = sec_axis(~ . / mult_cov, name = "")
  ) +
  theme_cowplot() +
  labs(x = "") +
  theme(
    text = element_text(size = 18), axis.text = element_text(size = 14),
    plot.margin = margin(0, 30, 10, 10)
  )

figura_unida <- plot_grid(p_riq, p_cov, ncol = 1, align = "v", rel_heights = c(1, 1), labels = "AUTO")

# Final drawing to add a shared label for the secondary y-axis
figura_final <- ggdraw(figura_unida) +
  draw_label(
    "Guanacos (ind/ha; dashed) & water stress (dotted)", 
    x = 0.98,          
    y = 0.5,           
    angle = -90,       
    size = 20          
  )
print(figura_final)



# 8. SOIL PROPERTIES ####

nut_ex <- read.csv("nut_ex.csv") |> 
  mutate(Plot = as.character(Plot),
         Fenced_Area = ifelse(Herbivory == "Control", 
                              paste0("fence_n_", Plot), 
                              paste0("fence_y_", Plot)))


# 8.1 Total Inorganic Nitrogen (TIN = NH4 + NO3) ####

mod <- glmmTMB(
  nh4 + no3 ~ Herbivory + (1 | Plot) + (1|Fenced_Area),
  dispformula = ~Herbivory,
  nut_ex,
  REML = T
)

mod0 <- update(mod, dispformula = ~1)
anova(mod, mod0)

mod1 <- update(mod0, . ~ . - (1 | Plot))
anova(mod0, mod1)

mod2 <- update(mod0, . ~ . - (1 | Fenced_Area))
anova(mod0, mod2)

mod1 <- update(mod2, . ~ . - (1 | Plot))
anova(mod2, mod1)

Anova(mod1)
check_model(mod1)
emmeans(mod1, pairwise ~ Herbivory)

tukey_ni <- data.frame(
  Herbivory = rep(c("Control", "Partial", "Exclosure"), 7),
  distance = c(
    2.6,
    1.65,
    1.6
  ),
  txt = c(
    "a",
    "b",
    "b"
  )
)

ggplot(
  nut_ex |>
    mutate(Herbivory = fct_relevel(Herbivory, "Control", "Partial")),
  aes(Herbivory, nh4 + no3, color = Herbivory)
) +
  stat_summary() +
  geom_text(
    data = tukey_ni |>
      mutate(Herbivory = fct_relevel(Herbivory, "Control", "Partial")),
    aes(Herbivory, distance, label = txt, color = Herbivory),
    inherit.aes = F,
    show.legend = FALSE,
    size = 7
  ) +
  theme_cowplot() +
  scale_color_brewer(
    palette = "Set1",
    breaks = c("Control", "Partial", "Exclosure"),
    labels = c("Large + small herbivores", "Large excluded", "All excluded")
  ) +
  labs(x = "", y = "Total inorganic\nnitrogen concentration", color = "") +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "white", color = "black"),
    legend.justification = "center",
    text = element_text(size = 18),
    axis.text = element_text(size = 14)
  ) -> ni


# 8.2 Organic matter ####

modelo0 <- glmmTMB(
  organic_matter ~ Herbivory + (1 | Plot) + (1|Fenced_Area),
  data = nut_ex,
  REML = T
)

modelo1 <- update(modelo0, . ~ . - (1 | Plot))
anova(modelo0, modelo1)
modelo2 <- update(modelo0, . ~ . - (1 | Fenced_Area))
anova(modelo0, modelo2)
modelo3 <- update(modelo1, . ~ . - (1 | Fenced_Area))
anova(modelo1, modelo3)

Anova(modelo3)
check_model(modelo3)
emmeans(modelo3, pairwise ~ Herbivory)

tukey_mo <- data.frame(
  Herbivory = rep(c("Control", "Partial", "Exclosure"), 7),
  distance = c(
    2.33,
    2.7,
    2.77
  ),
  txt = c(
    "a",
    "ab",
    "b"
  )
)

(ggplot(
  nut_ex |>
    mutate(Herbivory = fct_relevel(Herbivory, "Control", "Partial")),
  aes(Herbivory, organic_matter, color = Herbivory)
) +
    stat_summary() +
    geom_text(
      data = tukey_mo |>
        mutate(Herbivory = fct_relevel(Herbivory, "Control", "Partial")),
      aes(Herbivory, distance, label = txt, color = Herbivory),
      inherit.aes = F,
      show.legend = FALSE,
      size = 7
    ) +
    theme_cowplot() +
    scale_color_brewer(
      palette = "Set1",
      breaks = c("Control", "Partial", "Exclosure"),
      labels = c("Large + small herbivores", "Large excluded", "All excluded")
    ) +
    labs(x = "", y = "Organic matter content (%)", color = "") +
    theme(
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black"),
      legend.justification = "center",
      text = element_text(size = 18),
      axis.text = element_text(size = 14)
    ) -> om)


# Plot: Figure 5 (Multipane Soil properties)

plot_grid(
  get_plot_component(ni+
                       theme(
                         # legend.spacing = unit(1.5, "cm"),
                         legend.text = element_text(margin = margin(r = 25))), "guide-box-top"),
  plot_grid(
    ni + theme(legend.position = "none"),
    om + theme(legend.position = "none"),
    labels = "AUTO",
    ncol = 2,
    align = "vh", axis = "tblr",
    label_x = 0.20,  # Offset to adjust label position
    label_y = 0.95
  ),
  ncol = 1,
  rel_heights = c(.1, .9)
)



# 9. HARE ACTIVITY (Droppings Count) ####

hares <- read.csv("hares.csv") |> 
  mutate(Plot = as.character(Plot),
         yr = as.character(yr),
         Fenced_Area = ifelse(Herbivory == "Control", 
                              paste0("fence_n_", Plot), 
                              paste0("fence_y_", Plot)),
         uidplot = paste(Herbivory, Plot, sep = "_"))

# Poisson GLMM for count data including year interaction
mod01 <- glmmTMB(
  droppings ~ Herbivory * yr + (1 | Plot) + ar1(yr + 0 | uidplot) + (1 | Fenced_Area),
  REML = T,
  hares,
  dispformula = ~Herbivory,
  family = poisson()
)

mod00 <- update(mod01, dispformula = ~1)
anova(mod00, mod01)

mod000 <- update(mod00, . ~ . - ar1(fecha + 0 | uidplot))
anova(mod000, mod00)


mod0 <- update(mod000, . ~ . - (1 | Fenced_Area))
anova(mod0, mod000)
mod1 <- update(mod000, . ~ . - (1 | Plot))
anova(mod1, mod000)

mod2 <- update(mod0, . ~ . - (1 | Plot))
anova(mod2, mod0)


Anova(mod2)
check_model(mod2)
emmeans(mod2, pairwise ~ yr)

# Plot: Figure S8 (Hare activity)

ggplot(hares, aes(yr, droppings, color = Herbivory)) +
  stat_summary(size = 1.2) +
  theme_cowplot() +
  scale_color_brewer(
    palette = "Set1",
    breaks = c("Control", "Partial", "Exclosure"),
    labels = c("Large + small herbivores", "Large excluded", "All excluded")
  ) +
  theme(
    legend.position = "top",
    legend.justification = "center",
    text = element_text(size = 20),
    axis.text = element_text(size = 16)
  ) +
  labs(color = "", x = "", y = expression(Hare~activity~(droppings/50~m^2)))