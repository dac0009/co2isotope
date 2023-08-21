library(isocalcR)
library(tidyr)
library(dplyr)
library(cowplot)
library(ggplot2)
library(ggthemes)
library(gt)
library(Metrics)
library(car)

data("CO2data")
data("piru13C")
summary(CO2data)
summary(piru13C)

combined_data <- merge(piru13C, CO2data, by.x = "Year", by.y = "yr", all = TRUE)
combined_data$d13C_wood <- combined_data$wood.d13C / 1000

combined_data$iWUE_simple <- mapply(d13C.to.iWUE,
                                    d13C.plant = combined_data$wood.d13C,
                                    year = combined_data$Year,
                                    elevation = combined_data$Elevation_m,
                                    temp = combined_data$MGT_C,
                                    method = "simple",
                                    tissue = "wood")

combined_data$iWUE_photorespiration <- mapply(d13C.to.iWUE,
                                              d13C.plant = combined_data$wood.d13C,
                                              year = combined_data$Year,
                                              elevation = combined_data$Elevation_m,
                                              temp = combined_data$MGT_C,
                                              method = "photorespiration",
                                              frac = combined_data$frac)

combined_data$iWUE_mesophyll <- mapply(d13C.to.iWUE,
                                       d13C.plant = combined_data$wood.d13C,
                                       year = combined_data$Year,
                                       elevation = combined_data$Elevation_m,
                                       temp = combined_data$MGT_C,
                                       method = "mesophyll",
                                       frac = combined_data$frac)

data_long <- combined_data %>%
  select(Year, d13C_wood, iWUE_simple, iWUE_photorespiration, iWUE_mesophyll, Ca, d13C.atm) %>%
  pivot_longer(cols = c(d13C_wood, iWUE_simple, iWUE_photorespiration, iWUE_mesophyll),
               names_to = "Variable",
               values_to = "Value")

data_subset <- combined_data %>% 
  select(Year, wood.d13C, MGT_C, Elevation_m, Ca, iWUE_mesophyll, iWUE_photorespiration, iWUE_simple) %>%
  na.omit()


m1 <- lm(d13C.atm ~ Ca, data = CO2data) # Atmospheric CO2 concentration (ppm) & 13C sequestration
m1sum <- summary(m1)  

m2 <- lm(wood.d13C ~ iWUE_mesophyll + iWUE_photorespiration + iWUE_simple, data = data_subset) # 13C concentration & Intrinsic Water Use Effeciency
m2sum <- summary(m2)

m3 <- lm(Elevation_m ~ wood.d13C + iWUE_mesophyll + iWUE_photorespiration + iWUE_simple, data = data_subset) # Elevation & 13C concentration + iWUE fractions
m3sum <- summary(m3)

m4 <- lm(MGT_C ~ wood.d13C + iWUE_mesophyll + iWUE_photorespiration + iWUE_simple, data = data_subset) # Avg growth temperature & 13C concentration + iWUE fractions
m4sum <- summary(m4)


corCO2data <- cor(CO2data[,2:3], method = "pearson")
cor_piruCO2 <- cor(data_subset[,2:8], method = "pearson")



outm1 <- outlierTest(m1)
outm2 <- outlierTest(m2)
outm3 <- outlierTest(m3)
outm4 <- outlierTest(m4)



fit1 <- aov(m1)
fit2 <- aov(m2)
fit3 <- aov(m3)
fit4 <- aov(m4)


ncv_m1 <- ncvTest(m1)
ncv_m2 <- ncvTest(m2)
ncv_m3 <- ncvTest(m3)
ncv_m4 <- ncvTest(m4)

kru_m1 <- kruskal.test(d13C.atm ~ Ca, data = CO2data)
kru_m2 <- kruskal.test(wood.d13C ~ iWUE_simple, data = data_subset)
kru_m3 <- kruskal.test(Elevation_m ~ iWUE_simple, data = data_subset)
kru_m4 <- kruskal.test(MGT_C ~ iWUE_simple, data = data_subset)


ggplot(CO2data, aes(x = yr, y = Ca)) +
  geom_line(color = "#B40F20", lwd = 1.5) +
  theme_cowplot() +
  labs(x = "Year", y = "Atmospheric CO2 (ppm)") +
  theme(axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.position = "none") +
  scale_x_continuous(breaks = c(1940, 1960, 1980, 2000, 2020),
                     labels = c("1940", "1960", "1980", "2000", "2020"),
                     expand = c(0, 0),
                     limits = c(1940, 2020)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1))


ggplot(CO2data, aes(x = yr, y = d13C.atm)) +
  geom_line(color = "#003300", size = 1.5) +
  theme_cowplot() +
  labs(x = "Year", y = "Carbon Isotope Composition") +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.position = "none") +
  scale_x_continuous(breaks = c(1940, 1960, 1980, 2000, 2020),
                     labels = c("1940", "1960", "1980", "2000", "2020"),
                     expand = c(0, 0),
                     limits = c(1940, 2020)) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1))



ggplot(piru13C, aes(x = Year, y = wood.d13C)) +
  geom_line(color = "#003300") +
  theme_cowplot() +
  labs(x = "Year", y = "δ13CO2") +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
  scale_x_continuous(breaks = c(1940, 1960, 1980, 2000, 2020),
                     labels = c("1940", "1960", "1980", "2000", "2020"),
                     expand = c(0, 0),
                     limits = c(1940, 2020))



data_subset %>%
  gt() %>%
  tab_header(title = "Water Use Effeciency of δ13CO2 Plants") %>%
  cols_label(Year = "Year", wood.d13C = "δ13CO2 Wood", MGT_C = "Temperature (C)", Elevation_m = "Elevation (m)", Ca = "Atmosphere CO2 (ppm)", 
             iWUE_mesophyll = "WUE Mesophyll", iWUE_photorespiration = "WUE Photorespiration", iWUE_simple = "WUE Simple (No Frac)") %>%
  tab_spanner(label = "Calculated WUE Fraction", columns = c("iWUE_mesophyll", "iWUE_photorespiration", "iWUE_simple")) %>%
  tab_source_note(source_note = "Data Altered From Mathias & Thomas 2018; Belmecheri & Lavergne 2020") %>%
  tab_style(style = list(cell_text(align="center", weight="bold")),locations=cells_column_labels()) %>%
  tab_style(style = list(cell_text(align="center", weight="bold", color="#02401B")),locations = cells_title())%>%
  tab_style(style = list(cell_text(align="center", style="italic", color="#02401B")),locations = cells_column_spanners()) %>%
  tab_style(style = list(cell_text(align="center")), locations = cells_body()) %>%
  tab_options(table_body.hlines.color="#02401B", table.border.top.color="#0C1707", table.border.bottom.color = "#02401B",
              table_body.border.bottom.color =  "#0C1707",
              heading.border.bottom.color = "#0C1707",
              column_labels.border.top.color = "#02401B",
              column_labels.border.bottom.color = "#02401B",)



crPlots(m1)
crPlots(m3)
crPlots(m4)


ggplot(data_long, aes(x = Year, y = Value, color = Variable)) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(group = Variable), color = "#000000") +
  theme_bw() +
  labs(x = "Year", y = "iWUE") +
  scale_x_continuous(breaks = seq(1940, 2020, 10), limits = c(1940, 2020)) +
  scale_color_tableau() +
  facet_wrap(~Variable, scales = "free_y") +
  theme(legend.position = "bottom")



print(m1sum)
print(m2sum)
print(m3sum)
print(m4sum)



print(corCO2data)
print(cor_piruCO2)


print(outm1)
print(outm2)
print(outm3)
print(outm4)




print(fit1)
print(fit2)
print(fit3)
print(fit4)

print(ncv_m1)
print(ncv_m2)
print(ncv_m3)
print(ncv_m4)


print(kru_m1)
print(kru_m2)
print(kru_m3)
print(kru_m4)








