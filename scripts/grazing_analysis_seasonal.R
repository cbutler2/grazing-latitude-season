# fit models to urchin grazing rates, kelp carbon and nitrogen, and temperature data across seasonal gradient
# Claire Butler
# April 2024




# required packages
library(pacman)
p_load(MASS, tidyverse, mgcv, lubridate, geosphere)



# Seasonal data -----------------------------------------------------------

# Fit models to grazing data
model_data <- read_csv('data/seasonal_grazing.csv')

model_data <- model_data %>% 
  mutate(across(where(is.character), as.factor),
         Cage_ID = factor(Cage_ID)) 


# centro
model_data_centro <- model_data %>% 
  filter(Species == 'Centrostephanus rodgersii')

fit1 <- gam(list(resid ~ s(Gut_Index) + 
                   s(Gonad_Index) + s(daylength) + s(Mean_Logger_Temp) +
                   s(Cage_ID, bs = 're') + Duration_Exposure, ~s(dayofyear)), 
            family = gaulss(), 
            method = 'REML',
            data = model_data_centro)


# helio 
model_data_helio <- model_data %>% 
  filter(Species == 'Heliocidaris erythrogramma')

fit2 <- gam(list(resid ~ s(Gut_Index) + 
                   s(Mean_Logger_Temp) + te(Gonad_Index, Mean_Nitrogen) + 
                   s(Cage_ID, bs = 're') + Duration_Exposure, ~s(dayofyear)), 
            family = gaulss(), 
            method = 'REML',
            data = model_data_helio)


## Fit models to kelp CN data

nut1 <- read_csv('data/seasonal_cn.csv')

nut1 <- nut1 %>%
  mutate(Date = dmy(Date),
         time = as.numeric(Date)/1000)

# calculate mol ratio
nut1 <- nut1 %>%
  mutate(carbon_mol = Per_Carbon_Factored/12.01,
         nitrogen_mol = Per_Nitrogen_Factored/14.01,
         ratio_mol = carbon_mol/nitrogen_mol)

# % nitrogen
fit.nit <- gam(Per_Nitrogen_Factored ~ s(time, k = 6), method = 'REML', data = nut1)

# mol ratio
fit.cn <- gam(ratio_mol ~ s(time, k = 6), method = 'REML', data = nut1)

## Temperature data

# logger data
temp1 <- read_csv('data/seasonal_temp_logger.csv')

temp1 <- temp1 %>%
  mutate(time = as.numeric(datetime)/1000)

temp1_sum <- temp1 %>%
  group_by(month, year) %>%
  summarise(N    = length(temp),
            mean = mean(temp),
            sd   = sd(temp),
            se   = sd / sqrt(N)) %>%
  unite('date', year:month, sep = '-') %>%
  mutate(date = ym(date))



# satellite data

sst <- read_csv('data/latitudinal_satellite_sst.csv')

sst <- sst %>%
  mutate(Site = factor(Site)) %>%
  filter(Site == 'Fortescue') %>%
  filter(date > ymd('2021-09-01') & date < ymd('2023-03-30'))

## Gut and gonad index
gonad_indices <- model_data %>%
  group_by(Date_Deploy, Species) %>%
  summarise(N    = length(Gonad_Index),
            mean = mean(Gonad_Index),
            sd   = sd(Gonad_Index),
            se   = sd / sqrt(N)) %>%
  mutate(Index = 'Gonad Index')

gut_indices <- model_data %>%
  group_by(Date_Deploy, Species) %>%
  summarise(N    = length(Gut_Index),
            mean = mean(Gut_Index),
            sd   = sd(Gut_Index),
            se   = sd / sqrt(N)) %>%
  mutate(Index = 'Gut Index')

indices <- rbind(gonad_indices, gut_indices)

indices_centro <- indices %>%
  filter(Species == 'Centrostephanus rodgersii')

indices_helio <- indices %>%
  filter(Species == 'Heliocidaris erythrogramma')
