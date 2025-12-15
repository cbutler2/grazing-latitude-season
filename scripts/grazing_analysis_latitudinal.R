# fit models to urchin grazing rates, kelp carbon and nitrogen, and temperature data across latitude
# Claire Butler
# April 2024



# required packages
library(pacman)
p_load(MASS, tidyverse, mgcv, lubridate)



# Latitudinal data ---------------------------------------------------------


## Fit gaulss model to grazing data
model_data <- read_csv('data/latitudinal_grazing.csv')

model_data <- model_data %>% 
  mutate(Date_Deploy = dmy(Date_Deploy),
         Date_Retrieve = dmy(Date_Retrieve)) %>% 
  mutate(across(where(is.character), as.factor),
         Cage_ID = factor(Cage_ID)) 


# centros
model_data_centro <- filter(model_data, Species == 'Centrostephanus rodgersii')

fit1 <- gam(list(resid ~ s(Gonad_Index) + s(Gut_Index) + 
                   s(Latitude, k=4) + s(Cage_ID, bs = 're'), ~s(Latitude, k = 4)), 
            family = gaulss(), 
            method = 'REML', 
            data = model_data_centro)

# helios
model_data_helio <- filter(model_data, Species == 'Heliocidaris erythrogramma')

fit2 <- gam(list(resid ~ s(Gonad_Index) + s(Gut_Index) + 
                   s(Latitude, k=5) + s(Cage_ID, bs = 're'),  ~s(Latitude, k=5)), 
            family = gaulss(), 
            method = 'REML', 
            data = model_data_helio)



## Fit models to kelp CN data

nut2 <- read_csv('data/latitudinal_cn.csv')

# calculate molecular ratio
nut2 <- nut2 %>% 
  mutate(carbon_mol = Per_Carbon_Factored/12.01,
         nitrogen_mol = Per_Nitrogen_Factored/14.01,
         ratio_mol = carbon_mol/nitrogen_mol)

# % nitrogen
fit.nit <- gam(Per_Nitrogen_Factored ~ s(Latitude, k=4), method = 'REML', data = nut2)

# mol ratio
fit.cn <- gam(log(ratio_mol) ~ s(Latitude, k=4), method = 'REML', data = nut2)


## Fit models to urchin abundance data

rls_urchins <- read_csv('data/rls_urchins.csv')

# filter to southeast Australia
rls_urchins <- rls_urchins %>% 
  filter(longitude >= 145.92 & longitude <= 153.63)

### Summarise to get average per 50m2
rls_abundance <- rls_urchins %>% 
  group_by(country, area, ecoregion, realm, location, site_code, site_name, longitude,latitude, species_name, survey_id, block) %>% 
  summarise(abundance=sum(total)) %>% # sum counts of each species within each transect (=survey_id) and block - this sums counts across size classes
  group_by(country, area, ecoregion, realm, location, site_code, site_name, longitude,latitude,species_name, survey_id) %>% 
  summarise(abundance = mean(abundance)) %>% # average counts per transect (=survey_id), because RLS has 2 blocks per transect, and ATRC has only 1
  group_by(country, area, ecoregion, realm, location, site_code, site_name, longitude,latitude,species_name) %>% # group by site
  summarise(N = length(abundance),
            mean = mean(abundance),
            sd = sd(abundance),
            se = sd / sqrt(N)) %>% # average abundance per site across time
  filter(country == 'Australia' & (species_name=='Centrostephanus rodgersii'| species_name=='Heliocidaris erythrogramma')) %>%  # note there is more rows than for unique sites because both species exist at some sites???
  droplevels()

rls_abundance <- rls_abundance %>% 
  ungroup() %>% 
  dplyr::select(site_name, latitude, longitude, species_name, N, mean, sd, se) %>% 
  rename(Site = site_name,
         Latitude = latitude,
         Species = species_name) 

# centro model
rls_centros <- rls_abundance %>% 
  filter(Species== 'Centrostephanus rodgersii')

fit.rlsc <- glm.nb(log(mean) ~ Latitude + I(Latitude^2), data = rls_centros)

# helio model
rls_helios <- rls_abundance %>% 
  filter(Species== 'Heliocidaris erythrogramma')

fit.rlsh <- glm.nb(log(mean) ~ Latitude, data = rls_helios)



## Temperature data

# logger data
temp2 <- read_csv('data/latitudinal_temp_logger.csv')

temp2_sum <- temp2 %>% 
  group_by(Site, Latitude) %>% 
  summarise(N = length(Temp),
            mean = mean(Temp),
            sd = sd(Temp),
            se = sd / sqrt(N))


# satellite data
sst <- read_csv('data/latitudinal_satellite_sst.csv')

sst <- sst %>% 
  mutate(Month = format(date, "%m"))

sst <- sst %>% 
  mutate(Site = factor(Site),
         Month = factor(Month)) %>% 
  filter(Site == 'Sawtell' & Month == '02' | 
           Site == 'Forster' & Month == '02' |
           Site == 'Shellharbour' & Month == '02' |
           Site == 'Merimbula' & Month == '03' |
           Site == 'Fortescue' & Month == '03')

sst <- sst %>% 
  mutate(Year = format(date, "%Y"))

sst_sum <- sst %>%
  group_by(Site, Latitude, Longitude, Year, Month) %>%
  summarise(N = length(sst),
            mean = mean(sst),
            sd = sd(sst),
            se = sd / sqrt(N))


fit.sst <- gam(mean ~ s(Latitude, k=5), method = 'REML', data = sst_sum)


# Grazing pressure --------------------------------------------------------

## Centros

# create test data
testdata_centro = data.frame(Latitude=rep(seq(min(model_data_centro$Latitude),max(model_data_centro$Latitude) + 1, 0.5), 1),
                             Species = rep("Centrostephanus rodgersii", each = 28),
                             # Whole_Weight=mean(fit1$model$Whole_Weight),
                             Gonad_Index=mean(fit1$model$Gonad_Index),
                             Gut_Index=mean(fit1$model$Gut_Index),
                             Cage_ID = 4)

# perform bootstarp
n_boot <- 100 # define the number of bootstrap samples
cboot_gr <- model_data_centro # store original data
cboot_ab <- rls_centros # store original data
bootstrap_results_centro <- matrix(NA, nrow = n_boot, ncol = 28) # create a matrix to store the bootstrap results

set.seed(123)  # for reproducibility

for(i in 1:n_boot){
  
  # take random sample of rows from raw grazing data, with replacement
  igrazing <- sample(1:nrow(model_data_centro), nrow(model_data_centro),
                     replace = TRUE)
  
  # create new data set with these rows
  cboot_gr <- model_data_centro[igrazing,] 
  
  # Fit the model to the resampled data
  fit1_boot <- update(fit1, data = cboot_gr)
  
  # Do the same for abundance
  iabund <- sample(1:nrow(rls_centros), nrow(rls_centros),
                   replace = TRUE)
  
  cboot_ab <- rls_centros[iabund,]
  
  fit.rlsc_boot <- update(fit.rlsc, data = cboot_ab)
  
  # Predict and calculate total consumption
  predicted_grazing_boot <- predict(fit1_boot, newdata = testdata_centro)
  predicted_abundance_boot <- predict(fit.rlsc_boot, newdata = testdata_centro, type = 'response')
  total_consumption_boot <- (exp(predicted_abundance_boot)*200) * (predicted_grazing_boot/1000) # conversion from g/50m2 to kg/ha
  
  # Store the results
  bootstrap_results_centro[i, ] <- total_consumption_boot[, 1]
}


testdata_centro$mean <- apply(bootstrap_results_centro, 2, mean)
testdata_centro$sd <- apply(bootstrap_results_centro, 2, sd)
testdata_centro$se <- testdata_centro$sd/sqrt(50)

gp_centros <- testdata_centro %>% 
  select(Latitude, Species, mean, sd, se)


## Helios
# create test data 
testdata_helio = data.frame(Latitude=rep(seq(min(model_data_helio$Latitude),max(model_data_helio$Latitude) + 1, 0.5),1),
                            Species = rep("Heliocidaris erythrogramma", each = 28),
                            # Whole_Weight=mean(fit2$model$Whole_Weight),
                            Gonad_Index=mean(fit2$model$Gonad_Index),
                            Gut_Index=mean(fit2$model$Gut_Index),
                            Cage_ID = 4)



# Perform the bootstrap

n_boot <- 100 # define the number of bootstrap samples
hboot_gr <- model_data_helio # store original data
hboot_ab <- rls_helios # store original data
bootstrap_results_helio <- matrix(NA, nrow = n_boot, ncol = 28) # create a matrix to store the bootstrap results


set.seed(123)  # for reproducibility

for(i in 1:n_boot){
  
  # take random sample of rows from raw grazing data, with replacement
  igrazing <- sample(1:nrow(model_data_helio), nrow(model_data_helio),
                     replace = TRUE)
  
  hboot_gr <- model_data_helio[igrazing,] # create new data set with these rows
  
  # Fit the model to the resampled data
  fit2_boot <- update(fit2, data = hboot_gr)
  
  # Do the same for abundance
  iabund <- sample(1:nrow(rls_helios), nrow(rls_helios),
                   replace = TRUE)
  
  hboot_ab <- rls_helios[iabund,]
  
  fit.rlsh_boot <- update(fit.rlsh, data = hboot_ab)
  
  # Predict and calculate total consumption
  predicted_grazing_boot <- predict(fit2_boot, newdata = testdata_helio)
  predicted_abundance_boot <- predict(fit.rlsh_boot, newdata = testdata_helio, type = 'response')
  total_consumption_boot <- (exp(predicted_abundance_boot)*200) * (predicted_grazing_boot/1000)
  
  # Store the results
  bootstrap_results_helio[i, ] <- total_consumption_boot[, 1]
}


testdata_helio$mean <- apply(bootstrap_results_helio, 2, mean)
testdata_helio$sd <- apply(bootstrap_results_helio, 2, sd)
testdata_helio$se <- testdata_helio$sd/sqrt(50)

gp_helios <- testdata_helio %>% 
  select(Latitude, Species, mean, sd, se)




