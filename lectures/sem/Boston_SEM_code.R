# library(devtools)
# install_github("jslefche/piecewiseSEM")

library(piecewiseSEM)

# Read in data
keeley = read.csv("keeley.csv")

head(keeley)

### Create list of structured equations

modelList = list(
  
  lm(abiotic ~ distance, data = keeley),
  
  lm(hetero ~ distance, data = keeley),
  
  lm(rich ~ abiotic + hetero, data = keeley)
  
)

### Evaluate fit

sem.fit(modelList, data = keeley)

# Add significant path

modelList2 = list(
  
  lm(abiotic ~ distance, data = keeley),
  
  lm(hetero ~ distance, data = keeley),
  
  lm(rich ~ distance + abiotic + hetero, data = keeley)
  
)

sem.fit(modelList2, keeley)

### Evaluate individual model assumptions & fits

# Graphical evaluation
par(mfrow = c(2, 2))

lapply(modelList2, plot, which = 1)

sem.model.fits(modelList2)

### Get coefficients

sem.coefs(modelList2)

sem.coefs(modelList2, keeley, standardize = "scale")

### Compare to lavaan

library(lavaan)

modelList2.lavaan = sem.lavaan(modelList2, keeley)

summary(modelList2.lavaan, rsq = T, standardize = T)

### Correlated errors

modelList3 = list(
  
  lm(abiotic ~ distance, data = keeley),
  
  lm(hetero ~ distance, data = keeley),
  
  lm(rich ~ distance + hetero, data = keeley)
  
)

# Includes path in d-sep tests
sem.fit(modelList3, keeley)

sem.fit(modelList3, keeley, corr.errors = "rich ~~ abiotic")

# Conducts significance test
sem.coefs(modelList3, keeley, corr.errors = "rich ~~ abiotic")

sem.coefs(modelList3, keeley, standardize = "scale", corr.errors = "rich ~~ abiotic")

# Compare to lavaan
sem.lavaan(modelList3, keeley, corr.errors = "rich ~~ abiotic")

### Compare using AIC
modelList4 = list(
  
  lm(hetero ~ distance, data = keeley),
  
  lm(rich ~ distance + hetero, data = keeley)
  
)

sem.fit(modelList4, keeley, add.vars = "abiotic")

### Get partial correlation plot
dev.off()

partial.resid(rich ~ distance, modelList2, keeley)

####################################################################

#install.packages("gridExtra")
library(gridExtra)
# install.packages("nlme")
library(nlme)
# install.packages("lme4")
library(lme4)

library(piecewiseSEM)

### Shipley (2009)

# Load data
data(shipley2009)

# Create list of structured equations

modelList = list(
  
  lme(DD ~ lat, random = ~1|site/tree, na.action = na.omit, 
      data = shipley2009),
  
  lme(Date ~ DD, random = ~1|site/tree, na.action = na.omit, 
      data = shipley2009),
  
  lme(Growth ~ Date, random = ~1|site/tree, na.action = na.omit, 
      data = shipley2009),
  
  glmer(Live ~ Growth + (1|site) + (1|tree), 
        family=binomial(link = "logit"), data = shipley2009) 
  
)

# Evaluate fit of entire SEM
sem.fit(modelList, shipley2009)

# ...of individual models
sem.model.fits(modelList)

# Plot residuals vs fitted values
par(mfrow = c(2, 2))

do.call(grid.arrange, c(lapply(modelList, plot), nrow = 2))

# Get scaled coefficients

(coef.table = sem.coefs(modelList, shipley2009, standardize = "scale"))

# Get indirect effect of latitude on survival

prod(coef.table$estimate)

### Alternate model

# Generate list of structured equations 

modelList2 = list(
  
  lme(DD ~ lat, random = ~1|site/tree, na.action = na.omit, 
      data = shipley2009),
  
  lme(Date ~ DD, random = ~1|site/tree, na.action = na.omit, 
      data = shipley2009),
  
  lme(Growth ~ Date + DD, random = ~1|site/tree, na.action = na.omit, 
      data = shipley2009),
  
  glmer(Live ~ Growth + (1|site) + (1|tree), 
        family=binomial(link = "logit"), data = shipley2009) 
  
)

# Evaluate fit

sem.fit(modelList2, shipley2009)

### Phytoplankton example

# Load data
durocher = read.csv("durocher.csv")

head(durocher)

# Remove NAs
# durocher = durocher[complete.cases(durocher), ]

# Create model list
modelList3 = list(
  
  lme(CR ~ Std.Temp, random = ~ 1 | Pond.ID, na.action = na.omit, durocher),
  
  lme(Prich ~ Std.Temp, random = ~ 1 | Pond.ID, na.action = na.omit, durocher),
  
  lme(Pbio ~ Prich, random = ~ 1 | Pond.ID, na.action = na.omit, durocher),
  
  lme(GPP ~ Pbio, random = ~ 1 | Pond.ID, na.action = na.omit, durocher)
  
)

# Evaluate fit
sem.fit(modelList3, durocher)

# Get standardized coefficients
sem.coefs(modelList3, durocher, standardize = "scale")

# Get individual model fits
sem.model.fits(modelList3)
