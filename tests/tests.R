library(virtualspecies)
a <- matrix(rep(dnorm(1:100, 50, sd = 25)),
            nrow = 100, ncol = 100, byrow = TRUE)
env <- c(rast(a * dnorm(1:100, 50, sd = 25)),
         rast(a * 1:100))
names(env) <- c("variable1", "variable2")
plot(env) # Illustration of the variables

# Easy creation of the parameter list:
# see in real time the shape of the response functions
parameters <- formatFunctions(variable1 = c(fun = 'dnorm', mean = 1e-04,
                                             sd = 1e-04),
                              variable2 = c(fun = 'linearFun', a = 1, b = 0))

plotResponse(x = env, parameters = parameters, 
             approach = "response")

sp1 <- generateSpFromFun(env, parameters, plot = TRUE)

# If you provide env, then you can see the shape of response functions:
parameters <- formatFunctions(x = env,
                              variable1 = c(fun = 'dnorm', mean = 1e-04,
                                             sd = 1e-04),
                              variable2 = c(fun = 'linearFun', a = 1, b = 0))

sp1 <- generateSpFromFun(env, parameters, plot = TRUE)
sp1 <- generateSpFromFun(env, parameters, plot = TRUE,
                         formula = NULL,
                         species.type = "additive")
sp1 <- generateSpFromFun(env, parameters, plot = TRUE,
                         formula = NULL,
                         species.type = "multiplicative")
sp1 <- generateSpFromFun(env, parameters, plot = TRUE,
                         formula = NULL,
                         species.type = "multiplicative",
                         rescale = FALSE)
sp1 <- generateSpFromFun(env, parameters, plot = TRUE,
                         formula = NULL,
                         species.type = "multiplicative",
                         rescale = FALSE,
                         rescale.each.response = FALSE)
sp1 <- generateSpFromFun(env, parameters, plot = TRUE,
                         formula = "variable1 + variable2",
                         species.type = "multiplicative",
                         rescale = FALSE,
                         rescale.each.response = TRUE)
sp1 <- generateSpFromFun(env, parameters, plot = TRUE,
                         formula = "sqrt(variable1) + variable2 + 2 * variable2^2 + variable2^3",
                         species.type = "multiplicative",
                         rescale = FALSE,
                         rescale.each.response = TRUE)
sp1 <- generateSpFromFun(env, parameters, plot = TRUE,
                         formula = "sqrt(variable1) + variable2 + 2 * variable2^2 + variable2^3",
                         species.type = "multiplicative",
                         rescale = TRUE,
                         rescale.each.response = TRUE)


sp2 <- convertToPA(sp1,
                   PA.method = "threshold")
sp2 <- convertToPA(sp1,
                   PA.method = "threshold",
                   beta = 0.5)
sp2 <- convertToPA(sp1,
                   PA.method = "threshold",
                   beta = 0)
sp2 <- convertToPA(sp1,
                   PA.method = "threshold",
                   beta = 1)
sp2 <- convertToPA(sp1,
                   PA.method = "threshold",
                   beta = "random")
sp2 <- convertToPA(sp1,
                   PA.method = "threshold",
                   species.prevalence = .2)
sp2 <- convertToPA(sp1,
                   PA.method = "threshold",
                   species.prevalence = .02)
sp2 <- convertToPA(sp1,
                   PA.method = "threshold",
                   species.prevalence = .9)
sp2 <- convertToPA(sp1,
                   PA.method = "threshold",
                   species.prevalence = .99)
sp2 <- convertToPA(sp1,
                   PA.method = "probability",
                   prob.method = "logistic",
                   beta = "random",
                   a = NULL,
                   b = NULL,
                   species.prevalence = .99)
sp2 <- convertToPA(sp1,
                   PA.method = "probability",
                   prob.method = "logistic",
                   beta = "random",
                   a = NULL,
                   b = NULL,
                   species.prevalence = .5)
sp2 <- convertToPA(sp1,
                   PA.method = "probability",
                   prob.method = "logistic",
                   beta = .5,
                   a = NULL,
                   b = NULL)


sp3 <- generateRandomSp(env)

a <- matrix(rep(dnorm(1:100, 50, sd = 25)),
            nrow = 100, ncol = 100, byrow = TRUE)

env1 <- c(rast(a * dnorm(1:100, 50, sd = 25)),
          rast(a * 1:100),
          rast(a),
          rast(t(a)))
names(env1) <- c("var1", "var2", "var3", "var4")
b <- matrix(rep(dnorm(1:100, 25, sd = 50)),
            nrow = 100, ncol = 100, byrow = TRUE)

env2 <- c(rast(b * dnorm(1:100, 50, sd = 25)),
          rast(b * 1:100),
          rast(b),
          rast(t(b)))

names(env2) <- c("var1", "var2", "var3", "var4")

# Generating a species with the BCA
sp4 <- generateSpFromBCA(raster.stack.current = env1, raster.stack.future = env2)
plotResponse(sp4)

sp4 <- convertToPA(sp4)
plotSuitabilityToProba(sp4)

sp4 <- convertToPA(sp4,
                   PA.method = "threshold")
plotSuitabilityToProba(sp4)
sp4 <- convertToPA(sp4,
                   PA.method = "threshold",
                   beta = 0.5)
plotSuitabilityToProba(sp4)
sp4 <- convertToPA(sp4,
                   PA.method = "threshold",
                   beta = 0)
plotSuitabilityToProba(sp4)
sp4 <- convertToPA(sp4,
                   PA.method = "threshold",
                   beta = 1)
plotSuitabilityToProba(sp4)
sp4 <- convertToPA(sp4,
                   PA.method = "threshold",
                   beta = "random")
plotSuitabilityToProba(sp4)
sp4 <- convertToPA(sp4,
                   PA.method = "threshold",
                   species.prevalence = .2)
plotSuitabilityToProba(sp4)
sp4 <- convertToPA(sp4,
                   PA.method = "threshold",
                   species.prevalence = .02)
plotSuitabilityToProba(sp4)
sp4 <- convertToPA(sp4,
                   PA.method = "threshold",
                   species.prevalence = .9)
plotSuitabilityToProba(sp4)
sp4 <- convertToPA(sp4,
                   PA.method = "threshold",
                   species.prevalence = .99)
plotSuitabilityToProba(sp4)
sp4 <- convertToPA(sp4,
                   PA.method = "probability",
                   prob.method = "logistic",
                   beta = "random",
                   a = NULL,
                   b = NULL,
                   species.prevalence = .99)
plotSuitabilityToProba(sp4)
sp4 <- convertToPA(sp4,
                   PA.method = "probability",
                   prob.method = "logistic",
                   beta = "random",
                   a = NULL,
                   b = NULL,
                   species.prevalence = .5)
plotSuitabilityToProba(sp4)
sp4 <- convertToPA(sp4,
                   PA.method = "probability",
                   prob.method = "logistic",
                   beta = .5,
                   a = NULL,
                   b = NULL)
plotSuitabilityToProba(sp4)


for(i in 1:100){
  sp5 <- generateRandomSp(env)
}

env <- c(rast(a * dnorm(1:100, 50, sd = 25)),
         rast(a * 1:100),
         rast(a * logisticFun(1:100, alpha = 10, beta = 70)),
         rast(t(a)),
         rast(exp(a)),
         rast(log(a)))
names(env) <- paste("Var", 1:6, sep = "")

for(i in 1:100){
  sp5 <- generateRandomSp(env)
}

env <- c(rast(a * dnorm(1:100, 50, sd = 25)),
         rast(a * 1:100),
         rast(a * logisticFun(1:100, alpha = 10, beta = 70)),
         rast(t(a)),
         rast(exp(a)),
         rast(log(a)))
names(env) <- paste("Var", 1:6, sep = "")

sp5 <- generateRandomSp(env)

samp1 <- sampleOccurrences(sp5,
                           n = 50)
samp1 <- sampleOccurrences(sp5,
                           n = 50,
                           type = "presence-absence")
samp1 <- sampleOccurrences(sp5,
                           n = 50,
                           type = "presence-absence",
                           extract.probability = TRUE)
samp1 <- sampleOccurrences(sp5,
                           n = 50,
                           type = "presence-absence",
                           sample.prevalence = .9,
                           extract.probability = TRUE)
samp1 <- sampleOccurrences(sp5,
                           n = 50,
                           type = "presence-absence",
                           sample.prevalence = .1,
                           extract.probability = TRUE)

# 
# 
# worldclim <- geodata::worldclim_global(var = "bio", res = 10, path = tempdir())
# names(worldclim) <- paste0("bio", 1:19)
# 
# my.stack <- worldclim[[c("bio2", "bio5", "bio6", "bio12", "bio13", "bio14")]]
# random.sp <- generateSpFromPCA(my.stack,
#                                axes = 1:3,
#                                niche.breadth = "narrow")
# 
# random.sp <- convertToPA(random.sp)
# 
# worldmap <- rnaturalearth::ne_countries(returnclass = "sf")
# 
# samp1 <- sampleOccurrences(random.sp,
#                            n = 50)
# samp1 <- sampleOccurrences(random.sp,
#                            n = 50,
#                            sampling.area = "Morocco",
#                            error.probability = 0.1,
#                            detection.probability = .9)
# 
# samp1 <- sampleOccurrences(random.sp,
#                            n = 50,
#                            sampling.area = worldmap[worldmap$sovereignt == 
#                                                       "France", ],
#                            error.probability = 0.1,
#                            detection.probability = .9)
# 
# 
# 
# samp1 <- sampleOccurrences(random.sp,
#                            n = 50,
#                            sampling.area = ext(0, 180, 0, 90),
#                            error.probability = 0.1,
#                            detection.probability = .9)
# 
# 
# samp1 <- sampleOccurrences(random.sp,
#                            n = 50,
#                            sampling.area = ext(0, 180, 0, 90),
#                            bias = "extent",
#                            bias.strength = 50,
#                            error.probability = 0.1,
#                            detection.probability = .9)
# 
# 
# samp1 <- sampleOccurrences(random.sp,
#                            n = 50,
#                            sampling.area = ext(0, 180, 0, 90),
#                            bias = "country",
#                            bias.area = "Egypt",
#                            bias.strength = 50,
#                            error.probability = 0.1,
#                            detection.probability = .9)
# samp1 <- sampleOccurrences(random.sp,
#                            n = 50,
#                            sampling.area = ext(0, 180, 0, 90),
#                            bias = "region",
#                            bias.area = "Africa",
#                            bias.strength = 50,
#                            error.probability = 0.1,
#                            detection.probability = .9)
# 
# samp1 <- sampleOccurrences(random.sp,
#                            n = 50,
#                            sampling.area = ext(0, 180, 0, 90),
#                            bias = "continent",
#                            bias.area = "Africa",
#                            bias.strength = 50,
#                            error.probability = 0.1,
#                            detection.probability = .9)
# 
# samp1 <- sampleOccurrences(random.sp,
#                            n = 50,
#                            sampling.area = ext(0, 180, 0, 90),
#                            bias = "polygon",
#                            bias.area = worldmap[worldmap$sovereignt == 
#                                                   "Egypt", ],
#                            bias.strength = 200,
#                            error.probability = 0.1,
#                            detection.probability = .9)
# 
# 
# samp1 <- sampleOccurrences(random.sp,
#                            n = 50,
#                            sampling.area = ext(0, 180, 0, 90),
#                            bias = "polygon",
#                            bias.area = NULL,
#                            bias.strength = 200,
#                            error.probability = 0.1,
#                            detection.probability = .9)
# 
# 
# samp1 <- sampleOccurrences(random.sp,
#                            n = 50,
#                            sampling.area = ext(0, 180, 0, 90),
#                            bias = "extent",
#                            bias.area = NULL,
#                            bias.strength = 200,
#                            error.probability = 0.1,
#                            detection.probability = .9)
# 
# 
# samp1 <- sampleOccurrences(random.sp,
#                            type = "presence-absence",
#                            n = 50,
#                            sampling.area = ext(0, 180, 0, 90),
#                            bias = "extent",
#                            bias.area = NULL,
#                            bias.strength = 200,
#                            error.probability = 0.1,
#                            detection.probability = .9)
# 
# samp1 <- sampleOccurrences(random.sp,
#                            type = "presence-absence",
#                            n = 50,
#                            error.probability = 0.1,
#                            detection.probability = .9,
#                            correct.by.suitability = TRUE)
# 
# 
# samp1 <- sampleOccurrences(random.sp,
#                            type = "presence-absence",
#                            n = 50,
#                            error.probability = 0.1,
#                            detection.probability = .9,
#                            correct.by.suitability = TRUE,
#                            bias = "manual",
#                            weights = exp(random.sp$suitab.raster))
