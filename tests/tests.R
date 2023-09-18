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

sp4 <- convertToPA(sp4,
                   PA.method = "threshold")
sp4 <- convertToPA(sp4,
                   PA.method = "threshold",
                   beta = 0.5)
sp4 <- convertToPA(sp4,
                   PA.method = "threshold",
                   beta = 0)
sp4 <- convertToPA(sp4,
                   PA.method = "threshold",
                   beta = 1)
sp4 <- convertToPA(sp4,
                   PA.method = "threshold",
                   beta = "random")
sp4 <- convertToPA(sp4,
                   PA.method = "threshold",
                   species.prevalence = .2)
sp4 <- convertToPA(sp4,
                   PA.method = "threshold",
                   species.prevalence = .02)
sp4 <- convertToPA(sp4,
                   PA.method = "threshold",
                   species.prevalence = .9)
sp4 <- convertToPA(sp4,
                   PA.method = "threshold",
                   species.prevalence = .99)
sp4 <- convertToPA(sp4,
                   PA.method = "probability",
                   prob.method = "logistic",
                   beta = "random",
                   a = NULL,
                   b = NULL,
                   species.prevalence = .99)
sp4 <- convertToPA(sp4,
                   PA.method = "probability",
                   prob.method = "logistic",
                   beta = "random",
                   a = NULL,
                   b = NULL,
                   species.prevalence = .5)
sp4 <- convertToPA(sp4,
                   PA.method = "probability",
                   prob.method = "logistic",
                   beta = .5,
                   a = NULL,
                   b = NULL)
