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

plotResponse(x = env, parameters = parameters, rescale = TRUE, approach = "response")


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
