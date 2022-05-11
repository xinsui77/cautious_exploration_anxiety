library(ggplot2)
library(dplyr)
library(readr)

sim_collect = read.csv("/Users/xinsui/Google Drive/cpi/server/logs/collect_sim.csv", header = TRUE, sep = ",")
block_collect = read.csv("/Users/xinsui/Google Drive/cpi/server/logs/collect_block.csv", header = TRUE, sep = ",")
trial_collect = read.csv("/Users/xinsui/Google Drive/cpi/server/logs/collect_trial.csv", header = TRUE, sep = ",")

names(sim_collect) <- c("nsim", "participant", "threshold")
block_collect <- block_collect[,2:7]
trial_collect <- trial_collect[,2:10]

df <- merge(x=trial_collect, y=sim_collect[,c(1,3)], by="nsim")
df$y <- ifelse(df$arm==1, 1, 0) 
df$v <- df$mu1 - df$mu2
df$ru <- sqrt(df$var1) - sqrt(df$var2)
df$tu <- sqrt(df$var1 + df$var2)
df$v_tu <- df$v/df$tu

thresholds <- unique(df$threshold)
weights <- expand.grid(threshold=thresholds,
                       w0=0,
                       w_V=0,
                       w_RU=0,
                       w_VbyTU=0,
                       w0_sig=0)


test <- df[df$threshold==400,]
reg <- glm(y ~ v + ru + v_tu,
                  family = binomial(link = "probit"),
                  data = test)  

summary(reg)
coef(reg)

weights[9,2:5] <- coef(reg)





# for (t in thresholds){
#   haha <- df[df$threshold==t,]
#   probit_reg <- glm(y ~ v + ru + v_tu,
#                     family = binomial(link = "probit"),
#                     data = haha)
#   weights$w0[weights$threshold==t,] <- coef(probit_reg)[1]
#   weights$w_V[weights$threshold==t,] <- coef(probit_reg)[2]
#   weights$w_RU[weights$threshold==t,] <- coef(probit_reg)[3]
#   weights$w_VbyTU[weights$threshold==t,] <- coef(probit_reg)[4]
#   remove(probit_reg)
#   remove(haha)
# }

ggplot() +
  #geom_line(aes(color=threshold, group=threshold)) +
  geom_point(weights, mapping = aes(x=threshold, y=w_V, color = "V")) +
  geom_line(weights, mapping = aes(x=threshold, y=w_V, color = "V")) +
  geom_point(weights, mapping = aes(x=threshold, y=w_RU, color = "RU")) +
  geom_line(weights, mapping = aes(x=threshold, y=w_RU, color = "RU")) +
  geom_point(weights, mapping = aes(x=threshold, y=w_VbyTU, color = "V/TU")) +
  geom_line(weights, mapping = aes(x=threshold, y=w_VbyTU, color = "V/TU")) +
  labs(x = "Threshold", 
       y = "Probit Regression Coefficients",
       color = "Regressor") +
  theme_minimal(base_size = 16)











#### input data from server simulation output files
data = read.csv("/Users/xinsui/Google Drive/cpi/server/logs/finished/all_params_0205.csv", header = FALSE, sep = "")


### input data from MANY server simulation output files
data <- list.files(path="/Users/xinsui/Google Drive/cpi/server/results_187837/", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows 


### cleaning up data file
data = data[,2:4]
names(data) <-  c("sigma", "lambda", "ln_bestT")
data$bestT <- exp(data$ln_bestT)
data$lambda_c <- as.character(data$lambda)
data$lambda_rounded <- as.character(as.integer(data$lambda))
data$TU <- sqrt(2)*data$sigma

#### heatmap
ggplot(data, aes(x=lambda, y=sigma, fill= bestT)) + 
  geom_tile() +
  geom_text(aes(label = round(bestT, 1)), color = "white") +
  labs(fill = "Optimal Threshold", x = "Intolerance to Uncertainty (\u03BB)", y = "Uncertainty in the Environment (\u03C3)") +
  theme_minimal(base_size = 14) # + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#### input data from meanSD simulation output files
# data_meanSD = read.csv("/Users/xinsui/Google Drive/a_semester_3/cpi/server/logs/small/out/all_params_0215.csv", header = FALSE, sep = "")
# data_meanSD = data_meanSD[,2:4]
# names(data_meanSD) <-  c("sigma", "lambda", "log_best_t")
# data_meanSD$actual_best_t <- exp(data_meanSD$log_best_t)

# ggplot(data_meanSD, aes(x=lambda, y=sigma, fill= log_best_t)) + 
#   geom_tile() +
#   geom_text(aes(label = round(log_best_t, 2)), color = "white")


ggplot(data[data$lambda<10,], aes(x=sigma, y=bestT, shape=lambda_rounded, color=lambda_rounded)) + 
  geom_point() + geom_smooth(method="lm", se=FALSE) 

ggplot(data[data$sigma > 10,], aes(x=lambda, y=sigma, fill= exp_t)) + 
  geom_tile() +
  

#### correlation plots
library("ggpubr")
ggscatter(data, x = "bestT", y = "TU", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", #cor.coef.coord = c(12,300),
          xlab = "Optimal Threshold", ylab = "TU")

ggscatter(data, x = "lambda", y = "actual_best_t", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "lambda", ylab = "actual_best_t")

ggscatter(data[data$sigma <= 10,], x = "sigma", y = "actual_best_t", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "sigma", ylab = "actual_best_t")

ggscatter(data, x = "lambda", y = "sigma", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "lambda", ylab = "sigma")

hist(data$sigma, breaks = 20)

#### linear regression to check on interaction effect between sigma and lambda
lm_log_t = lm(log_best_t ~ sigma + lambda + sigma*lambda, data[])
summary(lm_log_t) 

install.packages("car")
library("car")
data("Ornstein")
mod <- lm(interlocks ~ log(assets), data=Ornstein)
newd <- data.frame(assets=exp(4:12))



tbl_predict <- data.frame(crossing(sigma = c(5L, 10L, 15L), lambda = c(2.5, 5, 7.5)))

predicted_data <- tibble::tibble(sigma = c(5,10,15), lambda = c())
tmp <- predict(lm_log_t, tbl_predict, interval = 'prediction')
tbl_predict$pred_y <- predict(lm_log_t, tbl_predict, interval = 'prediction')[,1]

ggplot(tbl_predict, aes(lambda, pred_y, group = sigma)) +
  geom_point(aes(color = sigma)) +
  geom_line(aes(color = sigma))

?predict

lm_actual_t = lm(actual_best_t ~ sigma + lambda + sigma*lambda, data)
summary(lm_actual_t) 


library(tidyverse)

predict(mod, newd, interval="prediction")
