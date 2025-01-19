##############################
#  Load Libraries and Setup  #
##############################
library(shiny)
library(ggplot2)
library(dplyr)
library(simputation)
library(patchwork)
library(gridExtra)

options(scipen = 999)  # Turn off scientific notation

##################################################
#  (1) Code for All Scenarios: MCAR, MAR, MNAR   #
##################################################

############
#  MCAR
############
set.seed(124)

# (A) Generate MCAR data
generate_mcar_dataset <- function(n = 1000) {
  X <- rnorm(n, mean = 30, sd = 10)
  Y <- abs(X ^ 4 + rnorm(n, mean = 0, sd = 200000))
  data.frame(X = X, Y = Y)
}

# Generate 10 datasets with MCAR
mcar_datasets <- lapply(1:10, function(i) generate_mcar_dataset())

# (B) Process MCAR datasets (introduce missingness, impute, etc.)
process_mcar_dataset <- function(dataset) {
  mcar <- rbinom(nrow(dataset), 1, 0.3)
  cX <- dataset$X
  dataset <- cbind(dataset, cX, mcar)
  dataset$X[mcar == 1] <- NA
  
  dataset$cXm  <- dataset$X
  dataset$Xcim <- dataset$X
  
  # Uncongenial
  dataset <- impute_lm(dataset, X ~ Y)
  # Congenial
  dataset <- impute_lm(dataset, Xcim ~ Y + I(Y ^ (1 / 4)))
  
  dataset
}
mcar_processed_datasets <- lapply(mcar_datasets, process_mcar_dataset)

# (C) Fit MCAR models
fit_mcar_models <- function(dataset) {
  lm0 <- lm(Y ~ cX + I(cX ^ 4), data = dataset)
  lm1 <- lm(Y ~ X + I(X ^ 4), data = dataset)
  lm2 <- lm(Y ~ cXm + I(cXm ^ 4), data = dataset, na.action = na.omit)
  lm3 <- lm(Y ~ Xcim + I(Xcim ^ 4), data = dataset)
  list(lm0 = lm0, lm1 = lm1, lm2 = lm2, lm3 = lm3)
}
mcar_models <- lapply(mcar_processed_datasets, fit_mcar_models)

# (D) Generate MCAR test sets
generate_mcar_testset <- function(n = 10000) {
  X <- rnorm(n, mean = 30, sd = 10)
  Y <- abs(X ^ 4 + rnorm(n, mean = 0, sd = 200000))
  cX <- cXm <- Xcim <- X
  data.frame(X = X, Y = Y, cX = cX, cXm = cXm, Xcim = Xcim)
}
mcar_testsets <- lapply(1:10, function(i) generate_mcar_testset())

# (E) Predict on MCAR test sets
predict_mcar <- function(testset, models) {
  testset$Y_true        <- predict(models$lm0, newdata = testset)
  testset$Y_uncon_impute <- predict(models$lm1, newdata = testset)
  testset$Y_cc          <- predict(models$lm2, newdata = testset)
  testset$Y_con_impute  <- predict(models$lm3, newdata = testset)
  testset
}
mcar_predicted_testsets <- Map(predict_mcar, mcar_testsets, mcar_models)

# (F) Evaluate MCAR
rmse <- function(y, y_hat) sqrt(mean((y - y_hat)^2))

calculate_rmse <- function(testset) {
  list(
    rmse_true         = rmse(testset$Y, testset$Y_true),
    rmse_uncon_impute = rmse(testset$Y, testset$Y_uncon_impute),
    rmse_cc           = rmse(testset$Y, testset$Y_cc),
    rmse_con_impute   = rmse(testset$Y, testset$Y_con_impute)
  )
}
mcar_rmse_results <- lapply(mcar_predicted_testsets, calculate_rmse)
mcar_rmse_summary <- do.call(rbind, lapply(mcar_rmse_results, as.data.frame))
rownames(mcar_rmse_summary) <- paste0("TestSet_", 1:10)

mcar_rmse_summary_df <- data.frame(
  method    = c("MCAR_True Model", "MCAR_Uncongenial Imputation", 
                "MCAR_Complete Case Analysis", "MCAR_Congenial Imputation"),
  mean_rmse = c(mean(mcar_rmse_summary$rmse_true),
                mean(mcar_rmse_summary$rmse_uncon_impute),
                mean(mcar_rmse_summary$rmse_cc),
                mean(mcar_rmse_summary$rmse_con_impute))
)

mcar_r_squared_df <- data.frame(
  method        = c("MCAR_True Model", "MCAR_Uncongenial Imputation", 
                    "MCAR_Complete Case Analysis", "MCAR_Congenial Imputation"),
  mean_r_squared = c(
    mean(sapply(mcar_models, function(x) summary(x$lm0)$r.squared)),
    mean(sapply(mcar_models, function(x) summary(x$lm1)$r.squared)),
    mean(sapply(mcar_models, function(x) summary(x$lm2)$r.squared)),
    mean(sapply(mcar_models, function(x) summary(x$lm3)$r.squared))
  )
)
mcar_results <- merge(mcar_rmse_summary_df, mcar_r_squared_df, by = "method")

# Example dataset for plotting in MCAR scenario
mcar_ex_dat      <- mcar_processed_datasets[[1]]
mcar_true_model  <- mcar_models[[1]]$lm0
mcar_uncon_model <- mcar_models[[1]]$lm1
mcar_cc_model    <- mcar_models[[1]]$lm2
mcar_con_model   <- mcar_models[[1]]$lm3


############
#  MAR
############
set.seed(124)
generate_mar_dataset <- function(n = 1000) {
  X <- rnorm(n, mean = 30, sd = 10)
  Y <- abs(X ^ 4 + rnorm(n, mean = 0, sd = 200000))
  data.frame(X = X, Y = Y)
}
mar_datasets <- lapply(1:10, function(i) generate_mar_dataset())

process_mar_dataset <- function(dataset) {
  rankY <- rank(-dataset$Y)
  mar_prob <- rankY / nrow(dataset)
  mar <- rbinom(nrow(dataset), 1, mar_prob * (0.3 / mean(mar_prob)))
  
  cX <- dataset$X
  dataset <- cbind(dataset, cX, mar)
  dataset$X[mar == 1] <- NA
  
  dataset$cXm  <- dataset$X
  dataset$Xcim <- dataset$X
  
  dataset <- impute_lm(dataset, X ~ Y)
  dataset <- impute_lm(dataset, Xcim ~ Y + I(Y ^ (1/4)))
  dataset
}
mar_processed_datasets <- lapply(mar_datasets, process_mar_dataset)

fit_mar_models <- function(dataset) {
  lm0 <- lm(Y ~ cX + I(cX ^ 4), data = dataset)
  lm1 <- lm(Y ~ X + I(X ^ 4), data = dataset)
  lm2 <- lm(Y ~ cXm + I(cXm ^ 4), data = dataset, na.action = na.omit)
  lm3 <- lm(Y ~ Xcim + I(Xcim ^ 4), data = dataset)
  list(lm0 = lm0, lm1 = lm1, lm2 = lm2, lm3 = lm3)
}
mar_models <- lapply(mar_processed_datasets, fit_mar_models)

generate_mar_testset <- function(n = 10000) {
  X <- rnorm(n, mean = 30, sd = 10)
  Y <- abs(X ^ 4 + rnorm(n, mean = 0, sd = 200000))
  cX <- cXm <- Xcim <- X
  data.frame(X = X, Y = Y, cX = cX, cXm = cXm, Xcim = Xcim)
}
mar_testsets <- lapply(1:10, function(i) generate_mar_testset())

predict_mar <- function(testset, models) {
  testset$Y_true         <- predict(models$lm0, newdata = testset)
  testset$Y_uncon_impute <- predict(models$lm1, newdata = testset)
  testset$Y_cc           <- predict(models$lm2, newdata = testset)
  testset$Y_con_impute   <- predict(models$lm3, newdata = testset)
  testset
}
mar_predicted_testsets <- Map(predict_mar, mar_testsets, mar_models)

mar_rmse_results <- lapply(mar_predicted_testsets, calculate_rmse)
mar_rmse_summary <- do.call(rbind, lapply(mar_rmse_results, as.data.frame))
rownames(mar_rmse_summary) <- paste0("TestSet_", 1:10)

mar_rmse_summary_df <- data.frame(
  method    = c("MAR_True Model", "MAR_Uncongenial Imputation", 
                "MAR_Complete Case Analysis", "MAR_Congenial Imputation"),
  mean_rmse = c(mean(mar_rmse_summary$rmse_true),
                mean(mar_rmse_summary$rmse_uncon_impute),
                mean(mar_rmse_summary$rmse_cc),
                mean(mar_rmse_summary$rmse_con_impute))
)

mar_r_squared_df <- data.frame(
  method        = c("MAR_True Model", "MAR_Uncongenial Imputation", 
                    "MAR_Complete Case Analysis", "MAR_Congenial Imputation"),
  mean_r_squared = c(
    mean(sapply(mar_models, function(x) summary(x$lm0)$r.squared)),
    mean(sapply(mar_models, function(x) summary(x$lm1)$r.squared)),
    mean(sapply(mar_models, function(x) summary(x$lm2)$r.squared)),
    mean(sapply(mar_models, function(x) summary(x$lm3)$r.squared))
  )
)
mar_results <- merge(mar_rmse_summary_df, mar_r_squared_df, by = "method")

mar_ex_dat      <- mar_processed_datasets[[1]]
mar_true_model  <- mar_models[[1]]$lm0
mar_uncon_model <- mar_models[[1]]$lm1
mar_cc_model    <- mar_models[[1]]$lm2
mar_con_model   <- mar_models[[1]]$lm3


############
#  MNAR
############
set.seed(124)
generate_mnar_dataset <- function(n = 1000) {
  X <- rnorm(n, mean = 30, sd = 10)
  Y <- abs(X ^ 4 + rnorm(n, mean = 0, sd = 200000))
  data.frame(X = X, Y = Y)
}
mnar_datasets <- lapply(1:10, function(i) generate_mnar_dataset())

process_mnar_dataset <- function(dataset) {
  rankX <- rank(dataset$X)
  mnar_prob <- rankX / nrow(dataset)
  mnar <- rbinom(nrow(dataset), 1, mnar_prob * (0.3 / mean(mnar_prob)))
  
  cX <- dataset$X
  dataset <- cbind(dataset, cX, mnar)
  dataset$X[mnar == 1] <- NA
  
  dataset$cXm  <- dataset$X
  dataset$Xcim <- dataset$X
  
  dataset <- impute_lm(dataset, X ~ Y)
  dataset <- impute_lm(dataset, Xcim ~ Y + I(Y ^ (1/4)))
  dataset
}
mnar_processed_datasets <- lapply(mnar_datasets, process_mnar_dataset)

fit_mnar_models <- function(dataset) {
  lm0 <- lm(Y ~ cX + I(cX ^ 4), data = dataset)
  lm1 <- lm(Y ~ X + I(X ^ 4), data = dataset)
  lm2 <- lm(Y ~ cXm + I(cXm ^ 4), data = dataset, na.action = na.omit)
  lm3 <- lm(Y ~ Xcim + I(Xcim ^ 4), data = dataset)
  list(lm0 = lm0, lm1 = lm1, lm2 = lm2, lm3 = lm3)
}
mnar_models <- lapply(mnar_processed_datasets, fit_mnar_models)

generate_mnar_testset <- function(n = 10000) {
  X <- rnorm(n, mean = 30, sd = 10)
  Y <- abs(X ^ 4 + rnorm(n, mean = 0, sd = 200000))
  cX <- cXm <- Xcim <- X
  data.frame(X = X, Y = Y, cX = cX, cXm = cXm, Xcim = Xcim)
}
mnar_testsets <- lapply(1:10, function(i) generate_mnar_testset())

predict_mnar <- function(testset, models) {
  testset$Y_true        <- predict(models$lm0, newdata = testset)
  testset$Y_uncon_impute <- predict(models$lm1, newdata = testset)
  testset$Y_cc          <- predict(models$lm2, newdata = testset)
  testset$Y_con_impute  <- predict(models$lm3, newdata = testset)
  testset
}
mnar_predicted_testsets <- Map(predict_mnar, mnar_testsets, mnar_models)

mnar_rmse_results <- lapply(mnar_predicted_testsets, calculate_rmse)
mnar_rmse_summary <- do.call(rbind, lapply(mnar_rmse_results, as.data.frame))
rownames(mnar_rmse_summary) <- paste0("TestSet_", 1:10)

mnar_rmse_summary_df <- data.frame(
  method    = c("MNAR_True Model", "MNAR_Uncongenial Imputation", 
                "MNAR_Complete Case Analysis", "MNAR_Congenial Imputation"),
  mean_rmse = c(mean(mnar_rmse_summary$rmse_true),
                mean(mnar_rmse_summary$rmse_uncon_impute),
                mean(mnar_rmse_summary$rmse_cc),
                mean(mnar_rmse_summary$rmse_con_impute))
)

mnar_r_squared_df <- data.frame(
  method        = c("MNAR_True Model", "MNAR_Uncongenial Imputation", 
                    "MNAR_Complete Case Analysis", "MNAR_Congenial Imputation"),
  mean_r_squared = c(
    mean(sapply(mnar_models, function(x) summary(x$lm0)$r.squared)),
    mean(sapply(mnar_models, function(x) summary(x$lm1)$r.squared)),
    mean(sapply(mnar_models, function(x) summary(x$lm2)$r.squared)),
    mean(sapply(mnar_models, function(x) summary(x$lm3)$r.squared))
  )
)
mnar_results <- merge(mnar_rmse_summary_df, mnar_r_squared_df, by = "method")

mnar_ex_dat      <- mnar_processed_datasets[[1]]
mnar_true_model  <- mnar_models[[1]]$lm0
mnar_uncon_model <- mnar_models[[1]]$lm1
mnar_cc_model    <- mnar_models[[1]]$lm2
mnar_con_model   <- mnar_models[[1]]$lm3


##########################################
#  (2) Shiny UI Definition
##########################################
ui <- fluidPage(
  titlePanel("Uncongeniality App"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      # Scenario selection
      selectInput(
        inputId = "scenario",
        label   = "Choose Scenario:",
        choices = c("MCAR", "MAR", "MNAR"),
        selected = "MCAR"
      ),
      
      # Plot/Table selection
      radioButtons(
        inputId = "outputType",
        label   = "What do you want to see?",
        choices = c(
          "Data Plot",
          "Prediction Models",
          "Calibration Lines",
          "RMSE & R2 Table"
        ),
        selected = "Data Plot"
      ),
      
      # Explanation text
      div(
        p("This app demonstrates the impact of uncongeniality between imputation models and substantive models in a simplified setting. The setting involves only 2 variables (X and Y) that have a non-linear relationship (Y = X^4 + noise). There are different scenarios of missingness mechanisms (MCAR, MAR, and MNAR)."),
        p("You can choose to look at 3 plots and one table for each scenario. 
           The Data Plot shows the relationship between X and Y, with missing values indicated by color. 
           The Prediction Models plot shows the trained prediction models for uncongenial and congenial imputations. 
           The Calibration Lines plot shows the relationship between the predicted values and the true values. 
           The RMSE & R2 Table shows the average RMSE and R^2 values for each model across 10 test sets. 
           For a more in-depth explanation, please read my report.")
      )
    ),
    mainPanel(
      width = 9,
      uiOutput("plotOrTableUI")
    )
  )
)


##########################################
#  (3) Shiny Server Logic
##########################################
server <- function(input, output, session) {
  
  # Helper function: dataPlot
  dataPlot <- function(scenario) {
    if (scenario == "MCAR") {
      ggplot(mcar_ex_dat, aes(x = cX, y = Y)) +
        geom_point(aes(color = factor(mcar))) +
        scale_color_manual(values = c("0" = "blue", "1" = "red")) +
        theme_classic() +
        labs(title = "MCAR Data", x = "X", y = "Y", color = "Missing Indicator") +
        theme(legend.position = "bottom")
    } else if (scenario == "MAR") {
      ggplot(mar_ex_dat, aes(x = cX, y = Y)) +
        geom_point(aes(color = factor(mar))) +
        scale_color_manual(values = c("0" = "blue", "1" = "red")) +
        theme_classic() +
        labs(title = "MAR Data", x = "X", y = "Y", color = "Missing Indicator") +
        theme(legend.position = "bottom")
    } else {
      ggplot(mnar_ex_dat, aes(x = cX, y = Y)) +
        geom_point(aes(color = factor(mnar))) +
        scale_color_manual(values = c("0" = "blue", "1" = "red")) +
        theme_classic() +
        labs(title = "MNAR Data", x = "X", y = "Y", color = "Missing Indicator") +
        theme(legend.position = "bottom")
    }
  }
  
  # Helper function: predictionModelsPlot
  predictionModelsPlot <- function(scenario) {
    if (scenario == "MCAR") {
      p1 <- ggplot(mcar_ex_dat, aes(x = X, y = Y)) +
        geom_point(aes(color = factor(mcar))) +
        scale_color_manual(values = c("0" = "blue", "1" = "red")) +
        stat_function(fun = function(x) {
          coef(mcar_uncon_model)[1] + coef(mcar_uncon_model)[2]*x + coef(mcar_uncon_model)[3]*x^4
        }, color = "black") +
        stat_function(fun = function(x) {
          coef(mcar_true_model)[1] + coef(mcar_true_model)[2]*x + coef(mcar_true_model)[3]*x^4
        }, color = "red", linetype = "dashed") +
        labs(title = "MCAR - Uncongenial vs. True", x = "X", y = "Y") +
        theme_classic() +
        theme(legend.position = "none")
      
      p2 <- ggplot(mcar_ex_dat, aes(x = Xcim, y = Y)) +
        geom_point(aes(color = factor(mcar))) +
        scale_color_manual(values = c("0" = "blue", "1" = "red")) +
        stat_function(fun = function(x) {
          coef(mcar_con_model)[1] + coef(mcar_con_model)[2]*x + coef(mcar_con_model)[3]*x^4
        }, color = "black") +
        stat_function(fun = function(x) {
          coef(mcar_true_model)[1] + coef(mcar_true_model)[2]*x + coef(mcar_true_model)[3]*x^4
        }, color = "red", linetype = "dashed") +
        labs(title = "MCAR - Congenial vs. True", x = "X", y = "Y") +
        theme_classic() +
        theme(legend.position = "none")
      
      p1 + p2
      
    } else if (scenario == "MAR") {
      p1 <- ggplot(mar_ex_dat, aes(x = X, y = Y)) +
        geom_point(aes(color = factor(mar))) +
        scale_color_manual(values = c("0" = "blue", "1" = "red")) +
        stat_function(fun = function(x) {
          coef(mar_uncon_model)[1] + coef(mar_uncon_model)[2]*x + coef(mar_uncon_model)[3]*x^4
        }, color = "black") +
        stat_function(fun = function(x) {
          coef(mar_true_model)[1] + coef(mar_true_model)[2]*x + coef(mar_true_model)[3]*x^4
        }, color = "red", linetype = "dashed") +
        labs(title = "MAR - Uncongenial vs. True", x = "X", y = "Y") +
        theme_classic() +
        theme(legend.position = "none")
      
      p2 <- ggplot(mar_ex_dat, aes(x = Xcim, y = Y)) +
        geom_point(aes(color = factor(mar))) +
        scale_color_manual(values = c("0" = "blue", "1" = "red")) +
        stat_function(fun = function(x) {
          coef(mar_con_model)[1] + coef(mar_con_model)[2]*x + coef(mar_con_model)[3]*x^4
        }, color = "black") +
        stat_function(fun = function(x) {
          coef(mar_true_model)[1] + coef(mar_true_model)[2]*x + coef(mar_true_model)[3]*x^4
        }, color = "red", linetype = "dashed") +
        labs(title = "MAR - Congenial vs. True", x = "X", y = "Y") +
        theme_classic() +
        theme(legend.position = "none")
      
      p1 + p2
      
    } else {
      # MNAR
      p1 <- ggplot(mnar_ex_dat, aes(x = X, y = Y)) +
        geom_point(aes(color = factor(mnar))) +
        scale_color_manual(values = c("0" = "blue", "1" = "red")) +
        stat_function(fun = function(x) {
          coef(mnar_uncon_model)[1] + coef(mnar_uncon_model)[2]*x + coef(mnar_uncon_model)[3]*x^4
        }, color = "black") +
        stat_function(fun = function(x) {
          coef(mnar_true_model)[1] + coef(mnar_true_model)[2]*x + coef(mnar_true_model)[3]*x^4
        }, color = "red", linetype = "dashed") +
        labs(title = "MNAR - Uncongenial vs. True", x = "X", y = "Y") +
        theme_classic() +
        theme(legend.position = "none")
      
      p2 <- ggplot(mnar_ex_dat, aes(x = Xcim, y = Y)) +
        geom_point(aes(color = factor(mnar))) +
        scale_color_manual(values = c("0" = "blue", "1" = "red")) +
        stat_function(fun = function(x) {
          coef(mnar_con_model)[1] + coef(mnar_con_model)[2]*x + coef(mnar_con_model)[3]*x^4
        }, color = "black") +
        stat_function(fun = function(x) {
          coef(mnar_true_model)[1] + coef(mnar_true_model)[2]*x + coef(mnar_true_model)[3]*x^4
        }, color = "red", linetype = "dashed") +
        labs(title = "MNAR - Congenial vs. True", x = "X", y = "Y") +
        theme_classic() +
        theme(legend.position = "none")
      
      p1 + p2
    }
  }
  
  # Helper: calibrationPlot
  calibrationPlot <- function(scenario) {
    create_calibration <- function(true_values, predicted_values, bins = 10) {
      bin <- cut(predicted_values, breaks = bins)
      data.frame(
        PredictedMean = tapply(predicted_values, bin, mean, na.rm = TRUE),
        ObservedMean  = tapply(true_values, bin, mean, na.rm = TRUE)
      )
    }
    
    if (scenario == "MCAR") {
      mcar_cc_list_uncon <- lapply(seq_len(10), function(i) {
        test_i <- as.data.frame(mcar_predicted_testsets[[i]])
        cc_i <- create_calibration(test_i$Y, test_i$Y_uncon_impute)
        cc_i$Model <- paste0("Testset ", i)
        cc_i$Type  <- "uncon"
        cc_i
      })
      mcar_cc_list_con <- lapply(seq_len(10), function(i) {
        test_i <- as.data.frame(mcar_predicted_testsets[[i]])
        cc_i <- create_calibration(test_i$Y, test_i$Y_con_impute)
        cc_i$Model <- paste0("Testset ", i)
        cc_i$Type  <- "con"
        cc_i
      })
      df_comb <- bind_rows(mcar_cc_list_uncon, mcar_cc_list_con)
      df_comb$Type <- factor(df_comb$Type, levels = c("uncon", "con"),
                             labels = c("uncongenial", "congenial"))
      
      ggplot(df_comb, aes(x = PredictedMean, y = ObservedMean, 
                          group = interaction(Model, Type), color = Type)) +
        geom_point() + geom_line() + 
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        labs(title = "Calibration Plots for MCAR", x = "Mean Predicted", y = "Mean Observed") +
        theme_classic() +
        scale_color_manual(values = c("uncongenial" = "red", "congenial" = "green")) +
        theme(legend.position = "bottom")
      
    } else if (scenario == "MAR") {
      mar_cc_list_uncon <- lapply(seq_len(10), function(i) {
        test_i <- as.data.frame(mar_predicted_testsets[[i]])
        cc_i <- create_calibration(test_i$Y, test_i$Y_uncon_impute)
        cc_i$Model <- paste0("Testset ", i)
        cc_i$Type  <- "uncon"
        cc_i
      })
      mar_cc_list_con <- lapply(seq_len(10), function(i) {
        test_i <- as.data.frame(mar_predicted_testsets[[i]])
        cc_i <- create_calibration(test_i$Y, test_i$Y_con_impute)
        cc_i$Model <- paste0("Testset ", i)
        cc_i$Type  <- "con"
        cc_i
      })
      df_comb <- bind_rows(mar_cc_list_uncon, mar_cc_list_con)
      df_comb$Type <- factor(df_comb$Type, levels = c("uncon", "con"),
                             labels = c("uncongenial", "congenial"))
      
      ggplot(df_comb, aes(x = PredictedMean, y = ObservedMean, 
                          group = interaction(Model, Type), color = Type)) +
        geom_point() + geom_line() + 
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        labs(title = "Calibration Plots for MAR", x = "Mean Predicted", y = "Mean Observed") +
        theme_classic() +
        scale_color_manual(values = c("uncongenial" = "red", "congenial" = "green")) +
        theme(legend.position = "bottom")
      
    } else {
      # MNAR
      mnar_cc_list_uncon <- lapply(seq_len(10), function(i) {
        test_i <- as.data.frame(mnar_predicted_testsets[[i]])
        cc_i <- create_calibration(test_i$Y, test_i$Y_uncon_impute)
        cc_i$Model <- paste0("Testset ", i)
        cc_i$Type  <- "uncon"
        cc_i
      })
      mnar_cc_list_con <- lapply(seq_len(10), function(i) {
        test_i <- as.data.frame(mnar_predicted_testsets[[i]])
        cc_i <- create_calibration(test_i$Y, test_i$Y_con_impute)
        cc_i$Model <- paste0("Testset ", i)
        cc_i$Type  <- "con"
        cc_i
      })
      df_comb <- bind_rows(mnar_cc_list_uncon, mnar_cc_list_con)
      df_comb$Type <- factor(df_comb$Type, levels = c("uncon", "con"),
                             labels = c("uncongenial", "congenial"))
      
      ggplot(df_comb, aes(x = PredictedMean, y = ObservedMean, 
                          group = interaction(Model, Type), color = Type)) +
        geom_point() + geom_line() + 
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        labs(title = "Calibration Plots for MNAR", x = "Mean Predicted", y = "Mean Observed") +
        theme_classic() +
        scale_color_manual(values = c("uncongenial" = "red", "congenial" = "green")) +
        theme(legend.position = "bottom")
    }
  }
  
  # Explanation text helper
  explanationText <- function(scenario, outType) {
    if (outType == "Data Plot") {
      switch(scenario,
             "MCAR" = "Notes: Missingness in X is completely random. Missings are spread all over the data.",
             "MAR"  = "Notes: Observations with a low Y-value are more likely to be missing. Missings are in the lower left corner.",
             "MNAR" = "Notes: Observations with a high X-value are more likely to be missing. Missings are in the top right corner."
      )
    } else if (outType == "Prediction Models") {
      switch(scenario,
             "MCAR" = "Notes: The uncongenial imputation (red points) can deviate significantly from the true model. Congenial is closer.",
             "MAR"  = "Notes: The uncongenial imputation is off mainly in the bottom left corner. Congenial is closer to the true model.",
             "MNAR" = "Notes: Uncongenial imputation is off in the top right corner. Congenial is closer to the true model."
      )
    } else if (outType == "Calibration Lines") {
      switch(scenario,
             "MCAR" = "Notes: Uncongenial predicted values are generally higher than observed. Congenial is closer to the 45-degree line.",
             "MAR"  = "Notes: Uncongenial and congenial are both fairly close to the 45-degree line.",
             "MNAR" = "Notes: Uncongenial is higher than observed. Congenial is closer to the 45-degree line."
      )
    } else if (outType == "RMSE & R2 Table") {
      switch(scenario,
             "MCAR" = "Notes: RMSE is higher, R2 lower for uncongenial compared to congenial. Congenial is more accurate.",
             "MAR"  = "Notes: Slightly higher RMSE, lower R2 for uncongenial vs. congenial. Difference is smaller than MCAR.",
             "MNAR" = "Notes: RMSE is higher, R2 is lower for uncongenial. Congenial is more accurate."
      )
    }
  }
  
  # metricsTable helper
  metricsTable <- function(scenario) {
    if (scenario == "MCAR") {
      mcar_results
    } else if (scenario == "MAR") {
      mar_results
    } else {
      mnar_results
    }
  }
  
  # Output: UI that either shows a plot or a table + explanation
  output$plotOrTableUI <- renderUI({
    scenarioChoice <- input$scenario
    outType        <- input$outputType
    
    mainOutput <- if (outType == "RMSE & R2 Table") {
      tableOutput("metricsTable")
    } else {
      plotOutput("scenarioPlot", height = "600px", width = "100%")
    }
    
    textBelow <- explanationText(scenarioChoice, outType)
    
    tagList(
      mainOutput,
      tags$br(),
      tags$p(style = "font-style: italic; color: gray;", textBelow)
    )
  })
  
  # Output: scenarioPlot
  output$scenarioPlot <- renderPlot({
    scenarioChoice <- input$scenario
    outType        <- input$outputType
    
    if (outType == "Data Plot") {
      dataPlot(scenarioChoice)
    } else if (outType == "Prediction Models") {
      predictionModelsPlot(scenarioChoice)
    } else if (outType == "Calibration Lines") {
      calibrationPlot(scenarioChoice)
    } else {
      NULL
    }
  })
  
  # Output: metricsTable
  output$metricsTable <- renderTable({
    scenarioChoice <- input$scenario
    metricsTable(scenarioChoice)
  })
}

##########################################
#  (4) Run the Shiny App
##########################################
shinyApp(ui = ui, server = server)
