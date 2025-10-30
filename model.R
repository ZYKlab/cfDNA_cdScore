# Colorectal Cancer 5-Marker Prediction Model
# CD Score-based Version

# 1. Load required packages --------------------------------------------------------------
library(tidyverse)
library(jsonlite)

# 2. Load trained model --------------------------------------------------------
load("fit_model.Rdata")

# 3. Define model parameters -----------------------------------------------------------
model_params <- list(
  model_name = "Colorectal Cancer 5-Marker Prediction Model",
  version = "1.0", 
  creation_date = Sys.Date(),
  markers = c("cg14015706", "cg20275528", "cg07589773", "cg19283840", "cg16935295"),
  coefficients = coef(fit.train),
  intercept = coef(fit.train)[1],
  marker_weights = coef(fit.train)[-1],
  model_type = "Logistic Regression",
  cd_threshold = 0.5
)

# 4. Define CD Score-based diagnostic functions -----------------------------------------------------------

#' CD Score-based diagnostic function
#' @param cg14015706 Methylation value of marker 1
#' @param cg20275528 Methylation value of marker 2  
#' @param cg07589773 Methylation value of marker 3
#' @param cg19283840 Methylation value of marker 4
#' @param cg16935295 Methylation value of marker 5
#' @return List containing CD Score and diagnostic results
diagnose_by_cd_score <- function(cg14015706, cg20275528, cg07589773, cg19283840, cg16935295) {
  # Calculate CD Score
  cd_score <- model_params$intercept +
    (cg14015706 * model_params$coefficients["cg14015706"]) +
    (cg20275528 * model_params$coefficients["cg20275528"]) +
    (cg07589773 * model_params$coefficients["cg07589773"]) +
    (cg19283840 * model_params$coefficients["cg19283840"]) +
    (cg16935295 * model_params$coefficients["cg16935295"])
  
  # Direct diagnosis
  diagnosis <- ifelse(cd_score >= model_params$cd_threshold, "Positive", "Negative")
  risk_level <- ifelse(cd_score >= model_params$cd_threshold, "High Risk", "Low Risk")
  
  return(list(
    cd_score = round(cd_score, 4),
    diagnosis = diagnosis,
    risk_level = risk_level,
    threshold = model_params$cd_threshold  # Using 0.5 threshold
  ))
}

#' Simplified CD Score calculation function
#' @param marker_values Vector of 5 marker values
calculate_cd_score <- function(marker_values) {
  weights <- model_params$coefficients
  cd_score <- weights[1] + sum(marker_values * weights[-1])
  
  return(round(cd_score, 4))
}

# 5. Save model files -----------------------------------------------------------

# 5.1 Save coefficients file
coefficients_df <- data.frame(
  variable = names(model_params$coefficients),
  coefficient = as.numeric(model_params$coefficients)
)
write.csv(coefficients_df, "model_coefficients.csv", row.names = FALSE)

# 5.2 Save CD Score parameters
cd_score_params <- data.frame(
  marker = c("intercept", model_params$markers),
  weight = as.numeric(model_params$coefficients),
  threshold = model_params$cd_threshold  # Using 0.5 threshold
)
write.csv(cd_score_params, "cd_score_parameters.csv", row.names = FALSE)

# 5.3 Save metadata
metadata_df <- data.frame(
  parameter = c("model_name", "version", "creation_date", "markers", "cd_threshold", "model_type"),
  value = c(
    model_params$model_name,
    model_params$version,
    as.character(model_params$creation_date),
    paste(model_params$markers, collapse = ", "),
    as.character(model_params$cd_threshold),
    model_params$model_type
  )
)
write.csv(metadata_df, "model_metadata.csv", row.names = FALSE)

# 5.4 Generate complete diagnostic script
diagnosis_script <- '
# Colorectal Cancer 5-Marker Prediction Model - CD Score-based Diagnostic Algorithm
# Model Version: 1.0
# Created: 2024

# Model coefficients and threshold
intercept <- %.6f
marker_weights <- c(
  cg14015706 = %.6f,
  cg20275528 = %.6f, 
  cg07589773 = %.6f,
  cg19283840 = %.6f,
  cg16935295 = %.6f
)
cd_threshold <- 0.5  # Diagnostic threshold set to 0.5

# CD Score-based diagnostic function
diagnose_crc <- function(cg14015706, cg20275528, cg07589773, cg19283840, cg16935295) {
  # Calculate CD Score (including intercept)
  cd_score <- intercept +
    (cg14015706 * marker_weights["cg14015706"]) +
    (cg20275528 * marker_weights["cg20275528"]) +
    (cg07589773 * marker_weights["cg07589773"]) +
    (cg19283840 * marker_weights["cg19283840"]) +
    (cg16935295 * marker_weights["cg16935295"])
  
  # Direct diagnosis (threshold = 0.5)
  diagnosis <- ifelse(cd_score >= cd_threshold, "Positive", "Negative")
  risk_level <- ifelse(cd_score >= cd_threshold, "High Risk", "Low Risk")
  
  return(list(
    cd_score = round(cd_score, 4),
    diagnosis = diagnosis,
    risk_level = risk_level,
    threshold = cd_threshold
  ))
}

# CD Score calculation only function
calculate_cd_score <- function(cg14015706, cg20275528, cg07589773, cg19283840, cg16935295) {
  cd_score <- intercept +
    (cg14015706 * marker_weights["cg14015706"]) +
    (cg20275528 * marker_weights["cg20275528"]) +
    (cg07589773 * marker_weights["cg07589773"]) +
    (cg19283840 * marker_weights["cg19283840"]) +
    (cg16935295 * marker_weights["cg16935295"])
  
  return(round(cd_score, 4))
}

# Batch diagnosis function
batch_diagnose <- function(marker_matrix) {
  # marker_matrix should be a dataframe or matrix, each row represents a sample, columns correspond to 5 markers
  results <- apply(marker_matrix, 1, function(row) {
    cd_score <- calculate_cd_score(row[1], row[2], row[3], row[4], row[5])
    diagnosis <- ifelse(cd_score >= cd_threshold, "Positive", "Negative")
    risk_level <- ifelse(cd_score >= cd_threshold, "High Risk", "Low Risk")
    return(c(cd_score = cd_score, diagnosis = diagnosis, risk_level = risk_level))
  })
  
  return(as.data.frame(t(results)))
}
'

# Format script
formatted_script <- sprintf(
  diagnosis_script,
  model_params$coefficients["(Intercept)"],
  model_params$coefficients["cg14015706"],
  model_params$coefficients["cg20275528"], 
  model_params$coefficients["cg07589773"],
  model_params$coefficients["cg19283840"],
  model_params$coefficients["cg16935295"]
)

writeLines(formatted_script, "crc_diagnosis_by_cd_score.R")

# 5.5 Save JSON file
write(toJSON(model_params, pretty = TRUE), "model_parameters.json")

# 6. Validation and testing -------------------------------------------------------------

# 6.1 Print model summary
cat("=== CD Score-based Diagnostic Model Parameter Summary ===\n")
cat("Model Name:", model_params$model_name, "\n")
cat("Version:", model_params$version, "\n") 
cat("Markers:", paste(model_params$markers, collapse = ", "), "\n")
cat("Intercept:", round(model_params$intercept, 6), "\n")
cat("CD Score Diagnostic Threshold:", model_params$cd_threshold, "\n")
cat("Marker Weights:\n")
for(i in 2:length(model_params$coefficients)) {
  cat("  ", names(model_params$coefficients)[i], ":", round(model_params$coefficients[i], 6), "\n")
}

# 6.2 Test example
cat("\n=== Example Test ===\n")
example_data <- c(0.5, 0.6, 0.7, 0.8, 0.9)
names(example_data) <- model_params$markers

example_result <- diagnose_by_cd_score(
  example_data[1], example_data[2], example_data[3], 
  example_data[4], example_data[5]
)

cat("Input Methylation Values:\n")
for(i in 1:5) {
  cat("  ", names(example_data)[i], ":", example_data[i], "\n")
}
cat("CD Score:", example_result$cd_score, "\n")
cat("Diagnostic Threshold:", example_result$threshold, "\n")
cat("Diagnosis Result:", example_result$diagnosis, "\n")
cat("Risk Level:", example_result$risk_level, "\n")

# Explain diagnostic logic
cat("\nDiagnostic Logic: CD Score >= 0.5 → Positive (High Risk)\n")
cat("                  CD Score < 0.5 → Negative (Low Risk)\n")

# 7. Generate documentation ---------------------------------------------------------------

documentation <- '
# Colorectal Cancer 5-Marker Prediction Model - CD Score-based Simplified Version

## Model Overview
This model uses the CD Score of 5 methylation markers to directly assess colorectal cancer risk, without calculating probability values.

## CD Score Calculation
CD Score = Intercept + (Marker1 × Weight1) + (Marker2 × Weight2) + (Marker3 × Weight3) + (Marker4 × Weight4) + (Marker5 × Weight5)

## Diagnostic Rules
- **CD Score ≥ 0.5**: Positive (High Risk)
- **CD Score < 0.5**: Negative (Low Risk)

## Threshold Explanation
The diagnostic threshold is set to 0.5, an empirical value that can be adjusted based on clinical needs.

## Model Files
- `model_coefficients.csv`: Complete model coefficients (including intercept)
- `cd_score_parameters.csv`: CD Score calculation weights and threshold
- `model_metadata.csv`: Model metadata
- `crc_diagnosis_by_cd_score.R`: Diagnostic algorithm script
- `model_parameters.json`: JSON format parameters

## Usage Instructions
See the diagnostic algorithm script for detailed function descriptions and usage examples.
'

writeLines(documentation, "README.md")

# 8. Completion information --------------------------------------------------------------

cat("\n=== CD Score-based Diagnostic Model Files Generation Completed ===\n")
cat("Generated the following files:\n")
cat("1. model_coefficients.csv - Model coefficients (including intercept)\n")
cat("2. cd_score_parameters.csv - CD Score parameters and threshold\n")
cat("3. model_metadata.csv - Model metadata\n")
cat("4. crc_diagnosis_by_cd_score.R - Diagnostic algorithm script\n")
cat("5. model_parameters.json - JSON format parameters\n")
cat("6. README.md - Usage documentation\n")
cat("\nCD Score-based simplified diagnostic model is ready for use.\n")
