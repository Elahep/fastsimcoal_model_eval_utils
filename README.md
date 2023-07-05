# fastsimcoal_model_eval_utils
Evaluating fastsimcoal models by calculating AIC, delta AIC, Akaike weights, and other statistical criteria in R.


### 1- Boxplots of model likelihoods
After running several repeats per model, we can visually evaluate variations in MaxEstLhood per repeat using boxplots or point plots. It is best to have non-overlapping variation across the models.

```
library(ggplot2)

# Read the tables into a list
tables <- list(
  model1 = read.table("1allbestruns.txt", header = TRUE),
  model2 = read.table("2allbestruns.txt", header = TRUE),
  model3 = read.table("3allbestruns.txt", header = TRUE),
  model4 = read.table("4allbestruns.txt", header = TRUE),
  model5 = read.table("5allbestruns.txt", header = TRUE),
  model6 = read.table("6allbestruns.txt", header = TRUE)
)

# Combine the data into a single data frame
data <- do.call(rbind, tables)

# Create a factor variable for model names
data$model <- rep(names(tables), sapply(tables, nrow))

# Plotting the box plots
ggplot(data, aes(x = model, y = MaxEstLhood, fill = model)) +
  geom_boxplot() +
  labs(x = "Model", y = "MaxEstLhood") +
  ggtitle("Variations in MaxEstLhood by Model") +
  theme_minimal() +
  theme(legend.position = "none")


## since boxplot doesn't show much variation (which is good! We don't
# want overlapping MaxEstLhood between our models), we'll use
# another graph to visually explore the tiny variations
# using a point plot to show the individual observations for each model. 
# This can help reveal any small variations in the data.

# Plotting the point plot
final <- ggplot(data, aes(x = model, y = MaxEstLhood, color = model)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5) +
  labs(x = "Model", y = "MaxEstLhood") +
  ggtitle("Variations in MaxEstLhood by Model") +
  theme_minimal() +
  theme(legend.position = "none")

final
ggsave("variations in MaxEstLhood per model.pdf", final, width = 8, height = 6, units = "in")
````

### 2- Calculating AIC and other statistics
Using the minimum MaxEstLhood of each model, we can calculate AIC to find the best model. The model with the smallest AIC will be chosen as the best fitting model.

```
# Define the number of parameters (k) for each model, looking at the .bestlhoods output file.
k <- 15

# Create a list of data frames
model1 = read.table("1allbestruns.txt", header = TRUE)
model2 = read.table("2allbestruns.txt", header = TRUE)
model3 = read.table("3allbestruns.txt", header = TRUE)
model4 = read.table("4allbestruns.txt", header = TRUE)
model5 = read.table("5allbestruns.txt", header = TRUE)
model6 = read.table("6allbestruns.txt", header = TRUE)

## calculate AIC per model for all runs and choose 
# the smallest from all runs

AIC1 <- min(2 * k - 2 * (model1$MaxEstLhood / log10(exp(1))))
AIC2 <- min(2 * k - 2 * (model2$MaxEstLhood / log10(exp(1))))
AIC3 <- min(2 * k - 2 * (model3$MaxEstLhood / log10(exp(1))))
AIC4 <- min(2 * k - 2 * (model4$MaxEstLhood / log10(exp(1))))
AIC5 <- min(2 * k - 2 * (model5$MaxEstLhood / log10(exp(1))))
AIC6 <- min(2 * k - 2 * (model6$MaxEstLhood / log10(exp(1))))

## also to create a table:

AIC_values <- c(
  min(2 * k - 2 * (model1$MaxEstLhood / log10(exp(1)))),
  min(2 * k - 2 * (model2$MaxEstLhood / log10(exp(1)))),
  min(2 * k - 2 * (model3$MaxEstLhood / log10(exp(1)))),
  min(2 * k - 2 * (model4$MaxEstLhood / log10(exp(1)))),
  min(2 * k - 2 * (model5$MaxEstLhood / log10(exp(1)))),
  min(2 * k - 2 * (model6$MaxEstLhood / log10(exp(1))))
)

AIC_results <- data.frame(Model = paste("Model", 1:6), AIC = AIC_values)

## also calculate delta likelihood by subtracting Max EstLhood and ObsLhood
## and create a table

MaxObsLhood <- -1338.704 #this is constant across models

DeltaLhood <- c(
  deltaL1 <- MaxObsLhood-min(model1$MaxEstLhood),
  deltaL2 <- MaxObsLhood-min(model2$MaxEstLhood),
  deltaL3 <- MaxObsLhood-min(model3$MaxEstLhood),
  deltaL4 <- MaxObsLhood-min(model4$MaxEstLhood),
  deltaL5 <- MaxObsLhood-min(model5$MaxEstLhood),
  deltaL6 <- MaxObsLhood-min(model6$MaxEstLhood))

DeltaLhood_results <- data.frame(Model = paste("Model", 1:6), DeltaLhood = DeltaLhood)

## now we can calculate AIC weights:
library(qpcR)

w = akaike.weights(x = AIC_values)
w$deltaAIC  ##delta AIC values
w$rel.LL  ##relative likelihoods
w$weights  ##the Akaike weights

deltaAIC_results <- data.frame(Model = paste("Model", 1:6), deltaAIC = w$deltaAIC)
relLL_results <- data.frame(Model = paste("Model", 1:6), relative_Lhoods = w$rel.LL)
Akaikeweights_results <- data.frame(Model = paste("Model", 1:6), Akaike_weights = w$weights)

## The ΔAIC values indicate the difference in information loss 
#between each model and the best model. Lower ΔAIC values indicate
#better model fit, with a difference of 2 or less being considered 
#weak evidence against the model with the lowest AIC. 
#Larger differences, such as ΔAIC > 4 or 7, provide stronger evidence
#against the model with higher ΔAIC. So the larger the difference, it
#means the better our selected model is.

##Here are some general guidelines:
#ΔAIC < 2: There is very weak evidence against the model compared to the best model. The model can be considered competitive and may provide a similar fit to the data.
#2 ≤ ΔAIC < 4: There is weak evidence against the model compared to the best model. The model may still be plausible, but it is less supported by the data.
#4 ≤ ΔAIC < 7: There is substantial evidence against the model compared to the best model. The model is less likely to be the best-fitting model.
#ΔAIC ≥ 7: There is strong evidence against the model compared to the best model. The model is highly unlikely to be the best-fitting model.

##Now we can create a table with all the necessary information
#we need to report in the manuscript.

#put all maxEstLogLhoods in a table:

maxEstlogLhoods <- c(
  min(model1$MaxEstLhood),
  min(model2$MaxEstLhood),
  min(model3$MaxEstLhood),
  min(model4$MaxEstLhood),
  min(model5$MaxEstLhood),
  min(model6$MaxEstLhood)
)

maxEstLogLhoods_results <- data.frame(Model = paste("Model", 1:6), maxEstlogLhoods = maxEstlogLhoods)


final <- cbind(maxEstLogLhoods_results,DeltaLhood_results$DeltaLhood,
               AIC_results$AIC,deltaAIC_results$deltaAIC,
               Akaikeweights_results$Akaike_weights,relLL_results$relative_Lhoods)

colnames(final) <- c("Model","maxEstLogLhoods",
                     names(DeltaLhood_results)[2],names(AIC_results)[2], 
                     names(deltaAIC_results)[2], names(Akaikeweights_results)[2],
                     names(relLL_results)[2])
write.table(final, "AIC_results.txt", sep = "\t")
```
