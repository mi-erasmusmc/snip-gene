require(caret)
require(doMC)
registerDoMC(cores = 3)

removalColumns = c("EKP_name", "EKP_ID", "Name", "Gene.start..bp.", "Gene.end..bp.", "Chromosome.scaffold.name")
input = features[,-which(colnames(features) %in% removalColumns)]

model = train(Class ~ ., data = input, method = "ranger", metric = "ROC",
              trControl=trainControl(method = "cv", number = 10,
                                     classProbs = TRUE, savePredictions = TRUE, 
                                     summaryFunction = twoClassSummary, sampling = "down"))
