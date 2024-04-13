
library(GEOquery)
okunan1=getGEO("GSE50161",GSEMatrix = TRUE)
class(okunan1)

eset1=okunan1[[1]]
dim(eset1)
kopya=eset1
kopya=(exprs(kopya))

colnames(pData(eset1))
durum=factor(pData(eset1)$characteristics_ch1)

levels(durum)[levels(durum)=="tissue: normal brain (cerebellum)"]="normal brain"
levels(durum)[levels(durum)=="tissue: normal brain (frontal lobe of cerebral cortex)"]="normal brain"





annotation(eset1)<-"hgu133plus2.db"
library(genefilter)
filtrelenmis=nsFilter(eset1,
                      var.cutoff=0.90)
dim(filtrelenmis$eset)

sonveri=data.frame(t(exprs(filtrelenmis$eset)))





set.seed(111)
sD1nD1r=floor(.80*length(durum)) 
print(sD1nD1r) 

ind=sample.int(n=length(durum),size=sD1nD1r,replace=F)
ind

veriegitim=sonveri[ind,]
sD1nD1fegitim=durum[ind]

veritest=sonveri[-ind,]
sD1nD1ftest=durum[-ind]

splitIndex <- createDataPartition(durum, p = 0.7, list = FALSE)
train_data <- sonveri[splitIndex, ]
test_data <- sonveri[-splitIndex, ]
train_labels <- durum[splitIndex]
test_labels <- durum[-splitIndex]

library(e1071)  
library(caret)

data_frame1 <- as.data.frame(sonveri)
data_frame1$durum=durum
dim(data_frame1)
data_frame1$durum <- as.factor(data_frame1$durum)
model <- train(durum ~ ., data = data_frame1, method = "rf")
feature_importance <- varImp(model)
print(feature_importance)

model <- train(x = train_data, y = train_labels, method = "rf")
predictions <- predict(model, newdata = test_data)
conf_matrix <- confusionMatrix(predictions, test_labels)
accuracy <- conf_matrix$overall["Accuracy"]
kappa <- conf_matrix$overall["Kappa"]
cat("Do??ruluk (Accuracy):", accuracy, "\n")
cat("Kappa De??eri:", kappa, "\n")





data_frame1 <- as.data.frame(sonveri)
data_frame1$durum=durum
dim(data_frame1)
data_frame1$durum <- as.factor(data_frame1$durum)
model <- train(durum ~ ., data = data_frame1, method = "svmLinear")
feature_importance <- varImp(model)
print(feature_importance)
model <- train(x = train_data, y = train_labels, method = "svmLinear")
predictions <- predict(model, newdata = test_data)
conf_matrix <- confusionMatrix(predictions, test_labels)
accuracy <- conf_matrix$overall["Accuracy"]
kappa <- conf_matrix$overall["Kappa"]
cat("Do??ruluk (Accuracy):", accuracy, "\n")
cat("Kappa De??eri:", kappa, "\n")


data_frame1 <- as.data.frame(sonveri)
data_frame1$durum=durum
dim(data_frame1)
data_frame1$durum <- as.factor(data_frame1$durum)
model <- train(durum ~ ., data = data_frame1, method = "nb")
feature_importance <- varImp(model)
print(feature_importance)

model <- train(x = train_data, y = train_labels, method = "nb")
predictions <- predict(model, newdata = test_data)
conf_matrix <- confusionMatrix(predictions, test_labels)
accuracy <- conf_matrix$overall["Accuracy"]
kappa <- conf_matrix$overall["Kappa"]
cat("Do??ruluk (Accuracy):", accuracy, "\n")
cat("Kappa De??eri:", kappa, "\n")










