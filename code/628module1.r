
remove(list = ls())
bodyfat <- read.csv("BodyFat.csv")
#Units change
bodyfat$WEIGHT <- 0.4536*bodyfat$WEIGHT
bodyfat$HEIGHT <- 2.54*bodyfat$HEIGHT
head(bodyfat,n=5) 

#####Data Cleaning####
par(mfrow=c(3,2)) 
par(mgp=c(1.8,.5,0), mar=c(3,3,1,1)) #"Beautifies" plots when creating multiple figures. Google this for more info.
hist(bodyfat$BODYFAT,breaks=20,
     main="Histogram of Body Fat %",xlab="Body Fat %")
hist(bodyfat$DENSITY,breaks=20,
     main="Histogram of Density",xlab="Density (gm/cm^3)")
hist(bodyfat$AGE,breaks=20,
     main="Histogram of Age",xlab="Age (yrs)")
hist(bodyfat$WEIGHT,breaks=30,
     main="Histogram of Weight",xlab="Weight (kg)")
hist(bodyfat$HEIGHT,breaks=30,
     main="Histogram of Height",xlab="Height (cm)")

bodyfat[bodyfat$BODYFAT>=45|bodyfat$BODYFAT<3,]

bodyfat[bodyfat$DENSITY<1|bodyfat$DENSITY>1.1,]

bodyfat[bodyfat$AGE>80,]

bodyfat[bodyfat$WEIGHT>150,]

bodyfat[bodyfat$HEIGHT<140,]

siri_data <- data.frame(ID=bodyfat$IDNO,BODYFAT=bodyfat$BODYFAT,DENSITY_INVERSE=1/bodyfat$DENSITY)
plot(BODYFAT~DENSITY_INVERSE,data=siri_data, pch = 16, cex = 1.3, col = "blue", main = "Body fat percentage against density", xlab = "1/Density (cm^3/gm)", ylab = "Body fat percentage (%)")
abline(lm(BODYFAT~DENSITY_INVERSE,data=siri_data))

siri.lm = lm(BODYFAT~DENSITY_INVERSE,data=siri_data) 
siri.res = resid(siri.lm)
bodyfat[abs(siri.res)>1.5,]

bodyfat<-bodyfat[-c(172,182,42,48,76,96),]
dim(bodyfat)

#### Data transformation####
bodyfat$NECK3 <- bodyfat$ADIPOSITY^3
bodyfat$CHEST3 <- bodyfat$CHEST^3
bodyfat$ABDOMEN3 <- bodyfat$ABDOMEN^3
bodyfat$HIP3 <- bodyfat$HIP^3
bodyfat$THIGH3 <- bodyfat$THIGH^3
bodyfat$KNEE3 <- bodyfat$KNEE^3
bodyfat$ANKLE3 <- bodyfat$ANKLE^3
bodyfat$BICEPS3 <- bodyfat$BICEPS^3
bodyfat$FOREARM3 <- bodyfat$FOREARM^3
bodyfat$WRIST3 <- bodyfat$WRIST^3

#### Variable selection ####
bodyfat$volume <- bodyfat$WEIGHT/bodyfat$DENSITY
volume_data <- bodyfat
volume_data$IDNO<-NULL;volume_data$BODYFAT<-NULL;volume_data$DENSITY<-NULL
step_model <- step(lm(volume~.,data=volume_data),direction = "both",k=log(nrow(volume_data)),trace = 0)
summary(step_model)

predict(step_model,newdata=data.frame(WEIGHT=70,HIP=95,WRIST=17,ABDOMEN3=85^3,HIP3=95^3),interval="predict")

confint(step_model)

bodyfat_density <- data.frame(BODYFAT=bodyfat$BODYFAT,DENSITY_INVERSE=1/bodyfat$DENSITY)
density_model <- lm(BODYFAT~DENSITY_INVERSE,bodyfat_density)
summary(density_model)

predict(density_model,newdata=data.frame(DENSITY_INVERSE=1/1.07),interval="predict")

predict(step_model,newdata=data.frame(WEIGHT=70,HIP=95,WRIST=17,ABDOMEN3=85^3,HIP3=95^3),interval="predict")/70*457-414

predict_bodyfat <- function(traindata,testdata){
  bodyfat_density_data_train <- data.frame(BODYFAT=traindata$BODYFAT,DENSITY_INVERSE=1/traindata$DENSITY)
  bodyfat_density_model <- lm(BODYFAT~DENSITY_INVERSE,bodyfat_density_data_train)
  volume_data_train <- traindata
  volume_data_train$volume <- volume_data_train$WEIGHT/volume_data_train$DENSITY
  volume_data_train$IDNO <- NULL;volume_data_train$BODYFAT <- NULL;volume_data_train$DENSITY <- NULL
  volume_model <- lm(volume~WEIGHT+HIP+WRIST+ABDOMEN3+HIP3,volume_data_train)
  volume_data_test <- testdata
  predicted_volume <- predict(volume_model,volume_data_test)
  bodyfat_density_data_test <- data.frame(BODYFAT=testdata$BODYFAT,DENSITY_INVERSE=predicted_volume/testdata$WEIGHT)
  predicted_bodyfat <- predict(bodyfat_density_model,bodyfat_density_data_test)
  return(predicted_bodyfat)
}

set.seed(2018)
MSE <- c()
for(i in 1:10){
  test_index <- sample(1:nrow(bodyfat),floor(nrow(bodyfat)/10))
  testdata <- bodyfat[test_index,];traindata <- bodyfat[-test_index,]
  predicted.values <- predict_bodyfat(traindata,testdata)
  MSE <- c(MSE,sum((predicted.values-testdata$BODYFAT)^2)/floor(nrow(bodyfat)/10))
}
mean(MSE)

plot(predict(step_model),rstandard(step_model),pch=23,bg="red",cex=1.2,
    xlab="Predicted Volume(dm^3)", ylab="Standardized Residuals",main="Standardized Residual Plot")
abline(a=0,b=0,col="black",lwd=3)

bodyfat[predict(step_model)>160,c(1,ncol(bodyfat),2:17)]

pii = hatvalues(step_model)
cooki = cooks.distance(step_model)
par(mfrow = c(2,1))
n = dim(bodyfat)[1]
plot(1:n,pii,type="p",pch=23,bg="red",cex=1.2,
     xlab="Index (Each Observation)",ylab="Pii",main="Influence Values (Pii)")
plot(1:n,cooki,type="p",pch=23,bg="red",cex=1.2,
     xlab="Index (Each Observation)",ylab="Cook's Distance",main="Influence Values (Cook's Distance)")

pii[39]

bodyfat[cooki>0.08,c(1,ncol(bodyfat),2:17)]

qqnorm(rstandard(step_model),pch=23,bg="red",cex=1.2)
abline(a=0,b=1,col="black",lwd=3)

additive <- data.frame(Volume=bodyfat$volume,Residual=rstandard(step_model),Weight=bodyfat$WEIGHT,Hip=bodyfat$HIP,
                      Wrist=bodyfat$WRIST,Abdomen3=bodyfat$ABDOMEN3,Hip3=bodyfat$HIP3)
pairs(additive)
