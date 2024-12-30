library(seasonal)
library(x13binary)
library(readxl)
library(zoo)
library(xts)

checkX13()
df_data=read_excel("C:\\Users\\user\\Desktop\\dataind_seasadj.csv")
#df_rail<-na.omit(df_rail)
df_data_ts<-ts(df_data$original,start=c(2008,1),end=c(2022,7),frequency=12)
susu
df_out<-xts(mod$data, order.by = as.Date(df_data$time, format='%Y:%m'))
write.zoo(df_out,"R_seasadj_iip.csv",sep=",",col.names=TRUE)

plot(mod)
plot(resid(mod), main = "Residuals")
qqnorm(resid(mod), main = "Residuals compared to Normal")
pacf(resid(mod), "Partial ACF of residuals")