
library("hydroGOF")
#library("zoo")

setwd("/data07/shh/Stars4Water/hbv/cali")

gauge <- c("12.70","12.171","12.178","12.188","12.192","12.193","12.197","12.207","12.209","12.212","12.215")   # sorted based on the gauge name not gauge id
version <- c(1,1,1,1,1,1,0,0,1,1,1)
gauge2 <- c("01200070","01200171","01200178","01200188","01200192","01200193","01200197","01200207","01200209","01200212","01200215")
Lperarea <- c(17.43,18.27,22.69,17.99,23.52,16.23,28.68,18.16,21.87,17.03,31.47)
w_Lperarea <- Lperarea/sum(Lperarea)
area <- c(568.5,79.4,311,4.67,74.77,51.54,184.7,269.88,554.1,11.02,119.7)
w_area <- area/sum(area)



for (i in 1:length(gauge)) {        

obs_file <-  paste("/data07/shh/Stars4Water/data/Q/",gauge[i],".0.1001.",version[i],sep="")

obs_o <- read.table(obs_file,  header=F,  na.strings = "-9999.000000")
obs_time <- as.Date(obs_o[,1], format="%Y%m%d/%H%M")
obs_year <- as.POSIXlt(obs_time)$year+1900

sim_o <- read.table(paste("results/hbv_",gauge2[i],".var",sep=""), header=F)
sim_time <-   strptime(sim_o[,1],"%Y%m%d/%H%M")
sim_year <-  as.POSIXlt(sim_time)$year+1900

syear <- min(sim_year)
eyear <- max(sim_year)


obs <- obs_o[which(obs_year>=syear&obs_year<=eyear),2]
sim <- sim_o[which(sim_year>=syear&sim_year<=eyear),2]

select <- which(is.na(obs)==FALSE&obs>0&sim>0)

 NSE <- round(1 - sum((obs[select]-sim[select])^2)/sum((obs[select]-mean(obs[select]))^2),2)
 kge <- KGE(sim[select],obs[select],method="2012")

bias <- (mean(sim[select])-mean(obs[select]))/mean(obs[select])
 LKGE <- KGE(log(sim[select]),log(obs[select]),method="2012")
 LNSE <- round(1 - sum((log(obs[select])-log(sim[select]))^2)/sum((log(obs[select])-mean(log(obs[select])))^2),2)

boxcoxKGE <- KGE((sim[select]^0.3-(0.01*mean(sim[select]))^0.3)/0.3,(obs[select]^0.3-(0.01*mean(obs[select]))^0.3)/0.3,method="2012")


#round(1 - sum((log(obs[select])-log(sim[select]))^2)/sum((log(obs[select])-mean(log(obs[select])))^2),2)

#### error term based
RMOV <- abs(mean(sim[select]-obs[select]))
RMS <- mean(abs(sim[select]-obs[select]))

highobsQ <- quantile(obs,0.75,na.rm=T)
select <- which(is.na(obs)==FALSE&obs>highobsQ&sim>0)
RMS_H1 <- mean(abs(sim[select]-obs[select]))

lowobsQ <- quantile(obs,0.25,na.rm=T)
select <- which(is.na(obs)==FALSE&obs<lowobsQ&obs>0&sim>0)
RMS_L1 <- mean(abs(sim[select]-obs[select]))
 #NSE <- RMOV +RMS+RMS_H1+RMS_L1      #### ?? RMS_H1>>RMS_L1 or only use RMS, the defaul OF of PEST


if(i==1) {
out <- c(kge,NSE,LNSE,bias, LKGE,RMOV,RMS,RMS_H1,RMS_L1,boxcoxKGE)
} else {
out <- rbind(out, c(kge,NSE,LNSE,bias, LKGE,RMOV,RMS,RMS_H1,RMS_L1,boxcoxKGE))
}
} # end loop i


out_all <- apply(abs(out),2,FUN = function(x) weighted.mean(x, w=w_area,na.rm = TRUE))
out_m <- data.frame(gauge=gauge,kge=out[,1],nse=out[,2],lnse=out[,3],bias=out[,4],lkge=out[,5],rmov=out[,6],rms=out[,7],rms_h1=out[,8],rms_l1=out[,9],boxkge=out[,10])

write(c("NSE",out_all), file = "nse.txt",
      ncolumns = 1,
      append = F, sep = " ")
write.table(out_m,"nse_all.txt",col.names=T,row.names=F)

