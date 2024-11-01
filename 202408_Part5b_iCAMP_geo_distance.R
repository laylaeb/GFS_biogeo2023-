iCAMP.dist <-read.csv("iCAMP_geo_distance.csv", head=T)

library(ggplot2)
p <-ggplot() + geom_point(data=iCAMP.dist, aes(geo.dist, HoS*100), col=rgb(0,0,1,0.25), cex=2)+geom_smooth(data=iCAMP.dist, aes(geo.dist, HoS*100),method="lm", col="blue")+theme_bw()+scale_x_continuous(trans='log10')+labs(x="geographic distance [m]", y="contribution [%]")
p <-p + geom_point(data=iCAMP.dist, aes(geo.dist, DL*100), col=rgb(1,0,0,0.25), cex=2)+geom_smooth(data=iCAMP.dist, aes(geo.dist, DL*100),method="lm", col="red")
p
