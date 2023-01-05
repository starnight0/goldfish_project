cc.cc<-read.table(file="out1.sclust.mafft.dS.CYPCAR.CYPCAR.txt", stringsAsFactors=F,header=T, comment.char="");
cc.zf<-read.table(file="out1.sclust.mafft.dS.CYPCAR.ENSDAR.txt", stringsAsFactors=F,header=T, comment.char="");
cc.gf<-read.table(file="out1.sclust.mafft.dS.CYPCAR.carAur.txt", stringsAsFactors=F,header=T, comment.char="");
gf.gf<-read.table(file="out1.sclust.mafft.dS.carAur.carAur.txt", stringsAsFactors=F,header=T, comment.char="");
gc.gf<-read.table(file="out1.sclust.mafft.dS.CTEIDE.carAur.txt", stringsAsFactors=F,header=T, comment.char="");
gc.zf<-read.table(file="out1.sclust.mafft.dS.CTEIDE.ENSDAR.txt", stringsAsFactors=F,header=T, comment.char="");
zf.gf<-read.table(file="out1.sclust.mafft.dS.ENSDAR.carAur.txt", stringsAsFactors=F,header=T, comment.char="");
cf.gf<-read.table(file="out1.sclust.mafft.dS.ENSAMX.carAur.txt", stringsAsFactors=F,header=T, comment.char="");
cf.zf<-read.table(file="out1.sclust.mafft.dS.ENSAMX.ENSDAR.txt", stringsAsFactors=F,header=T, comment.char="");
cc.cf<-read.table(file="out1.sclust.mafft.dS.CYPCAR.ENSAMX.txt", stringsAsFactors=F,header=T, comment.char="");
m=2;

cc.gf.h<-hist(cc.gf[cc.gf[,'dS']<m,'dS'], breaks=100*m, plot=F)
cc.cc.h<-hist(cc.cc[cc.cc[,'dS']<m,'dS'], breaks=100*m, plot=F)
gf.gf.h<-hist(gf.gf[gf.gf[,'dS']<m,'dS'], breaks=100*m, plot=F)
gc.zf.h<-hist(gc.zf[gc.zf[,'dS']<m,'dS'], breaks=100*m, plot=F)
gc.gf.h<-hist(gc.gf[gc.gf[,'dS']<m,'dS'], breaks=100*m, plot=F)
cf.gf.h<-hist(cf.gf[cf.gf[,'dS']<m,'dS'], breaks=100*m, plot=F)
zf.gf.h<-hist(zf.gf[zf.gf[,'dS']<m,'dS'], breaks=100*m, plot=F)
cc.zf.h<-hist(cc.zf[cc.zf[,'dS']<m,'dS'], breaks=100*m, plot=F)
cf.zf.h<-hist(cf.zf[cf.zf[,'dS']<m,'dS'], breaks=100*m, plot=F)
cc.cf.h<-hist(cc.cf[cc.cf[,'dS']<m,'dS'], breaks=100*m, plot=F)

gc.zf.y <- gc.zf.h$density/sum(gc.zf.h$density);
gc.gf.y <- gc.gf.h$density/sum(gc.gf.h$density);
cc.cf.y <- cc.cf.h$density/sum(cc.cf.h$density);
cc.zf.y <- cc.zf.h$density/sum(cc.zf.h$density);
cc.cc.y <- cc.cc.h$density/sum(cc.cc.h$density);
cc.gf.y <- cc.gf.h$density/sum(cc.gf.h$density);
gf.gf.y <- gf.gf.h$density/sum(gf.gf.h$density);
cf.gf.y <- cf.gf.h$density/sum(cf.gf.h$density);
cf.zf.y <- cf.zf.h$density/sum(cf.zf.h$density);
zf.gf.y <- zf.gf.h$density/sum(zf.gf.h$density);
y.max <- max( c(cc.gf.y, cc.cc.y, gf.gf.y, gc.gf.y, gc.zf.y, zf.gf.y, cf.zf.y) );

gc.zf.max.i <- max.col(t(gc.zf.y))
gc.gf.max.i <- max.col(t(gc.gf.y))
cc.cf.max.i <- max.col(t(cc.cf.y))
cc.zf.max.i <- max.col(t(cc.zf.y))
cc.cc.max.i <- max.col(t(cc.cc.y))
cc.gf.max.i <- max.col(t(cc.gf.y))
cf.gf.max.i <- max.col(t(cf.gf.y))
cf.zf.max.i <- max.col(t(cf.zf.y))
zf.gf.max.i <- max.col(t(zf.gf.y))
gf.gf.max.i <- max.col(t(gf.gf.y))

gc.zf.max.x <- gc.zf.h$mid[gc.zf.max.i];
gc.gf.max.x <- gc.gf.h$mid[gc.gf.max.i];
cc.cf.max.x <- cc.cf.h$mid[cc.cf.max.i];
cc.zf.max.x <- cc.zf.h$mid[cc.zf.max.i];
cc.cc.max.x <- cc.cc.h$mid[cc.cc.max.i];
cc.gf.max.x <- cc.gf.h$mid[cc.gf.max.i];
cf.gf.max.x <- cf.gf.h$mid[cf.gf.max.i];
cf.zf.max.x <- cf.zf.h$mid[cf.zf.max.i];
zf.gf.max.x <- zf.gf.h$mid[zf.gf.max.i];
gf.gf.max.x <- gf.gf.h$mid[gf.gf.max.i];

gc.zf.max.y <- gc.zf.y[gc.zf.max.i];
gc.gf.max.y <- gc.gf.y[gc.gf.max.i];
cc.cf.max.y <- cc.cf.y[cc.cf.max.i];
cc.zf.max.y <- cc.zf.y[cc.zf.max.i];
cc.cc.max.y <- cc.cc.y[cc.cc.max.i];
cc.gf.max.y <- cc.gf.y[cc.gf.max.i];
cf.gf.max.y <- cf.gf.y[cf.gf.max.i];
cf.zf.max.y <- cf.zf.y[cf.zf.max.i];
zf.gf.max.y <- zf.gf.y[zf.gf.max.i];
gf.gf.max.y <- gf.gf.y[gf.gf.max.i];

max.y <- list(cc.gf=cc.gf.max.y, cc.cc=cc.cc.max.y, gf.gf=gf.gf.max.y, gc.gf=gc.gf.max.y, gc.zf=gc.zf.max.y, zf.gf=zf.gf.max.y, cf.zf=cf.zf.max.y)
max.x <- list(cc.gf=cc.gf.max.x, cc.cc=cc.cc.max.x, gf.gf=gf.gf.max.x, gc.gf=gc.gf.max.x, gc.zf=gc.zf.max.x, zf.gf=zf.gf.max.x, cf.zf=cf.zf.max.x)
h <- list(cc.gf=cc.gf.h, cc.cc=cc.cc.h, gf.gf=gf.gf.h, gc.gf=gc.gf.h, gc.zf=gc.zf.h, zf.gf=zf.gf.h, cf.zf=cf.zf.h)
y <- list(cc.gf=cc.gf.y, cc.cc=cc.cc.y, gf.gf=gf.gf.y, gc.gf=gc.gf.y, gc.zf=gc.zf.y, zf.gf=zf.gf.y, cf.zf=cf.zf.y)
col <- c(rgb(1,0,0),  rgb(0,0.8,0.1),  rgb(0,0,1),  rgb(0.5,0.2,0.7),  rgb(1,1,0.2),  rgb(0,0.8,0.9),  rgb(1,0,0.8))
alpha=0.3
cola <- c(rgb(1,0,0,alpha),  rgb(0,0.8,0.1,alpha),  rgb(0,0,1,alpha),  rgb(0.5,0.2,0.7,alpha),  rgb(1,1,0.2,alpha),  rgb(0,0.8,0.9,alpha),  rgb(1,0,0.8,alpha))
text <- c("CC vs. GF", "CC vs. CC", "GF vs. GF", "GC vs. GF", "GC vs. ZF", "Zf vs. GF", "CF vs. ZF")

tiff("dS_distribution.tiff", w=3000,h=2000,res=300,compress="lzw+p");
plot(h[[1]]$mid, y[[1]], type='l', col=col[1], lwd=2, ylim=c(0,y.max*1.1), xlab="dS", ylab="Percentage", main="dS distribution", axes=F)
axis(1,lwd=2);
axis(2,lwd=2);
legend(x=1.5, y=y.max, legend=text[1:6], text.col=col[1:6], fill=col[1:6],border=F,bg="lightgray", box.lwd=0);
for (i in 1:6) {
	lines(h[[i]]$mid, y[[i]], type='l', col=col[i], lwd=2);
	polygon(c(0, h[[i]]$mid, m), c(0, y[[i]], 0), col = cola[i], border = NA);
}
for (i in 1:6) {
	points(max.x[[i]], max.y[[i]], col=col[i], pch=1);
	text(max.x[[i]], max.y[[i]]+y.max*0.02, label=max.x[[i]], col=col[i]);
}
dev.off()

