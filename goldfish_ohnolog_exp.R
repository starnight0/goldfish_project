library('gplots');
library('tsne');
library('entropy');

###############################
# Sub Functions
###############################
#{{{
plot_tree <- function(h, lab, col, lab.col, cex=1, lab.cex=1, pch=NA, dLeaf=NULL, horiz=F )
{
    h1d <- as.dendrogram(h)
    local( {
        lab1 <<- function(x) {
            if (is.leaf(x)) {
                i <<- i+1

                attr(x, "label") <- lab[h$order[i]]
                attr(x, "nodePar") <- c(attr(x, "nodePar"), list(col=col[h$order[i]], lab.col=lab.col[h$order[i]], lab.cex=lab.cex,  cex=cex, pch=pch) ); 
            }
            x
        }
        i <- 0
    } )
    h1da <- dendrapply(h1d, lab1);
    plot(h1da, dLeaf=dLeaf, horiz=horiz);
    return(h1da);
}

cor_dist <- function(x) { as.dist(1-cor(t(x))) }
cor_dist1 <- function(x) { as.dist(1-cor(t(x), method="spearman")) }

czl_cutree <- function(h, min_height, min_size)
{
#    hc <- cutree(h,h=min_height);
    h.n <- length(h$order);
#    hc.size <- tapply(rep(1,h.n), hc, sum);
    size1<-rep(0,h.n-1);
    is_root <- rep(0, h.n-1);
    parent0 <- rep(0,h.n);
    parent <- rep(0,h.n-1);
    id=0;
    for (i in 1:(h.n-1)) {
        i1 = h$merge[i,1]
        i2 = h$merge[i,2]
        if (i1<0) {
            parent0[-i1] = i;
            if (i2<0) { size1[i]=2; parent0[-i2]=i;} else { size1[i]=size1[i2]+1; parent[i2]=i;}
        } else {
            parent[i1] = i;
            if (i2<0) { size1[i]=size1[i1]+1; parent0[-i2]=i; } else { size1[i]=size1[i1]+size1[i2]; parent[i2]=i; }
        }
        if (h$height[i]>=min_height && size1[i]>=min_size) {
            if (h$height[i1]<min_height && size1[i1]>=min_size && h$height[i2]<min_height && size1[i2]>=min_size) {
                id=id+1; is_root[i1] = id;
                id=id+1; is_root[i2] = id;
            } else if (h$height[i1]>min_height && size1[i1]<min_size || h$height[i2]>min_height && size1[i2]>=min_size) {
                if (is_root[i1]==0 && is_root[i2]==0) { id=id+1; is_root[i]=id; }
                else {
                    if (is_root[i1]==0) { id=id+1; is_root[i1]=id; }
                    if (is_root[i2]==0) { id=id+1; is_root[i2]=id; }
                }
            }
        }
    }
    hc <- rep(0,h.n);
    for (i in 1:h.n) {
        j=parent0[i];
        while (is_root[j]==0) { j=parent[j]; }
        hc[i]=is_root[j];
    }
    return(hc);
}

czl_cutree2 <- function(h, min_size)
{
    h.n <- length(h$order);
    size1<-rep(0,h.n-1);
    is_root <- rep(0, h.n-1);
    parent0 <- rep(0,h.n);
    parent <- rep(0,h.n-1);
    id=0;
    for (i in 1:(h.n-1)) {
        i1 = h$merge[i,1]
        i2 = h$merge[i,2]
        if (i1<0) {
            parent0[-i1] = i;
            if (i2<0) { size1[i]=2; parent0[-i2]=i;} else { size1[i]=size1[i2]+1; parent[i2]=i;}
        } else {
            parent[i1] = i;
            if (i2<0) { size1[i]=size1[i1]+1; parent0[-i2]=i; } else { size1[i]=size1[i1]+size1[i2]; parent[i2]=i; }
        }
    }
    to_break<-c(h.n-1);
    while(length(to_break)>0) {
        to_break1<-c();
        for (ii in 1:length(to_break)) {
            i <- to_break[ii];  i1 <- h$merge[i,1];  i2 <- h$merge[i,2];
            if (size1[i1]<min_size || size1[i]<min_size*2) {
                if (size1[i2]<min_size || size1[i]<min_size*2) {
                    id=id+1; is_root[i]=id;
                } else {
                    id=id+1; is_root[i1]=id;
                    to_break1 <- c(to_break1, i2);
                }
            } else {
                if (size1[i2]<min_size) {
                    id=id+1; is_root[i2]=id;
                    to_break1 <- c(to_break1, i1);
                } else {
                    to_break1 <- c(to_break1, i1);
                    to_break1 <- c(to_break1, i2);
                }
            }
        }
        to_break<-to_break1;
    }
    hc <- rep(0,h.n);
    for (i in 1:h.n) {
        j=parent0[i];
        while (is_root[j]==0) { j=parent[j]; }
        hc[i]=is_root[j];
    }
    hc1<-hc[h$order];
    hc1d<-diff(hc1);
    hc1d[hc1d!=0]=1;
    hc1 <- cumsum(c(1,hc1d));
    hc[ h$order ] <- hc1;

    return(hc);
}

czl_get_leaf <- function(h, root)
{
    leaf<-c();
    for (i in 1:2) {
        if (h$merge[root,i]<0) {
            leaf = c(leaf,-h$merge[root,i]);
        } else {
            leaf = c(leaf, czl_get_leaf(h,h$merge[root,i]));
        }
    }
    return(leaf);
}

czl_hclust1<-function(dist)
{
    return(hclust(dist, method="complete"));
}

#}}}
################################################################

# red, green, blue, cyan, magenta, yellow, gray
color8 <- c("#FF0000", "#00FF00", "#0000FF", "#00FFFF", "#FF00FF", "#FFFF00", "#808080", "#FFC0CB");
color8a <- paste(color8, "60", sep="");
color <- c("#FF0000", "#00FF00", "#0000FF", "#00FFFF", "#FF00FF", 
		"#FFFF00", "#808080", "#FFC0CB", "#D00000", "#00D000", 
		"#000080", "#008080", "#800080", "#808000", "#7FFFE0", 
		"#FFE0D0", "#FFA010", "#A03030", "#DEB888", "#609AA0",
		"#80FF00");
# color <- rainbow(32);

##########################
# Load GO annotation
##########################
# {{{
id_to_gs0<-list();
id_to_gs<-list();
gs <- list();
for (f in c("GO2_MF", "GO2_CC", "GO2_BP")) {
    id_to_gs0[[ f ]] <- read.table(paste("~/data/goldfish/11549472/sergey_canu70x/arrow/carAur01/big/ips/GF/gene.ips.id_to_", f, ".txt", sep=""), head=F, col.names=c("gene_id", "annot_id", "annot_name", "level"), sep="\t", quote="");
    id_to_gs[[f]] <- split(id_to_gs0[[f]][,2], id_to_gs0[[f]][,1]);
    gs[[f]] <- split(id_to_gs0[[f]][,1], id_to_gs0[[f]][,2]);
}
for (f in c("Reactome", "KEGG") ) {
    id_to_gs0[[ f ]] <- read.table(paste("~/data/goldfish/11549472/sergey_canu70x/arrow/carAur01/big/ips/GF/gene.ips.id_to_", f, ".txt", sep=""), head=F, col.names=c("gene_id", "annot_id", "annot_name"), sep="\t", quote="");
    id_to_gs0[[f]] <- cbind(id_to_gs0[[f]], rep(0, dim(id_to_gs0[[f]])[1]));
    colnames(id_to_gs0[[f]]) = c("gene_id", "annot_id", "annot_name", "level");
    id_to_gs1[[f]][,4]=0;
    id_to_gs[[f]] <- split(id_to_gs0[[f]][,2], id_to_gs0[[f]][,1]);
    gs[[f]] <- split(id_to_gs0[[f]][,1], id_to_gs0[[f]][,2]);
}
gs_id_name <- list();
for (f in c("GO2_MF", "GO2_CC", "GO2_BP", "Reactome", "KEGG")) {
    gs_id_name[[f]] <- unique(id_to_gs0[[f]][,2:4]);
    rownames(gs_id_name[[f]]) <- gs_id_name[[f]][,1];
}
# }}}
########################

###########################
# Load gene annotation
###########################
bgp <- read.table("~/data/goldfish/11549472/sergey_canu70x/arrow/carAur03/big/carAur03.noM.gene.bgp", stringsAsFactor=F, head=F, sep="\t");
rownames(bgp) <- bgp[,18];
bgp.exon_len <- strsplit(bgp[,11], ",", fixed=T);
bgp.exon_len <- lapply(bgp.exon_len, as.numeric);
names(bgp.exon_len) <- rownames(bgp);
bgp.exon_num <- bgp[,10];
names(bgp.exon_num) <- rownames(bgp);
bgp.tran_len <- sapply(bgp.exon_len, sum);

#############################

tpm_thres=1;


#ohno_pair <- read.table('~/data/goldfish/11549472/sergey_canu70x/arrow/ohnolog/fish4.cluster/czl_ohno_syn.out3/rescue_m1.6.cluster.pair_extend.norm_coverage.txt', head=F);
ohno_pair <- read.table('~/data/goldfish/11549472/sergey_canu70x/arrow/ohnolog/fish4.cluster/czl_ohno_syn.out3/pair_from_chainnet.extended.txt', stringsAsFactor=F, head=F, sep="\t");
tpm <- read.table('~/data/goldfish/11549472/sergey_canu70x/arrow/RNA_star_run4/out/RSEM_out/mat/gene.tpm.ungroup.mat', stringsAsFactor=F);
ltpm <- log2(tpm[,2:(m+1)]+1);
#tpm21 <- read.table('~/data/goldfish/11549472/sergey_canu70x/arrow/RNA_star_run2/out/RSEM_out1/mat/gene.tpm.ungroup.mat')
#tpm22 <- read.table('~/data/goldfish/11549472/sergey_canu70x/arrow/RNA_star_run2/out/RSEM_out2/mat/gene.tpm.ungroup.mat')
n<-dim(ohno_pair)[1]
m<-dim(tpm)[2]-1;
colnames(ohno_pair) <- c('chr1','gid1','name1','chr2','gid2','name2', 'identity', 'score', 'cov1', 'cov2')
max1 <- apply(tpm[as.character(ohno_pair[,2]),2:(m+1)],1,max);
max2 <- apply(tpm[as.character(ohno_pair[,5]),2:(m+1)],1,max);
f <- (max1>=tpm_thres | max2>=tpm_thres) & bgp.exon_num[ohno_pair[,2]]>=2 & bgp.exon_num[ohno_pair[,5]]>=2  & bgp.tran_len[ohno_pair[,2]]>=90 & bgp.tran_len[ohno_pair[,5]]>=90
ohno_pair_f <- ohno_pair[f,]
n1 <- dim(ohno_pair_f)[1]
tpmA0 <- tpm[as.character(ohno_pair_f[,2]),];   tpmA <- tpmA0[,2:(m+1)]; ltpmA <- log2(tpmA+1);
tpmB0 <- tpm[as.character(ohno_pair_f[,5]),];   tpmB <- tpmB0[,2:(m+1)]; ltpmB <- log2(tpmB+1);
m <- dim(tpmA)[2];
sm_names <- colnames(tpmA);
sm_names <- colnames(tpmB);

a1 <- grep("Brain", colnames(tpmA));
a2 <- grep("Gill", colnames(tpmA));
tpm7A <- cbind( apply(tpmA[,a1],1,max), tpmA[,2], apply(tpmA[,a2],1,max), tpmA[,5], tpmA[,8:10]);
colnames(tpm7A) <- c("Brain", "Eye", "Gill", "Bone", "Heart", "Muscle", "TailFin");
ltpm7A <- log2(tpm7A+1);
tpm7B <- cbind( apply(tpmB[,a1],1,max), tpmB[,2], apply(tpmB[,a2],1,max), tpmB[,5], tpmB[,8:10]);
colnames(tpm7B) <- c("Brain", "Eye", "Gill", "Bone", "Heart", "Muscle", "TailFin");
ltpm7B <- log2(tpm7B+1);
names <- paste(rownames(tpmA), rownames(tpmB), sep="__"); 
gene_names<- as.character(tpmA0[,1]);
gene_names[ gene_names==rownames(tpmA) ] = names[ gene_names==rownames(tpmA) ];
gene_names[ gene_names=="_" ] = names[ gene_names=="_" ];
gene_names2 <- c( paste(gene_names, rownames(tpmA), sep="|A|"), paste(gene_names, rownames(tpmB), sep="|B|") );

exp_sm_count <- apply(tpm7A+tpm7B>=tpm_thres, 1, sum);
exp_sm_countA <- apply(tpm7A>=tpm_thres, 1, sum);
exp_sm_countB <- apply(tpm7B>=tpm_thres, 1, sum);
max_norm_ltpm_div <- apply( abs(ltpm7A-ltpm7B), 1, max );

# get tissue specific-expressed gene-pairs
f3<-(tpm7A>=4 & tpm7B<0.5) | (tpm7B>=4 & tpm7A<0.5);
specific_count <- apply(f3, 1, sum )
f4<-f3[ specific_count>0, ];
specific_sm<-rep("",length(specific_count));
for (i in 1:length(specific_count)) {
    if (sum(f3[i,])>0) { specific_sm[i] <- paste(sm_names[f3[i,]], sep=",", collapse=";") }
}
a<-cbind(rep(gene_names[specific_count>0],2), rbind(tpm7A[specific_count>0,], tpm7B[specific_count>0,]), as.character(rep(specific_sm[specific_count>0],2)));




##################################
# cluster A test
##################################
# {{{
n2<-50;
ltpm_cor_dist2 <- cor_dist1( as.matrix(rbind(ltpmA[1:n2,],ltpmB[1:n2,])) );
data2 <- as.matrix(rbind(ltpmA[1:n2,2:m],ltpmB[1:n2,2:m]));
data2.rowmean <- rowMeans(data2)
data2.rowsd <- apply(data2, 1, sd);
norm_data2 <- (data2-data2.rowmean)/data2.rowsd;
rownames2 <- rep(gene_names[1:n2],2);
h2 <- hclust( ltpm_cor_dist2,  method="complete" );
#h1 <- hclust( tpm_cor_dist,  method="ward.2D" );
max_z <- ceiling(max(max(norm_data2), -min(norm_data2)));
hc21 <- cutree(h2, h=0.5);
hc21.size <- tapply(rep(1,(n2*2)), hc21, sum);
hc24 <- czl_cutree2(h2, min_size=20)
hc24.size <- tapply(rep(1,(n2*2)), hc24, sum);
png("test.dendrogram.png", height=1536, width=1024, res=150);
plot_tree(h2, rownames2, hc24, lab.cex=0.5, horiz=T);
dev.off();

    d <- norm_data2[hc24==1,h1_sm$order];
    g <- rownames2[hc24==1];
    n4 <- sum(hc24==1);
    heatmap.2(d, labRow=g, distfun=cor_dist1, hclust=czl_hclust1, trace="none", scale="none", col=redblue(255), lwid=c(100+n4,40*m+200), lhei=c(200,20*n4+200), margins=c(9,9), dendrogram="row");
# }}}


##################################
# cluster
##################################
# {{{
f5<-apply(ltpm7A>=1,1,sum)>0 & apply(ltpm7B>=1,1,sum)>0;
tpm7A1  <- tpm7A[f5,];
ltpm7A1 <- ltpm7A[f5,];
tpm7B1  <- tpm7B[f5,];
ltpm7B1 <- ltpm7B[f5,];
data1 <- as.data.frame(rbind(ltpm7A1,ltpm7B1));
data1.n <- dim(data1)[1];
data1.m <- dim(data1)[2];
data1.n1 <- dim(data1)[1]/2;
data1.rowmean <- apply(data1, 1, mean);
data1.rowsd   <- apply(data1, 1, sd);
norm_data1_1 <- (data1-data1.rowmean);
norm_data1 <- (data1-data1.rowmean)/data1.rowsd;
rm(data1.rowmean);
rm(data1.rowsd);
data1.gene_name1 <- gene_names[f5];
data1.rownames1 <- c( paste(gene_names[f5], "1", sep="|"),paste(gene_names[f5], "2", sep="|") );
norm_data1.max_z <- ceiling(max(max(norm_data1), -min(norm_data1)));
data1.max_z <- ceiling(max(max(data1), -min(data1)));
norm_data1_1.max_z <- ceiling(max(max(norm_data1_1), -min(norm_data1_1)));
ltpm_cor_dist1 <- cor_dist1( norm_data1 );
ltpm_p1_dist   <- dist(data1, method="minkowski", p=1);

###### cluster by sample ######
ltpm_cor_dist_sm <- cor_dist( t(norm_data1) );
h1c_sm <- hclust( ltpm_cor_dist_sm,  method="ward.D2" );
ord_sm <- h1c_sm$order;
#ltpm_euc_dist_sm <- dist( t(data1) );
#h1e_sm <- hclust( ltpm_euc_dist_sm,  method="ward.D2" );
######


h<-list();
hc_v<- list();
same_v <- list();
h.methods <- c("complete", "average", "mcquitty", "ward.D2");
h.methods2 <- c("complete P", "average P", "mcquitty P", "ward.D2 P");
for (hi in 1:8) {
    if (hi<=4) { h[[hi]] <- hclust( ltpm_cor_dist,  method=h.methods[hi] ); }
    else  { h[[hi]] <- hclust( ltpm_cor_dist,  method=h.methods[hi-4] ); }
    hc_v[[hi]] <- matrix(0, 100, n1*2);
    same_v[[hi]] <- rep(0,100);
    for (k in 1:100) { 
        hc_v[[hi]][k,] <- cutree(h[[hi]], k=k);
        same_v[[hi]][k] = sum(hc_v[[hi]][k,1:n1]==hc_v[[hi]][k,n1+(1:n1)]);  
    }
}
par(mar=c(4,4,4,10)); plot(1:100, same_v[[1]], cex=0.5, col=1);
for (i in 2:8) { points(1:100, same_v[[i]], col=i, cex=0.5); }
legend(x=par()$usr[2], y=par()$usr[4], legend=c(h.methods,h.methods2), text.col=1:8, xpd=T, border=F);

#--------------------------------------
# clutering using euclidean distance
#--------------------------------------
# {{{
ltpm_euc_dist  <- dist(data1, method="euclidean");
h1e <- hclust(ltpm_euc_dist, method="ward.D2");
rm(ltpm_euc_dist);
v <- rep(0,100);
for (nc in 1:100) {
	cid<<-cutree(h1e, k=nc);  csz <- tapply(rep(1,(data1.n)), cid, sum);
	same <- cid[1:data1.n1]==cid[data1.n1+(1:data1.n1)];
	pair_same<-matrix(0, nc,nc);
	for (i in 1:data1.n1) {
		c1 <- cid[i];
		c2 <- cid[i+data1.n1];
		pair_same[c1,c2] <- pair_same[c1,c2] +1;
		pair_same[c2,c1] <- pair_same[c2,c1] +1;
	}
	v[nc] <- sum(diag(pair_same));
}
h1e.nc <- 20;
h1e.cid <- cutree(h1e, k=h1e.nc);  h1e.csz <- tapply(rep(1,(data1.n)), h1e.cid, sum);
h1e.same <- h1e.cid[1:data1.n1]==h1e.cid[data1.n1+(1:data1.n1)];
h1e.pair_same<-matrix(0, h1e.nc,h1e.nc);
h1e.cord<- unique(h1e.cid[h1e$order])
h1e.ord_csz <- h1e.csz[ h1e.cord ];
for (i in 1:data1.n1) {
    c1 <- h1e.cid[i];
    c2 <- h1e.cid[i+data1.n1];
    h1e.pair_same[c1,c2] <- h1e.pair_same[c1,c2] +1;
    h1e.pair_same[c2,c1] <- h1e.pair_same[c2,c1] +1;
}
h1e.pair_same_frac<-matrix(0, h1e.nc,h1e.nc);
for (i in 1:h1e.nc) {
for (j in 1:h1e.nc) {
    if (i==j) {
        h1e.pair_same_frac[i,j] <- h1e.pair_same[i,j]/h1e.csz[i];
    } else {
        h1e.pair_same_frac[i,j] <- h1e.pair_same[i,j]/(h1e.csz[i]-h1e.pair_same[i,i] + h1e.csz[j]-h1e.pair_same[j,j]);
    }
}
}

#---------------------------------
# heatmap
#---------------------------------
#{{{
data1.z <- max(apply(data1,2,function(x) { quantile(x,(0:50)/50)[50]; }))

local( {
w = 1024;
max_z <<-ceiling(data1.z);
d<-data1[h1e$order,ord_sm];
d[ which(d>max_z,arr.ind=T) ] = max_z; 
png(paste("~/tmp/log_tpm_euc_heatmap/K20/all.small.png", sep=""), width=w, height=w*5, res=150);
#par(mar=c(1,1,1,1)+0.1);
#layout(matrix(c(1,2),1,2), widths=c(1,1))
#par(mar=c(1,1,1,2)+0.1);
#plot_tree(h1e, rep(".",data1.n), color[h1e.cid], color[h1e.cid], lab.cex=4, horiz=T);
#par(mar=c(1,1,1,2)+0.1);
layout(matrix(c(1,2),2,1), height=c(1,19))
par(mar=c(1,1,0,1)+0.1);
image(seq(0,max_z,length.out=256), 0:1, matrix(seq(0,max_z,length.out=256),256,1), col=bluered(255), axes=F, zlim=c(0, max_z), xlab="", ylab="");
axis(1, at=seq(0,max_z,1) )
par(mar=c(0,1,2,1)+0.1);
image(0:data1.m, 0:data1.n, t(d), col=bluered(255), axes=F, zlim=c(0, max_z), xlab="", ylab="");
hline=cumsum(h1e.ord_csz);
for (i in 1:(h1e.nc-1)) {
	lines(c(-1,data1.m+1), c(hline[i], hline[i]), col='gray2', xpd=T);
}
#    image(0:data1.m, 0:n4, t(d), col=rgb((0:255)/255,0,0), axes=F, zlim=c(0, max_z), xlab="", ylab="");
dev.off();
} );

for (j in 1:h1e.nc) {
	i<- h1e.cord[j];
	max_z=data1.z;
    d <- data1[h1e$order,ord_sm] [h1e.cid[h1e$order]==i,];
	d[ which(d>max_z,arr.ind=T) ] = max_z; 
    g <- data1.rownames1[h1e$order][h1e.cid[h1e$order]==i];
    n4 <- sum(h1e.cid[h1e$order]==i);
    
#    png(paste("norm_log_tpm_heatmap/part", sprintf("%03d",i) , ".large.png", sep=""), width=100+n4+40*m+200, height=200+20*n4+200, res=150);
#    heatmap.2(d, labRow=g, distfun=cor_dist1, hclust=czl_hclust1, trace="none", scale="none", col=redblue(255), lwid=c(100+n4,40*m+200), lhei=c(200,20*n4+200), margins=c(9,9), dendrogram="row");
#    layout( matrix(1:6, 2,3), width=c(100+n4, 40+m, 200), height=c(20*n4,200))
#    image(0:(m-1), 0:n4, t(d), col=redblue(255), axes=F, zlim=c(-max_z, max_z), xlab="", ylab="");
#    dev.off();

    w = 800;
    png(paste("~/tmp/log_tpm_euc_heatmap/K20/C", sprintf("%03d",i) , ".small.png", sep=""), width=w, height=w*2, res=150);
    par(mar=c(0,0,0,0)+0.1);
    image(0:data1.m, 0:n4, t(d), col=bluered(255), axes=F, zlim=c(0, max_z), xlab="", ylab="");
#    image(0:data1.m, 0:n4, t(d), col=rgb((0:255)/255,0,0), axes=F, zlim=c(0, max_z), xlab="", ylab="");
	dev.off();

	w=40*data1.m;
    png(paste("~/tmp/log_tpm_euc_heatmap/K20/C", sprintf("%03d",i) , ".large.png", sep=""), width=w, height=8*n4, res=150);
    par(mar=c(0,0,0,0)+0.1);
    image(0:data1.m, 0:n4, t(d), col=bluered(255), axes=F, zlim=c(0, max_z), xlab="", ylab="");
#    image(0:data1.m, 0:n4, t(d), col=rgb((0:255)/255,0,0), axes=F, zlim=c(0, max_z), xlab="", ylab="");
	dev.off();

    png(paste("log_tpm_euc_heatmap/K20/C", sprintf("%03d",i) , ".small.png", sep=""), width=w, height=w*2, res=150);
    par(mar=c(10,3,3,3));
    image(0:data1.m, 0:n4, t(d), col=redblue(255), axes=F, zlim=c(-max_z, max_z), xlab="", ylab="");
    text(0:(m-2), rep(-5,m-1), sm_names[h1_sm$order], xpd=T, srt=-45, pos=4);
    dev.off();
}
#}}}


#---------------------------------
# shared gene pairs between clusters
#---------------------------------
local( {
par(mar=c(0,0,0,0)+0.1);
c<- unique(h1e.cid[h1e$order])
sz<-h1e.csz[c];
d<-h1e.pair_same_frac[c,c]
d[ which(d>0.4, arr.ind=T) ] = 0.4;
png(paste("~/tmp/dendrogram.euc.k",h1e.nc, ".pair_frac.png", sep=""), height=1024, width=1024, res=150);
x<-c(0,cumsum(sz));
image(x, x, t(d), col=bluered(255), zlim=c(0,0.4), axes=F, xlab="", ylab="");
for (i in 1:h1e.nc) {
for (j in 1:h1e.nc) {
	x1<-h1e.pair_same[c,c][i,j];
	x2<-h1e.pair_same_frac[c,c][i,j];
	if (x2>=0.05 && x1>=20) {
		if (sz[i]>data1.n/30 && sz[j]>data1.n/30) {
			if (i>j) { text(x[i]-0.5+sz[i]/2, x[j]-0.5+sz[j]/2, sprintf("%2.1f", x2*100)); }
			else { text(x[i]-0.5+sz[i]/2, x[j]-0.5+sz[j]/2, x1);}
		} else {
#			text(x[i]-0.5+sz[i]/2, x[j]-0.5+sz[j]/2, sprintf("%2.1f", x1*100), cex=0.3);
		}
	}
}
}
dev.off();
} );
par(mar=c(1,0,1,4)+0.1);
image(0:1, x, matrix(diag(h1e.pair_same_frac)[c],1,h1e.nc), col=bluered(255), axes=F, xlab="", ylab="");
for (j in 1:h1e.nc) {
	x1<-h1e.pair_same_frac[c,c][j,j];
	if (h1e.ord_csz[j]<data1.n/50) {
		text(1.6, x[j]+h1e.ord_csz[j]/2, sprintf("%2.1f", x1*100), xpd=T);
	} else {
		text(0.5, x[j]+h1e.ord_csz[j]/2, sprintf("%2.1f", x1*100));
	}
}
# }}}

#------------------------------
# cluster use cor dist
#------------------------------
# {{{
ltpm_cor_dist  <- cor_dist( norm_data1 );
h1c <- hclust(ltpm_cor_dist, method="ward.D2");
rm(ltpm_cor_dist);
h1c.nc <- 20;
h1c.cid <- cutree(h1c, k=h1c.nc);  h1c.csz <- tapply(rep(1,(data1.n)), h1c.cid, sum);
h1c.same <- h1c.cid[1:data1.n1]==h1c.cid[data1.n1+(1:data1.n1)];
h1c.pair_same<-matrix(0, h1c.nc,h1c.nc);
for (i in 1:data1.n1) {
    c1 <- h1c.cid[i];
    c2 <- h1c.cid[i+data1.n1];
    h1c.pair_same[c1,c2] <- h1c.pair_same[c1,c2] +1;
    h1c.pair_same[c2,c1] <- h1c.pair_same[c2,c1] +1;
}
h1c.pair_same_frac<-matrix(0, h1c.nc,h1c.nc);
for (i in 1:h1c.nc) {
for (j in 1:h1c.nc) {
    if (i==j) {
        h1c.pair_same_frac[i,j] <- h1c.pair_same[i,j]/h1c.csz[i];
    } else {
        h1c.pair_same_frac[i,j] <- h1c.pair_same[i,j]/(h1c.csz[i]-h1c.pair_same[i,i] + h1c.csz[j]-h1c.pair_same[j,j]);
    }
}
}

norm_data1.z <- max(abs(norm_data1));
for (i in 1:h1c.nc) {
	max_z=norm_data1.z;
    d <- norm_data1[h1c$order,ord_sm] [h1c.cid[h1c$order]==i,];
#	d[ which(d>max_z,arr.ind=T) ] = max_z; 
    g <- data1.rownames1[h1c$order][h1c.cid[h1c$order]==i];
    n4 <- sum(h1c.cid[h1c$order]==i);
    
    w = 800;
    png(paste("~/tmp/log_tpm_cor_heatmap/K20/C", sprintf("%03d",i) , ".small.png", sep=""), width=w, height=w*2, res=150);
    par(mar=c(0,0,0,0)+0.1);
    image(0:data1.m, 0:n4, t(d), col=bluered(255), axes=F, zlim=c(-max_z, max_z), xlab="", ylab="");
#    image(0:data1.m, 0:n4, t(d), col=rgb((0:255)/255,0,0), axes=F, zlim=c(0, max_z), xlab="", ylab="");
	dev.off();

	w=40*data1.m;
    png(paste("~/tmp/log_tpm_cor_heatmap/K20/C", sprintf("%03d",i) , ".large.png", sep=""), width=w, height=8*n4, res=150);
    par(mar=c(0,0,0,0)+0.1);
    image(0:data1.m, 0:n4, t(d), col=bluered(255), axes=F, zlim=c(-max_z, max_z), xlab="", ylab="");
#    image(0:data1.m, 0:n4, t(d), col=rgb((0:255)/255,0,0), axes=F, zlim=c(0, max_z), xlab="", ylab="");
	dev.off();

#    png(paste("log_tpm_cor_heatmap/K20/C", sprintf("%03d",i) , ".small.png", sep=""), width=w, height=w*2, res=150);
#    par(mar=c(10,3,3,3));
#    image(0:data1.m, 0:n4, t(d), col=redblue(255), axes=F, zlim=c(-max_z, max_z), xlab="", ylab="");
#    text(0:(m-2), rep(-5,m-1), sm_names[h1_sm$order], xpd=T, srt=-45, pos=4);
#    dev.off();
}
# }}}

# }}}

# plot trees
png(paste("~/tmp/dendrogram.euc.k",h1e.nc, ".png", sep=""), height=2048, width=data1.n*8, res=150);
plot_tree(h1e, data1.rownames1, color[h1e.cid], color[h1e.cid], cex=1, lab.cex=0.4);
dev.off();
png(paste("~/tmp/dendrogram.euc.k",h1e.nc, ".small.png", sep=""), height=2048, width=4096, res=150);
par(mar=c(3,0,0,0)+0.1);
plot_tree(h1e, rep(".",data1.n), color[h1e.cid], color[h1e.cid], lab.cex=4);
dev.off();
png(paste("~/tmp/dendrogram.cor.k",h1c.nc, ".small.png", sep=""), height=1536, width=2048, res=150);
plot_tree(h1c, rep(".",data1.n), color[h1e.cid], color[h1c.cid], pch=NA, cex=1, lab.cex=4);
dev.off();



# {{{
ps <- matrix(0, nc1,nc1);
for (i in 1:nc1) {
    a1<-(hc1[1:n1]==i) & !(hc1[n1+(1:n1)]==i) | !(hc1[1:n1]==i) & (hc1[n1+(1:n1)]==i);
    b1 <- hc1==i & rep(a1,2);
    c1 <- (1:(n1*2))[b1];
    m1 <- length(c1);
for (j in 1:nc1) {
    if (i==j) { next; }
    a2<-(hc1[1:n1]==j) & !(hc1[n1+(1:n1)]==j) | !(hc1[1:n1]==j) & (hc1[n1+(1:n1)]==j);
    b2 <- hc1==j & rep(a2,2);
    c2 <- (1:(n1*2))[b2];
    m2 <- length(c2);
    nn <- 0;
    for (run in 1:10000) {
        c <- sample(c(c1,c2), replace=F);
        c3 <- c[1:m1]%%n1;  c4<-c[(m1+1):(m1+m2)]%%n1;
        d <- sum(duplicated(c3,c4));
        if (d>pair_same[i,j]) { nn=nn+1;}
    }
    ps[i,j]<-nn/10000;
}
}
# }}}

jpeg("norm_log_tpm.heatmap.jpg", width=20*m, height=14*n1, res=150);
image(t(norm_data1[h1$order,]), col=redblue(255), axes=F, zlim=c(-max_z, max_z));
dev.off();

library(gplots);
n2=n1;
jpeg("heatmap.jpg", width=100+m*40+n2/4, height=400+n2*10*2, res=150);
hm.2<-heatmap.2(as.matrix(rbind(tpmA[1:n2,2:m],tpmB[1:n2,2:m])), distfun=cor_dist, labRow=rep(gene_names[1:n2],2), scale='row', col=bluered(255), trace='none', lmat=rbind(c(4,3),c(2,1)), lwid=c(100+n2/4,40*m), lhei=c(300,10*n2*2), cexCol=0.8, margin=c(16,16) );
dev.off();

#---------------------------------
# t-SNE or PCA test
#---------------------------------
tsne1e <- tsne(dist(data1))
pca1e <- princomp(data1)
pca1c <- princomp(norm_data1, cor=T)



# }}}

##############################################
# Gene Set Analysis
##############################################
# {{{
id_to_gs10 <- list();  id_to_gs1 <- list();  gs1<-list();
for (f in c("GO2_MF", "GO2_BP", "GO2_CC", "Reactome", "KEGG") ) {
    id_to_gs10[[f]] <- id_to_gs0[[f]][ (id_to_gs0[[f]][,1] %in% rownames(data1)), ];
    id_to_gs1[[f]] <- split(id_to_gs10[[f]][,2], id_to_gs10[[f]][,1]);
    gs1[[f]] <- split(id_to_gs10[[f]][,1], id_to_gs10[[f]][,2]);
}
corr_all_gs <- c();
group_names <- c();
group_means <- c();
ps1 <- c();
ps2 <- c();
for (i in 1:length(gs1[[f]])) {
    aid <- names(gs1[[f]])[i];
    aname <- paste( gs_id_name[[f]][aid,3], as.character( gs_id_name[[f]][aid,2]), sep="|" );
    sz <- length(gs1[[f]][[i]]);
    if (sz>0) {
        a1 <- (1:n1)[ rownames(tpmA) %in% gs1[[f]][[i]] ]
        a2 <- (1:n1)[ rownames(tpmB) %in% gs1[[f]][[i]] ]
        a <- unique(sort( c(a1,a2) ));
        if (length(a)>=20) {
            tmp.corr <-  corr1[a];
            p1 <- wilcox.test(corr1, tmp.corr, alternative="less")$p.value;
            p2 <- wilcox.test(corr1, tmp.corr, alternative="greater")$p.value;
            corr_all_gs <- rbind(corr_all_gs, cbind(rep(i, length(tmp.corr)), tmp.corr) );
            group_names <- c(group_names, aname);
            group_means <- c(group_means, mean(tmp.corr));
            ps1 <- c(ps1, p1);
            ps2 <- c(ps2, p2);
        }
    }
}
colnames(corr_all_gs) <- c("annot_idx",  "corr");
######## Boxplot ##############
n4 <- length(group_names);
png( paste("corr_by_", f, ".png", sep=""), w=n4*30, h=800, res=150 );
col <- rep(NA, length(group_names));
col[ ps1<0.01 ] = "pink";
col[ ps1<0.001 ] = "red";
col[ ps2<0.01 ] = "lightblue";
col[ ps2<0.001 ] = "blue";
par(mar=c(12,2,2,2));  boxplot(corr~annot_idx, corr_all_gs[,], names=as.character(group_names), axes=F, col=col); axis(2);
text( (1:length(group_names))-0.5, par()$usr[3], group_names, xpd=T, srt=-45, pos=4, cex=0.8, col=col);
dev.off();
# }}}


############################################
# Looking for relation between 'corr' and 'iden', expressed tissue counts, exon lost, CNE lost
############################################
# {{{
dataC <- data.frame(iden=ohno_pair_f[,7], div=max_norm_ltpm_div, specfic_cout=specific_count);
#noise <- matrix(runif(r1*m,0,0.01), n1,m);
dataC$corr <- rep(0,n1);
dataC$corr1 <- rep(0,n1);
for (i in 1:n1) {
	tmp.var <- var(t(ltpm7A[i,]),t(ltpm7B[i,]));
	if (tmp.var>0) {
		dataC$corr[i]  = cor(t(ltpm7A[i,]),t(ltpm7B[i,])); 
		dataC$corr1[i] = cor(t(ltpm7A[i,]),t(ltpm7B[i,]), method="spearman");
	}
}
#var1 <- c(); for (i in 1:n1) { var1 <- c(var1, var(t(ltpmA[i,]),t(ltpmB[i,]))); }
dataC$diff_p <- rep(1,n1);
dataC$corr_p <- rep(1,n1);
for (i in 1:n1) {
	dataC$diff_p[i] <- t.test(as.numeric(ltpm7A[i,]), as.numeric(ltpm7B[i,]) )$p.value;
	dataC$corr_p[i] <- cor.test(as.numeric(ltpm7A[i,]), as.numeric(ltpm7B[i,]), alternative="greater" )$p.value;
}
dataC$diff_p[ is.na(dataC$diff_p) ] = 1;
dataC$corr_p[ is.na(dataC$corr_p) ] = 1;
dataC$comb_p <- apply(cbind(dataC$diff_p, 1-dataC$corr_p), 1, min);
dataC$comb_logp <- -log10(dataC$comb_p);
dataC$comb_logp[ is.infinite(dataC$comb_logp) ] = 12;
dataC$comb_logp[ dataC$comb_logp>12 ] = 12;
dataC$euc = sqrt(apply( (ltpmA-ltpmB)^2, 1, sum));
max_euc = ceiling(max(dataC$euc));
#dataC$exp_sd1 <- apply(cbind(dataC$exp_sdA, dataC$exp_sdB),1,min);
#dataC$exp_sd2 <- apply(cbind(dataC$exp_sdA, dataC$exp_sdB),1,max);
#dataC$exp_mean1 <- apply(cbind(dataC$exp_meanA, dataC$exp_meanB),1,min);
#dataC$exp_mean2 <- apply(cbind(dataC$exp_meanA, dataC$exp_meanB),1,max);


dataAB <- data.frame( iden=rep(ohno_pair_f[,7],2), 
		div       = rep(max_norm_ltpm_div,2), 
		corr      = rep(dataC$corr,2), 
		comb_p    = rep(dataC$comb_p,2), 
		comb_logp = rep(dataC$comb_logp,2), 
		euc       = rep(dataC$euc,2),
		exp_sd    = c(apply(ltpm7A, 1, sd), apply(ltpm7B, 1, sd)),
		exp_mean  = c(apply(ltpm7A, 1, mean), apply(ltpm7B, 1, mean))
		);
dataAB$exp_sd_div_mean <- dataAB$exp_sd/dataAB$exp_mean;

#-----------------------------------------------
# plot identity histogram for each corr group
#-----------------------------------------------
# {{{
iden_group <- ceiling( (dataC$iden-86)/2 );
iden_group[ iden_group<=0 ] = 0;
iden_group[ iden_group>=5 ] = 5;
dataC$iden_group <- iden_group;
dataAB$iden_group <- rep(iden_group,2);
xlab=c("<86", paste( 86+(0:3)*2, 86+(1:4)*2, sep="-"), ">94");
png("/Users/chenz11/Documents/goldfish_project/tables and figures/WGD/iden_corr.boxplot.png", width=1500, height=1000, res=150)
boxplot( corr ~ iden_group, dataC, names=xlab);
dev.off();
iden_corr.lm <- lm(corr~iden, dataC);
summary(iden_corr.lm);

corr_hist <- list();
corr_hist.ymax <- 0;
for (i in 1:6) {
	corr_hist[[i]] <- hist(corr[iden_group==i-1], breaks=c(-1,0,0.2,0.4,0.6,0.8,1.0), plot=F);
	max1 <- max(corr_hist[[i]]$density); 
	if (max1 > corr_hist.ymax) { corr_hist.ymax=max1; }
}
corr_hist.percent <- corr_hist[[1]]$counts*100/sum(corr_hist[[1]]$counts);
for (i in 2:6) {
	corr_hist.percent <- cbind(corr_hist.percent, corr_hist[[i]]$counts*100/sum(corr_hist[[i]]$counts) );
}
colnames(corr_hist.percent)<-xlab;
rownames(corr_hist.percent)<- c("<0", "0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0");

euc_hist <- list();
euc_hist.ymax <- 0;
for (i in 1:6) {
	euc_hist[[i]] <- hist(dataC$euc[dataC$iden_group==i-1], breaks=c(seq(0,8,2),max_euc), plot=F);
	max1 <- max(euc_hist[[i]]$density); 
	if (max1 > euc_hist.ymax) { euc_hist.ymax=max1; }
}
euc_hist.percent <- euc_hist[[1]]$counts*100/sum(euc_hist[[1]]$counts);
for (i in 2:6) {
	euc_hist.percent <- cbind(euc_hist.percent, euc_hist[[i]]$counts*100/sum(euc_hist[[i]]$counts) );
}
colnames(euc_hist.percent)<-xlab;
rownames(euc_hist.percent)<-c("0-2", "2-4", "4-6", "6-8", ">=8");

par(mar=c(3,3,3,10));
plot(euc_hist[[1]]$mids, euc_hist[[1]]$density, ylim=c(0,max1), type='l', col=color[1]);
k1 <- length(euc_hist[[1]]$mids);
for (i in 1:6) {
	polygon( c(euc_hist[[i]]$mids[1],euc_hist[[i]]$mids, euc_hist[[i]]$mids[k1]), c(0,euc_hist[[i]]$density,0), col=color8a[i], border=F );
}
for (i in 1:6) {
	lines(euc_hist[[i]]$mids, euc_hist[[i]]$density, col=color8[i]);
}
legend(x=max_euc+1, y=max1, legend=xlab, text.col=2:7, fill=2:7, box.lwd=0, border=2:7, xpd=T)


plot(corr_hist[[1]]$mids, corr_hist[[1]]$density, ylim=c(0,corr_hist.ymax), type='l', col=greenred(11)[1]);
for (i in 2:5) {
    lines(corr_hist[[i]]$mids, corr_hist[[i]]$density, type='l', col=greenred(11)[i]);
}
# }}}

#-----------------------------------------------
# number of expressed tissue vs correlation
#-----------------------------------------------
# result: weak negative correlation, lm.p=0.0259, R^2=0.00038
# {{{
dataAB$exp_sm_count <- c(exp_sm_countA, exp_sm_countB);
dataC$exp_sm_count <- exp_sm_count
exp_corr_hist <- list();
exp_corr_hist.ymax <- 0;
for (i in 1:7) {
	exp_corr_hist[[i]] <- hist(dataC$corr[dataC$exp_sm_count==i], breaks=c(-1,0,0.2,0.4,0.6,0.8,1.0), plot=F);
	max1 <- max(exp_corr_hist[[i]]$density); 
	if (max1 > exp_corr_hist.ymax) { exp_corr_hist.ymax=max1; }
}
exp_corr_hist.density <- exp_corr_hist[[1]]$density;
exp_corr_hist.percent <- exp_corr_hist[[1]]$counts*100/sum(exp_corr_hist[[1]]$counts);
for (i in 2:7) {
	exp_corr_hist.density <- cbind(exp_corr_hist.density, exp_corr_hist[[i]]$density);
	exp_corr_hist.percent <- cbind(exp_corr_hist.percent, exp_corr_hist[[i]]$counts*100/sum(exp_corr_hist[[i]]$counts) );
}
colnames(exp_corr_hist.density)<-1:7;
rownames(exp_corr_hist.density)<- c("<0", "0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0");
colnames(exp_corr_hist.percent)<-1:7;
rownames(exp_corr_hist.percent)<- c("<0", "0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0");

exp_euc_hist <- list();
exp_euc_hist.ymax <- 0;
max_euc = ceiling(max(dataC$euc))
for (i in 1:7) {
	exp_euc_hist[[i]] <- hist(dataC$euc[dataC$exp_sm_count==i], breaks=c(0:5,max_euc), plot=F);
	max1 <- max(exp_euc_hist[[i]]$density); 
	if (max1 > exp_euc_hist.ymax) { exp_euc_hist.ymax=max1; }
}
exp_euc_hist.percent <- exp_euc_hist[[1]]$counts*100/sum(exp_euc_hist[[1]]$counts);
for (i in 2:7) {
	exp_euc_hist.percent <- cbind(exp_euc_hist.percent, exp_euc_hist[[i]]$counts*100/sum(exp_euc_hist[[i]]$counts) );
}
colnames(exp_corr_hist.percent)<-1:7;
rownames(exp_corr_hist.percent)<- c( paste(0:4,1:5,sep="-"), ">5");
# }}}

#-----------------------------------------------
# CNE lost ~~ corr or euc
#-----------------------------------------------
# {{{
# 1: gene ID
# 2-6: exon counts: all, not CC not ZF, CC not ZF, ZF not CC, CC ZF
# 7-11: lost bp
# 12-16: lost element
# 17-31: same for CNE
lost_count <- read.table("~//data/goldfish/11549472/sergey_canu70x/arrow/carAur03/big/WGD/exon_CNE.by_chain.AB.tab37.gene_count.5K.txt", head=F, sep="\t");
rownames(lost_count)<-as.character(lost_count[,1])
lost_exon_countA <- matrix(0,n1,15);
lost_CNE_countA <- matrix(0,n1,15);
rownames(lost_exon_countA) <- as.character( rownames(tpm7A) );
rownames(lost_CNE_countA) <- as.character( rownames(tpm7A) );
a <- (1:n1)[ rownames(lost_exon_countA) %in% rownames(lost_count) ];
lost_exon_countA[a, ] <- as.matrix(lost_count[ as.character(rownames(tpm7A)[a]),2:16]);
lost_CNE_countA[a, ]  <- as.matrix(lost_count[ as.character(rownames(tpm7A)[a]),22:36]);
lost_exonA <- rep("",n1*5);  dim(lost_exonA) <- c(n1,5);
lost_exonA[a,] <- as.matrix( lost_count[ as.character(rownames(tpm7A)[a]), 17:21] );
lost_CNEA  <- rep("",n1*5);  dim(lost_CNEA) <- c(n1,5);
lost_CNEA[a,]  <- as.matrix( lost_count[ as.character(rownames(tpm7A)[a]), 37:41] );

lost_exon_countB <- matrix(0,n1,15);
lost_CNE_countB  <- matrix(0,n1,15);
rownames(lost_exon_countB) <- as.character( rownames(tpm7B) );
rownames(lost_CNE_countB)  <- as.character( rownames(tpm7B) );
a <- (1:n1)[ rownames(lost_exon_countB) %in% rownames(lost_count) ];
lost_exon_countB[a, ] <- as.matrix(lost_count[ as.character(rownames(tpm7B)[a]),2:16]);
lost_CNE_countB[a, ]  <- as.matrix(lost_count[ as.character(rownames(tpm7B)[a]),22:36]);
lost_exonB <- rep("",n1*5);  dim(lost_exonB) <- c(n1,5);
lost_exonB[a,] <- as.matrix( lost_count[ as.character(rownames(tpm7B)[a]), 17:21] );
lost_CNEB  <- rep("",n1*5);  dim(lost_CNEB) <- c(n1,5);
lost_CNEB[a,]  <- as.matrix( lost_count[ as.character(rownames(tpm7B)[a]), 37:41] );

dataC$"exon_lost_n" <- lost_exon_countA[,1] + lost_exon_countB[,1] ;
dataC$"exon_lost_w" <- lost_exon_countA[,11]+ lost_exon_countB[,11] ;
dataC$"CNE_lost_n"  <- lost_CNE_countA[,1]  + lost_CNE_countB[,1];
dataC$"CNE_lost_w"  <- lost_CNE_countA[,11] + lost_CNE_countB[,11];
dataC$"exon_lost_nn" <- dataC$"exon_lost_n";
dataC$"exon_lost_nn"[dataC$"exon_lost_nn">3] = 3;
dataC$"CNE_lost_nn" <- dataC$"CNE_lost_n";
dataC$"CNE_lost_nn"[dataC$"CNE_lost_nn">3] = 3;

dataC$"exon_lost_n2" <- apply(lost_exon_countA[,4:5]   + lost_exon_countB[,4:5],1,sum) ;
dataC$"exon_lost_w2" <- apply(lost_exon_countA[,14:15] + lost_exon_countB[,14:15],1,sum) ;
dataC$"CNE_lost_n2"  <- apply(lost_CNE_countA[,4:5]    + lost_CNE_countB[,4:5],1,sum) ;
dataC$"CNE_lost_w2"  <- apply(lost_CNE_countA[,14:15]  + lost_CNE_countB[,14:15],1,sum) ;

dataAB$exon_lost_n <- c(lost_exon_countA[,1] , lost_exon_countB[,1]);
dataAB$exon_lost_w <- c(lost_exon_countA[,11], lost_exon_countB[,11]);
dataAB$CNE_lost_n  <- c(lost_CNE_countA[,1]  , lost_CNE_countB[,1]);
dataAB$CNE_lost_w  <- c(lost_CNE_countA[,11] , lost_CNE_countB[,11]);
dataAB$"exon_lost_nn" <- dataAB$"exon_lost_n";
dataAB$"exon_lost_nn"[dataAB$"exon_lost_nn">3] = 3;
dataAB$"CNE_lost_nn" <- dataAB$"CNE_lost_n";
dataAB$"CNE_lost_nn"[dataAB$"CNE_lost_nn">3] = 3;

# boxplot euc ~ exon_lost_nn
# {{{
local({
# W = 1521300, p-value = 5.871e-07
wilcox.test(dataC$euc[dataC$exon_lost_nn==0], dataC$euc[dataC$exon_lost_nn==1], alternative="less");
# W = 12147, p-value = 0.0001122
wilcox.test(dataC$euc[dataC$exon_lost_nn==1], dataC$euc[dataC$exon_lost_nn==2], alternative="less");
# W = 5716, p-value = 0.09654
wilcox.test(dataC$euc[dataC$exon_lost_nn==2], dataC$euc[dataC$exon_lost_nn==3], alternative="less");
png("exon_lost-euc.png", width=900, height=900, res=150);
par(mar=c(4,4,5,3));
boxplot(euc ~ exon_lost_nn, dataC, names=c(0, 1, 2, ">=3"), border=c("black",rep('red',3)), axes=F); axis(2); text(1:4,rep(-2,4),c(0, 1, 2, ">=3"), xpd=T); 
lines(matrix(c(1.05,1.05,1.95,1.95, 28, 29, 29, 25),4,2), type='l', xpd=T); text(1.5,30.5, paste("",5.87e-06,sep=""), xpd=T); text(1.5,28,"***", col='red', xpd=T);
lines(matrix(c(2.05,2.05,2.95,2.95, 24, 27, 27, 23),4,2), type='l', xpd=T); text(2.5,28.5, paste("",0.00011,sep=""), xpd=T); text(2.5,26,"***", col='red', xpd=T);
lines(matrix(c(3.05,3.05,3.95,3.95, 23, 25, 25, 20),4,2), type='l', xpd=T); text(3.5,26.5, paste("",0.09,sep=""), xpd=T);
dev.off();
});
# }}}

# boxplot corr ~ exon_lost_nn
# {{{
local({
# W = 1521300, p-value = 5.871e-07
p0<<-wilcox.test(dataC$corr[dataC$exon_lost_nn==0], dataC$corr[dataC$exon_lost_nn==1], alternative="less");
# W = 12147, p-value = 0.0001122
p1<<-wilcox.test(dataC$corr[dataC$exon_lost_nn==1], dataC$corr[dataC$exon_lost_nn==2], alternative="less");
# W = 5716, p-value = 0.09654
p2<<-wilcox.test(dataC$corr[dataC$exon_lost_nn==2], dataC$corr[dataC$exon_lost_nn==3], alternative="less");
png("exon_lost-corr.png", width=900, height=900, res=150);
par(mar=c(4,4,5,3));
boxplot(corr ~ exon_lost_nn, dataC, names=c(0, 1, 2, ">=3"), border=c("black",rep('red',3)), axes=F); axis(2); text(1:4,rep(-2,4),c(0, 1, 2, ">=3"), xpd=T); 
#lines(matrix(c(1.05,1.05,1.95,1.95, 28, 29, 29, 25),4,2), type='l', xpd=T); text(1.5,30.5, paste("",5.87e-06,sep=""), xpd=T); text(1.5,28,"***", col='red', xpd=T);
#lines(matrix(c(2.05,2.05,2.95,2.95, 24, 27, 27, 23),4,2), type='l', xpd=T); text(2.5,28.5, paste("",0.00011,sep=""), xpd=T); text(2.5,26,"***", col='red', xpd=T);
#lines(matrix(c(3.05,3.05,3.95,3.95, 23, 25, 25, 20),4,2), type='l', xpd=T); text(3.5,26.5, paste("",0.09,sep=""), xpd=T);
dev.off();
});
# }}}

boxplot(corr ~ exon_lost_nn, dataC, names=c(0, 1, 2, ">=3") )
a1<- cbind(lost_exon_countA+lost_exon_countB, corr, corr1, max_norm_ltpm_div, norm_ltpm_p);
a2<- cbind(lost_CNE_countA+lost_CNE_countB  , corr, corr1, max_norm_ltpm_div, norm_ltpm_p);
colnames(a1) <- c(paste("n",1:5,sep=""),  paste("l",1:5,sep=""),  paste("w",1:5,sep=""),  "corr", "corr1", "div", "wp");
colnames(a2) <- c(paste("n",1:5,sep=""),  paste("l",1:5,sep=""),  paste("w",1:5,sep=""),  "corr", "corr1", "div", "wp");
ew0 <- a1[,"w1"]
cw0 <- a2[,"w1"]
cw2 <- a2[,"w4"] + a2[,"w5"]
en0 <- a1[,"n1"]
cn0 <- a2[,"n1"]
cn2 <- a2[,"n4"] + a2[,"n5"]

ew11<-sum(corr[en0==0 & cn0==0]>=0.5)
ew12<-sum(corr[en0==0 & cn0==0]<0.5)
ew21<-sum(corr[ew0>0 & cn0==0]>=0.5)
ew22<-sum(corr[ew0>0 & cn0==0]<0.5)
fisher.test(matrix(c(ew11,ew12,ew21,ew22),2,2))

ew11<-sum(en0==0 & cn0==0 & dataC$comb_p>0.5)
ew12<-sum(en0> 0 & cn0==0 & dataC$comb_p>0.5)
ew21<-sum(en0==0 & cn0==0 & (dataC$comb_p<0.1))
ew22<-sum(en0> 0 & cn0==0 & (dataC$comb_p<0.1))
fisher.test(matrix(c(ew11,ew12,ew21,ew22),2,2))

cw11<-sum(cn0==0 & en0==0 & dataC$comb_p>0.5)
cw12<-sum(cn0> 0 & en0==0 & dataC$comb_p>0.5)
cw21<-sum(cn0==0 & en0==0 & (dataC$comb_p<0.1))
cw22<-sum(cn0> 0 & en0==0 & (dataC$comb_p<0.1))
fisher.test(matrix(c(cw11,cw12,cw21,cw22),2,2))

cw11<-sum(corr[cn0==0 & en0==0 & specific_count>0 & specific_count<=1]>=0.8)
cw12<-sum(corr[cn0==0 & en0==0 & specific_count>0 & specific_count<=1]<0.8)
cw21<-sum(corr[cw0/cn0==2 & en0==0 & specific_count>0 & specific_count<=1]>=0.8)
cw22<-sum(corr[cw0/cn0==2 & en0==0 & specific_count>0 & specific_count<=1]<0.8)
fisher.test(matrix(c(cw11,cw12,cw21,cw22),2,2))

cn11<-sum(cn0==0 & en0==0 & specific_count==0)
cn12<-sum(cn0>0 & en0==0 & specific_count==0)
cn21<-sum(cn0==0 & en0==0 & specific_count>1)
cn22<-sum(cn0>0 & en0==0 & specific_count>1)
fisher.test(matrix(c(cn11,cn12,cn21,cn22),2,2))

cn11<-sum(cn0==0 & en0==0 & max_norm_ltpm_div<=0.25)
cn12<-sum(cn0>0  & en0==0 & max_norm_ltpm_div<=0.25)
cn21<-sum(cn0==0 & en0==0 & max_norm_ltpm_div>=4)
cn22<-sum(cn0>0  & en0==0 & max_norm_ltpm_div>=4)
fisher.test(matrix(c(cn11,cn12,cn21,cn22),2,2))


cw11<-sum(corr[cn0==0 & en0==0 & ohno_pair_f[,7]<90 ]>=0.5)
cw12<-sum(corr[cn0==0 & en0==0 & ohno_pair_f[,7]<90 ]<0.5)
cw21<-sum(corr[cw0>1 & en0==0 & ohno_pair_f[,7]<90 ]>=0.5)
cw22<-sum(corr[cw0>1 & en0==0 & ohno_pair_f[,7]<90 ]<0.5)
fisher.test(matrix(c(cw11,cw12,cw21,cw22),2,2))
boxplot(corr[b2==0 & a1[,1]==0], corr[b2>0 & a1[,1]==0]) 


b <- a[ specific_count>0, ]
boxplot( corr ~ CNE_lost, b);


x1 <- lost_countA[,4]>0 & f2 & specific_count>0;
x2 <- lost_countB[,4]>0 & f2 & specific_count>0;
a1 <- cbind(rownames(tpmA)[x1], rownames(tpmB)[x1], lost_CNE_A[x1,1])
a2 <- cbind(rownames(tpmA)[x2], rownames(tpmB)[x2], lost_CNE_B[x2,1])
write.table(rbind(a1,a2), "pairs.txt", sep="\t", quote=F, col.names=F, row.names=F)
# }}}

#-----------------------------------------------
# tissue exp SD ~ CNE lost
#-----------------------------------------------
tmp.f1 <- dataAB$CNE_lost_n[1:n1]==0 & dataAB$CNE_lost_n[n1+(1:n1)]>0;
tmp.f2 <- dataAB$CNE_lost_n[1:n1]>0 & dataAB$CNE_lost_n[n1+(1:n1)]==0;
# W = 570110, p-value = 0.0003266
wilcox.test(c(dataAB$exp_sd[1:n1][tmp.f1], dataAB$exp_sd[n1+(1:n1)][tmp.f2]), c(dataAB$exp_sd[n1+(1:n1)][tmp.f1],dataAB$exp_sd[1:n1][tmp.f2]))
# W = 555320, p-value = 5.063e-06
wilcox.test(c(dataAB$exp_mean[1:n1][tmp.f1], dataAB$exp_mean[n1+(1:n1)][tmp.f2]), c(dataAB$exp_mean[n1+(1:n1)][tmp.f1],dataAB$exp_mean[1:n1][tmp.f2]))
# W = 647160, p-value = 0.03583
wilcox.test(c(dataAB$exp_sd_div_mean[1:n1][tmp.f1], dataAB$exp_sd_div_mean[n1+(1:n1)][tmp.f2]), c(dataAB$exp_sd_div_mean[n1+(1:n1)][tmp.f1],dataAB$exp_sd_div_mean[1:n1][tmp.f2]))


# }}}

##########################################
# For each cluster analysis
##########################################
# {{{
cdata_list <- list();
for (i in 1:h1e.nc) {
	tmp.f <- h1e.cid==i;
	d <- data.frame(same=rep(h1e.same[tmp.f],2));
	for (name in names(dataAB) ) {
		d[,name] <- dataAB[f5,][tmp.f,name];
	}
	cdata_list[[i]] <-  d;
}

cdata <- data.frame(same_frac=diag(h1e.pair_same_frac), csz=h1e.csz);
cdata$exp_mean <- rep(0,h1e.nc);
for (i in 1:h1e.nc) {
	d <- data1[ h1e.cid==i, ];
	cdata$exp_mean[i] <- sum(d)/(h1e.csz[i]*data1.m);
	cdata$mean_exp_mean[i] <- mean(apply(d,1,mean));
	cdata$mean_exp_sd[i] <- mean(apply(d,1,sd));
	cdata$ngene_lost_exon[i] <- sum(cdata_list[[i]]$exon_lost_n>0);
	cdata$ngene_lost_CNE[i]  <- sum(cdata_list[[i]]$CNE_lost_n>0);
	cdata$ngene_lost_exon_or_CNE[i]  <- sum(cdata_list[[i]]$exon_lost_n+cdata_list[[i]]$CNE_lost_n>0);
	cdata$fgene_lost_exon[i] <- cdata$ngene_lost_exon[i] / h1e.csz[i];
	cdata$fgene_lost_CNE[i] <- cdata$ngene_lost_CNE[i] / h1e.csz[i];
	cdata$fgene_lost_exon_or_CNE[i] <- cdata$ngene_lost_exon_or_CNE[i] / h1e.csz[i];
	cdata$mean_iden[i]       <- mean(cdata_list[[i]]$iden);
	cdata$mean_corr[i]       <- mean(cdata_list[[i]]$corr);
	cdata$mean_comb_logp[i]  <- mean(cdata_list[[i]]$comb_logp);
}
plot(mean_exp_sd ~ mean_exp_mean, cdata, col="red");
text(cdata$mean_exp_mean, cdata$mean_exp_sd, 1:h1e.nc);
a <- cbind(
		(cdata$mean_exp_mean - mean(cdata$mean_exp_mean)) / sd(cdata$mean_exp_mean),
		(cdata$mean_iden - mean(cdata$mean_iden)) / sd(cdata$mean_iden),
		(cdata$fgene_lost_exon - mean(cdata$fgene_lost_exon)) / sd(cdata$fgene_lost_exon),
		(cdata$fgene_lost_CNE - mean(cdata$fgene_lost_CNE)) / sd(cdata$fgene_lost_CNE),
		(cdata$mean_exp_sd - mean(cdata$mean_exp_sd)) / sd(cdata$mean_exp_sd)
		);
# }}}

##############################################
# singleton gene 
##############################################
lost_count.exon_num<-bgp.exon_num[ rownames(lost_count) ]



