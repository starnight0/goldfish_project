require('Rcpp');
require('dendextend');
###############################
# Sub Functions
###############################
#{{{
remove_id_version = function(v) { v=gsub("\\.[0-9]+","",v); }

cppFunction('NumericVector pair_apply(NumericMatrix x, NumericMatrix y, Function func) {
	int nrow = x.nrow();
	NumericVector out(nrow, 0.0);
	for (int i=0; i<nrow; i++) {
		NumericMatrix::Row row1 = x( i, _);
		NumericMatrix::Row row2 = y( i, _);
		out[i] = as<double>(func(row1, row2));
	}
	return out;
}');

cppFunction('NumericVector pair_euc_dist(NumericMatrix x, NumericMatrix y) {
	int nrow = x.nrow();
	int nc = x.ncol();
	NumericVector out(nrow, 0.0);
	for (int i=0; i<nrow; i++) {
		double d=0;
		for (int j=0; j<nc; j++) {
			d += pow(x(i,j)-y(i,j),2);
		}
		d = sqrt(d);
		out[i] = d;
	}
	return out;
}');

cppFunction('NumericVector pair_cor(NumericMatrix x, NumericMatrix y) {
	int nrow = x.nrow();
	int nc = x.ncol();
	NumericVector out(nrow, 0.0);
	for (int i=0; i<nrow; i++) {
		double m1=0, m2=0, var1=0, var2=0;
		double cov=0;
		for (int j=0; j<nc; j++) {
			m1 += x(i,j);
			m2 += y(i,j);
		}
		m1/=nc;
		m2/=nc;
		for (int j=0; j<nc; j++) {
			double x1 = x(i,j)-m1;
			double y1 = y(i,j)-m2;
			cov += x1*y1;
			var1 += pow(x1, 2);
			var2 += pow(y1, 2);
		}
		cov/=nc-1;
		var1/=nc-1;
		var2/=nc-1;
		if (nc==1 || var1==0 || var2==0) { out[i]=0; }
		else { out[i] = cov/sqrt(var1*var2); }
	}
	return out;
}');

#cppFunction('NumericVector pair_cor_p(NumericMatrix x, NumericMatrix y) {
#	// Obtain environment containing function
#	Environment base("package:stats"); 
#
#	// Make function callable from C++
#	Function cor_test = base["cor.test"];    
#
#	int nrow = x.nrow();
#	NumericVector out(nrow, 0.0);
#	for (int i=0; i<nrow; i++) {
#		NumericMatrix::Row x1 = x(i,_);
#		NumericMatrix::Row y1 = y(i,_);
#		// Call the function and receive its list output
#		List test_out = cor_test(_["x"]=x1, _["y"]=y1);
#		out[i] = test_out["p.value"];
#	//	if (out[i]==NA_REAL) { out[i]=1.0; }
#	}
#	return out;
#}');

cppFunction('NumericVector pair_cor_p(NumericMatrix x, NumericMatrix y) {
	// Obtain environment containing function
//	Environment base("package:stats"); 

	// Make function callable from C++
//	Function t_test = base["t.test"];    

	int nrow = x.nrow();
	int nc = x.ncol();
	NumericVector out(nrow, 0.0);
	for (int i=0; i<nrow; i++) {
		double m1=0, m2=0, var1=0, var2=0;
		double cov=0;
		for (int j=0; j<nc; j++) {
			m1 += x(i,j);
			m2 += y(i,j);
		}
		m1/=nc;
		m2/=nc;
		for (int j=0; j<nc; j++) {
			double x1 = x(i,j)-m1;
			double y1 = y(i,j)-m2;
			cov += x1*y1;
			var1 += pow(x1, 2);
			var2 += pow(y1, 2);
		}
		cov/=nc-1;
		var1/=nc-1;
		var2/=nc-1;
		if (var1==0 || var2==0) { out[i]=0.0; }
		else {
			double r = cov/sqrt(var1*var2);
			if (r>=1 || r<= -1) { out[i]=0.0; }
			else {
				double t = r*sqrt(nc-2)/sqrt(1-r*r);
				if (t<0) t=-t;
				out[i] = (1.0 - R::pt(t, nc-2, NA_REAL, false))*2;
			}
		}
//		NumericMatrix::Row x1 = x(i,_);
//		NumericMatrix::Row y1 = y(i,_);
		// Call the function and receive its list output
//		List test_out = t_test(_["x"]=x1, _["y"]=y1);
//		out[i] = test_out["p.value"];
	//	if (out[i]==NA_REAL) { out[i]=1.0; }
	}
	return out;
}');

cppFunction('NumericVector pair_t_p(NumericMatrix x, NumericMatrix y) {
	// Obtain environment containing function
//	Environment base("package:stats"); 

	// Make function callable from C++
//	Function t_test = base["t.test"];    

	int nrow = x.nrow();
	int nc = x.ncol();
	NumericVector out(nrow, 0.0);
	for (int i=0; i<nrow; i++) {
		double m=0, var=0;
		for (int j=0; j<nc; j++) {
			m += x(i,j)-y(i,j);
		}
		m/=nc;
		for (int j=0; j<nc; j++) {
			double z = x(i,j)-y(i,j)-m;
			var += pow(z, 2);
		}
		var/=nc-1;
		if (var==0) {
			if (m==0) out[i]=0.0;
			else out[i]=1.0;
		} else {
			double t = m*sqrt(nc)/sqrt(var);
			if (t<0) t=-t;
			out[i] = (1.0 - R::pt(t, nc-1, NA_REAL, false))*2;
		}
	}
	return out;
}');


pair_apply2 = function(x,y, FUN)
{
	if (ncol(x)!=ncol(y)) { return(NA); }
	n=ncol(x);
	r=apply(cbind(x,y), 1, FUN);
}

pair_cor2 = function(x, y)
{
	n=ncol(x);
	res=apply(cbind(x, y), 1, function(z) { m=length(z)/2; cor(z[1:m],z[(m+1):(2*m)]); });
	res[ is.na(res) ] = 0;
	res
}

pair_cor_p2 = function(x, y)
{
	res=apply(cbind(x, y), 1, function(z) { n=length(z)/2; cor.test(z[1:n],z[(n+1):(2*n)])$p.value; });
	res[ is.na(res) ] = 1.0;
	res
}

pair_t_p2 = function(x, y)
{
	res = apply(cbind(x, y), 1, function(z) { n=length(z)/2; t.test(z[1:n],z[(n+1):(2*n)],paired=T)$p.value; });
	res[ is.na(res) ] = 1.0;
	res
}

arrange_data = function(data, nts, names, colnames) {
	m=nrow(data); n=ncol(data); 
	x=matrix(as.vector(t(data)), m*n/nts, nts, byrow=T);
	rownames(x) = unlist(t(names));
	colnames(x) = colnames
	x
} 

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

fate_classify = function(data, nts, i0, i1, i2, thres, thres1, thres2, fc, C1, C3, cor_p, tp) {
	a0=data[,(nts*(i0-1)+1):(nts*i0)];
	a1=data[,(nts*(i1-1)+1):(nts*i1)];
	a2=data[,(nts*(i2-1)+1):(nts*i2)];
	a3=a1+a2;
# log2(fpkm+1)
	la0 = log2(a0+1);  la1 = log2(a1+1);  la2 = log2(a2+1);  la3 = log2(a3+1);
	max0 = apply(a0, 1, max);   max1 = apply(a1, 1, max);   max2 = apply(a2, 1, max);
	mean0 = apply(a0,1,mean);   mean1 = apply(a1,1,mean);   mean2=apply(a2,1,mean);   mean3 = apply(a3,1,mean);
	f4 = max0>=thres & max1>=thres & max2>=thres;

	r01 = pair_cor(la0,la1);   r01.p = pair_cor_p(la0,la1);   tp01 = pair_t_p(la0, la1);
	r02 = pair_cor(la0,la2);   r02.p = pair_cor_p(la0,la2);   tp02 = pair_t_p(la0, la2);
	r03 = pair_cor(la0,la3);   r03.p = pair_cor_p(la0,la3);   tp03 = pair_t_p(la0, la3);
	r12 = pair_cor(la1,la2);   r12.p = pair_cor_p(la1,la2);   tp12 = pair_t_p(la1, la2);

	r01_hc = r01>=C3 & r01.p<cor_p;
	r01_mc = r01>=C1;
	r01_lc = r01<C1;
	r01_de = tp01<tp;
	r02_hc = r02>=C3 & r02.p<cor_p;
	r02_mc = r02>=C1;
	r02_lc = r02<C1;
	r02_de = tp02<tp;
	r12_hc = r12>=C3 & r12.p<cor_p;
	r12_mc = r12>=C1;
	r12_lc = r12<C1;
	r12_de = tp12<tp;
	r03_hc = r03>=C3 & r03.p<cor_p;
	r03_mc = r03>=C1;
	r03_lc = r03<C1;
	r03_de = tp03<tp;

# on-off
	a_on0_off1 = a1<thres1 & a0>=thres2 & (a1==0 | a0/a1>=fc);
	a_on0_off2 = a2<thres1 & a0>=thres2 & (a2==0 | a0/a2>=fc);
	a_on0_off3 = a3<thres1 & a0>=thres2 & (a3==0 | a0/a3>=fc);
	a_on1_off0 = a0<thres1 & a1>=thres2 & (a0==0 | a1/a0>=fc);
	a_on1_off2 = a2<thres1 & a1>=thres2 & (a2==0 | a1/a2>=fc);
	a_on1_off3 = a3<thres1 & a1>=thres2 & (a3==0 | a1/a3>=fc);
	a_on2_off0 = a0<thres1 & a2>=thres2 & (a0==0 | a2/a0>=fc);
	a_on2_off1 = a1<thres1 & a2>=thres2 & (a1==0 | a2/a1>=fc);
	a_on2_off3 = a3<thres1 & a2>=thres2 & (a3==0 | a2/a3>=fc);
	a_on3_off0 = a0<thres1 & a3>=thres2 & (a0==0 | a3/a0>=fc);
	a_on3_off1 = a1<thres1 & a3>=thres2 & (a1==0 | a3/a1>=fc);
	a_on3_off2 = a2<thres1 & a3>=thres2 & (a2==0 | a3/a2>=fc);


	n_on0_off1 = apply( a_on0_off1, 1, sum);
	n_on0_off2 = apply( a_on0_off2, 1, sum);
	n_on0_off3 = apply( a_on0_off3, 1, sum);
	n_on1_off0 = apply( a_on1_off0, 1, sum);
	n_on1_off2 = apply( a_on1_off2, 1, sum);
	n_on1_off3 = apply( a_on1_off3, 1, sum);
	n_on2_off0 = apply( a_on2_off0, 1, sum);
	n_on2_off1 = apply( a_on2_off1, 1, sum);
	n_on2_off3 = apply( a_on2_off3, 1, sum);
	n_on3_off0 = apply( a_on3_off0, 1, sum);
	n_on3_off1 = apply( a_on3_off1, 1, sum);
	n_on3_off2 = apply( a_on3_off2, 1, sum);

# greater
	gt01 = a0>=thres2 & a0>=a1*fc;
	gt02 = a0>=thres2 & a0>=a2*fc;
	gt03 = a0>=thres2 & a0>=a3*fc;
	gt10 = a1>=thres2 & a1>=a0*fc;
	gt12 = a1>=thres2 & a1>=a2*fc;
	gt13 = a1>=thres2 & a1>=a3*fc;
	gt20 = a2>=thres2 & a2>=a0*fc;
	gt21 = a2>=thres2 & a2>=a1*fc;
	gt23 = a2>=thres2 & a2>=a3*fc;
	gt30 = a3>=thres2 & a3>=a0*fc;
	gt31 = a3>=thres2 & a3>=a1*fc;
	gt32 = a3>=thres2 & a3>=a2*fc;

	n_gt01 = apply(a0>=thres2 & a0>=a1*fc, 1, sum);
	n_gt02 = apply(a0>=thres2 & a0>=a2*fc, 1, sum);
	n_gt03 = apply(a0>=thres2 & a0>=a3*fc, 1, sum);
	n_gt10 = apply(a1>=thres2 & a1>=a0*fc, 1, sum);
	n_gt12 = apply(a1>=thres2 & a1>=a2*fc, 1, sum);
	n_gt13 = apply(a1>=thres2 & a1>=a3*fc, 1, sum);
	n_gt20 = apply(a2>=thres2 & a2>=a0*fc, 1, sum);
	n_gt21 = apply(a2>=thres2 & a2>=a1*fc, 1, sum);
	n_gt23 = apply(a2>=thres2 & a2>=a3*fc, 1, sum);
	n_gt30 = apply(a3>=thres2 & a3>=a0*fc, 1, sum);
	n_gt31 = apply(a3>=thres2 & a3>=a1*fc, 1, sum);
	n_gt32 = apply(a3>=thres2 & a3>=a2*fc, 1, sum);

	n = apply(a0>=thres2,1,sum);

	on_off_conserved1  = f4 & n_gt01==0 & n_gt10==0;
	on_off_conserved2  = f4 & n_gt02==0 & n_gt20==0;
	on_off_conserved3  = f4 & n_gt03==0 & n_gt30==0;

	coexpressed1 = f4 & n_on0_off1==0 & n_on1_off0==0;
	coexpressed2 = f4 & n_on0_off2==0 & n_on2_off0==0;
	dosage_coexpressed = f4 & n_on0_off3==0 & n_on3_off0==0;
	double_coexpressed  = coexpressed1 & coexpressed2;
	single_coexpressed  = coexpressed1 | coexpressed2;
	any_coexpressed    = dosage_coexpressed | single_coexpressed;
	ab_coexpressed = f4 & n_on1_off2==0 & n_on2_off1==0;

	on_off_subfunc = f4 & apply(a_on0_off1 & a_on2_off1,1,sum)>0 & apply(a_on0_off2 & a_on1_off2,1,sum)>0 & n_on1_off0==0 & n_on2_off0==0;

	on_off_neofunc1 = f4 & n_on1_off0>0 & n_on0_off1==0 & n_on0_off2==0;
	on_off_neofunc2 = f4 & n_on2_off0>0 & n_on0_off1==0 & n_on0_off2==0;
	on_off_neofunc  = f4 & (on_off_neofunc1 | on_off_neofunc2);

	on_off_neofunc_WGD1 = f4 & apply(a_on1_off0 & a_on1_off2, 1, sum)>0 & n_on2_off1==0 & n_on0_off1==0 & n_on0_off2==0;
	on_off_neofunc_WGD2 = f4 & apply(a_on2_off0 & a_on2_off1, 1, sum)>0 & n_on1_off2==0 & n_on0_off1==0 & n_on0_off2==0;
	on_off_neofunc_WGD  = f4 & (on_off_neofunc_WGD1 | on_off_neofunc_WGD2);

	nonfunc1= max0>=thres2 & max1<thres1 & n_on0_off1>0;
	nonfunc2= max0>=thres2 & max2<thres1 & n_on0_off2>0;
	nonfunc = max0>=thres2 & (max1<thres1 & n_on0_off1>0 | max2<thres1 & n_on0_off2>0);

#	nonfunc_WGD1= max0>=thres2 & max2>=thres2 & max1<thres1 & apply(a_on0_off1& a_on2_off1,1,sum)>0;
#	nonfunc_WGD2= max0>=thres2 & max1>=thres2 & max2<thres1 & apply(a_on0_off2& a_on1_off2,1,sum)>0;
	nonfunc_WGD1= max0>=thres2 & max2>=thres2 & max1<thres1 & n_on0_off1>0 & n_on2_off1>0;
	nonfunc_WGD2= max0>=thres2 & max1>=thres2 & max2<thres1 & n_on0_off2>0 & n_on1_off2>0;
	nonfunc_WGD = nonfunc_WGD1 | nonfunc_WGD2;

	part_nonfunc1 = n_on0_off1>0 & apply(!a_on0_off1 & (a_on1_off0 | a_on0_off2 | a_on2_off0), 1, sum)==0 & !nonfunc;
	part_nonfunc2 = n_on0_off2>0 & apply(!a_on0_off2 & (a_on2_off0 | a_on0_off1 | a_on1_off0), 1, sum)==0 & !nonfunc;
	part_nonfunc = part_nonfunc1 | part_nonfunc2; 

	on_off_part_subfunc = f4 & apply(a_on0_off1 & a_on2_off1,1,sum)>0 & apply(a_on0_off2 & a_on1_off2,1,sum)>0 & !on_off_subfunc & !nonfunc;
	
	on_off_part_neofunc1 = f4 & n_on1_off0>0 & !on_off_neofunc & !nonfunc & !on_off_neofunc;
	on_off_part_neofunc2 = f4 & n_on2_off0>0 & !on_off_neofunc & !nonfunc & !on_off_neofunc;
	on_off_part_neofunc = f4 & (on_off_part_neofunc1 | on_off_part_neofunc2) & !on_off_neofunc & !nonfunc;

	part_nonfunc_WGD1 = apply(a_on0_off1 & a_on2_off1 & !a_on0_off2,1,sum)>0 & apply(!(a_on0_off1 & a_on2_off1 & !a_on0_off2) & (a_on1_off0 | a_on0_off2 | a_on2_off0), 1, sum)==0 & !nonfunc;
	part_nonfunc_WGD2 = apply(a_on0_off2 & a_on1_off2 & !a_on0_off1,1,sum)>0 & apply(!(a_on0_off2 & a_on1_off2 & !a_on0_off1) & (a_on2_off0 | a_on0_off1 | a_on1_off0), 1, sum)==0 & !nonfunc;
	part_nonfunc_WGD = part_nonfunc_WGD1 | part_nonfunc_WGD2; 

	on_off_part_neofunc_WGD1 = f4 & apply(a_on1_off0 & a_on1_off2, 1, sum)>0 & !on_off_neofunc & !nonfunc & !on_off_neofunc_WGD;
	on_off_part_neofunc_WGD2 = f4 & apply(a_on2_off0 & a_on2_off1, 1, sum)>0 & !on_off_neofunc & !nonfunc & !on_off_neofunc_WGD;
	on_off_part_neofunc_WGD  = f4 & (on_off_part_neofunc_WGD1 | on_off_part_neofunc_WGD2) & !on_off_neofunc & !nonfunc;


	res=data.frame(
		ab_high_corr     = f4 & r12_hc,
		ab_mid_corr      = f4 & r12_mc,
		ab_low_corr      = f4 & r12_lc,
		ab_diff          = f4 & tp12<tp,
		ab_high_conserved     = f4 & r12_hc & !r12_de,
		ab_mid_conserved      = f4 & r12_mc & !r12_de,
		ab_not_conserved      = f4 & (r12_lc | r12_de),

		ab_coexpressed      = ab_coexpressed,
		ab_on_off_conserved = f4 & n_gt12==0 & n_gt21==0,

		ab_all_greater1 = tp12<tp & apply(a2<a1,1,sum)==nts,
		ab_all_greater2 = tp12<tp & apply(a1<a2,1,sum)==nts,

		single_high_corr = f4 & (r01_hc | r02_hc),
		single_mid_corr  = f4 & (r01_mc | r02_mc),
		double_high_corr = f4 & (r01_hc & r02_hc),
		double_mid_corr  = f4 & (r01_mc & r02_mc),
		high_dosage_corr = f4 & r03_hc,
		mid_dosage_corr  = f4 & r03_mc,
		not_dosage_corr  = f4 & r03_lc,
		any_high_corr    = f4 & ( (r01_hc | r02_hc) | r03_hc ),
		any_mid_corr     = f4 & ( (r01_mc | r02_mc) | r03_mc ),

		single_high_conserved = f4 & (r01_hc & !r01_de | r02_hc & !r02_de),
		single_mid_conserved  = f4 & (r01_mc & !r01_de | r02_mc & !r02_de),
		double_high_conserved = f4 & (r01_hc & !r01_de & r02_hc & !r02_de),
		double_mid_conserved  = f4 & (r01_mc & !r01_de & r02_mc & !r02_de),
		high_dosage_balance   = f4 & r03_hc & !r03_de,
		mid_dosage_balance    = f4 & r03_mc & !r03_de,
		not_dosage_balance    = f4 & (r03_lc | r03_de),
		any_high_conserved    = f4 & ( r01_hc & !r01_de | r02_hc & !r02_de | r03_hc & !r03_de ),
		any_mid_conserved     = f4 & ( r01_mc & !r01_de | r02_mc & !r02_de | r03_mc & !r03_de ),
		double_diverged = f4 & ((r01_lc | r01_de) & (r02_lc | r02_de)),
		single_diverged = f4 & ((r01_lc | r01_de) | (r02_lc | r02_de)),

		on_off_conserved1 = on_off_conserved1,
		on_off_conserved2 = on_off_conserved2,
		on_off_double_conserved  = on_off_conserved1 & on_off_conserved2,
		on_off_single_conserved  = on_off_conserved1 | on_off_conserved2,
		on_off_dosage_conserved  = on_off_conserved3,
		on_off_any_conserved = on_off_conserved1 | on_off_conserved2 | on_off_conserved3,

		coexpressed1 = coexpressed1, 
		coexpressed2 = coexpressed2, 
		dosage_coexpressed = dosage_coexpressed,
		double_coexpressed = double_coexpressed, 
		single_coexpressed = single_coexpressed,
		any_coexpressed    = any_coexpressed   ,

		cor_subfunc       = f4 & r12_lc & r03_hc & r01<r03 & r02<r03,
		cor_weak1_subfunc = f4 &          r03_hc & r01<r03 & r02<r03,
		cor_weak2_subfunc = f4 & r12_lc & r03_mc & r01<r03 & r02<r03,
		cor_weak3_subfunc = f4 &          r03_mc & r01<r03 & r02<r03,
		on_off_subfunc = on_off_subfunc,
		on_off_part_subfunc = on_off_part_subfunc,

		cor_neofunc1 = f4 & r12_lc & (r02_hc & r01_lc),
		cor_neofunc2 = f4 & r12_lc & (r01_hc & r02_lc),
		cor_neofunc = f4 & r12_lc & (r01_hc & r02_lc | r02_hc & r01_lc),
		cor_weak_neofunc1 = f4 & r12_lc & (r02_mc & r01_lc & n_gt10>0),
		cor_weak_neofunc2 = f4 & r12_lc & (r01_mc & r02_lc & n_gt20>0),
		cor_weak_neofunc  = f4 & r12_lc & (r01_mc & r02_lc & n_gt20>0 | r02_mc & r01_lc & n_gt10>0),
		on_off_neofunc1 = on_off_neofunc1,
		on_off_neofunc2 = on_off_neofunc2,
		on_off_neofunc  = on_off_neofunc ,
		on_off_part_neofunc1 = on_off_part_neofunc1, 
		on_off_part_neofunc2 = on_off_part_neofunc2, 
		on_off_part_neofunc  = on_off_part_neofunc , 
		on_off_neofunc_WGD1 = on_off_neofunc_WGD1,
		on_off_neofunc_WGD2 = on_off_neofunc_WGD2,
		on_off_neofunc_WGD  = on_off_neofunc_WGD ,
		on_off_part_neofunc_WGD1 = on_off_part_neofunc_WGD1,
		on_off_part_neofunc_WGD2 = on_off_part_neofunc_WGD2,
		on_off_part_neofunc_WGD  = on_off_part_neofunc_WGD ,

		nonfunc1= nonfunc1,
		nonfunc2= nonfunc2,
		nonfunc = nonfunc ,
		nonfunc_WGD1= nonfunc_WGD1,
		nonfunc_WGD2= nonfunc_WGD2,
		nonfunc_WGD = nonfunc_WGD ,
		part_nonfunc1 = part_nonfunc1,
		part_nonfunc2 = part_nonfunc2, 
		part_nonfunc  = part_nonfunc ,
		part_nonfunc_WGD1 = part_nonfunc_WGD1,
		part_nonfunc_WGD2 = part_nonfunc_WGD2, 
		part_nonfunc_WGD  = part_nonfunc_WGD ,

		greater1 = tp01<tp & mean1>mean0,
		greater2 = tp02<tp & mean2>mean0,
		greater3 = tp03<tp & mean3>mean0,
		less1    = tp01<tp & mean1<mean0,
		less2    = tp02<tp & mean2<mean0,
		less3    = tp03<tp & mean3<mean0,
		on_off_greater1 = n_gt01==0 & n_gt10>0,
		on_off_greater2 = n_gt02==0 & n_gt20>0,
		on_off_greater3 = n_gt03==0 & n_gt30>0,
		on_off_less1    = n_gt01>0 & n_gt10==0,
		on_off_less2    = n_gt02>0 & n_gt20==0,
		on_off_less3    = n_gt03>0 & n_gt30==0,

		all_greater1 = tp01<tp & apply(a0<a1,1,sum)==nts,
		all_greater2 = tp02<tp & apply(a0<a2,1,sum)==nts,
		all_greater3 = tp03<tp & apply(a0<a3,1,sum)==nts,
		all_less1 = tp01<tp & apply(a0>a1,1,sum)==nts,
		all_less2 = tp02<tp & apply(a0>a2,1,sum)==nts,
		all_less3 = tp03<tp & apply(a0>a3,1,sum)==nts
	);
	res
}

is_strict_neofunc = function(data1, i0, i1, i2, j0, j1, j2) {
	a0=data1[,(nts*(i0-1)+1):(nts*i0)];
	a1=data1[,(nts*(i1-1)+1):(nts*i1)];
	a2=data1[,(nts*(i2-1)+1):(nts*i2)];
	b0=data1[,(nts*(j0-1)+1):(nts*j0)];
	b1=data1[,(nts*(j1-1)+1):(nts*j1)];
	b2=data1[,(nts*(j2-1)+1):(nts*j2)];
	a_max12 = a1;
	a_max12[ which(a1<a2, arr=T) ] = a2[which(a1<a2, arr=T)];
	b_max12 = b1;
	b_max12[ which(b1<b2, arr=T) ] = b2[which(b1<b2, arr=T)];
	a_max0 = apply(a0, 1, max);
	a_max1 = apply(a1, 1, max);
	a_max2 = apply(a2, 1, max);
	b_max0 = apply(b0, 1, max);
	b_max1 = apply(b1, 1, max);
	b_max2 = apply(b2, 1, max);
	min_max_ab = apply(cbind(a_max0,a_max1,a_max2,b_max0, b_max1,b_max2), 1, min);
	min_max_ab>=1 & apply(a0<0.1 & b0<0.1 & (a_max12>=1 & b_max12>=1), 1, sum)>0
}

is_strict_subfunc = function(data1, i0, i1, i2, j0, j1, j2) {
	a0=data1[,(nts*(i0-1)+1):(nts*i0)];
	a1=data1[,(nts*(i1-1)+1):(nts*i1)];
	a2=data1[,(nts*(i2-1)+1):(nts*i2)];
	b0=data1[,(nts*(j0-1)+1):(nts*j0)];
	b1=data1[,(nts*(j1-1)+1):(nts*j1)];
	b2=data1[,(nts*(j2-1)+1):(nts*j2)];
	a_max0 = apply(a0, 1, max);
	a_max1 = apply(a1, 1, max);
	a_max2 = apply(a2, 1, max);
	b_max0 = apply(b0, 1, max);
	b_max1 = apply(b1, 1, max);
	b_max2 = apply(b2, 1, max);
	min_max_ab = apply(cbind(a_max0,a_max1,a_max2,b_max0, b_max1,b_max2), 1, min);
	n1 = apply(a0>=1 & b0>=1 & a1>=1 & a2<1 & a1-a2>=1 & b1>=1 & b2<1 & b1-b2>=1, 1, sum)
	n2 = apply(a0>=1 & b0>=1 & a1>=1 & a2<1 & a1-a2>=1 & b2>=1 & b1<1 & b2-b1>=1, 1, sum)
	n3 = apply(a0>=1 & b0>=1 & a2>=1 & a1<1 & a2-a1>=1 & b1>=1 & b2<1 & b1-b2>=1, 1, sum)
	n4 = apply(a0>=1 & b0>=1 & a2>=1 & a1<1 & a2-a1>=1 & b2>=1 & b1<1 & b2-b1>=1, 1, sum)
	min_max_ab>=1 & (n1>0 & n4>0 | n2>0 & n3>0)
}


subfunc_score = function(data1, i0, i1, i2, thres) {
	a0=data1[,(nts*(i0-1)+1):(nts*i0)];
	a1=data1[,(nts*(i1-1)+1):(nts*i1)];
	a2=data1[,(nts*(i2-1)+1):(nts*i2)];
	a3 = a1+a2;
	r01 = pair_cor(a0,a1);   r01.p = pair_cor_p(a0,a1);
	r02 = pair_cor(a0,a2);   r02.p = pair_cor_p(a0,a2);
	r03 = pair_cor(a0,a3);   # r03.p = pair_cor_p(a0,a3);
	r12 = pair_cor(a1,a2);   # r12.p = pair_cor_p(a1,a2);
	max0 = apply(a0, 1, max);
	max1 = apply(a1, 1, max);
	max2 = apply(a2, 1, max);
	min_max12 = apply(cbind(max1,max2),1,min);
#	r0max = max(r01,r02);
#	r0min = min(r01,r02);
	r0min = apply(cbind(r01,r02), 1, min);
	r = apply(cbind(r0min,r12),1,max);
	x = rep(0,length(r03));
	x = sqrt(abs(r03*(r03-r)/2));
	x[r03<0 | r03<r] = -x[r03<0 | r03<r];
	x
}

neofunc_score = function(data, i0, i1, i2) {
	a0=data[,(nts*(i0-1)+1):(nts*i0)];
	a1=data[,(nts*(i1-1)+1):(nts*i1)];
	a2=data[,(nts*(i2-1)+1):(nts*i2)];
#	a3=data[,(nts*(i3-1)+1):(nts*i3)];
	r01 = pair_cor(a0,a1);   # r01.p = pair_cor_p(a0,a1);
	r02 = pair_cor(a0,a2);   # r02.p = pair_cor_p(a0,a2);
#	r12 = pair_cor(a1,a2);   # r12.p = pair_cor_p(a1,a2);
#	f12= r12<C1
	r0max = apply(cbind(r01,r02),1,max);
	r0min = apply(cbind(r01,r02),1,min);
	x = rep(0,length(r01));
	x = sqrt(abs(r0max*(r0max-r0min)/2));
	x[r0max<0] = -x[r0max<0];
	x
}


is_tissue_subfunc = function(data, i0, i1, i2, thres, thres1, thres2, delta) {
	n=tot_copy+tot_sp; 
	a0=data[,(nts*(i0-1)+1):(nts*i0)];
	a1=data[,(nts*(i1-1)+1):(nts*i1)];
	a2=data[,(nts*(i2-1)+1):(nts*i2)];
	a3=a1+a2;
	max0 = apply(a0, 1, max);
	max1 = apply(a1, 1, max);
	max2 = apply(a2, 1, max);
	f4 = max0>=thres & max1>=thres & max2>=thres;
	f1 = apply( a1<thres1 & a0>=thres2 & a2>=thres2 & a0-a1>=delta & a2-a1>=delta, 1, sum) > 0;
	f2 = apply( a2<thres1 & a0>=thres2 & a1>=thres2 & a0-a2>=delta & a1-a2>=delta, 1, sum) > 0;
	list(res=f4&f1&f2, res1=f4&f1, res2=f4&f2);
}

is_tissue_complete_subfunc = function(data, i0, i1, i2, i3, thres, thres1, thres2, delta) {
	n=tot_copy+tot_sp; 
	a0=data[,(nts*(i0-1)+1):(nts*i0)];
	a1=data[,(nts*(i1-1)+1):(nts*i1)];
	a2=data[,(nts*(i2-1)+1):(nts*i2)];
	a3=data[,(nts*(i3-1)+1):(nts*i3)];
	max0 = apply(a0, 1, max);
	max1 = apply(a1, 1, max);
	max2 = apply(a2, 1, max);
	f4 = max0>=thres & max1>=thres & max2>=thres;
	n1 = apply( a1<thres1 & a0>=thres2 & a2>=thres2 & a0-a1>=delta & a2-a1>=delta, 1, sum);
	n2 = apply( a2<thres1 & a0>=thres2 & a1>=thres2 & a0-a2>=delta & a1-a2>=delta, 1, sum);
	n  = apply(a0>=thres2, 1, sum);
	f = n1>0 & n2>0 & (n1+n2)==n;
	list(res=f4&f);
}

is_tissue_neofunc = function(data, i0, i1, i2, i3, thres=1, thres1=1, thres2=1, delta=1) {
	a0=data[,(nts*(i0-1)+1):(nts*i0)];
	a1=data[,(nts*(i1-1)+1):(nts*i1)];
	a2=data[,(nts*(i2-1)+1):(nts*i2)];
	a3=data[,(nts*(i3-1)+1):(nts*i3)];
	max0 = apply(a0, 1, max);
	max1 = apply(a1, 1, max);
	max2 = apply(a2, 1, max);
	f4 = max0>=thres & max1>=thres & max2>=thres;
	f1 = apply( a0<thres1 & a2<thres1 & a1>=thres2 & a1-a0>=delta & a1-a2>=delta, 1, sum) > 0;
#	f1b = apply(a0>=thres1 & a2>=thres1, 1, sum) == apply(a0<=thres1, 1, sum);
	f2 = apply( a0<thres1 & a1<thres1 & a2>=thres2 & a2-a0>=delta & a2-a1>=delta, 1, sum) > 0;
#	f2b = apply(a0>=thres1 & a1>=thres1, 1, sum) == apply(a0<=thres1, 1, sum);
	res1 = f4 & f1;
	res2 = f4 & f2;
	list( res = res1|res2, res1=res1, res2=res2 );
}

is_nonfunc = function(data, i0, i1, i2, i3, thres=1, thres1=1, delta=1) {
	a0=data[,(nts*(i0-1)+1):(nts*i0)];
	a1=data[,(nts*(i1-1)+1):(nts*i1)];
	a2=data[,(nts*(i2-1)+1):(nts*i2)];
	a3=data[,(nts*(i3-1)+1):(nts*i3)];
	max0 = apply(a0, 1, max);
	max1 = apply(a1, 1, max);
	max2 = apply(a2, 1, max);
	f1 = max0>=thres & max1<thres1 & max0-max1>=delta;
	f2 = max0>=thres & max2<thres1 & max0-max2>=delta;
	res2= f2;
	res1= f1;
	list(res=res1|res2, res1=res1, res2=res2);
}

is_dosage_balance = function(data, i0, i1, i2, i3, thres, tp, C1=0.6) {
	a0=data[,(nts*(i0-1)+1):(nts*i0)];
	a3=data[,(nts*(i3-1)+1):(nts*i3)];
	r03 = pair_cor(a0,a3);   r03.p = pair_cor_p(a0,a3);
	t03 = pair_t_p(a0, a3);
	max0 = apply(data[,(nts*(i0-1)+1):(nts*i0)], 1, max);
	max1 = apply(data[,(nts*(i1-1)+1):(nts*i1)], 1, max);
	max2 = apply(data[,(nts*(i2-1)+1):(nts*i2)], 1, max);
	f3 = r03 >=C1 & t03>=tp; # ZF--(CC1+CC2) >= 0.6
	list(res=f3 & max0>thres &  max1>thres & max2>thres);
}

# count sub-func
# Definition in Braasch 2014
# C1=0.75, C2=0.75
# cor(ZF,(CC1+CC2))>=C2
# cor(CC1,CC2)<C1
# ( cor(ZF,CC1)<cor(ZF,(CC1+CC2)) AND cor(ZF,CC2)<cor(ZF,(CC1+CC2)) )
# 
# My definiation:
# C1=0.6, C2=0.75
# cor(ZF,(CC1+CC2))>=C2, 
# cor(CC1,CC2) < C1,
# ( cor(ZF,CC1)<cor(ZF,(CC1+CC2)) AND cor(ZF,CC2)<cor(ZF,(CC1+CC2)) )
# # cor(ZF,CC1)<0.6 & cor(ZF,CC2)<0.6
# # cor(CC1,CC2)<cor(ZF,CC1), cor(CC1,CC2)<cor(ZF,CC2)
# # cor(ZF,CC1)<cor(ZF,(CC1+CC2))
# # cor(ZF,CC2)<cor(ZF,(CC1+CC2))
# # cor(CC1,CC2) < cor(ZF,(CC1+CC2))
# # all ohnologs is expressed in at least one tissue
is_subfunc= function(data1, i0, i1, i2, i3, thres1=1, C1=0.6, C2=0.75) {
	a0=data1[,(nts*(i0-1)+1):(nts*i0)];
	a1=data1[,(nts*(i1-1)+1):(nts*i1)];
	a2=data1[,(nts*(i2-1)+1):(nts*i2)];
	a3=data1[,(nts*(i3-1)+1):(nts*i3)];
	r01 = pair_cor(a0,a1);   r01.p = pair_cor_p(a0,a1);
	r02 = pair_cor(a0,a2);   r02.p = pair_cor_p(a0,a2);
	r03 = pair_cor(a0,a3);   r03.p = pair_cor_p(a0,a3);
	r12 = pair_cor(a1,a2);   r12.p = pair_cor_p(a1,a2);
	max0 = apply(a0, 1, max);
	max1 = apply(a1, 1, max);
	max2 = apply(a2, 1, max);
	f3 = r03 >=C2 & r03.p<0.1; # ZF--(CC1+CC2) >= C2
	f12= r12<C1;  # cor(CC1,CC2) <C1
	f0 = r01<r03 & r02<r03; # cor(ZF,CC1)<cor(ZF,(CC1+CC2)) AND cor(ZF,CC2)<cor(ZF,(CC1+CC2))
	res= f3 & f12 & f0 & max0>=thres & max1>=thres & max2>=thres;
	list(res=res);
}


# Definition in Braasch 2014
# C1=0.75, C2=0.75
# cor(CC1,CC2)<C1
# (cor(ZF,CC1)<C1 & cor(ZF,CC2)>=C1) or (cor(ZF,CC2)<C1 & cor(ZF,CC1)>=C1)
is_neofunc = function(data, i0, i1, i2, i3, thres, C1,C2) {
	n=tot_copy+tot_sp; 
	a0=data[,(nts*(i0-1)+1):(nts*i0)];
	a1=data[,(nts*(i1-1)+1):(nts*i1)];
	a2=data[,(nts*(i2-1)+1):(nts*i2)];
	a3=data[,(nts*(i3-1)+1):(nts*i3)];
	r01 = pair_cor(a0,a1);   r01.p = pair_cor_p(a0,a1);
	r02 = pair_cor(a0,a2);   r02.p = pair_cor_p(a0,a2);
	r03 = pair_cor(a0,a3);   #r03.p = pair_cor_p(a0,a3);
	r12 = pair_cor(a1,a2);   #r12.p = pair_cor_p(a1,a2);
	max0 = apply(a0, 1, max);
	max1 = apply(a1, 1, max);
	max2 = apply(a2, 1, max);
	f4 = max0>=thres & max1>=thres & max2>=thres;
	f12= r12<C1;
	f1 = r02>=C2 & r02.p<0.1 & r01<C1;
	f2 = r01>=C2 & r01.p<0.1 & r02<C1;
	res1=f4&f12&f1;
	res2=f4&f12&f2;
	list(res=res1|res2, res1=res1, res2=res2);
}

# go_test
go_test = function(dataset, selected, all)
{
	n = length(all);
	m = length(selected);
	p = sapply(dataset, FUN=function(x) {
			n2=sum(x %in% all)
			if (n2<10 || n2>500) {
				1.0;
			} else {
				n11 = sum(x %in% selected);
				n12 = n2 - n11;
				a = matrix(c(n11,m-n11, n12, n-m-n12), 2,2);
				fisher.test(a)$p.value;
			}
		});
	r = sapply(dataset, FUN=function(x) {
			n2=sum(x %in% all)
			if (n2<10 || n2>500) {
				1.0;
		   	} else {
				n11 = sum(x %in% selected);
				n12 = n2 - n11;
				a = matrix(c(n11,m-n11, n12, n-m-n12), 2,2);
				fisher.test(a)$estimate;
			}
		});
	fdr = p.adjust(p, method='fdr');
	cbind(p,ratio=r, fdr=fdr);
}

random_data=function(data, by_row=F, replace=F)
{
	m=nrow(data);
	n=ncol(data);
	
	res = data;
	if (by_row) {
		for (i in 1:m) {
			res[i,]=sample(data[i,], n, replace=replace);
		}
	} else {
		for (j in 1:n) {
			res[,j]=sample(data[,j], m, replace=replace);
		}
	}
	res
}

random_data2=function(data, win, by_row=F, replace=F)
{
	m=nrow(data);
	n=ncol(data);
	
	res = data;
	if (by_row) {
		for (i in 1:(m/win)) {
			idx = sample(1:n, n, replace=replace);
			res[(win*(i-1)+1):(win*i),]=data[(win*(i-1)+1):(win*i), idx];
		}
	} else {
		for (j in 1:(n/win)) {
			idx = sample(1:m, m, replace=replace);
			res[,(win*(j-1)+1):(win*j)]=data[idx,(win*(j-1)+1):(win*j)];
		}
	}
	res
}

ohno_sum=function(data, tissues, sps)
{
	nts = length(tissues);
	tot_sp = length(levels(as.factor(sps)));
	sps1 = unique(sps);
	sum_data = matrix(0, nrow(data), tot_sp*nts);
	colnames(sum_data) = paste0("sum_of_", rep(sps1,each=nts), '_', rep(tissues,tot_sp));
	for (i in 1:tot_sp) {
		js = (1:length(sps))[sps==sps1[i]];
		if (length(js)==1) {
			j0 = js[1];
			sum_data[,(nts*(i-1)+1):(i*nts)] = data[,(nts*(j0-1)+1):(j0*nts)];
		} else {
			for (j in 1:nts) {
				sum_data[,nts*(i-1)+j] = apply(data[,(nts*(js-1)+j)], 1,sum);
			}
		}
	}
	sum_data
}

#! barplot with sd bar
barplot_sd <- function(height, sd, bar_len_inch=0.2, ...)
# {{{
{
	h1 = height+sd;
	h0 = height-sd;
	center<-barplot(height, ...);
	arrows(center, h0, center, h1, code=3, angle=90, length=bar_len_inch/2, xpd=T );
}
# }}}

#! barplot with sd bar
barplot_sd1 <- function(height, sd, bar_len_inch=0.2, ...)
# {{{
{
	h1 = height+sd;
	h0 = height-sd;
	center<-barplot(height, ...);
	arrows(center, height, center, h1, code=2, angle=90, length=bar_len_inch/2, xpd=T );
}
# }}}

plot_mod_pattern <- function(data, genes, sm_col='black')
# {{{
{
	n <- ncol(data);
	x <- data[genes,];
	if (is.vector(x)) { x <- matrix(x, 1, n); }
	m <- nrow(x);
	x.sm_mean <- apply(x, 2, mean);
	y0 <- min( apply(x, 2, function(x) {quantile(x,0.01)} ) );
	y1 <- max( apply(x, 2, function(x) {quantile(x,0.99)} ) );
	plot(1:n, x.sm_mean, col=sm_col, ylim=c(y0,y1), cex=0.5, xlab="Sample", ylab="Expression", axes=F);
	axis(2);
	for (i in 1:m) {
		lines(1:n, x[i,], col=labels2colors(1:20)[(i %% 20)+1]);
	}
	points(1:n, x.sm_mean, col=sm_col, ylim=c(y0,y1));
}
# }}}

evalue_cluster <- function(data, cid)
# {{{
{
    n1 <- nrow(data);   n2 <- ncol(data);
    csz <- tapply(rep(1,length(cid)), cid, sum);
    k <- length(csz); # number of cluster
    dist_DI <- 0; # Dunn Index with delta=mean distance between mean and all points
	corr_DI <- 0;
    WSS <- 0;   BSS <- 0;   TSS <- sum( apply(data,2,var)*(n1-1) );
    mean_max_within_dist <- 0;
    mean_min_within_corr <- 0;
    min_within_corr <- rep(0,k);
    max_within_dist <- rep(0,k);
    means0 <- apply(data, 2, mean);
    means <- matrix(0, k, ncol(data2));
    wsss <- matrix(0, k, n2);
    bsss <- matrix(0, k, n2);
    for (kk in 1:k) {
        for (m in 1:n2) {
            means[kk,m] <- mean(data2[cid==kk,m]);
            if (csz[kk]==1) { wsss[kk,m] <- 0; }
            else { wsss[kk,m] <- (csz[kk]-1)*var(data2[cid==kk,m]); }
            bsss[kk,m] <- csz[kk]*(means[kk,m]-means0[m])^2;
        }
    }
    WSS <- sum( apply(wsss, 2, sum) );
    BSS <- sum( apply(bsss, 2, sum) );

    for (kk in 1:k) {
        min_within_corr[kk] <- min( cor(means[kk,], t(data)) );
        max_within_dist[kk] <- max( apply( (t(data[cid==k])-means[kk,])^2, 1, sum ) );
    }
    mean_min_within_corr <- mean( min_within_corr );
    mean_max_within_dist <- mean( max_within_dist );

    a1 <- c();   a2 <- c();
    for (k1 in 1:(k-1)) {
        for (k2 in (k1+1):k) {
            a1 = c(a1, sum( (means[k1,]-means[k2,])^2 ) );
            a2 = c(a2, cor(means[k1,],means[k2,]) );
        }
    }
    min_between_dist <- min(a1);
    max_between_corr <- max(a2);
    dist_DI = min_between_dist / mean_max_within_dist;
    corr_DI = max_between_corr / mean_min_within_corr;

    lm.aics <- rep(0, n2);   lm.bics <- rep(0, n2);
    for (m in 1:n2) {
        lm.res <- lm(data[,m] ~ factor(cid));
        lm.aics[m] <- AIC(lm.res);
        lm.bics[m] <- BIC(lm.res);
    }
    lm.aic <- mean(lm.aics);
    lm.bic <- mean(lm.bics);

    list(WSS=WSS, BSS=BSS, TSS=TSS, 
        mean_min_within_corr=mean_min_within_corr, 
        mean_max_within_dist=mean_max_within_dist,
        min_between_dist = min_between_dist,   max_between_corr = max_between_corr,
        aic=lm.aic,  bic=lm.bic,  dist_DI = dist_DI, corr_DI=corr_DI);
}
# }}}

# max_nc: maximal number of cluster
choose_cluster_num_WSS <- function(data, h, max_nc, delta=0.01, frac=0.01)
# {{{
{
	m = nrow(data);
	n = ncol(data);
	TSS <- sum( apply(data,2,var)*(m-1) );
	WSS <- rep(0,max_nc);
	WSS[1] = 1.0;
	delta2 = 0;
	nc0=0;
	for (nc in 2:max_nc) {
		cid = cutree(h, k=nc);
		wss = 0;
		for (j in 1:n) {
			lm.res <- lm(data[,j] ~ factor(cid));
			aov.res <- aov(lm.res);
			wss = c(wss, sum(aov.res$residuals^2));
		}
		WSS[nc] = sum(wss)/TSS;
		if (nc==2) {delta2=WSS[nc-1]-WSS[nc];}
		if (WSS[nc-1]-WSS[nc] < delta && WSS[nc-1]-WSS[nc]<delta2*frac ) { if (nc0==0) {nc0 = nc;} }
	}
	plot(2:nc, WSS[2:nc], xlab="NC", ylab="WSS");
	gc();
	list(nc=nc0, WSS=WSS);
}
# }}}

choose_cluster_num_BIC<- function(data, h, max_nc)
# {{{
{
	m = nrow(data);
	n = ncol(data);
	BIC <- rep(0,max_nc);
	nc0=0;
	for (nc in 2:max_nc) {
		cid = cutree(h, k=nc);
		bic = c();
		for (j in 1:n) {
			lm.res <- lm(data[,j] ~ 0+factor(cid));
			bic = c(bic, BIC(lm.res) );
		}
		BIC[nc] = sum(bic);
	}
	nc0 = max.col(t(-BIC[2:max_nc]))+1;
	plot(2:nc, BIC[2:nc], xlab="NC", ylab="BIC");
	gc();
	list(BIC=BIC, nc=nc0);
}
# }}}

choose_cluster_num_silhouette<- function(data, h, max_nc)
# {{{
{
	m = nrow(data);
	n = ncol(data);
	v <- rep(0,max_nc);
	nc0=0;
	for (nc in 2:max_nc) {
		cid = cutree(h, k=nc);
		sil = silhouette(cid, (1-cor(t(data)))/2);
		v[nc] = min(summary(sil)$avg.width);
	}
	nc0 = max.col(t(v[2:max_nc]))+1;
	plot(2:nc, v[2:nc], xlab="NC", ylab="mean silhouette");
	gc();
	nc0;
}
# }}}

plot_gene_exp <- function(data.group_mean, data.group_sd, id, col='grey', ...)
# {{{
{
	ymax = max(data.group_mean[id,] + data.group_sd[id,]);
	ymin = min(c(0,data.group_mean[id,] - data.group_sd[id,]));
	center<-barplot(data.group_mean[id,], ylim=c(ymin,ymax), col=col, ...);
	arrows(center, data.group_mean[id,]-data.group_sd[id,], center, data.group_mean[id,]+data.group_sd[id,], code=3, angle=90, length=0.3, xpd=T );
}
# }}}

plot_mod_pattern <- function(data, genes, sm_col='black')
# {{{
{
	n <- ncol(data);
	x <- data[genes,];
	if (is.vector(x)) { x <- matrix(x, 1, n); }
	m <- nrow(x);
	x.sm_mean <- apply(x, 2, mean);
	y0 <- min( apply(x, 2, function(x) {quantile(x,0.01)} ) );
	y1 <- max( apply(x, 2, function(x) {quantile(x,0.99)} ) );
	plot(1:n, x.sm_mean, col=sm_col, ylim=c(y0,y1), cex=0.5, xlab="Sample", ylab="Expression", axes=F);
	axis(2);
	for (i in 1:m) {
		lines(1:n, x[i,], col=labels2colors(1:20)[(i %% 20)+1]);
	}
	points(1:n, x.sm_mean, col=sm_col, ylim=c(y0,y1));
}
# }}}

mat_degree <- function(mat)
# {{{
{
	n <- ncol(mat);
	m <- nrow(mat);
	if (n!=m) { return(c()); }
	for (i in 1:n) { mat[i,i]=0; }
	degree <- apply(mat>0,1,sum);
	degree;
}
# }}}

mod_jaccard <- function(data, cid1, cid2)
# {{{
{
	nc1 = max(cid1);
	nc2 = max(cid2);
	csz1 = table(cid1);
	csz2 = table(cid2);
	mat = matrix(0,nc1,nc2);
	for (i in 1:length(cid1)) {
		i1=cid1[i];
		i2=cid2[i];
		mat[i1,i2] = mat[i1,i2] + 1;
	}
	for (i1 in 1:nc1) {
		for (i2 in 1:nc2) {
			mat[i1,i2] = mat[i1,i2] / (csz1[i1]+csz2[i2]-mat[i1,i2]);
		}
	}
	mat
}
# }}}

mod_frac <- function(data, cid1, cid2)
# {{{
{
	nc1 = max(cid1);
	nc2 = max(cid2);
	csz1 = table(cid1);
	csz2 = table(cid2);
	mat = matrix(0,nc1,nc2);
	for (i in 1:length(cid1)) {
		i1=cid1[i];
		i2=cid2[i];
		mat[i1,i2] = mat[i1,i2] + 1;
	}
	for (i1 in 1:nc1) {
		for (i2 in 1:nc2) {
			mat[i1,i2] = mat[i1,i2] / min(csz1[i1],csz2[i2]);
		}
	}
	mat
}
# }}}

mod_share <- function(data, cid1, cid2)
{
	nc1 = max(cid1);
	nc2 = max(cid2);
	mat = matrix(0,nc1,nc2);
	apply( cbind(cid1,cid2), 1, function(x) { mat[x]= mat[x]+1; } );
	mat
}

mod_share <- function(data, cid1, cid2)
{
	nc1 = max(cid1);
	nc2 = max(cid2);
	mat = matrix(0,nc1,nc2);
	apply( cbind(cid1,cid2), 1, function(x) { mat[x]= mat[x]+1; } );
	mat
}

mod_cor1 <- function(data, cid1, cid2)
{
	nc1 = max(cid1);
	nc2 = max(cid2);
	data.cor = cor(t(data));
	cor1 = matrix(0,nc1, nc2);
	for (i1 in 1:nc1) {
		for (i2 in 1:nc2) {
			cor1[i1,i2] = mean( data.cor[cid1==i1,cid2==i2] );
		}
	}
	cor1
}

mod_cor2 <- function(data, cid1, cid2)
{
	nc1 = max(cid1);
	nc2 = max(cid2);
	cor1 = matrix(0,nc1, nc2);
	for (i1 in 1:nc1) {
		for (i2 in 1:nc2) {
			x1 = data[cid1==i1];
			x2 = data[cid2==i2];
			c1 = apply(x1, 2, mean);
			c2 = apply(x2, 2, mean);
			cor1[i1,i2] = cor(c1,c2);
		}
	}
	cor1
}


# data: expression matrix, each row is a gene. gene name is entrez ID
# data_eg_id: vector, Entrez ID for each gene in the data, same order as data
# gene_kegg: a vector, gene to kegg mapping , names is gene entrez ID
plot_kegg_heatmap <- function(data, log2FC_data, PValue_data, data_eg_id=NULL, gene_kegg, path_id, out, gene_symbol=NULL, gene_lab_col='black', species='dre')
{
	if (is.null(data_eg_id)) { data_eg_id = rownames(data); }
	eg_ids = names(gene_kegg)[gene_kegg==path_id];
	eg_ids = eg_ids[ eg_ids %in% data_eg_id ];
	f = data_eg_id %in% eg_ids & apply(!is.na(log2FC_data) & log2FC_data!=0,1,sum)>1;
	f1 = f;
# for multiple entrez ID, choose the first
	c1 = rep(0,length(eg_ids));
	names(c1) = eg_ids;
	for (i in 1:length(f)) {
		if (!f[i]) {next;}
		if (c1[ data_eg_id[i] ]==0) { c1[ data_eg_id[i] ]=1; }
		else { f1[ i ] = F; }
	}
	eg_ids = data_eg_id[f];
	eg_ids1 = data_eg_id[f1];
##
	if (is.null(gene_symbol)) {
		ylab = paste0( rownames(data)[f], '|',  mapIds(org.Dr.eg.db, keys=eg_ids, keytype='ENTREZID', column=c('SYMBOL'), multiVals='first') );
	} else {
		ylab = gene_symbol[ f ];
	}
	x2 = log2FC_data[f1,];  rownames(x2) = eg_ids1;
	idtype='KEGG';
	same_layer=T;
	if (species=='ko') { idtype='KEGG'; }
	else {idtyype="entrez"; same_layer=F;}
	pathview(gene.data = x2, pathway.id = path_id, species = species, gene.idtype=idtype, kegg.native=T, same.layer=same_layer, map.symbol=T, limit=list(gene=c(-2,2), cpd=c(-2,2)), bins=list(gene=20, cpd=20));
	file=paste0(path_id, ".png");   file.remove(file);
	file=paste0(path_id, ".xml");   file.remove(file);
	file=paste0(path_id, ".pathview.multi.png");
	if (file.exists(file))  {
		file.rename( file, paste0(out, path_id, '.multi.png') );
	}
	file=paste0(path_id, "pathview.png");
	if (file.exists(file))  {
		file.rename( file, paste0(out, path_id, '.png') );
	}

	pdf( paste0(out, path_id, '.heatmap.pdf'), wi=10, he=12+0.05*sum(f));
	x = data[f,];
	z.max = max(x,-x);
#	labeledHeatmap(x, xLabels=colnames(x), yLabels=ylab, colors=bluered(255), zlim=c(-z.max, z.max), x.adj.lab.y = 0, main='Expression' );
	m=c(1.5,3)/par('csi');
	heatmap.2(x, trace='no', col=bluered(255), zlim=c(-z.max,z.max), labRow = ylab, margins=m, Colv='none', dendrogram='row', dist=cor_dist, main='Normalized Expression', lhei=c(2,10), colRow=gene_lab_col[f])
	x = log2FC_data[f,];
	x[ which(is.na(x), arr=T) ] = 0;
	z.max = max(x,-x);
#	labeledHeatmap(x, xLabels=colnames(x), yLabels=ylab, colors=bluered(255), zlim=c(-z.max, z.max), x.adj.lab.y = 0, main='log2FC');
	if (!missing(PValue_data) && !is.null(PValue_data)) {
		p = PValue_data[f,];
		p[ which(is.na(p), arr=T) ] = 1;
		p1=p;
		p1[ which(p>=0.05, arr=T) ] = '';
		p1[ which(p<0.05, arr=T) ] = '*';
		p1[ which(p<0.01, arr=T) ] = '**';
		p1[ which(p<0.001, arr=T) ] = '***';
		heatmap.2(x, trace='no', col=bluered(255), zlim=c(-z.max,z.max), labRow = ylab, margins=m, Colv='none', dendrogram='row', dist=cor_dist, main='Normalized Expression', lhei=c(1,10), colRow=gene_lab_col[f], cellnote=p1, notecex=0.8, notecol='black' )
	} else {
		heatmap.2(x, trace='no', col=bluered(255), zlim=c(-z.max,z.max), labRow = ylab, margins=m, Colv='none', dendrogram='row', dist=cor_dist, main='Normalized Expression', lhei=c(1,10), colRow=gene_lab_col[f] )
	}
	dev.off();
	0;
}


# mat: a correlation or distance matrix
# h  : hclust results
# cid: cluster assignment
# select: selected data rows (columns), sorted.
plot_image <- function(mat, h, cid, select, col=heat.colors(255))
# {{{
{
	n = ncol(mat);
	n1 = length(select);
	cid1 = cid[h$order][select];
	grid1 = (1:n1)[ diff(cid1)!=0 ];
	x = mat[h$order,h$order][select,select];
	image(0:n1, 0:n1, x, col=col);
	for (i in 1:length(grid1)) {
		y = grid1[i];
		abline(h=y, col='blue');
		abline(v=y, col='blue');
	}
}
# }}}

plot_link <- function(x0, y0, x1, y1, r1, r2, side, ...)
{
	tmp_f=function(x, cx, r1) { a1=1-((x-cx)/r1)^2; a1[a1<0]=0; a1 }
	if (side==1) { # bottom
		if (2*r1>x1-x0) { r1=(x1-x0)/2; }
		if (r2>y1-y0) { r2=(y1-y0); }
		lines(c(x0+r1,x1-r1), c(y0,y0), ...);
		lines(c(x0,x0), c(y0+r2,y1), ...);
		lines(c(x1,x1), c(y0+r2,y1), ...);
		cx = x0+r1;
		cy = y0+r2;
		if (r1>0 && r2>0) {
			curve( cy-r2*sqrt(1-tmp_f(x,cx,r1)), xlim=c(x0,x0+r1), add=T, ...);
			cx = x1-r1;
			curve( cy-r2*sqrt(1-tmp_f(x,cx,r1)), xlim=c(x1-r1,x1), add=T, ...);
		}
	} else if (side==2) { # left
		if (r1>x1-x0) { r1=x1-x0; }
		if (2*r2>y1-y0) { r2=(y1-y0)/2; }
		lines(c(x0+r1,x1), c(y0,y0), ...);
		lines(c(x0+r1,x1), c(y1,y1), ...);
		lines(c(x0,x0), c(y0+r2,y1-r2), ...);
		if (r1>0 && r2>0) {
			cx = x0+r1;
			cy = y0+r2;
			curve( cy-r2*tmp_f(x,cx,r1), xlim=c(x0,x0+r1), add=T, ...);
			cy = y1-r2;
			curve( cy+r2*tmp_f(x,cx,r1), xlim=c(x0,x0+r1), add=T, ...);
		}
	} else if (side==3) { # top
		if (2*r1>x1-x0) { r1=(x1-x0)/2; }
		if (r2>y1-y0) { r2=(y1-y0); }
		lines(c(x0+r1,x1-r1), c(y1,y1), ...);
		lines(c(x0,x0), c(y0,y1-r2), ...);
		lines(c(x1,x1), c(y0,y1-r2), ...);
		if (r1>0 && r2>0) {
			cx = x0+r1;
			cy = y1-r2;
			curve( cy+r2*tmp_f(x,cx,r1), xlim=c(x0,x0+r1), add=T, ...);
			cx = x1-r1;
			curve( cy+r2*tmp_f(x,cx,r1), xlim=c(x1-r1,x1), add=T, ...);
		}
	} else if (side==4) { # right
		if (r1>x1-x0) { r1=x1-x0; }
		if (2*r2>y1-y0) { r2=(y1-y0)/2; }
		lines(c(x0,x1-r1), c(y0,y0), ...);
		lines(c(x0,x1-r1), c(y1,y1), ...);
		lines(c(x1,x1), c(y0+r2,y1-r2), ...);
		if (r1>0 && r2>0) {
			cx = x1-r1;
			cy = y0+r2;
			curve( cy-r2*tmp_f(x,cx,r1), xlim=c(x1-r1,x1), add=T, ...);
			cy = y1-r2;
			curve( cy+r2*tmp_f(x,cx,r1), xlim=c(x1-r1,x1), add=T, ...);
		}
	}
}

# hclust/dendrogram help function
# {{{
# only available for continue leaves
sub_hclust <- function(h, idx)
{
	out = h;
	
	n = nrow(h$merge)+1;
	out$order = 1:length(idx);
	tmp_rank = rep(0,n);
	tmp_rank[h$order] = (1:n);
	j0 = min(tmp_rank[idx]);
	j1 = max(tmp_rank[idx]);
	is = (1:(n-1))[ (h$merge[,1] %in% -idx) | (h$merge[,2] %in% -idx) ];
	is0 = c();
	while(length(is)>1) {
		is0 = c(is0,is);
		is = (1:(n-1))[ ((h$merge[,1] %in% is) | (h$merge[,2] %in% is))];
		is = is[!(is%in%is0)]
	}
	rev_is = rep(0,n-1);
	rev_is[is0] = 1:length(is0);
	out$merge = h$merge[is0,];
	out$merge[out$merge[,1]>0,1] = rev_is[out$merge[out$merge[,1]>0,1]];
	out$merge[out$merge[,2]>0,2] = rev_is[out$merge[out$merge[,2]>0,2]];
	out$merge[out$merge[,1]<0,1] = -(tmp_rank[-out$merge[out$merge[,1]<0,1]]-j0+1);
	out$merge[out$merge[,2]<0,2] = -(tmp_rank[-out$merge[out$merge[,2]<0,2]]-j0+1);
	out$height = h$height[is0];
	out$labels = h$labels[h$order[idx]];
	out$old_to_new = tmp_rank-j0+1;

	out
}

sub_hclust2 <- function(h, idx)
{
	h1=h;
	n=length(h$order);
	n1 = length(idx);
	to_new = rep(0,2*n-1);
	to_new[idx] = 1:n1;
	h1$order = to_new[h$order];
	h1$order = h1$order[h1$order>0];
	h1$labels = h$labels[idx];

	h1$merge = NULL;
	h1$height = rep(0, n1-1);
	m1 = n1;
	back = to_new;
	for (i in 1:(n-1)) {
		j1=h$merge[i,1];
		j2=h$merge[i,2];
		if (j1<0) { j1= -j1; }
		else if (j1>0) { j1= n+j1; }
		if (j2<0) { j2= -j2; }
		else if (j2>0) { j2= n+j2; }
		j3 = back[j1];
		j4 = back[j2];
		if (j3>0) {
			if (j4>0) {
				m1 = m1+1;
				to_new[n+i] = m1;
				back[n+i] = m1;
				if (j3<=n1) { j5=-j3; } else { j5=j3-n1; }
				if (j4<=n1) { j6=-j4; } else { j6=j4-n1; }
				h1$merge = rbind(h1$merge, c(j5,j6));
				h1$height[m1-n1] = h$height[i];
			} else {
				back[n+i] = j3;
			}
		} else {
			if (j4>0) {
				back[n+i] = j4;
			} else {
			}
		}
	}
	
	h1
}

hclust_more <- function(data, h, nc)
{
	res = list(h=h, nc=nc);
	res$cid = cutree(h, k=nc);
	n = res$n = length(res$cid);
	res$csz = tapply( rep(1,length(res$cid)), res$cid, sum);
	res$order = h$order;
	res$rank=h$order;    res$rank[h$order] = 1:n;
	res$dnd = as.dendrogram(h);

	res$cid = c(res$cid, rep(0,n-1));
	res$cnode = rep(0, nc);
	res$cdist = matrix(0, nc, nc);
	nodes=list();
	for (i in 1:(n-1)) {
		j1 = h$merge[i,1];
		j2 = h$merge[i,2];
		if (j1<0) { j1= -j1; }
		else if (j1>0) { j1= n+j1; }
		if (j2<0) { j2= -j2; }
		else if (j2>0) { j2= n+j2; }
		cid1 = res$cid[j1];
		cid2 = res$cid[j2];
		if (cid1==cid2 && cid1!=0) { res$cid[n+i] = cid1; res$cnode[cid1]=n+i; }
		else {
			if (cid1>0) {
				if (cid2>0) {
					res$cdist[cid1,cid2] = h$height[i];
					res$cdist[cid2,cid1] = h$height[i];
					nodes[[i]] = c(cid1,cid2);
				} else {
					cid2s = nodes[[j2-n]];
					res$cdist[cid1,cid2s] = 2*h$height[i] - h$height[res$cnode[cid1]-n] - h$height[res$cnode[cid2s]-n]; 
					res$cdist[cid2s,cid1] = res$cdist[cid1,cid2s];
					nodes[[i]] = c(cid1, cid2s);
				}
			} else {
				if (cid2>0) {
					cid1s = nodes[[j1-n]];
					res$cdist[cid1s,cid2] = 2*h$height[i] - h$height[res$cnode[cid1s]-n] - h$height[res$cnode[cid2]-n]; 
					res$cdist[cid2,cid1s] = res$cdist[cid1s,cid2];
					nodes[[i]] = c(cid1s, cid2);
				} else {
					cid1s = nodes[[j1-n]];
					cid2s = nodes[[j2-n]];
					pair = cbind(rep(cid1s,each=length(cid2s)), rep(cid2s,length(cid1s)));
					res$cdist[pair] = 2*h$height[i] - h$height[res$cnode[pair[,1]]-n] - h$height[res$cnode[pair[,2]]-n]; 
					res$cdist[pair[,2:1]] = res$cdist[pair];
					nodes[[i]] = c(cid1s, cid2s);
				}
			}
		}
	}

	res$parent = rep(0, n*2-1);
	res$height = rep(0, n*2-1);
	for (k in 1:2) {
		f=h$merge[,k]<0;
		res$parent[ -h$merge[f,k] ] = n+(1:(n-1))[f];
		f=h$merge[,k]>0;
		res$parent[ n+h$merge[f,k] ] = n+(1:(n-1))[f];
	}
	res$height[(n+1):(2*n-1)] = h$height;

	res$c_mean = matrix(0, nc, ncol(data));
	for (i in 1:ncol(data)) {
		res$c_mean[,i] = tapply(data[,i], res$cid[1:res$n], mean);
	}

	res;
}
# }}}

# draw chromosome region
plot_chr_lfpkm <- function(genes, lfpkm, nts, color)
{
#	f = genes[,'chr']=='chr' & genes[,'begin']>=begin & genes[,'end']<=end;
#	genes1 = genes[f,];
	begin0 = min(genes[,'begin']);
	end0   = max(genes[,'end']);
	ng = nrow(genes1);
	par(mar=c(1,1,1,1)+0.1)
	plot(c(0,100*ng), c(0,100), col='white', xpd=T, ann=F, axes=F);
	gene_y0 = 80;
	gene_thick = 10;
	gene_y1 = gene_y0+gene_thick;
	exp_y0 = 20;
	exp_thick = 10;
	for (i in 1:ng) {
		gid = genes[i,'id'];
		s   = genes[i,'strand'];
		beg = genes[i,'begin'];
		end = genes[i,'end'];
		x0 = 100*i-90;
		x1 = 100*i-10;
		rect(x0, gene_y0, x1, gene_y1, col=color[i]);
		m = ncol(lfpkm)/nts;
		if (!(gid %in% rownames(lfpkm))) { next; }
		z = NULL;
		for (j in 1:m) {
			z = rbind(z,lfpkm[gid,(nts*(j-1)+1):(j*nts)]);
		}
		z = z[m:1,];
		image(seq(x0,x1,length.out=nts+1), exp_y0+(0:m)*exp_thick, t(z), add=T, col=rgb(1,(256:0)/256,(256:0)/256));
	}
}























