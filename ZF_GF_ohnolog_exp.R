library('gplots'); # hist2d, heatmap.2
library('tsne'); # tsne
library('entropy');
library('DESeq2');
library('edgeR');
library('limma');
library('fastcluster');
library('GOstats');
library("extrafont");
library("plot3D"); # scatter3D
library("rgl");
library("car"); # scatter3d
library(RColorBrewer)
font_import();
loadfonts();
par(family="Arial");

source('/Users/chenz11/my_program3/src/goldfish_project/ZF_GF_CC_ohnolog_exp.func.R');

options(stringsAsFactors=FALSE);
wd="~/data/goldfish/11549472/sergey_canu70x/arrow/carAur03/big/WGD2/";
setwd(wd);

tpm_thres=1;
tissues=c("Brain", "Eye", "Heart", "Gill", "Muscle", "Tail");
tissues1=c('B', 'E', 'H', 'G', 'M', 'T');
nts = length(tissues);
sps=c('ZF', 'GF'); # DONOT change the order
sps2=c('ZF', 'GF', 'GF'); # DONOT change the order
sps3=c('ZF', 'GFa', 'GFb'); # DONOT change the order
copies=c(1,2,2);
names(copies)=sps;
tot_copy=3;
tot_sp = length(sps);
fate_names=c('double_high_corr', 'double_high_conserved', 'high_dosage_balance', 'mid_dosage_balance','cor_subfunc','cor_neofunc','nonfunc', 'on_off_subfunc', 'on_off_neofunc');
fate_names1=c('high double correlated', 'double conserved', 'high dosage balance', 'dosage balance','sub-func','neo-func','non-func', 'sub-functionalized', 'neo-functionalized');
names(fate_names1) = fate_names;

on_off_fate_names = c('double_coexpressed', 'on_off_subfunc', 'on_off_neofunc_WGD', 'nonfunc_WGD', 'on_off_part_subfunc', 'on_off_part_neofunc_WGD', 'part_nonfunc_WGD');
on_off_disp_fate_names = c('co-express', 'sub-F', 'neo-F', 'non-F', 'partial sub-F', 'partial neo-F', 'partial non-F')
names(on_off_disp_fate_names) = on_off_fate_names;

# red, green, blue, cyan, magenta, yellow, gray
color8 <- c("#FF0000", "#00FF00", "#0000FF", "#00FFFF", "#FF00FF", "#FFFF00", "#808080", "#FFC0CB");
color8a <- paste(color8, "60", sep="");
color <- c("#FF0000", "#00FF00", "#0000FF", "#00FFFF", "#FF00FF", 
		"#FFFF00", "#808080", "#FFC0CB", "#D00000", "#00D000", 
		"#000080", "#008080", "#800080", "#808000", "#7FFFE0", 
		"#FFE0D0", "#FFA010", "#A03030", "#DEB888", "#609AA0",
		"#80FF00", "#9E0142", "#D53E4F", "#F46D43", "#FDAE61",
	   	"#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5",
	   	"#3288BD", "#5E4FA2");
blueyellow255=colorpanel(255, low='blue',mid='white',high='yellow');
yellowblue255=colorpanel(255, low='yellow',mid='white',high='blue');
# color <- rainbow(32);


##########################
# Load GO, Interpro  annotation
##########################
# {{{
if (file.exists('Rdata/func_annot.Rdata')) {
	load('Rdata/func_annot.Rdata');
} else {
	id_to_gs0<-list();
	id_to_gs<-list();
	gs <- list();
	for (f in c("GO2_MF", "GO2_CC", "GO2_BP")) {
		id_to_gs0[[f]] = list();
		id_to_gs[[f]] = list();
		gs[[f]]=list();
		for (sp in c("ZF","CC","GF")) {
			fname = paste0('input/',sp,'/gene.ips.id_to_', f, ".txt");
			if (file.exists(fname)) {
				id_to_gs0[[f]][[sp]] <- read.table(fname, head=F, col.names=c("gene_id", "annot_id", "annot_name", "level"), sep="\t", quote="");
				id_to_gs[[f]][[sp]] <- split(id_to_gs0[[f]][[sp]][,2], id_to_gs0[[f]][[sp]][,1]);
				gs[[f]][[sp]] <- split(id_to_gs0[[f]][[sp]][,1], id_to_gs0[[f]][[sp]][,2]);
			}
		}
	}
	for (f in c("Reactome", "KEGG", "Interpro") ) {
	for (sp in c("ZF","CC","GF")) {
		fname = paste0('input/',sp,'/gene.ips.id_to_', f, ".txt");
		if (file.exists(fname)) {
			id_to_gs0[[f]][[sp]] <- read.table(fname, head=F, col.names=c("gene_id", "annot_id", "annot_name"), sep="\t", quote="");
			id_to_gs0[[f]][[sp]] <- cbind(id_to_gs0[[f]][[sp]], rep(0, dim(id_to_gs0[[f]][[sp]])[1]));
			colnames(id_to_gs0[[f]][[sp]]) = c("gene_id", "annot_id", "annot_name", "level");
			id_to_gs0[[f]][[sp]][,4]=0;
			id_to_gs[[f]][[sp]] <- split(id_to_gs0[[f]][[sp]][,2], id_to_gs0[[f]][[sp]][,1]);
			gs[[f]][[sp]] <- split(id_to_gs0[[f]][[sp]][,1], id_to_gs0[[f]][[sp]][,2]);
		}
	}
	}
	gs_id_name <- list();
	for (f in c("GO2_MF", "GO2_CC", "GO2_BP", "Reactome", "KEGG", "Interpro")) {
	for (sp in c("ZF","CC","GF")) {
		if (!is.null(id_to_gs0[[f]][[sp]])) {
			gs_id_name[[f]] <- rbind(gs_id_name[[f]], unique(id_to_gs0[[f]][[sp]][,2:4]));
		}
	}
		gs_id_name[[f]] = unique(gs_id_name[[f]]);
		rownames(gs_id_name[[f]]) <- gs_id_name[[f]][,1];
	}

	child_to_parent = list();
	fname = paste0('input/Interpro.66.0.child_to_parent.txt');
	child_to_parent$Interpro = read.table(fname, head=F, sep="\t");
	rownames(child_to_parent$Interpro) = child_to_parent$Interpro[,1];
	save(id_to_gs0, id_to_gs, gs, gs_id_name, child_to_parent, file='Rdata/func_annot.Rdata');
}
## }}}
########################

###########################
# Load gene annotation
###########################
# {{{
#bgp <- list(
#		GF=>read.table("~/data/goldfish/11549472/sergey_canu70x/arrow/carAur03/big/carAur03.noM.gene.bgp", stringsAsFactor=F, head=F, sep="\t"),
#		CC=>read.table("~/data/common_carp/NCBI_Cyprinus_carpio/ref_common_carp_genome_top_level.rename_chr.re_id.gff3.bgp"),
#		ZF=>read.table("~/data/ensembl85/bgp/Danio_rerio.GRCz10.85.bgp"));

if (file.exists("Rdata/annot.Rdata")) {
	load("Rdata/annot.Rdata");
} else {
	bgp <- list(
			GF=read.table(paste0(wd,'/input/GF.gene.bgp'), stringsAsFactor=F, head=F, sep="\t"),
			CC=read.table(paste0(wd,'/input/CC.gene.bgp'), stringsAsFactor=F, head=F, sep="\t"),
			ZF=read.table(paste0(wd,'/input/ZF.gene.bgp'), stringsAsFactor=F, head=F, sep="\t") );

	annot = list();
	genes = list();
	for (sp in c('GF','CC','ZF')) {
		rownames(bgp[[sp]]) = bgp[[sp]][,4]; 
		a1 = list();
		a1$gene_names = bgp[[sp]][,18];
		a1$exon_len=strsplit(bgp[[sp]][,11], ",", fixed=T);
		a1$exon_len <- lapply(a1$exon_len, as.numeric); names(a1$exon_len) = bgp[[sp]][,4];
		a1$exon_num <- bgp[[sp]][,10]; names(a1$exon_num) = bgp[[sp]][,4];
		a1$tran_len <- sapply(a1$exon_len, sum); names(a1$tran_len) = bgp[[sp]][,4];
		annot[[sp]] = a1;
		f=bgp[[sp]][,1]!='chrMT';
		a2 = unique(bgp[[sp]][f,c(1,6,18,19)]);
		colnames(a2) = c('chr', 'strand', 'id', 'name');
		rownames(a2) = a2[,'id'];
		gids = bgp[[sp]][f,18];
		begs = tapply(bgp[[sp]][f,2], gids, min)[rownames(a2)];
		ends = tapply(bgp[[sp]][f,3], gids, max)[rownames(a2)];
		genes[[sp]] = cbind(a2, begin=begs, end=ends);
	}
	save(bgp, annot, genes, file="Rdata/annot.Rdata")
}
# }}}
#############################


###########################
# Load clusters
###########################
## {{{
#dfile=paste0(wd,"/Rdata/anchor.Rdata");
#if (file.exists(dfile)) {
#	load("Rdata/anchor.Rdata");
#} else {
#	file=paste0(wd, "/input/anchor.fix");
#	anchor=read.table(file, head=T, sep='\t');
##	tmp.clust = read.table("input/rescue_m1.6.cluster.txt", sep='\t', head=F);
##	f.ZF = grep('\\(ZF', tmp.clust[,4]);
##	f.CC = grep('\\(CC', tmp.clust[,4]);
##	f.GF = grep('\\(GF', tmp.clust[,4]);
#	
#	anchor[,'anchor_tid'] = gsub("\\.[0-9]+", "", anchor[,'anchor_tid'] );
#	anchor[,'anchor_gid'] = gsub("\\.[0-9]+", "", anchor[,'anchor_gid'] );
##	anchor = anchor[anchor[,"clust_id"]>=0 | anchor[,'to_dup']!=".",];
#	rownames(anchor) = paste0(anchor[,'species'],'|',anchor[,'anchor_gid']); # name is gene id
#	anchor[,'clust_id'] = as.numeric(anchor[,'clust_id']);
#	ncluster = max(anchor[,'clust_id'])+1;
##	cluster_is_dup = rep(F,ncluster); 
##	cluster_is_dup[f.ZF] = T;
##	cluster_is_dup[f.CC] = T;
##	cluster_is_dup[f.GF] = T;
#	cluster = list();
#	length(cluster)=ncluster;
#	cluster_sizes=matrix(0, ncluster, 5);
#	colnames(cluster_sizes) = c('total', 'ZF', 'GC', 'CC', 'GF');
#	for (i in 1:nrow(anchor)) {
#		cid=anchor[i,'clust_id']+1;
#		if (cid==0) { next; }
#		sp = anchor[i,'species'];
#		if (is.null(cluster[[cid]])) cluster[[cid]]=list();
#		if (anchor[i,'to_dup']==".") {
#			cluster[[cid]][['total']] = c(cluster[[cid]][['total']], i);
#			cluster[[cid]][[sp]] = c(cluster[[cid]][[sp]], i);
#			cluster_sizes[cid, sp]=cluster_sizes[cid, sp]+1;
#			cluster_sizes[cid, 'total']=cluster_sizes[cid, 'total']+1;
#		} 
#	}
#	anchor_dup = list();
#	length(anchor_dup)=nrow(anchor);
#	names(anchor_dup) = rownames(anchor);
#	f=anchor[,'to_dup']!=".";
#	dup = paste0(anchor[,'species'],'|',anchor[,'to_dup']);
#	for (i in (1:nrow(anchor))[f]) {
#		anchor_dup[[dup[i]]] = c(anchor_dup[[dup[i]]], i); 
#	}
#	for (i in 1:nrow(anchor)) {
#		if (length(anchor_dup[[i]])>0) {
#			names(anchor_dup[[i]]) = rownames(anchor)[ anchor_dup[[i]] ];
#		}
#	}
##	cluster = split(1:nrow(anchor), anchor[,'clust_id']);
#
#	file=paste0(wd, "/input/qC50.loci.txt");
#	anchor2=read.table(file, head=T, sep='\t');
#	anchor2[,'anchor_gid'] = gsub("\\.[0-9]+", "", anchor2[,'anchor_gid'] );
#	anchor2 = anchor2[anchor2[,"clust_id"]>=0,];
#	rownames(anchor2) = paste0(anchor2[,'species'],'|',anchor2[,'anchor_gid']); # name is gene id
#	anchor2[,'clust_id'] = as.numeric(anchor2[,'clust_id']);
#	ncluster2 = max(anchor2[,'clust_id'])+1;
#	cluster2 = list();
#	length(cluster2)=ncluster2;
#	cluster_sizes2=matrix(0, ncluster2, 5);
#	colnames(cluster_sizes2) = c('total', 'ZF', 'GC', 'CC', 'GF');
#	for (i in 1:nrow(anchor2)) {
#		cid=anchor2[i,'clust_id']+1;
#		if (cid==0) { next; }
#		sp = anchor2[i,'species'];
#		if (is.null(cluster2[[cid]])) cluster2[[cid]]=list();
#		cluster2[[cid]][[sp]] = c(cluster2[[cid]][[sp]], i);
#		cluster_sizes2[cid, sp]=cluster_sizes2[cid, sp]+1;
#		cluster_sizes2[cid, 'total']=cluster_sizes2[cid, 'total']+1;
#	}
#
#	cluster_is_changed = rep(F,ncluster); # 0: no change, 1: changed by anchor2
#	for (cid in 1:ncluster) {
#		if (is.null(cluster2[[cid]])) { cluster_is_changed[cid]=T; next; }
#		if (cid>ncluster2) { cluster_is_changed[cid]=T; next; }
#		m=0;
#		for (j in 1:tot_sp) {
#			if (cluster_sizes[cid,sps[j]]==cluster_sizes2[cid,sps[j]]) {m=m+1; }
#			else {break;}
#		}
#		if (m<tot_sp) cluster_is_changed[cid]=T;
#	}
#	rm(anchor2, cluster2, ncluster2, cluster_sizes2);
#
#	edge = list(prot=read.table('input/blastp.ZF_GC_CC_GF.gene.f3.m6.gz', head=F),
#			nucl= read.table('input/blastn.ZF_GC_CC_GF.gene.f3.m6.gz', head=F));
#	rownames(edge$prot) = paste0(edge$prot[,1], "-", edge$prot[,2]);
#	rownames(edge$nucl) = paste0(edge$nucl[,1], "-", edge$nucl[,2]);
#	cluster_is_dup=rep(F, ncluster);
#	cluster_is_dup3=matrix(rep(F, ncluster*tot_sp), ncluster, tot_sp);
#	colnames(cluster_is_dup3) = sps;
#	for (i in 1:ncluster) {
#		if (cluster_sizes[i,'ZF']==0) { next; }
#		zf_idx = cluster[[i]][['ZF']][1];
#		if (length(anchor_dup[[zf_idx]])>0) {
#			cluster_is_dup[i]=T;
#			cluster_is_dup3[i,'ZF']=T;
#			next;
#		}
#		for (sp in c('CC','GF')) {
#			bad=F;
#			idxs = cluster[[i]][[sp]];
#			for (idx in idxs) {
#				if (length(anchor_dup[[idx]])>0) {
#					id_pair = paste0(rownames(anchor)[zf_idx], '-', rownames(anchor)[idx]);
#					cov0 = edge$nucl[id_pair,18];
#					if (is.na(cov0) || is.null(cov0)) { bad=T; }
#					else {
#						id_pairs = paste0(rownames(anchor)[zf_idx], '-', rownames(anchor)[anchor_dup[[idx]]]);
#						cov = max(edge$nucl[id_pairs,18]);
#						if (!is.na(cov) && cov>25 && cov>=cov0/2) { bad=T; }
#					}
#					if (bad) { break;}
#				}
#			}
#			cluster_is_dup3[i,sp]=bad;
#		}
#		if (cluster_is_dup3[i,'CC'] || cluster_is_dup3[i,'GF'] ) { cluster_is_dup[i]=T; }
#	}
#
#	cluster_is_syn = rep(F,ncluster);
#	tmp_idxs=unlist(sapply(1:ncluster,function(x) {cluster[[x]][['total']]}));
#	tmp_cids=unlist(sapply(1:ncluster,function(x) {rep(x,cluster_sizes[x,'total'])}));
#	cluster_is_syn[unique(tmp_cids)] = tapply(tmp_idxs, tmp_cids, function(x) {sum(anchor[x,'to_block_each_species']!="")>=3}, default=F);
#	save(anchor,anchor_dup, ncluster,cluster,cluster_sizes, cluster_is_changed, cluster_is_syn, cluster_is_dup, edge, file="Rdata/anchor.Rdata");
#}
## }}}
###########################

###########################
# Load ZF-GF-GF Gene triplets
###########################
if ( file.exists("Rdata/ZF_GF.data.Rdata") ) {
	load("Rdata/ZF_GF_data.Rdata")
} else {
	data1 = list();
	a = read.table('input/dcs_lost_run3/long20000.lost.out.gene.triplet.stat.txt', head=T, comment.char='', check.names=F);
	colnames(a)[1] = gsub('^#', '', colnames(a)[1])
	rownames(a) = a[,1];
	data1$gene_ids2 = a[,1:3];
	colnames(a) = gsub('^.*([01][01][01])_.*$', '\\1', colnames(a));
	colnames(data1$gene_ids2) = c('ZF','GFa','GFb');
	data1$exon_loss_num  = a[,3+1:(8*1)];
	data1$exon_loss_len  = a[,3+(8*1+1):(8*2)];
	mm=apply(data1$exon_loss_len,1,sum);
	data1$exon_loss_perc = data1$exon_loss_len/mm*100;
	data1$exon_loss_perc[mm==0,1] = 100;
	data1$exon_loss_perc[mm==0,2:8] = 0;

	data1$CNE_loss_num  = a[,3+(8*2+1):(8*3)];
	data1$CNE_loss_len  = a[,3+(8*3+1):(8*4)];
	mm=apply(data1$CNE_loss_len,1,sum);
	data1$CNE_loss_perc = data1$CNE_loss_len/mm*100;
	data1$CNE_loss_perc[mm==0,1] = 100;
	data1$CNE_loss_perc[mm==0,2:8] = 0;

	data1$CNE5k_loss_num  = a[,3+(8*4+1):(8*5)];
	data1$CNE5k_loss_len  = a[,3+(8*5+1):(8*6)];
	mm=apply(data1$CNE5k_loss_len,1,sum);
	data1$CNE5k_loss_perc = data1$CNE5k_loss_len/mm*100;
	data1$CNE5k_loss_perc[mm==0,1] = 100;
	data1$CNE5k_loss_perc[mm==0,2:8] = 0;

	data1$CNE10k_loss_num  = a[,3+(8*6+1):(8*7)];
	data1$CNE10k_loss_len  = a[,3+(8*7+1):(8*8)];
	mm=apply(data1$CNE10k_loss_len,1,sum);
	data1$CNE10k_loss_perc = data1$CNE10k_loss_len/mm*100;
	data1$CNE10k_loss_perc[mm==0,] = c(100,rep(0,7));
	data1$CNE10k_loss_perc[mm==0,1] = 100;
	data1$CNE10k_loss_perc[mm==0,2:8] = 0;

	a1 = read.table('input/dcs_lost_run3/long20000.lost.out.gene.triplet.blastn.align_stat.txt', head=F, comment='');
	data1$nucl_align=list();
	names=c('size1', 'size2', 'cov1', 'cov2', 'iden', 'alen', 'bit', 'score');
	for (i in 1:8) {
		name = names[i];
		data1$nucl_align[[name]] = a1[,c(3+(i+(1:9-1)*8))];
		rownames(data1$nucl_align[[name]]) = rownames(data1$gene_ids2);
	}

	a2 = read.table('input/dcs_lost_run3/long20000.lost.out.gene.triplet.blastp.align_stat.txt', head=F, comment='');
	data1$prot_align=list();
	names=c('size1', 'size2', 'cov1', 'cov2', 'iden', 'alen', 'bit', 'score');
	for (i in 1:8) {
		name = names[i];
		data1$prot_align[[name]] = a2[,c(3+(i+(1:9-1)*8))];
		rownames(data1$prot_align[[name]]) = rownames(data1$gene_ids2);
	}

	data1$gene_names2 = NULL;
	for (i in 1:tot_copy) {
		sp = sps2[i];
		a=unique(bgp[[sp]][,18:19]);
		rownames(a) = a[,1];
		data1$gene_names2 = cbind(data1$gene_names2, a[ data1$gene_ids2[,i], 2]);
	}

###########################
# Load TPM
# {{{
	exp = list();
	exp[['tpm']] = list(thres=1);
	exp[['fpkm']] = list(thres=1);
	exp[['count']] = list(thres=5);
	for (type in c('tpm', 'fpkm', 'count')) {
		tpm_file=list(
				ZF=paste0(wd,'/input/ZF/gene.', type, '.mat'),
				CC=paste0(wd,'/input/CC/gene.', type, '.mat'),
				GF=paste0(wd,'/input/GF/gene.', type, '.mat')
				);
		tpm_sp=c();
		tpm=c();
		for (sp in sps) {
			a=read.table(tpm_file[[sp]], sep="\t", head=T, row.names=1);
			a=a[,tissues];
			tpm=rbind(tpm,a);
			tpm_sp = c(tpm_sp, rep(sp,nrow(a)));
		}
		rownames(tpm) = remove_id_version(rownames(tpm));
		rownames(tpm) = paste0(tpm_sp,'|',rownames(tpm));

		exp[[type]][['raw']] = tpm;
		exp[[type]][['raw_sp']] = tpm_sp;
	}

	for (type in c('tpm', 'fpkm', 'count')) {
		sp_gids0 = rownames(exp[[type]][['raw']]);
		tpm2 = matrix(0, nrow(data1$gene_ids2), tot_copy*nts);
		rownames(tpm2) = rownames(data1$gene_ids2);
		k0 = 0;
		j=0;
		for (sp in sps) {
			for (i in 1:copies[sp]) {
				j=j+1;
				gids0 = data1$gene_ids2[,1];
				gids = data1$gene_ids2[,j];
				sp_gids = paste0(sp, '|', gids);
				f =  sp_gids %in% sp_gids0 ;
				
				tpm2[gids0[f],(k0+1):(k0+nts)] = as.matrix(exp[[type]][['raw']][sp_gids[f],]);
				k0 = k0+nts;
			}
		}
		colnames(tpm2)=c( paste0(rep('ZF',nts), '.', tissues), 
			paste0(rep('GFa',nts), '.', tissues), paste0(rep('GFb',nts), '.', tissues) );
		sum_tpm2 = ohno_sum(tpm2, tissues, sps2);
		data1[[type]] = cbind(tpm2,sum_tpm2);
		ltype = paste0('log2_',type);
		data1[[ltype]] = log2(data1[[type]]+1);
# sum of two paralog

#		qnorm_ltpm2 = normalizeBetweenArrays(ltpm2[,1:(nts*tot_copy)]);
		low = apply(data1[[type]][,1:(tot_copy*nts)]>=thres, 1, sum)<3;
		data1[[ paste0('low_',type) ]] = low;
		data1[[ paste0('low2_',type) ]] = NULL;
		for (i in 1:(tot_copy+tot_sp)) {
			low = apply(data1[[type]][,(nts*(i-1)+1):(nts*i)]>=thres, 1, sum)==0;
			data1[[ paste0('low2_',type) ]] = cbind(data1[[ paste0('low2_',type) ]], low);
		}
	}
# }}}

# set interpros for genes
	data1$gene_Interpros2 = list();
	for (i in 1:tot_copy) {
		sp = sps2[i];
		data1$gene_Interpros2[[i]] = id_to_gs$Interpro[[sp]][ data1$gene_ids2[,i] ];
		names(data1$gene_Interpros2[[i]]) = data1$gene_ids2[,i] ;
	}
##

	thres=1;

##############################################
# compute score between ortholog-ohnolog or ohnolog-ohnolog
##############################################
# {{{
	data1$exp = list();
	for (etype in c('fpkm', 'tpm')) {
		data1$exp[[etype]] = list();
	}

	for (etype in c('fpkm', 'tpm')) {
		fpkm  = data1[[etype]];
		lfpkm = data1$exp[[etype]]$log2 = log2(data1[[etype]]+1);

		m=nrow(lfpkm);
		n=ncol(lfpkm)/nts;
		data1$exp[[etype]]$median2 = matrix(0, m, tot_copy+tot_sp);
		data1$exp[[etype]]$mean2   = matrix(0, m, tot_copy+tot_sp);
		data1$exp[[etype]]$sd2     = matrix(0, m, tot_copy+tot_sp);
		data1$exp[[etype]]$max2    = matrix(0, m, tot_copy+tot_sp);
		for (i in 1:n) {
			a = lfpkm[, (nts*(i-1)+1):(nts*i)];
			data1$exp[[etype]]$median2[,i] = apply(a, 1, median);
			data1$exp[[etype]]$mean2[,i]   = apply(a, 1, mean);
			data1$exp[[etype]]$sd2[,i]     = apply(a, 1, sd);
			data1$exp[[etype]]$max2[,i]    = apply(a, 1, max);
		}
		data1$exp[[etype]]$norm2 = lfpkm;
		data1$exp[[etype]]$zscore2 = lfpkm;
		for (i in 1:n) {
			data1$exp[[etype]]$norm2[,(nts*(i-1)+1):(nts*i)] = lfpkm[,(nts*(i-1)+1):(nts*i)] - data1$exp[[etype]]$mean2[,i];
			sd1 = data1$exp[[etype]]$sd2[,i];
			sd1[is.na(sd1)] = 1;
			data1$exp[[etype]]$zscore2[,(nts*(i-1)+1):(nts*i)] = data1$exp[[etype]]$norm2[,(nts*(i-1)+1):(nts*i)] / sd1;
		}
		data1$exp[[etype]]$exp[[etype]]$median1 = apply(lfpkm[,1:(nts*tot_copy)], 1, median);
		data1$exp[[etype]]$mean1   = apply(lfpkm[,1:(nts*tot_copy)], 1, mean);
		data1$exp[[etype]]$sd1     = apply(lfpkm[,1:(nts*tot_copy)], 1, sd);
		data1$exp[[etype]]$zscore1 = (lfpkm-data1$exp[[etype]]$mean1)/data1$exp[[etype]]$sd1;

# count number of expressed samples
		x <- matrix(0, nrow(lfpkm), tot_copy+tot_sp);
		rownames(x) = rownames(lfpkm);
		colnames(x) = c(sps3, paste0("sum_of_", sps));
		for (i in 1:ncol(x)) {
			a=lfpkm[, ((nts)*(i-1)+1):(nts*i)];
			x[,i] <- apply(a>=thres, 1, sum);
		}
		data1$exp[[etype]]$exp_sm_count = x;
# tissue-specific score (<=1, near 1 is single tissue-specific)
		x = matrix(0, nrow(lfpkm), tot_copy+tot_sp);
		y = matrix(0, nrow(lfpkm), tot_copy+tot_sp);
		rownames(x) = rownames(lfpkm);
		colnames(x) = c(sps3, paste0("sum_of_", sps));
		for (i in 1:ncol(x)) {
			a = lfpkm[,((i-1)*nts+1):(i*nts)];
			x[,i] = apply(a,1,max)/(apply(a,1,sum)+thres/10);
			f = x[,i]>=0.4;
			y[f,i] = max.col(a[f,]);

		}
		data1$exp[[etype]]$specific_score = x;
		data1$exp[[etype]]$specific_tissue = y;
# entropy (>=0, near 0 is more tissue-specific)
		x = matrix(0, nrow(lfpkm), tot_copy+tot_sp);
		rownames(x) = rownames(lfpkm);
		colnames(x) = c(sps3, paste0("sum_of_", sps));
		for (i in 1:ncol(x)) {
			a = lfpkm[,((i-1)*nts+1):(i*nts)];
			x[,i] = apply(a,1,function(x) {entropy(x+thres/100,unit='log2')});
		}
		data1$exp[[etype]]$entropy_score = x;
# gene pair correlation
		n=ncol(lfpkm)/nts;
		x=matrix(0, nrow(lfpkm), n*n);
		rownames(x) = rownames(lfpkm);
		y=matrix(0, nrow(lfpkm), n*n);
		rownames(y) = rownames(lfpkm);
		for (i in 1:n) {
			for (j in 1:n) {
				if (i==j) {next;}
				x[,(i-1)*n+j] = pair_cor(lfpkm[,(nts*(i-1)+1):(nts*i)], lfpkm[,(nts*(j-1)+1):(nts*j)]);
				f=is.na(x[,(i-1)*n+j]);
				x[f,(i-1)*n+j] = 0;
				y[,(i-1)*n+j] = pair_cor_p(lfpkm[,(nts*(i-1)+1):(nts*i)], lfpkm[,(nts*(j-1)+1):(nts*j)]);
			}
		}
		data1$exp[[etype]]$pair_cor = x;
		data1$exp[[etype]]$pair_cor_p = y;
# gene pair euclidean distance
		n=ncol(lfpkm)/nts;
		data1$exp[[etype]]$pair_euc_dist=matrix(0, nrow(lfpkm), n*n);
		for (i in 1:n) {
			for (j in 1:n) {
				if (i==j) {next;}
				data1$exp[[etype]]$pair_euc_dist[,(i-1)*n+j] = pair_euc_dist(lfpkm[,(nts*(i-1)+1):(nts*i)], lfpkm[,(nts*(j-1)+1):(nts*j)]);
			}
		}
# t.test, paired
		n=ncol(lfpkm)/nts;
		data1$exp[[etype]]$tp=matrix(0, m, n*n);
		for (i in 1:n) {
			for (j in 1:n) {
				if (i==j) {next;}
				data1$exp[[etype]]$tp[,(i-1)*n+j] = pair_t_p(lfpkm[,(nts*(i-1)+1):(nts*i)], lfpkm[,(nts*(j-1)+1):(nts*j)]);
				f=is.na(data1$exp[[etype]]$tp[,(i-1)*n+j]);
				data1$exp[[etype]]$tp[f,(i-1)*n+j] = 1.0;
			}
		}
# gene pair maximal distance
		data1$exp[[etype]]$pair_max_dist=matrix(0, m, n*n);
		for (i in 1:n) {
			for (j in 1:n) {
				if (i==j) {next;}
				data1$exp[[etype]]$pair_max_dist[,(i-1)*n+j] = apply(lfpkm[,c( (nts*(i-1)+1):(nts*i), (nts*(j-1)+1):(nts*j))], 1, function(x) {max(abs(x[1:nts]-x[(nts+1):(nts*2)]))});
			}
		}

# on-off
		x = matrix(0, nrow(fpkm), (tot_sp+tot_copy)^2*nts);
		rownames(x) = rownames(fpkm);
		for (k in c(2,3,5,8)) {
			i=floor((k-1)/(tot_sp+tot_copy))+1;
			j=(k-1)%%(tot_sp+tot_copy)+1;
			a1=fpkm[,(nts*(i-1)+1):(nts*i)];
			a2=fpkm[,(nts*(j-1)+1):(nts*j)];
			x[,(nts*(k-1)+1):(k*nts)] = a1<thres & a2>=thres & (a1==0 | a2/a1>=2);
			k2 = (j-1)*(tot_sp+tot_copy)+i;
			x[,(nts*(k2-1)+1):(k2*nts)] = a2<thres & a1>=thres & (a2==0 | a1/a2>=2);
		}
		data1$exp[[etype]]$on_off = x;

# number of co-expressed tissue
		m=nrow(lfpkm);
		n=ncol(lfpkm)/nts;
		x=matrix(0,m,n*n);
		y=matrix(0,m,n*n);
		z=matrix(0,m,n*n);
		for (i in 1:n) {
			for (j in 1:n) {
				if (i==j) {next;}
				a1=fpkm[,(nts*(i-1)+1):(nts*i)];
				a2=fpkm[,(nts*(j-1)+1):(nts*j)];
				a = a1>=thres & a2>=thres;
#			a = a1>=thres & (a2>=thres | a1-a2>thres*0.5) | a2>=thres & (a1>=thres | a2-a1>thres*0.5);
				x[,(i-1)*n+j] = apply(a, 1, sum);
				a = a1<thres & a2<thres;
				y[,(i-1)*n+j] = apply(a, 1, sum);
				a = (a1<0.1& a2>=thres & (a1==0 | a2/a1>=2)) | (a2<0.1& a1>=thres & (a2==0 | a1/a2>=2));
				z[,(i-1)*n+j] = apply(a, 1, sum);
			}
		}
		data1$exp[[etype]]$pair_coexp_tissue_num = x;
		data1$exp[[etype]]$pair_cosil_tissue_num = y;
		data1$exp[[etype]]$pair_on_off_tissue_num = z;
	}

# compute jacard index for Interpro domain shared between genes
# {{{
	m=nrow(data1$fpkm);
	n=tot_copy;
	x=matrix(0,m,n*n);
	y=matrix(0,m,n*n);
	rownames(x) = rownames(lfpkm);
	rownames(y) = rownames(lfpkm);
	for (i in 1:(n-1)) {
		for (j in (i+1):n) {
			a1 = data1$gene_Interpros2[[i]];
			n1 = sapply(a1, length);
			a2 = data1$gene_Interpros2[[j]];
			n2 = sapply(a2, length);
			x[,(i-1)*n+j] = sapply(1:m, function(x) {
				if(length(a1[[x]])==0 || length(a2[[x]])==0) {nn=0; }
				else {
					x1=child_to_parent$Interpro[a2[[x]],2];
					if (length(x1)>0) {x2=child_to_parent$Interpro[x1,2];} else {x2=c();}
					a=c(a2[[x]],x1,x2);
					nn1 = sum(a1[[x]] %in% a);
					x1=child_to_parent$Interpro[a1[[x]],2];
					if (length(x1)>0) {x2=child_to_parent$Interpro[x1,2];} else {x2=c();}
					a=c(a1[[x]],x1,x2);
					nn2 = sum(a2[[x]] %in% a);
					nn = sum(a1[[x]] %in% a2[[x]]);
					nn=nn1+nn2-nn;
				}
				nn } );
			x[,(j-1)*n+i] = x[,(i-1)*n+j];
			y[,(i-1)*n+j] = x[,(i-1)*n+j]/(n1+n2-x[,(i-1)*n+j]);
			y[,(j-1)*n+i] = y[,(i-1)*n+j];
		}
	}
	data1$pair_interpor_shared_num = x;
	data1$pair_interpor_jaccard = y;
# }}}

# compute exon/CNE difference length, percentage
# {{{
	for (feat in c('exon', 'CNE5k')) {
		for (type in c('num', 'len', 'perc')) {
			name = paste0(feat, '_loss_', type);
			name2 = paste0(feat, '_', type, '_diff');

			m=nrow(data1$fpkm);
			n=tot_copy;
			x=matrix(0,m,n*n);
			x[,2] = data1[[name]]$`010` + data1[[name]]$`100` + data1[[name]]$`011` + data1[[name]]$`101` ; 
			x[,3] = data1[[name]]$`001` + data1[[name]]$`100` + data1[[name]]$`011` + data1[[name]]$`110` ; 
			x[,6] = data1[[name]]$`010` + data1[[name]]$`001` + data1[[name]]$`110` + data1[[name]]$`101` ; 
			rownames(x) = data1$gene_ids2[,1];
			data1[[name2]] = x;
		}
	}
# }}}

# set gene in each chromsome
	data1$chr_gids = list();
	tmp$gids = data1$gene_ids2[,'ZF'];
	tmp$a = unique(bgp$ZF[,c(1,18)]);
	f = regexpr('^[0-9]+', tmp$a[,1])>0 & (tmp$a[,2] %in% tmp$gids);
	data1$chr_gids$ZF = split(tmp$a[f,2], as.factor(tmp$a[f,1]));

	tmp$gids = unlist(data1$gene_ids2[,2:3]);
	tmp$a = unique(bgp$GF[,c(1,18)]);
	f = regexpr('^LG[0-9]', tmp$a[,1])>0 & (tmp$a[,2] %in% tmp$gids);
	data1$chr_gids$GF = split(tmp$a[f,2], as.factor(tmp$a[f,1]));

#}}}

	data1$f = !data1$low_fpkm & !data1$low_tpm & data1$nucl_align$iden[,2]<98 & data1$nucl_align$iden[,3]<98 & data1$nucl_align$iden[,6]<98 


	save(data1, file='Rdata/ZF_GF_data.Rdata');
}

thres = 1;
f0 = data1$f & data1$exp$fpkm$max2[,1]>=thres & data1$exp$tpm$max2[,1]>=thres
#############################################
# stat of exon loss
#############################################
# {{{
tmp=list();
tmp$GF1_GF2_reciprocal_triplet_num  = sum( apply(data1$exon_loss_num[,c('001','101')],1,sum)>0 & apply(data1$exon_loss_num[,c('010','110')],1,sum)>0 );
tmp$ZF_gain_GF1_GF2_reciprocal_triplet_num  = sum( data1$exon_loss_num[,'101']>0 & data1$exon_loss_num[,'110']>0 );

tmp$x = as.matrix(data1$exon_loss_num);
tmp$tot_num  = sum(apply(tmp$x[f0,], 1, sum));
tmp$config_num  = apply(tmp$x[f0,], 2, sum);
tmp$ZF_gain_GF1_loss_num    = sum(tmp$x[f0,c('010','011')]%*%c(1,1));
tmp$ZF_gain_GF2_loss_num    = sum(tmp$x[f0,c('001','011')]%*%c(1,1));
tmp$ZF_gain_GF_loss_num     = tmp$ZF_gain_GF1_loss_num + tmp$ZF_gain_GF2_loss_num ;
tmp$ZF_gain_GF_loss_num_WGD = sum(apply(tmp$x[f0,c('010','001')], 1, sum));
tmp$ZF_loss_GF_gain_num     = sum(tmp$x[f0,c('110','101','100')]%*%c(1,1,2));
tmp$ZF_loss_GF_gain_num_WGD = sum(apply(tmp$x[f0,c('110','101')], 1, sum));

tmp$x = as.matrix(data1$exon_loss_len);
tmp$tot_len  = sum(apply(tmp$x[f0,], 1, sum));
tmp$config_len  = apply(tmp$x[f0,], 2, sum);
tmp$ZF_gain_GF1_loss_len    = sum(tmp$x[f0,c('010','011')]%*%c(1,1));
tmp$ZF_gain_GF2_loss_len    = sum(tmp$x[f0,c('001','011')]%*%c(1,1));
tmp$ZF_gain_GF_loss_len     = tmp$ZF_gain_GF1_loss_len + tmp$ZF_gain_GF2_loss_len ;
tmp$ZF_gain_GF_loss_len_WGD = sum(apply(tmp$x[f0,c('010','001')], 1, sum));
tmp$ZF_loss_GF_gain_len     = sum(tmp$x[f0,c('110','101','100')]%*%c(1,1,2));
tmp$ZF_loss_GF_gain_len_WGD = sum(apply(tmp$x[f0,c('110','101')], 1, sum));

tmp$config_perc  = tmp$config_len*100/tmp$tot_len;
tmp$ZF_gain_GF1_loss_perc    = tmp$ZF_gain_GF1_loss_len   /tmp$tot_len*100;
tmp$ZF_gain_GF2_loss_perc    = tmp$ZF_gain_GF2_loss_len   /tmp$tot_len*100;
tmp$ZF_gain_GF_loss_perc     = tmp$ZF_gain_GF_loss_len    /2/tmp$tot_len*100;
tmp$ZF_gain_GF_loss_perc_WGD = tmp$ZF_gain_GF_loss_len_WGD/2/tmp$tot_len*100;
tmp$ZF_loss_GF_gain_perc     = tmp$ZF_loss_GF_gain_len    /2/tmp$tot_len*100;
tmp$ZF_loss_GF_gain_perc_WGD = tmp$ZF_loss_GF_gain_len_WGD/2/tmp$tot_len*100;

tmp$y = data.frame(
		Count = c(tmp$tot_num, tmp$config_num, tmp$ZF_gain_GF1_loss_num, tmp$ZF_gain_GF2_loss_num, tmp$ZF_gain_GF_loss_num, tmp$ZF_gain_GF_loss_num_WGD, tmp$ZF_loss_GF_gain_num, tmp$ZF_loss_GF_gain_num_WGD), 
		Length= c(tmp$tot_len, tmp$config_len, tmp$ZF_gain_GF1_loss_len, tmp$ZF_gain_GF2_loss_len, tmp$ZF_gain_GF_loss_len, tmp$ZF_gain_GF_loss_len_WGD, tmp$ZF_loss_GF_gain_len, tmp$ZF_loss_GF_gain_len_WGD), 
		Percent = sprintf('%.2f', c(100, tmp$config_perc, tmp$ZF_gain_GF1_loss_perc, tmp$ZF_gain_GF2_loss_perc, tmp$ZF_gain_GF_loss_perc, tmp$ZF_gain_GF_loss_perc_WGD, tmp$ZF_loss_GF_gain_perc, tmp$ZF_loss_GF_gain_perc_WGD) )
		);
rownames(tmp$y) = c('Total', names(tmp$config_num), 'ZF gain GF1 loss', 'ZF gain GF2 loss', 'ZF gain GF loss', 'ZF gain GF (singleton gain)', 'ZF loss GF gain', 'ZF loss GF gain (singleton gain)');
write.table(tmp$y, file=paste0("output/ZF_GF.triplet_exon_loss_stat.txt"), sep="\t", quote=F);
# }}}

#############################################
# stat of CNE loss
#############################################
# {{{
tmp=list();
tmp$GF1_GF2_reciprocal_triplet_num  = sum( apply(data1$CNE5k_loss_num[,c('001','101')],1,sum)>0 & apply(data1$CNE5k_loss_num[,c('010','110')],1,sum)>0 );
tmp$ZF_gain_GF1_GF2_reciprocal_triplet_num  = sum( data1$CNE5k_loss_num[,'101']>0 & data1$CNE5k_loss_num[,'110']>0 );

tmp$x = as.matrix(data1$CNE5k_loss_num);
tmp$tot_num  = sum(apply(tmp$x[f0,], 1, sum));
tmp$config_num  = apply(tmp$x[f0,], 2, sum);
tmp$ZF_gain_GF1_loss_num    = sum(tmp$x[f0,c('010','011')]%*%c(1,1));
tmp$ZF_gain_GF2_loss_num    = sum(tmp$x[f0,c('001','011')]%*%c(1,1));
tmp$ZF_gain_GF_loss_num     = tmp$ZF_gain_GF1_loss_num + tmp$ZF_gain_GF2_loss_num ;
tmp$ZF_gain_GF_loss_num_WGD = sum(apply(tmp$x[f0,c('010','001')], 1, sum));
tmp$ZF_loss_GF_gain_num     = sum(tmp$x[f0,c('110','101','100')]%*%c(1,1,2));
tmp$ZF_loss_GF_gain_num_WGD = sum(apply(tmp$x[f0,c('110','101')], 1, sum));

tmp$x = as.matrix(data1$CNE5k_loss_len);
tmp$tot_len  = sum(apply(tmp$x[f0,], 1, sum));
tmp$config_len  = apply(tmp$x[f0,], 2, sum);
tmp$ZF_gain_GF1_loss_len    = sum(tmp$x[f0,c('010','011')]%*%c(1,1));
tmp$ZF_gain_GF2_loss_len    = sum(tmp$x[f0,c('001','011')]%*%c(1,1));
tmp$ZF_gain_GF_loss_len     = tmp$ZF_gain_GF1_loss_len + tmp$ZF_gain_GF2_loss_len ;
tmp$ZF_gain_GF_loss_len_WGD = sum(apply(tmp$x[f0,c('010','001')], 1, sum));
tmp$ZF_loss_GF_gain_len     = sum(tmp$x[f0,c('110','101','100')]%*%c(1,1,2));
tmp$ZF_loss_GF_gain_len_WGD = sum(apply(tmp$x[f0,c('110','101')], 1, sum));

tmp$config_perc  = tmp$config_len*100/tmp$tot_len;
tmp$ZF_gain_GF1_loss_perc    = tmp$ZF_gain_GF1_loss_len   /tmp$tot_len*100;
tmp$ZF_gain_GF2_loss_perc    = tmp$ZF_gain_GF2_loss_len   /tmp$tot_len*100;
tmp$ZF_gain_GF_loss_perc     = tmp$ZF_gain_GF_loss_len    /tmp$tot_len/2*100;
tmp$ZF_gain_GF_loss_perc_WGD = tmp$ZF_gain_GF_loss_len_WGD/tmp$tot_len/2*100;
tmp$ZF_loss_GF_gain_perc     = tmp$ZF_loss_GF_gain_len    /tmp$tot_len/2*100;
tmp$ZF_loss_GF_gain_perc_WGD = tmp$ZF_loss_GF_gain_len_WGD/tmp$tot_len/2*100;

tmp$y = data.frame(
		Count = c(tmp$tot_num, tmp$config_num, tmp$ZF_gain_GF1_loss_num, tmp$ZF_gain_GF2_loss_num, tmp$ZF_gain_GF_loss_num, tmp$ZF_gain_GF_loss_num_WGD, tmp$ZF_loss_GF_gain_num, tmp$ZF_loss_GF_gain_num_WGD), 
		Length= c(tmp$tot_len, tmp$config_len, tmp$ZF_gain_GF1_loss_len, tmp$ZF_gain_GF2_loss_len, tmp$ZF_gain_GF_loss_len, tmp$ZF_gain_GF_loss_len_WGD, tmp$ZF_loss_GF_gain_len, tmp$ZF_loss_GF_gain_len_WGD), 
		Percent = sprintf('%.2f', c(100, tmp$config_perc, tmp$ZF_gain_GF1_loss_perc, tmp$ZF_gain_GF2_loss_perc, tmp$ZF_gain_GF_loss_perc, tmp$ZF_gain_GF_loss_perc_WGD, tmp$ZF_loss_GF_gain_perc, tmp$ZF_loss_GF_gain_perc_WGD) )
		);
rownames(tmp$y) = c('Total', names(tmp$config_num), 'ZF gain GF1 loss', 'ZF gain GF2 loss', 'ZF gain GF loss', 'ZF gain GF (singleton gain)', 'ZF loss GF gain', 'ZF loss GF gain (singleton gain)');
write.table(tmp$y, file=paste0("output/ZF_GF.triplet_CNE5k_loss_stat.txt"), sep="\t", quote=F);
# }}}


#############################################
# plot histogram of sequence identity
#############################################
# {{{
tmp=list();
jpeg("output/plot/ZF_GF.nucl_iden_histogram.jpg", res=res, w=res*6, h=res*4);
tmp$a1  = data1$nucl_align$iden[,2];
tmp$a2  = data1$nucl_align$iden[,3];
tmp$a12 = data1$nucl_align$iden[,6];
tmp$b1  = apply(data1$nucl_align$iden[,2:3],1,min);
tmp$b2  = apply(data1$nucl_align$iden[,2:3],1,max);
res=300;
xlab='Nucleotide identity (%)';
ylab='Frequency (%)';
col=paste0(color[1:4], '60');
hist(c(tmp$a12[tmp$a12>0]), col=col[1], border=F, breaks=seq(60,100,1), freq=F, xlab=xlab, ylab=ylab, main='');
#hist(c(tmp$a1[tmp$a1>0],tmp$a2[tmp$a2>0]), col=col[2], border=F, breaks=seq(60,100,1), add=T, freq=F);
hist(c(tmp$b1[tmp$b1>0]), col=col[3], border=F, breaks=seq(60,100,1), add=T, freq=F);
hist(c(tmp$b2[tmp$b2>0]), col=col[4], border=F, breaks=seq(60,100,1), add=T, freq=F);
cx=par('cxy')[1];  cy=par('cxy')[2];
legend(x=par('usr')[1]+2*cx, y=par('usr')[4]-2*cy, legend=c('GF1-GF2', 'min(ZF-GF1,ZF-GF2)', 'max(ZF-GF1,ZF-GF2)'), col=col, border=col, fill=col, bg='grey', box.col='grey', xpd=T);
dev.off();

jpeg("output/plot/ZF_GF.exon_perc_diff_histogram.jpg", res=res, w=res*6, h=res*4);
tmp$a1  = data1$exon_perc_diff[,2];
tmp$a2  = data1$exon_perc_diff[,3];
tmp$a   = c(tmp$a1,tmp$a2);
tmp$b1  = apply(data1$exon_perc_diff[,2:3], 1, min);
tmp$b2  = apply(data1$exon_perc_diff[,2:3], 1, max);
tmp$a12 = data1$exon_perc_diff[,6]; # between GF1-GF2
tmp$a3  = apply(data1$exon_loss_perc[,c('011', '100')], 1, sum); # between ZF-WGD
col=paste0(color[1:4], '60');
xlab='Exon gain/loss (%)';
ylab='Frequency (%)';
#tmp$a012 = data1$exon_loss_perc[,c('001')] + data1$exon_loss_perc[,c('010')];
#hist(c(tmp$a1[tmp$a1>=0],tmp$a2[tmp$a2>=0]), col='#0000FF80', border=F, breaks=seq(0,100,2), freq=F);
hist(c(tmp$a12[tmp$a12>=0]), col=col[1], border=F, breaks=seq(0,100,2), freq=F, xlab=xlab, ylab=ylab, main='');
hist(c(tmp$b1), col=col[2], border=F, breaks=seq(0,100,2), freq=F, add=T);
hist(c(tmp$b2), col=col[3], border=F, breaks=seq(0,100,2), freq=F, add=T);
#hist(c(tmp$a[tmp$a>=0]), col=col[2], border=F, breaks=seq(0,100,2), add=T, freq=F);
#hist(c(tmp$a3[tmp$a3>=0]), col='#0000FF80', border=F, breaks=seq(0,100,2), add=T, freq=F);
#hist(c(tmp$a012[tmp$a012>=0]), col='#FF000080', border=F, breaks=seq(0,100,2), add=T, freq=F);
legend(x=par('usr')[2]-2*cx, y=par('usr')[4]-2*cy, legend=c('GF1-GF2', 'min(ZF-GF1,ZF-GF2)', 'max(ZF-GF1,ZF-GF2)'), col=col, border=col, fill=col, bg='grey', box.col='grey', xpd=T, xjust=1.0);
dev.off();

jpeg("output/plot/ZF_GF.exon_perc_diff_histogram.ZF_always_exists.jpg", res=res, w=res*6, h=res*4);
tmp$tot = apply(data1$exon_loss_perc[,c('000','001','010','011')], 1, sum);
tmp$a1  = apply(data1$exon_loss_perc[,c('010', '011')], 1, sum);
tmp$a2  = apply(data1$exon_loss_perc[,c('001', '011')], 1, sum);
tmp$a12 = apply(data1$exon_loss_perc[,c('001', '010')], 1, sum);
tmp$b1  = apply(cbind(tmp$a1,tmp$a2), 1, min);
tmp$b2  = apply(cbind(tmp$a1,tmp$a2), 1, max);
col=paste0(color[1:4], '60');
xlab='Exon gain/loss (%)';
ylab='Frequency (%)';
hist(c(tmp$a12[tmp$a12>=0]), col=col[1], border=F, breaks=seq(0,100,2), freq=F, xlab=xlab, ylab=ylab, main='');
hist(c(tmp$b1), col=col[2], border=F, breaks=seq(0,100,2), freq=F, add=T);
hist(c(tmp$b2), col=col[3], border=F, breaks=seq(0,100,2), freq=F, add=T);
legend(x=par('usr')[2]-2*cx, y=par('usr')[4]-2*cy, legend=c('GF1-GF2', 'min(ZF-GF1,ZF-GF2)', 'max(ZF-GF1,ZF-GF2)'), col=col, border=col, fill=col, bg='grey', box.col='grey', xpd=T, xjust=1.0);
dev.off();

jpeg("output/plot/ZF_GF.exon_num_diff_histogram.jpg", res=res, w=res*6, h=res*41);
tmp$a1  = data1$exon_num_diff[,2];
tmp$a2  = data1$exon_num_diff[,3];
tmp$a12 = data1$exon_num_diff[,6];
tmp$max = 20;
hist(c(tmp$a1[tmp$a1>=0 & tmp$a1<=20],tmp$a2[tmp$a2>=0 & tmp$a2<=20]), col='#0000FF80', border=F, breaks=seq(0,tmp$max,1), freq=F);
hist(c(tmp$a12[tmp$a12>=0 & tmp$a12<=20]), col='#FF000080', border=F, breaks=seq(0,tmp$max,1), add=T, freq=F);
dev.off();

jpeg("output/plot/ZF_GF.CNE5k_perc_diff_histogram.jpg", res=res, w=res*6, h=res*4);
tmp=list();
tmp$a1  = data1$CNE5k_perc_diff[,2];
tmp$a2  = data1$CNE5k_perc_diff[,3];
tmp$a   = c(tmp$a1,tmp$a2);
tmp$b1  = apply(data1$CNE5k_perc_diff[,2:3], 1, min);
tmp$b2  = apply(data1$CNE5k_perc_diff[,2:3], 1, max);
tmp$a12 = data1$CNE5k_perc_diff[,6]; # between GF1-GF2
tmp$a3  = apply(data1$CNE5k_loss_perc[,c('011', '100')], 1, sum); # between ZF-WGD
col=paste0(color[1:4], '60');
xlab='CNE gain/loss (%)';
ylab='Frequency (%)';
#tmp$a012 = data1$CNE5k_loss_perc[,c('001')] + data1$CNE5k_loss_perc[,c('010')];
#hist(c(tmp$a1[tmp$a1>=0],tmp$a2[tmp$a2>=0]), col='#0000FF80', border=F, breaks=seq(0,100,2), freq=F);
hist(c(tmp$a12[tmp$a12>=0]), col=col[1], border=F, breaks=seq(0,100,2), freq=F, xlab=xlab, ylab=ylab, main='');
hist(c(tmp$b1), col=col[2], border=F, breaks=seq(0,100,2), freq=F, add=T);
hist(c(tmp$b2), col=col[3], border=F, breaks=seq(0,100,2), freq=F, add=T);
#hist(c(tmp$a[tmp$a>=0]), col=col[2], border=F, breaks=seq(0,100,2), add=T, freq=F);
#hist(c(tmp$a3[tmp$a3>=0]), col='#0000FF80', border=F, breaks=seq(0,100,2), add=T, freq=F);
#hist(c(tmp$a012[tmp$a012>=0]), col='#FF000080', border=F, breaks=seq(0,100,2), add=T, freq=F);
legend(x=par('usr')[2]-2*cx, y=par('usr')[4]-2*cy, legend=c('GF1-GF2', 'min(ZF-GF1,ZF-GF2)', 'max(ZF-GF1,ZF-GF2)'), col=col, border=col, fill=col, bg='grey', box.col='grey', xpd=T, xjust=1.0);
dev.off();

jpeg("output/plot/ZF_GF.CNE5k_perc_diff_histogram.ZF_always_exists.jpg", res=res, w=res*6, h=res*4);
tmp=list();
tmp$tot = apply(data1$CNE5k_loss_perc[,c('000','001','010','011')], 1, sum);
tmp$a1  = apply(data1$CNE5k_loss_perc[,c('010', '011')], 1, sum);
tmp$a2  = apply(data1$CNE5k_loss_perc[,c('001', '011')], 1, sum);
tmp$a12 = apply(data1$CNE5k_loss_perc[,c('001', '010')], 1, sum);
tmp$b1  = apply(cbind(tmp$a1,tmp$a2), 1, min);
tmp$b2  = apply(cbind(tmp$a1,tmp$a2), 1, max);
col=paste0(color[1:4], '60');
xlab='CNE gain/loss (%)';
ylab='Frequency (%)';
hist(c(tmp$a12[tmp$a12>=0]), col=col[1], border=F, breaks=seq(0,100,2), freq=F, xlab=xlab, ylab=ylab, main='');
hist(c(tmp$b1), col=col[2], border=F, breaks=seq(0,100,2), freq=F, add=T);
hist(c(tmp$b2), col=col[3], border=F, breaks=seq(0,100,2), freq=F, add=T);
legend(x=par('usr')[2]-2*cx, y=par('usr')[4]-2*cy, legend=c('GF1-GF2', 'min(ZF-GF1,ZF-GF2)', 'max(ZF-GF1,ZF-GF2)'), col=col, border=col, fill=col, bg='grey', box.col='grey', xpd=T, xjust=1.0);
dev.off();
# }}}


# count ohnolog cluster for each fate
# {{{
thres2=2;
f0 = data1$f & data1$exp$fpkm$max2[,1]>=thres & data1$exp$tpm$max2[,1]>=thres
for (etype in c('fpkm','tpm')) {
	fpkm  = data1[[etype]];
	lfpkm = data1$exp[[etype]]$log2;
	data1$exp[[etype]]$fate = list(
		GF = fate_classify(fpkm, nts, 1,2,3, thres, thres, thres2, fc=2, C1=0.6, C3=0.75, cor_p=0.1, tp=0.01)
		);
	all_fates = names(data1$exp[[etype]]$fate$GF);

	tot=sum(f0);
	data1$exp[[etype]]$fate_n = list();
	sp='GF';
	for (ft in names(data1$exp[[etype]]$fate[[sp]])) {
		data1$exp[[etype]]$fate_n[[sp]][[ft]] = sum(data1$exp[[etype]]$fate[[sp]][[ft]] & f0);
	}

	a=data1$exp[[etype]]$fate_n;
	x=matrix(0, length(a$GF), 4);
	rownames(x)=names(a$GF);
	for (ft in names(a$GF)) {
		x[ft,1] = a$GF[ft];
	}
	x
	x/tot

	y=matrix(0, length(all_fates), 4);
	colnames(y) = c("Triplets", "Percentage", "Mean identity (ZF~GF)", "Mean Identity (GF1~GF2)");
	rownames(y) = all_fates;
	for (ft in all_fates) {
		y[ft,1] = data1$exp[[etype]]$fate_n$GF[ft];
		y[ft,2] = sprintf("%.2f", data1$exp[[etype]]$fate_n$GF[ft]*100/tot);
		ft1=paste0(ft,'1');
		ft2=paste0(ft,'1');
		if (ft=='nonfunc' || ft=='cor_neofunc') {
			a = c(data1$nucl_align$iden[data1$exp[[etype]]$fate$GF[[ft1]] & f0, 2],   data1$nucl_align$iden[data1$exp[[etype]]$fate$GF[[ft2]] & f0, 3])
		} else if (regexpr('1$', ft)>0){
			a = data1$nucl_align$iden[data1$exp[[etype]]$fate$GF[[ft]] & f0, 2];
		} else if (regexpr('2$', ft)>0){
			a = data1$nucl_align$iden[data1$exp[[etype]]$fate$GF[[ft]] & f0, 3];
		} else {
			a = apply(data1$nucl_align$iden[data1$exp[[etype]]$fate$GF[[ft]] & f0, 2:3], 1, min);
		}
		a[a=0]=70;
		y[ft,3] = sprintf("%.2f", mean( a ) );
#		a1 = apply(data1$nucl_align$cov1[data1$fate$GF[[ft]] & data1$f, 2:3], 1, min);
#		a2 = apply(data1$nucl_align$cov2[data1$fate$GF[[ft]] & data1$f, 2:3], 1, min);
#		a = apply(cbind(a1,a2),1,min);
#		a[a=0]=50;
#		y[ft,4] = sprintf("%.2f", mean(a) );
		a = c(data1$nucl_align$iden[data1$exp[[etype]]$fate$GF[[ft]] & f0, 6]);
		a[a=0]=50;
		y[ft,4] = sprintf("%.2f", mean(a) );
	}
	write.table(y, file=paste0("output/ZF_GF.",etype,".fate_identity.txt"), sep="\t", quote=F);
}


# }}}

# bootstrap 1000 1
# {{{
lfpkm1000 = list();
for (k in 1:1000) {
	lfpkm1000[[k]] = random_data1(data1$lfpkm[,1:(nts*tot_copy)], replace=F);
	sum_tpm2 = matrix(0, data1$n, tot_sp*nts);
	colnames(sum_tpm2) = paste0("sum_of_", rep(sps,each=nts), '.', rep(tissues,tot_sp));
	j0 = 0;
	for (i in 1:tot_sp) {
		for (j in 1:copies[sps[i]]) {
			j0 = j0+1;
			sum_tpm2[,(nts*(i-1)+1):(i*nts)] = sum_tpm2[,(nts*(i-1)+1):(i*nts)] + lfpkm1000[[k]][,(nts*(j0-1)+1):(j0*nts)];
		}
	}
	lfpkm1000[[k]] = cbind(lfpkm1000[[k]], sum_tpm2);
}
fate1000 = list(
	CC_double_conserved=list(),
	GF_double_conserved=list(),
	CC_subfunc=list(),
	GF_subfunc=list(),
	CC_neofunc=list(),
	GF_neofunc=list()
		);
for (k in 1:1000) {
	tmp_lfpkm = random_data(data$lfpkm[,1:(nts*tot_copy)], replace=F);
	tmp_sum = ohno_sum(tmp_lfpkm, tissues, sps2);
	tmp_lfpkm = cbind(tmp_lfpkm, tmp_sum);
	C1=0.6; C2=0.75;
	fate1000$CC_double_conserved[[k]] = is_double_conserved(lfpkm1000[[k]], 1,2,3,7,lthres, C1);
	fate1000$GF_double_conserved[[k]] = is_double_conserved(lfpkm1000[[k]], 1,4,5,8,lthres, C1);
	
	fate1000$CC_subfunc[[k]] = is_subfunc(lfpkm1000[[k]], 1,2,3,7,lthres, C1,C2);
	fate1000$GF_subfunc[[k]] = is_subfunc(lfpkm1000[[k]], 1,4,5,8,lthres, C1,C2);

	fate1000$CC_neofunc[[k]] = is_neofunc(lfpkm1000[[k]], 1,2,3,7,lthres, C1,C2);
	fate1000$GF_neofunc[[k]] = is_neofunc(lfpkm1000[[k]], 1,4,5,8,lthres, C1,C2);

	if (k%%100==0) { print(k); }
}

fate1000_n = list();
for (sp in c("CC", "GF")) {
for (ft in c("double_conserved", "subfunc", "neofunc")) {
	name = paste0(sp,'_',ft);
	fate1000_n[[name]] = sapply(fate1000[[name]], function(z) { sum(z$res) });
}
}

data$fate_p = c();
for (sp in c("CC", "GF")) {
for (ft in c("double_conserved", "subfunc", "neofunc")) {
	name = paste0(sp,'_',ft);
	data$fate_p[name] = sum(fate1000_n[[name]]>=data$fate_n[[name]])/1000;
}}
# }}}

# bootstrap 1000 2
# {{{
fate1000_2 = list(
	CC_double_conserved=list(),
	GF_double_conserved=list(),
	CC_subfunc=list(),
	GF_subfunc=list(),
	CC_neofunc=list(),
	GF_neofunc=list()
		);
for (k in 1:1000) {
	C1=0.6; C2=0.75;
	tmp_lfpkm = random_data2(data$lfpkm[,1:(nts*tot_copy)], win=nts, replace=F);
	tmp_sum = ohno_sum(tmp_lfpkm, tissues, sps2);
	tmp_lfpkm = cbind(tmp_lfpkm, tmp_sum);

	fate1000_2$CC_double_conserved[[k]] = is_double_conserved(tmp_lfpkm, 1,2,3,7,lthres, C1);
	fate1000_2$GF_double_conserved[[k]] = is_double_conserved(tmp_lfpkm, 1,4,5,8,lthres, C1);
	fate1000_2$CC_subfunc[[k]] = is_subfunc(tmp_lfpkm, 1,2,3,7,lthres, C1,C2);
	fate1000_2$GF_subfunc[[k]] = is_subfunc(tmp_lfpkm, 1,4,5,8,lthres, C1,C2);
	fate1000_2$CC_neofunc[[k]] = is_neofunc(tmp_lfpkm, 1,2,3,7,lthres, C1,C2);
	fate1000_2$GF_neofunc[[k]] = is_neofunc(tmp_lfpkm, 1,4,5,8,lthres, C1,C2);

	if (k%%100==0) { print(k); }
}
fate1000_2_n = list();
for (sp in c("CC", "GF")) {
for (ft in c("double_conserved", "subfunc", "neofunc")) {
	name = paste0(sp,'_',ft);
	fate1000_2_n[[name]] = sapply(fate1000_2[[name]], function(z) { sum(z$res) });
}
}

data$fate_p2 = c();
for (sp in c("CC", "GF")) {
for (ft in c("double_conserved", "subfunc", "neofunc")) {
	name = paste0(sp,'_',ft);
	data$fate_p2[name] = sum(!is.na(fate1000_2_n[[name]]) & fate1000_2_n[[name]]>=data$fate_n[[name]])/(1000-sum(is.na(fate1000_2_n[[name]])));
}}
# }}}

save(data1, file='Rdata/ZF_GF_data.Rdata');

############################################################
# draw 2D-histogram between: ZF GF CC
############################################################
# {{{
# 1. cor(ZF,CC1) and cor(ZF,CC2)
# 2. cor(ZF,GF1) and cor(ZF,GF2)
# 3. cor(ZF,CC1) and cor(ZF,GF1)
# 4. cor(ZF,CC2) and cor(ZF,GF2)
# 5. cor(ZF,CC1+CC2) and cor(ZF,GF1+GF2)
# 6. cor(CC1,CC2) and cor(GF1,GF2)
res=300;
jpeg("output/plot/cor_hist2d.jpg", res=res, w=res*9, h=res*11);
par(mfrow=c(3,2)); par(mar=c(5,5,1,6)+0.1);
for (k in 1:6) {
	if (k==1) {
		x = data$pair_cor[,2];   y = data$pair_cor[,3];
		xlab = "Correlation of ZF and CC1";
		ylab = "Correlation of ZF and CC2";
	} else if (k==2) {
		x = data$pair_cor[,4];   y = data$pair_cor[,5];
		xlab = "Correlation of ZF and GF1";
		ylab = "Correlation of ZF and GF2";
	} else if (k==3) {
		x = data$pair_cor[,2];   y = data$pair_cor[,4];
		xlab = "Correlation of ZF and CC1";
		ylab = "Correlation of ZF and GF1";
	} else if (k==4) {
		x = data$pair_cor[,3];   y = data$pair_cor[,5];
		xlab = "Correlation of ZF and CC2";
		ylab = "Correlation of ZF and GF2";
	} else if (k==5) {
		x = data$pair_cor[,7];   y = data$pair_cor[,8];
		xlab = "Correlation of ZF and CC1+CC2";
		ylab = "Correlation of ZF and GF1+GF2";
	} else if (k==6) {
		x = data$pair_cor[,(2-1)*8+3];   y = data$pair_cor[,(4-1)*8+5];
		xlab = "Correlation of CC1 and CC2";
		ylab = "Correlation of GF1 and GF2";
	}
	n = length(x);
	n11 = sum(x<0.6 & y<0.6);    p11 = sprintf("%.2f", n11*100/n);
	n12 = sum(x<0.6 & y>=0.6);   p12 = sprintf("%.2f", n12*100/n);
	n21 = sum(x>=0.6 & y<0.6);   p21 = sprintf("%.2f", n21*100/n);
	n22 = sum(x>=0.6 & y>=0.6);  p22 = sprintf("%.2f", n22*100/n);
	x = c(x,-1);  y = c(y,-1);
	a = hist2d(x,y, nbin=40, col=blueyellow255, xlab=xlab, ylab=ylab, same.scale=T, cex=1.3, cex.axis=1.3, cex.lab=1.6);
	for (i in 1:255) {
		rect(1.1,-1+1.0*(i-1)/256,1.2,-1+1.0*i/256, col=blueyellow255[i], border=blueyellow255[i], xpd=T);
	}
	text(1.22, -1  , sprintf("%.1f%%",min(a$counts)*100/data$n), xpd=T, adj=0, cex=1.4);
	text(1.22, -0.5, sprintf("%.1f%%",(min(a$counts)+max(a$counts))*100/2/data$n), xpd=T, adj=0, cex=1.4);
	text(1.22, 0   , sprintf("%.1f%%",max(a$counts)*100/data$n), xpd=T, adj=0, cex=1.4);
#	image(c(1.05,1.15), -1+0.5*(0:256)/256, matrix(0:255,1,256), col=blueyellow255, xpd=T, add=T);
	abline(h=0.6, col='gray');   abline(v=0.6, col='gray');
	text(-0.2,-0.2, paste0(n11, "\n(", p11, "%)"), xpd=T, cex=1.5, col='orange');
	text(-0.2,0.8, paste0(n12, "\n(", p12, "%)"), xpd=T, cex=1.5, col='orange');
	text(0.8,-0.2, paste0(n21, "\n(", p21, "%)"), xpd=T, cex=1.5, col='orange');
	text(1.25,0.8, paste0(n22, "\n(", p22, "%)"), xpd=T, cex=1.5, col='orange');
	text(0.7,-0.9, 0.6, xpd=T, col='gray', cex=1.4);   text(-0.9,0.66, 0.6, xpd=T, col='gray', cex=1.4);
}
dev.off();
#x1 = data$pair_cor[,2];
#x2 = data$pair_cor[,3];
#y1 = data$pair_cor[,4];
#y2 = data$pair_cor[,5];
#dxy1 = abs(x1-y1)+abs(x2-y2);
#dxy2 = abs(x1-y2)+abs(x2-y1);
#x = c(x1,x2);
#y3 = y1;
#y4 = y2;
#y3[dxy2<dxy1] = y2[dxy2<dxy1];
#y4[dxy2<dxy1] = y1[dxy2<dxy1];
#y = c(y3,y4);
#hist2d(x,y, nbin=20, col=rgb(1,1,(255:0)/255))
# }}}

#-----------------------------------------------
# 20190104 plot corr ~ identity , corr ~ coverage
#-----------------------------------------------
# {{{
tmp=list();
for (k in 1:2) {
	if (k==1) { # corr ~ iden(GF1,GF2)
		tmp$sep = c(0,seq(80,98,2));
		tmp$x = data1$nucl_align$iden[,6];
		tmp$y = data1$exp$fpkm$pair_cor[,8];
		tmp$f = f0 & data1$exp$fpkm$sd2[,2]>0 & data1$exp$fpkm$sd2[,3]>0;
		jpeg_name = "output/plot/ZF_GF.iden_cor_GF1_GF2.jpg";
		jpeg_name2 = "output/plot/ZF_GF.iden_cor_GF1_GF2.boxplot.jpg";
		xlab = 'Nuclieotide Identity between GF pair';
	} else if (k==2) { # corr ~ iden(ZF,GF)
		tmp$sep = c(0,seq(70,98,2));
		tmp$x = c(data1$nucl_align$iden[,2], data1$nucl_align$iden[,3]);
		tmp$y = c(data1$exp$fpkm$pair_cor[,2], data1$exp$fpkm$pair_cor[,3]);
		tmp$f = c(f0 & data1$exp$fpkm$sd2[1]>0 & data1$exp$fpkm$sd2[2]>0, f0 & data1$exp$fpkm$sd2[1]>0 & data1$exp$fpkm$sd2[3]>0);
		jpeg_name = "output/plot/ZF_GF.iden_cor_ZF_GF.jpg";
		jpeg_name2 = "output/plot/ZF_GF.iden_cor_ZF_GF.boxplot.jpg";
		xlab = 'Nuclieotide Identity between ZF and GF';
#	} else if (k==3) { # corr ~ cov(GF1,GF2)
#		tmp$sep = c(0,seq(30,100,5));
#		tmp$x1 = data1$nucl_align$cov1[,6];
#		tmp$x2 = data1$nucl_align$cov2[,6];
#		tmp$x = apply(cbind(tmp$x1,tmp$x2), 1, min);
#		tmp$y = data1$pair_cor[,8];
#		tmp$f = data1$f & data1$sd2[,2]>0 & data1$sd2[,3]>0;
#		jpeg_name = "output/plot/ZF_GF.cov_cor_ZF_GF.jpg";
#		jpeg_name2 = "output/plot/ZF_GF.cov_cor_ZF_GF.boxplot.jpg";
#		xlab = 'Nuclieotide Coverage between ZF and GF';
#	} else if (k==4) { # corr ~ cov(ZF,GF)
#		tmp$sep = c(0,seq(30,100,5));
#		tmp$x1 = c(data1$nucl_align$cov1[,2], data1$nucl_align$cov1[,3]);
#		tmp$x2 = c(data1$nucl_align$cov2[,2], data1$nucl_align$cov2[,3]);
#		tmp$x = apply(cbind(tmp$x1,tmp$x2), 1, min);
#		tmp$y = c(data1$pair_cor[,2], data1$pair_cor[,3]);
#		tmp$f = c(data1$f & data1$sd2[1]>0 & data1$sd2[2]>0, data1$f & data1$sd2[1]>0 & data1$sd2[3]>0);
#		jpeg_name = "output/plot/ZF_GF.cov_cor_GF1_GF2.jpg";
#		jpeg_name2 = "output/plot/ZF_GF.cov_cor_GF1_GF2.boxplot.jpg";
#		xlab = 'Nuclieotide Coverage between GF pair';
	}
	tmp$m = length(tmp$sep)-1;
	tmp$xx = matrix(0, tmp$m, 4);
	tmp$x1 = NULL;
	for (i in 1:tmp$m) {
		tmp$f1 = tmp$f & tmp$x>0 & tmp$x>tmp$sep[i] & tmp$x<=tmp$sep[i+1];
		tmp$xx[i,] = c( mean(tmp$x[tmp$f1]),  mean(tmp$y[tmp$f1]),   
				sd(tmp$y[tmp$f1])/sqrt(sum(tmp$f1)),    sum(tmp$f1));
		tmp$x1 = rbind(tmp$x1, cbind(i,tmp$y[tmp$f1]));
	}
	if (k==1 || k==2) { colnames(tmp$x1) = c("Iden", "Corr"); } 
	if (k==3 || k==4) { colnames(tmp$x1) = c("Cov", "Corr"); } 

	ylab='Average correlation';
	ylab2='Correlation';
	xtick = c(paste0('<',tmp$sep[2]), tmp$sep[3:(tmp$m+1)]);
	jpeg(jpeg_name, res=res, width=6*res, height=4*res);
	plot(1:tmp$m, tmp$xx[,2], type='l', ylim=c(min(tmp$xx[,2]-0.05),max(tmp$xx[,2]+0.05)), axes=F, xlab=xlab, ylab=ylab);
#axis(1,at=c(1:tmp$m), lab=c('<80', paste0(tmp$idens1[2:tmp$m],'-',tmp$idens1[3:(tmp$m+1)])) );
	axis(1,at=c(1,seq(2,tmp$m,2),tmp$m), lab=c(paste0('<',tmp$sep[2]), tmp$sep[seq(3,tmp$m+1,2)], 100) );
	axis(2);
	points(1:tmp$m, tmp$xx[,2], cex=log2(tmp$xx[,4])/4, pch=21, bg='#0000FF80', col='#0000FF');
	dev.off();

	jpeg(jpeg_name2, res=res, width=6*res, height=4*res);
	if (k==1 || k==2) {
		colnames(tmp$x1) = c("Iden", "Corr");
		boxplot( Corr ~ Iden, tmp$x1, names=xtick, ylim=c(0.4,1), col='#30C080', xlab=xlab, ylab=ylab2);
	} 
	if (k==3 || k==4) {
		colnames(tmp$x1) = c("Cov", "Corr");
		boxplot( Corr ~ Cov, tmp$x1, names=xtick, ylim=c(0.4,1), col='#30C080', xlab=xlab, ylab=ylab2);
	}
	dev.off()
}

#iden_group <- ceiling( (dataC$iden-86)/2 );
#iden_group[ iden_group<=0 ] = 0;
#iden_group[ iden_group>=5 ] = 5;
#dataC$iden_group <- iden_group;
#dataAB$iden_group <- rep(iden_group,2);
#xlab=c("<86", paste( 86+(0:3)*2, 86+(1:4)*2, sep="-"), ">94");
#png("/Users/chenz11/Documents/goldfish_project/tables and figures/WGD/iden_corr.boxplot.png", width=1500, height=1000, res=150)
#boxplot( corr ~ iden_group, dataC, names=xlab);
#dev.off();
#iden_corr.lm <- lm(corr~iden, dataC);
#summary(iden_corr.lm);
#
#corr_hist <- list();
#corr_hist.ymax <- 0;
#for (i in 1:6) {
#	corr_hist[[i]] <- hist(corr[iden_group==i-1], breaks=c(-1,0,0.2,0.4,0.6,0.8,1.0), plot=F);
#	max1 <- max(corr_hist[[i]]$density); 
#	if (max1 > corr_hist.ymax) { corr_hist.ymax=max1; }
#}
#corr_hist.percent <- corr_hist[[1]]$counts*100/sum(corr_hist[[1]]$counts);
#for (i in 2:6) {
#	corr_hist.percent <- cbind(corr_hist.percent, corr_hist[[i]]$counts*100/sum(corr_hist[[i]]$counts) );
#}
#colnames(corr_hist.percent)<-xlab;
#rownames(corr_hist.percent)<- c("<0", "0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0");
#
#euc_hist <- list();
#euc_hist.ymax <- 0;
#for (i in 1:6) {
#	euc_hist[[i]] <- hist(dataC$euc[dataC$iden_group==i-1], breaks=c(seq(0,8,2),max_euc), plot=F);
#	max1 <- max(euc_hist[[i]]$density); 
#	if (max1 > euc_hist.ymax) { euc_hist.ymax=max1; }
#}
#euc_hist.percent <- euc_hist[[1]]$counts*100/sum(euc_hist[[1]]$counts);
#for (i in 2:6) {
#	euc_hist.percent <- cbind(euc_hist.percent, euc_hist[[i]]$counts*100/sum(euc_hist[[i]]$counts) );
#}
#colnames(euc_hist.percent)<-xlab;
#rownames(euc_hist.percent)<-c("0-2", "2-4", "4-6", "6-8", ">=8");
#
#par(mar=c(3,3,3,10));
#plot(euc_hist[[1]]$mids, euc_hist[[1]]$density, ylim=c(0,max1), type='l', col=color[1]);
#k1 <- length(euc_hist[[1]]$mids);
#for (i in 1:6) {
#	polygon( c(euc_hist[[i]]$mids[1],euc_hist[[i]]$mids, euc_hist[[i]]$mids[k1]), c(0,euc_hist[[i]]$density,0), col=color8a[i], border=F );
#}
#for (i in 1:6) {
#	lines(euc_hist[[i]]$mids, euc_hist[[i]]$density, col=color8[i]);
#}
#legend(x=max_euc+1, y=max1, legend=xlab, text.col=2:7, fill=2:7, box.lwd=0, border=2:7, xpd=T)
#
#
#plot(corr_hist[[1]]$mids, corr_hist[[1]]$density, ylim=c(0,corr_hist.ymax), type='l', col=greenred(11)[1]);
#for (i in 2:5) {
#    lines(corr_hist[[i]]$mids, corr_hist[[i]]$density, type='l', col=greenred(11)[i]);
#}
# }}}

#-----------------------------------------------
# 20190104 corr ~ exon loss, corr ~ CNE loss
#-----------------------------------------------
# {{{
# exon difference by length percentage, cor(GF1,GF2)
# {{{
for (etype in c('fpkm','tpm')) {
	tmp=list();
	tmp$f = f0 & data1$exp[[etype]]$sd2[,2]>0 & data1$exp[[etype]]$sd2[,3]>0;
	tmp$pair_cor = data1$exp[[etype]]$pair_cor;
	tmp$exon_diff = data1$exon_perc_diff[,6];
	tmp$sep = c(-1, seq(0,50,5), 100);
	tmp$x = tmp$exon_diff;
	tmp$y = tmp$pair_cor[,8];
	tmp$m = length(tmp$sep)-1;
	tmp$xx = matrix(0, tmp$m, 4);
	tmp$x1 = NULL;
	for (i in 1:tmp$m) {
		tmp$f1 = tmp$f & tmp$x>tmp$sep[i] & tmp$x<=tmp$sep[i+1];
		tmp$xx[i,] = c( mean(tmp$x[tmp$f1]),  mean(tmp$y[tmp$f1]),   
				sd(tmp$y[tmp$f1])/sqrt(sum(tmp$f1)),    sum(tmp$f1));
		tmp$x1 = rbind(tmp$x1, cbind(i,tmp$y[tmp$f1]));
	}
	colnames(tmp$x1) = c('ExonDiff', 'Corr');

	jpeg_name = paste0("output/plot/ZF_GF.",etype,".exon_diff_cor_GF1_GF2.jpg");
	jpeg_name2 = paste0("output/plot/ZF_GF.",etype,".exon_diff_cor_GF1_GF2.boxplot.jpg");
	xlab='Exon difference between GF gene pairs (%)';
	ylab='Average correlation';
	ylab2='Correlation';
	xtick = c(paste0('=',tmp$sep[2]), tmp$sep[3:(tmp$m+1)]);
# line plot with size
	jpeg(jpeg_name, res=res, width=6*res, height=4*res);
	plot(1:tmp$m, tmp$xx[,2], type='l', ylim=c(min(tmp$xx[,2]-0.05),max(tmp$xx[,2]+0.05)), axes=F, xlab=xlab, ylab=ylab);
#axis(1,at=c(1:tmp$m), lab=c('<80', paste0(tmp$idens1[2:tmp$m],'-',tmp$idens1[3:(tmp$m+1)])) );
	axis(1,at=c(1:tmp$m), lab=xtick );
	axis(2);
	points(1:tmp$m, tmp$xx[,2], cex=log2(tmp$xx[,4])/4, pch=21, bg='#0000FF80', col='#0000FF');
	dev.off();
# boxplot
	jpeg(jpeg_name2, res=res, width=6*res, height=4*res);
	colnames(tmp$x1) = c("ExonDiff", "Corr");
	boxplot( Corr ~ ExonDiff, tmp$x1, names=xtick, ylim=c(0.4,1), col='#30C080', xlab=xlab, ylab=ylab2);
	dev.off();
}
# }}}

# exon difference by length percentage, cor(ZF,GF)
# {{{
for (etype in c('fpkm','tpm')) {
	tmp=list();
	tmp$f = c( f0 & data1$exp[[etype]]$sd2[,1]>0 & data1$exp[[etype]]$sd2[,2]>0, f0 & data1$exp[[etype]]$sd2[,1]>0 & data1$exp[[etype]]$sd2[,3]>0);
	tmp$pair_cor = data1$exp[[etype]]$pair_cor;
	tmp$sep = c(-1, seq(0,50,5), 100);
	tmp$x1 = data1$exon_perc_diff[,2];
	tmp$x2 = data1$exon_perc_diff[,3];
	tmp$x = c(tmp$x1,tmp$x2);
	tmp$y = c(tmp$pair_cor[,2], tmp$pair_cor[,3]);
	tmp$m = length(tmp$sep)-1;
	tmp$xx = matrix(0, tmp$m, 4);
	tmp$x1 = NULL;
	for (i in 1:tmp$m) {
		tmp$f1 = tmp$f & tmp$x>tmp$sep[i] & tmp$x<=tmp$sep[i+1];
		tmp$xx[i,] = c( mean(tmp$x[tmp$f1]),  mean(tmp$y[tmp$f1]),   
				sd(tmp$y[tmp$f1])/sqrt(sum(tmp$f1)),    sum(tmp$f1));
		tmp$x1 = rbind(tmp$x1, cbind(i,tmp$y[tmp$f1]));
	}
	colnames(tmp$x1) = c('ExonDiff', 'Corr');

	jpeg_name = paste0("output/plot/ZF_GF.",etype, ".exon_diff_cor_ZF_GF.jpg");
	jpeg_name2 = paste0("output/plot/ZF_GF.",etype,".exon_diff_cor_ZF_GF.boxplot.jpg");
	xlab='Exon difference between ZF and GF (%)';
	ylab='Average correlation';
	ylab2='Correlation';
	xtick = c(paste0('=',tmp$sep[2]), tmp$sep[3:(tmp$m+1)]);
# line plot with size
	jpeg(jpeg_name, res=res, width=6*res, height=4*res);
	plot(1:tmp$m, tmp$xx[,2], type='l', ylim=c(min(tmp$xx[,2]-0.05),max(tmp$xx[,2]+0.05)), axes=F, xlab=xlab, ylab=ylab);
#axis(1,at=c(1:tmp$m), lab=c('<80', paste0(tmp$idens1[2:tmp$m],'-',tmp$idens1[3:(tmp$m+1)])) );
	axis(1,at=c(1:tmp$m), lab=c(paste0('<',tmp$sep[2]), tmp$sep[3:(tmp$m+1)]) );
	axis(2);
	points(1:tmp$m, tmp$xx[,2], cex=log2(tmp$xx[,4])/4, pch=21, bg='#0000FF80', col='#0000FF');
	dev.off();
# boxplot
	jpeg(jpeg_name2, res=res, width=6*res, height=4*res);
	colnames(tmp$x1) = c("ExonDiff", "Corr");
	boxplot( Corr ~ ExonDiff, tmp$x1, names=xtick, ylim=c(0.3,1), col='#30C080', xlab=xlab, ylab=ylab2, axes=F);
	axis(1,at=c(1:tmp$m), lab=xtick, xpd=T );
	axis(2);
	dev.off();
}
# }}}

# CNE difference by length percentage, cor(GF1,GF2)
# {{{
for (etype in c('fpkm','tpm')) {
	fpkm = data1[[etype]];
	tmp=list();
	tmp$exon_diff = data1$exon_perc_diff[,6];
	tmp$f = f0 & data1$exp[[etype]]$sd2[,2]>0 & data1$exp[[etype]]$sd2[,3]>0 & tmp$exon_diff==0 & data1$nucl_align$iden[,6]>=94 & data1$nucl_align$iden[,6]<97 & apply(fpkm[,7:18]>=thres, 1, sum)>=2
	tmp$sep = c(-1, seq(0,50,10), 100);
	tmp$x = data1$CNE5k_perc_diff[,6];
	tmp$y = data1$exp[[etype]]$pair_cor[,8];
	tmp$y1 = apply(data1$exp[[etype]]$sd2[,2:3],1,mean);
	tmp$y2 = data1$exp[[etype]]$pair_on_off_tissue_num[,8];
	tmp$m = length(tmp$sep)-1;
	tmp$xx = matrix(0, tmp$m, 6);
	tmp$x1 = NULL;
	for (i in 1:tmp$m) {
		tmp$f1 = tmp$f & tmp$x>tmp$sep[i] & tmp$x<=tmp$sep[i+1];
		tmp$xx[i,] = c( mean(tmp$x[tmp$f1]),  mean(tmp$y[tmp$f1]),   
				sd(tmp$y[tmp$f1])/sqrt(sum(tmp$f1)),    sum(tmp$f1),
				mean(tmp$y1[tmp$f1]), mean(tmp$y2[tmp$f1])
				);
		tmp$x1 = rbind(tmp$x1, cbind(i,tmp$y[tmp$f1], tmp$y1[tmp$f1], tmp$y2[tmp$f1]));
	}
	colnames(tmp$x1) = c('CNEDiff', 'Corr', 'SD', 'N_tissue');

	wilcox.test(tmp$x1[tmp$x1[,1]<=2,2], tmp$x1[tmp$x1[,1]>2,2], alternative='greater');
	wilcox.test(tmp$x1[tmp$x1[,1]<=2,3], tmp$x1[tmp$x1[,1]>2,3], alternative='less');

	jpeg_name  = paste0("output/plot/ZF_GF.",etype,".CNE5k_diff_cor_GF1_GF2.jpg");
	jpeg_name2 = paste0("output/plot/ZF_GF.",etype,".CNE5k_diff_cor_GF1_GF2.boxplot.jpg");
	xlab='CNE difference between GF gene pairs (%)';
	ylab='Average correlation';
	ylab2='Correlation';
	xtick = c(paste0('=',tmp$sep[2]), tmp$sep[3:(tmp$m+1)]);
# line plot with size
	jpeg(jpeg_name, res=res, width=6*res, height=4*res);
	plot(1:tmp$m, tmp$xx[,2], type='l', ylim=c(min(tmp$xx[,2]-0.05),max(tmp$xx[,2]+0.05)), axes=F, xlab=xlab, ylab=ylab);
#axis(1,at=c(1:tmp$m), lab=c('<80', paste0(tmp$idens1[2:tmp$m],'-',tmp$idens1[3:(tmp$m+1)])) );
	axis(1,at=c(1:tmp$m), lab=xtick );
	axis(2);
	points(1:tmp$m, tmp$xx[,2], cex=log2(tmp$xx[,4])/4, pch=21, bg='#0000FF80', col='#0000FF');
	dev.off();
# boxplot
	jpeg(jpeg_name2, res=res, width=6*res, height=4*res);
	colnames(tmp$x1) = c('CNEDiff', 'Corr', 'SD', 'N_tissue');
	boxplot( Corr ~ CNEDiff, tmp$x1, names=xtick, ylim=c(0.4,1), col='#30C080', xlab=xlab, ylab=ylab2);
	dev.off();
}
# }}}

# CNE5k difference by length percentage, cor(ZF,GF)
# {{{
for (etype in c('fpkm','tpm')) {
	fpkm = data1[[etype]];
	tmp=list();
	tmp$x1 = data1$exon_perc_diff[,2];
	tmp$x2 = data1$exon_perc_diff[,3];
	tmp$exon_diff = tmp$x = c(tmp$x1,tmp$x2);
	tmp$mean_lfpkm1 = rep(data1$exp[[etype]]$mean2[,1],2);
	tmp$mean_lfpkm2 = unlist(data1$exp[[etype]]$mean2[,2:3]);
	tmp$f = c( f0 & data1$exp[[etype]]$sd2[,1]>0 & data1$exp[[etype]]$sd2[,2]>0, f0 & data1$exp[[etype]]$sd2[,1]>0 & data1$exp[[etype]]$sd2[,3]>0) & tmp$exon_diff==0 & unlist(data1$nucl_align$iden[,2:3])>85 & unlist(data1$nucl_align$iden[,2:3])<97;
	tmp$sep = c(-1, seq(0,50,10), 100);
	tmp$x1 = data1$CNE5k_perc_diff[,2];
	tmp$x2 = data1$CNE5k_perc_diff[,3];
	tmp$x = c(tmp$x1,tmp$x2);
	tmp$y = c(data1$exp[[etype]]$pair_cor[,2], data1$exp[[etype]]$pair_cor[,3]);
	tmp$dist = c(data1$exp[[etype]]$pair_max_dist[,2], data1$exp[[etype]]$pair_max_dist[,3]);
	tmp$m = length(tmp$sep)-1;
	tmp$xx = matrix(0, tmp$m, 7);
	tmp$x1 = NULL;
	for (i in 1:tmp$m) {
		tmp$f1 = tmp$f & tmp$x>tmp$sep[i] & tmp$x<=tmp$sep[i+1];
		tmp$xx[i,] = c( mean(tmp$x[tmp$f1]),  mean(tmp$y[tmp$f1]),   
				sd(tmp$y[tmp$f1])/sqrt(sum(tmp$f1)),    sum(tmp$f1),
				mean_fpkm1=mean(tmp$mean_lfpkm1[tmp$f1]),
				mean_fpkm2=mean(tmp$mean_lfpkm2[tmp$f1]),
				mean(tmp$dist[tmp$f1])
				);
		tmp$x1 = rbind(tmp$x1, cbind(i,tmp$y[tmp$f1], tmp$mean_lfpkm1[tmp$f1], tmp$mean_lfpkm2[tmp$f1], tmp$dist[tmp$f1]));
	}
	colnames(tmp$x1) = c('CNE5k_Diff', 'Corr', 'mean_ZF_lfpkm', 'mean_GF_lfpkm', 'log2FC');

	jpeg_name  = paste0("output/plot/ZF_GF.",etype,".CNE5k_diff_cor_ZF_GF.jpg");
	jpeg_name2 = paste0("output/plot/ZF_GF.",etype,".CNE5k_diff_cor_ZF_GF.boxplot.jpg");
	xlab='CNE5k difference between ZF and GF (%)';
	ylab='Average correlation';
	ylab2='Correlation';
	xtick = c(paste0('=',tmp$sep[2]), tmp$sep[3:(tmp$m+1)]);
# line plot with size
	jpeg(jpeg_name, res=res, width=6*res, height=4*res);
	plot(1:tmp$m, tmp$xx[,2], type='l', ylim=c(min(tmp$xx[,2]-0.05),max(tmp$xx[,2]+0.05)), axes=F, xlab=xlab, ylab=ylab);
#axis(1,at=c(1:tmp$m), lab=c('<80', paste0(tmp$idens1[2:tmp$m],'-',tmp$idens1[3:(tmp$m+1)])) );
	axis(1,at=c(1:tmp$m), lab=xtick );
	axis(2);
	points(1:tmp$m, tmp$xx[,2], cex=log2(tmp$xx[,4])/4, pch=21, bg='#0000FF80', col='#0000FF');
	dev.off();
# boxplot
	jpeg(jpeg_name2, res=res, width=6*res, height=4*res);
	boxplot( Corr ~ CNE5k_Diff, tmp$x1, names=xtick, ylim=c(0.3,1), col='#30C080', xlab=xlab, ylab=ylab2, axes=F);
	axis(1,at=c(1:tmp$m), lab=xtick, xpd=T );
	axis(2);
	dev.off();
#
	jpeg_name  = paste0("output/plot/ZF_GF.",etype,".CNE5k_diff_max_dist_ZF_GF.jpg");
	jpeg_name2 = paste0("output/plot/ZF_GF.",etype,".CNE5k_diff_max_dist_ZF_GF.boxplot.jpg");
	xlab='CNE5k difference between ZF and GF (%)';
	ylab='Average correlation';
	ylab2='Correlation';
	xtick = c(paste0('=',tmp$sep[2]), tmp$sep[3:(tmp$m+1)]);
# line plot with size
	jpeg(jpeg_name, res=res, width=6*res, height=4*res);
	plot(1:tmp$m, tmp$xx[,7], type='l', ylim=c(min(tmp$xx[,7]-0.05),max(tmp$xx[,7]+0.05)), axes=F, xlab=xlab, ylab=ylab);
#axis(1,at=c(1:tmp$m), lab=c('<80', paste0(tmp$idens1[2:tmp$m],'-',tmp$idens1[3:(tmp$m+1)])) );
	axis(1,at=c(1:tmp$m), lab=xtick );
	axis(2);
	points(1:tmp$m, tmp$xx[,7], cex=log2(tmp$xx[,4])/4, pch=21, bg='#0000FF80', col='#0000FF');
	dev.off();
# boxplot
	jpeg(jpeg_name2, res=res, width=6*res, height=4*res);
	boxplot( log2FC ~ CNE5k_Diff, tmp$x1, names=xtick, ylim=c(0.3,1), col='#30C080', xlab=xlab, ylab=ylab2, axes=F);
	axis(1,at=c(1:tmp$m), lab=xtick, xpd=T );
	axis(2);
	dev.off();
}

# }}}

tmp$domain_diff = 1-data1$pair_interpor_jaccard;
# }}}

# output table: number of co-expressed tissues ~ number of triplets
# {{{
for (etype in c('fpkm','tpm')) {
	x = matrix(0, nts+1, 4);
	colnames(x) = c('ZF-GF1', 'ZF-GF2', 'ZF-GF_pair', 'GF1-GF2');
	rownames(x) = 0:6;
	xf = matrix(0, nts+1, 4);
	colnames(xf) = c('ZF-GF1', 'ZF-GF2', 'ZF-GF_pair', 'GF1-GF2');
	rownames(xf) = 0:6;
	for (i in 1:7) {
		x[i,1] = sum(data1$exp[[etype]]$pair_coexp_tissue_num[f0,2]==i-1);
		x[i,2] = sum(data1$exp[[etype]]$pair_coexp_tissue_num[f0,3]==i-1);
		x[i,3] = sum(data1$exp[[etype]]$pair_coexp_tissue_num[f0,5]==i-1);
		x[i,4] = sum(data1$exp[[etype]]$pair_coexp_tissue_num[f0,8]==i-1);
		xf[i,] = x[i,]/sum(f0);
	}
	tmp_out = xf;
	for (i in 1:7) { tmp_out[i,] = paste0(x[i,], ' (', sprintf('%.2f', xf[i,]*100), '%)'); }
	write.table(tmp_out, file=paste0('output/ZF_GF.',etype,'.triplet_count_coexp_tissue.txt'), sep="\t", quote=F);
	rm(tmp_out, x, xf);
}
# }}}

# plot correlation between ZF, GF
# {{{
tmp_cor_mean = c();
tmp_cor_sd = c();
tmp_cor_se = c();
# GF
f=data1$f & data1$max2[,1]>=thres;
x=x4=apply(data1$pair_cor[f,2:3], 1, max);
tmp_cor_mean["maximum of \nZF~GF1\nZF~GF2"] = mean(x);
tmp_cor_sd["maximum of \nZF~GF1\nZF~GF2"] = sd(x);
tmp_cor_se["maximum of \nZF~GF1\nZF~GF2"] = sd(x)/sqrt(length(x));
x=x5=apply(data1$pair_cor[f,2:3], 1, min);
tmp_cor_mean["minimum of \nZF~GF1\nZF~GF2"] = mean(x);
tmp_cor_sd["minimum of \nZF~GF1\nZF~GF2"] = sd(x);
tmp_cor_se["minimum of \nZF~GF1\nZF~GF2"] = sd(x)/sqrt(length(x));
x=x1=c(data1$pair_cor[f,2], data1$pair_cor[f,3]);
tmp_cor_mean['ZF~GF'] = mean(x);
tmp_cor_sd['ZF~GF'] = sd(x);
tmp_cor_se['ZF~GF'] = sd(x)/sqrt(length(x));
x=x2=data1$pair_cor[f,5];
tmp_cor_mean['ZF~GF_pair'] = mean(x);
tmp_cor_sd['ZF~GF_pair'] = sd(x);
tmp_cor_se['ZF~GF_pair'] = sd(x)/sqrt(length(x));
x=x3=data1$pair_cor[f,5*(2-1)+3];
tmp_cor_mean['GF1~GF2'] = mean(x);
tmp_cor_sd['GF1~GF2'] = sd(x);
tmp_cor_se['GF1~GF2'] = sd(x)/sqrt(length(x));
p1=wilcox.test(x1,x2)$p.value; # 1.986e-15
p2=wilcox.test(x2,x3)$p.value; # < 2.2e-16

res=300;
jpeg("output/plot/ZF_GF.cor_bar_se.jpg", res=res, w=res*7, h=res*6);
par(mar=c(4,5,5,0.5)+0.1);
par(cex=1.2, cex.lab=1.2, cex.axis=1.1);
barplot_sd1(tmp_cor_mean, tmp_cor_se, width=0.9, space=c(rep(0.111,5)), col=color[1:5], ylim=c(0.4,0.8), border=color[1:5], xpd=F, ylab="Average Correlation", names.arg=F);
cx=par('cxy')[1];
cy=par('cxy')[2];
text(c(1:5)*1-0.45, 0.4-cy, names(tmp_cor_mean), xpd=T, cex=0.9, adj=c(0.5,1));
axis(2);
y0=tmp_cor_mean[3]+tmp_cor_se[3]+cy*0.5;
y1=tmp_cor_mean[4]+tmp_cor_se[4]+cy*0.5;
y2=max(y0,y1)+cy;
lines(c(2.55,2.55), c(y0,y2));
lines(c(2.55,3.55), c(y2,y2));
lines(c(3.55,3.55), c(y1,y2));
text(3, y2+cy, sprintf("%.2e",p1));
y0=tmp_cor_mean[4]+tmp_cor_se[4]+cy*0.5;
y1=tmp_cor_mean[5]+tmp_cor_se[5]+cy*0.5;
y2=max(y0,y1)+cy;
y1=y0=max(y0,y1);
lines(c(3.55,3.55), c(y0,y2), xpd=T);
lines(c(3.55,4.55), c(y2,y2), xpd=T);
lines(c(4.55,4.55), c(y1,y2), xpd=T);
text(4, y2+cy, sprintf("%.2e",p2), xpd=T);

dev.off();
# }}}

# plot maximum lfpkm diffenence (across tissues) between ZF, GF
# {{{
tmp_dist_mean = c();
tmp_dist_sd = c();
tmp_dist_se = c();
f0 = data1$f & data1$max2[,1]>=thres
# GF
f=f0;
x=x4=apply(data1$pair_max_dist[f,2:3], 1, max);
tmp_dist_mean["maximum of \nZF~GF1\nZF~GF2"] = mean(x);
tmp_dist_sd["maximum of \nZF~GF1\nZF~GF2"] = sd(x);
tmp_dist_se["maximum of \nZF~GF1\nZF~GF2"] = sd(x)/sqrt(length(x));
x=x5=apply(data1$pair_max_dist[f,2:3], 1, min);
tmp_dist_mean["minimum of \nZF~GF1\nZF~GF2"] = mean(x);
tmp_dist_sd["minimum of \nZF~GF1\nZF~GF2"] = sd(x);
tmp_dist_se["minimum of \nZF~GF1\nZF~GF2"] = sd(x)/sqrt(length(x));
x=x1=c(data1$pair_max_dist[f,2], data1$pair_max_dist[f,3]);
tmp_dist_mean['ZF~GF'] = mean(x);
tmp_dist_sd['ZF~GF'] = sd(x);
tmp_dist_se['ZF~GF'] = sd(x)/sqrt(length(x));
x=x2=data1$pair_max_dist[f,5];
tmp_dist_mean['ZF~GF_pair'] = mean(x);
tmp_dist_sd['ZF~GF_pair'] = sd(x);
tmp_dist_se['ZF~GF_pair'] = sd(x)/sqrt(length(x));
x=x3=data1$pair_max_dist[f,5*(2-1)+3];
tmp_dist_mean['GF1~GF2'] = mean(x);
tmp_dist_sd['GF1~GF2'] = sd(x);
tmp_dist_se['GF1~GF2'] = sd(x)/sqrt(length(x));
p1=wilcox.test(x1,x2)$p.value; # 1.986e-15
p2=wilcox.test(x2,x3)$p.value; # < 2.2e-16

res=300;
jpeg("output/plot/ZF_GF.max_log2FC_bar_se.jpg", res=res, w=res*7, h=res*6);
par(mar=c(4,5,5,0.5)+0.1);
par(cex=1.2, cex.lab=1.2, cex.axis=1.1);
barplot_sd1(tmp_dist_mean, tmp_dist_se, width=0.9, space=c(rep(0.111,5)), col=color[1:5], ylim=c(min(tmp_dist_mean)-1, max(tmp_dist_mean+1)), border=color[1:5], xpd=F, ylab="Maximum log2FC across tissues", names.arg=F);
cx=par('cxy')[1];
cy=par('cxy')[2];
text(c(1:5)*1-0.45, min(tmp_dist_mean)-1-cy, names(tmp_dist_mean), xpd=T, cex=0.9, adj=c(0.5,1));
axis(2);
y0=tmp_dist_mean[3]+tmp_dist_se[3]+cy*0.5;
y1=tmp_dist_mean[4]+tmp_dist_se[4]+cy*0.5;
y2=max(y0,y1)+cy;
lines(c(2.55,2.55), c(y0,y2));
lines(c(2.55,3.55), c(y2,y2));
lines(c(3.55,3.55), c(y1,y2));
text(3, y2+cy, sprintf("%.2e",p1));
y0=tmp_dist_mean[4]+tmp_dist_se[4]+cy*0.5;
y1=tmp_dist_mean[5]+tmp_dist_se[5]+cy*0.5;
y2=max(y0,y1)+cy;
y1=y0=max(y0,y1);
lines(c(3.55,3.55), c(y0,y2), xpd=T);
lines(c(3.55,4.55), c(y2,y2), xpd=T);
lines(c(4.55,4.55), c(y1,y2), xpd=T);
text(4, y2+cy, sprintf("%.2e",p2), xpd=T);

dev.off();
# }}}

# draw subfunc/neofunc expression
# {{{
b = fate_classify(data1$fpkm, nts, 1,2,3, thres, 0.5, 2, fc=2, C1=0.6, C3=0.75, cor_p=0.1, tp=0.01)
res=300;
for (i in (1:nrow(data1$fpkm))[data1$exp$fpkm$fate$GF$on_off_subfunc & f0]) {
	type='';
	if (b$on_off_subfunc[i]) { type='C'; } else { type='D'; }
	gid = data1$gene_ids2[i,1];
	chr = genes$ZF[gid,'chr'];
	pos = floor((genes$ZF[gid,'begin'])/1000);
	gname = genes$ZF[gid,'name'];
	disp_name = paste0(chr,"_",sprintf('%06i', pos),'_',gid,'_',gname, '_', type);
	jpeg(paste0("output/plot/subfunc/",disp_name,".jpg"), res=res, w=res*4, h=res*4);
	y1 = as.vector(t(matrix(data1$exp$fpkm$log2[i,c(1:(3*nts))], nts, 3)));
	col = rep(c('gray','red','blue'), 6);
	w=0.9;
	space=rep(c(0.11*3,0,0),6);
	par(family="Arial");
	par(mar=c(2.5,3,4.5,2)+0.1);
	barplot(y1, width=w, col=col, space=space, border=col, main=disp_name);
	text((1:6)*3-1.5, -par("cxy")[2], tissues, xpd=T);
	dev.off();
}
for (i in (1:nrow(data1$fpkm))[data1$exp$fpkm$fate$GF$on_off_neofunc_WGD & f0]) {
	type='';
	if (b$on_off_neofunc[i]) { type='C'; } else { type='D'; }
	gid = data1$gene_ids2[i,1];
	chr = genes$ZF[gid,'chr'];
	pos = floor((genes$ZF[gid,'begin'])/1000);
	gname = genes$ZF[gid,'name'];
	disp_name = paste0(chr,"_",sprintf('%06i', pos),'_',gid,'_',gname, '_', type);
	jpeg(paste0("output/plot/neofunc/",disp_name,".jpg"), res=res, w=res*4, h=res*4);
	y1 = as.vector(t(matrix(data1$exp$fpkm$log2[i,c(1:(3*nts))], nts, 3)));
	col = rep(c('gray','red','blue'), 6);
	w=0.9;
	space=rep(c(0.11*3,0,0),6);
	par(family="Arial");
	par(mar=c(2.5,3,4.5,2)+0.1);
	barplot(y1, width=w, col=col, space=space, border=col, main=disp_name);
	text((1:6)*3-1.5, -par("cxy")[2], tissues, xpd=T);
	dev.off();
}
# }}}

# select neo-func gene for paper
# f2=data1$fate$GF$on_off_neofunc & data1$nucl_align$iden[2]>=87 & data1$nucl_align$iden[3]>=87 & apply(data1$exon_loss_per[,c('000','010','011')],1,sum)>=95

# draw chromosome region
f = genes$ZF[,'chr']=='1' & genes$ZF[,'begin']>=11300000 & genes$ZF[,'end']<=11800000;
genes1 = genes$ZF[f,];
gids1 = genes1[,'id'];
color1 = rep('grey', nrow(genes1));
names(color1) = gids1;
a1=data1$fate$GF$on_off_subfunc
names(a1) = data1$gene_ids2[,1];
color1[ a1[ gids1 ] ] = color[3];
a2=data1$fate$GF$on_off_neofunc
names(a2) = data1$gene_ids2[,1];
color1[ a2[ gids1 ] ] = color[4];
plot_chr_lfpkm(genes1, lfpkm[,1:18], 6, color1);


# example of sub-functionalization:
# gene: ENSDARG00000101585, tecra, data.i=3928
# gene: data.i=73, "si:ch211-213o11.11"
# {{{
# i select from subfunc_score>0.5 & is_subfunc in both CC and GF
res=300;
jpeg("output/plot/subfunc_si_ch211-213o11.11.jpg", res=res, w=res*8, h=res*4);
i=73;
y1 = as.vector(t(matrix(data$lfpkm[i,c(1:(3*nts))], nts, 3)));
col = rep(c('gray','red','blue'), 6);
w=0.9;
space=rep(c(0.11*3,0,0),6);
par(family="Arial");
par(mfrow=c(1,2), mar=c(5,3,1,1)+0.1);
barplot(y1, width=w, col=col, space=space, border=col, ylim=c(0,4));
text((1:6)*3-1.5, -par("cxy")[2], tissues, xpd=T);
y2 = as.vector(t(matrix(data$lfpkm[i,c(1:nts,(3*nts+1):(5*nts))], nts, 3)));
barplot(y2, width=w, col=col, space=space, border=col, ylim=c(0,4));
text((1:6)*3-1.5, -par("cxy")[2], tissues, xpd=T);
dev.off();
# }}}

# example of neo-functionalization:
# gene: ENSDARG00000078797, dennd3a, data.i=1844
# gene: ENSDARG00000045705, metg1, data.i=12
# gene:nrn1la,data.i=1711,gene_ids="ENSDARG00000069368","109065766","109109317","CA00029509","CA00029499"
# {{{
res=300;
jpeg("output/plot/neofunc_nrn1la.jpg", res=res, w=res*8, h=res*4);
i=1711
par(family="Arial");
par(mfrow=c(1,2), mar=c(5,3,1,1)+0.1);
y1 = as.vector(t(matrix(data$lfpkm[i,c(1:(3*nts))], nts, 3)));
col = rep(c('gray','red','blue'), 6);
w=0.9;
space=rep(c(0.111*3,0,0),6);
barplot(y1, width=w, col=col, space=space, border=col, ylim=c(0,4));
text((1:6)*3-1.5, -par("cxy")[2], tissues, xpd=T);
y2 = as.vector(t(matrix(data$lfpkm[i,c(1:nts,(3*nts+1):(5*nts))], nts, 3)));
barplot(y2, width=w, col=col, space=space, border=col, ylim=c(0,4));
text((1:6)*3-1.5, -par("cxy")[2], tissues, xpd=T);
dev.off();
# }}}

############################################################
# draw 2D-histogram between: ZF GF
############################################################
# {{{
# 1,1  dist(GF1,GF2) and cor(GF1,GF2);
# 1,2  cor(ZF,GF1) and cor(ZF,GF1+GF2);
# 2,1  cor(ZF,GF1+GF2) and cor(ZF,GF2);
# 2,2  cor(ZF,GF1) and cor(ZF,GF2)
res=300;
jpeg("output/plot/ZF_GF.cor_hist2d.jpg", res=res, w=res*10, h=res*8.5);
par(mfrow=c(2,2)); par(mar=c(5,5.7,1,5.5)+0.1);
for (k in 1:4) {
	f = data1$f;     tot = sum(f);
	if (k==1) {
		x = data1$pair_cor[f,8];   y = data1$pair_max_dist[f,8];
		y[y>4] = 4;
		xlab = "Correlation of GF1 and GF2";
		ylab = "Maximum log2 fold change \nbetween GF1 and GF2";
	} else if (k==2) {
		x = data1$pair_cor[f,2];   y = data1$pair_cor[f,5];
		xlab = "Correlation of ZF and GF1";
		ylab = "Correlation of ZF and GF pairs";
	} else if (k==3) {
		x = data1$pair_cor[f,5];   y = data1$pair_cor[f,3];
		xlab = "Correlation of ZF and GF pairs";
		ylab = "Correlation of ZF and GF2";
	} else if (k==4) {
		x = data1$pair_cor[f,2];   y = data1$pair_cor[f,3];
		xlab = "Correlation of ZF and GF1";
		ylab = "Correlation of ZF and GF2";
	}
	n = length(x);
	x_cut = 0.6;      y_cut = 0.6;
	if (k==1) { y_cut=1; }
	n11 = sum(x<0.6 & y<y_cut);    p11 = sprintf("%.2f", n11*100/n);
	n12 = sum(x<0.6 & y>=y_cut);   p12 = sprintf("%.2f", n12*100/n);
	n21 = sum(x>=0.6 & y<y_cut);   p21 = sprintf("%.2f", n21*100/n);
	n22 = sum(x>=0.6 & y>=y_cut);  p22 = sprintf("%.2f", n22*100/n);
	x = c(x,-1);  if (k!=1) { y = c(y,-1); } else { y=c(y,0); }
	a = hist2d(cbind(x,y), nbin=20, col=blueyellow255, xlab=xlab, ylab=ylab, same.scale=F, cex=1.3, cex.axis=1.3, cex.lab=1.6);
	cx=par('cxy')[1];   cy=par('cxy')[2];
#	image(c(1.05,1.15), -1+0.5*(0:256)/256, matrix(0:255,1,256), col=blueyellow255, xpd=T, add=T);
	if (k==1) {
		for (i in 1:255) { rect(1.1,max(y)/2*(i-1)/256,1.2,max(y)/2*i/256, col=blueyellow255[i], border=blueyellow255[i], xpd=T); }
		text(1.22, 0  , sprintf("%.1f%%",min(a$counts)*100/tot), xpd=T, adj=0, cex=1.4);
		text(1.22, max(y)/4, sprintf("%.1f%%",(min(a$counts)+max(a$counts))*100/2/tot), xpd=T, adj=0, cex=1.4);
		text(1.22, max(y)/2, sprintf("%.1f%%",max(a$counts)*100/tot), xpd=T, adj=0, cex=1.4);
		abline(h=y_cut, col='gray');   abline(v=0.6, col='gray');
		text(-0.2,(max(y)+y_cut)*0.5, paste0(n11, "\n(", p11, "%)"), xpd=T, cex=1.5, col='orange');
		text(-0.2,y_cut*0.6  , paste0(n12, "\n(", p12, "%)"), xpd=T, cex=1.5, col='orange');
		text(0.8 ,(max(y)+y_cut)*0.5, paste0(n21, "\n(", p21, "%)"), xpd=T, cex=1.5, col='orange');
		text(0.8 ,y_cut*0.6  , paste0(n22, "\n(", p22, "%)"), xpd=T, cex=1.5, col='orange');
		text(0.6+cx,min(y)+cy*0.6, 0.6, xpd=T, col='gray', cex=1.4);   text(-0.9,y_cut+cy*0.6, y_cut, xpd=T, col='gray', cex=1.4);
	} else {
		for (i in 1:255) { rect(1.1,-1+1.0*(i-1)/256,1.2,-1+1.0*i/256, col=blueyellow255[i], border=blueyellow255[i], xpd=T); }
		text(1.22, -1  , sprintf("%.1f%%",min(a$counts)*100/tot), xpd=T, adj=0, cex=1.4);
		text(1.22, -0.5, sprintf("%.1f%%",(min(a$counts)+max(a$counts))*100/2/tot), xpd=T, adj=0, cex=1.4);
		text(1.22, 0   , sprintf("%.1f%%",max(a$counts)*100/tot), xpd=T, adj=0, cex=1.4);
		abline(h=y_cut, col='gray');   abline(v=0.6, col='gray');
		text(-0.2,-0.2, paste0(n11, "\n(", p11, "%)"), xpd=T, cex=1.5, col='orange');
		text(-0.2,0.8, paste0(n12, "\n(", p12, "%)"), xpd=T, cex=1.5, col='orange');
		text(0.8,-0.2, paste0(n21, "\n(", p21, "%)"), xpd=T, cex=1.5, col='orange');
		text(1.25,0.8, paste0(n22, "\n(", p22, "%)"), xpd=T, cex=1.5, col='orange');
		text(0.6+cx,-0.9, 0.6, xpd=T, col='gray', cex=1.4);   text(-0.9,0.6+cy*0.6, 0.6, xpd=T, col='gray', cex=1.4);
	}
}
dev.off();
# }}}

# ZF--GF1--GF2 Triplet: iden ~ number of category
# {{{
tmp_fate_names=on_off_fate_names;
tmp_disp_fate_names=on_off_disp_fate_names;
for (etype in c('fpkm','tpm')) {
	fpkm = data1[[etype]];
	data1$exp[[etype]]$fate_iden_count = list();
	data1$exp[[etype]]$fate_iden = list();
	xx=75:95;
	lxx=length(xx);
	fate=data1$exp[[etype]]$fate;
	for (ft in tmp_fate_names) {
		x=matrix(c(xx,rep(0,lxx*2)), lxx, 3);
		colnames(x) = c("iden", "ZF_GF", "GF1_GF2");
		if (regexpr("coexpress|conserved|corr",ft)>0 || regexpr("dosage",ft)>0 || regexpr('subfunc', ft)>0 ) {
			a1 = c(data1$nucl_align$iden[,2],data1$nucl_align$iden[,3]);
			f2 = a1>=0 & rep(f0,2);
			b1 = rep(fate$GF[[ft]],2);
		} else {
			a1 = apply(data1$nucl_align$iden[,2:3],1,min);
			ft1 = paste0(ft, '1');    ft2 = paste0(ft, '2');
			f=fate$GF[[ft1]] & !fate$GF[[ft2]];
			a1[f] = data1$nucl_align$iden[f,2];
			f=fate$GF[[ft2]] & !fate$GF[[ft1]];
			a1[f] = data1$nucl_align$iden[f,3];
			f2 = a1>=0 & f0;
			b1 = fate$GF[[ft]];
		}

		a = data1$nucl_align$iden[,6];

		y=list();
		y$ZF_GF = a1[fate$GF[[ft]] & f2];
		y$GF1_GF2 = a[fate$GF[[ft]] & a>=0 & f0];
		for (i in 1:lxx) {
			f = a1<=xx[i] & a1>=0 & f0;
			x[i,2] = sum(b1 & f)/sum(b1[f2]);

			f = a<=xx[i] & a>=0 & f0;
			x[i,3] = sum(fate$GF[[ft]] & f)/sum(fate$GF[[ft]][a>=0 & f0]);
		}
		data1$exp[[etype]]$fate_iden_count[[ft]]=x;
		data1$exp[[etype]]$fate_iden[[ft]]=y;
	}
}
for (etype in c('fpkm','tpm')) {
	res=300; cex=1.3;
	jpeg(paste0("output/plot/ZF_GF.",etype,".fate_ZF_GF_iden_count.jpg"), res=res, w=res*8, h=res*4);
	par('mar'=c(4,5.5,1.5,15)+0.1);
	xlab="ZF-GF Nucleic Identity (%)";
	ylab="Cumulative sums of triplets (%)";
	plot(cbind(c(75,95),c(0,100)), col='white', axes=F, xlab=xlab, ylab=ylab, cex=cex, cex.lab=cex, cex.axis=cex, xpd=T);
	axis(1, cex.axis=1.3, line=0.6); axis(2, cex.axis=1.3, line=0.6);
	i=0;
	cx=par('cxy')[1];   cy=par('cxy')[2];
	y = par('usr')[4]-cy;
#for (ft in fate_names[c(2,4,5:9)])
	for (ft in tmp_fate_names[1:4])
	{
		x=data1$exp[[etype]]$fate_iden_count[[ft]];
		i=i+1;
		lines(x[,1],x[,2]*100, col=color[i], lwd=1.5);
#	lines(x[,1],log10(x[,2]*100), col=color[i], lwd=1.5);
		points(95+cx,y, pch=16, xpd=T, col=color[i]);
		text(95+cx*2,y, paste0(tmp_disp_fate_names[ft]), xpd=T, col=color[i], adj=c(0,0.5), cex=1.6);
		y=y-cy*1.8;
	}
	dev.off();

	res=300; cex=1.3;
	jpeg(paste0("output/plot/ZF_GF.",etype,".fate_GF_pair_iden_count.jpg"), res=res, w=res*8, h=res*4);
	par('mar'=c(4,5.5,1.5,16)+0.1);
	xlab="GF1-GF2 Nucleic Identity (%)";
	ylab="Cumulative sums of triplets (%)";
	plot(cbind(c(75,95),c(0,100)), col='white', axes=F, xlab=xlab, ylab=ylab, cex=cex, cex.lab=cex, cex.axis=cex, xpd=T);
	axis(1, cex.axis=1.3, line=0.6); axis(2, cex.axis=1.3, line=0.6);
	i=0;
	cx=par('cxy')[1];    cy=par('cxy')[2];
	y = par('usr')[4]-cy;
#for (ft in fate_names[c(2,4,5:9)]) 
	for (ft in tmp_fate_names[1:4])
	{
		x=data1$exp[[etype]]$fate_iden_count[[ft]];
		i=i+1;
		lines(x[,1],x[,3]*100, col=color[i], lwd=1.5);
#	lines(x[,1],log10(x[,3]*100), col=color[i], lwd=1.5);
		points(95+cx,y,pch=20, xpd=T, col=color[i]);
		text(95+cx*2,y, paste0(tmp_disp_fate_names[ft]), xpd=T, col=color[i], adj=c(0,0.5), cex=1.6);
		y=y-cy*1.8;
	}
	dev.off();
}
# }}}

# ZF--GF1--GF2 Triplet: exon_dff ~ number of category
# {{{
tmp_fate_names=on_off_fate_names;
tmp_disp_fate_names=on_off_disp_fate_names;
for (etype in c('fpkm','tpm')) {
	fpkm = data1[[etype]];
	data1$exp[[etype]]$fate_exon_diff_count = list();
	data1$exp[[etype]]$fate_exon_diff = list();
	xx=0:100;
	lxx=length(xx);
	fate=data1$exp[[etype]]$fate;
	for (ft in tmp_fate_names) {
		x=matrix(c(xx,rep(0,lxx*2)), lxx, 3);
		colnames(x) = c("exon_diff", "ZF_GF", "GF1_GF2");
		if (regexpr("double_coexpress|conserve|corr|dosage",ft)>0 || regexpr('subfunc', ft)>0 ) {
			a1 = c(data1$exon_perc_diff[,2],data1$exon_perc_diff[,3]);
			f2 = rep(f0,2);
			b1 = rep(fate$GF[[ft]],2);
		} else {
			a1 = apply(data1$exon_perc_diff[,2:3],1,max);
			ft1 = paste0(ft, '1');    ft2 = paste0(ft, '2');
			f=fate$GF[[ft1]] & !fate$GF[[ft2]];
			a1[f] = data1$exon_perc_diff[f,2];
			f=fate$GF[[ft2]] & !fate$GF[[ft1]];
			a1[f] = data1$exon_perc_diff[f,3];
			f2 = f0;
			b1 = fate$GF[[ft]];
		}
		y=list();
		y$ZF_GF = a1[b1];
		a = data1$exon_perc_diff[,6];
		y$GF1_GF2 = a[ fate$GF[[ft]] ];
		for (i in 1:lxx) {
			f = a1<=xx[i] & f0;
			x[i,2] = sum(b1 & f)/sum(b1[f2]);

			f = a<=xx[i] & a>=0 & f0;
			x[i,3] = sum(fate$GF[[ft]] & f)/sum(fate$GF[[ft]][f0]);
		}
		data1$exp[[etype]]$fate_exon_diff_count[[ft]]=x;
		data1$exp[[etype]]$fate_exon_diff[[ft]]=y;
	}
}
for (etype in c('fpkm','tpm')) {
	res=300;  cex=1.3;
	jpeg(paste0("output/plot/ZF_GF.",etype,".fate_ZF_GF_exon_diff_count.jpg"), res=res, w=res*8, h=res*4);
	par('mar'=c(4,5.5,1.5,15)+0.1);
	xlab="ZF-GF exon gain/loss (%)";
	ylab="Cumulative sums of triplets (%)";
	plot(cbind(c(0,100),c(0,100)), col='white', axes=F, xlab=xlab, ylab=ylab, cex=cex, cex.lab=cex, cex.axis=cex, xpd=T);
	axis(1, cex.axis=1.3, line=0.6); axis(2, cex.axis=1.3, line=0.6);
	i=0;
	cx=par('cxy')[1];   cy=par('cxy')[2];
	y = par('usr')[4]-cy;
#for (ft in fate_names[c(2,4,5:9)])
	for (ft in tmp_fate_names[1:4])
	{
		x=data1$exp[[etype]]$fate_exon_diff_count[[ft]];
		i=i+1;
		lines(x[,1],x[,2]*100, col=color[i], lwd=1.5);
		points(100+cx, y, pch=16, xpd=T, col=color[i]);
		text(100+cx*2, y, paste0(tmp_disp_fate_names[ft]), xpd=T, col=color[i], adj=c(0,0.5), cex=1.5);
		y=y-cy*1.3;
	}
	dev.off();

	res=300;  cex=1.2
	jpeg(paste0("output/plot/ZF_GF.",etype,".fate_GF_pair_exon_diff_count.jpg"), res=res, w=res*9, h=res*5);
	par('mar'=c(4,5.5,1.5,15)+0.1);
	xlab="GF1-GF2 exon gain/loss (%)";
	ylab="Cumulative sums of triplets (%)";
	plot(cbind(c(0,100),c(0,100)), col='white', axes=F, xlab=xlab, ylab=ylab, cex=cex, cex.lab=cex, cex.axis=cex, xpd=T);
	axis(1, cex.axis=1.3, line=0.6); axis(2, cex.axis=1.3, line=0.6);
	i=0;
	cx=par('cxy')[1];
	cy=par('cxy')[2];
	y = par('usr')[4]-cy;
	for (ft in tmp_fate_names[1:4]) {
		x=data1$exp[[etype]]$fate_exon_diff_count[[ft]];
		i=i+1;
		lines(x[,1],x[,3]*100, col=color[i], lwd=1.5);
		points(100+cx, y, pch=20, xpd=T, col=color[i]);
		text(100+cx*2, y, paste0(tmp_disp_fate_names[ft]), xpd=T, col=color[i], adj=c(0,0.5), cex=1.5);
		y=y-cy*1.3;
	}
	dev.off();
}
# }}}

# ZF--GF1--GF2 Triplet: CNE5k_dff ~ number of category
# {{{
tmp_fate_names=on_off_fate_names;
tmp_disp_fate_names=on_off_disp_fate_names;
for (etype in c('fpkm','tpm')) {
	fpkm = data1[[etype]];
	fate = data1$exp[[etype]]$fate;
	data1$exp[[etype]]$fate_CNE5k_diff_count = list();
	data1$exp[[etype]]$fate_CNE5k_diff = list();
	xx=0:100;
	lxx=length(xx);
	for (ft in tmp_fate_names) {
		x=matrix(c(xx,rep(0,lxx*2)), lxx, 3);
		colnames(x) = c("CNE5k_diff", "ZF_GF", "GF1_GF2");
		if (regexpr("double_coexp|conserve|corr|dosage",ft)>0 || regexpr('subfunc', ft)>0 ) {
			a1 = c(data1$CNE5k_perc_diff[,2],data1$CNE5k_perc_diff[,3]);
			f2 = rep(f0,2);
			b1 = rep(fate$GF[[ft]],2);
		} else {
			a1 = apply(data1$CNE5k_perc_diff[,2:3],1,max);
			ft1 = paste0(ft, '1');    ft2 = paste0(ft, '2');
			f=fate$GF[[ft1]] & !fate$GF[[ft2]];
			a1[f] = data1$CNE5k_perc_diff[f,2];
			f=fate$GF[[ft2]] & !fate$GF[[ft1]];
			a1[f] = data1$CNE5k_perc_diff[f,3];
			f2 = f0;
			b1 = fate$GF[[ft]];
		}
		a = data1$CNE5k_perc_diff[,6];
		y=list();
		y$ZF_GF = a1[b1 & f2];
		y$GF1_GF2 = a[fate$GF[[ft]] & f0];
		for (i in 1:lxx) {
			f = a1<=xx[i] & a1>=0 & f2;
			x[i,2] = sum(b1 & f)/sum(b1[f2]);

			f = a<=xx[i] & a>=0 & f0;
			x[i,3] = sum(fate$GF[[ft]] & f)/sum(fate$GF[[ft]][f0]);
		}
		data1$exp[[etype]]$fate_CNE5k_diff_count[[ft]]=x;
		data1$exp[[etype]]$fate_CNE5k_diff[[ft]]=y;
	}
}
res=600; cex=1.3;
for (etype in c('fpkm','tpm')) {
	tiff(paste0("output/plot/ZF_GF.",etype,".fate_ZF_GF_CNE5k_diff_count.tif"), res=res, w=res*6, h=res*3.5, compress='lzw');
	par('mar'=c(4,4.5,1.5,10)+0.1);
	xlab="ZF-GF CNE gain/loss (%)";
	ylab="Cumulative sums of triplets (%)";
	plot(cbind(c(0,100),c(0,100)), col='white', axes=F, xlab=xlab, ylab=ylab, cex=cex, cex.lab=cex, cex.axis=cex, xpd=T);
	axis(1, cex.axis=1.3, line=0.6); axis(2, cex.axis=1.3, line=0.6);
	i=0;
	cx=par('cxy')[1];   cy=par('cxy')[2];
	y = par('usr')[4]-cy;
#for (ft in fate_names[c(2,4,5:9)])
	for (ft in tmp_fate_names[1:4])
	{
		x=data1$exp[[etype]]$fate_CNE5k_diff_count[[ft]];
		i=i+1;
		lines(x[,1],x[,2]*100, col=color[i], lwd=1.5);
#	lines(x[,1],log10(x[,2]*100), col=color[i], lwd=1.5);
		points(100+cx, y, pch=16, xpd=T, col=color[i]);
		text(100+cx*2, y, paste0(tmp_disp_fate_names[ft]), xpd=T, col=color[i], adj=c(0,0.5), cex=1.8);
		y=y-cy*1.3;
	}
	dev.off();

	tiff(paste0("output/plot/ZF_GF.",etype,".fate_GF_pair_CNE5k_diff_count.tif"), res=res, w=res*6, h=res*3.5, compress='lzw');
	par('mar'=c(4,4.5,1.5,10)+0.1);
	xlab="GF1-GF2 CNE gain/loss (%)";
	ylab="Cumulative sums of triplets (%)";
	plot(cbind(c(0,100),c(0,100)), col='white', axes=F, xlab=xlab, ylab=ylab, cex=cex, cex.lab=cex, cex.axis=cex, xpd=T);
	axis(1, cex.axis=1.3, line=0.6); axis(2, cex.axis=1.3, line=0.6);
	i=0;
	cx=par('cxy')[1];
	cy=par('cxy')[2];
	y = par('usr')[4]-cy;
	for (ft in tmp_fate_names[1:4]) {
		x=data1$exp[[etype]]$fate_CNE5k_diff_count[[ft]];
		i=i+1;
		lines(x[,1],x[,3]*100, col=color[i], lwd=1.5);
		points(100+cx, y, pch=20, xpd=T, col=color[i]);
		text(100+cx*2, y, paste0(tmp_disp_fate_names[ft]), xpd=T, col=color[i], adj=c(0,0.5), cex=1.8);
		y=y-cy*1.3;
	}
	dev.off();
}
# }}}

# ZF--GF1--GF2 Triplet: maximal lfpkm ~ number of category
# {{{
tmp_fate_names=on_off_fate_names;
tmp_disp_fate_names=on_off_disp_fate_names;
for (etype in c('fpkm','tpm')) {
	data1$exp[[etype]]$fate_exp_count = list();
	data1$exp[[etype]]$fate_exp = list();
	fate = data1$exp[[etype]]$fate$GF;
	xx = seq(1,10,0.5);
	lxx = length(xx);
	for (ft in tmp_fate_names) {
		y=list();
		x=matrix(c(xx,rep(0,lxx*2)), lxx, 3);
		colnames(x) = c("lfpkm", "ZF", "GF");
		a1 = data1$exp[[etype]]$max2[,1];
		a2 = data1$exp[[etype]]$max2[,5];
		y[[1]] = a1[fate[[ft]] & f0];
		y[[2]] = a2[fate[[ft]] & f0];
#	y[[2]] = a = apply(data1$mean2[,2:3],1,sum);
		for (i in 1:lxx) {
			f = a1<=xx[i] & f0;
			x[i,2] = sum(fate[[ft]] & f)/sum(fate[[ft]] & f0);
			f = a2<=xx[i] & f0;
			x[i,3] = sum(fate[[ft]] & f)/sum(fate[[ft]] & f0);
		}
		data1$exp[[etype]]$fate_exp_count[[ft]]=x;
		names(y) = c("ZF", "GF");
		data1$exp[[etype]]$fate_exp[[ft]]=y;
	}

	res=600;
	for (sp in c('ZF','GF') ) {
	tiff(paste0("output/plot/ZF_GF.",etype,".fate_",sp,"_max_exp_count.tif"), res=res, w=res*6, h=res*3.5, compress='lzw');
	par('mar'=c(4,4.5,1.5,10)+0.1);
	xlab=paste0("log2(",etype,"+1)");
	ylab="Cumulative sums of triplets (%)";
	plot(cbind(c(0,10),c(0,100)), col='white', axes=F, xlab=xlab, ylab=ylab, cex=cex, cex.lab=cex, cex.axis=cex, xpd=T);
	axis(1, cex.axis=1.3, line=0.6); axis(2, cex.axis=1.3, line=0.6);
	i=0;
	cx=par('cxy')[1];
	cy=par('cxy')[2];
	y = par('usr')[4]-cy;
	for (ft in tmp_fate_names[1:4]) {
		x=data1$exp[[etype]]$fate_exp_count[[ft]];
		i=i+1;
		if (sp=='ZF') { lines(x[,1], x[,2]*100, col=color[i], lwd=1.5); }
		else if (sp=='GF') { lines(x[,1], x[,3]*100, col=color[i], lwd=1.5); }
		points(par('usr')[2]+1.5*cx,y,pch=16, xpd=T, col=color[i]);
		text(par('usr')[2]+cx*2.5, y, paste0(tmp_disp_fate_names[ft]), xpd=T, col=color[i], adj=c(0,0.5), cex=1.8);
		y=y-cy*1.5;
	}
	dev.off();
	}
}
# }}}

# ZF--GF1--GF2 Triplet: expressed tissues ~ number of category
# {{{
tmp_fate_names=on_off_fate_names;
tmp_disp_fate_names=on_off_disp_fate_names;
for (etype in c('fpkm','tpm')) {
	data1$exp[[etype]]$fate_exp_sm_count = list();
	data1$exp[[etype]]$fate_exp_sm = list();
	fate = data1$exp[[etype]]$fate$GF;
	xx = seq(1,nts,1);
	lxx = length(xx);
	for (ft in tmp_fate_names) {
		y=list();
		x=matrix(c(xx,rep(0,lxx*2)), lxx, 3);
		colnames(x) = c("tissue_num", "ZF", "GF");
		a1 = data1$exp[[etype]]$exp_sm_count[,1];
		y[[1]] = a1[fate[[ft]] & f0];
		a2 = data1$exp[[etype]]$exp_sm_count[,5];
		y[[2]] = a2[fate[[ft]] & f0];
#	y[[2]] = a = apply(data1$exp_sm_count[,2:3],1,max);
		for (i in 1:lxx) {
			f = a1<=xx[i] & f0;
			x[i,2] = sum(fate[[ft]] & f)/sum(fate[[ft]] & f0);

			f = a2<=xx[i] & f0;
			x[i,3] = sum(fate[[ft]] & f)/sum(fate[[ft]] & f0);
		}
		data1$exp[[etype]]$fate_exp_sm_count[[ft]]=x;
		names(y) = c("ZF", "GF");
		data1$exp[[etype]]$fate_exp_sm[[ft]]=y;
	}

	res=600;
	for (sp in c("ZF","GF") ) {
		tiff(paste0("output/plot/ZF_GF.",etype,".fate_exp_tissue_num_",sp,".tif"), res=res, w=res*6, h=res*3.5, compress='lzw');
		par('mar'=c(4,4.5,1.5,10)+0.1);
		xlab="No. expressed tissues";
		ylab="Cumulative sums of triplets (%)";
		plot(cbind(c(1,6),c(0,100)), col='white', axes=F, xlab=xlab, ylab=ylab, cex=cex, cex.lab=cex, cex.axis=cex, xpd=T);
		axis(1, cex.axis=1.3, line=0.6); axis(2, cex.axis=1.3, line=0.6);
		i=0;
		cx=par('cxy')[1];
		cy=par('cxy')[2];
		y = par('usr')[4]-cy;
		for (ft in tmp_fate_names[1:4]) {
			x=data1$exp[[etype]]$fate_exp_sm_count[[ft]];
			i=i+1;
			if (sp=="ZF") { lines(x[,1], x[,2]*100, col=color[i], lwd=1.5); }
			else if (sp=="GF") { lines(x[,1], x[,3]*100, col=color[i], lwd=1.5); }
			points(par('usr')[2]+cx,y,pch=16, xpd=T, col=color[i]);
			text(par('usr')[2]+2*cx, y, paste0(tmp_disp_fate_names[ft]), xpd=T, col=color[i], adj=c(0,0.5), cex=1.5);
			y=y-cy*1.5;
		}
		dev.off();
	}
}
# }}}

# Pairwise Test identity, exon diff, CNE diff on ZF-GF1 vs. ZF-GF2 in each fate
# {{{
for (etype in c('fpkm', 'tpm')) {
	tmp = list();
	fate = data1$exp[[etype]]$fate$GF;
	res = NULL;
	res_names=NULL;
	for (ft in on_off_fate_names[3:4]) {
		tmp[[ft]]=list();
		if (regexpr("double_coexp|conserve|corr|dosage",ft)>0 || regexpr('subfunc', ft)>0 ) {
			tmp$f = f0 & fate[[ft]];
			tmp[[ft]] = list( iden1 = apply(data1$nucl_align$iden[tmp$f,2:3],1,min),
					iden2 = apply(data1$nucl_align$iden[tmp$f,2:3],1,max),
					exon_diff1 = apply(data1$exon_perc_diff[tmp$f,2:3],1,min),
					exon_diff2 = apply(data1$exon_perc_diff[tmp$f,2:3],1,max),
					CNE5k_diff1 = apply(data1$CNE5k_perc_diff[tmp$f,2:3],1,min),
					CNE5k_diff2 = apply(data1$CNE5k_perc_diff[tmp$f,2:3],1,max)
					);
		} else {
			ft1 = paste0(ft, '1');    ft2 = paste0(ft, '2');
			f1 = fate[[ft1]]&f0;   f2=fate[[ft2]]&f0;
			a1 = data1$nucl_align$iden[,2];   a2 = data1$nucl_align$iden[,3];;
			tmp[[ft]]$iden1 = c(a2[f1],a1[f2]);
			tmp[[ft]]$iden2 = c(a1[f1],a2[f2]);
			a1 = data1$exon_perc_diff[,2];   a2 = data1$exon_perc_diff[,3];;
			tmp[[ft]]$exon_diff1 = c(a2[f1],a1[f2]);
			tmp[[ft]]$exon_diff2 = c(a1[f1],a2[f2]);
			a1 = data1$CNE5k_perc_diff[,2];   a2 = data1$CNE5k_perc_diff[,3];;
			tmp[[ft]]$CNE5k_diff1 = c(a2[f1],a1[f2]);
			tmp[[ft]]$CNE5k_diff2 = c(a1[f1],a2[f2]);
		}
		mean = NULL;
		for (name in names(tmp[[ft]])) { mean=c(mean, mean(tmp[[ft]][[name]])); }
		names(mean) = names(tmp[[ft]]);
		for (name in c('iden', 'exon_diff', 'CNE5k_diff') ) {
			name1 = paste0(name, '1'); name2 = paste0(name, '2');
			a1 = tmp[[ft]][[name1]];   a2 = tmp[[ft]][[name2]];
			p1 = wilcox.test(a1,a2, pair=T)$p.value;
			p2 = t.test(a1,a2, pair=T)$p.value;
			res = rbind(res, c(mean1=sprintf("%.2f",mean[name1]), sprintf("%.2f",mean2=mean[name2]), diff=sprintf("%.2f",mean[name1]-mean[name2]), sprintf("%.4f",p1), sprintf("%.4f",p2)) );
			res_names = c(res_names, paste0(ft, '_', name));
		}
	}
	rownames(res) = res_names;
	colnames(res) = c('mean1', 'mean2', 'difference', 'P (Wilcoxon rank test)', 'P (t-test)');
	write.table(res, paste0("output/ZF_GF.",etype,".fate_pair_test.txt"), sep="\t", quote=F);
}
# }}}

# Test ZF-GF (or GF1,GF2) identity, exon diff, CNE diff between different fate groups
for (etype in c('fpkm', 'tpm')) {
	res = NULL;
	for (feat in c('iden', 'exon_diff', 'CNE5k_diff') ) {
		for (pair in c('ZF_GF', 'GF1_GF2')) {
			for (i1 in 1:3) {
				ft1 = on_off_fate_names[i1];
				for (i2 in (i1+1):4) {
					ft2 = on_off_fate_names[i2];
					name1 = paste0('fate_',feat);
					name2 = paste0('fate_',feat);
					a1 = data1$exp[[etype]][[name1]][[ft1]][[pair]];
					a2 = data1$exp[[etype]][[name2]][[ft2]][[pair]];
					p = wilcox.test(a1,a2)$p.value;
					res = rbind(res, c(feat, pair, ft1, ft2, sprintf('%.2f',mean(a1)), sprintf('%.2f',mean(a2)), sprintf('%.4f',p)));
				}
			}
		}
	}
	colnames(res) = c('Feature', 'species pair', 'Fate1', 'Fate2', 'Mean1', 'Mean2', 'P value');
	write.table(res, paste0("output/ZF_GF.",etype,".fate_fate_pair_test.txt"), row.names=F, sep="\t", quote=F);
}

############################################################
# ZF--GF1--GF2 Triplet: GO over-representation
############################################################
# {{{
for (etype in c('fpkm', 'tpm')) {
	tmp_fate_names=on_off_fate_names;
	data1$exp[[etype]]$fate_enrich_func=list();
	fate=data1$exp$fpkm$fate;
	for (domain in c('GO2_BP', 'GO2_MF', 'GO2_CC', 'KEGG')) {
		print(domain);
		data1$exp[[etype]]$fate_enrich_func[[domain]]=list();
		sp='GF';
		data1$exp[[etype]]$fate_enrich_func[[domain]][[sp]]=list();
		for (ft in tmp_fate_names) {
			print(ft);
			data1$exp[[etype]]$fate_enrich_func[[domain]][[sp]][[ft]] = go_test(gs[[domain]]$ZF, rownames(data1$fpkm)[fate[[sp]][[ft]] & f0], rownames(data1$fpkm)[f0]);
		}
	}
	res=300;
	p = 0.05;
	fdr=0.25;
	for (domain in c('GO2_BP', 'GO2_MF', 'GO2_CC', 'KEGG')) {
		for (ft in tmp_fate_names[1:4]) {
			fname=paste0("output/plot/go_test/ZF_GF.",etype,".fate_",domain,"_",ft,".jpg");
			jpeg(fname, res=res, w=res*12, h=res*8);
			for (sp in c("GF")) {
				par(mar=c(25,5,1,12)+0.1);
				cy = par('cxy')[2];
				a = data1$exp[[etype]]$fate_enrich_func[[domain]][[sp]][[ft]];
				a = a[gs_id_name[[domain]][,'annot_id'],];
				a = a[ order(a[,'p']), ];
				lv = gs_id_name[[domain]][ rownames(a), 'level' ];
				if (domain=='KEGG') {
					a_up = a[ a[,'ratio']>1, ];
					a_dn = a[ a[,'ratio']<1, ];
				} else {
					a_up = a[ a[,'ratio']>1 & lv>=2 & lv<=6, ];
					a_dn = a[ a[,'ratio']<1 & lv>=2 & lv<=6, ];
				}
				if (sum(a_up[,'p']<p & a_up[,'fdr']<fdr)<5) { up_term = rownames(a_up)[1:(min(length(a_up[,1]),5))]; } else { up_term = rownames(a_up)[a_up[,'p']<p & a_up[,'fdr']<fdr]; }
				if (sum(a_dn[,'p']<p & a_dn[,'fdr']<fdr)<5) { dn_term = rownames(a_dn)[1:(min(length(a_dn[,1]),5))]; } else { dn_term = rownames(a_dn)[a_dn[,'p']<p & a_dn[,'fdr']<fdr]; }
				if (length(up_term)>25) { up_term=up_term[1:25]; }
				if (length(dn_term)>25) { dn_term=dn_term[1:25]; }
				bar_h = c(-log10( a[up_term, 'p'] ),  rev(-log10( a[dn_term, 'p'] )) );
				if (length(bar_h)==0) { next; }
				names(bar_h) = gs_id_name[[domain]][ c(up_term,rev(dn_term)), 'annot_name' ];
				col=c(rep('red', length(up_term)), rep('blue', length(dn_term)));
				barplot(bar_h, col=col, border=col, width=0.8, space=0.25, names.arg=F, ylab="log10 P-value", cex.lab=2, cex.axis=1.5);
				text(1:length(bar_h)-0.5, -cy*2, names(bar_h), srt=-45, adj=0, cex=1.5, xpd=T)
			}
			dev.off();
		}
	}
}
# print 
for (ft in tmp_fate_names) {
	print(ft);
	p=0.05; fdr=0.05;
	for (domain in c('GO2_BP', 'GO2_MF', 'GO2_CC', 'KEGG')) {
		print(domain);
		a = data1$exp[[etype]]$fate_enrich_func[[domain]][["GF"]][[ft]];
		a = a[ a[,'ratio']>1, ];
		a = a[ order(a[,'p']), ];
		f = a[,'p']<p & a[,'fdr']<fdr;
		if (sum(f)<5) { a = a[1:5,]; } else { a = a[ a[,'p']<p & a[,'fdr']<fdr, ]; }
		name=gs_id_name[[domain]][rownames(a),c('annot_name','level')];
		a = cbind(name,a);
		print(a);
	}
}
# low and high identity GF1--GF2
f0 = data1$f & data1$max2[,1]>=thres;
fh = data1$nucl_align$iden[,6]>=quantile(data1$nucl_align$iden[,6], 0.75);
fl = data1$nucl_align$iden[,6]<=quantile(data1$nucl_align$iden[,6], 0.25);
data1$high_iden_go_test = list();
data1$low_iden_go_test = list();
for (domain in c('GO2_BP', 'GO2_MF', 'KEGG')) {
	data1$high_iden_go_test[[domain]] = list();
	data1$high_iden_go_test[[domain]]$GF1_GF2 = go_test(gs[[domain]]$ZF, rownames(data1$fpkm)[fh & f0], rownames(data1$fpkm)[f0]);
	data1$low_iden_go_test[[domain]] = list();
	data1$low_iden_go_test[[domain]]$GF1_GF2 = go_test(gs[[domain]]$ZF, rownames(data1$fpkm)[fl & f0], rownames(data1$fpkm)[f0]);
}
for (domain in c('GO2_BP', 'GO2_MF', 'GO2_CC', 'KEGG')) {
	print('low identity');
	a = data1$low_iden_go_test[[domain]]$GF1_GF2;
	a = a[ a[,'ratio']>1, ];
	a = a[ order(a[,'p']), ];
	f = a[,'p']<0.05 & a[,'fdr']<0.1;
	if (sum(f)<5) { a = a[1:5,]; } else { a = a[ a[,'p']<0.05 & a[,'fdr']<0.1, ]; }
	name=gs_id_name[[domain]][rownames(a),c('annot_name','level')];
	a = cbind(name,a);
	print(a);

	print("\nhigh identity");
	a = data1$high_iden_go_test[[domain]]$GF1_GF2;
	a = a[ a[,'ratio']>1, ];
	a = a[ order(a[,'p']), ];
	f = a[,'p']<0.05 & a[,'fdr']<0.1;
	if (sum(f)<5) { a = a[1:5,]; } else { a = a[ a[,'p']<0.05 & a[,'fdr']<0.1, ]; }
	name=gs_id_name[[domain]][rownames(a),c('annot_name','level')];
	a = cbind(name,a);
	print(a);

}
# }}}

# test expression ~ nucleotide identity
tmp=list();
tmp$f1= data1$exp$fpkm$mean2[,3]-data1$exp$fpkm$mean2[,2]>=1 ;
tmp$f2= data1$exp$fpkm$mean2[,2]-data1$exp$fpkm$mean2[,3]>=1 ;
tmp$iden_low_exp  = c(data1$nucl_align$iden[tmp$f1,2], data1$nucl_align$iden[tmp$f2,3]);
tmp$iden_high_exp = c(data1$nucl_align$iden[tmp$f1,3], data1$nucl_align$iden[tmp$f2,2]);
wilcox.test(tmp$iden_low_exp, tmp$iden_high_exp, pair=T); # 0.03, less


############################################################
# ZF--GF1--GF2 Triplet: by chromosome
############################################################
corr_p = c();
for (chr in 1:25) {
	chr1 = paste0('LG', chr);
	chr2 = paste0('LG', 25+chr);
	a1 = unique(bgp$GF[bgp$GF[,1]==chr1,c(1,18:19)]);
	rownames(a1) = a1[,2];
	a2 = unique(bgp$GF[bgp$GF[,1]==chr2,c(1,18:19)]);
	rownames(a2) = a2[,2];
	f1 = data1$gene_ids2[,2] %in% a1[,2];
	f2 = data1$gene_ids2[,3] %in% a2[,2];
	ff1 = f1 & f2;
	f1 = data1$gene_ids2[,2] %in% a2[,2];
	f2 = data1$gene_ids2[,3] %in% a1[,2];
	ff2 = f1 & f2;
	pair_cor = rbind(data1$pair_cor[ff1,2:3], data1$pair_cor[ff2,3:2]);
	print(chr);
	res=wilcox.test(pair_cor[,1], pair_cor[,2], paired=T);
	corr_p = c(corr_p, res$p.value);
	print(res);

	lfpkm_pair = rbind(lfpkm[ff1,(nts*1+1):(nts*3)], lfpkm[ff2,c((nts*2+1):(nts*3), (nts*1+1):(nts*2))] );
	lfpkm_mean = cbind(apply(lfpkm_pair[,1:nts],1,mean), apply(lfpkm_pair[,(nts+1):(nts*2)],1,mean)); 
	lfpkm_pair_tp = pair_t_p(lfpkm_pair[,1:nts], lfpkm_pair[,(nts+1):(nts*2)]);
	lfpkm_pair_fdr = p.adjust(lfpkm_pair_tp, method='fdr');
}

# iden ~ number of category
# {{{
data$fate_iden_count = list();
data$fate_iden - list();
xx=70:100;
lxx=length(xx);
for (ft in fate_names) {
	x=matrix(c(xx,rep(0,lxx*4)), lxx, 5);
	colnames(x) = c("iden", "ZF_CC","ZF_GF", "CC1_CC2", "GF1_GF2");
	if (ft=="double_conserved" || ft=="dosage_balance" || ft=="subfunc") {
		a1 = apply(data$nucl_iden[,2:3],1,min);
		a2 = apply(data$nucl_iden[,4:5],1,min);
	} else {
		a1 = apply(data$nucl_iden[,2:3],1,min);
		ft1 = paste0(ft, '1');    ft2 = paste0(ft, '2');
		f=data$fate$CC[[ft1]] & !data$fate$CC[[ft2]];
		a1[f] = data$nucl_iden[f,2];
		f=data$fate$CC[[ft2]] & !data$fate$CC[[ft1]];
		a1[f] = data$nucl_iden[f,3];

		a2 = apply(data$nucl_iden[,4:5],1,min);
		f=data$fate$GF[[ft1]] & !data$fate$GF[[ft2]];
		a2[f] = data$nucl_iden[f,4];
		f=data$fate$GF[[ft2]] & !data$fate$GF[[ft1]];
		a2[f] = data$nucl_iden[f,5];
	}
	y=list();
for (i in 1:lxx) {
	f = a1<=xx[i] & a1>=0;   f2= a1>=0;
	y$ZF_CC = a1[data$fate$CC[[ft]]];
	x[i,2] = sum(data$fate$CC[[ft]] & f)/sum(data$fate$CC[[ft]][f2]);
	f = a2<=xx[i] & a2>=0;   f2= a2>=0;
	y$ZF_GF = a2[data$fate$GF[[ft]]];
	x[i,3] = sum(data$fate$GF[[ft]] & f)/sum(data$fate$GF[[ft]][f2]);

	a = data$nucl_iden[,8];
	f = a<=xx[i] & a>=0;   f2= a>=0;
	y$CC1_CC2 = a[data$fate$CC[[ft]]];
	x[i,4] = sum(data$fate$CC[[ft]] & f)/sum(data$fate$CC[[ft]][f2]);

	a = data$nucl_iden[,20];
	f = a<=xx[i] & a>=0;   f2= a>=0;
	y$GF1_GF2 = a[data$fate$GF[[ft]]];
	x[i,5] = sum(data$fate$CC[[ft]] & f)/sum(data$fate$CC[[ft]][f2]);
}
	data$fate_iden_count[[ft]]=x;
	data$fate_iden[[ft]]=y;
}

res=300;
jpeg("output/plot/ZF_iden_count.jpg", res=res, w=res*8, h=res*6);
par('mar'=c(5,6,1,2)+0.1);
plot(cbind(c(70,100),c(0,1)), col='white', axes=F, xlab="Nucleic Identity (%)", ylab="Number of ohnolog clusters\n with identity less than x", cex=1.5, cex.lab=1.5, cex.axis=1.5);
axis(1, cex.axis=1.5);
axis(2, cex.axis=1.5);
i=0;
cx=par('cxy')[1];   cy=par('cxy')[2];
y = 0.4;
for (ft in fate_names) {
	x=data$fate_iden_count[[ft]];
	i=i+1;
#	lines(x[,c(1,2)], col=color[i], lty='dashed', lwd=1.5);
	lines(x[,c(1,3)], col=color[i], lwd=1.5);
#	rect(70,0.8,74+cx*7*0.9,0.95, col='#E0E0E0', border='lightgray');
#	lines(c(71,73),c(0.9,0.9));
#	text(74,0.9, "ZF vs. GF", adj=0, xpd=T);
#	lines(c(71,73),c(0.85,0.85), lty='dashed');
#	text(74,0.85, "ZF vs. CC", adj=0, xpd=T);
	points(88+cx,y,pch=16, xpd=T, col=color[i]);
	text(88+cx*2, y, paste0(fate_names1[ft]), xpd=T, col=color[i], adj=c(0,0.5), cex=1.8);
	y=y-cy*1.5;
}
dev.off();

res=300;
jpeg("output/plot/pair_iden_count.jpg", res=res, w=res*8, h=res*6);
par('mar'=c(5,5,1,10)+0.1);
plot(cbind(c(70,100),c(0,1)), col='white', axes=F, xlab="Nucleic Identity (%)", ylab="Number of ohnolog clusters\n with identity less than x");
axis(1);
axis(2);
i=0;
cx=par('cxy')[1];
cy=par('cxy')[2];
y = 1;
for (ft in fate_names) {
	x=data$fate_iden_count[[ft]];
	i=i+1;
	lines(x[,c(1,4)], col=color[i], lty='dashed');
	lines(x[,c(1,5)], col=color[i]);
	points(100+cx,y,pch=20, xpd=T, col=color[i]);
	text(100+cx*2, y, paste0(ft), xpd=T, col=color[i], adj=0);
	y=y-cy;
}
dev.off();
# }}}

# iden X coverage ~ number of category
# {{{
data$fate_iden_cov_count = list();
xx=20:100;
lxx=length(xx);
for (ft in fate_names) {
	x=matrix(c(xx,rep(0,lxx*4)), lxx, 5);
	colnames(x) = c("iden x coverage", "ZF_CC","ZF_GF", "CC1_CC2", "GF1_GF2");
	if (ft=="double_conserved" || ft=="dosage_balance" || ft=="subfunc") {
		a1 = apply(data$nucl_iden[,2:3]*data$nucl_cov[,2:3]/100,1,min);
		a2 = apply(data$nucl_iden[,4:5]*data$nucl_cov[,4:5]/100,1,min);
	} else {
		a1 = apply(data$nucl_iden[,2:3]*data$nucl_cov[,2:3]/100,1,min);
		ft1 = paste0(ft, '1');
		ft2 = paste0(ft, '2');
		f=data$fate$CC[[ft1]] & !data$fate$CC[[ft2]];
		a1[f] = data$nucl_iden[f,2]*data$nucl_cov[f,2];
		f=data$fate$CC[[ft2]] & !data$fate$CC[[ft1]];
		a1[f] = data$nucl_iden[f,3]*data$nucl_cov[f,3];

		a2 = apply(data$nucl_iden[,4:5]*data$nucl_cov[,4:5]/100,1,min);
		f=data$fate$GF[[ft1]] & !data$fate$GF[[ft2]];
		a2[f] = data$nucl_iden[f,4];
		f=data$fate$GF[[ft2]] & !data$fate$GF[[ft1]];
		a2[f] = data$nucl_iden[f,5];
	}
for (i in 1:lxx) {
	f = a1<=xx[i] & a1>0;
	f2= a1>0;
	x[i,2] = sum(data$fate$CC[[ft]] & f)/sum(data$fate$CC[[ft]][f2]);
	f = a2<=xx[i] & a2>0;
	f2= a2>0;
	x[i,3] = sum(data$fate$GF[[ft]] & f)/sum(data$fate$GF[[ft]][f2]);

	a = data$nucl_iden[,8]*data$nucl_cov[,8];
	f = a<=xx[i] & a>0;
	f2= a>0;
	x[i,4] = sum(data$fate$CC[[ft]] & f)/sum(data$fate$CC[[ft]][f2]);
	a = data$nucl_iden[,20]*data$nucl_cov[,20];
	f = a<=xx[i] & a>0;
	f2= a>0;
	x[i,5] = sum(data$fate$CC[[ft]] & f)/sum(data$fate$CC[[ft]][f2]);
}
	data$fate_iden_cov_count[[ft]]=x;
}

res=300;
jpeg("output/plot/ZF_iden_cov_count.jpg", res=res, w=res*8, h=res*6);
par('mar'=c(5,6,1,2)+0.1);
plot(cbind(c(xx[1],xx[lxx]),c(0,1)), col='white', axes=F, xlab="Nucleic Identity x Coverage (%)", ylab="Number of ohnolog clusters\n with identity less than x", cex=1.5, cex.lab=1.5, cex.axis=1.5);
axis(1, cex.axis=1.5);
axis(2, cex.axis=1.5);
i=0;
cx=par('cxy')[1];
cy=par('cxy')[2];
y = 0.4;
for (ft in fate_names) {
	x=data$fate_iden_cov_count[[ft]];
	i=i+1;
	lines(x[,c(1,2)], col=color[i], lty='dashed', lwd=1.5);
	lines(x[,c(1,3)], col=color[i], lwd=1.5);
	rect(xx[1],0.8,74+cx*7*0.9,0.95, col='#E0E0E0', border='lightgray');
	lines(xx[1]+c(1,3),c(0.9,0.9));
	text(xx[1]+4,0.9, "ZF vs. GF", adj=0, xpd=T);
	lines(xx[1]+c(1,3),c(0.85,0.85), lty='dashed');
	text(xx[1]+4,0.85, "ZF vs. CC", adj=0, xpd=T);
	points(xx[lxx]-cx*10,y,pch=16, xpd=T, col=color[i]);
	text(xx[lxx]-cx*8, y, paste0(fate_names1[ft]), xpd=T, col=color[i], adj=c(0,0.5), cex=1.8);
	y=y-cy*1.5;
}
dev.off();

res=300;
jpeg("output/plot/pair_iden_count.jpg", res=res, w=res*8, h=res*6);
par('mar'=c(5,5,1,10)+0.1);
plot(cbind(c(70,100),c(0,1)), col='white', axes=F, xlab="Nucleic Identity (%)", ylab="Number of ohnolog clusters\n with identity less than x");
axis(1);
axis(2);
i=0;
cx=par('cxy')[1];
cy=par('cxy')[2];
y = 1;
for (ft in fate_names) {
	x=data$fate_iden_count[[ft]];
	i=i+1;
	lines(x[,c(1,4)], col=color[i], lty='dashed');
	lines(x[,c(1,5)], col=color[i]);
	points(100+cx,y,pch=20, xpd=T, col=color[i]);
	text(100+cx*2, y, paste0(ft), xpd=T, col=color[i], adj=0);
	y=y-cy;
}
dev.off();
# }}}

# lfpkm ~ number of category
# {{{
data$fate_lfpkm_count = list();
data$fate_lfpkm = list();
xx = seq(0,12,0.1);
lxx = length(xx);
for (ft in fate_names) {
	y=list();
	x=matrix(c(xx,rep(0,lxx*4)), lxx, 5);
	colnames(x) = c("mean_lfpkm", "ZF_CC", "ZF_GF", "CC", "GF");
for (i in 1:lxx) {
	y[[1]] = a = data$mean2[,1];
	f = a<=xx[i];
	x[i,2] = sum(data$fate$CC[[ft]] & f)/sum(data$fate$CC[[ft]]);

	y[[2]] = a = data$mean2[,1];
	f = a<=xx[i];
	x[i,3] = sum(data$fate$GF[[ft]] & f)/sum(data$fate$GF[[ft]]);

	y[[3]] = a = apply(data$mean2[,2:3],1,max);
	f = a<=xx[i];
	x[i,4] = sum(data$fate$CC[[ft]] & f)/sum(data$fate$CC[[ft]]);

	y[[4]] = a = apply(data$mean2[,4:5],1,max);
	f = a<=xx[i];
	x[i,5] = sum(data$fate$GF[[ft]] & f)/sum(data$fate$GF[[ft]]);
}
	data$fate_lfpkm_count[[ft]]=x;
	names(y) = c("ZF_CC", "ZF_GF", "CC", "GF");
	data$fate_lfpkm[[ft]]=y;
}

res=300;
for (sp in c("CC","GF") ) {
jpeg(paste0("output/plot/ZF_",sp,"_lfpkm_count.jpg"), res=res, w=res*8, h=res*6);
par('mar'=c(5,6,1,2)+0.1);
plot(cbind(c(xx[1],xx[lxx]),c(0,1)), col='white', axes=F, xlab="log2(fpkm+1)", ylab="Number of ohnolog clusters\n with zebrafish log2(fpkm+1) less than x", cex=1.5, cex.lab=1.5, cex.axis=1.5);
axis(1, cex.axis=1.5);
axis(2, cex.axis=1.5);
i=0;
cx=par('cxy')[1];
cy=par('cxy')[2];
y = 0.4;
for (ft in fate_names) {
	x=data$fate_lfpkm_count[[ft]];
	i=i+1;
	if (sp=="CC") { lines(x[,c(1,2)], col=color[i], lwd=1.5); }
	else if (sp=="GF") { lines(x[,c(1,3)], col=color[i], lwd=1.5); }
#	lines(x[,c(1,3)], col=color[i], lty='dashed', lwd=1.5);
#	lines(x[,c(1,4)], col=color[i], lty='dotted', lwd=1.5);
#	rect(70,0.8,74+cx*7*0.9,0.95, col='#E0E0E0', border='lightgray');
#	lines(c(71,73),c(0.9,0.9));
#	text(74,0.9, "ZF vs. GF", adj=0, xpd=T);
#	lines(c(71,73),c(0.85,0.85), lty='dashed');
#	text(74,0.85, "ZF vs. CC", adj=0, xpd=T);
	points(6+cx,y,pch=16, xpd=T, col=color[i]);
	text(6+cx*2, y, paste0(fate_names1[ft]), xpd=T, col=color[i], adj=c(0,0.5), cex=1.8);
	y=y-cy*1.5;
}
dev.off();
}
# }}}

# expressed tissues ~ number of category
# {{{
data$fate_exp_sm_count = list();
data$fate_exp_sm = list();
xx = seq(1,nts,1);
lxx = length(xx);
for (ft in fate_names) {
	y=list();
	x=matrix(c(xx,rep(0,lxx*4)), lxx, 5);
	colnames(x) = c("tissue_num", "ZF_CC", "ZF_GF", "CC", "GF");
for (i in 1:lxx) {
	y[[1]] = a = data$exp_sm_count[,1];
	f = a<=xx[i];
	x[i,2] = sum(data$fate$CC[[ft]] & f)/sum(data$fate$CC[[ft]]);

	y[[2]] = a = data$exp_sm_count[,1];
	f = a<=xx[i];
	x[i,3] = sum(data$fate$GF[[ft]] & f)/sum(data$fate$GF[[ft]]);

	y[[3]] = a = apply(data$exp_sm_count[,2:3],1,max);
	f = a<=xx[i];
	x[i,4] = sum(data$fate$CC[[ft]] & f)/sum(data$fate$CC[[ft]]);

	y[[4]] = a = apply(data$exp_sm_count[,4:5],1,max);
	f = a<=xx[i];
	x[i,5] = sum(data$fate$GF[[ft]] & f)/sum(data$fate$GF[[ft]]);
}
	data$fate_exp_sm_count[[ft]]=x;
	names(y) = c("ZF_CC", "ZF_GF", "CC", "GF");
	data$fate_exp_sm[[ft]]=y;
}

res=300;
for (sp in c("CC","GF") ) {
jpeg(paste0("output/plot/ZF_",sp,"_exp_tissue.jpg"), res=res, w=res*8, h=res*6);
par('mar'=c(5,6,1,2)+0.1);
plot(cbind(c(xx[1],xx[lxx]),c(0,1)), col='white', axes=F, xlab="No. zebrafish expressed tissues", ylab="Number of ohnolog clusters\n with No. zebrafish expressed tissues less than x", cex=1.5, cex.lab=1.5, cex.axis=1.5);
axis(1, cex.axis=1.5);
axis(2, cex.axis=1.5);
i=0;
cx=par('cxy')[1];
cy=par('cxy')[2];
y = 1.0;
for (ft in fate_names) {
	x=data$fate_exp_sm_count[[ft]];
	i=i+1;
	if (sp=="CC") { lines(x[,c(1,2)], col=color[i], lwd=1.5); }
	else if (sp=="GF") { lines(x[,c(1,3)], col=color[i], lwd=1.5); }
#	lines(x[,c(1,3)], col=color[i], lty='dashed', lwd=1.5);
#	lines(x[,c(1,4)], col=color[i], lty='dotted', lwd=1.5);
#	rect(70,0.8,74+cx*7*0.9,0.95, col='#E0E0E0', border='lightgray');
#	lines(c(71,73),c(0.9,0.9));
#	text(74,0.9, "ZF vs. GF", adj=0, xpd=T);
#	lines(c(71,73),c(0.85,0.85), lty='dashed');
#	text(74,0.85, "ZF vs. CC", adj=0, xpd=T);
	points(1,y,pch=16, xpd=T, col=color[i]);
	text(1+cx, y, paste0(fate_names1[ft]), xpd=T, col=color[i], adj=c(0,0.5), cex=1.8);
	y=y-cy*1.5;
}
dev.off();
}
# }}}


############################################################
# test iden, bit
############################################################
# {{{
fate_names=c('double_conserved', 'dosage_balance','subfunc','neofunc','nonfunc')
for (k in 1:4) {
	p =matrix(1.0, 2, length(fate_names));
	rownames(p) = c("CC","GF");
	colnames(p) = fate_names;
	x=list();
	y=list();
	for (sp in c("CC","GF")) {
		if (k==1 || k==2) {
			idx=2:3;
			if (sp=="GF") {idx=4:5; }
		} else if (k==3 || k==4) {
			idx=(2-1)*5+3;
			if (sp=="GF") {idx=(4-1)*5+5; }
		}
		x[[sp]]=list();
		y[[sp]]=list();
		if (k==1) {
			a = apply(data$nucl_iden[,idx], 1, min);
		} else if (k==2) {
			a = apply(data$nucl_bit[,idx], 1, min);
		} else if (k==3) {
			a = data$nucl_iden[,idx];
		} else if (k==4) {
			a = data$nucl_bit[,idx];
		}
		x0= a[ data$fate[[sp]][['double_conserved']] & a>0 ];
		for (ft in c('double_conserved', 'dosage_balance','subfunc','neofunc','nonfunc')) {
			x1= a[ data$fate[[sp]][[ft]] & a>0 ];
#		y1= a[ !data$fate[[sp]][[ft]] & a>0 ];
			if (length(x1)>=5) {
				data$fate_ZF_iden_p[sp,ft] = wilcox.test(x0,x1)$p.value;
			}
			x[[sp]][[ft]] = x1;
			y[[sp]][[ft]] = y1;
		}
	}
	if (k==1) { data$fate_ZF_iden_p = p; }
	else if (k==2) { data$fate_ZF_bit_p = p; }
	else if (k==3) { data$fate_pair_iden_p = p; }
	else if (k==4) { data$fate_pair_bit_p = p; }
	res=300;
	if (k==1) { jpeg("output/plot/ZF_iden_wilcox.jpg", res=res, w=res*9, h=res*5); }
	else if (k==2) { jpeg("output/plot/ZF_bit_wilcox.jpg", res=res, w=res*9, h=res*5); }
	else if (k==3) { jpeg("output/plot/pair_iden_wilcox.jpg", res=res, w=res*9, h=res*5); }
	else if (k==4) { jpeg("output/plot/pair_bit_wilcox.jpg", res=res, w=res*9, h=res*5); }
	par(mfrow=c(1,2), mar=c(4,3,6,1)+0.1);
	for (sp in c('CC','GF')) {
		col = c("red", "red", "blue", "blue", "gray");
		boxplot(x[[sp]], axes=F, col=col, border="darkgray");
		axis(2);
		y0 = par('usr')[3];
		y1 = par('usr')[4];
		dy = 0.1*(y1-y0)/par('pin')[2];
		cy = par('cxy')[2];
		lines(c(1,1), c(y1,y1+dy), xpd=T);
		lines(c(1,3), c(y1+dy,y1+dy), xpd=T);
		lines(c(3,3), c(y1,y1+dy), xpd=T);
		text(2, y1+dy+cy/2, sprintf("%.2e",data$fate_ZF_iden_p[sp,3]), xpd=T);
		y1 = y1+dy+cy;
		lines(c(1,1), c(y1,y1+dy), xpd=T);
		lines(c(1,4), c(y1+dy,y1+dy), xpd=T);
		lines(c(4,4), c(y1,y1+dy), xpd=T);
		text(2, y1+dy+cy/2, sprintf("%.2e",data$fate_ZF_iden_p[sp,4]), xpd=T);
		y1 = y1+dy+cy;
		lines(c(1,1), c(y1,y1+dy), xpd=T);
		lines(c(1,5), c(y1+dy,y1+dy), xpd=T);
		lines(c(5,5), c(y1,y1+dy), xpd=T);
		text(2, y1+dy+cy/2, sprintf("%.2e",data$fate_ZF_iden_p[sp,5]), xpd=T);
		abline(h=y0, lwd=2);
		abline(v=par('usr')[1], lwd=2);
		text(1:5, rep(c(y0-cy/2,y0-cy-cy/2),5)[1:5], c("conserved", "dosage balance", "sub-func", 'neo-func', 'non-func'), xpd=T, col=col);
	}
	dev.off();
}


#CC.a1=as.vector( (data.n.iden*data.n.cov/10000)[data.fate$CC.double_conserved,2:3] );
#CC.a2=as.vector( (data.n.iden*data.n.cov/10000)[data.fate$CC.dosage_balance,2:3] );
#CC.a12=as.vector( (data.n.iden*data.n.cov/10000)[data.fate$CC.double_conserved,8] );
#CC.a0=as.vector( (data.n.iden*data.n.cov/10000)[data.fate$CC.double_conserved|data.fate$CC.dosage_balance,2:3] );
#CC.a012=as.vector( (data.n.iden*data.n.cov/10000)[data.fate$CC.double_conserved|data.fate$CC.dosage_balance,8] );
#GF.a1=as.vector( (data.n.iden*data.n.cov/10000)[data.fate$GF.double_conserved,4:5] );
#GF.a2=as.vector( (data.n.iden*data.n.cov/10000)[data.fate$GF.dosage_balance,4:5] );
#GF.a12=as.vector( (data.n.iden*data.n.cov/10000)[data.fate$GF.dosage_balance,5*(4-1)+5] );
#GF.a0=as.vector( (data.n.iden*data.n.cov/10000)[data.fate$GF.double_conserved|data.fate$GF.dosage_balance,4:5] );
#GF.a012=as.vector( (data.n.iden*data.n.cov/10000)[data.fate$GF.double_conserved|data.fate$GF.dosage_balance,5*(4-1)+5] );
## dosage balance ~~ conserved expressed
#wilcox.test(CC.a1[CC.a1>0],CC.a2[CC.a2>0]); # p=0.4673
#wilcox.test(GF.a1[GF.a1>0],GF.a2[GF.a2>0]); # p=0.041
## non-func ~~ dosage balance or conserved expressed
#CC.b= c( (data.n.iden*data.n.cov/10000)[data.fate$CC.nonfunc1,2],(data.n.iden*data.n.cov/10000)[data.fate$CC.nonfunc2,3] );
#CC.b2= c( (data.n.iden*data.n.cov/10000)[data.fate$CC.nonfunc1,3],(data.n.iden*data.n.cov/10000)[data.fate$CC.nonfunc2,2] );
#CC.b12=(data.n.iden*data.n.cov/10000)[data.fate$CC.nonfunc,8];
#wilcox.test(CC.a0[CC.a0>0],CC.b[CC.b>0], alternative='greater'); # p<2.2e-16
#wilcox.test(CC.b2[CC.b2>0],CC.b[CC.b>0], alternative='greater'); # p<2.2e-16
#wilcox.test(CC.a012[CC.a012>0], CC.b12[CC.b12>0],alternative='greater'); # p=5.015e-14
#GF.b= c( (data.n.iden*data.n.cov/10000)[data.fate$GF.nonfunc1,4],(data.n.iden*data.n.cov/10000)[data.fate$GF.nonfunc2,5] );
#GF.b2= c( (data.n.iden*data.n.cov/10000)[data.fate$GF.nonfunc1,5],(data.n.iden*data.n.cov/10000)[data.fate$GF.nonfunc2,4] );
#GF.b12=(data.n.iden*data.n.cov/10000)[data.fate$GF.nonfunc,5*(4-1)+5];
#wilcox.test(GF.a0[GF.a0>0],GF.b[GF.b>0], alternative='greater'); # p=5.287e-14
#wilcox.test(GF.b2[GF.b2>0],GF.b[GF.b>0], alternative='greater'); # p<2.2e-16
#wilcox.test(GF.a012[GF.a012>0], GF.b12[GF.b12>0],alternative='greater'); # p<2.2e-16
## sub-function ~~ 
#CC.b1= as.vector( (data.n.iden*data.n.cov/10000)[data.fate$CC.subfunc,2]);
#CC.b2= as.vector( (data.n.iden*data.n.cov/10000)[data.fate$CC.subfunc,3]);
#CC.b= as.vector( (data.n.iden*data.n.cov/10000)[data.fate$CC.subfunc,2:3]);
#CC.b12= as.vector( (data.n.iden*data.n.cov/10000)[data.fate$CC.subfunc,8]);
#wilcox.test(CC.a0[CC.a0>0],CC.b[CC.b>0], alternative='greater'); # p=0.011
#wilcox.test(CC.b1[CC.b1>0],CC.b2[CC.b2>0]);
#wilcox.test(CC.a012[CC.a012>0],CC.b12[CC.b12>0], alternative='greater'); # p=0.002
#GF.b1= as.vector( (data.n.iden*data.n.cov/10000)[data.fate$GF.subfunc,4]);
#GF.b2= as.vector( (data.n.iden*data.n.cov/10000)[data.fate$GF.subfunc,5]);
#GF.b= as.vector( (data.n.iden*data.n.cov/10000)[data.fate$GF.subfunc,4:5]);
#GF.b12= as.vector( (data.n.iden*data.n.cov/10000)[data.fate$GF.subfunc,20]);
#wilcox.test(GF.a0[GF.a0>0],GF.b[GF.b>0], alternative='greater'); # p=2.509e-05
#wilcox.test(GF.b1[GF.b1>0],GF.b2[GF.b2>0]);
#wilcox.test(GF.a012[GF.a012>0],GF.b12[GF.b12>0], alternative='greater'); # p=0.004
## neo-function ~~
#CC.b1= as.vector( (data.n.iden*data.n.cov/10000)[data.fate$CC.neofunc1,2]);
#CC.b2= as.vector( (data.n.iden*data.n.cov/10000)[data.fate$CC.neofunc2,3]);
#CC.b= c(CC.b1,CC.b2);
#CC.b12= as.vector( (data.n.iden*data.n.cov/10000)[data.fate$CC.neofunc,8]);
#wilcox.test(CC.a0[CC.a0>0],CC.b[CC.b>0], alternative='greater'); # p=4.696e-05
#wilcox.test(CC.a012[CC.a012>0],CC.b12[CC.b12>0], alternative='greater'); # p=0.001
#GF.b1= as.vector( (data.n.iden*data.n.cov/10000)[data.fate$GF.neofunc1,4]);
#GF.b2= as.vector( (data.n.iden*data.n.cov/10000)[data.fate$GF.neofunc2,5]);
#GF.b= c(GF.b1,GF.b2);
#GF.b12= as.vector( (data.n.iden*data.n.cov/10000)[data.fate$GF.neofunc,20]);
#wilcox.test(GF.a0[GF.a0>0],GF.b[GF.b>0], alternative='greater'); # p=6.385e-05
#wilcox.test(GF.a012[GF.a012>0],GF.b12[GF.b12>0], alternative='greater'); # p=0.002
# }}}

############################################################
# GO over-representation
############################################################
# {{{
tmp_fates = fate_names[c(2,4,7:9)];
data$fate_enrich_func=list();
for (domain in c('GO2_BP', 'GO2_MF', 'KEGG')) {
	data$fate_enrich_func[[domain]]=list();
	for (sp in c("CC", "GF")) {
		data$fate_enrich_func[[domain]][[sp]]=list();
		for (ft in tmp_fate_names) {
			data$fate_enrich_func[[domain]][[sp]][[ft]] = go_test(gs[[domain]], rownames(data$lfpkm)[data$fate[[sp]][[ft]]], rownames(data$lfpkm));
		}
	}
}
res=300;
for (domain in c('GO2_BP', 'GO2_MF')) {
	for (ft in tmp_fate_names) {
		fname=paste0("output/plot/",domain,"_",ft,".jpg");
		jpeg(fname, res=res, w=res*12, h=res*15);
		par(mfrow=c(2,1));
	for (sp in c("CC", "GF")) {
		par(mar=c(18,4,1,12)+0.1);
		cy = par('cxy')[2];
		a = data$fate_enrich_func[[domain]][[sp]][[ft]];
		a = a[gs_id_name[[domain]][,'annot_id'],];
		a = a[ order(a[,'p']), ];
		lv = gs_id_name[[domain]][ rownames(a), 'level' ];
		a_up = a[ a[,'ratio']>1 & lv>=2 & lv<=6, ];
		a_dn = a[ a[,'ratio']<1 & lv>=2 & lv<=6, ];
		if (sum(a_up[,'p']<0.1 & a_up[,'fdr']<0.25)<5) { up_term = rownames(a_up)[1:5]; } else { up_term = rownames(a_up)[a_up[,'p']<0.1 & a_up[,'fdr']<0.25]; }
		if (sum(a_dn[,'p']<0.1 & a_dn[,'fdr']<0.25)<5) { dn_term = rownames(a_dn)[1:5]; } else { dn_term = rownames(a_dn)[a_dn[,'p']<0.1 & a_dn[,'fdr']<0.25]; }
		if (length(up_term)>20) { up_term=up_term[1:20]; }
		if (length(dn_term)>20) { dn_term=dn_term[1:20]; }
		bar_h = c(-log10( a[up_term, 'p'] ),  rev(-log10( a[dn_term, 'p'] )) );
		if (length(bar_h)==0) { next; }
		names(bar_h) = gs_id_name[[domain]][ c(up_term,rev(dn_term)), 'annot_name' ];
		col=c(rep('red', length(up_term)), rep('blue', length(dn_term)));
		barplot(bar_h, col=col, border=col, width=0.8, space=0.25, names.arg=F);
		text(1:length(bar_h)-0.5, -cy*0.5, names(bar_h), srt=-45, adj=0, xpd=T)
	}
		dev.off();
	}
}
# low and high identity GF1--GF2
fh = data$nucl_iden[,20]>=quantile(data$nucl_iden[,20], 0.75);
fl = data$nucl_iden[,20]<=quantile(data$nucl_iden[,20], 0.25);
data$high_iden_go_test = list();
data$low_iden_go_test = list();
for (domain in c('GO2_BP', 'GO2_MF', 'KEGG')) {
	data$high_iden_go_test[[domain]] = list();
	data$high_iden_go_test[[domain]]$GF1_GF2 = go_test(gs[[domain]], rownames(data$lfpkm)[fh], rownames(data$lfpkm));
	data$low_iden_go_test[[domain]] = list();
	data$low_iden_go_test[[domain]]$GF1_GF2 = go_test(gs[[domain]], rownames(data$lfpkm)[fl], rownames(data$lfpkm));
}
# }}}

##########################################
# expression ~ identity
##########################################
# {{{
tmp=list();
tmp$a1 = data$max2[,4];
tmp$a2 = data$max2[,5];
tmp$a3 = data$mean2[,4];
tmp$a4 = data$mean2[,5];
tmp$iden1 = data$nucl_iden[,4];
tmp$iden2 = data$nucl_iden[,5];
tmp$a = data.frame(max_lfpkm=c(tmp$a1,tmp$a2),  mean_lfpkm=c(tmp$a3,tmp$a4), iden=c(tmp$iden1,tmp$iden2));
tmp$xx = c(60,seq(76,92,4),100);
tmp$lxx = length(tmp$xx);
tmp$max_lfpkm_by_iden=list();
tmp$mean_lfpkm_by_iden=list();
for (i in 2:tmp$lxx) {
	iden0 = tmp$xx[i-1];
	iden1 = tmp$xx[i];
	tmp$max_lfpkm_by_iden[[i-1]] = tmp$a$max_lfpkm[tmp$a$iden>iden0 & tmp$a$iden<=iden1 ];
	tmp$mean_lfpkm_by_iden[[i-1]] = tmp$a$mean_lfpkm[tmp$a$iden>iden0 & tmp$a$iden<=iden1 ];
}
names(tmp$max_lfpkm_by_iden) = c("<=76", paste0(tmp$xx[2:(tmp$lxx-1)],"-",tmp$xx[3:tmp$lxx]));

tmp$xx = c(0,seq(1,6),100);
tmp$lxx = length(tmp$xx);
tmp$iden_by_max_lfpkm=list();
tmp$iden_by_mean_lfpkm=list();
for (i in 2:tmp$lxx) {
	e0 = tmp$xx[i-1];
	e1 = tmp$xx[i];
	tmp$iden_by_max_lfpkm[[i-1]] = tmp$a$iden[tmp$a$max_lfpkm>e0& tmp$a$max_lfpkm<=e1 & tmp$a$iden>0];
	tmp$iden_by_mean_lfpkm[[i-1]] = tmp$a$iden[tmp$a$mean_lfpkm>e0& tmp$a$mean_lfpkm<=e1 & tmp$a$iden>0];
}
names(tmp$iden_by_max_lfpkm) = c("<=1", paste0(tmp$xx[2:(tmp$lxx-2)],"-",tmp$xx[3:(tmp$lxx-1)]), ">=6");
names(tmp$iden_by_mean_lfpkm) = c("<=1", paste0(tmp$xx[2:(tmp$lxx-2)],"-",tmp$xx[3:(tmp$lxx-1)]), ">=6");

res=300;
jpeg(paste0("output/plot/ZF_GF_max_lfpkm_by_iden_groups.jpg"), res=res, w=res*8, h=res*6);
boxplot(tmp$max_lfpkm_by_iden, xlab="Identity groups (%)", ylab="Maximal log2(FPKM+1) across six tissues", cex.lab=1.2);
dev.off();
jpeg(paste0("output/plot/ZF_GF_iden_by_max_lfpkm_groups.jpg"), res=res, w=res*8, h=res*6);
boxplot(tmp$iden_by_max_lfpkm, ylab="Identity (%)", xlab="Maximal log2(FPKM+1) across six tissues", cex.lab=1.2);
dev.off();
jpeg(paste0("output/plot/ZF_GF_iden_by_average_lfpkm_groups.jpg"), res=res, w=res*8, h=res*6);
boxplot(tmp$iden_by_mean_lfpkm, ylab="Identity (%)", xlab="Average log2(FPKM+1) across six tissues", cex.lab=1.2);
dev.off();
#boxplot(tmp$mean_lfpkm_by_iden);
tmp$f_high_iden = tmp$a$iden>=88;
tmp$f_low_iden = tmp$a$iden<85;
wilcox.test(tmp$a$max_lfpkm[tmp$f_high_iden], tmp$a$max_lfpkm[tmp$f_low_iden]);
tmp$f_high_exp = tmp$a$max_lfpkm>=3;
tmp$f_low_exp = tmp$a$max_lfpkm<=2;
wilcox.test(tmp$a$iden[tmp$f_high_exp], tmp$a$iden[tmp$f_low_exp]);
rm(tmp);

tmp$mat = matrix(0,2,2);
tmp$mat[1,1] = sum(data$max2[,5]-data$max2[,4]>=0.5 & data$max2[,4]<1 & data$nucl_iden[,4]<data$nucl_iden[,5]);
tmp$mat[1,2] = sum(data$max2[,5]-data$max2[,4]>=0.5 & data$max2[,4]<1 & data$nucl_iden[,4]>data$nucl_iden[,5]);
tmp$mat[2,1] = sum(data$max2[,4]-data$max2[,5]>=0.5 & data$max2[,5]<1 & data$nucl_iden[,4]<data$nucl_iden[,5]);
tmp$mat[2,2] = sum(data$max2[,4]-data$max2[,5]>=0.5 & data$max2[,5]<1 & data$nucl_iden[,4]>data$nucl_iden[,5]);
fisher.test(tmp$mat);

tmp$mat[1,1] = sum(data$pair_cor[,5]-data$pair_cor[,4]>=0.2 & data$nucl_iden[,4]<data$nucl_iden[,5] & data$nucl_iden[,4]<85);
tmp$mat[1,2] = sum(data$pair_cor[,5]-data$pair_cor[,4]>=0.2 & data$nucl_iden[,4]>data$nucl_iden[,5] & data$nucl_iden[,5]<85);
tmp$mat[2,1] = sum(data$pair_cor[,4]-data$pair_cor[,5]>=0.2 & data$nucl_iden[,4]<data$nucl_iden[,5] & data$nucl_iden[,4]<85);
tmp$mat[2,2] = sum(data$pair_cor[,4]-data$pair_cor[,5]>=0.2 & data$nucl_iden[,4]>data$nucl_iden[,5] & data$nucl_iden[,5]<85);
fisher.test(tmp$mat);

tmp$c11 = data$nucl_iden[data$max2[,5]-data$max2[,4]>=0.5 & data$max2[,4]<1,4]
tmp$c12 = data$nucl_iden[data$max2[,4]-data$max2[,5]>=0.5 & data$max2[,5]<1,5]
tmp$c1 = c(tmp$c11,tmp$c12);
tmp$c21 = data$nucl_iden[data$max2[,5]-data$max2[,4]>=0.5 & data$max2[,4]<1,5]
tmp$c22 = data$nucl_iden[data$max2[,4]-data$max2[,5]>=0.5 & data$max2[,5]<1,4]
tmp$c2 = c(tmp$c21,tmp$c22);
wilcox.test(tmp$c1,tmp$c2, paired = T)

tmp$c11 = data$nucl_iden[data$max2[,5]-data$max2[,4]>=0.5,4]
tmp$c12 = data$nucl_iden[data$max2[,4]-data$max2[,5]>=0.5,5]
tmp$c1 = c(tmp$c11,tmp$c12);
tmp$c21 = data$nucl_iden[data$max2[,5]-data$max2[,4]>=0.5,5]
tmp$c22 = data$nucl_iden[data$max2[,4]-data$max2[,5]>=0.5,4]
tmp$c2 = c(tmp$c21,tmp$c22);
wilcox.test(tmp$c1,tmp$c2, paired = T)
# }}}


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
# cluster individual genes
##################################
# {{{
tmp=list();
lfpkm = data1$log2_fpkm;
m=nrow(lfpkm);
tmp$names2 = as.vector(t(matrix(paste0( rep(data1$gene_ids2[,1],3), rep(c('','_1','_2'),each=m) ), m, 3)));
tmp$f_div = data1$f & !data1$fate$GF$double_high_conserved;
#tmp$f_div = data1$pair_cor[,2]<0.75 | data1$pair_cor[,3]<0.75 | data1$pair_cor[,8]<0.75 ; 
tmp$lfpkm6 = arrange_data(lfpkm[tmp$f_div,1:18], 6, data1$gene_ids2[tmp$f_div,1:3], tissues);
tmp$f = !as.vector(t(data1$low2_fpkm[tmp$f_div,1:3]));
tmp$lfpkm = tmp$lfpkm6[tmp$f,];
data1$genes2_for_h1 = tmp$f_div;
data1$genes_for_h1 = rep(0, m*3);
data1$genes_for_h1[rep(tmp$f_div,3)] = tmp$f;
h1=list();
h1$gene_ids2 = data1$gene_ids2[tmp$f_div,];
h1$data2=lfpkm[tmp$f_div,];
h1$data = tmp$lfpkm;
h1$cor_dist = cor_dist(h1$data);
h1$n = nrow(h1$data);
h1$n2 = nrow(h1$gene_ids2);
h1$h = hclust(h1$cor_dist, method='ward.D2');
h1$wss1 = choose_cluster_num_WSS(h1$data, h1$h, 200);
h1$nc = 20
res = hclust_more(h1$data, h1$h, h1$nc);
for (name in names(res)) { h1[[name]] = res[[name]]; }
rm(res);

h1$idx2 = matrix(0, h1$n2, 3);
rownames(h1$idx2) = rownames(h1$data2);
a = 1:h1$n;
names(a) = rownames(h1$data);
for (i in 1:3) {
	f = h1$gene_ids2[,i] %in% names(a);
	gid1 = h1$gene_ids2[f,1];
	gid2 = h1$gene_ids2[f,i];
	h1$idx2[gid1,i] = a[gid2];
	h1$idx2[!f,i] = 0;
}
rm(a,gid1,gid2);

##########################################
# For each cluster analysis
##########################################
# {{{
# compute number of shared homolog
h1$cid2 = matrix(0, nrow(h1$gene_ids2), 3);
h1$csz2 = matrix(0, h1$nc, 3);
rownames(h1$cid2) = h1$gene_ids2[,1];
for (i in 1:tot_copy) {
	f=h1$gene_ids2[,i] %in% rownames(h1$data);
	tmp$gids2 = h1$gene_ids2[f,1];
	tmp$gids1 = h1$gene_ids2[f,i];
	h1$cid2[tmp$gids2,i] = h1$cid[tmp$gids1];
	f1 = h1$cid2[,i]>0;
	a = tapply(rep(1,h1$n2)[f1], h1$cid2[f1,i], sum);
	h1$csz2[,i] = a;
}
h1$count = tapply(rep(1,nrow(h1$gene_ids2)), apply(h1$cid2[,1:3], 1, function(x) { paste(x,collapse='_');}), sum);
h1$shared_count = list();
h1$shared_frac  = list();
for (i in 1:(tot_copy-1)) {
for (j in (i+1):tot_copy) {
	mat = matrix(0, h1$nc, h1$nc);
	for (k in 1:nrow(h1$gene_ids2)) {
		mat[ h1$cid2[k,i], h1$cid2[k,j] ] =mat[ h1$cid2[k,i], h1$cid2[k,j] ] +1;
	}
	min = matrix(0, h1$nc, h1$nc);
	for (k1 in 1:h1$nc) { for (k2 in 1:h1$nc) {min[k1,k2] = min(h1$csz2[k1,i]+h1$csz2[k1,j], h1$csz2[k2,i]+h1$csz2[k2,j])} } 
	h1$shared_count[[(i-1)*tot_copy+j]] = mat;
	h1$shared_count[[(j-1)*tot_copy+i]] = mat;
	h1$shared_frac [[(i-1)*tot_copy+j]] = mat/min;
	h1$shared_frac [[(j-1)*tot_copy+i]] = mat/min;
}
}
h1$shared_count_p = list();
for (nn in 1:1000) {
	cid2 = h1$cid2;
	for (i in 1:tot_copy) { cid2[,i] = sample(cid2[,i], h1$n2, replace=F); }

	for (i in 1:(tot_copy-1)) {
	for (j in (i+1):tot_copy) {
		if (nn==1) { h1$shared_count_p[[(i-1)*tot_copy+j]] = matrix(0, h1$nc, h1$nc); }
		mat = matrix(0, h1$nc, h1$nc);
		for (k in 1:nrow(h1$gene_ids2)) {
			mat[ cid2[k,i], cid2[k,j] ] =mat[ cid2[k,i], cid2[k,j] ] +1;
		}
		h1$shared_count_p[[(i-1)*tot_copy+j]] = h1$shared_count_p[[(i-1)*tot_copy+j]] + (mat+t(mat)>(h1$shared_count[[(i-1)*tot_copy+j]])+t(h1$shared_count[[(i-1)*tot_copy+j]]));
	}
	}
	if (nn%%100==0) { print(nn); }
}
for (i in 1:(tot_copy-1)) {
for (j in (i+1):tot_copy) {
	h1$shared_count_p[[(j-1)*tot_copy+i]] = h1$shared_count_p[[(i-1)*tot_copy+j]] = h1$shared_count_p[[(i-1)*tot_copy+j]] /1000;
}
}

h1$near_cid_pair = NULL;
a1 = h1$parent[h1$cnode];
b1 = sort(a1);
b1 = b1[duplicated(b1)];
for (i in b1) {
	h1$near_cid_pair = rbind(h1$near_cid_pair, (1:h1$nc)[a1==i]);
}
h1$GF_same_count = sum(diag(h1$shared_count[[6]]));
h1$GF_near_count = sum(h1$shared_count[[6]][h1$near_cid_pair]) + sum(h1$shared_count[[6]][h1$near_cid_pair[,2:1]]);

h1$cnode_superc = rep(0, length(h1$cnode));
for (i in 1:h1$nc) {
	j = h1$cnode[i];
	while (h1$parent[j]!=h1$n*2-1) { j = h1$parent[j];}
	h1$cnode_superc[i]  = j;
}
h1$GF_superc_same_count=0;
h1$GF_superc_diff_count=0;
a = h1$shared_count[[6]];
for (i in 1:h1$nc) {
for (j in 1:h1$nc) {
	if (h1$cnode_superc[i]==h1$cnode_superc[j]) { h1$GF_superc_same_count=h1$GF_superc_same_count+a[i,j]; }
	else { h1$GF_superc_diff_count=h1$GF_superc_diff_count+a[i,j]; }
}
}


# compute mean of log2 FPKM cross cluster for each sample
h1$lfpkm_mean18 = matrix(0, h1$nc, tot_copy*nts);
for (i in 1:tot_copy) {
	cid2=h1$cid2[,i];
	for (j in 1:nts) {
		h1$lfpkm_mean18[,(i-1)*nts+j] = tapply(h1$data2[cid2!=0,(i-1)*nts+j], h1$cid2[cid2!=0,i], mean);
	}
}

h1$lfpkm_sd18 = matrix(0, h1$nc, tot_copy*nts);
for (i in 1:tot_copy) {
	cid2=h1$cid2[,i];
	for (j in 1:nts) {
		h1$lfpkm_sd18[,(i-1)*nts+j] = tapply(h1$data2[cid2!=0,(i-1)*nts+j], h1$cid2[cid2!=0,i], sd);
	}
}

# }}}

a = h1$shared_count[[6]]+t(h1$shared_count[[6]]);
draw_pair1 = which(a>=20,arr=T)
draw_pair1 = draw_pair1[draw_pair1[,1]<draw_pair1[,2],];
draw_pair = which(a>=10 & h1$shared_count_p[[6]]<0.01 ,arr=T)
draw_pair = draw_pair[draw_pair[,1]<draw_pair[,2],];
# random select 1% of goldfish genes for plot
# {{{
	frac = 200/h1$n2;
	h1s = list();
	h1$GF_idx400 = c();
	csz = h1$csz2[,2]+h1$csz2[,3];
	for (i in 1:h1$nc) {
		if (csz[i]<=5) {
			a = (1:h1$n2)[h1$cid2[,2]==i | h1$cid2[,3]==i ];
			h1$GF_idx400 = c(h1$GF_idx400, a);
		} else {
			a = sample((1:h1$n2)[h1$cid2[,2]==i], ceiling(h1$csz2[i,2]*frac), replace=F);
			b = sample((1:h1$n2)[h1$cid2[,3]==i], ceiling(h1$csz2[i,3]*frac), replace=F);
			h1$GF_idx400 = c(h1$GF_idx400, unique(c(a,b)));
		}
	}
	n0 = sum(h1$cid2[,2]==0);
	a = sample((1:h1$n2)[h1$cid2[,2]==0], ceiling(n0*frac), replace=F);
	n0 = sum(h1$cid2[,3]==0);
	b = sample((1:h1$n2)[h1$cid2[,3]==0], ceiling(n0*frac), replace=F);
	h1$GF_idx400 = c(h1$GF_idx400,unique(c(a,b)));
	h1$GF_idx400 = unique(sort(h1$GF_idx400));
	idx   = sort(as.vector(h1$idx2[h1$GF_idx400,2:3]));
	idx   = idx[idx>0];
	subtree = sub_hclust2(h1$h, idx);

	jpeg(paste0("output/plot/ZF_GF.GF_individual_gene.lfpkm_zscore1_subtree.400.k",h1$nc,".jpg"), res=res, w=res*12, h=res*6);
	plot(subtree);
	dev.off();


	h1s$gene_ids2 = h1$gene_ids2[h1$GF_idx400,];
	h1s$n2 = length(h1$GF_idx400);
	h1s$data2 = h1$data2[h1$GF_idx400,];
	h1s$data  = h1$data [idx,];
	h1s$n     = nrow(h1s$data);
	h1s$order = order(h1$rank[idx]);
	h1s$rank=h1s$order;
	h1s$rank[h1s$order] = 1:h1s$n;
	h1s$cid = h1$cid[idx];
	h1s$csz = tapply(rep(1,h1s$n), h1s$cid, sum);
	mean = apply(h1s$data, 1, mean);
	sd   = apply(h1s$data, 1, sd  );
	z    =(h1s$data-mean)/sd;
	z    = z[h1s$order,]; 
	max_z <- ceiling(max(max(z), -min(z)));

	tmp=list();
	tmp$c_ord = unique(h1s$cid[h1s$order]);
	tmp$c_rnk = tmp$c_ord;   tmp$c_rnk[tmp$c_ord] = 1:h1$nc;
	cl_y = rep(0, h1$nc+1);
	for (ci in tmp$c_ord) {
		cl_y[tmp$c_rnk[ci]+1]=cl_y[tmp$c_rnk[ci]]+h1s$csz[ci];
	}
	res=300;
	jpeg(paste0("output/plot/ZF_GF.GF_individual_gene.lfpkm_zscore1_heatmap.400.k",h1$nc,".jpg"), res=res, w=res*9, h=res*12);
	par(mar=c(0.5,8,0.5,10)+0.1);
	image(0:ncol(z), 0:nrow(z), t(z), col=bluered(255), axes=F, zlim=c(-max_z, max_z), xlab="", ylab="");
	cx = par('cxy')[1];  cy = par('cxy')[2];
	abline( v=(1:(tot_copy+tot_sp-1))*nts);
	abline( v=tot_copy*nts, lwd=3, col='black');
	y=0;
	lines(c(0,nts), c(0,0), col='darkgray', xpd=T);
	for (ci in unique(h1s$cid[h1s$order])) {
		y1=y+h1s$csz[ci];
#		rect(49,y,50,y1,col=color[ci],border=F,xpd=T);
		rect(-2*cx,y,-cx,y1,col=color[ci],border=F,xpd=T);
		lines(c(0,nts), c(y1,y1), col='darkgray', xpd=T);
#		lines(c(49,100), c(y1,y1), col='darkgray', xpd=T);
		y=y1;
	}

# draw links between clusters
	have_link = matrix(rep(0, h1$nc), h1$nc,1);
	ord = order(abs(draw_pair[,2]-draw_pair[,1]));
	for (i in ord) {
		j0=draw_pair[i,1];
		j1=draw_pair[i,2];
		mm = h1$shared_count[[6]][j0,j1];
		ff = h1$shared_frac[[6]][j0,j1];
		if (abs(j1-j0)>1) {
			for (k in 1:ncol(have_link)) {
				ovl = sum(have_link[(j0+1):(j1-1),k]);
				if (ovl==0) { break;}
			}
			if (ovl!=0) {
				have_link = cbind(have_link, rep(0,h1$nc));
				k=k+1;
			}
		} else {
			k=1;
		}
		if (abs(j1-j0)>1) { have_link[(j0+1):(j1-1),k] = 1; }
		y0 = (cl_y[j0]+cl_y[j0+1])/2;
		y1 = (cl_y[j1]+cl_y[j1+1])/2;
		dy = min(cy/2, (y1-y0)/4);
		y0 = y0+dy;
		y1 = y1-dy;
		plot_link(nts+cx*0.5, y0, nts+cx*0.5*(k+3), y1, cx/3, cy/4, 4, lwd=mm/4, xpd=T, col=color[1+(i-1)%%30]);
	}
	dev.off();
# }}}

# random select 1% genes for plot
# {{{
	frac = 400/h1$n2;
	h2s = list();
	h1$idx400 = c();
	for (i in 1:h1$nc) {
		if (h1$csz2[i,1]<=2) {
			a = (1:h1$n2)[h1$cid2[,1]==i];
			h1$idx400 = c(h1$idx400, a);
		} else {
			a = sample((1:h1$n2)[h1$cid2[,1]==i], ceiling(h1$csz2[i,1]*frac), replace=F);
			h1$idx400 = c(h1$idx400, a);
		}
	}
	n0 = sum(h1$cid2[,1]==0);
	a = sample((1:h1$n2)[h1$cid2[,1]==0], ceiling(n0*frac), replace=F);
	h1$idx400 = c(h1$idx400,a);
	h1$idx400 = sort(h1$idx400);
	h2s$gene_ids2 = h1$gene_ids2[h1$idx400,];
	h2s$n2 = length(h1$idx400);
	h2s$data2 = h1$data2[h1$idx400,];
	idx   = as.vector(h1$idx2[h1$idx400,]);
	idx   = idx[idx>0];
	h2s$data  = h1$data [idx,];
	h2s$n     = nrow(h2s$data);
	h2s$order = order(h1$rank[idx]);
	h2s$rank=h2s$order;
	h2s$rank[h2s$order] = 1:h2s$n;
	h2s$cid = h1$cid[idx];
	h2s$csz = tapply(rep(1,h2s$n), h2s$cid, sum);
	mean = apply(h2s$data, 1, mean);
	sd   = apply(h2s$data, 1, sd  );
	z    =(h2s$data-mean)/sd;
	z    = z[h2s$order,]; 
	max_z <- ceiling(max(max(z), -min(z)));

	tmp=list();
	tmp$c_ord = unique(h2s$cid[h2s$order]);
	tmp$c_rnk = tmp$c_ord;   tmp$c_rnk[tmp$c_ord] = 1:h1$nc;
	cl_y = rep(0, h1$nc+1);
	for (ci in tmp$c_ord) {
		cl_y[tmp$c_rnk[ci]+1]=cl_y[tmp$c_rnk[ci]]+h2s$csz[ci];
	}
	res=300;
	jpeg(paste0("output/plot/ZF_GF.individual_gene.lfpkm_zscore1_heatmap.400.k",h1$nc,".jpg"), res=res, w=res*9, h=res*12);
	par(mar=c(0.5,8,0.5,10)+0.1);
	image(0:ncol(z), 0:nrow(z), t(z), col=bluered(255), axes=F, zlim=c(-max_z, max_z), xlab="", ylab="");
	cx = par('cxy')[1];  cy = par('cxy')[2];
	abline( v=(1:(tot_copy+tot_sp-1))*nts);
	abline( v=tot_copy*nts, lwd=3, col='black');
	y=0;
	lines(c(0,nts), c(0,0), col='darkgray', xpd=T);
	for (ci in unique(h2s$cid[h2s$order])) {
		y1=y+h2s$csz[ci];
#		rect(49,y,50,y1,col=color[ci],border=F,xpd=T);
		rect(-2*cx,y,-cx,y1,col=color[ci],border=F,xpd=T);
		lines(c(0,nts), c(y1,y1), col='darkgray', xpd=T);
#		lines(c(49,100), c(y1,y1), col='darkgray', xpd=T);
		y=y1;
	}
#	k=0;
#	dx = min(0.2,cx*0.2);
#	lx = min(0.5,cx*0.5);
#	x0 = 50+dx;
#	for (ft in fate_names) {
#		k=k+1;
#		for (sp in c("CC","GF")) {
#			d = data$fate[[sp]][[ft]][h1$idx500][h2s$order];
#			x1 = x0+lx;
#			for (i in (1:h2s$n)[d]) {
#				rect(x0,i+0.2,x1,i+0.98,col=color[k],border=color[k],xpd=T);
#			}
#			x0 = x1+dx;
#		}
#	}

# draw expressed tissue names for each cluster
#	y=0;
#	for (ci in unique(h2s$cid[h2s$order])) {
#		y1=y+h2s$csz[ci];
#		a=h1$exp_tissue_gene_frac[[1]][ci,];
#		ts1=tissues1;
#		ts1[a<0.4]='-';
#		text( (-(nts+1)+(1:nts))*cx-2, rep((y+y1)/2,nts), ts1, xpd=T, adj=c(0.5,0.5), family='Arial', cex=1.2);
#		y=y1;
#	}
#	dev.off();

# draw links between clusters
	have_link = matrix(rep(0, h1$nc), h1$nc,1);
	ord = order(abs(draw_pair[,2]-draw_pair[,1]));
	for (i in ord) {
		j0=draw_pair[i,1];
		j1=draw_pair[i,2];
		mm = h1$shared_count[[6]][j0,j1];
		ff = h1$shared_frac[[6]][j0,j1];
		if (abs(j1-j0)>1) {
			for (k in 1:ncol(have_link)) {
				ovl = sum(have_link[(j0+1):(j1-1),k]);
				if (ovl==0) { break;}
			}
			if (ovl!=0) {
				have_link = cbind(have_link, rep(0,h1$nc));
				k=k+1;
			}
		} else {
			k=1;
		}
		if (abs(j1-j0)>1) { have_link[(j0+1):(j1-1),k] = 1; }
		y0 = (cl_y[j0]+cl_y[j0+1])/2;
		y1 = (cl_y[j1]+cl_y[j1+1])/2;
		dy = min(cy/2, (y1-y0)/4);
		y0 = y0+dy;
		y1 = y1-dy;
		plot_link(nts+cx*0.5, y0, nts+cx*0.5*(k+3), y1, cx/3, cy/4, 4, lwd=mm/10, xpd=T, col=color[1+(i-1)%%30]);
	}
	dev.off();
# }}}

data1$h1 = h1;
data1$h1_GF400 = h1s;
data1$h1_400 = h2s;
# }}}

save(data1, file='Rdata/ZF_GF_data.Rdata');

##################################
# cluster triplets
##################################
# {{{
lfpkm = data1$exp$fpkm$log2;
h=list();
h$f=f0;
lfpkm1 = lfpkm[h$f,];
data1$pca1 = princomp(lfpkm1[,1:(nts*tot_copy)], scores=T); # 4: 79.2% variance, 7: 90.3% variance
#data1$tsne1_k3 = tsne(lfpkm1, k=3);
#data1$tsne1_k2 = tsne(lfpkm1, k=2);
data1$clust = list();
#h$cor_dist = cor_dist(lfpkm1[,1:(nts*tot_copy)]);
h$lfpkm    = lfpkm1;
h$h = hclust(cor_dist(lfpkm1[,1:(nts*tot_copy)]), method="ward.D2");
h$wss1 = choose_cluster_num_WSS(h$lfpkm[,1:(nts*tot_copy)], h$h, 200);
res=hclust_more(h$lfpkm, h$h, nc=20);
for (name in names(res)) {h[[name]] = res[[name]]; }
# compute number of expressed genes in each cluster for each tissue
x=matrix(0, h$nc, (tot_copy+tot_sp)*nts);
colnames(x) = colnames(data1$fpkm);
for (i in 1:(tot_copy+tot_sp)) {
	j0 = nts*(i-1);
	for (j in 1:nts) {
		x[,j0+j] = tapply(h$lfpkm[,j0+j]>=data1$max2[h$f,i]-2 & data1$fpkm[h$f,j0+j]>=thres | h$lfpkm[,j0+j]>=3, h$cid[1:h$n], sum);
	}
}
h$c_exp_gene_count = x;
h$c_exp_gene_frac = h$c_exp_gene_count ;
h$c_exp_gene_frac = h$c_exp_gene_count / matrix(rep(h$csz,nts*(tot_copy+tot_sp)),h$nc,nts*(tot_copy+tot_sp));

# count ohno number/portion of fate categories in each cluster
all_fates = names(data1$fate$GF);
x = matrix(0, h$nc, length(all_fates));
colnames(x) = all_fates;
for (j in 1:length(all_fates)) {
	ft=all_fates[j];
	x[,j] = tapply(data1$fate$GF[[ft]][data1$f], h$cid[1:h$n], sum);
}
h$c_fate_gene_count = x;
h$c_fate_gene_frac = x/matrix(rep(h$csz,length(all_fates)),h$nc,length(all_fates));

# random select 400 genes for plot
# {{{
h$idx400 = c();
frac = 400/h$n;
for (i in 1:length(h$csz)) {
	if (h$csz[i]<=2) { h$idx400 = c(h$idx400, (1:h$n)[h$cid[1:h$n]==i]); }
	else { h$idx400 = c(h$idx400, sample((1:h$n)[h$cid[1:h$n]==i], ceiling(h$csz[i]*frac), replace=F)); }
}
hs = list();
hs$n = length(h$idx400);
hs$order = order(h$rank[h$idx400]);
hs$rank=hs$order;
hs$rank[hs$order] = 1:hs$n;
hs$cid = h$cid[h$idx400];
hs$csz = tapply(rep(1,hs$n), hs$cid, sum);

subtree = sub_hclust2(h$h, h$idx400);
jpeg(paste0("output/plot/ZF_GF.triplet.exp_zscore1_subtree.400.k",h$nc,".jpg"), res=res, w=res*10, h=res*5);
par(mar=c(0.5,0.5,0.5,0.5)+0.1);
plot(subtree);
dev.off();

swap = data1$pair_cor[h$f,2] < data1$pair_cor[h$f,3];
swap = swap[h$idx400];
z=data1$zscore1[h$f,][h$idx400,][hs$order,];
z[swap,(nts+1):(nts*3)] = z[swap, c((nts*2+1):(nts*3),(nts*1+1):(nts*2))]
max_z <- ceiling(max(max(z), -min(z)));

# key
image(seq(-max_z, max_z,length.out=255), 0:1, matrix(seq(-max_z, max_z,length.out=255),255,1), col=bluered(255), ann=F, axes=F);
axis(1);

res=300;
jpeg(paste0("output/plot/ZF_GF.triplet.exp_zscore1_heatmap.400.k",h$nc,".jpg"), res=res, w=res*9, h=res*10);
par(mar=c(0.5,8,0.5,0.5)+0.1);
image(0:ncol(z), 0:nrow(z), t(z), col=bluered(255), axes=F, zlim=c(-max_z, max_z), xlab="", ylab="");
cx = par('cxy')[1];
cy = par('cxy')[2];
abline( v=(1:(tot_copy+tot_sp-1))*nts);
abline( v=tot_copy*nts, lwd=3, col='black');
y=0;
lines(c(0,ncol(h$lfpkm)), c(0,0), col='darkgray', xpd=T);
#lines(c(nrow(h$lfpkm)+1,100), c(0,0), col='darkgray', xpd=T);
for (ci in unique(hs$cid[hs$order])) {
	y1=y+hs$csz[ci];
#		rect(49,y,50,y1,col=color[ci],border=F,xpd=T);
	rect(-2,y,-1,y1,col=color[ci],border=F,xpd=T);
	lines(c(0,ncol(h$lfpkm)), c(y1,y1), col='darkgray', xpd=T);
	y=y1;
}

# draw expressed tissue names for each cluster
y=0;
for (ci in unique(hs$cid[hs$order])) {
	y1=y+hs$csz[ci];
	a=h$c_exp_gene_frac[ci,1:nts];
	ts1=tissues1;
	ts1[a<0.4]='-';
	text( (-(nts+1)+(1:nts))*cx-2, rep((y+y1)/2,nts), ts1, xpd=T, adj=c(0.5,0.5), family='Arial', cex=1.2);
	y=y1;
}
dev.off();
# }}}
data1$h3 = h;

# }}}

save(data1, file='Rdata/ZF_GF_data.Rdata');



##############################
# old pairs
##############################
old=list();
old$ohno_pair <- read.table('~/data/goldfish/11549472/sergey_canu70x/arrow/ohnolog/fish4.cluster/czl_ohno_syn.out3/pair_from_chainnet.extended.txt', stringsAsFactor=F, head=F, sep="\t");
old$tpm <- read.table('~/data/goldfish/11549472/sergey_canu70x/arrow/RNA_star_run4/out/RSEM_out/mat/gene.tpm.ungroup.mat', stringsAsFactor=F);
#old$fpkm <- read.table('~/data/goldfish/11549472/sergey_canu70x/arrow/RNA_star_run4/out/RSEM_out/mat/gene.fpkm.ungroup.mat', stringsAsFactor=F);

tpm_thres=1;
n<-dim(old$ohno_pair)[1]
m<-dim(old$tpm)[2]-1;
colnames(old$ohno_pair) <- c('chr1','gid1','name1','chr2','gid2','name2', 'identity', 'score', 'cov1', 'cov2')
max1 <- apply(old$tpm[as.character(old$ohno_pair[,2]),2:(m+1)],1,max);
max2 <- apply(old$tpm[as.character(old$ohno_pair[,5]),2:(m+1)],1,max);
f <- (max1>=tpm_thres | max2>=tpm_thres) & annot$GF$exon_num[paste0(old$ohno_pair[,2],'_R1')]>=2 & annot$GF$exon_num[paste0(old$ohno_pair[,5], '_R1')]>=2  & annot$GF$tran_len[paste0(old$ohno_pair[,2],'_R1')]>=90 & annot$GF$tran_len[paste0(old$ohno_pair[,5],'_R1')]>=90
old$ohno_pair_f <- old$ohno_pair[f,]
n1 <- dim(old$ohno_pair_f)[1]
old$tpmA0 <- old$tpm[as.character(old$ohno_pair_f[,2]),];   old$tpmA <- old$tpmA0[,2:(m+1)]; old$ltpmA <- log2(old$tpmA+1);
old$tpmB0 <- old$tpm[as.character(old$ohno_pair_f[,5]),];   old$tpmB <- old$tpmB0[,2:(m+1)]; old$ltpmB <- log2(old$tpmB+1);
m <- dim(old$tpmA)[2];
sm_names <- colnames(old$tpmA);
sm_names <- colnames(old$tpmB);

a1 <- grep("Brain", colnames(old$tpmA));
a2 <- grep("Gill", colnames(old$tpmA));
old$tpm7A <- cbind( apply(old$tpmA[,a1],1,max), old$tpmA[,2], apply(old$tpmA[,a2],1,max), old$tpmA[,5], old$tpmA[,8:10]);
colnames(old$tpm7A) <- c("Brain", "Eye", "Gill", "Bone", "Heart", "Muscle", "TailFin");
old$ltpm7A <- log2(old$tpm7A+1);
old$tpm7B <- cbind( apply(old$tpmB[,a1],1,max), old$tpmB[,2], apply(old$tpmB[,a2],1,max), old$tpmB[,5], old$tpmB[,8:10]);
colnames(old$tpm7B) <- c("Brain", "Eye", "Gill", "Bone", "Heart", "Muscle", "TailFin");
old$ltpm7B <- log2(old$tpm7B+1);
names <- paste(rownames(old$tpmA), rownames(old$tpmB), sep="__"); 
old$gene_names<- as.character(old$tpmA0[,1]);
old$gene_names[ old$gene_names==rownames(old$tpmA) ] = names[ old$gene_names==rownames(old$tpmA) ];
old$gene_names[ old$gene_names=="_" ] = names[ old$gene_names=="_" ];
old$gene_names2 <- c( paste(old$gene_names, rownames(old$tpmA), sep="|A|"), paste(old$gene_names, rownames(old$tpmB), sep="|B|") );

# count number of co-express tissue
a1=old$tpm7A;
a2=old$tpm7B;
old$pair_coexp_tissue_num = apply(a1>=tpm_thres & a2>=tpm_thres, 1, sum);
n1 = nrow(a1);
tapply(rep(1,n1), old$pair_coexp_tissue_num, sum)
old$pair_cor = pair_cor(as.matrix(old$ltpm7A), as.matrix(old$ltpm7B));
old$pair_cor_p = pair_cor_p(as.matrix(old$ltpm7A), as.matrix(old$ltpm7B));





##################################
# select divergence genes for cluster
##################################
#tmp_lfpkm_max = apply(data$lfpkm[,1:30], 1, max);
#tmp_min_cor = apply(data$pair_cor[,c(2:5,1*8+3,3*8+5)],  1,min);
#data$filter2 = tmp_lfpkm_max>=lthres & tmp_min_cor<C2;
#data2 = list();
#for (name in c("lfpkm", "gene_ids2", "gene_names2", "zscore1", "zscore2")) {
#	data2[[name]] = data[[name]][data$filter2,];
#}
#

#data$sm_cor_dist = cor_dist(t(data$lfpkm[,1:30]));
#data$sm_h1 = hclust(data$sm_cor_dist, method="ward.D2");
lfpkm1 = lfpkm[data1$f,];
data$pca1 = princomp(data1$lfpkm[data1$f,1:30], scores=T); # 4: 78.8% variance, 9: 90.18% variance
data$tsne1_k4 = tsne(data$lfpkm[,1:30], k=4);
data$tsne1_k3 = tsne(data$lfpkm[,1:30], k=3);
data$tsne1_k2 = tsne(data$lfpkm[,1:30], k=2);
data$kmc1 = kmeans(data$lfpkm[,1:30], 1000);
data$cor_dist = cor_dist(data$lfpkm[,1:30]);
data$clust = list();
h1=list();
h1$h = hclust(data$cor_dist, method="ward.D2");
#data$h2 = hclust(data$cor_dist, method="average");
tmp_wss1 = choose_cluster_num_WSS(data$lfpkm[,1:30], h1$h, 200);
#tmp_wss2 = choose_cluster_num_WSS(data$lfpkm, data$h2, 100);

process_cluster = function(h1, nc)
# {{{
{
	h1$nc = nc;
	h1$cid = cutree(h1$h, k=h1$nc);
	h1$csz = tapply(rep(1,nrow(data$lfpkm)), h1$cid, sum);
	h1$dnd = as.dendrogram(h1$h);
	h1$centers = matrix(0, h1$nc, ncol(data$lfpkm));
	h1$order = h1$h$order;
	h1$rank=h1$h$order;    h1$rank[h1$h$order] = 1:data$n;
	colnames(h1$centers) = colnames(data$lfpkm);
	for (i in 1:ncol(data$lfpkm)) { h1$centers[,i] = tapply(data$lfpkm[,i], h1$cid, mean); }
# count number of ohno in each cluster that express in each tissue
	x=list();
	for (i in 1:(tot_copy+tot_sp)) {
		x[[i]] = matrix(0, h1$nc, nts);
		colnames(x[[i]]) = tissues;
		j0 = nts*(i-1);
	for (j in 1:nts) {
		x[[i]][,j] = tapply(data$lfpkm[,j0+j]>=data$max2[,i]-2 & data$lfpkm[,j0+j]>=data$max2[,i]/2 & data$lfpkm[,j0+j]>=lthres | data$lfpkm[,j0+j]>=4, h1$cid, sum);
	}
	}
	h1$exp_tissue_gene_count = x;
	h1$exp_tissue_gene_frac = h1$exp_tissue_gene_count ;
	for (i in 1:(tot_copy+tot_sp)) {
		h1$exp_tissue_gene_frac[[i]] = x[[i]]/matrix(rep(h1$csz,nts),h1$nc,nts);
	}

# count ohno number/portion of fate categories in each cluster
	x=list();
	x$CC=x$GF=matrix(0, h1$nc, length(fate_names));
	colnames(x$CC) = fate_names;
	colnames(x$GF) = fate_names;
	for (j in 1:length(fate_names)) {
		ft=fate_names[j];
		x$CC[,j] = tapply(data$fate$CC[[ft]], h1$cid, sum);
		x$GF[,j] = tapply(data$fate$GF[[ft]], h1$cid, sum);
	}
	h1$fate_gene_count = x;
	h1$fate_gene_frac = list();
	h1$fate_gene_frac$CC = x$CC/matrix(rep(h1$csz,length(fate_names)),h1$nc,length(fate_names));
	h1$fate_gene_frac$GF = x$GF/matrix(rep(h1$csz,length(fate_names)),h1$nc,length(fate_names));
}
# }}}

# random select 500 genes for plot
# {{{
	h1s = list();
	h1$idx500 = c();
	for (i in 1:length(h1$csz)) {
		h1$idx500 = c(h1$idx500, sample((1:data$n)[h1$cid==i], ceiling(h1$csz[i]*0.1), replace=F));
	}
	h1s$n = length(h1$idx500);
	h1s$order = order(h1$rank[h1$idx500]);
	h1s$rank=h1s$order;
	h1s$rank[h1s$order] = 1:h1s$n;
	h1s$cid = h1$cid[h1$idx500];
	h1s$csz = tapply(rep(1,h1s$n), h1s$cid, sum);
	z=data$zscore1[h1$idx500,][h1s$order,];
	max_z <- ceiling(max(max(z), -min(z)));

	res=300;
	jpeg(paste0("output/plot/exp_zscore1_heatmap.k",h1$nc,".jpg"), res=res, w=res*12, h=res*9);
	par(mar=c(0.5,8,0.5,0.5)+0.1);
	image(0:ncol(z), 0:nrow(z), t(z), col=bluered(255), axes=F, zlim=c(-max_z, max_z), xlab="", ylab="");
	cx = par('cxy')[1];
	cy = par('cxy')[2];
	abline( v=(1:(tot_copy+tot_sp-1))*nts);
	abline( v=tot_copy*nts, lwd=3, col='black');
	y=0;
	lines(c(0,48), c(0,0), col='darkgray', xpd=T);
	lines(c(49,100), c(0,0), col='darkgray', xpd=T);
	for (ci in unique(h1s$cid[h1s$order])) {
		y1=y+h1s$csz[ci];
#		rect(49,y,50,y1,col=color[ci],border=F,xpd=T);
		rect(-2,y,-1,y1,col=color[ci],border=F,xpd=T);
		lines(c(0,48), c(y1,y1), col='darkgray', xpd=T);
#		lines(c(49,100), c(y1,y1), col='darkgray', xpd=T);
		y=y1;
	}
#	k=0;
#	dx = min(0.2,cx*0.2);
#	lx = min(0.5,cx*0.5);
#	x0 = 50+dx;
#	for (ft in fate_names) {
#		k=k+1;
#		for (sp in c("CC","GF")) {
#			d = data$fate[[sp]][[ft]][h1$idx500][h1s$order];
#			x1 = x0+lx;
#			for (i in (1:h1s$n)[d]) {
#				rect(x0,i+0.2,x1,i+0.98,col=color[k],border=color[k],xpd=T);
#			}
#			x0 = x1+dx;
#		}
#	}

# draw expressed tissue names for each cluster
	y=0;
	for (ci in unique(h1s$cid[h1s$order])) {
		y1=y+h1s$csz[ci];
		a=h1$exp_tissue_gene_frac[[1]][ci,];
		ts1=tissues1;
		ts1[a<0.4]='-';
		text( (-(nts+1)+(1:nts))*cx-2, rep((y+y1)/2,nts), ts1, xpd=T, adj=c(0.5,0.5), family='Arial', cex=1.2);
		y=y1;
	}
	dev.off();
# }}}

data$clust[[1]] = h1;


heatmap.2()
f=h1$cid==1;
z=data$zscore1[,];
max_z <- ceiling(max(max(z), -min(z)));
image(0:ncol(z), 0:nrow(z), t(z), col=bluered(255), axes=F, zlim=c(-max_z, max_z), xlab="", ylab="");
abline( v=(1:(tot_copy+tot_sp-1))*nts);



h2=list();
h2$h = hclust(dist(data$lfpkm,method="maximum"), method="ward.D2");
tmp_wss2 = choose_cluster_num_WSS(data$lfpkm, data$h2, 100);
h2$nc = 22;
h2$cid= cutree(h2$h, k=h2$nc);
h2$csz= tapply(rep(1,nrow(data2$lfpkm)), h2$cid, sum);
data2$clust[[2]] = h2;
save(data2, file="Rdata/data2.Rdata");

###### cluster by sample ######
ltpm_cor_dist_sm <- cor_dist( t(ltpm1) );
h1c_sm <- hclust( ltpm_cor_dist_sm,  method="ward.D2" );
ord_sm <- h1c_sm$order;
#ltpm_euc_dist_sm <- dist( t(data1) );
#h1e_sm <- hclust( ltpm_euc_dist_sm,  method="ward.D2" );
######

##################################
# cluster
##################################
# {{{
f.ZF=max_tpm1[(1:(tot_gene/5))*5-4]>=tpm_thres; # ZF express
f.CC1=max_tpm1[(1:(tot_gene/5))*5-3]>=tpm_thres; # CC1 express
f.CC2=max_tpm1[(1:(tot_gene/5))*5-2]>=tpm_thres; # CC2 express
f.GF1=max_tpm1[(1:(tot_gene/5))*5-1]>=tpm_thres; # GF1 express
f.GF2=max_tpm1[(1:(tot_gene/5))*5]>=tpm_thres; # GF2 express

# most gene lost one copy expression in GF and(or) CC are shared, 
# i.e. they lost the expression in their common ancestry
a=matrix(0,2,2); # 543,933,684,7756, P-value ***
a[1,1]=sum(f.ZF & !(f.CC1 & f.CC2) & !(f.GF1 & f.GF2));
a[1,2]=sum(f.ZF & !(f.CC1 & f.CC2) & (f.GF1 & f.GF2));
a[2,1]=sum(f.ZF & (f.CC1 & f.CC2) & !(f.GF1 & f.GF2));
a[2,2]=sum(f.ZF & (f.CC1 & f.CC2) & (f.GF1 & f.GF2));
#### 

f = rep(f.ZF,each=5) & apply(tpm1>=tpm_thres,1,sum)>0
data1 = as.data.frame(norm_ltpm1[f,]);
data1.gene_name <- rownames(tpm1)[f];
data1.sp <- tpm_sp[f];
data1.n <- dim(data1)[1];
data1.m <- dim(data1)[2];
data1.rowmean <- apply(data1, 1, mean);
data1.rowsd   <- apply(data1, 1, sd);
norm_data1 <- (data1-data1.rowmean)/data1.rowsd;

#data1.cov = cov(t(norm_data1));
#ltpm_cor_dist <- cor_dist( norm_data1 );

norm_data1.princomp=princomp(norm_data1)
norm_data1.kms2000=kmeans(norm_data1, centers=2000)
norm_data1.kms4000=kmeans(norm_data1, centers=4000)
centers1 = norm_data1.kms4000$centers;
centers1.cor_dist=cor_dist(centers1);
centers1.h1 = hclust(centers1.cor_dist, method='average');
centers1.h2 = hclust(centers1.cor_dist, method='ward.D2');

save(data1, norm_data1, norm_data1.princomp, norm_data1.kms4000, 
		file=paste0(wd,'/Rdata/cluster.Rdata'));

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



