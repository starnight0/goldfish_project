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
wd="~/data/goldfish/11549472/sergey_canu70x/arrow/carAur03/big/WGD2";
setwd(wd);

tpm_thres=1;
tissues=c("Brain", "Eye", "Heart", "Gill", "Muscle", "Tail");
tissues1=c('B', 'E', 'H', 'G', 'M', 'T');
nts = length(tissues);
sps=c('ZF', 'CC', 'GF'); # DONOT change the order
sps2=c('ZF', 'CC', 'CC', 'GF', 'GF'); # DONOT change the order
sps3=c('ZF', 'CC1', 'CC2', 'GF1', 'GF2'); # DONOT change the order
copies=c(1,2,2);
names(copies)=sps;
tot_copy=5;
tot_sp = length(sps);
fate_names=c('double_conserved', 'dosage_balance','subfunc','neofunc','nonfunc')
fate_names1=c('double-conserved', 'dosage-balance','sub-func','neo-func','non-func')
names(fate_names1) = fate_names;

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
	for (sp in c('GF','CC','ZF')) {
		rownames(bgp[[sp]]) = bgp[[sp]][,4]; 
		a1 = list();
		a1$gene_names = bgp[[sp]][,18];
		a1$exon_len=strsplit(bgp[[sp]][,11], ",", fixed=T);
		a1$exon_len <- lapply(a1$exon_len, as.numeric);
		a1$exon_num <- bgp[[sp]][,10];
		a1$tran_len <- sapply(a1$exon_len, sum);
		annot[[sp]] = a1;
	}
	save(bgp, annot, file="Rdata/annot.Rdata")
}
# }}}
#############################


###########################
# Load clusters
###########################
# {{{
dfile=paste0(wd,"/Rdata/anchor.Rdata");
if (file.exists(dfile)) {
	load("Rdata/anchor.Rdata");
} else {
	file=paste0(wd, "/input/anchor.fix");
	anchor=read.table(file, head=T, sep='\t');
#	tmp.clust = read.table("input/rescue_m1.6.cluster.txt", sep='\t', head=F);
#	f.ZF = grep('\\(ZF', tmp.clust[,4]);
#	f.CC = grep('\\(CC', tmp.clust[,4]);
#	f.GF = grep('\\(GF', tmp.clust[,4]);
	
	anchor[,'anchor_tid'] = gsub("\\.[0-9]+", "", anchor[,'anchor_tid'] );
	anchor[,'anchor_gid'] = gsub("\\.[0-9]+", "", anchor[,'anchor_gid'] );
#	anchor = anchor[anchor[,"clust_id"]>=0 | anchor[,'to_dup']!=".",];
	rownames(anchor) = paste0(anchor[,'species'],'|',anchor[,'anchor_gid']); # name is gene id
	anchor[,'clust_id'] = as.numeric(anchor[,'clust_id']);
	ncluster = max(anchor[,'clust_id'])+1;
#	cluster_is_dup = rep(F,ncluster); 
#	cluster_is_dup[f.ZF] = T;
#	cluster_is_dup[f.CC] = T;
#	cluster_is_dup[f.GF] = T;
	cluster = list();
	length(cluster)=ncluster;
	cluster_sizes=matrix(0, ncluster, 5);
	colnames(cluster_sizes) = c('total', 'ZF', 'GC', 'CC', 'GF');
	for (i in 1:nrow(anchor)) {
		cid=anchor[i,'clust_id']+1;
		if (cid==0) { next; }
		sp = anchor[i,'species'];
		if (is.null(cluster[[cid]])) cluster[[cid]]=list();
		if (anchor[i,'to_dup']==".") {
			cluster[[cid]][['total']] = c(cluster[[cid]][['total']], i);
			cluster[[cid]][[sp]] = c(cluster[[cid]][[sp]], i);
			cluster_sizes[cid, sp]=cluster_sizes[cid, sp]+1;
			cluster_sizes[cid, 'total']=cluster_sizes[cid, 'total']+1;
		} 
	}
	anchor_dup = list();
	length(anchor_dup)=nrow(anchor);
	names(anchor_dup) = rownames(anchor);
	f=anchor[,'to_dup']!=".";
	dup = paste0(anchor[,'species'],'|',anchor[,'to_dup']);
	for (i in (1:nrow(anchor))[f]) {
		anchor_dup[[dup[i]]] = c(anchor_dup[[dup[i]]], i); 
	}
	for (i in 1:nrow(anchor)) {
		if (length(anchor_dup[[i]])>0) {
			names(anchor_dup[[i]]) = rownames(anchor)[ anchor_dup[[i]] ];
		}
	}
#	cluster = split(1:nrow(anchor), anchor[,'clust_id']);

	file=paste0(wd, "/input/qC50.loci.txt");
	anchor2=read.table(file, head=T, sep='\t');
	anchor2[,'anchor_gid'] = gsub("\\.[0-9]+", "", anchor2[,'anchor_gid'] );
	anchor2 = anchor2[anchor2[,"clust_id"]>=0,];
	rownames(anchor2) = paste0(anchor2[,'species'],'|',anchor2[,'anchor_gid']); # name is gene id
	anchor2[,'clust_id'] = as.numeric(anchor2[,'clust_id']);
	ncluster2 = max(anchor2[,'clust_id'])+1;
	cluster2 = list();
	length(cluster2)=ncluster2;
	cluster_sizes2=matrix(0, ncluster2, 5);
	colnames(cluster_sizes2) = c('total', 'ZF', 'GC', 'CC', 'GF');
	for (i in 1:nrow(anchor2)) {
		cid=anchor2[i,'clust_id']+1;
		if (cid==0) { next; }
		sp = anchor2[i,'species'];
		if (is.null(cluster2[[cid]])) cluster2[[cid]]=list();
		cluster2[[cid]][[sp]] = c(cluster2[[cid]][[sp]], i);
		cluster_sizes2[cid, sp]=cluster_sizes2[cid, sp]+1;
		cluster_sizes2[cid, 'total']=cluster_sizes2[cid, 'total']+1;
	}

	cluster_is_changed = rep(F,ncluster); # 0: no change, 1: changed by anchor2
	for (cid in 1:ncluster) {
		if (is.null(cluster2[[cid]])) { cluster_is_changed[cid]=T; next; }
		if (cid>ncluster2) { cluster_is_changed[cid]=T; next; }
		m=0;
		for (j in 1:tot_sp) {
			if (cluster_sizes[cid,sps[j]]==cluster_sizes2[cid,sps[j]]) {m=m+1; }
			else {break;}
		}
		if (m<tot_sp) cluster_is_changed[cid]=T;
	}
	rm(anchor2, cluster2, ncluster2, cluster_sizes2);

	edge = list(prot=read.table('input/blastp.ZF_GC_CC_GF.gene.f3.m6.gz', head=F),
			nucl= read.table('input/blastn.ZF_GC_CC_GF.gene.f3.m6.gz', head=F));
	rownames(edge$prot) = paste0(edge$prot[,1], "-", edge$prot[,2]);
	rownames(edge$nucl) = paste0(edge$nucl[,1], "-", edge$nucl[,2]);
	cluster_is_dup=rep(F, ncluster);
	cluster_is_dup3=matrix(rep(F, ncluster*tot_sp), ncluster, tot_sp);
	colnames(cluster_is_dup3) = sps;
	for (i in 1:ncluster) {
		if (cluster_sizes[i,'ZF']==0) { next; }
		zf_idx = cluster[[i]][['ZF']][1];
		if (length(anchor_dup[[zf_idx]])>0) {
			cluster_is_dup[i]=T;
			cluster_is_dup3[i,'ZF']=T;
			next;
		}
		for (sp in c('CC','GF')) {
			bad=F;
			idxs = cluster[[i]][[sp]];
			for (idx in idxs) {
				if (length(anchor_dup[[idx]])>0) {
					id_pair = paste0(rownames(anchor)[zf_idx], '-', rownames(anchor)[idx]);
					cov0 = edge$nucl[id_pair,18];
					if (is.na(cov0) || is.null(cov0)) { bad=T; }
					else {
						id_pairs = paste0(rownames(anchor)[zf_idx], '-', rownames(anchor)[anchor_dup[[idx]]]);
						cov = max(edge$nucl[id_pairs,18]);
						if (!is.na(cov) && cov>25 && cov>=cov0/2) { bad=T; }
					}
					if (bad) { break;}
				}
			}
			cluster_is_dup3[i,sp]=bad;
		}
		if (cluster_is_dup3[i,'CC'] || cluster_is_dup3[i,'GF'] ) { cluster_is_dup[i]=T; }
	}

	cluster_is_syn = rep(F,ncluster);
	tmp_idxs=unlist(sapply(1:ncluster,function(x) {cluster[[x]][['total']]}));
	tmp_cids=unlist(sapply(1:ncluster,function(x) {rep(x,cluster_sizes[x,'total'])}));
	cluster_is_syn[unique(tmp_cids)] = tapply(tmp_idxs, tmp_cids, function(x) {sum(anchor[x,'to_block_each_species']!="")>=3}, default=F);
	save(anchor,anchor_dup, ncluster,cluster,cluster_sizes, cluster_is_changed, cluster_is_syn, cluster_is_dup, edge, file="Rdata/anchor.Rdata");
}
# }}}
###########################


###########################
# Load TPM
# {{{
dfile=paste0(wd,'/Rdata/tpm.Rdata');
if ( file.exists(dfile) ) {
	load("Rdata/tpm.Rdata");
} else {
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

	cluster_is_exp_bad=rep(F,ncluster); # low or NA
	tmp_idxs=unlist(sapply(1:ncluster,function(x) {cluster[[x]][['total']]}));
	tmp_cids=unlist(sapply(1:ncluster,function(x) {rep(x,cluster_sizes[x,'total'])}));
	tmp_sps  = anchor[tmp_idxs,'species'];
	tmp_gids = anchor[tmp_idxs,'anchor_gid'];
	a=grep('^GC',tmp_sp_gids, invert=T);
	tmp_sps = tmp_sps[a];
	tmp_gids = tmp_gids[a];
	tmp_idxs = tmp_idxs[a];
	tmp_cids = tmp_cids[a];
	tmp_sp_gids=paste0(tmp_sps, '|', tmp_gids);
	tmp_fpkm_max = apply(exp$fpkm$raw[match(tmp_sp_gids, rownames(exp$fpkm$raw)),], 1, max);
	tmp_tpm_max = apply(exp$tpm$raw[match(tmp_sp_gids, rownames(exp$tpm$raw)),], 1, max);
	tmp_max = apply(cbind(tmp_fpkm_max, tmp_tpm_max), 1, min);
	cluster_is_exp_bad[unique(tmp_cids)] = tapply(tmp_max, tmp_cids, function(x) {sum(is.na(x))>0 | sum(x>=exp$fpkm$thres)<2}, default=T);


	cluster_is_good_double2=!cluster_is_exp_bad & !cluster_is_changed & !cluster_is_dup & cluster_sizes[,'ZF']==1 & cluster_sizes[,'CC']==2 & cluster_sizes[,'GF']==2;
	cluster_is_good_syn_double2=cluster_is_good_double2 & cluster_is_syn;
	cluster_is_good_CC0GF1 =!cluster_is_exp_bad & !cluster_is_changed & !cluster_is_dup & cluster_sizes[,'ZF']==1 & cluster_sizes[,'CC']==0 & cluster_sizes[,'GF']==1;
	cluster_is_good_CC0GF2 =!cluster_is_exp_bad & !cluster_is_changed & !cluster_is_dup & cluster_sizes[,'ZF']==1 & cluster_sizes[,'CC']==0 & cluster_sizes[,'GF']==2;
	cluster_is_good_CC1GF0 =!cluster_is_exp_bad & !cluster_is_changed & !cluster_is_dup & cluster_sizes[,'ZF']==1 & cluster_sizes[,'CC']==1 & cluster_sizes[,'GF']==0;
	cluster_is_good_CC1GF1 =!cluster_is_exp_bad & !cluster_is_changed & !cluster_is_dup & cluster_sizes[,'ZF']==1 & cluster_sizes[,'CC']==1 & cluster_sizes[,'GF']==1;
	cluster_is_good_CC1GF2 =!cluster_is_exp_bad & !cluster_is_changed & !cluster_is_dup & cluster_sizes[,'ZF']==1 & cluster_sizes[,'CC']==1 & cluster_sizes[,'GF']==2;
	cluster_is_good_CC2GF0 =!cluster_is_exp_bad & !cluster_is_changed & !cluster_is_dup & cluster_sizes[,'ZF']==1 & cluster_sizes[,'CC']==2 & cluster_sizes[,'GF']==0;
	cluster_is_good_CC2GF1 =!cluster_is_exp_bad & !cluster_is_changed & !cluster_is_dup & cluster_sizes[,'ZF']==1 & cluster_sizes[,'CC']==2 & cluster_sizes[,'GF']==1;

	cid_lists = list(
			CC2GF2=(1:ncluster)[cluster_is_good_double2],
			syn_CC2GF2=(1:ncluster)[cluster_is_good_syn_double2],
			CC0GF1=(1:ncluster)[cluster_is_good_CC0GF1],
			CC0GF2=(1:ncluster)[cluster_is_good_CC0GF2],
			CC1GF0=(1:ncluster)[cluster_is_good_CC1GF0],
			CC1GF1=(1:ncluster)[cluster_is_good_CC1GF1],
			CC1GF2=(1:ncluster)[cluster_is_good_CC1GF2],
			CC2GF0=(1:ncluster)[cluster_is_good_CC2GF0],
			CC2GF1=(1:ncluster)[cluster_is_good_CC2GF1]
			);
		
	for (type in c('tpm', 'fpkm', 'count')) {
		type1="syn_CC2GF2";
		cids = cid_lists[[type1]];
		tpm2 = matrix(0, length(cids), tot_copy*nts);
		gene_ids2 = c();
		k0 = 0;
		for (sp in sps) {
			for (i in 1:copies[sp]) {
				tmp_idxs=unlist(sapply(cids,function(x) {cluster[[x]][[sp]][i]}));
				tmp_cids=cids;
				tmp_gids = anchor[tmp_idxs,'anchor_gid'];
				sp_gids=paste0(sp, '|', tmp_gids);
				tpm2[,(k0+1):(k0+nts)] = as.matrix(exp[[type]][['raw']][sp_gids,]);
				gene_ids2 = cbind(gene_ids2, tmp_gids);
				k0 = k0+nts;
			}
		}
		rownames(tpm2) = gene_ids2[,1];
		colnames(tpm2)=c( paste0(rep('ZF',nts), '.', tissues),
			paste0(rep('CC1',nts), '.', tissues), paste0(rep('CC2',nts), '.', tissues),
			paste0(rep('GF1',nts), '.', tissues), paste0(rep('GF2',nts), '.', tissues) );
		ltpm2 = log2(tpm2+tpm_thres);

# sum of two paralog
		sum_tpm2 = matrix(0, nrow(tpm2), tot_sp*nts);
		colnames(sum_tpm2) = paste0("sum_of_", rep(sps,each=nts), '.', rep(tissues,tot_sp));
		j0 = 0;
		for (i in 1:tot_sp) {
			for (j in 1:copies[sps[i]]) {
				j0 = j0+1;
				sum_tpm2[,(nts*(i-1)+1):(i*nts)] = sum_tpm2[,(nts*(i-1)+1):(i*nts)] + tpm2[,(nts*(j0-1)+1):(j0*nts)];
			}
		}
		tpm2 = cbind(tpm2[,1:(nts*tot_copy)], sum_tpm2);
		ltpm2 = log2(tpm2+tpm_thres);

		qnorm_ltpm2 = normalizeBetweenArrays(ltpm2[,1:(nts*tot_copy)]);
		sum_tpm2 = matrix(0, nrow(tpm2), tot_sp*nts);
		colnames(sum_tpm2) = paste0("sum_of_", rep(sps,each=nts), '.', rep(tissues,tot_sp));
		j0 = 0;
		for (i in 1:tot_sp) {
			for (j in 1:copies[sps[i]]) {
				j0 = j0+1;
				sum_tpm2[,(nts*(i-1)+1):(i*nts)] = sum_tpm2[,(nts*(i-1)+1):(i*nts)] + 2^qnorm_ltpm2[,(nts*(j0-1)+1):(j0*nts)]-tpm_thres;
			}
		}
		qnorm_ltpm2 = cbind(qnorm_ltpm2[,1:(nts*tot_copy)], log2(sum_tpm2+tpm_thres));

		exp[[type]][['raw2']] = tpm2;
		exp[[type]][['lraw2']] = ltpm2;
		exp[[type]][['gene_ids2']] = gene_ids2;
		exp[[type]][['qnorm_lraw2']] = qnorm_ltpm2;
	}

#		tpm1 = matrix(0, length(cids)*tot_copy, nts);
#		ltpm1 = log2(tpm1+tpm_thres);
#		for (i in 1:nts) {
#			tpm1[,i] = as.vector(t(tpm2[, (1:tot_copy)*nts-nts+i]));
#		}
#		tpm1.gene_ids = as.vector(t(gene_ids2));
#		colnames(tpm1) = tissues;
#		rownames(tpm1) = paste0(rep(sps2,nrow(tpm2)),'|', tpm1.gene_ids);

###########################
# compute scores by species
# {{{
	for (type in c('tpm', 'fpkm', 'count')) {
		ltpm2 = exp[[type]][['lraw2']];
		thres = exp[[type]][['thres']];
	
# count expressed sample for each gene, layout: cluster by ortholog/paralog
		exp_sm_count <- matrix(0, nrow(ltpm2), tot_copy+tot_sp);
		for (i in 1:ncol(exp_sm_count)) {
			k0=(nts)*(i-1)+1;
			k1=(nts)*i;
			exp_sm_count[,i] <- apply(ltpm2[,k0:k1]>=thres, 1, sum);
		}
		exp[[type]][['exp_sm_count']] = exp_sm_count;
# tissue-specific score (<=1, near 1 is single tissue-specific)
		specific_score = matrix(0, nrow(ltpm2), tot_copy+tot_sp);
		specific_tissue= matrix(0, nrow(ltpm2), tot_copy+tot_sp);
		rownames(specific_score) = rownames(ltpm2);
		colnames(specific_score) = c(sps3, paste0("sum_of_", sps));
		for (i in 1:ncol(specific_score)) {
			a = ltpm2[,((i-1)*nts+1):(i*nts)];
			specific_score[,i] = apply(a,1,max)/(apply(a,1,sum)+thres);
		}
		exp[[type]][['specific_score']] = specific_score;
# entropy (>=0, near 0 is more tissue-specific)
		entropy_score = matrix(0, nrow(ltpm2), tot_copy+tot_sp);
		rownames(entropy_score) = rownames(ltpm2);
		colnames(entropy_score) = c(sps3, paste0("sum_of_", sps));
		for (i in 1:ncol(entropy_score)) {
			a = ltpm2[,((i-1)*nts+1):(i*nts)];
			entropy_score[,i] = apply(a,1,function(x) {entropy(x+0.001,unit='log2')});
		}
		exp[[type]][['entropy_score']] = entropy_score;
# mean, SD
		ltpm2.mean = matrix(0, nrow(ltpm2), tot_copy+tot_sp);
		rownames(ltpm2.mean) = rownames(ltpm2);
		colnames(ltpm2.mean) = c(sps3, paste0("sum_of_", sps));
		for (i in 1:ncol(ltpm2.mean)) {
			a = ltpm2[,((i-1)*nts+1):(i*nts)];
			ltpm2.mean[,i] = apply(a,1,mean);
		}
		ltpm2.sd = matrix(0, nrow(ltpm2), tot_copy+tot_sp);
		rownames(ltpm2.sd) = rownames(ltpm2);
		colnames(ltpm2.sd) = c(sps3, paste0("sum_of_", sps));
		for (i in 1:ncol(ltpm2.sd)) {
			a = ltpm2[,((i-1)*nts+1):(i*nts)];
			ltpm2.sd[,i] = apply(a,1,sd);
		}

		ltpm2.norm1 = ltpm2;
		ltpm2.norm2 = ltpm2;
		for (i in 1:ncol(ltpm2.sd)) {
			sd1 = ltpm2.sd[,i];
			sd1[sd1==0] = 1;
			a = ltpm2[,((i-1)*nts+1):(i*nts)];
			ltpm2.norm1[,((i-1)*nts+1):(i*nts)] = (a-ltpm2.mean[,i]);
			ltpm2.norm2[,((i-1)*nts+1):(i*nts)] = (a-ltpm2.mean[,i])/sd1;
		}

#		norm_ltpm1=ltpm1;
#		for (i in 1:tot_copy) {
#			norm_ltpm1[(1:nrow(norm_ltpm2))*tot_copy-tot_copy+i,] = norm_ltpm2[,nts*(i-1)+(1:nts)];
#		}
		exp[[type]][['lraw2.mean']] = ltpm2.mean;
		exp[[type]][['lraw2.sd']] = ltpm2.sd;
		exp[[type]][['lraw2.center_norm']] = ltpm2.norm1;
		exp[[type]][['lraw2.zscore_norm']] = ltpm2.norm2;
	}
# }}}
###########################

###########################
# compute score of ortholog/paralog-pairs
###########################
# {{{
	for (type in c('tpm', 'fpkm', 'count')) {
		ltpm2 = exp[[type]][['lraw2']];
		thres = exp[[type]][['thres']];
# correlation
		m=nrow(ltpm2);
		n=ncol(ltpm2)/nts;
		ltpm2.pair_cor=matrix(0, m, n*n);
		for (i in 1:n) {
			for (j in 1:n) {
				if (i==j) {next;}
				ltpm2.pair_cor[,(i-1)*n+j] = apply(ltpm2[,c( (nts*(i-1)+1):(nts*i), (nts*(j-1)+1):(nts*j))], 1, function(x) {cor(x[1:nts],x[(nts+1):(nts*2)])});
				f=is.na(ltpm2.pair_cor[,(i-1)*n+j]);
				ltpm2.pair_cor[f,(i-1)*n+j] = 0;
			}
		}
		exp[[type]][['lraw2.pair_cor']] = ltpm2.pair_cor;

		ltpm2.norm2.pair_cor=matrix(0, m, n*n);
		for (i in 1:n) {
			for (j in 1:n) {
				if (i==j) {next;}
				ltpm2.norm2.pair_cor[,(i-1)*n+j] = apply(ltpm2.norm2[,c( (nts*(i-1)+1):(nts*i), (nts*(j-1)+1):(nts*j))], 1, function(x) {cor(x[1:nts],x[(nts+1):(nts*2)])});
				f=is.na(ltpm2.norm2.pair_cor[,(i-1)*n+j]);
				ltpm2.norm2.pair_cor[f,(i-1)*n+j] = 0;
			}
		}
		exp[[type]][['lraw2.zscore_norm.pair_cor']] = ltpm2.norm2.pair_cor;
# t.test
		m=nrow(ltpm2);
		n=ncol(ltpm2)/nts;
		ltpm2.tp=matrix(0, m, n*n);
		for (i in 1:n) {
			for (j in 1:n) {
				if (i==j) {next;}
				ltpm2.tp[,(i-1)*n+j] = apply(ltpm2[,c( (nts*(i-1)+1):(nts*i), (nts*(j-1)+1):(nts*j))], 1, function(x) {t.test(x[1:nts],x[(nts+1):(nts*2)])$p.value});
				f=is.na(ltpm2.tp[,(i-1)*n+j]);
				ltpm2.tp[f,(i-1)*n+j] = 1.0;
			}
		}
		exp[[type]][['lraw2.pair_tp']] = ltpm2.tp;
	}
# }}}


	save(exp, cluster_is_exp_bad, cid_lists, file=paste0(wd,"/Rdata/tpm.Rdata"));
}
# }}}
###########################


###########################
# gene cDNA sequence similarity
###########################
#if (file.exists('Rdata/pairs.Rdata')) {
#	load('Rdata/pairs.Rdata');
#} else {
#	a = edge = read.table('input/in3.edge.gz', head=F);
#	rownames(a) = paste0(a[,1],'+',a[,2]);
#	b1=rep(sps2,each=5);
#	b2=rep(sps2,5);
#	tmp=c();
#	for (i in 1:tot_copy) { tmp=cbind(tmp,gene_ids2); }
#	gene_pairs=cbind(paste0(rep(b1,ncluster1),'|',rep(gene_ids1,each=5)), paste0(rep(b2,ncluster1), '|', as.vector(t(tmp)) ) );
#	gene_pairs = merge(gene_pairs, a[,c(1:2,4:5)], by=1:2, all.x=T, all.y=F, sort=F); 
#	rownames(gene_pairs) = paste0(gene_pairs[,1],'+', gene_pairs[,2]);
#	colnames(gene_pairs)=c('id1', 'id2', 'score', 'iden_perc');
#	save(gene_pairs, file='Rdata/pairs.Rdata');
#}

# initialize data (lfpkm, gene_ids2, gene_names2)
# {{{
ZF_idx  = 1:nts;
CC_idx1 = (nts+1):(2*nts);
CC_idx2 = (2*nts+1):(3*nts);
GF_idx1 = (3*nts+1):(4*nts);
GF_idx2 = (4*nts+1):(5*nts);
ZF_i = 1;
CC_i1 = 2; CC_i2 = 3;
GF_i1 = 4; GF_i2 = 5;
thres = 1.0;
thres1= 1.0;
thres2= 1.0;
lthres = log2(thres+1);
lthres1 = log2(thres1+1);
lthres2 = log2(thres2+1);
data =list();
thres = data$thres = 1;
data$fpkm = exp$fpkm$raw2[,1:(tot_copy*nts)];
data$gene_ids2 = exp$fpkm$gene_ids2;
data$gene_names2 = data$gene_ids2;
for (i in 1:tot_copy) {
	a=matrix(data$gene_ids2[,i],nrow(data$fpkm),1);
	sp = sps2[i];
	data$gene_names2[,i] = merge(a, anchor[ anchor[,'species']==sp, c('anchor_gid', 'anchor_name') ], by=1, sort=F)[,2];
}
# exchange ohno-1 and ohno-2 if the sum of both correlation is larger after exchange
tmp_cor1=pair_cor(data$fpkm[,CC_idx1],data$fpkm[,GF_idx1]) + pair_cor(data$fpkm[,CC_idx2],data$fpkm[,GF_idx2]);
tmp_cor2=pair_cor(data$fpkm[,CC_idx1],data$fpkm[,GF_idx2]) + pair_cor(data$fpkm[,CC_idx2],data$fpkm[,GF_idx1]);
f = tmp_cor2>tmp_cor1;
data$fpkm[f,c(GF_idx1,GF_idx2)] = data$fpkm[f,c(GF_idx2,GF_idx1)];
data$gene_ids2[f,4:5] = data$gene_ids2[f,5:4];
data$gene_names2[f,4:5] = data$gene_names2[f,5:4];
# exchange ohno-1 and ohno-2 if 
#   max(cor(ZF,CC-ohno-1),cor(ZF,GF-ohno-1)) is lower than max(cor(ZF,CC-ohno-2),cor(ZF,GF-ohno-2)) 
tmp_cor11=pair_cor(data$fpkm[,ZF_idx],data$fpkm[,CC_idx1]);
tmp_cor12=pair_cor(data$fpkm[,ZF_idx],data$fpkm[,GF_idx1]);
tmp_cor1 = apply(cbind(tmp_cor11,tmp_cor12), 1, max);
tmp_cor21=pair_cor(data$fpkm[,ZF_idx],data$fpkm[,CC_idx2]);
tmp_cor22=pair_cor(data$fpkm[,ZF_idx],data$fpkm[,GF_idx2]);
tmp_cor2 = apply(cbind(tmp_cor21,tmp_cor22), 1, max);
f = tmp_cor2>tmp_cor1;
data$fpkm[f,c(CC_idx1,CC_idx2)] = data$fpkm[f,c(CC_idx2,CC_idx1)];
data$fpkm[f,c(GF_idx1,GF_idx2)] = data$fpkm[f,c(GF_idx2,GF_idx1)];
data$gene_ids2[f,2:3] = data$gene_ids2[f,3:2];
data$gene_ids2[f,4:5] = data$gene_ids2[f,5:4];
data$gene_names2[f,2:3] = data$gene_names2[f,3:2];
data$gene_names2[f,4:5] = data$gene_names2[f,5:4];

tmp_sum = ohno_sum(data$fpkm, tissues, sps2);
data$fpkm = cbind(data$fpkm, tmp_sum);
data$lfpkm = log2(data$fpkm+1);

data$gene_Interpros2 = list();
for (i in 1:tot_copy) {
	sp = sps2[i];
	data$gene_Interpros2[[i]] = id_to_gs$Interpro[[sp]][ data$gene_ids2[,i] ];
	names(data$gene_Interpros2[[i]]) = data$gene_ids2[,i] ;
}
# }}}

###########################
# load gene pair identity, bit, ...
###########################
# {{{

data$n=nrow(data$lfpkm);
data$prot_iden = matrix(0, data$n, tot_copy*tot_copy);
data$prot_bit  = matrix(0, data$n, tot_copy*tot_copy);
data$prot_cov  = matrix(0, data$n, tot_copy*tot_copy);
for (i in 1:tot_copy) {
	for (j in 1:tot_copy) {
		sp1 = sps2[i];
		sp2 = sps2[j];
		a = merge( cbind(paste0(sp1, '|', data$gene_ids2[,i]), paste0(sp2, '|', data$gene_ids2[,j])),   edge$prot, by=1:2, all.x=T, all.y=F, sort=F);
		a = merge( cbind(paste0(sp1, '|', data$gene_ids2[,i]), paste0(sp2, '|', data$gene_ids2[,j])), a, by=1:2, sort=F);
		data$prot_iden[,tot_copy*(i-1)+j] = a[,3];
		data$prot_bit [,tot_copy*(i-1)+j] = a[,12];
		if (sp1=='ZF') {
			data$prot_cov [,tot_copy*(i-1)+j] = a[,19];
		} else {
			data$prot_cov [,tot_copy*(i-1)+j] = apply(a[,18:19],1,min);
		}
	}
}
data$prot_bit [which(is.na(data$prot_bit),arr=T) ] = 0;
data$prot_iden[which(is.na(data$prot_iden),arr=T)] = 0;
data$prot_cov [which(is.na(data$prot_cov),arr=T) ] = 0;
for (i in 1:tot_copy) {   for (j in 1:tot_copy) {
		aij=data$prot_iden[,tot_copy*(i-1)+j];
		aji=data$prot_iden[,tot_copy*(j-1)+i];
		data$prot_iden[,tot_copy*(i-1)+j] = apply(cbind(aij,aji),1,max);
		data$prot_iden[,tot_copy*(j-1)+i] = data$prot_iden[,tot_copy*(i-1)+j];
} }
for (i in 1:tot_copy) {   for (j in 1:tot_copy) {
		aij=data$prot_bit[,tot_copy*(i-1)+j];
		aji=data$prot_bit[,tot_copy*(j-1)+i];
		data$prot_bit[,tot_copy*(i-1)+j] = apply(cbind(aij,aji),1,max);
		data$prot_bit[,tot_copy*(j-1)+i] = data$prot_bit[,tot_copy*(i-1)+j];
} }
for (i in 1:tot_copy) {   for (j in 1:tot_copy) {
		aij=data$prot_cov[,tot_copy*(i-1)+j];
		aji=data$prot_cov[,tot_copy*(j-1)+i];
		data$prot_cov[,tot_copy*(i-1)+j] = apply(cbind(aij,aji),1,max);
		data$prot_cov[,tot_copy*(j-1)+i] = data$prot_cov[,tot_copy*(i-1)+j];
} }

data$nucl_iden = matrix(0, data$n, tot_copy*tot_copy);
data$nucl_bit  = matrix(0, data$n, tot_copy*tot_copy);
data$nucl_cov  = matrix(0, data$n, tot_copy*tot_copy);
for (i in 1:tot_copy) {
	for (j in 1:tot_copy) {
		sp1 = sps2[i];
		sp2 = sps2[j];
		a = merge(  cbind(paste0(sp1, '|', data$gene_ids2[,i]), paste0(sp2, '|', data$gene_ids2[,j])),   edge$nucl, by=1:2, all.x=T, all.y=F, sort=F);
		a = merge( cbind(paste0(sp1, '|', data$gene_ids2[,i]), paste0(sp2, '|', data$gene_ids2[,j])), a, by=1:2, sort=F);
		data$nucl_iden[,tot_copy*(i-1)+j] = a[,3];
		data$nucl_bit [,tot_copy*(i-1)+j] = a[,12];
		data$nucl_cov [,tot_copy*(i-1)+j] = apply(a[,18:19],1,min);
	}
}
data$nucl_bit[which(is.na(data$nucl_bit),arr=T)] = 0;
data$nucl_iden[which(is.na(data$nucl_iden),arr=T)] = 0;
data$nucl_cov[which(is.na(data$nucl_cov),arr=T)] = 0;
for (i in 1:tot_copy) {   for (j in 1:tot_copy) {
		aij=data$nucl_iden[,tot_copy*(i-1)+j];
		aji=data$nucl_iden[,tot_copy*(j-1)+i];
		data$nucl_iden[,tot_copy*(i-1)+j] = apply(cbind(aij,aji),1,max);
		data$nucl_iden[,tot_copy*(j-1)+i] = data$nucl_iden[,tot_copy*(i-1)+j];
} }
for (i in 1:tot_copy) {   for (j in 1:tot_copy) {
		aij=data$nucl_bit[,tot_copy*(i-1)+j];
		aji=data$nucl_bit[,tot_copy*(j-1)+i];
		data$nucl_bit[,tot_copy*(i-1)+j] = apply(cbind(aij,aji),1,max);
		data$nucl_bit[,tot_copy*(j-1)+i] = data$nucl_bit[,tot_copy*(i-1)+j];
} }
for (i in 1:tot_copy) {   for (j in 1:tot_copy) {
		aij=data$nucl_cov[,tot_copy*(i-1)+j];
		aji=data$nucl_cov[,tot_copy*(j-1)+i];
		data$nucl_cov[,tot_copy*(i-1)+j] = apply(cbind(aij,aji),1,max);
		data$nucl_cov[,tot_copy*(j-1)+i] = data$nucl_cov[,tot_copy*(i-1)+j];
} }
save(data, file='Rdata/data.Rdata');
# }}}

##############################################
# compute score between ortholog-ohnolog or ohnolog-ohnolog
##############################################
# {{{
m=nrow(data$lfpkm);
n=ncol(data$lfpkm)/nts;
data$median2 = matrix(0, m, tot_copy+tot_sp);
data$mean2   = matrix(0, m, tot_copy+tot_sp);
data$sd2     = matrix(0, m, tot_copy+tot_sp);
data$max2    = matrix(0, m, tot_copy+tot_sp);
for (i in 1:n) {
	a = data$lfpkm[, (nts*(i-1)+1):(nts*i)];
	data$median2[,i] = apply(a, 1, median);
	data$mean2[,i]   = apply(a, 1, mean);
	data$sd2[,i]     = apply(a, 1, sd);
	data$max2[,i]    = apply(a, 1, max);
}
data$norm2 = data$lfpkm;
data$zscore2 = data$lfpkm;
for (i in 1:n) {
	data$norm2[,(nts*(i-1)+1):(nts*i)] = data$lfpkm[,(nts*(i-1)+1):(nts*i)] - data$mean2[,i];
	sd1 = data$sd2[,i];
	sd1[is.na(sd1)] = 1;
	data$zscore2[,(nts*(i-1)+1):(nts*i)] = data$norm2[,(nts*(i-1)+1):(nts*i)] / sd1;
}
data$median1 = apply(data$lfpkm[,1:(nts*tot_copy)], 1, median);
data$mean1 = apply(data$lfpkm[,1:(nts*tot_copy)], 1, mean);
data$sd1   = apply(data$lfpkm[,1:(nts*tot_copy)], 1, sd);
data$zscore1 = (data$lfpkm-data$mean1)/data$sd1;

# count number of expressed samples
data$exp_sm_count <- matrix(0, nrow(data$lfpkm), tot_copy+tot_sp);
rownames(data$exp_sm_count) = rownames(data$lfpkm);
colnames(data$exp_sm_count) = c(sps3, paste0("sum_of_", sps));
for (i in 1:ncol(data$exp_sm_count)) {
	a=data$lfpkm[, ((nts)*(i-1)+1):(nts*i)];
	data$exp_sm_count[,i] <- apply(a>=thres, 1, sum);
}
# tissue-specific score (<=1, near 1 is single tissue-specific)
data$specific_score = matrix(0, nrow(data$lfpkm), tot_copy+tot_sp);
data$specific_tissue = matrix(0, nrow(data$lfpkm), tot_copy+tot_sp);
rownames(data$specific_score) = rownames(data$lfpkm);
colnames(data$specific_score) = c(sps3, paste0("sum_of_", sps));
for (i in 1:ncol(data$specific_score)) {
	a = data$lfpkm[,((i-1)*nts+1):(i*nts)];
	data$specific_score[,i] = apply(a,1,max)/(apply(a,1,sum)+thres/10);
	f = data$specific_score[,i]>=0.4;
	data$specific_tissue[f,i] = max.col(a[f,]);

}
# entropy (>=0, near 0 is more tissue-specific)
data$entropy_score = matrix(0, nrow(data$lfpkm), tot_copy+tot_sp);
rownames(data$entropy_score) = rownames(data$lfpkm);
colnames(data$entropy_score) = c(sps3, paste0("sum_of_", sps));
for (i in 1:ncol(data$entropy_score)) {
	a = data$lfpkm[,((i-1)*nts+1):(i*nts)];
	data$entropy_score[,i] = apply(a,1,function(x) {entropy(x+thres/100,unit='log2')});
}
# gene pair correlation
n=ncol(data$lfpkm)/nts;
data$pair_cor=matrix(0, data$n, n*n);
data$pair_cor_p=matrix(0, data$n, n*n);
for (i in 1:n) {
	for (j in 1:n) {
		if (i==j) {next;}
		data$pair_cor[,(i-1)*n+j] = pair_cor(data$lfpkm[,(nts*(i-1)+1):(nts*i)], data$lfpkm[,(nts*(j-1)+1):(nts*j)]);
		f=is.na(data$pair_cor[,(i-1)*n+j]);
		data$pair_cor[f,(i-1)*n+j] = 0;
		data$pair_cor_p[,(i-1)*n+j] = pair_cor_p(data$lfpkm[,(nts*(i-1)+1):(nts*i)], data$lfpkm[,(nts*(j-1)+1):(nts*j)]);
	}
}
# gene pair euclidean distance
n=ncol(data$lfpkm)/nts;
data$pair_euc_dist=matrix(0, data$n, n*n);
for (i in 1:n) {
	for (j in 1:n) {
		if (i==j) {next;}
		data$pair_euc_dist[,(i-1)*n+j] = pair_euc_dist(data$lfpkm[,(nts*(i-1)+1):(nts*i)], data$lfpkm[,(nts*(j-1)+1):(nts*j)]);
	}
}
# t.test, paired
n=ncol(data$lfpkm)/nts;
data$tp=matrix(0, m, n*n);
for (i in 1:n) {
	for (j in 1:n) {
		if (i==j) {next;}
		data$tp[,(i-1)*n+j] = pair_t_p(data$lfpkm[,(nts*(i-1)+1):(nts*i)], data$lfpkm[,(nts*(j-1)+1):(nts*j)]);
		f=is.na(data$tp[,(i-1)*n+j]);
		data$tp[f,(i-1)*n+j] = 1.0;
	}
}
# gene pair maximal distance
data$pair_max_dist=matrix(0, m, n*n);
for (i in 1:n) {
	for (j in 1:n) {
		if (i==j) {next;}
		data$pair_max_dist[,(i-1)*n+j] = apply(data$lfpkm[,c( (nts*(i-1)+1):(nts*i), (nts*(j-1)+1):(nts*j))], 1, function(x) {max(abs(x[1:nts]-x[(nts+1):(nts*2)]))});
		data$pair_max_dist[f,(i-1)*n+j] = 0;
	}
}

# number of co-expressed tissue
m=data$n;
n=ncol(data$lfpkm)/nts;
x=matrix(0,m,n*n);
y=matrix(0,m,n*n);
for (i in 1:n) {
	for (j in 1:n) {
		if (i==j) {next;}
		a1=data$fpkm[,(nts*(i-1)+1):(nts*i)];
		a2=data$fpkm[,(nts*(j-1)+1):(nts*j)];
		a = a1>=thres & (a2>=thres | a1-a2>thres*0.5) | a2>=thres & (a1>=thres | a2-a1>thres*0.5);
		x[,(i-1)*n+j] = apply(a, 1, sum);
		a = a1<thres & a2<thres & abs(a2-a1)<0.5*thres;
		y[,(i-1)*n+j] = apply(a, 1, sum);
	}
}
data$pair_coexp_tissue_num = x;
data$pair_cosil_tissue_num = y;

m=data$n;
n=tot_copy;
x=matrix(0,m,n*n);
y=matrix(0,m,n*n);
rownames(x) = rownames(data$lfpkm);
rownames(y) = rownames(data$lfpkm);
for (i in 1:(n-1)) {
	for (j in (i+1):n) {
		a1 = data$gene_Interpros2[[i]];
		n1 = sapply(a1, length);
		a2 = data$gene_Interpros2[[j]];
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
data$pair_interpor_shared_num = x;
data$pair_interpor_jaccard = y;


save(data, file='Rdata/data.Rdata');
#}}}

# count ohnolog cluster for each fate
# {{{
data$fate = list(
	CC = fate_classify(data$fpkm, 1,2,3, thres, thres, thres2, fc=2.0, C1=0.6, C3=0.75, cor_p=0.1, tp=0.1),
	GF = fate_classify(data$fpkm, 1,4,5, thres, thres, thres2, fc=2.0, C1=0.6, C3=0.75, cor_p=0.1, tp=0.1)
	);

data$fate_n = list();
for (sp in c("CC", "GF")) {
for (ft in names(data$fate[[sp]])) {
	data$fate_n[[sp]][[ft]] = sum(data$fate[[sp]][[ft]]);
}}
data$neo_score = cbind(
		CC=neofunc_score(data$lfpkm, 1, 2, 3, 7),
		GF=neofunc_score(data$lfpkm, 1, 4, 5, 8)
		);
data$sub_score = cbind(
		CC=subfunc_score(data$lfpkm, 1, 2, 3, 7),
		GF=subfunc_score(data$lfpkm, 1, 4, 5, 8)
		);
x=matrix(0, length(data$fate_n$CC), 4);
rownames(x) = names(data$fate_n$CC);
for (ft in names(data$fate_n$CC)) {
	x[ft,1] = data$fate_n$CC[ft];
	x[ft,2] = data$fate_n$GF[ft];
	x[ft,3] = sum(data$fate$CC[[ft]] | data$fate$GF[[ft]]);
	x[ft,4] = sum(data$fate$CC[[ft]] & data$fate$GF[[ft]]);
}
a = rep(0,4);
a[1] = sum(data$fate$CC$double_conserved | data$fate$CC$dosage_balance);
a[2] = sum(data$fate$GF$double_conserved | data$fate$GF$dosage_balance);
a[3] = sum( (data$fate$CC$double_conserved | data$fate$CC$dosage_balance) | (data$fate$GF$double_conserved | data$fate$GF$dosage_balance) );
a[4] = sum( (data$fate$CC$double_conserved | data$fate$CC$dosage_balance) & (data$fate$GF$double_conserved | data$fate$GF$dosage_balance) );
x = rbind(x, conserved_or_balance=a);
x
x/data$n;

# }}}

# bootstrap 1000 1
# {{{
lfpkm1000 = list();
for (k in 1:1000) {
	lfpkm1000[[k]] = random_data(data$lfpkm[,1:(nts*tot_copy)], replace=F);
	sum_tpm2 = matrix(0, data$n, tot_sp*nts);
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

save(data, file='Rdata/data.Rdata');

# plot correlation between ZF, CC and GF
# {{{
tmp_cor_mean = c();
tmp_cor_sd = c();
tmp_cor_se = c();
x=x1=c(data$pair_cor[,2], data$pair_cor[,3]);
# CC
tmp_cor_mean['ZF~CC'] = mean(x);
tmp_cor_sd['ZF~CC'] = sd(x);
tmp_cor_se['ZF~CC'] = sd(x)/sqrt(length(x));
x=x2=data$pair_cor[,7];
tmp_cor_mean['ZF~CC_pair'] = mean(x);
tmp_cor_sd['ZF~CC_pair'] = sd(x);
tmp_cor_se['ZF~CC_pair'] = sd(x)/sqrt(length(x));
x=x3=data$pair_cor[,8*(2-1)+3];
tmp_cor_mean['CC1~CC2'] = mean(x);
tmp_cor_sd['CC1~CC2'] = sd(x);
tmp_cor_se['CC1~cc2'] = sd(x)/sqrt(length(x));
p1=wilcox.test(x1,x2)$p.value; # 3.226e-08
p2=wilcox.test(x2,x3)$p.value; # < 2.2e-16
# GF
x=x1=c(data$pair_cor[,4], data$pair_cor[,5]);
tmp_cor_mean['ZF~GF'] = mean(x);
tmp_cor_sd['ZF~GF'] = sd(x);
tmp_cor_se['ZF~GF'] = sd(x)/sqrt(length(x));
x=x2=data$pair_cor[,8];
tmp_cor_mean['ZF~GF_pair'] = mean(x);
tmp_cor_sd['ZF~GF_pair'] = sd(x);
tmp_cor_se['ZF~GF_pair'] = sd(x)/sqrt(length(x));
x=x3=data$pair_cor[,8*(4-1)+5];
tmp_cor_mean['GF1~GF2'] = mean(x);
tmp_cor_sd['GF1~GF2'] = sd(x);
tmp_cor_se['GF1~cc2'] = sd(x)/sqrt(length(x));
p3=wilcox.test(x1,x2)$p.value; # 1.986e-15
p4=wilcox.test(x2,x3)$p.value; # < 2.2e-16

res=300;
jpeg("output/plot/cor_bar_se.jpg", res=res, w=res*9, h=res*6);
par(mar=c(3.5,5,5,0.5)+0.1);
par(cex=1.2, cex.lab=1.2, cex.axis=1.1);
barplot_sd1(tmp_cor_mean, tmp_cor_se, width=0.9, space=c(0.111,0.111,0.111,1.222,0.111,0.111), col=color[1:3], ylim=c(0.4,0.8), border=color[1:3], xpd=F, ylab="Average Correlation", names.arg=F);
cx=par('cxy')[1];
cy=par('cxy')[2];
text(c(1:3,5:7)*1-0.48, 0.4-cy, names(tmp_cor_mean), xpd=T, cex=0.9);
axis(2);
y0=tmp_cor_mean[1]+tmp_cor_se[1]+cy*0.5;
y1=tmp_cor_mean[2]+tmp_cor_se[2]+cy*0.5;
y2=max(y0,y1)+cy;
lines(c(0.55,0.55), c(y0,y2));
lines(c(0.55,1.55), c(y2,y2));
lines(c(1.55,1.55), c(y1,y2));
text(1, y2+cy, sprintf("%.2e",p1));
y0=tmp_cor_mean[2]+tmp_cor_se[2]+cy*0.5;
y1=tmp_cor_mean[3]+tmp_cor_se[3]+cy*0.5;
y2=max(y0,y1)+cy;
y1=y0=max(y0,y1);
lines(c(1.55,1.55), c(y0,y2), xpd=T);
lines(c(1.55,2.55), c(y2,y2), xpd=T);
lines(c(2.55,2.55), c(y1,y2), xpd=T);
text(2, y2+cy, sprintf("%.2e",p2), xpd=T);

y0=tmp_cor_mean[4]+tmp_cor_se[4]+cy*0.5;
y1=tmp_cor_mean[5]+tmp_cor_se[5]+cy*0.5;
y2=max(y0,y1)+cy;
lines(c(4.55,4.55), c(y0,y2), xpd=T);
lines(c(4.55,5.55), c(y2,y2), xpd=T);
lines(c(5.55,5.55), c(y1,y2), xpd=T);
text(5, y2+cy, sprintf("%.2e",p3), xpd=T);
y0=tmp_cor_mean[5]+tmp_cor_se[5]+cy*0.5;
y1=tmp_cor_mean[6]+tmp_cor_se[6]+cy*0.5;
y2=max(y0,y1)+cy;
y1=y0=max(y0,y1);
lines(c(5.55,5.55), c(y0,y2), xpd=T);
lines(c(5.55,6.55), c(y2,y2), xpd=T);
lines(c(6.55,6.55), c(y1,y2), xpd=T);
text(6, y2+cy, sprintf("%.2e",p4), xpd=T);
dev.off();
# }}}

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
# draw 2D-histogram between:
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
tmp_fate_names=c('double_conserved', 'dosage_balance', 'single_conserved', 'diverged', 'dosage_imbalance', 'subfunc','neofunc','nonfunc')
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
##########################################
# expression ~ fate
##########################################


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

data$sm_cor_dist = cor_dist(t(data$lfpkm[,1:30]));
data$sm_h1 = hclust(data$sm_cor_dist, method="ward.D2");
data$pca1 = princomp(data$lfpkm[,1:30], scores=T); # 4: 78.8% variance, 9: 90.18% variance
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

# random select 500 genes for plot
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



