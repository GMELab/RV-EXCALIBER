#!/usr/env/R

#=======================================================================================================================
# Title: get_gCF_association_test_rvexcaliber.R
#=======================================================================================================================
#    This file is part of RV-EXCALIBER.
#
#    rvexcaliber is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    rvexcaliber is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with rvexcaliber.  If not, see <https://www.gnu.org/licenses/>.

# Authors: Ricky Lali, Michael Chong, Arghavan Omidi, Pedrum Mohammadi-Shemirani, Ann Le, Edward Cui, and Guillaum Pare
#=======================================================================================================================




# Script start


# Load packages


if (!"ggplot2" %in% installed.packages()) {

  stop("Please ensure the 'ggplot2' R package is installed")

} else if (!"dplyr" %in% installed.packages()) {

  stop("Please ensure the 'dplyr' R package is installed")

} else if (!"data.table" %in% installed.packages()) {

  stop("Please ensure the 'data.table' R package is installed")

}

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))

args                                    <- commandArgs(TRUE)
outdir                                  <- args[1]
internal_study                          <- args[2]
MAF_MCAP                                <- args[3]
eth                                     <- args[4]
coverage                                 <- args[5]
adjust_type                             <- args[6]
testing_RVBurdenMatrix_pruned           <- args[7]
gnomAD_RVBurdenMatrix_pruned_iCFadjust  <- args[8]
testing_SummaryAssociations_df          <- args[9]
ranking_SummaryAssociations_df          <- args[10]


# Reset the 'adjust_type'


adjust_type <-
  "iCFgCFadjust"

assoc_mod <-
  "rvexcaliber_base"

filter_names_mod_txt <-
  paste0(MAF_MCAP,"_",eth,"_",coverage,"_",assoc_mod,".txt")

filter_names_txt <-
  paste0(MAF_MCAP,"_",eth,"_",coverage,"_",adjust_type,"_",assoc_mod,".txt")

filter_names_png <-
  paste0(MAF_MCAP,"_",eth,"_",coverage,"_",adjust_type,"_",assoc_mod,".png")


# Define functions


fast_read <- function(mat,colID) {

    read_mat <-
      transpose(
        fread(
          file=
            mat,
          header=
            TRUE
        ),
        keep.names=
          colID,
        make.names=
          colID
      )

    samples <-
      as.character(
        read_mat$Gene
      )

    read_mat_sub <-
      data.matrix(
        read_mat[,2:ncol(read_mat)]
      ); rownames(read_mat_sub) <- samples

    rm(read_mat)
    rm(samples)

    return(read_mat_sub)

}


array_to_df <- function(firstHeader, secondHeader, input) {

    input$firstHeader <-
      rownames(input)

    rownames(input) <-
      c()

    colnames(input) <-
      c(secondHeader, firstHeader)

    return(input)

}


convert_to_write_ready <- function(mat) {

  tmat <-
    t(
      mat
    )

  write_mat <-
    cbind.data.frame(
      Gene=
        rownames(
          tmat
        ),
      data.frame(
        tmat,
        row.names=
          NULL
      )
    )

  return(write_mat)

}


matrix_to_df <- function(input_matrix, obs_or_exp) {

  gene_count_matrix <-
    as.matrix(
      apply(
        input_matrix,
        2,
        sum
      )
    )

  if (obs_or_exp == "test") {

      colID <-
        "test_allele_count"

  } else if (obs_or_exp == "gnomAD") {

      colID <-
        paste0("test_gnomAD_allele_count_",adjust_type)

  }

  colnames(gene_count_matrix) <-
    colID

  Gene <-
    rownames(gene_count_matrix)

  gene_count_df <-
    cbind.data.frame(Gene, gene_count_matrix)

  rownames(gene_count_df) <-
    c()

  return(gene_count_df)

}

rvexcaliber_base <- function(p_list, obs_list) {

    delta_count_mean <-
      mean(
        obs_list - p_list
      )

	  delta_count_sd <-
      sd(
        obs_list - p_list
      )

    z_score <-
      delta_count_mean / (delta_count_sd/(sqrt(length(obs_list))))

	  finP <-
      1-pnorm(z_score)

    return(finP)
}

gco_prop <- function(df) {

  original_gco <-
    nrow(df)

  testing_pdiff <-
    0.55

  new_gco <-
    original_gco*((1-(testing_pdiff/2)))/(1+(testing_pdiff/2))

  # Calculate 99% CI of new_gco proportion using binomial exact test

  z <-
    qnorm(0.995)

  gene_prop <-
    new_gco / original_gco

  lb_prop <-
    gene_prop - (z * sqrt(gene_prop*(1-gene_prop)/original_gco))

  ub_prop <-
    gene_prop + (z * sqrt(gene_prop*(1-gene_prop)/original_gco))

  lb_gco <-
    lb_prop * original_gco

  ub_gco <-
    ub_prop * original_gco

  prop_out <-
    c(new_gco, lb_gco, ub_gco)

  return(prop_out)

}

find_allele_filt <- function(df, lb_prop, ub_prop) {

		fin <-
      0

		allele_filt <-
      0

		while(fin < lb_prop | fin > ub_prop) {

		    df_filt <-
          df[df[,3]>=allele_filt,]

				fin <-
          nrow(df_filt)

				allele_filt <-
          allele_filt + 0.01

		}

    if (allele_filt > 10) {

        allele_filt <-
          10

    }

		return(allele_filt)

}


# Note that the following 'estlambda' function was written in the GenABEL
# package (https://github.com/cran/GenABEL/blob/master/R/estlambda.R)


"estlambda" <-
    function(data,plot = FALSE, proportion=1.0, method="regression", filter=TRUE, df=1,... ) {
	data <- data[which(!is.na(data))]
	if (proportion>1.0 || proportion<=0)
		stop("proportion argument should be greater then zero and less than or equal to one")
	ntp <- round(proportion*length(data))
	if (ntp<1) stop("no valid measurments")
	if (ntp==1) {
		warning(paste("One measurment, Lambda = 1 returned"))
		return(list(estimate=1.0,se=999.99))
	}
	if (ntp<10) warning(paste("number of points is too small:",ntp))
	if (min(data)<0) stop("data argument has values <0")
	if (max(data)<=1) {
#		lt16 <- (data < 1.e-16)
#		if (any(lt16)) {
#			warning(paste("Some probabilities < 1.e-16; set to 1.e-16"))
#			data[lt16] <- 1.e-16
#		}
		data <- qchisq(data,1,lower.tail=FALSE)
	}
	if (filter)
	{
		data[which(abs(data)<1e-8)] <- NA
	}
	data <- sort(data)
	ppoi <- ppoints(data)
	ppoi <- sort(qchisq(ppoi,df=df,lower.tail=FALSE))
	data <- data[1:ntp]
	ppoi <- ppoi[1:ntp]
#	s <- summary(lm(data~offset(ppoi)))$coeff
# 	bug fix thanks to Franz Quehenberger

	out <- list()
	if (method=="regression") {
		s <- summary(lm(data~0+ppoi))$coeff
		out$estimate <- s[1,1]
		out$se <- s[1,2]
	} else if (method=="median") {
		out$estimate <- median(data,na.rm=TRUE)/qchisq(0.5,df)
		out$se <- NA
	} else if (method=="KS") {
		limits <- c(0.5,100)
		out$estimate <- estLambdaKS(data,limits=limits,df=df)
		if ( abs(out$estimate-limits[1])<1e-4 || abs(out$estimate-limits[2])<1e-4 )
			warning("using method='KS' lambda too close to limits, use other method")
		out$se <- NA
	} else {
		stop("'method' should be either 'regression' or 'median'!")
	}

	if (plot) {
		lim <- c(0,max(data,ppoi,na.rm=T))
#		plot(ppoi,data,xlim=lim,ylim=lim,xlab="Expected",ylab="Observed", ...)
		plot(ppoi,data,xlab="Expected",ylab="Observed", ...)
		abline(a=0,b=1)
		abline(a=0,b=out$estimate,col="red")
	}

	out
}


QQplot <- function(Pvals, cqq, cci, titleC) {

    observed <-
      -log10(sort(Pvals,decreasing=F))

    expected <-
      -log10(1:length(observed)/length(observed))

    conf95 <-
      rep(
        0,
        length(
          expected
        )
      )

    conf05 <-
      rep(
        0,
        length(
          expected
        )
      )

    for (i in 1:length(expected)) {

        conf95[i] <-
          qbeta(
            0.95,
            i,
            length(
              expected
            )-i+1
          )

        conf05[i] <-
          qbeta(
            0.05,
            i,
            length(
              expected
            )-i+1
          )

    }

    dfqq <-
      data.frame(
        observed=observed,
        expected=expected,
        conf05=conf05,
        conf95=conf95
      )

    id_max <-
        max(
          -log10(conf05),
          observed
        )

    ggplot(dfqq, aes(x=expected,y=observed)) +
        geom_point(shape=16,size=1.2,colour=cqq) +
        geom_line(aes(x=expected,y=expected), colour=cci) +
        geom_ribbon(aes(ymax=-log10(conf05), ymin=-log10(conf95)), fill=cci, alpha=0.2) +
        ggtitle(titleC) +
        scale_x_continuous(limits = c(0, max(expected))) +
        scale_y_continuous(limits = c(0, id_max)) +
        ylab(expression(Observed~~-log[10](italic(p)))) +
        xlab(expression(Theoretical~~-log[10](italic(p)))) +
        theme_light() +
        theme(
          plot.title=element_text(size=12,hjust=0.5),
          axis.title.y=element_text(size=20),
          axis.title.x=element_text(size=20),
          axis.text.y=element_text(size=18,colour="black"),
          axis.text.x=element_text(size=18,colour="black"),
          panel.background=element_rect(colour="black",size=1),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()
        )
    ggsave(
      paste0("rvexcaliber_testing_",internal_study,"_QQPlot_allele_filter_",filter_names_png)
    )

}


Mountainplot <- function(df, titleC) {

  ggplot(df, aes(x=Gene_index,y=cumulative_distribution,colour=Model)) +
    geom_point(size=0.2,colour="grey50") +
    geom_smooth(method="loess",se=FALSE,lwd=1.5) +
    ggtitle(titleC) +
    xlab("Gene Index") +
    ylab("Cumulative sum\nof delta allele count") +
    theme_light() +
    theme(
      plot.title=element_text(size=18,hjust=0.5),
      axis.title.y=element_text(size=20),
      axis.title.x=element_text(size=20),
      axis.text.y=element_text(size=18,colour="black"),
      axis.text.x=element_text(size=18,colour="black"),
      legend.title=element_text(size=14,colour="black"),
      legend.text=element_text(size=14,colour="black"),
      panel.background=element_rect(colour="black",size=1),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank()
    )
    ggsave(
      paste0("rvexcaliber_testing_",internal_study,"_MountainPlot_",filter_names_png)
    )

}


# End of functions


# Load in testing/ranking and gnomAD rare variant burden matrices

ranking_SummaryAssociations_df <-
  read.delim(
    file=
      ranking_SummaryAssociations_df,
    header=
      TRUE,
    sep=
      "\t"
  )

testing_SummaryAssociations_df <-
  read.delim(
    file=
      testing_SummaryAssociations_df,
    header=
      TRUE,
    sep=
      "\t"
  )

testing_RVBurdenMatrix_pruned <-
  fast_read(
    testing_RVBurdenMatrix_pruned,
    "Gene"
  )

gnomAD_RVBurdenMatrix_pruned_iCFadjust <-
  fast_read(
    gnomAD_RVBurdenMatrix_pruned_iCFadjust,
    "Gene"
  )


# Merge discovery dataframes (1 for P value ordering, the other for calculation of gCF)


rank_test_df <-
  merge(
    ranking_SummaryAssociations_df,
    testing_SummaryAssociations_df,
    by="Gene"
  )


# Set colnames for clarity


colnames(rank_test_df) <-
  c(
    "Gene",
    "rank_allele_count",
    "rank_gnomAD_allele_count_iCFadjust",
    "rank_P_rvexcaliber_base_iCFadjust",
    "test_allele_count",
    "test_gnomAD_allele_count_iCFadjust",
    "test_P_rvexcaliber_base_iCFadjust"
  )


# The merge will alphabetically re-order according the 'merge by field', which
# is "Gene" in the case; so we will re-order by ranking P-value


rank_test_df_P_order <-
  rank_test_df[order(rank_test_df$rank_P_rvexcaliber_base_iCFadjust),] # Can we deprecate? It just gets re-ordered again below


# Assign bin structure according to expected allele count and P-value


rank_test_df$rank_gnomAD_allele_count_iCFadjust_quantile <-
  ntile(
    rank_test_df$rank_gnomAD_allele_count_iCFadjust,
    5
  )


# Variable re-assignment


rank_test_df <-
  do.call(
    'rbind',
    lapply(
      split(
        rank_test_df,
        rank_test_df$rank_gnomAD_allele_count_iCFadjust_quantile
      ),
      function(x) {
        x$rank_P_rvexcaliber_base_iCFadjust_quantile <-
          ntile(
            x[,which(colnames(rank_test_df)=="rank_P_rvexcaliber_base_iCFadjust")],
            10
          )
      return(x)
      }
    )
  ); rownames(rank_test_df) <- c()


# Assign index of gCF


rank_test_df$gCF_index <-
  as.factor(
    paste0(
      rank_test_df$rank_gnomAD_allele_count_iCFadjust_quantile,
      ".",
      rank_test_df$rank_P_rvexcaliber_base_iCFadjust_quantile
    )
  )


# Merge the allele counts for the internal testing dataset and gnomAD based
# on their gCF index

  gCF_merge <-
    merge(
      array_to_df(
        "gCF_index",
        "test_allele_count_sum",
        data.frame(
          tapply(
            rank_test_df$test_allele_count,
            rank_test_df$gCF_index,
            sum
          )
        )
      ),
      test_gnomAD_allele_count_iCFadjust_gCF <-
        array_to_df(
          "gCF_index",
          "test_gnomAD_allele_count_iCFadjust_sum",
          data.frame(
            tapply(
              rank_test_df$test_gnomAD_allele_count_iCFadjust,
              rank_test_df$gCF_index,
              sum
            )
          )
        ),
      by="gCF_index"
    )


# Calculate gCF


gCF_merge$gCF <-
  gCF_merge$test_allele_count_sum / gCF_merge$test_gnomAD_allele_count_iCFadjust_sum


# Merge back to genes based on the gCF index


rank_test_gCF_df <-
  merge(
    rank_test_df,
    gCF_merge[,c(1,4)],
    by="gCF_index"
  )


# Prepare for element-wise matrix multiplication


gCF_matrix <-
  do.call(
    "cbind",
    lapply(
      lapply(
        as.list(
          colnames(
            gnomAD_RVBurdenMatrix_pruned_iCFadjust
          )
        ),
        function(x)
        rank_test_gCF_df[rank_test_gCF_df$Gene %in% x,]$gCF
      ),
      function(y)
      rep(
        y,
        nrow(
          gnomAD_RVBurdenMatrix_pruned_iCFadjust
        )
      )
    )
  )

rownames(gCF_matrix) <-
  as.character(
    rownames(
      gnomAD_RVBurdenMatrix_pruned_iCFadjust
    )
  )

colnames(gCF_matrix) <-
  as.character(
    colnames(
      gnomAD_RVBurdenMatrix_pruned_iCFadjust
    )
  )


# Apply the gCF to the iCF-adjusted gnomAD matrix


gnomAD_RVBurdenMatrix_pruned_iCFgCFadjust <-
  gnomAD_RVBurdenMatrix_pruned_iCFadjust * gCF_matrix


# Generate the delta count matrix (i.e. element wise matrix substraction)


delta_RVBurdenMatrix_pruned_iCFgCFadjust <-
  testing_RVBurdenMatrix_pruned - gnomAD_RVBurdenMatrix_pruned_iCFgCFadjust


# Write matrices


setwd(outdir)

write.table(
  convert_to_write_ready(
    gnomAD_RVBurdenMatrix_pruned_iCFgCFadjust
  ),
  paste0("rvexcaliber_gnomAD_RVBurdenMatrix_pruned_",filter_names_txt),
  quote=
    FALSE,
  row.names=
    FALSE,
  sep=
    "\t"
)


# Distribute matrix into a list where each element represents column (i.e. a single gene) from the matrix


p_list <-
  split(
    gnomAD_RVBurdenMatrix_pruned_iCFgCFadjust,
    col(
      gnomAD_RVBurdenMatrix_pruned_iCFgCFadjust
    )
  )

obs_list <-
  split(
    testing_RVBurdenMatrix_pruned,
    col(
      testing_RVBurdenMatrix_pruned
    )
  )


# Calculate gene-based P-values using rvexcaliber_base


gene_P <-
  cbind.data.frame(
    Gene=
      colnames(
        gnomAD_RVBurdenMatrix_pruned_iCFgCFadjust
      ),
    P_rvexcaliber_base=
      unlist(
        mapply(
          rvexcaliber_base,
          p_list,
          obs_list,
          SIMPLIFY=
            FALSE
        )
      )
  )

names(gene_P)[2] <-
  "test_P_rvexcaliber_base_iCFgCFadjust"


testing_allele_count_df <-
  matrix_to_df(
    testing_RVBurdenMatrix_pruned,
    "test"
  )

gnomAD_allele_count_df <-
  matrix_to_df(
    gnomAD_RVBurdenMatrix_pruned_iCFgCFadjust,
    "gnomAD"
  )


# Merge gene-based observed counts, expected counts, and P-values


rvexcaliber_fin_df <-
  merge(
    merge(
      testing_allele_count_df,
      gnomAD_allele_count_df,
      by=
        "Gene"
    ),
    gene_P,
    by=
      "Gene"
  )


# Order by significant gene-based P-value


rvexcaliber_fin_df_P_order <-
  rvexcaliber_fin_df[order(rvexcaliber_fin_df[,4]),]

allele_filt_prop_stats <-
  gco_prop(
    rvexcaliber_fin_df_P_order
  )

allele_filt <-
  find_allele_filt(
    rvexcaliber_fin_df_P_order,
    allele_filt_prop_stats[2],
    allele_filt_prop_stats[3]
  )

rvexcaliber_fin_df_P_order_allele_filt <-
  rvexcaliber_fin_df_P_order[rvexcaliber_fin_df_P_order[,3]>=allele_filt,]


# Write Summary Association files


write.table(
  rvexcaliber_fin_df_P_order,
  paste0("rvexcaliber_testing_",internal_study,"_SummaryAssociations_",filter_names_txt),
  quote=
    FALSE,
  row.names=
    FALSE,
  sep=
    "\t"
)

write.table(
  rvexcaliber_fin_df_P_order_allele_filt,
  paste0("rvexcaliber_testing_",internal_study,"_SummaryAssociations_allele_filter_",filter_names_txt),
  quote=
    FALSE,
  row.names=
    FALSE,
  sep=
    "\t"
)


# Calculate Genomin Inflation at the median of gene-based P-values


GenomicInflationFactor_median <-
  estlambda(
    data=
      rvexcaliber_fin_df_P_order_allele_filt$test_P_rvexcaliber_base_iCFgCFadjust,
    plot=
      FALSE,
    proportion=
      1,
    method=
      "median",
    filter=
      TRUE,
    df=
      1
  )$estimate

lambda <-
  data.frame(
    assoc_mod=
      assoc_mod,
    adjust_type=
      adjust_type,
    allele_filt=
      allele_filt,
    GenomicInfationFactor_median=
      GenomicInflationFactor_median
  )


# Write Genomic Inflation file


write.table(
  lambda,
  paste0("rvexcaliber_testing_",internal_study,"_GenomicInflationFactor_allele_filter_",filter_names_txt),
  quote=
    FALSE,
  row.names=
    FALSE,
  sep=
    "\t"
)


# Generate quantile-quantile plot to visualize gene-based P-value distribution


nb <-
  nrow(
    rvexcaliber_fin_df_P_order_allele_filt
  )

suppressMessages(
  QQplot(
    rvexcaliber_fin_df_P_order_allele_filt[,4],
    "black",
    "red",
    paste0("internal_study = ", internal_study, "; association model = ",assoc_mod,"\n","adjustment type = ",adjust_type, "; gnomAD count threshold = ",round(allele_filt,3),"\n","GenomicInfationFactor_median = ",round(GenomicInflationFactor_median,2),"; MAF_MCAP = ", MAF_MCAP,"; n genes = ", nb)
  )
)

# Generate mountain plots for 'iCFadjust' and 'iCFgCFadjust' models


mountain <-
  merge(
    rank_test_df_P_order,
    rvexcaliber_fin_df_P_order,
    by="Gene"
  )

mountain_cumulative_sums <-
  rbind(
    data.frame(
      Gene_index=
        1:nrow(mountain),
      cumulative_distribution=
        cumsum(
          mountain[order(mountain$rank_P_rvexcaliber_base_iCFadjust),]$test_allele_count.y - mountain[order(mountain$rank_P_rvexcaliber_base_iCFadjust),]$test_gnomAD_allele_count_iCFadjust
        ),
      Model=
        rep(
          "iCF",
          nrow(
            mountain
          )
        )
    ),
    data.frame(
      Gene_index=
        1:nrow(mountain),
      cumulative_distribution=
        cumsum(
          mountain[order(mountain$rank_P_rvexcaliber_base_iCFadjust),]$test_allele_count.y - mountain[order(mountain$rank_P_rvexcaliber_base_iCFadjust),]$test_gnomAD_allele_count_iCFgCFadjust
        ),
      Model=
        rep(
          "iCF and gCF",
          nrow(
            mountain
          )
        )
    )
  )

suppressMessages(
  Mountainplot(
    mountain_cumulative_sums,
    paste0("internal_study = ", internal_study,"\n","MAF_MCAP = ", MAF_MCAP)
  )
)

# FIN
