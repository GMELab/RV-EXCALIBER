#!/usr/env/R

#=======================================================================================================================
# Title: get_iCF_association_test_rvexcaliber.R
#=======================================================================================================================
#    This file is part of RV-EXCALIBER.
#
#    RV-EXCALIBER is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    RV-EXCALIBER is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with RV-EXCALIBER.  If not, see <https://www.gnu.org/licenses/>.

# Authors: Ricky Lali, Michael Chong, Arghavan Omidi, Pedrum Mohammadi-Shemirani, Ann Le, Edward Cui, and Guillaum Pare
#=======================================================================================================================




# Script start


# Load packages


if (!"ggplot2" %in% installed.packages()) {

  stop("Please ensure the 'ggplot2' R package is installed")

} else if (!"data.table" %in% installed.packages()) {

  stop("Please ensure the 'data.table' R package is installed")

}

suppressMessages(library(ggplot2))
suppressMessages(library(data.table))

args <-
  commandArgs(TRUE)

path_to_rvexcaliber <-
  args[1]

outdir <-
  args[2]

internal_testing_study <-
  args[3]

internal_ranking_study <-
  args[4]

stage <-
  args[5]

MAF_MCAP <-
  args[6]

eth <-
  args[7]

coverage <-
  args[8]

adjust_type <-
  args[9]

testing_RVBurdenMatrix <-
  args[10]

ranking_RVBurdenMatrix <-
  args[11]

gnomAD_RVBurdenMatrix  <-
  args[12]

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

  if (obs_or_exp == "obs") {

      if (stage == "stage0") {

          colID <-
            "rank_obs_count"

      } else if (stage == "stage1" | stage == "stage2") {

          colID <-
            "test_obs_count"

      }


  } else if (obs_or_exp == "gnomAD") {

      if (stage == "stage0") {

          colID <-
            paste0("rank_gnomAD_count_",adjust_type)

      } else if (stage == "stage1" | stage == "stage2") {

        colID <-
          paste0("test_gnomAD_count_",adjust_type)

      }

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


# End of functions


# Load in testing/ranking and gnomAD rare variant burden matrices


setwd(outdir)

default_ranking_SummaryAssociations_df <-
  read.delim(
    file=paste0(path_to_rvexcaliber,"/default_ranking_dataset/rvexcaliber_ranking_MIGen_exome_controls_0.001_0.025_nfe_high_coverage_coding_iCFadjust_rvexcaliber_base.txt"),
    header=
      TRUE,
    sep=
      "\t"
  )

if (stage == "stage0") {

  obs_RVBurdenMatrix <-
    fast_read(
      ranking_RVBurdenMatrix,
      "Gene"
    )

  gnomAD_RVBurdenMatrix <-
    fast_read(
      gnomAD_RVBurdenMatrix,
      "Gene"
    )

} else if (stage == "stage1") {

  obs_RVBurdenMatrix <-
    fast_read(
      testing_RVBurdenMatrix,
      "Gene"
    )

  gnomAD_RVBurdenMatrix <-
    fast_read(
      gnomAD_RVBurdenMatrix,
      "Gene"
    )

} else if (stage == "stage2") {

  obs_RVBurdenMatrix <-
    fast_read(
      testing_RVBurdenMatrix,
      "Gene"
    )

  ranking_SummaryAssociations_df <-
    read.delim(
      file=paste0("rvexcaliber_ranking_",internal_ranking_study,"_SummaryAssociations_",filter_names_txt),
      header=
        TRUE,
      sep=
        "\t"
    )

  gnomAD_RVBurdenMatrix <-
    fast_read(
      gnomAD_RVBurdenMatrix,
      "Gene"
    )

}


# Check dimentionality using the 'nrow' of the gnomAD rare variant burden matrix as the dimention standard
# (i.e. the testing rare variant burden matrix rows should equal rows of the gnomAD rare varient burden matrix)


if (nrow(obs_RVBurdenMatrix) != nrow(gnomAD_RVBurdenMatrix) & ncol(obs_RVBurdenMatrix) != nrow(gnomAD_RVBurdenMatrix)) {

    stop("Number of samples in both the testing and gnomAD rare variant burden matrices should be the same")

}


# Prune out genes that are under-powerered to inform the iCF (i.e. genes with observed or expected counts <1)


if (stage == "stage0") {

    keep_genes <-
        intersect(
          as.character(
            names(
              which(
                apply(
                  obs_RVBurdenMatrix,
                  2,
                  sum
                )>=1
              )
            )
          ),
          as.character(
            names(
              which(
                apply(
                  gnomAD_RVBurdenMatrix,
                  2,
                  sum
                )>=1
              )
            )
          )
        )

} else if (stage == "stage1" | stage == "stage2") {

    if (stage == "stage1") {

    ranking_SummaryAssociations_df <-
      default_ranking_SummaryAssociations_df


    } else if (stage == "stage2") {

      DN <-
        "Do Nothing"

      rm(DN)

    }

    keep_genes <-
      intersect(
        intersect(
          as.character(
            names(
              which(
                apply(
                  obs_RVBurdenMatrix,
                  2,
                  sum
                )>=1
              )
            )
          ),
          as.character(
            names(
              which(
                apply(
                  gnomAD_RVBurdenMatrix,
                  2,
                  sum
                )>=1
              )
            )
          )
        ), as.character(ranking_SummaryAssociations_df$Gene)
      )

} else {

    stop("Internal stop (should not appear)")

}

obs_RVBurdenMatrix_pruned <-
  obs_RVBurdenMatrix[,keep_genes]

gnomAD_RVBurdenMatrix_pruned <-
  gnomAD_RVBurdenMatrix[,keep_genes]


# Check dimensionality of pruned internal testing or ranking rare variant burden matrix
# with the pruned gnomAD rare variant burden matrix


if ((nrow(obs_RVBurdenMatrix_pruned) == nrow(gnomAD_RVBurdenMatrix_pruned)) & (ncol(obs_RVBurdenMatrix_pruned) == ncol(gnomAD_RVBurdenMatrix_pruned))) {

     DN <-
       "Do Nothing"

     rm(DN)

} else {

	   stop("The dimensions of the testing/ranking rare variant burden matrix does not match the gnomAD rare variant burden matrix.")

}


# Calculate per-individual iCF via element-wise matrix division


iCF <-
  as.list(
    apply(
      obs_RVBurdenMatrix_pruned,
      1,
       sum
     ) /
     apply(
       gnomAD_RVBurdenMatrix_pruned,
       1,
       sum
     )
  )

iCF_matrix <-
  do.call(
    'rbind',
    lapply(
      iCF,
      function(x)
      rep(
        x,
        ncol(
          obs_RVBurdenMatrix_pruned
        )
      )
    )
  )

rownames(iCF_matrix) <-
  as.character(
    rownames(
      obs_RVBurdenMatrix_pruned
    )
  )

colnames(iCF_matrix) <-
  rep(
    "iCF",
    ncol(
      obs_RVBurdenMatrix_pruned
    )
  )

if (adjust_type == "noadjust") {

   gnomAD_RVBurdenMatrix_pruned_iCFadjust <-
     gnomAD_RVBurdenMatrix_pruned * 1


} else if (adjust_type == "iCFadjust") {

  # Apply the iCF to the expected count matrix

   gnomAD_RVBurdenMatrix_pruned_iCFadjust <-
     gnomAD_RVBurdenMatrix_pruned * iCF_matrix

}


# Generate the delta count matrix (i.e. element wise matrix substraction)


delta_RVBurdenMatrix <-
  obs_RVBurdenMatrix_pruned - gnomAD_RVBurdenMatrix_pruned_iCFadjust


# Write matrices


if (stage == "stage1" | stage == "stage2") {

    write.table(
      convert_to_write_ready(
        obs_RVBurdenMatrix_pruned
      ),
      paste0("rvexcaliber_testing_",internal_testing_study,"_RVBurdenMatrix_pruned_",filter_names_mod_txt),
      quote=
        FALSE,
      row.names=
        FALSE,
      sep=
        "\t"
    )

    write.table(
      convert_to_write_ready(
        gnomAD_RVBurdenMatrix_pruned_iCFadjust
      ),
      paste0("rvexcaliber_gnomAD_RVBurdenMatrix_pruned_",internal_testing_study,"_",filter_names_txt),
      quote=
        FALSE,
      row.names=
        FALSE,
      sep=
        "\t"
    )

}


# Distribute matrix into a list where each element represents column (i.e. a single gene) from the matrix


p_list <-
  split(
    gnomAD_RVBurdenMatrix_pruned_iCFadjust,
    col(
      gnomAD_RVBurdenMatrix_pruned_iCFadjust
    )
  )

obs_list <-
  split(
    obs_RVBurdenMatrix_pruned,
    col(
      gnomAD_RVBurdenMatrix_pruned_iCFadjust
    )
  )


# Calculate gene-based P-values using rvexcaliber_base


gene_P <-
  cbind.data.frame(
    Gene=
      colnames(
        gnomAD_RVBurdenMatrix_pruned_iCFadjust
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

if (stage == "stage0") {

  names(gene_P)[2] <-
      paste0("rank_rvexcaliber_base_",adjust_type)

} else if (stage == "stage1" | stage == "stage2") {

  names(gene_P)[2] <-
    paste0("test_P_rvexcaliber_base_",adjust_type)

} else {

  stop("Internal stop (should not appear).")

}

obs_allele_count_df <-
  matrix_to_df(
    obs_RVBurdenMatrix_pruned,
    "obs"
  )

gnomAD_allele_count_df <-
  matrix_to_df(
    gnomAD_RVBurdenMatrix_pruned_iCFadjust,
    "gnomAD"
  )


# Merge gene-based observed counts, expected counts, and P-values


rvexcaliber_fin_df <-
  merge(
    merge(
      obs_allele_count_df,
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


# Write Summary Association files


if (stage == "stage0") {

    test_rank <-
      paste0("ranking_",internal_ranking_study)

} else if (stage == "stage1" | stage == "stage2") {

    test_rank <-
      paste0("testing_",internal_testing_study)

}

write.table(
  rvexcaliber_fin_df_P_order,
  paste0("rvexcaliber_",test_rank,"_SummaryAssociations_",filter_names_txt),
  quote=
    FALSE,
  row.names=
    FALSE,
  sep=
    "\t"
)


# FIN
