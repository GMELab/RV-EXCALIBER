#!/usr/env/R

#===============================================================================================================================
# Title: sub_generate_RVBurdenMatrix_gnomAD_rvexcaliber.R
#===============================================================================================================================
#    This file is part of rvexcaliber.
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

# Copyright 2021 Ricky Lali, Michael Chong, Arghavan Omidi, Pedrum Mohammadi-Shemirani, Ann Le, Edward Cui, and Guillaume Pare
#===============================================================================================================================




# Script start


# Load packages


if (!"data.table" %in% installed.packages()) {

  stop("Please ensure the 'data.table' R package is installed")

}

library(data.table)

args                          <- commandArgs(TRUE)
outdir                        <- args[1]
internal_study                <- args[2]
MAF_MCAP                      <- args[3]
eth                           <- args[4]
coverage                      <- args[5]
internal_study_RVBurdenMatrix <- args[6]


# Define functions


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


# End of functions


# Generate list of the variant list files


setwd(outdir)

chr <-
  as.list(1:22)

gnomAD_extracted.l <-
  lapply(
    chr,
    function(x)
    paste0(x,"_gnomAD_variant_list_extracted_",MAF_MCAP,"_",eth,"_",coverage,".txt")
  )


# Read in the gnomAD variant lists


gnomAD_extracted <-
  lapply(
    gnomAD_extracted.l,
    read.delim,
    header=
      FALSE,
    sep=
      " "
  )


# Remove rows of variant list with the same gnomAD allele count
# and gene ID


extracted_df_sub <-
  lapply(
    gnomAD_extracted,
    function(x)
    x[!duplicated(x[,c(4,6)]),][,c(4,6)]
  )


# Generate a nested list of CMAC elements (daughter list) per chromosmome
# (parent list)


CMAC <-
  lapply(
    extracted_df_sub,
    function(x)
    as.list(
      as.numeric(
        x[,2]
      )
    )
  )


# Load in sample IDs from the internal rare variant burden matrix

samples <-
  as.character(
    colnames(
      fread(
        file=
          internal_study_RVBurdenMatrix,
        header=
          TRUE
      )
    )
  )[-1]


# Generate a rare variant burden matrix that corresponds to the internal
# testing and/or ranking dataset


CMAC_rep <-
  do.call(
    "cbind",
    lapply(
      CMAC,
      function(x)
      do.call(
        "cbind",
        lapply(
          x,
          function(x)
          rep(
            x,
            length(
              samples
            )
          )
        )
      )
    )
  )

colnames(CMAC_rep) <-
  as.character(
    do.call(
      "rbind.data.frame",
      extracted_df_sub
    )[,1]
  )

rownames(CMAC_rep) <-
  as.character(
    samples
  )


# Write matrix


CMAC_rep_write <-
  convert_to_write_ready(
    CMAC_rep
  )

write.table(
  CMAC_rep_write,
  paste0("gnomAD_RVBurdenMatrix_",internal_study,"_",MAF_MCAP,"_",eth,"_",coverage,".txt"),
  quote=
    FALSE,
  row.names=
    FALSE,
  sep=
    " "
)
