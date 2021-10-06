#!/usr/env/R

#=======================================================================================================================
# Title: sub_combine_RVBurdenMatrix_internal_rvexcaliber.R
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


if (!"data.table" %in% installed.packages()) {

  stop("Please ensure the 'data.table' R package is installed")

}

library(data.table)

args <-
  commandArgs(TRUE)

outdir <-
  args[1]

combined_internal_dataset <-
  args[2]

MAF_MCAP <-
  args[3]

eth <-
  args[4]

coverage <-
  args[5]

internal_RVBurdenMatrix_filelist <-
  args[6]

filter_names_mod <-
  paste(MAF_MCAP, eth, coverage, sep="_")


# Load in gnomAD chromosome files as list elements


internal_RVBurdenMatrix_filelist <-
  read.delim(
    file=
      internal_RVBurdenMatrix_filelist,
    header=
      FALSE
  )


# Split indivual rare variant burden matrices into list elements


internal_RVBurdenMatrix.l <-
  split(
    unique(
      as.character(
        internal_RVBurdenMatrix_filelist[which(grepl(filter_names_mod,internal_RVBurdenMatrix_filelist$V1)),]
      )
    ),
    1:length(
      unique(
        as.character(
          internal_RVBurdenMatrix_filelist[unique(which(grepl(filter_names_mod,as.character(internal_RVBurdenMatrix_filelist$V1)))),]
        )
      )
    )
  )


# Read in the internal rare variant burden matrices


internal_RVBurdenMatrix <-
  lapply(
    internal_RVBurdenMatrix.l,
    function(x)
    fread(
      x,
      header=
        TRUE,
      data.table=
        FALSE
    )
  )


# Write combined matrix


setwd(outdir)

internal_RVBurdenMatrix_combined <-
  do.call(
    "rbind.data.frame",
    internal_RVBurdenMatrix
  )

write.table(
  internal_RVBurdenMatrix_combined,
  paste0(combined_internal_dataset,"_RVBurdenMatrix_",MAF_MCAP,"_",eth,"_",coverage,".txt"),
  quote=
    FALSE,
  row.names=
    FALSE,
  sep=
    " "
)
