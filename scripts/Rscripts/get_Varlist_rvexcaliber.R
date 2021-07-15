#!/usr/env/R

#=======================================================================================================================
# Title: get_Varlist_rvexcaliber.R
# Generate gene-based assoc_modiation statistics based on iCF-adjusted expected allele counts in gnomAD
#=======================================================================================================================
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

# Authors: Ricky Lali, Michael Chong, Arghavan Omidi, Pedrum Mohammadi-Shemirani, Ann Le, Edward Cui, and Guillaum Pare
#=======================================================================================================================




# Script start


# Load packages


if (!"dplyr" %in% installed.packages()) {

  stop("Please ensure the 'dplyr' R package is installed")

}

suppressMessages(library(dplyr))


# Define functions


getethbreak <- function(ethf, dfMAF2f=dfMAF2) {

    maineth <-
      unlist(
        strsplit(
          ethf,
          ""
        )
      )

    ind_perc <-
      which(
        !is.na(
          suppressWarnings(
            as.numeric(
              unlist(
                strsplit(
                  ethf,
                  ""
                )
              )
            )
          )
        )
      )

    perc <-
      as.numeric(
        paste(
          maineth[ind_perc],
          collapse=
            ""
        )
      )

    char_ind <-
      which(
        is.na(
          suppressWarnings(
            as.numeric(
              unlist(
                strsplit(
                  ethf,
                  ""
                )
              )
            )
          )
        )
      )

    char <-
      paste(
        maineth[char_ind],
        collapse=
          ""
      )

    gnomAD_eth_indf <-
      which(
        grepl(
          char,
          colnames(
            dfMAF2f
          )
        )
      )

    out <-
      c(
        gnomAD_eth_indf,
        perc
      )

    return(out)

}

getWeightedeth <- function(gnomAD_eth_indf, perc, dfMAF2f=dfMAF2) {

    f1 <- function(x,y) {
        z <-
          dfMAF2f[,x]*(y/100)
        return(z)
    }

    dfMAF2f$AF_weighted <-
      apply(
        do.call(
          "cbind",
          mapply(
            f1,
            gnomAD_eth_indf,
            perc,
            SIMPLIFY=
              FALSE
            )
        ),
        1,
        sum
      )

    return(dfMAF2f)

}

getBinary <- function(df, thresholds, filter, int_MAF=internal_MAF_threshold) {

    if (filter == "MAF") {

        if (ncol(df) > 5) {


            if (length(which(apply(df,1,function(x) sum(is.na(x))==5)))>=1) {


                # Obtain index of variants that are not observed in gnomAD,
                # but present within your case dataset

                no_MAF <-
                  which(
                    apply(
                      df,
                      1,
                      function(x)
                      sum(
                        is.na(
                          x
                        )
                      )==5
                    )
                  )

                dft <-
                  as.data.frame(
                    t(
                      df
                    )
                  )

                dft[,no_MAF] <-
                  rep(
                    df[no_MAF,]$AF_int,
                    each=
                      ncol(
                        df
                      )
                  )

                df <-
                  as.data.frame(
                    t(
                      dft
                    )
                  )

                f1 <- function(x, y) {

                        indexed_rownames <-
                          rownames(df)

                        outa <-
                          transform(
                            df[-no_MAF,],
                            binary=
                              ifelse(
                                df[-no_MAF,][,x]<y & df[-no_MAF,][,ncol(df)]<int_MAF,
                                1,
                                0
                              )
                          ); rownames(outa) <- indexed_rownames[-no_MAF]

                        outb <-
                          transform(
                            df[no_MAF,],
                            binary=
                              ifelse(
                                apply(
                                  df[no_MAF,],
                                  1,
                                  function(x)
                                  all(
                                    x<int_MAF
                                  )
                                ),
                                1,
                                0
                              )
                          ); rownames(outb) <- no_MAF

                         out <-
                           rbind(
                             outa,
                             outb
                           )

                         suppressWarnings(out <-
                           out[order(as.numeric(rownames(out))),])

                    return(out)

                }

            } else {


                f1 <- function(x, y) {

                        out <-
                          transform(
                            df,
                            binary=
                              ifelse(
                                df[,x]<y & df[,ncol(df)]<int_MAF,
                                1,
                                0
                              )
                          )

                 return(out)

                 }

            }

            df <-
              df[,-6]

        } else if (ncol(df) == 5) {


            f1 <- function(x, y) {

                    out <-
                      transform(
                        df,
                        binary=
                          ifelse(
                            df[,x]<y,
                            1,
                            0
                          )
                      )

                return(out)

            }

        }

    suffix <-
      c(
        "nfe",
        "afr",
        "sas",
        "eas",
        "amr"
      )

    } else if (filter == "MCAP") {

        f1 <- function(x, y) {

            out <-
              transform(
                df,
                binary=
                  ifelse(
                    df[,x]>y,
                    1,
                    0
                  )
              )

            return(out)

        }

    suffix <-
      "all"

    } else {

        stop("Incorrent 'filter' input")

    }

    ind <-
      as.list(
        1:ncol(df)
      )

    thresholds_l <-
      as.list(
        rep(
          thresholds,
          1,
          each=
            ncol(
              df
            )
        )
      )

    df_binary <-
      mapply(
        f1,
        ind,
        thresholds_l,
        SIMPLIFY=
          FALSE
      )

    df_binary_comb <-
      do.call(
        "cbind",
        lapply(
          df_binary,
          function(x)
          x[,ncol(df)+1]
        )
      )

    rownames(df_binary_comb) <-
      ID$ID

    colnames(df_binary_comb) <-
      paste0(paste0(filter,"_",rep(thresholds,each=ncol(df)),"_",rep(suffix,length(thresholds))))

    return(df_binary_comb)

}

getExtract <- function(thresholds_a, thresholds_b, mat_a, mat_b) {

    if (!is.numeric(thresholds_a) | !is.numeric(thresholds_b)) {

        stop("Thresholds must be numeric vectors")

    }

    thresholds_a_l <-
      as.list(
        thresholds_a
      )

    thresholds_b_l <-
      as.list(
        thresholds_b
      )

    MAF_out <-
      lapply(
        thresholds_a_l,
        function(x)
        names(
          which(
            apply(
              mat_a[,which(grepl(x,colnames(mat_a)))],1,sum
            ) == 5
          )
        )
      )

    if (ncol(dfMCAP_binary)==1) {

        MCAP_out <-
        lapply(
          thresholds_b_l,
          function(x)
          names(
            mat_b[mat_b[,which(grepl(paste0("MCAP_",x,"_all"),colnames(mat_b)))]==1,]
          )
        )

    } else {

        MCAP_out <-
          lapply(
            thresholds_b_l,
            function(x)
            rownames(
              mat_b[mat_b[,which(grepl(paste0("MCAP_",x,"_all"),colnames(mat_b)))]==1,]
            )
          )

    }

    pairwise_int <-
    lapply(
      do.call(
        c,
        lapply(
          MAF_out,
          function(x)
          lapply(
            MCAP_out,
            function(y)
            intersect(x,y)
          )
        )
      ),
      function(x)
      as.character(
        x
      )
    )

    # Use pairwise_int to subset ID dataframe and generate gnomAD cumulative minor allele counts (CMAC)

    pairwise_int_a <-
      lapply(
        pairwise_int,
        function(x)
        ID[ID$ID %in% x,]
      )

    CMAC_f            <- function(n) {
        z <- sum(n)*2
        return(z)
    }

    CMAC <-
      lapply(
        pairwise_int,
        function(x)
        data.frame(
          Gene=
            names(
              tapply(
                ID[ID$ID %in% x,]$gnomAD_MAF,ID[ID$ID %in% x,]$Gene,CMAC_f
              )
            ),
          gnomAD_CMAC=
            as.numeric(
              tapply(
                ID[ID$ID %in% x,]$gnomAD_MAF,ID[ID$ID %in% x,]$Gene,CMAC_f
              )
            )
        )
     )

    # Merge on Gene

    f1_merge <- function (x,y,u) { z <- merge(x,y,by=u); return(z) }

    out <-
      mapply(
        f1_merge,
        pairwise_int_a,
        CMAC,
        "Gene",
        SIMPLIFY=
          FALSE
     )

    out_fin <-
      lapply(
        out,
        function(x)
        x[,-1]
      )

    # Save names as global list variable

    file_names <<-
      do.call(
        c,
        lapply(
          as.list(
            thresholds_a
          ),
          function(x)
          lapply(
            as.list(
              thresholds_b
            ),
            function(y)
            paste0(x,"_",y)
          )
        )
      )

    return(out_fin)

}

getWrite <- function(x,y) {

    out <-
      write.table(
        x,
        file=
          y,
        quote=
          F,
        row.names=
          F,
        col.names=
          F,
        sep=
          " ",
      )

    return(out)

}


# End of functions


# Load in testing/ranking and gnomAD rare variant burden matrices


args              <- commandArgs(TRUE)

outdir            <- args[1]
study             <- args[2]
input             <- args[3]

if (grepl(",",args[4])) {

    gnomAD_MAF_threshold <-
      as.numeric(
        unlist(
          strsplit(
            args[4],
            ","
          )
        )
      )

} else {

    gnomAD_MAF_threshold <-
      as.numeric(
        args[4]
      )

}

internal_MAF_threshold <- args[5]

if (grepl(",",args[6])) {

MCAP_threshold <-
  as.numeric(
    unlist(
      strsplit(
        args[6],
        ","
      )
    )
  )

} else {

MCAP_threshold   <- as.numeric(args[6])

}

if (grepl("_",args[7])) {

eth <-
  as.character(
    unlist(
      strsplit(
        args[7],
        "_"
      )
    )
  )

} else {

eth              <- args[7]

}

coverage         <- args[8]


dfAll <-
  read.table(
    file=
      input,
    header=
      TRUE,
    sep=
      "\t"
  )


# Index gnomAD MAF, internal MAF, and MCAP fields


MAF_fields <-
  which(
    names(
      dfAll
    ) %in%
    c(
      "AF_nfe",
      "AF_afr",
      "AF_sas",
      "AF_eas",
      "AF_amr",
      "AF_int"
    )
  )


MCAP_fields <-
  which(
    names(
      dfAll
    ) == "M.CAP_score"
  )

  dfMAF <-
    dfAll[,MAF_fields]

  dfMCAP <-
    dfAll[,MCAP_fields]


# N.B. Due to the potenial presence of "." in annoation fields, all numerical fields
# will automatically be read in as factor variables. So, must convert to numeric
# prior to proceeding to getBinary function


dfMAF <-
  suppressWarnings( # 'as.numeric' will intoduce 'NA' by coersion, this is normal
    mutate_all(     # behaviour and is the reason for the warning
      dfMAF,
      function(x)
      as.numeric(
        as.character(
          x
        )
      )
    )
  )


#


dfMAF2 <-
  suppressWarnings( # 'as.numeric' will intoduce 'NA' by coersion, this is normal
    mutate_all(     # behaviour and is the reason for the warning
      dfMAF,
      function(x)
      as.numeric(
        as.character(
          x
        )
      )
    )
  )


# Set 'NA' fields to 0 in dfMAF2 so that missing gnomAD allele frequencies
# will not be weighted if multiple 'eth' variables are provided
# dfMAF will be left as is, since missing gnomAD allele frequencies will be
# set to 'internal_MAF_threshold'


dfMAF2[is.na(dfMAF2)] <-
  0

dfMCAP <-
  suppressWarnings( # 'as.numeric' will introduce 'NA' by coersion, this is normal
    data.frame(     # this is normal behaviour and is the reason for the warning
      M.CAP_score=
        as.numeric(
          as.character(
            dfMCAP
          )
        )
    )
  )


# Set 'NA's in dfMCAP to 0 for now (these are variants that are not
# nonsynonymous SNVs or those that have a MAF >0.01)


dfMCAP[is.na(dfMCAP)] <-
  0


# Index gnomAD ancestry


if (length(eth) > 1) {

   eth_l <-
     as.list(
       eth
     )

    eth_ind_break_l <-
      as.list(
        unlist(
          lapply(
            eth_l,
            function(x)
            getethbreak(x)
          )
        )[seq(1,length(unlist(lapply(eth_l, function(x) getethbreak(x)))),2)]
      )
    perc_break_l <-
      as.list(
        unlist(
          lapply(
            eth_l,
            function(x)
            getethbreak(x)
          )
        )[seq(2,length(unlist(lapply(eth_l, function(x) getethbreak(x)))),2)]
      )

    dfMAF2 <-
      getWeightedeth(
        eth_ind_break_l,
        perc_break_l
      )

    gnomAD_eth_ind <-
      ncol(dfMAF2)

} else if (length(eth)==1) {

    gnomAD_eth_ind <-
      which(
        grepl(
          eth[[1]],
          names(
            dfMAF2
          )
        )
      )

}


# Establish variant field IDs


if (grepl("chr",dfAll$Chr)[1]) {

    suppressWarnings( # 'as.numeric' will be defulated to 'NA' by coersion if
      ID <-           # if gnomAD AF is missing; this is normal behaviour and
        data.frame(   # is the reason for the warning
          ID=
            paste0(
              gsub("chr","",dfAll$Chr),":",dfAll$Pos,":",dfAll$Ref,":",dfAll$Alt," ",dfAll$Alt2," ",1," ",dfAll$Gene.refGene),
          Gene=
            dfAll$Gene.refGene,
          gnomAD_MAF=
            as.numeric(
              as.character(
                dfMAF2[,gnomAD_eth_ind]
              )
            )
        )
    )

} else {

    suppressWarnings( # 'as.numeric' will be defulated to 'NA' by coersion if
      ID <-           # if gnomAD AF is missing; this is normal behaviour and
        data.frame(   # is the reason for the warning
          ID=
          paste0(
            dfAll$Chr,":",dfAll$Pos,":",dfAll$Ref,":",dfAll$Alt," ",dfAll$Alt2," ",1," ",dfAll$Gene.refGene),
          Gene=
            dfAll$Gene.refGene,
          gnomAD_MAF=
            as.numeric(
              as.character(
                dfMAF2[,gnomAD_eth_ind]
              )
            )
        )
    )

}


# Get binary matrices


dfMAF_binary <-
  getBinary(
    dfMAF,
    gnomAD_MAF_threshold,
    "MAF"
  )

dfMCAP_binary <-
  getBinary(
    dfMCAP,
    MCAP_threshold,
    "MCAP"
  )

# All variants that are 'disrputive' will automatically get a 1 for each
# MCAP category


disruptive <-
  c(
    "splicing",
    "stoploss",
    "stopgain",
    "frameshift_deletion",
    "frameshift_insertion"
  )

dfMCAP_binary[which(dfAll$Func.refGene %in% disruptive | dfAll$ExonicFunc.refGene %in% disruptive),] <-
  1


# Convert 'eth' back to being underscore delimited


if (length(eth) > 1) {

    eth <-
      paste(
        eth,
        collapse=
          "_"
      )

}


# Extract variants that intersect


setwd(outdir)

set <-
  mapply(
    getWrite,
    getExtract(
      gnomAD_MAF_threshold,
      MCAP_threshold,
      dfMAF_binary,
      dfMCAP_binary
    ),
    paste0(study,"_","variant_list_extracted_",file_names,"_",eth,"_",coverage,".txt"),
    SIMPLIFY=
      FALSE
 )
