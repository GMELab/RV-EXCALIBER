#!/bin/bash

#===============================================================================================================================
# Title: get_RVBurdenMatrix_internal_rvexcaliber.clean.sh
# Generate a rare variant burden matrix for your internal testing and/or ranking dataset using user defined gnomAD MAF, internal
# MAF, and MCAP thresholds
#===============================================================================================================================
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

# Copyright 2021 Ricky Lali, Michael Chong, Arghavan Omidi, Pedrum Mohammadi-Shemirani, Ann Le, Edward Cui, and Guillaume Pare
#===============================================================================================================================
# ESSENTIAL DEPENDENCIES:

# 1. ANNOVAR (version release: 2019-10-24
# reference: Wang K, Li M, Hakonarson H. ANNOVAR: functional annotation of
#            genetic variants from high-throughput sequencing data. Nucleic Acids Res.
#            2010;38(16):e164. doi:10.1093/nar/gkq603

# 2. plink version 1.9;
# reference: Chang CC, Chow CC, Tellier LC, Vattikuti S, Purcell SM, Lee JJ.
#            Second-generation PLINK: rising to the challenge of larger and
#            richer datasets. Gigascience. 2015;4(1):7.
#            doi:10.1186/s13742-015-0047-8

# 3. bedtools version 2.25.0;
# reference: Quinlan AR and Hall IM, 2010. BEDTools: a flexible suite of
#            utilities for comparing genomic features. Bioinformatics. 26, 6,
#            pp. 841–842)
#===============================================================================================================================




#-------------------------------------------------------------------------------------------------------------------------------
#	EXPLANATION OF INPUT ARGUMENTS:
#-------------------------------------------------------------------------------------------------------------------------------

# indir

# Description:     full path to the directory containing the input binary plink files (.bed, .fam, .bim)

# Expected input:  character input

#-------------------------------------------------------------------------------------------------------------------------------

# outdir

# Description:     full path to the directory to where all output files will be directed

# Expected input:  character input

#-------------------------------------------------------------------------------------------------------------------------------

# internal_dataset

# Description:     name of the binary plink files for your internal ** testing ** or ** ranking ** dataset (i.e. identical to
#                  what is specified with the plink --bfile or --out commands)

# Expected input:  character input

# Notes:           all output files will be named using the 'internal_dataset' variable

#-------------------------------------------------------------------------------------------------------------------------------

# gnomAD_MAF_threshold

# Description:     upper-limit minor allele frequency threshold from ** gnomAD ** used to filter variants

# Expected input:  comma separated numeric value for multiple thresholds (e.g. 0.001,0.01)

#                  OR

#                  single numeric value for a single threshold (e.g. 0.001)

#-------------------------------------------------------------------------------------------------------------------------------

# internal_MAF_threshold

# Description:     upper-limit minor allele frequency threshold from your testing or ranking dataset used to filter variants

# Expected input:  single numeric value for a single threshold (e.g. 0.001)

# Notes:           if your internal dataset has a small sample size (typically < 1000), the 'internal_MAF_threshold' should
#                  be higher than the 'gnomAD_MAF_threshold'

#-------------------------------------------------------------------------------------------------------------------------------

# MCAP_threshold

# Description:     lower-limit threshold MCAP score used to filter nonsynonymous SNVs

# Expected input:  comma separated numeric value for multiple thresholds (e.g. 0.009,0.025)

#                  OR

#                  single numeric value for a single threshold (e.g. 0.009)

#-------------------------------------------------------------------------------------------------------------------------------

# eth

# Description:     ethnic group of the comparator gnomAD frequency you selected in the
#                  get_RVBurdenMatrix_internal_rvexcaliber.clean.sh script

# Expected input:  character input, can be one of < "nfe", "afr", "sas", "eas", "amr" >
#                  "nfe" = non-Finnish Euroepean, "afr" = African, "sas" = South Asian, "eas" = East Asian; "amr" = Latino

#                  OR

#                  If you selected multiple 'eth' variables to generate a weighted comparator gnomAD allele frequency,
#                  please separate the variables and their percentage weights with with an underscore (e.g. "nfe90_afr10")

# Notes:           gene-based gnomAD allele counts will be calculated based on the selected comparator ethnic group in gnomAD

#-------------------------------------------------------------------------------------------------------------------------------

# coverage

# Description:     conduct site-based intersection of input internal_dataset with high coverage coding (hcc) regions in gnomAD

# Expected input:  can be one of < "hcc", "rcc" >
#                  "hcc" = high-coverage coding; "rcc" = regular-coverage coding

# Notes:           high coverage coding regions are defined as sites in gnomAD that have >20X coverage in >90% of individuals

#-------------------------------------------------------------------------------------------------------------------------------

# job_interval

# Description:     number of genes that are applied in parallel to generate final rare variant burden matrix

# Expected input:  whole-number, numeric integer (no decimal) input

# Notes:           number of cpu cores being used will match the job interval

#-------------------------------------------------------------------------------------------------------------------------------




#-------------------------------------------------------------------------------------------------------------------------------
#	PREPARING THE SCRIPT: 4 entries must be entered by the user and are explained below
#-------------------------------------------------------------------------------------------------------------------------------


# ** < Please enter the full path to the cloned '/rvexcaliber' directory (example below): > **


path_to_rvexcaliber="/torrent_vol_2018/lalir_dbGaP_repositories/temp_out/github_RV-EXCALIBER_ADSP/RV-EXCALIBER"


# Example: "/genetics/rvexcaliber"


# ** < Please enter the full path the directory containing the ANNOVAR perl scripts (example below): > **


annovar="/genetics_work/lalir/programs/annovar"


# Example: "/genetics/tools/annovar"


# The following scripts were used to obtain the ** necessary ** input files from ANNOVAR used by rvexcaliber

# 1. The refGene database:
# perl ${annovar}/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene ${annovar}/humandb/

# 2. The gnomAD (version 2.1.1) allele frequency exome annotations
# perl ${annovar}/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad211_exome  ${annovar}/humandb/

# 3. The in silico protein-scoring annotations from dbnsfp
# perl ${annovar}/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp35c  ${annovar}/humandb/


# ** < Please enter the full path the plink (version 1.9) command-line executable (example below): > **


plink="/genetics/PAREG/common/plink_1.9/plink"


# Example: "/genetics/tools/plink_1.9/plink"


# ** < Please enter the full path to the 'intersectBed' command-line executable from bedtools (example below) > **
# Note: this path can be left blank if the 'coverage' is set to "rcc"


bedtools="/genetics_work/lalir/programs/bedtools2/bin/intersectBed"


# Example: "/genetics/tools/bedtools2/bin/intersectBed"


#-------------------------------------------------------------------------------------------------------------------------------




#-------------------------------------------------------------------------------------------------------------------------------
# Check script preparation
#-------------------------------------------------------------------------------------------------------------------------------


if [[ ! -d ${path_to_rvexcaliber} ]]; then

    echo -e "\n"
    echo "Error: The full path to the /RV-EXCALIBER directory that was downloaded locally does not exist"
    echo "Please ensure that this path is correctly defined above, under the "PREPARING THE SCRIPT" section"
    exit 1

elif [[ ! -d ${annovar} ]]; then

    echo -e "\n"
    echo "Error: The full path to the /annovar directory does not exist"
    echo "Please ensure that this path is correctly defined above, under the "PREPARING THE SCRIPT" section"
    exit 1

elif [[ ! -f ${plink} ]]; then

    echo -e "\n"
    echo "Error: The full path to the plink v1.9 command-line executable does not exist"
    echo "Please ensure that this path is correctly defined above, under the "PREPARING THE SCRIPT" section"
    exit 1

elif [[ ! -f ${bedtools} ]]; then

    echo -e "\n"
    echo "Error: The full path to the full path to the bedtools' 'intersectBed' command-line executable does not exist"
    echo "Please ensure that this path is correctly defined above, under the "PREPARING THE SCRIPT" section"
    exit 1

fi


#-------------------------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------------------
# Defaulted directories are provided below (do not change):
#-------------------------------------------------------------------------------------------------------------------------------


scripts="${path_to_rvexcaliber}/scripts"
gnomAD="${path_to_rvexcaliber}/gnomAD_211_exomes_hg19"
gnomAD_filter="${path_to_rvexcaliber}/gnomAD_211_exomes_hg19_filter"
gnomAD_coverage="${path_to_rvexcaliber}/gnomAD_211_exomes_hg19_coverage"


#-------------------------------------------------------------------------------------------------------------------------------




# Directory and naming variables


indir=$1
outdir=$2
internal_dataset=$3


# Filter variables


gnomAD_MAF_threshold=$4
internal_MAF_threshold=$5
MCAP_threshold=$6
eth=$7
coverage=$8


# Miscellaneous variables


job_interval=$9


# Default changes and variables


indir=$(echo ${indir} | sed 's/\/$//g')

outdir=$(echo ${outdir} | sed 's/\/$//g')

reference="hg19"

function stop_if_error {

    exit_status=$?

    if [ ${exit_status} -ne 0 ]; then

        echo -e "\n"

        echo "Please review error"

        exit "${exit_status}"

    fi

}


# Script start


# Check 'user entered' input variables


if [[ $# = 9 ]]; then

    if [[ -f ${indir}/${internal_dataset}.bed ]] && [[ -f ${indir}/${internal_dataset}.bim ]] && [[ -f ${indir}/${internal_dataset}.fam ]]; then

        dig='^[0-9]+(\.[0-9]+)?$'
        digw='^[0-9]+$'

        # Check variables for gnomAD MAF thresholds

        gnomAD_MAF_threshold_count=$(echo ${gnomAD_MAF_threshold} | awk -F"," '{for(i=1;i<=NF;i++) print $i}' | awk 'END { print NR }')

        gnomAD_MAF_threshold_dec=$(echo ${gnomAD_MAF_threshold} | awk -F"," '{for(i=1;i<=NF;i++) print $i}' | grep -wc ".")

        gnomAD_MAF_threshold_check=$(echo ${gnomAD_MAF_threshold} | sed -e 's/,//g' -e 's/\.//g')

        gnomAD_MAF_ge_1=$(echo ${gnomAD_MAF_threshold} | awk -F"," '{for(i=1;i<=NF;i++) print $i}' | awk '$1>=1' | awk 'END { print NR}')

       # Check variables for internal MAF thresholds

        internal_MAF_threshold_dec=$(echo ${internal_MAF_threshold} | grep -wc ".")

        internal_MAF_threshold_check=$(echo ${internal_MAF_threshold} | sed -e 's/,//g' -e 's/\.//g')

        internal_MAF_ge_1=$(echo ${internal_MAF_threshold} | awk -F"," '{for(i=1;i<=NF;i++) print $i}' | awk '$1>=1' | awk 'END { print NR}')

        # Check variables for MCAP thresholds

        MCAP_threshold_count=$(echo ${MCAP_threshold} | awk -F"," '{for(i=1;i<=NF;i++) print $i}' | awk 'END { print NR }')

        MCAP_threshold_dec=$(echo ${MCAP_threshold} | awk -F"," '{for(i=1;i<=NF;i++) print $i}' | grep -wc ".")

        MCAP_threshold_check=$(echo ${MCAP_threshold} | sed -e 's/,//g' -e 's/\.//g')

        MCAP_ge_1=$(echo ${MCAP_threshold} | awk -F"," '{for(i=1;i<=NF;i++) print $i}' | awk '$1>=1' | awk 'END { print NR}')

        # Check variables for gnomAD ethnicity

        eth_count=$(echo ${eth} | awk -F"_" '{for(i=1;i<=NF;i++) print $i}' | awk 'END { print NR }')

        eth_check=$(echo ${eth} | awk -F"_" '{for(i=1;i<=NF;i++) print $i}' | sed 's/[^a-z]*//g' | grep -Ewc "nfe|afr|sas|eas|amr")

        eth_num_check=$(echo ${eth} | awk -F"_" '{for(i=1;i<=NF;i++) print $i}' | sed 's/[^0-9.]*//g' | awk 'NF>0' | awk 'END { print NR}')

        eth_num_sum=$(echo ${eth} | awk -F"_" '{for(i=1;i<=NF;i++) print $i}' | sed 's/[^0-9.]*//g' | awk '{ sum += $1} END { print sum}')

        if [[ ! ${gnomAD_MAF_threshold_check} =~ ${dig} ]] || [[ ${gnomAD_MAF_threshold_count} -ne ${gnomAD_MAF_threshold_dec} ]]; then

            echo -e "\n"
            echo "Error: Ensure that the gnomAD MAF thresholds (i.e. 'gnomAD_MAF_threshold') are numeric and *comma separated* if using multiple MAF thresholds"
            echo "See 'Explanation of input arguments' for further information"
            exit 1

        elif [[ ${gnomAD_MAF_ge_1} -ge 1 ]]; then

            echo -e "\n"
            echo "Error: Ensure all gnomAD MAF thresholds (i.e. 'gnomAD_MAF_threshold') are less than 1"
            echo "See 'Explanation of input arguments' for further information"
            exit 1

        elif [[ ! ${internal_MAF_threshold_check} =~ ${dig} ]] || [[ ${internal_MAF_threshold_dec} -ne 1 ]]; then

            echo -e "\n"
            echo "Error: Ensure that the internal MAF threshold (i.e. 'internal_MAF_threshold) is a single numeric value"
            echo "See 'Explanation of input arguments' for further information"
            exit 1

        elif [[ ${internal_MAF_ge_1} -ge 1 ]]; then

            echo -e "\n"
            echo "Error: Ensure all internal MAF threshold (i.e. 'internal_MAF_threshold') is less than 1"
            echo "See 'Explanation of input arguments' for further information"
            exit 1

        elif [[ ! ${MCAP_threshold_check} =~ ${dig} ]] || [[ ${MCAP_threshold_count} -ne ${MCAP_threshold_dec} ]]; then

            echo -e "\n"
            echo "Error: Ensure that the MCAP thresholds (i.e. 'MCAP_threshold') are numeric and *comma separated* if using multiple MCAP thresholds"
            echo "See 'Explanation of input arguments' for further information"
            exit 1

        elif [[ ${MCAP_ge_1} -ge 1 ]]; then

            echo -e "\n"
            echo "Error: Ensure all MCAP thresholds (i.e. 'MCAP_threshold') are less than 1"
            echo "See 'Explanation of input arguments' for further information"
            exit 1

        elif [[ ${eth_count} -ne ${eth_check} ]]; then

            echo -e "\n"
            echo "Error: Ensure that the comparator gnomAD ethnicity variable (i.e. 'eth') is one of < "nfe", "afr", "sas", "eas", "amr" >"
            echo "or"
            echo "If you wish to generate a weighted comparator gnomAD allele frequency based on 2 or more individual 'eth' variables, please ensure they are separated by an underscore (_)"
            echo "See 'Explanation of input arguments' for further information"
            exit 1

        elif [[ ${eth_count} -gt 1 ]] && [[ ${eth_check} -ne ${eth_num_check} ]]; then

            echo -e "\n"
            echo "If you wish to generate a weighted comparator gnomAD allele frequency based on 2 or more individual comparator gnomAD ethnicity variables (i.e. 'eth'), please ensure there is a percentage weight is assigned to each 'eth' variable"
            echo "See 'Explanation of input arguments' for further information"
            exit 1

        elif [[ ${eth_count} -gt 1 ]] && [[ ${eth_num_sum} -ne 100 ]]; then

            echo -e "\n"
            echo "Error: If you wish to generate a weighted comparator gnomAD allele frequency based on 2 or more individual comparator gnomAD ethnicity variables (i.e. 'eth'), please ensure that their corresponding percentage weights add to 100"
            echo "See 'Explanation of input arguments' for further information"
            exit 1

        elif [[ ! ${coverage} = "hcc" ]] && [[ ! ${coverage} = "rcc" ]]; then

            echo -e "\n"
            echo "Error: Ensure that the 'high coverage coding  intersection variable' (i.e. 'coverage') is one of < "hcc" or "rcc" >"
            echo "See 'Explanation of input arguments' for further information"
            exit 1

        elif [[ ! ${job_interval} =~ ${digw} ]]; then

            echo -e "\n"
            echo "Error: Ensure that the number of parallel genes to be added to the rare variant burdent matrix (i.e. 'job_interval') is a single numeric integer"
            echo "See 'Explanation of input arguments' for further information"
            exit 1

        elif [[ ${job_interval} -lt 1 ]]; then

            echo -e "\n"
            echo "Error: Ensure that the number of parallel genes to be added to the rare variant burdent matrix (i.e. 'job_interval') is at least 1"
            echo "See 'Explanation of input arguments' for further information"
            exit 1

        elif [[ ! -d ${outdir} ]]; then

            echo -e "\n"
            echo "Error: output directory (i.e. 'outdir' does not exist)"
            echo "See 'Explanation of input arguments' for further information"
            exit 1

        else

            echo -e "\n"
            echo "Input variables successfully entered. Proceeding..."
            sleep 2s

        fi


        # Post check variable assignments


        function reassign {

            MAF_MCAP=$(for MAF in $(echo ${gnomAD_MAF_threshold} |
            awk -F"," '{ for(i=1;i<=NF;i++) print $i }'); do
            for MCAP in $(echo ${MCAP_threshold} |
            awk -F"," '{ for(i=1;i<=NF;i++) print $i }'); do
            echo $MAF"_"$MCAP;
            done;
            done)

       }
       reassign

        if [[ ${eth_count} -eq 1 ]]; then

            eth=$(echo ${eth} | sed 's/[^a-z]*//g')

        fi


        # Start script


        # Re-assign ID column in .bim file for the internal testing or ranking
        # dataset


        function convert_bim {(

          set -e

          awk 'BEGIN { FS=OFS="\t" } {

              $2=$1":"$4":"$6":"$5

          }1' ${indir}/${internal_dataset}.bim > ${outdir}/${internal_dataset}_re.bim

          cp -r ${indir}/${internal_dataset}.bed ${outdir}/${internal_dataset}_re.bed

          cp -r ${indir}/${internal_dataset}.fam ${outdir}/${internal_dataset}_re.fam

        )}
        convert_bim
        stop_if_error


        # For all variants that are mutually observed in the internal
        # testing/ranking dataset and gnomAD, remove those that do not have a
        # FILTER status of PASS in gnomAD


        function gnomAD_PASS_crossref {(

            set -e

            echo -e "\n"
            echo "Cross-referencing variants with gnomAD PASS sites"

            rm -f ${outdir}/${internal_dataset}_preprocessed.bed

            rm -f ${outdir}/${internal_dataset}_preprocessed.bim

            rm -f ${outdir}/${internal_dataset}_preprocessed.fam

            # Print ID column and remove asterisk alleles (if any)

            awk '{

                if ($5!="*")

                    print $2

            }' ${outdir}/${internal_dataset}_re.bim > ${outdir}/${internal_dataset}_sites.txt

            awk 'NR==FNR {

                 FILTER[$1]=$2; next

            }

            {

                print $0, ($1 in FILTER ? FILTER[$1] : "Not_obs")

            }' ${gnomAD_filter}/ALL_CHROM_gnomAD_filter.txt ${outdir}/${internal_dataset}_sites.txt |

            grep -Ew "PASS|Not_obs" |

            awk '{

                print $1

            }' > ${outdir}/${internal_dataset}_gnomAD_PASS_extract

            ${plink} --noweb \
            \
            --silent \
            \
            --bfile ${outdir}/${internal_dataset}_re \
            \
            --extract ${outdir}/${internal_dataset}_gnomAD_PASS_extract \
            \
            --keep-allele-order \
            \
            --make-bed \
            \
            --out ${outdir}/${internal_dataset}_preprocessed

            rm -f ${outdir}/${internal_dataset}_sites.txt

            rm -f ${outdir}/${internal_dataset}_re.bed

            rm -f ${outdir}/${internal_dataset}_re.bim

            rm -f ${outdir}/${internal_dataset}_re.fam

            rm -f ${outdir}/${internal_dataset}_gnomAD_PASS_extract

        )}
        gnomAD_PASS_crossref
        stop_if_error


        # Intersect regions in the internal testing or ranking dataset with the
        # high coverage coding regions in gnomAD


        function coverage {(

            set -e

            if [[ ! ${coverage} = "hcc" ]]; then

                echo -e "\n"
                echo "No coverage intersection selected. Will not intersect ${internal_dataset} with gnomAD high coverage coding regions"

                sleep 2s

            else

                echo -e "\n"
                echo "Intersecting ${internal_dataset} with gnomAD high coverage coding regions"

                sleep 2s

                awk '{

                    print $2

                }' ${outdir}/${internal_dataset}_preprocessed.bim |

                awk -F":" '{

                    print $1"\t"$2"\t"$2"\t"$0

                }' > ${outdir}/${internal_dataset}_preprocessed_hold.bed

                ${bedtools} \
                \
                -a ${gnomAD_coverage}/ALL_CHROM_gnomad_refGene_exons_highcov20X.mod.bed \
                \
                -b ${outdir}/${internal_dataset}_preprocessed_hold.bed \
                \
                -wb |

                cut -f 10 - |

                awk -F":" '!seen[$1,$2,$3,$4]++' > ${outdir}/${internal_dataset}_gnomAD_coverage_extract

                ${plink} --noweb \
                \
                --silent \
                \
                --bfile ${outdir}/${internal_dataset}_preprocessed \
                \
                --extract ${outdir}/${internal_dataset}_gnomAD_coverage_extract \
                \
                --keep-allele-order \
                \
                --make-bed \
                \
                --out ${outdir}/${internal_dataset}_preprocessed_coverage_complete

                rm -f ${outdir}/${internal_dataset}_preprocessed_hold.bed

                rm -f ${outdir}/${internal_dataset}_gnomAD_coverage_extract

                # *_preprocessed_coverage_complete.bed, _preprocessed_coverage_complete.bim,
                # and _preprocessed_coverage_complete.fam plink files will be renamed to
                # _preprocessed.bed _preprocessed.bim _preprocessed.fam, respectively

                mv ${outdir}/${internal_dataset}_preprocessed_coverage_complete.bed ${outdir}/${internal_dataset}_preprocessed.bed

                mv ${outdir}/${internal_dataset}_preprocessed_coverage_complete.bim ${outdir}/${internal_dataset}_preprocessed.bim

                mv ${outdir}/${internal_dataset}_preprocessed_coverage_complete.fam ${outdir}/${internal_dataset}_preprocessed.fam

                rm -f ${outdir}/${internal_dataset}_*preprocessed_coverage_complete*

            fi

        )}
        coverage
        stop_if_error


        # Generate ANNOVAR-ready file


        function prep_for_annovar {(

            set -e

            rm -f ${outdir}/${internal_dataset}.annovarInput*

            echo -e "\n"
            echo "Initiating refGene, gnomAD_211, and dbNSFP M-CAP annotations for ${internal_dataset}"

            rm -f ${outdir}/${internal_dataset}.annovarInput

            awk 'BEGIN {FS=OFS="\t"} {

                print $1, $4, $4, $6, $5

            }' ${outdir}/${internal_dataset}_preprocessed.bim |

            awk 'BEGIN {FS=OFS="\t"} {

              if  (length($4) > length($5)) $3=$3+(length($4)-length($5))

            }1' > ${outdir}/${internal_dataset}.annovarInput

        )}
        prep_for_annovar
        stop_if_error


        # Annotate variants in the internal testing or ranking set with refGene,
        # gnomAD 211, and dbNSFP


        {

        function annovar_annotation {(

            set -e

            rm -f ${outdir}/${internal_dataset}.${reference}_multianno*

            perl ${annovar}/table_annovar.pl \
            \
            -buildver ${reference} \
            \
            ${outdir}/${internal_dataset}.annovarInput \
            \
            -remove \
            \
            -protocol refGene,gnomad211_exome,dbnsfp35c \
            \
            -operation g,f,f \
            \
            -nastring . \
            \
            -polish \
            \
            ${annovar}/humandb/ \
            \
            -out ${outdir}/${internal_dataset}

            gzip -f ${outdir}/${internal_dataset}.annovarInput

        )}

        annovar_annotation
        stop_if_error

        } &> ${outdir}/${internal_dataset}.annovar.output


        # Extract allele frequencies from the inernal testing or ranking
        # dataset


        function get_internal_AF {(

            set -e

            ${plink} --noweb \
            \
            --silent \
            \
            --bfile ${outdir}/${internal_dataset}_preprocessed \
            \
            --freq \
            \
            --keep-allele-order \
            \
            --out ${outdir}/${internal_dataset}_preprocessed

            awk 'OFS="\t" {$1=$1}1' ${outdir}/${internal_dataset}_preprocessed.frq |

            awk -F"[\t:]" 'BEGIN {

                print "Chr", "Pos", "Ref", "Alt", "AF_int"

            } NR>1 {

                print $2, $3, $4, $5, $8

            }' |

            awk 'OFS="\t" {$1=$1}1 ' > ${outdir}/${internal_dataset}_preprocessed.mod.frq

            rm -f ${outdir}/${internal_dataset}_preprocessed.frq

        )}
        get_internal_AF
        stop_if_error


        # Generate a re-formatted per-variant annotation file


        function prep_for_R_input {(

            set -e

            rm -f ${outdir}/${internal_dataset}_pruned_annotation_for_R_input_${coverage}*

            echo -e "\n"
            echo "Re-arranging annovar annotation file to a format used for R-script filtering"

            egrep -w 'refGene|exonic|splicing' ${outdir}/${internal_dataset}.${reference}_multianno.txt |

            awk 'BEGIN {

                FS=OFS="\t"

            }

            NR==1 {

                for (i=1;i<=NF;i++) id[$i]=i;

                    next
            }

            {

                print $(id["Chr"]), $(id["Start"]), $(id["Ref"]), $(id["Alt"]), $(id["Alt"]),

                $(id["Func.refGene"]), $(id["Gene.refGene"]), $(id["ExonicFunc.refGene"]),

                $(id["AF_nfe"]), $(id["AF_afr"]), $(id["AF_sas"]),

                $(id["AF_eas"]), $(id["AF_amr"]), $(id["M-CAP_score"])

            }' |

            awk 'BEGIN {

                OFS="\t"

            } {$1=$1}1' |

            sed -e 's/nonsynonymous\tSNV/nonsynonymous_SNV/g' \
                \
                -e 's/synonymous\tSNV/synonymous_SNV/g' \
                \
                -e 's/frameshift\tsubstitution/frameshift_substitution/g' \
                \
                -e 's/nonframeshift\tsubstitution/nonframeshift_substitution/g' |

            awk 'OFS="\t" {

                if ((length($3) > length($4)) && ($8=="frameshift_substitution"))

                    $8="frameshift_deletion";

                else if ((length($3) < length($4)) && ($8=="frameshift_substitution"))

                    $8="frameshift_insertion";

                else if ((length($3) > length($4)) && ($8=="nonframeshift_substitution"))

                    $8="nonframeshift_deletion";

                else if ((length($3) < length($4)) && ($8=="nonframeshift_substitution"))

                    $8="nonframeshift_insertion"

            }1' |

            # Remove benign variants from file

            egrep -vw 'nonframeshift_insertion|nonframeshift_deletion|synonymous_SNV|unknown' |

            # Re-header

            awk 'BEGIN {

                FS=OFS="\t";

                print "Chr", "Pos", "Ref", "Alt", "Alt2",

                "Func.refGene", "Gene.refGene", "ExonicFunc.refGene",

                "AF_nfe", "AF_afr", "AF_sas",

                "AF_eas", "AF_amr", "M-CAP_score"

            }1' |

            # Add in internal allele frequencies

            awk 'BEGIN {

                FS=OFS="\t"

            }

            NR==FNR {

                AF[$1,$2,$3,$4]=$5; next

            }

            {

                print $0, (($1,$2,$3,$4) in AF ? AF[$1,$2,$3,$4] : "NA")

            }' ${outdir}/${internal_dataset}_preprocessed.mod.frq - |

            # Re-arrange fields

            awk 'BEGIN {

                FS=OFS="\t"

            }

            {

                print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $15, $14

            }' > ${outdir}/${internal_dataset}_pruned_annotation_for_R_input_${coverage}.txt

            gzip -f ${outdir}/${internal_dataset}.${reference}_multianno.txt

            rm -f ${outdir}/${internal_dataset}_preprocessed.frq

            rm -f ${outdir}/${internal_dataset}_preprocessed.mod.frq

        )}
        prep_for_R_input
        stop_if_error


        # Rename variant IDs for bim and R_input files


        function rename_ID {

            item=$1
            file=$2

            if [[ ${item} = "bim" ]]; then

                outdir_func=${outdir}

                awk 'BEGIN {

                    FS=OFS="\t"

                }

                {

                    print $2

                }' ${outdir_func}/${file} > ${outdir_func}/${file}.tmp

                file2=${file}.tmp

            elif [[ ${item} = "R_input" ]]; then

                outdir_func=${outdir}

                awk 'NR > 1 {

                    print $1":"$2":"$3":"$4

                }' ${outdir_func}/${file} > ${outdir_func}/${file}.tmp

                file2=${file}.tmp

                awk 'NR==1' ${outdir_func}/${file} > ${outdir_func}/${file2}_1

                awk 'NR>1' ${outdir_func}/${file} > ${outdir_func}/${file}_0

                mv ${outdir_func}/${file}_0 ${outdir_func}/${file}

            fi

            sed 's/:/\t/g' ${outdir_func}/${file2} |

            awk 'BEGIN {

                FS=OFS="\t"

            }

            {

                if (length($4)>length($3)) $5="I";

                else if (length($4)<length($3)) $5="D";

                else $5=$3

            }1' |

            awk 'BEGIN { FS=OFS="\t" } {

                print $0, length($4)-length($3)

            }' |

            awk 'BEGIN {

                FS=OFS="\t"

            }

            {

                if ($6==0) $6=$4

            }1' |

            awk 'BEGIN {

                FS=OFS="\t"

            }

            {

                if ($6<0) $6=($6*-1)

            }1' |

            awk '{

                print $1":"$2":"$5":"$6

            }' |

            awk 'dup_id[$1]++ {

                $1=$1"_var_"dup_id[$1]-1

            }1' |

            paste -d"\t" - ${outdir_func}/${file} > ${outdir_func}/${file2}_0

            if [[ ${item} = "bim" ]]; then

                awk 'BEGIN {

                    FS=OFS="\t"

                }

                {

                    print $2, $1, $4, $5, $6, $7

                }' ${outdir_func}/${file2}_0 > ${outdir_func}/${file2}_1

                mv ${outdir_func}/${file2}_1 ${outdir_func}/${file}

                rm -f ${outdir_func}/${file2}*

            elif [[ ${item} = "R_input" ]]; then

                sed 's/:/\t/g' ${outdir_func}/${file2}_0 |

                cut -f 1-4,9- >> ${outdir_func}/${file2}_1

                mv ${outdir_func}/${file2}_1 ${outdir_func}/${file}

                gzip -f ${outdir_func}/${file}

                rm -f ${outdir_func}/${file2}*

            fi

        }


        # Execute the rename_ID function, which re-names the indels in the
        # ID column of the .bim files to abbreviated versions


        function execute_rename_ID_bim {(

            set -e

            rename_ID bim ${internal_dataset}_preprocessed.bim

        )}
        execute_rename_ID_bim
        stop_if_error


        # Execute the rename_ID function, which re-names the indels in the
        # ID column of the R_input files to abbreviated versions


        function execute_rename_ID_R_input {(

            set -e

            rename_ID R_input ${internal_dataset}_pruned_annotation_for_R_input_${coverage}.txt

        )}
        execute_rename_ID_R_input
        stop_if_error


        # Generate a list of variants that meet the user-defined gnomAD MAF
        # threshold (i.e. 'gnomAD_MAF_theshold'), MCAP threshold
        # (i.e. 'MCAP_threshold'); all variants will also meet the user-defined
        # internal MAF threshold (i.e. 'internal_MAF_threshold')


        function get_Varlist {(

            set -e

            echo -e "\n"
            echo -e "Filtering:\ngnomAD MAF:${gnomAD_MAF_threshold}\ninternal MAF:${internal_MAF_threshold}\nMCAP:${MCAP_threshold}"

            sleep 2s

            echo -e "\n"
            echo "Filtering will be conducted on ${internal_dataset}"

            sleep 2s

            rm -f ${outdir}/${internal_dataset}_variant_list_extracted_${MAF_MCAP}_${coverage}*

            Rscript ${scripts}/Rscripts/get_Varlist_rvexcaliber.R ${outdir} ${internal_dataset} ${outdir}/${internal_dataset}_pruned_annotation_for_R_input_${coverage}.txt.gz ${gnomAD_MAF_threshold} ${internal_MAF_threshold} ${MCAP_threshold} ${eth} ${coverage} "N"

        )}
        get_Varlist
        stop_if_error


        # Leverage the plink v1.9 --score function to generate per-individual
        # observed allele counts in the internal testing or ranking dataset


        # Calculates per-gene burden score based on the observe allele counts in
        # the internal testing or ranking dataset and generates rare variant
        # burden matrix


        function sub_generate_RVBurdenMatrix_internal {(

            set -e

            echo -e "\n"
            echo "Generating rare-variant gene matrix for ${internal_dataset} and initializing rare variant burden scores per-gene"

            sleep 2s

            for MAF_MCAP in ${MAF_MCAP}; do

                rm -f ${outdir}/${internal_dataset}_RVBurdenMatrix_${MAF_MCAP}_${eth}_${coverage}*

                bash ${scripts}/shell/sub_generate_RVBurdenMatrix_internal_rvexcaliber.clean.sh ${indir} ${outdir} ${internal_dataset} ${MAF_MCAP} ${eth} ${coverage} ${job_interval} ${plink}

                gzip -f ${outdir}/${internal_dataset}_variant_list_extracted_${MAF_MCAP}_${eth}_${coverage}.txt

                gzip -f ${outdir}/${internal_dataset}_RVBurdenMatrix_${MAF_MCAP}_${eth}_${coverage}.txt

            done
            reassign

        )}
        sub_generate_RVBurdenMatrix_internal
        stop_if_error

        for MAF_MCAP in ${MAF_MCAP}; do

            if  [[ -f ${outdir}/${internal_dataset}_RVBurdenMatrix_${MAF_MCAP}_${eth}_${coverage}.txt.gz ]]; then

                echo -e "\n"
                echo "A rare variant burden matrix has been generated for ${internal_dataset} at a MAF_MCAP: ${MAF_MCAP}, eth: ${eth}, and coverage: ${coverage}"
                echo -e "\n"
                echo "** Successful completion **"

            else

                echo -e "\n"
                echo "Error: A rare variant burden matrix has not been generated for ${internal_dataset} at a MAF_MCAP: ${MAF_MCAP}, eth: ${eth}, and coverage: ${coverage}"

            fi

        done
        reassign

        rm -f ${outdir}/${internal_dataset}_preprocessed.bed

        rm -f ${outdir}/${internal_dataset}_preprocessed.bim

        rm -f ${outdir}/${internal_dataset}_preprocessed.fam

        rm -f ${outdir}/${internal_dataset}_preprocessed.log

        rm -f ${outdir}/${internal_dataset}_preprocessed.nosex

    else

        echo -e "\n"
        echo "Error: Ensure all binary plink files are present in the input directory and that the input directory is correct."

    fi

else

    echo "Error: Incorrect number of inputs"

fi
