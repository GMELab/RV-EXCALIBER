#!/bin/bash

#================================================================================================================================
# Title: get_SummaryAssociations_iCF_gCF_rvexcaliber.clean.sh
# Conducts a gene-based rare variant association test using RV-EXCALIBER
#================================================================================================================================
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

# The R Project for Statistical Computing
# libraries within R that are reqruired
#    ggplot2
#    dplyr
#    data.table
#===============================================================================================================================




#-------------------------------------------------------------------------------------------------------------------------------
# EXPLANATION OF INPUT ARGUMENTS
#-------------------------------------------------------------------------------------------------------------------------------

# outdir

# Description:     full path to the directory to where all output files will be directed

# Expected input:  character input

#-------------------------------------------------------------------------------------------------------------------------------

# internal_testing_dataset

# Description:     name of your  internal ** testing ** dataset (i.e. identical to what was specified for the 'internal_study'
#                  variable in the 'get_RVBurdenMatrix_internal_rvexcaliber.clean.sh' script) or the name of the
#                  'combined_internal_study' variable in the 'get_combined_RVBurdenMatrix_internal_rvexcaliber.clean.sh' script
#                  if you had to combine your internal testing dataset across chromosomes

# Expected input:  character input

# Notes:           all output files will be named using the 'study' variable

#-------------------------------------------------------------------------------------------------------------------------------

# internal_ranking_dataset

# Description:     name of your  internal ** ranking ** dataset (i.e. identical to what was specified for the 'internal_study'
#                  variable in the 'get_RVBurdenMatrix_internal_rvexcaliber.clean.sh' script) or the name of the
#                  'combined_internal_study' variable in the 'get_combined_RVBurdenMatrix_internal_rvexcaliber.clean.sh' script
#                  if you had to combine your internal ranking dataset across chromosomes

# Expected input:  character input, can be the full path to ranking rare variant burden matrix or NA if you do not have a
#                  dedicated internal ranking dataset

# Notes:           all output files will be named using the 'study' variable

#-------------------------------------------------------------------------------------------------------------------------------

# gnomAD_MAF_threshold

# Description:     upper-limit minor allele frequency threshold from ** gnomAD ** used to filter variants

# Expected input:  comma separated numeric value for multiple thresholds (e.g. 0.001,0.01)

#                  OR

#                  single numeric value for a single threshold (e.g. 0.001)

# Notes:           The gnomAD minor allele frequency thresholds must be greater than 1e-05 and less than 1

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

# Description:     conduct site-based intersection of input internal_study with high coverage coding (hcc) regions in gnomAD

# Expected input:  can be one of < "hrr", "rcc" >
#                  "hcc" = high-coverage coding; "rcc" = regular-coverage coding

# Notes:           high coverage coding regions are defined as sites in gnomAD that have >20X coverage in >90% of individuals

#-------------------------------------------------------------------------------------------------------------------------------

# adjust_type:

# Description:     type of adjustment to conduct on gnomAD cumulative minor allele counts (i.e. the expected counts)

# Expected input:  character input, can be on of < "noadjust", "fulladjust" >
#                  "noadjust"  = no adjustments will be performed on the gnomAD allele counts,
#                  "fulladjust" = gnomAD counts will be adjusted for according to the individual and gene correction factors
#                  (i.e. iCF and gCF, respectively)

#-------------------------------------------------------------------------------------------------------------------------------




#-------------------------------------------------------------------------------------------------------------------------------
#	PREPARING THE SCRIPT: 1 entry must be entered by the user and is explained below
#-------------------------------------------------------------------------------------------------------------------------------


# ** < Please enter the full path to the cloned '/rvexcaliber' directory (example below): > **


path_to_rvexcaliber=""


# Example: "/genetics/rvexcaliber"


#-------------------------------------------------------------------------------------------------------------------------------




#-------------------------------------------------------------------------------------------------------------------------------
# Check script preparation
#-------------------------------------------------------------------------------------------------------------------------------


if [[ ! -d ${path_to_rvexcaliber} ]]; then

    echo -e "\n"
    echo "Error: The full path to the /RV-EXCALIBER directory that was downloaded locally does not exist"
    echo "Please ensure that this path is correctly defined above, under the "PREPARING THE SCRIPT" section"
    exit 1

fi


#-------------------------------------------------------------------------------------------------------------------------------




#-------------------------------------------------------------------------------------------------------------------------------
# Defaulted directories are provided below (do not change)
#-------------------------------------------------------------------------------------------------------------------------------


scripts="${path_to_rvexcaliber}/scripts"


#-------------------------------------------------------------------------------------------------------------------------------




# Directory and name variables


outdir=$1
internal_testing_dataset=$2
internal_ranking_dataset=$3


# Filter variables


gnomAD_MAF_threshold=$4
MCAP_threshold=$5
eth=$6
coverage=$7
adjust_type=$8


# Default changes and variables


outdir=$(echo ${outdir} | sed 's/\/$//g')

function stop_if_error {

    exit_status=$?

    if [ ${exit_status} -ne 0 ]; then

        echo -e "\n"

        echo "Please review error"

        exit "${exit_status}"

    fi

}

assoc_mod="rvexcaliber_base"




# Script start


# Check 'user entered' input variables


if [[ $# = 8 ]]; then

    dig='^[0-9]+(\.[0-9]+)?$'
    digw='^[0-9]+$'

    # Check variables for gnomAD MAF thresholds

    gnomAD_MAF_threshold_count=$(echo ${gnomAD_MAF_threshold} | awk -F"," '{for(i=1;i<=NF;i++) print $i}' | awk 'END { print NR }')

    gnomAD_MAF_threshold_dec=$(echo ${gnomAD_MAF_threshold} | awk -F"," '{for(i=1;i<=NF;i++) print $i}' | grep -wc ".")

    gnomAD_MAF_threshold_check=$(echo ${gnomAD_MAF_threshold} | sed -e 's/,//g' -e 's/\.//g')

    gnomAD_MAF_ge_lim=$(echo ${gnomAD_MAF_threshold} | awk -F"," '{for(i=1;i<=NF;i++) print $i}' | awk '$1>=1' | awk 'END { print NR}')

    gnomAD_MAF_le_lim=$(echo ${gnomAD_MAF_threshold} | awk -F"," '{for(i=1;i<=NF;i++) print $i}' | awk '$1<1e-05' | awk 'END { print NR}')

    # Check variables for MCAP thresholds

    MCAP_threshold_count=$(echo ${MCAP_threshold} | awk -F"," '{for(i=1;i<=NF;i++) print $i}' | awk 'END { print NR }')

    MCAP_threshold_dec=$(echo ${MCAP_threshold} | awk -F"," '{for(i=1;i<=NF;i++) print $i}' | grep -wc ".")

    MCAP_threshold_check=$(echo ${MCAP_threshold} | sed -e 's/,//g' -e 's/\.//g')

    MCAP_ge_lim=$(echo ${MCAP_threshold} | awk -F"," '{for(i=1;i<=NF;i++) print $i}' | awk '$1>=1' | awk 'END { print NR}')

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

    elif [[ ${gnomAD_MAF_ge_lim} -ge 1 ]]; then

        echo -e "\n"
        echo "Error: Ensure all gnomAD MAF thresholds (i.e. 'gnomAD_MAF_threshold') are less than 1"
        echo "See 'Explanation of input arguments' for further information"
        exit 1

    elif [[ ${gnomAD_MAF_le_lim} -ge 1 ]]; then

        echo -e "\n"
        echo "Error: Ensure all gnomAD MAF thresholds (i.e. 'gnomAD_MAF_threshold') are greater than 1e-05"
        echo "See 'Explanation of input arguments' for further information"
        exit 1

    elif [[ ! ${MCAP_threshold_check} =~ ${dig} ]] || [[ ${MCAP_threshold_count} -ne ${MCAP_threshold_dec} ]]; then

        echo -e "\n"
        echo "Error: Ensure that the MCAP thresholds (i.e. 'MCAP_threshold') are numeric and *comma separated* if using multiple MCAP thresholds"
        echo "See 'Explanation of input arguments' for further information"
        exit 1

    elif [[ ${MCAP_ge_lim} -ge 1 ]]; then

        echo -e "\n"
        echo "Error: Ensure all MCAP thresholds (i.e. 'MCAP_threshold') are less than 1"
        echo "See 'Explanation of input arguments' for further information"
        exit 1

    elif [[ ${eth_count} -ne ${eth_check} ]]; then

        echo -e "\n"
        echo "Error: Ensure that the comparator gnomAD ethnicity variable (i.e. 'eth') is one of < "nfe", "afr", "sas", "eas", "amr" >; if there are more than one 'eth' variable, please ensure they are **comma separated**"
        echo "See 'Explanation of input arguments' for further information"
        exit 1

    elif [[ ${eth_count} -ne ${eth_check} ]]; then

        echo -e "\n"
        echo "Error: Ensure that the comparator gnomAD ethnicity variable (i.e. 'eth') is one of < "nfe", "afr", "sas", "eas", "amr" >"
        echo "or"
        echo "If you wish ot generate a weighted comparator gnomAD allele frequency based on 2 or more individual 'eth' variables, please ensure they are separated by an underscore (_)"
        echo "See 'Explanation of input arguments' for further information"
        exit 1

    elif [[ ${eth_count} -gt 1 ]] && [[ ${eth_check} -ne ${eth_num_check} ]]; then

        echo -e "\n"
        echo "Error: If you wish ot generate a weighted comparator gnomAD allele frequency based on 2 or more individual comparator gnomAD ethnicity variables (i.e. 'eth'), please ensure there is a percentage weight is assigned to each 'eth' variable"
        echo "See 'Explanation of input arguments' for further information"
        exit 1

    elif [[ ${eth_count} -gt 1 ]] && [[ ${eth_num_sum} -ne 100 ]]; then

        echo -e "\n"
        echo "Error: If you wish ot generate a weighted comparator gnomAD allele frequency based on 2 or more individual comparator gnomAD ethnicity variables (i.e. 'eth'), please ensure that their corresponding percentage weights add to 100"
        echo "See 'Explanation of input arguments' for further information"
        exit 1

    elif [[ ! ${coverage} = "hcc" ]] && [[ ! ${coverage} = "rcc" ]]; then

        echo -e "\n"
        echo "Error: Ensure that the 'site-based intersection' variable (i.e. 'coverage') is one of < "hcc" or "rcc" >"
        echo "See 'Explanation of input arguments' for further information"
        exit 1

    elif [[ ! ${adjust_type} = "noadjust" ]] && [[ ! ${adjust_type} = "fulladjust" ]]; then

        echo -e "\n"
        echo "Error: Ensure that the 'type of adjustment' variable (i.e. adjust_type)  is one of < "noadjust" or "fulladjust" >"
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


    # Check if rare variant burden matrices exist based on filtering criteria


    for MAF_MCAP in ${MAF_MCAP}; do

        internal_testing_RVBurdenMatrix=${outdir}/${internal_testing_dataset}_RVBurdenMatrix_${MAF_MCAP}_${eth}_${coverage}.txt.gz
        gnomAD_internal_testing_RVBurdenMatrix=${outdir}/gnomAD_RVBurdenMatrix_${internal_testing_dataset}_${MAF_MCAP}_${eth}_${coverage}.txt.gz

        if [[ ! ${internal_ranking_dataset} = NA ]]; then

            internal_ranking_RVBurdenMatrix=${outdir}/${internal_ranking_dataset}_RVBurdenMatrix_${MAF_MCAP}_${eth}_${coverage}.txt.gz
            gnomAD_internal_ranking_RVBurdenMatrix=${outdir}/gnomAD_RVBurdenMatrix_${internal_ranking_dataset}_${MAF_MCAP}_${eth}_${coverage}.txt.gz

        elif [[ ${internal_ranking_dataset} = NA ]]; then

            internal_ranking_RVBurdenMatrix=NA
            gnomAD_internal_ranking_RVBurdenMatrix=NA

        fi

        if [[ ! -f ${internal_testing_RVBurdenMatrix} ]]; then

            echo -e "\n"
            echo "Error: Cannot locate the rare variant burden matrix for the internal testing dataset"
            echo -e "\n"
            echo "Please ensure that:"
            echo -e "\n"
            echo "1) the specified naming criteria (i.e. 'internal_testing_dataset') and all filtering criteria (i.e. 'gnomAD_MAF_theshold', 'MCAP_threshold', 'eth', and 'coverage') are consistent with what was used in the 'get_RVBurdenMatrix_internal_rvexcaliber.clean.sh' script"
            echo "and"
            echo "2) the specified output directory (i.e. 'outdir') contains the rare variant burden matrix for the internal testing dataset that was generated in the 'get_RVBurdenMatrix_internal_rvexcaliber.clean.sh' script"
            exit 1

        elif [[ ! -f ${gnomAD_internal_testing_RVBurdenMatrix} ]]; then

            echo -e "\n"
            echo "Error: Cannot locate the gnomAD rare variant burden matrix that corresponds to the internal testing dataset"
            echo -e "\n"
            echo "Please ensure that:"
            echo -e "\n"
            echo "1) the specified naming criteria (i.e. 'internal_testing_dataset') and all filtering criteria (i.e. 'gnomAD_MAF_theshold', 'MCAP_threshold', 'eth', and 'coverage') are consistent with what was used in the 'get_RVBurdenMatrix_internal_rvexcaliber.clean.sh'  and 'get_RVBurdenMatrix_gnomAD_rvexcaliber.clean.sh' scripts"
            echo "and"
            echo "2) the specified output directory (i.e. 'outdir') contains the rare variant burden matrix for the internal testing dataset that was generated in the 'get_RVBurdenMatrix_internal_rvexcaliber.clean.sh' and 'get_RVBurdenMatrix_internal_rvexcaliber.clean.sh' scripts"
            exit 1

        elif [[ ! ${internal_ranking_dataset} = NA ]]; then

            if [[ ! -f ${internal_ranking_RVBurdenMatrix} ]]; then

                echo -e "\n"
                echo "Error: Cannot locate the rare variant burden matrix for the internal ranking dataset"
                echo -e "\n"
                echo "If you do not have an internal ranking dataset, please assign 'internal_ranking_dataset' to NA"
                echo -e "\n"
                echo "If you do have an internal ranking dataset, please ensure that:"
                echo -e "\n"
                echo "1) the specified naming criteria (i.e. 'internal_ranking_dataset') and all filtering criteria (i.e. 'gnomAD_MAF_theshold', 'MCAP_threshold', 'eth', and 'coverage') are consistent with what was used in the 'get_RVBurdenMatrix_internal_rvexcaliber.clean.sh' script"
                echo "and"
                echo "2) the specified output directory (i.e. 'outdir') contains the rare variant burden matrix for the internal testing dataset that was generated in the 'get_RVBurdenMatrix_internal_rvexcaliber.clean.sh' script"
                exit 1

            elif [[ ! -f ${gnomAD_internal_ranking_RVBurdenMatrix} ]]; then

                echo -e "\n"
                echo "Error: Cannot locate the gnomAD rare variant burden matrix that corresponds to the internal ranking dataset"
                echo -e "\n"
                echo "Please ensure that:"
                echo -e "\n"
                echo "1) the specified naming criteria (i.e. 'internal_ranking_dataset') and all filtering criteria (i.e. 'gnomAD_MAF_theshold', 'MCAP_threshold', 'eth', and 'coverage') are consistent with what was used in the 'get_RVBurdenMatrix_internal_rvexcaliber.clean.sh'  and 'get_RVBurdenMatrix_gnomAD_rvexcaliber.clean.sh' scripts"
                echo "and"
                echo "2) the specified output directory (i.e. 'outdir') contains the rare variant burden matrix for the internal testing dataset that was generated in the 'get_RVBurdenMatrix_internal_rvexcaliber.clean.sh' and 'get_RVBurdenMatrix_internal_rvexcaliber.clean.sh' scripts"
                exit 1

            fi

        fi


        # Generate Summary Associations for the user-defined internal testing or ranking
        # datasets using rvexcaliber


        function rvexcaliber_association_test {(

            set -e

            if [[ ${adjust_type} = "fulladjust" ]]; then

                adjust_type="iCFadjust"

            elif [[ ${adjust_type} = "noadjust" ]]; then

                :

            fi

            filter_names_mod="${MAF_MCAP}_${eth}_${coverage}_${assoc_mod}.txt"

            filter_names="${MAF_MCAP}_${eth}_${coverage}_${adjust_type}_${assoc_mod}.txt"

            function cond_state {

                if [[ ${adjust_type} = "iCFadjust" ]]; then

                    echo "Adjusting gnomAD allele counts with the iCF"

                elif [[ ${adjust_type} = "noadjust" ]]; then

                    echo "No adjustment on gnomAD allele counts"

                fi

            }

            if [[ -f ${internal_testing_RVBurdenMatrix} ]] && [[ ${internal_ranking_RVBurdenMatrix} = NA ]]; then

                echo -e "\n"
                echo "Using summary associations from the default ranking dataset"
                echo -e "\n"
                echo "and generating summary associations from user-defined testing dataset"
                echo -e "\n"
                cond_state

                stage="stage1"

                gnomAD_RVBurdenMatrix=${gnomAD_internal_testing_RVBurdenMatrix}

                Rscript ${scripts}/Rscripts/get_iCF_association_test_rvexcaliber.R ${path_to_rvexcaliber} ${outdir} ${internal_testing_dataset} NA ${stage} ${MAF_MCAP} ${eth} ${coverage} ${adjust_type} ${internal_testing_RVBurdenMatrix} ${internal_ranking_RVBurdenMatrix} ${gnomAD_RVBurdenMatrix}

                gzip -f ${outdir}/rvexcaliber_testing_${internal_testing_dataset}_RVBurdenMatrix_pruned_${filter_names_mod}

                gzip -f ${outdir}/rvexcaliber_gnomAD_RVBurdenMatrix_pruned_${internal_testing_dataset}_${filter_names}

                if [[ ! -f ${outdir}/rvexcaliber_testing_${internal_testing_dataset}_SummaryAssociations_${filter_names} ]]; then

                    echo -e "\n"
                    echo "Error: The Summmary Association files for the internal testing dataset: ${internal_testing_dataset}, using adjustment type: ${adjust_type} have not been generated"
                    exit 1

                elif [[ -f ${outdir}/rvexcaliber_testing_${internal_testing_dataset}_SummaryAssociations_${filter_names} ]]; then

                    echo -e "\n"
                    echo "Summary Assocciations have been generated for the internal testing dataset: ${internal_testing_dataset}, using adjustment type: ${adjust_type}"
                    echo -e "\n"
                    echo "** Successful completion **"

                fi

            elif [[ -f ${internal_testing_RVBurdenMatrix} ]] && [[ ! ${internal_ranking_RVBurdenMatrix} = NA ]] && [[ -f ${internal_ranking_RVBurdenMatrix} ]]; then

                echo -e "\n"
                echo "Generating summary associations from user-defined ranking dataset"

                echo -e "\n"
                cond_state

                stage="stage0"

                gnomAD_RVBurdenMatrix=${gnomAD_internal_ranking_RVBurdenMatrix}

                Rscript ${scripts}/Rscripts/get_iCF_association_test_rvexcaliber.R ${path_to_rvexcaliber} ${outdir} ${internal_testing_dataset} ${internal_ranking_dataset} ${stage} ${MAF_MCAP} ${eth} ${coverage} ${adjust_type} ${internal_testing_RVBurdenMatrix} ${internal_ranking_RVBurdenMatrix} ${gnomAD_RVBurdenMatrix}

                echo -e "\n"
                echo "Generating summary associations from user-defined testing dataset"

                echo -e "\n"
                cond_state

                stage="stage2"

                gnomAD_RVBurdenMatrix=${gnomAD_internal_testing_RVBurdenMatrix}

                Rscript ${scripts}/Rscripts/get_iCF_association_test_rvexcaliber.R ${path_to_rvexcaliber} ${outdir} ${internal_testing_dataset} ${internal_ranking_dataset} ${stage} ${MAF_MCAP} ${eth} ${coverage} ${adjust_type} ${internal_testing_RVBurdenMatrix} ${internal_ranking_RVBurdenMatrix} ${gnomAD_RVBurdenMatrix}

                gzip -f ${outdir}/rvexcaliber_testing_${internal_testing_dataset}_RVBurdenMatrix_pruned_${filter_names_mod}

                gzip -f ${outdir}/rvexcaliber_gnomAD_RVBurdenMatrix_pruned_${internal_testing_dataset}_${filter_names}

                if [[ -f ${outdir}/rvexcaliber_testing_${internal_testing_dataset}_SummaryAssociations_${filter_names} ]]; then

                    echo -e "\n"
                    echo "Summary Assocciations have been generated for the internal testing dataset: ${internal_testing_dataset}, using adjustment type: ${adjust_type}"
                    echo -e "\n"
                    echo "** Successful completion **"

                else

                    echo -e "\n"
                    echo "Error: The Summmary Association files for the internal testing dataset: ${internal_testing_dataset}, using adjustment type: ${adjust_type} have not been generated"
                    exit 1

                fi

            fi

            if [[ ${adjust_type} = "iCFadjust" ]]; then

                echo -e "\n"
                echo "Now adjusting gnomAD allele counts with the gCF"

                if [[ ${internal_ranking_RVBurdenMatrix} = NA ]]; then

                    ranking_dataset=${path_to_rvexcaliber}/default_ranking_dataset/rvexcaliber_ranking_MIGen_exome_controls_0.001_0.025_nfe_high_coverage_coding_iCFadjust_rvexcaliber_base.txt

                elif [[ ! ${internal_ranking_RVBurdenMatrix} = NA ]]; then

                    ranking_dataset=${outdir}/rvexcaliber_ranking_${internal_ranking_dataset}_SummaryAssociations_${filter_names}

                fi

                Rscript ${scripts}/Rscripts/get_gCF_association_test_rvexcaliber.R ${outdir} ${internal_testing_dataset} ${MAF_MCAP} ${eth} ${coverage} ${adjust_type} ${outdir}/rvexcaliber_testing_${internal_testing_dataset}_RVBurdenMatrix_pruned_${filter_names_mod}.gz ${outdir}/rvexcaliber_gnomAD_RVBurdenMatrix_pruned_${internal_testing_dataset}_${filter_names}.gz ${outdir}/rvexcaliber_testing_${internal_testing_dataset}_SummaryAssociations_${filter_names} ${ranking_dataset}

                rm -f ${outdir}/rvexcaliber_ranking_${internal_ranking_dataset}_SummaryAssociations_${filter_names}

                rm -f ${outdir}/rvexcaliber_testing_${internal_testing_dataset}_SummaryAssociations_${filter_names}

                adjust_type="iCFgCFadjust"

                filter_names="${MAF_MCAP}_${eth}_${coverage}_${adjust_type}_${assoc_mod}.txt"

                gzip -f ${outdir}/rvexcaliber_gnomAD_RVBurdenMatrix_pruned_${filter_names}

                if [[ -f ${outdir}/rvexcaliber_testing_${internal_testing_dataset}_SummaryAssociations_${filter_names} ]] && [[ -f ${outdir}/rvexcaliber_testing_${internal_testing_dataset}_SummaryAssociations_allele_filter_${filter_names} ]]; then

                    echo -e "\n"
                    echo "Summary Assocciations have been generated for the internal testing dataset: ${internal_testing_dataset}, using adjustment type: ${adjust_type}"
                    echo -e "\n"
                    echo  "** Successful completion **"

                else

                    echo -e "\n"
                    echo "Error: The summary Assocciations have been generated for the internal testing dataset: ${internal_testing_dataset}, using adjustment type: ${adjust_type}"
                    exit 1

                fi

            else

                :

            fi

        )}
        rvexcaliber_association_test
        stop_if_error

    done
    reassign

else

    echo -e "\n"
    echo "Error: Incorrect number of input arguments"

fi
