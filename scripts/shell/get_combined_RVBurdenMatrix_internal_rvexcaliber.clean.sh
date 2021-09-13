#!/bin/bash

#===============================================================================================================================
# Title: get_combined_RVBurdenMatrix_internal_rvexcaliber.clean.sh
# Generate a combined rare variant burden matrix for your internal testing and/or ranking dataset
#===============================================================================================================================
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

# Copyright 2021 Ricky Lali, Michael Chong, Arghavan Omidi, Pedrum Mohammadi-Shemirani, Ann Le, Edward Cui, and Guillaume Pare
#===============================================================================================================================




#-------------------------------------------------------------------------------------------------------------------------------
#	EXPLANATION OF INPUT ARGUMENTS:
#-------------------------------------------------------------------------------------------------------------------------------

# outdir

# Description:     full path to the directory to where all output files will be directed

# Expected input:  character input

#-------------------------------------------------------------------------------------------------------------------------------

# combined_combined_internal_dataset

# Description:     desired name of your combined internal ** testing ** or ** ranking ** dataset

# Expected input:  character input

#-------------------------------------------------------------------------------------------------------------------------------

# gnomAD_MAF_threshold

# Description:     upper-limit minor allele frequency threshold from ** gnomAD ** used to filter variants

# Expected input:  comma seperated numeric value for multiple thresholds (e.g. 0.001,0.01)
#                  OR
#                  single numeric value for a single threshold (e.g. 0.001)

#-------------------------------------------------------------------------------------------------------------------------------

# MCAP_threshold

# Description:     lower-limit threshold MCAP score used to filter nonsynonymous SNVs

# Expected input:  comma seperated numeric value for multiple thresholds (e.g. 0.009,0.025)
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

# Description:     conduct site-based intersection of input combined_internal_dataset with high coverage coding (hcc) regions in gnomAD

# Expected input:  can be one of < "hcc", "rcc" >
#                  "hcc" = high-coverage coding; "rcc" = regular-coverage coding

# Notes:           high coverage coding regions are defined as sites in gnomAD that have >20X coverage in >90% of individuals

#-------------------------------------------------------------------------------------------------------------------------------

# internal_RVBurdenMatrix_filelist

# Description:     a concatenated list of the full paths to the per-chromosome rare variant burden matrices for your internal
#                  testing and/or ranking dataset

# Expected input:  can be one of < "Y", "N" >
#                  "Y" = yes; "N" = no

# Notes:           high coverage coding regions are defined as sites in gnomAD that have >20X coverage in >90% of individuals

#-------------------------------------------------------------------------------------------------------------------------------




#-------------------------------------------------------------------------------------------------------------------------------
#	PREPARING THE SCRIPT: 1 entry must be entered by the user and are explained below
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




# Directory and name variables


outdir=$1
combined_internal_dataset=$2


# Fiter variables


gnomAD_MAF_threshold=$3
MCAP_threshold=$4
eth=$5
coverage=$6


# File list variables


internal_RVBurdenMatrix_filelist=$7


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


# Check 'user entered' input variables


if [[ $# = 7 ]]; then

    dig='^[0-9]+(\.[0-9]+)?$'
    digw='^[0-9]+$'

    # Check variables for gnomAD MAF thresholds

    gnomAD_MAF_threshold_count=$(echo ${gnomAD_MAF_threshold} | awk -F"," '{for(i=1;i<=NF;i++) print $i}' | awk 'END { print NR }')

    gnomAD_MAF_threshold_dec=$(echo ${gnomAD_MAF_threshold} | awk -F"," '{for(i=1;i<=NF;i++) print $i}' | grep -wc ".")

    gnomAD_MAF_threshold_check=$(echo ${gnomAD_MAF_threshold} | sed -e 's/,//g' -e 's/\.//g')

    gnomAD_MAF_ge_1=$(echo ${gnomAD_MAF_threshold} | awk -F"," '{for(i=1;i<=NF;i++) print $i}' | awk '$1>=1' | awk 'END { print NR}')

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
        echo "Error:  If you wish ot generate a weighted comparator gnomAD allele frequency based on 2 or more individual comparator gnomAD ethnicity variables (i.e. 'eth'), please ensure that their corresponding percentage weights add to 100"
        echo "See 'Explanation of input arguments' for further information"
        exit 1

    elif [[ ! ${coverage} = "hcc" ]] && [[ ! ${coverage} = "rcc" ]]; then

        echo -e "\n"
        echo "Error: Ensure that the 'high coverage coding  intersection variable' (i.e. 'coverage') is one of < "hcc" or "rcc" >"
        echo "See 'Explanation of input arguments' for further information"
        exit 1

    elif [[ ! -f ${internal_RVBurdenMatrix_filelist} ]]; then

       echo -e "\n"
       echo "Error: Ensure that the concatenated list of the full paths to the per-chromosome rare variant burden matrices for your internal testing and/or ranking dataset (i.e. 'internal_RVBurdenMatrix_filelist') is present"
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
       sleep 5s

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


    for MAF_MCAP in ${MAF_MCAP}; do

        function sub_combine_RVBurdenMatrix_internal {(

            set -e

            rm -f ${outdir}/${combined_internal_dataset}_RVBurdenMatrix_${MAF_MCAP}_${eth}_${coverage}*

            Rscript ${path_to_rvexcaliber}/scripts/Rscripts/sub_combine_RVBurdenMatrix_internal_rvexcaliber.R ${outdir} ${combined_internal_dataset} ${MAF_MCAP} ${eth} ${coverage} ${internal_RVBurdenMatrix_filelist}

            gzip -f ${outdir}/${combined_internal_dataset}_RVBurdenMatrix_${MAF_MCAP}_${eth}_${coverage}.txt

        )}
        sub_combine_RVBurdenMatrix_internal
        stop_if_error

        if [[ -f ${outdir}/${combined_internal_dataset}_RVBurdenMatrix_${MAF_MCAP}_${eth}_${coverage}.txt.gz ]]; then

          echo -e "\n"
          echo "A combined rare variant burden matrix has been generated for ${combined_internal_dataset} at a MAF_MCAP: ${MAF_MCAP} and coverage: ${coverage}"
          echo -e "\n"
          echo "** Successful completion **"

        else

          echo -e "\n"
          echo "A combined rare variant burden matrix has not been generated for ${combined_internal_dataset} at a MAF_MCAP: ${MAF_MCAP} and coverage: ${coverage}"

        fi

    done
    reassign

else

  echo -e "\n"
  echo "Error: Incorrect number of input arguments"

fi
