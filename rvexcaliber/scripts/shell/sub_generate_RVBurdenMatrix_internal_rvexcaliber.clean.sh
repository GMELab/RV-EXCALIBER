#!/bin/bash


#===============================================================================================================================
# Title: sub_generate_RVBurdenMatrix_internal_rvexcaliber.clean.sh
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





# Directory and name variables


indir=$1
outdir=$2
internal_dataset=$3


# Filter variables


MAF_MCAP=$4
eth=$5
coverage=$6


# Miscellaneous variables


job_interval=$7
plink=$8


# Script start


# Check 'user entered' input variables


if [[ $# = 8 ]]; then

    function cleanup {

        rm -rf ${outdir}/temp_gene_out_${internal_dataset}_${MAF_MCAP}

        rm -rf ${outdir}/logfiles_${internal_dataset}_${MAF_MCAP}

    }
    cleanup &
    wait

    function create_jobs_log {

        startpoint=0

        rm -rf ${outdir}/logfiles_${internal_dataset}_${MAF_MCAP}

        mkdir -p ${outdir}/logfiles_${internal_dataset}_${MAF_MCAP}

        geneCount=$(awk 'BEGIN {FS=OFS=" "} { # switch to tab delim
              print $(NF-2)
        }' ${outdir}/${internal_dataset}_variant_list_extracted_${MAF_MCAP}_${eth}_${coverage}.txt |
        sort -u |
        awk 'END {
            print NR
        }')

        for ((a=${startpoint};a<=${geneCount};a+=${job_interval})); do

            echo -e ${a}"\t"$((${a}+${job_interval})) >> ${outdir}/logfiles_${internal_dataset}_${MAF_MCAP}/${internal_dataset}_${MAF_MCAP}_jobs.log

        done

        ll=$(awk 'END { print NR}' ${outdir}/logfiles_${internal_dataset}_${MAF_MCAP}/${internal_dataset}_${MAF_MCAP}_jobs.log)

        awk -v ll="$ll" -v geneCount="$geneCount" 'BEGIN { OFS="\t" } {

            if (NR==ll) $2=geneCount

        }1' ${outdir}/logfiles_${internal_dataset}_${MAF_MCAP}/${internal_dataset}_${MAF_MCAP}_jobs.log |

        awk 'BEGIN {FS=OFS="\t"} {

            if ($1!=$2)

                print $0

        }' > ${outdir}/logfiles_${internal_dataset}_${MAF_MCAP}/${internal_dataset}_${MAF_MCAP}_jobs_fin.log
    }
    create_jobs_log &
    wait

    function getGeneMatrix {

        if [[ -d ${outdir}/temp_gene_out_${internal_dataset}_${MAF_MCAP} ]]; then

            echo "Directory made"
            #rm -rf ${outdir}/temp_gene_out_${internal_dataset}_${MAF_MCAP}/${internal_dataset}_${MAF_MCAP}_genelist.txt

        else

            mkdir ${outdir}/temp_gene_out_${internal_dataset}_${MAF_MCAP}

        fi

        ll=$(awk 'END { print NR }' ${outdir}/logfiles_${internal_dataset}_${MAF_MCAP}/${internal_dataset}_${MAF_MCAP}_jobs_fin.log)

        for j in $(seq 1 1 ${ll}); do

        f1=$(awk 'NR=="'${j}'" { print $1 }' ${outdir}/logfiles_${internal_dataset}_${MAF_MCAP}/${internal_dataset}_${MAF_MCAP}_jobs_fin.log)

        f2=$(awk 'NR=="'${j}'" { print $2 }' ${outdir}/logfiles_${internal_dataset}_${MAF_MCAP}/${internal_dataset}_${MAF_MCAP}_jobs_fin.log)

        for i in $(awk '{ print $(NF-2) }' ${outdir}/${internal_dataset}_variant_list_extracted_${MAF_MCAP}_${eth}_${coverage}.txt |
            sort -u |
            awk -v f1="$f1" -v f2="$f2" '{ if (NR>f1 && NR<=f2)
            print $0
            }'); do

            echo ${i} >> ${outdir}/temp_gene_out_${internal_dataset}_${MAF_MCAP}/${internal_dataset}_${MAF_MCAP}_genelist.txt

            awk '$(NF-2)=="'${i}'"' ${outdir}/${internal_dataset}_variant_list_extracted_${MAF_MCAP}_${eth}_${coverage}.txt > ${outdir}/temp_gene_out_${internal_dataset}_${MAF_MCAP}/${i}_${MAF_MCAP}_scoring.txt &

        done
        wait

            for i in $(awk '{ print $(NF-2) }' ${outdir}/${internal_dataset}_variant_list_extracted_${MAF_MCAP}_${eth}_${coverage}.txt |
                sort -u |
                awk -v f1="$f1" -v f2="$f2" '{ if (NR>f1 && NR<=f2)
                print $0
                }'); do

                ${plink} --noweb \
                \
                --silent \
                \
                --bfile ${outdir}/${internal_dataset}_preprocessed \
                \
                --keep-allele-order \
                \
                --score ${outdir}/temp_gene_out_${internal_dataset}_${MAF_MCAP}/${i}_${MAF_MCAP}_scoring.txt 1 2 3 sum \
                \
                --out ${outdir}/temp_gene_out_${internal_dataset}_${MAF_MCAP}/${i}_${MAF_MCAP} &

            done

        done

    }
    getGeneMatrix &
    wait

    function combineMatrix {

    awk 'BEGIN {FS=OFS=" "; print "Gene"} { print $1 }' ${indir}/${internal_dataset}.fam | paste -s -d' ' > ${outdir}/${internal_dataset}_RVBurdenMatrix_${MAF_MCAP}_${eth}_${coverage}.txt

    for Gene in $(cat ${outdir}/temp_gene_out_${internal_dataset}_${MAF_MCAP}/${internal_dataset}_${MAF_MCAP}_genelist.txt); do

        if [[ -f ${outdir}/temp_gene_out_${internal_dataset}_${MAF_MCAP}/${Gene}_${MAF_MCAP}.profile ]]; then

            awk 'BEGIN {FS=OFS=" "; print "'${Gene}'"} NR>1{

                print $5

            }' ${outdir}/temp_gene_out_${internal_dataset}_${MAF_MCAP}/${Gene}_${MAF_MCAP}.profile |

            paste -s -d' ' >> ${outdir}/${internal_dataset}_RVBurdenMatrix_${MAF_MCAP}_${eth}_${coverage}.txt

        else

            echo "${Gene} allele scoring not available"

        fi

    done

    }
    combineMatrix &
    wait

    function cleanup {

        rm -rf ${outdir}/temp_gene_out_${internal_dataset}_${MAF_MCAP}

        rm -rf ${outdir}/logfiles_${internal_dataset}_${MAF_MCAP}

    }
    cleanup &
    wait

else
    echo "Error: Incorrect number of input arguments"
fi
