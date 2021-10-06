#!/bin/bash

#===============================================================================================================================
# Title: get_ANNOVAR_annotations.clean.sh
# Extract a variants that meet both MAF and MCAP thresholds
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
#===============================================================================================================================




#-------------------------------------------------------------------------------------------------------------------------------
#	PREPARING THE SCRIPT: 1 entry must be entered by the user and are explained below
#-------------------------------------------------------------------------------------------------------------------------------


# ** < Please enter the full path the directory containing the ANNOVAR perl scripts (example below): > **


annovar=""


# Example: "/genetics/tools/annovar"


# The following scripts were used to obtain the ** necessary ** input files from ANNOVAR used by RV-EXCALIBER

# 1. The refGene database:
# perl ${annovar}/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene ${annovar}/humandb/

# 2. The gnomAD (version 2.1.1) allele frequency exome annotations
# perl ${annovar}/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad211_exome  ${annovar}/humandb/

# 3. The in silico protein-scoring annotations from dbnsfp
# perl ${annovar}/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp35c  ${annovar}/humandb/


#-------------------------------------------------------------------------------------------------------------------------------




#-------------------------------------------------------------------------------------------------------------------------------
# Check script preparation
#-------------------------------------------------------------------------------------------------------------------------------


if [[ ! -d ${annovar} ]]; then

    echo -e "\n"
    echo "Error: The full path to the /annovar directory does not exist"
    echo "Please ensure that this path is correctly defined above, under the "PREPARING THE SCRIPT" section"
    exit 1

fi


#-------------------------------------------------------------------------------------------------------------------------------




# 1. The refGene database:

perl ${annovar}/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene ${annovar}/humandb/

# 2. The gnomAD (version 2.1.1) allele frequency exome annotations

perl ${annovar}/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad211_exome  ${annovar}/humandb/

# 3. The in silico protein-scoring annotations from dbNSFP

perl ${annovar}/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp35c  ${annovar}/humandb/
