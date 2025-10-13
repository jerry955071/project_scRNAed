# This workflow performs orthologous clustering of 15 plant species
# The 15 plant species include:
#   P. trichocarpa (v4.1; Phytozome v13), 
#   E. grandis (v2.0; Phytozome v13), 
#   Physcomitrium patens (v6.1; Phytozome v13), 
#   Selaginella moellendorffii (v1.0; Phytozome v13), 
#   T. aralioides (https://doi.org/10.1093/gigascience/giz136), 
#   L. chinense (v1.0, TreeGenes),
#   Marchantia polymorpha (MpTak1v5.1; MarpolBase), 
#   Pinus taeda (v2.01; TreeGenes), 
#   Gnetum montanum (A genome for gnetophytes and early evolution of seed plants; Dryad Digital Repository), 
#   Amborella trichopoda (Ensembl Genomes, release 49), 
#   Oryza sativa (Ensembl Genomes, release 49), 
#   Arabidopsis thaliana (Araport11; The Arabidopsis Information Resource),
#   Coffea canephora (Coffee Genome Hub), and 
#   Solanum lycopersicum (ITAG4.0; Sol Genomics Network).

#   C. lanceolata:
#   "Haplotype-resolved de novo genome assemblies of four coniferous tree species", Journal of Forest Research, 2023;:1-7
#   "Chinese fir genome and the evolution of gymnosperms"


# ==== 1. Manually downloading plant genomes ====
# Phytozome Batch 1: Ptr, Egr, Ppa, and Smo
rule phytozome_batch1:
    output:
        directory("outputs/Orthologous_clustering/references/Phytozome_batch1")
    log:
        "logs/Orthologous_clustering/phytozome_batch1.log"
    shell:
        """
        bash scripts/phytozome-batch1.sh {output} 1> {log} 2> {log}
        """

# T. aralioides (https://doi.org/10.1093/gigascience/giz136)
rule tar:
    output:
        directory("outputs/Orthologous_clustering/references/Tar")
    log:
        "logs/Orthologous_clustering/tar.log"
    shell:
        """
        wd=$(pwd)
        mkdir -p {output}
        cd {output}
        wget https://s3.ap-northeast-1.wasabisys.com/gigadb-datasets/live/pub/10.5524/100001_101000/100657/Analyses/Structure/Trochodendronb_aralioides_chromosomes_pasa2.longest.filter.pep \
            1> $wd/{log} \
            2> $wd/{log}
        """

# L. chinense (v1.0, TreeGenes)
rule lch:
    output:
        directory("outputs/Orthologous_clustering/references/Lch")
    log:
        "logs/Orthologous_clustering/lch.log"
    shell:
        """
        wd=$(pwd)
        mkdir -p {output}
        cd {output}
        wget https://treegenesdb.org/FTP/Genomes/Lich/v1.0/annotation/Lich.1_0.pep.fa.gz \
            1> $wd/{log} \
            2> $wd/{log}
        """

# Marchantia polymorpha (MpTak1v5.1; MarpolBase)
rule mpo:
    output:
        directory("outputs/Orthologous_clustering/references/Mpo")
    log:
        "logs/Orthologous_clustering/mpo.log"
    shell:
        """
        wd=$(pwd)
        mkdir -p {output}
        cd {output}
        wget https://marchantia.info/data/MpTak_v7.1_standard_genome/MpTak_v7.1.protein.primary.fa \
            1> $wd/{log} \
            2> $wd/{log}
        """

# Pinus taeda (v2.01; TreeGenes)
rule pta:
    output:
        directory("outputs/Orthologous_clustering/references/Pta")
    log:
        "logs/Orthologous_clustering/pta.log"
    shell:
        """
        wd=$(pwd)
        mkdir -p {output}
        cd {output}
        wget https://treegenesdb.org/FTP/Genomes/Pita/v2.01/annotation/Pita.2_01.pep.fa.gz \
            1> $wd/{log} \
            2> $wd/{log}
        """

# Gnetum montanum (A genome for gnetophytes and early evolution of seed plants; Dryad Digital Repository)
rule gmo:
    output:
        directory("outputs/Orthologous_clustering/references/Gmo")
    log:
        "logs/Orthologous_clustering/gmo.log"
    shell:
        """
        mkdir -p {output}
        echo 'Direct download from Dryad using web browser.' > {log}
        """

#   Amborella trichopoda (Ensembl Genomes, release 61), 
#   Oryza sativa (Ensembl Genomes, release 61), 
rule ensembl_batch1:
    output:
        atr=directory("outputs/Orthologous_clustering/references/Atr"),
        osa=directory("outputs/Orthologous_clustering/references/Osa")
    log:
        "logs/Orthologous_clustering/ensembl_batch1.log"
    shell:
        """
        wd=$(pwd)
        mkdir -p {output.atr} {output.osa}
        
        # download Atr peptides
        cd {output.atr}
        wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/fasta/amborella_trichopoda/pep/Amborella_trichopoda.AMTR1.0.pep.all.fa.gz \
            1> $wd/{log} \
            2> $wd/{log}

        # download Osa peptides
        cd $wd/{output.osa}
        wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/fasta/oryza_sativa/pep/Oryza_sativa.IRGSP-1.0.pep.all.fa.gz \
            1>> $wd/{log} \
            2>> $wd/{log}
        """

# Arabidopsis thaliana
rule ath:
    output:
        directory("outputs/Orthologous_clustering/references/Ath")
    log:
        "logs/Orthologous_clustering/ath.log"
    shell:
        """
        wd=$(pwd)
        mkdir -p {output}
        cd {output}
        wget https://www.arabidopsis.org/api/download-files/download?filePath=Proteins/Araport11_protein_lists/Araport11_pep_20250411.gz \
            1> $wd/{log} \
            2> $wd/{log}
        """

# Coffea canephora
rule cca:
    output:
        directory("outputs/Orthologous_clustering/references/Cca")
    log:
        "logs/Orthologous_clustering/cca.log"
    shell:
        """
        mkdir -p {output}
        echo 'Direct download from Coffee Genome Hub using web browser.' > {log}
        """

#  Solanum lycopersicum (ITAG4.0; Sol Genomics Network).
rule sly:
    output:
        directory("outputs/Orthologous_clustering/references/Sly")
    log:
        "logs/Orthologous_clustering/sly.log"
    shell:
        """
        wd=$(pwd)
        mkdir -p {output}
        cd {output}
        wget https://solgenomics.net/ftp/tomato_genome/annotation/ITAG4.0_release/ITAG4.0_proteins.fasta \
            1> $wd/{log} \
            2> $wd/{log}
        """

# C. lanceolata
rule cla:
    output:
        directory("outputs/Orthologous_clustering/references/Cla")
    log:
        "logs/Orthologous_clustering/cla.log"
    shell:
        """
        # prepare output folder
        wd=$(pwd)
        mkdir -p {output}
        cd {output}
        
        # download
        $wd/src/ncbi-datasets/datasets download \
            genome accession GCA_027924645.1 \
                --dehydrated \
                --include gff3,rna,cds,protein,genome,seq-report \
                --filename ncbi_dataset.zip
            2> $wd/{log} \
            1> $wd/{log}
        
        # unzip 
        unzip ncbi_dataset.zip

        # rehydrate
        $wd/src/ncbi-datasets/datasets rehydrate \
            --directory . \
            2> $wd/{log} \
            1> $wd/{log}
        """


# get longest protein
rule get_longest_protein:        
    input:
        directory("outputs/Orthologous_clustering/references/Lch"),
    output:
        
# DOWNLOADED
# 48354 outputs/Orthologous_clustering/references/Ath/Araport11_pep_20250411
# 27313 outputs/Orthologous_clustering/references/Atr/Amborella_trichopoda.AMTR1.0.pep.all.fa
# 25574 outputs/Orthologous_clustering/references/Cca/coffea_canephora.protein.faa
# 36349 outputs/Orthologous_clustering/references/Phytozome_batch1/Phytozome/PhytozomeV12/Egrandis/annotation/Egrandis_297_v2.0.protein_primaryTranscriptOnly.fa
# 27491 outputs/Orthologous_clustering/references/Gmo/Gmm.final.pep
# 35269 outputs/Orthologous_clustering/references/Lch/Lich.1_0.pep.fa
# 18007 outputs/Orthologous_clustering/references/Mpo/MpTak_v7.1.protein.primary.fa

# 35328 outputs/Orthologous_clustering/references/Tar/Trochodendronb_aralioides_chromosomes_pasa2.longest.filter.pep
# 22273 outputs/Orthologous_clustering/references/Phytozome_batch1/Phytozome/PhytozomeV9/Smoellendorffii/annotation/Smoellendorffii_91_protein_primaryTranscriptOnly.fa
# 34699 outputs/Orthologous_clustering/references/Phytozome_batch1/Phytozome/PhytozomeV13/Ptrichocarpa/v4.1/annotation/Ptrichocarpa_533_v4.1.protein_primaryTranscriptOnly.fa
# 36043 outputs/Orthologous_clustering/references/Phytozome_batch1/Phytozome/PhytozomeV13/Ppatens/v6.1/annotation/Ppatens_870_v6.1.protein_primaryTranscriptOnly.fa
# 51751 outputs/Orthologous_clustering/references/Pta/Pita.2_01.pep.fa
# 42582 outputs/Orthologous_clustering/references/Osa/Oryza_sativa.IRGSP-1.0.pep.all.fa



# NOTE: Number of records in each protein fasta (from previous study)
# 27654 Athaliana_all.protein.fa
# 27313 Atrichopoda_all_protein.fa
# 25574 Ccanephora_all_protein.fa
# 37225 Clanceolata_all_protein.fa
# 36349 Egrandis_all_protein.fa
# 27491 Gmontanum_all_protein.fa
# 35269 Lchinense_all_protein.fa
# 24750 Mpolymorpha_all_protein.fa
# 35775 Osativa_all_protein.fa
# 32926 Ppatens_all_protein.fa
# 51751 Ptaeda_all_protein.fa
# 34699 Ptrichocarpa_all_protein.fa
# 34075 Slycopersicu_all_proteins.fa
# 22285 Smoellendorffii_all_protein.fa
# 35328 Taralioides_longest_protein.fa
