# Script_Bash_IFB 
## Collection de scripts bash pour analyse RNAseq sur cluster
adaptés de scripts initialement créés par [Pierre Poulain](https://github.com/pierrepo)  dans le cadre du [DU Omiques](https://odf.u-paris.fr/fr/offre-de-formation/diplome-d-universite-1/sciences-technologies-sante-STS/du-creation-analyse-et-valorisation-de-donnees-biologiques-omiques-DUSCAVD_117.html).


**00-script_cluster_genome_Index.sh** : _script pour indexation génome_

**00-script_concatenate-fastaq.sh** :  _script pour rassembler fichiers FastaQ_

**01-script_cluster_paired_Trim.sh** : _script pour analyse fichiers FastaQ en paired end avec filtrage_

**01-script_cluster_paired_noTrim.sh** : _script pour analyse fichiers FastaQ en paired end sans filtrage_

**01-script_cluster_single_Trim.sh** : _script pour analyse fichiers FastaQ en single end avec filtrage_

**01-script_cluster_single_noTrim.sh** : _script pour analyse fichiers FastaQ en single end sans filtrage_

**02-script_cluster_Count_aggreg.sh** : _script pour rassembler les tables de comptage_

**02-script_cluster_multiQC.sh** : _script pour rassembler les analyses fastqc_

**TruSeq3-SE.fa** : fichier nécessaire pour filtrage des séquences d'adaptateurs avec trimmomatic ([source](https://github.com/timflutre/trimmomatic/tree/master/adapters))

### Licence

![](CC-BY-SA.png)

Ce contenu est mis à disposition selon les termes de la licence [Creative Commons Attribution - Partage dans les Mêmes Conditions 4.0 International](https://creativecommons.org/licenses/by-sa/4.0/deed.fr) (CC BY-SA 4.0). Consultez le fichier [LICENSE](LICENSE) pour plus de détails.

This content is released under the [Creative Commons Attribution-ShareAlike 4.0 ](https://creativecommons.org/licenses/by-sa/4.0/deed.en) (CC BY-SA 4.0) license. See the bundled [LICENSE](LICENSE) file for details.
