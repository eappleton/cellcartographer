#nProc=8
#timestamp() {  STAMP=`date +%Y_%m_%d_%H_%M_%S`; }
#timestamp
#out="results_"$STAMP
#mkdir ${out}
#set -e
cd tba/model_training

################## FULL EXAMPLE ##################

###  extracts the sequences of of open chromatin regions. (minutes) ###
# creates ./output/neutrophil.fasta
mkdir ${out}
echo $CT
#find ./ -name '*'
python extract_sequences.py ../../resources/data/bed/${STAMP}/${CT}.bed ${genome} ./${out}/${CT}.fasta

###  generates GC matched background coordinates (minutes) ###
# creates./ouput/background.bed and ./output/background.fasta
python generate_background_coordinates.py ../../resources/data/bed/${STAMP}/${CT}.bed ${genome} ./${out}

### rename background files (seconds) ###
mv ./${out}/background.bed ./${out}/${CT}_background.bed
mv ./${out}/background.fasta ./${out}/${CT}_background.fasta

### calculates motif scores for open chromatin regions and background regions ###
# creates several files (you just need 2 - I use the others for other analyses):
# ./output/neutrophil_motif_scores.tsv 
# ./output/ neutrophil_motif_starts.tsv 
# ./output/ neutrophil_background_motif_scores.tsv 
# ./output/ neutrophil_background_motif_starts.tsv  
# ./output/ neutrophil_combined_features.tsv  <--- you'll need this one
# ./output/ neutrophil_labels.txt <-- you'll need this one
# this is a bit of a bottle neck step depending on how much computing resources you can throw
# for 100000 open chromatin regions (+100000 background sites) this should take 1 min for each motif
# there are 200 motifs and so the process will take 200 minutes if using just one processor core
# I parallelized the code and you can specific how many cores you'd like to use using the -num_proc parameter. I used 40 and so this process took about 5 minutes
python create_features.py ./${out}/${CT}.fasta ./${out}/${CT}_background.fasta ./${out} ./../default_motifs/*motif -num_proc ${nProc}

### trains model. ###
# This step will tell you motifs positively correlated and negatively correlated with open chromatin
# this step takes ~30 seconds per cross validation sets. The more cross validation sets the more robust the model outputs are - the default is 5
# some of the code I use for cross validation is out of date and so this may break if the versions don't match up. This would take me a few hours to fix and I intend to get to it soon
# creates:
# ./ouput/coefficients.tsv, model coefficients for each cross validation set
# ./output/performance.tsv, aucROC and precision for each cross validation set 
python train_classifier.py ./${out}/${CT}_combined_features.tsv ./${out}/${CT}_labels.txt ./${out}/

# rename model output
mv ./${out}/coefficients.tsv ./${out}/${CT}_coefficient.tsv
mv ./${out}/performance.tsv ./${out}/${CT}_performance.tsv
 
### Runs likelihood ratio test for each motif ###
# Tests how much removing a motif from the model hurts the model performance 
# This is basically an insilico mutagenesis screen. Reports unnormalized p-values
# takes 30 seconds to test each motif and so 100 minutes
# this code can be parallelized and I plan to - this would take me a couple of days
# creates:
# ./output/significance.tsv
python calc_feature_significance.py -num_procs ${nProc} ./${out}/${CT}_combined_features.tsv ./${out}/${CT}_labels.txt ./${out}/

### rename results of likelihood_ratio test ###
mv ./${out}/significance.tsv ./${out}/${CT}_significance.tsv

#Rename motifs in coefficients and significance files to gene names
#python map_motif_to_gene.py ./${out}/${CT}_coefficient.tsv ./${out}/${CT}_coefficient_gene.tsv ./../default_motifs/*jaspar
#python map_motif_to_gene.py ./${out}/${CT}_significance.tsv ./${out}/${CT}_significance_gene.tsv ./../default_motifs/*jaspar
python annotate_results_with_genes.py ./${out}/${CT}_coefficient.tsv ./${out}/${CT}_coefficient_gene.tsv
python annotate_results_with_genes.py ./${out}/${CT}_significance.tsv ./${out}/${CT}_significance_gene.tsv

python rank_order.py ./${out}/${CT}_coefficient_gene.tsv ./${out}/${CT}_significance_gene.tsv ./${out}/${CT}_gene_rank.tsv

#Move the results to the resources folder
mkdir -p ./../../resources/data/tba/${out}
mv ./${out}/* ./../../resources/data/tba/${out}
