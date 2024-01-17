### header ###
__author__ = "Evan Appleton"
__email__ = "evan.appleton@gmail.com"

# Import dependencies
import requests, json, os
#import json, os

# Function for building search terms
# Inputs cell type, sequencing technology, and cell source
# Outputs a search string
def buildENCODESearchString(CT, seqTech, source, org):
	
	#searchString="https://www.encodeproject.org/search/?frame=embedded&format=json&type=experiment&limit=all"
	searchString="https://www.encodeproject.org/search/?frame=embedded&format=json&type=Experiment&status=released&perturbed=false"

	# First define list of all types of organisms, assay categories, assays, and sources
	assayCategories = ["DNA binding", "Transcription"]
	assays = ["ChIP-seq","DNase-seq","polyA plus RNA-seq","shRNA RNA-seq","total RNA-seq","eCLIP","RNA microarray","DNAme array","WGBS","RRBS","small RNA-seq","microRNA-seq","ATAC-seq","RAMPAGE","RNA Bind-n-Seq","genotyping array","CAGE","single cell RNA-seq","Repli-seq","microRNA counts","siRNA-seq","CRISPRi RNA-seq","RIP-seq","MRE-seq","Repli-chip","MeDIP-seq","CRISPR RNA-seq","genetic modification DNase-seq","ChIA-PET","FAIRE-seq","Hi-C","PAS-seq","RIP-chip","ployA depleted RNA-seq","RNA-PET","whole genome sequencing assay","genotyping HTS","MS-MS","5C","TAB-seq","iCLIP","DNA-PET","5%27 RLm RACE","MNase-seq,Switchgear,Circulome-seq"]
	sources = ["immortalized cell line","tissue","whole organism","primary cell","in vitro differentiated cells","stem cell","cell-free sample","induced pluripotent stem cell line"]
	organisms = ["Homo sapiens","Mus musculus","Drosophila melanogaster","Caenorhabditis elegans","Drosophila pseudoobscura","Drosophila simulans","Drosophila mojavensis","Drosophila ananassae","Drosophila virilis","Drosophila yakuba"]

	# First determine if the assay listed is an assay or a caterory
	if seqTech.lower() in [x.lower() for x in assayCategories]:		

		#If it is, append this filter to the search URL
		lAssayCategories=[x.lower() for x in assayCategories]
		ind=lAssayCategories.index(seqTech.lower())
		#searchString+=("&assay_slims=" + assayCategories[ind].replace(" ", "+"))

	# Filter by assay
	if seqTech.lower() in [x.lower() for x in assays]:

                #If it is, append this filter to the search URL
                lAssays=[x.lower() for x in assays]
                ind=lAssays.index(seqTech.lower())
                searchString+=("&assay_title=" + assays[ind].replace(" ", "+"))

	# Filter by organism
	if org.lower() in [x.lower() for x in organisms]:

                #If it is, append this filter to the search URL
                lOrgs=[x.lower() for x in organisms]
                ind=lOrgs.index(org.lower())
                searchString+=("&replicates.library.biosample.donor.organism.scientific_name=" + organisms[ind].replace(" ", "+"))

	# Filter by source
	if source.lower() in [x.lower() for x in sources]:

                #If it is, append this filter to the search URL
                lSources=[x.lower() for x in sources]
                ind=lSources.index(source.lower())
                #searchString+=("&biosample_ontology.classification=" + sources[ind].replace(" ", "+"))

	# Search cell type
	if bool(CT):
		searchString+=("&biosample_ontology.term_name=" + CT.replace(" ", "+"))

	print(searchString)
	return searchString

# Get search extension string for designating files downloaded from search
def getSearchExt(CT, seqTech, source, org):

	searchExt="_"
#	if bool(CT):
#		searchExt+=(CT.replace(" ", "_") + "_")
	if bool(seqTech):
                searchExt+=(seqTech.replace(" ", "_").replace(",", "").replace(".", "") + "_")
	if bool(source):
                searchExt+=(source.replace(" ", "_").replace(",", "").replace(".", "") + "_")	
	if bool(org):
                searchExt+=(org.replace(" ", "_").replace(",", "").replace(".", "") + "_")

	return searchExt

# Main function for performing the ENCODE search
if __name__ == "__main__":

	execPath = os.path.dirname(os.path.realpath(__file__))
	
	# Force return from the server in JSON format
	HEADERS = {'accept': 'application/json'}
	CT_file = open('search_results_ct.txt', 'w')
	size_file = open('search_results_size.txt', 'w')
	URLs_file = open('search_results.txt', 'w')
	prefix_file = open('search_results_prefixes.txt', 'w') 
	searchExt_file = open('search_results_extensions.txt', 'w')
	
	# Make the output directory
	outDir=os.getenv('outDir')
	if not os.path.exists(outDir):
		os.makedirs(outDir)

	# Retrieve key terms from bash environment to build search strings if the direct search strings are not provided
	keyTerms=os.getenv('keyTerms').split('|')
	for key in keyTerms:
		
		i=0
		print(key)
		tokens=key.split(';')
		CT=tokens[0]
		seqTech=tokens[1]
		source=tokens[2]
		org=tokens[3]
		URL=buildENCODESearchString(CT, seqTech, source, org)
		search_ext=getSearchExt(CT, seqTech, source, org)	

		# GET the search result and extract the JSON response as a python dict
		#response = requests.get(URL, headers=HEADERS)
		response = requests.get(URL)
		response_json_dict = response.json()
		#print(response_json_dict)
		
		# Get all the files in this response
		graph = response_json_dict['@graph']
		files = [d['files'] for d in graph]
		
		# Loop through all the files and collect all of the fastq files, and print the hrefs to a file
		for exptSet in files:
			for f in exptSet:
				if f.get('file_type') == 'fastq':
					if f.get('replicate').get('experiment').get('status') == 'released':
	
						#If the experiment is a ChIP-seq experiment, only consider the targets (H3K27Ac, H3K4me1, H3K4me3)
						if "ChIP-seq" in seqTech:
							target=f.get('replicate').get('experiment').get('target').get('label')
							allowedTargets=["H3K27Ac","H3K4me1","H3K4me3"]
							if target not in allowedTargets:
								continue

						#Make a directory for each cell type				
						#ct=f.get('replicate').get('experiment').get('biosample_summary').replace(" ", "_").replace(",", "").replace(".", "").replace("(","").replace(")","")
						ct=CT.replace(" ", "_").replace(",", "").replace(".", "").replace("(","").replace(")","")
						ctDir=outDir+"/"+ct
						if not os.path.exists(ctDir):
        	        				os.makedirs(ctDir)

						i+=1
						size=f.get('file_size')
						size_file.write(str(f.get('file_size'))+'\n')
						CT_file.write(ct+'\n')
						URLs_file.write(f.get('href')+'\n')
						prefix_file.write(f.get('href').split('/')[-1].split('.')[0]+'\n')
						searchExt_file.write("_"+ct+search_ext+str(i)+".fastq.gz"+'\n')

	#Make a new directory for the search results and move file there
	os.rename("search_results.txt", outDir + "/search_results.txt")
	os.rename("search_results_size.txt", outDir + "/search_results_size.txt")
	os.rename("search_results_ct.txt", outDir + "/search_results_ct.txt")
	os.rename("search_results_prefixes.txt", outDir + "/search_results_prefixes.txt")
	os.rename("search_results_extensions.txt", outDir + "/search_results_extensions.txt")
