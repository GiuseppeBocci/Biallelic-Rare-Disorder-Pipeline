import json
import sys

json_file = sys.argv[1]
annotation_sources = ["MIM_morbid"] #, "Orphanet"] # , "Cancer_Gene_Census"]
selected_freq = 5/10_000

#[36,"transcript_consequences",3,"phenotypes",58,"source"]       "Cancer_Gene_Census"
#[36,"transcript_consequences",3,"phenotypes",59,"source"]       "MIM_morbid"
#[36,"transcript_consequences",3,"phenotypes",60,"source"]       "Orphanet"

with open(json_file) as file:
	res = json.load(file)

# Trim json

sel = []

for snp in res:
	is_relevant = False
	freq = '-'
	id_rs = '-'
	phenotypes = []
	hgvsp = []
	gene_ids = []
	terms = []
	biotypes = []
	if "transcript_consequences" in snp:
		for trc in snp["transcript_consequences"]:
			if trc["impact"] == "HIGH" or trc["impact"] == "MODERATE":
				if "sift_score" in trc and "polypen_prediction" in trc:
					if not(trc["sift_score"] < 0.05 and trc["polypen_prediction"] > 0.446):
						continue
				is_relevant = True
				if "biotype" in trc:
					biotypes.append(trc["biotype"])
				else:
					biotypes.append('-')
				# print(trc)
				if "phenotypes" in trc:
					phenotypes.append([phe["phenotype"] for phe in trc["phenotypes"] if phe["source"] in annotation_sources])
				else:
					phenotypes.append([])
				hgvsp.append(trc["hgvsp"] if "hgvsp" in trc else '-')
				gene_ids.append(trc["gene_id"] if "gene_id" in trc else '-')
				terms.append(trc["consequence_terms"])
	if is_relevant and "colocated_variants" in snp:
		for cv in snp["colocated_variants"]:
			if "frequencies" in cv:
				freqs = [cv["frequencies"][all]["gnomade"] for all in cv["frequencies"]]
				if all([all_f <= selected_freq for all_f in freqs]):
					freq = ''.join(map(str, freqs))
					id_rs = cv['id']
					# print(cv)
				else:
					is_relevant = False
				# print(cv)

	if is_relevant:
		for i in range(len(gene_ids)):
			sel.append(f"{snp['input']}\t{snp['allele_string']}\t{snp['strand']}\t{id_rs}\t{freq}\t{' '.join(terms[i])}\t{biotypes[i]}\t{hgvsp[i]}@{gene_ids[i]}\t{'; '.join(phenotypes[i])}")

print("INPUT\tALLELES\tSTRAND\tID\tFREQUENCIES\tCONSEQUENCES\tBIOTYPE\tHGSV@GENE_ID\tPHENOTYPES")
for s in sel:
	print(s)
