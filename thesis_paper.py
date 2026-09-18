# Umutcan Unaldi
# Saarland University
# Bioinformatics MSc
# Master Thesis

import os
import pandas as pd
import openpyxl
import warnings
import subprocess
import shutil
import requests
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
import re
from collections import Counter
import urllib.request

from Bio.PDB import PDBParser, PDBIO, NeighborSearch, Superimposer
from Bio.PDB.Model import Model
from Bio.PDB.Polypeptide import is_aa
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices

from Bio.SeqUtils import seq1, seq3

from modeller import Environ, log
from modeller.automodel import AutoModel, assess, autosched, refine
from modeller.parallel import job, local_slave

# Get rid of unknown extension warning
warnings.simplefilter("ignore", UserWarning)

# Uses the folder containing this script as the working directory
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# Links for the supplementary tables of the original study [38][39]
mmc2_url = (
    "https://ars.els-cdn.com/content/image/"
    "1-s2.0-S0092867416300435-mmc2.xlsx"
)
mmc3_url = (
    "https://ars.els-cdn.com/content/image/"
    "1-s2.0-S0092867416300435-mmc3.xlsx"
)

# Downloads the tables only when they are not already present
if not os.path.exists("mmc2.xlsx"):
    urllib.request.urlretrieve(mmc2_url, "mmc2.xlsx")

if not os.path.exists("mmc3.xlsx"):
    urllib.request.urlretrieve(mmc3_url, "mmc3.xlsx")

# mmc3.xlsx consists of all the interactors and isoforms
df_mmc3 = pd.read_excel("mmc3.xlsx", sheet_name="2B-Isoform PPIs")

# Two SEPT9 rows were stored as Excel dates instead of gene names
# These were also absent from the original analysis
df_mmc3 = df_mmc3[
    df_mmc3["Gene_Symbol"].apply(
        lambda gene: isinstance(gene, str)
    )
].reset_index(drop=True)

# Extract interactors
inters = df_mmc3.iloc[:, 6].unique()

# Unconventional naming by the authors is fixed
inters[inters == "LOC100288797"] = "TMEM239"

def search(gene_name):
    ''' Searches for inter proteins in UniProt database
    and returns information like accession number, gene name,
    synonyms of all found isos of that inter protein [1].'''

    # Filters: Human, and Swiss-Prot reviewed
    query = f'(gene_exact:{gene_name})+(organism_id:9606+AND+reviewed:true)'
    # Output format option is set to JSON with included isoform information
    url = (f"https://rest.uniprot.org/uniprotkb/search?query={query}"
           f"&format=json&fields=accession,gene_names,protein_name&includeIsoform=true")
    response = requests.get(url)

    return response.json().get("results", [])

def parse(entry):
    ''' Parses the entry retrieved from search function [2].'''
    all_names = []
    gene_data = entry.get('genes', [])
    for gene in gene_data:
        if 'geneName' in gene:
            all_names.append(gene['geneName']['value'])
        for syn in gene.get('synonyms', []):
            all_names.append(syn['value'])

    # Get the recommended protein name
    protein_name = (entry.get('proteinDescription', {}).
                    get('recommendedName', {}).
                    get('fullName', {}).
                    get('value', ''))
    return all_names, protein_name

def json(accession, fasta):
    ''' Fetches the results in JSON format for an inter protein
    and also returns the FASTA file for an inter protein.'''
    if fasta == False:
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    else:
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text if fasta else response.json()
    else:
        print(f"Couldn't fetch {'FASTA' if fasta else 'JSON'} for {accession}")
        return None

all_data = list()
# Cache is to prevent duplicates
cache = dict()
# Keeps record of the non-retrieved interactors
failed_genes = list()

if not os.path.exists("fasta_interactors"):
    os.makedirs("fasta_interactors")
    # Handles the interactor names used in research, and retrieves the information
    for research_name in inters:
        # Removes isoform information suffices (_ORF1, _ORF2, etc.)
        prefix_inter = research_name.split('_')[0]
        # Uses it from cache if it is already processed
        if prefix_inter in cache:
            entries = cache[prefix_inter]
        # Adds to cache
        else:
            entries = search(prefix_inter)
            cache[prefix_inter] = entries

        # Each isoform has different accession numbers (Q15645-1, Q15645-2, etc.)
        accession_numbers = list()
        for entry in entries:
            accession_numbers.append(entry['primaryAccession'])
        accession_numbers.sort()

        # Targets isoforms
        if '_ORF' in research_name:
            # It is failed if it doesn't find multiple accession numbers from the UniProt
            if len(accession_numbers) == 1 or len(accession_numbers) == 0:
                failed_genes.append(research_name)
                print(f"{research_name} failed because there are no isos in UniProt")
                continue

            # Extracts the last letter, which indicates the
            # isoform number (1 for _ORF1, 2 for _ORF2, etc.)
            iso_num = int(research_name[-1])
            # Assigns the suffices for isoforms
            # Q15645 and Q15645-1 always returns the first isoform
            if iso_num == 1:
                accession = accession_numbers[int(iso_num - 1)] + str(-1)
            else:
                accession = accession_numbers[int(iso_num - 1)]
        # It is failed if UniProt search didn't return any accession numbers
        else:
            if len(accession_numbers) == 0:
                failed_genes.append(research_name)
                continue
            accession = accession_numbers[0]
        print(f"Selected UniProt entry for {research_name}: {accession}")

        full_entry = json(accession, False)
        fasta = json(accession, True)

        synonyms, protein_name = parse(full_entry)

        if fasta:
            lines = fasta.strip().splitlines()

            # Checks if first line is a proper FASTA header
            if not lines or not lines[0].startswith(">"):
                print(f"{research_name} doesn't have a proper FASTA header")
                continue

            # Separates the header and the sequence
            header = lines[0]
            # Also, removes line breaks and spaces
            sequence = ''.join(lines[1:]).replace(" ", "")

            # Formats the header
            synonym_str = "_".join(synonyms)
            prot_name_clean = protein_name.replace(" ", "_")
            # New header includes all essential information
            customised_header = (f">{research_name}:"
                                 f"{prefix_inter}:"
                                 f"{accession}:"
                                 f"Synonyms-{synonym_str}:"
                                 f"{prot_name_clean}")

            # Ensures not to repeat writing the files for each run of the code (debugging)


            # Writes to files for each isoform that is retrieved
            output_path = os.path.join("fasta_interactors", f"{research_name}.fasta")
            print(f"Output path is: {output_path}")
            with open(output_path, "w") as f:
                f.write(customised_header + "\n")
                f.write(sequence + "\n")

        else:
            print(f"{research_name}, no FASTA information")

if failed_genes:
    print(f"Failed genes: {failed_genes}")

def write_iso_fasta(output_dir):
    '''This function writes fasta files for isos
    that are retrieved from the research paper.'''

    # mmc2.xlsx has all isoform FASTA information
    path = "mmc2.xlsx"
    sheet_name = "1B-Ref+Alt ORFs"
    df = pd.read_excel(path, sheet_name=sheet_name)

    # Extracts ORF column from the Excel and zip it with Iso ID
    iso_sequence_map = dict(zip(df["Isoform_ID"],
                                    df["Isoform__Open_Reading_Frame_Sequence"]))

    # Ensures not to repeat writing the files for each run of the code (debugging)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

        # Creates the FASTA files for each isoform
        for iso_name, sequence in iso_sequence_map.items():
            filename = os.path.join(output_dir, f"{iso_name}.fasta")
            with open(filename, "w") as fasta_file:
                fasta_file.write(f">{iso_name}\n{sequence}\n")

write_iso_fasta("fasta_isoforms")

input_fasta_isos = "fasta_isoforms"
output_translated_fasta = "translated_isoforms"

# Ensures not to repeat writing the files for each run of the code (debugging)
if not os.path.exists(output_translated_fasta):
    os.makedirs(output_translated_fasta)

    # Need to extract each frame
    frames = [1, 2, 3, -1, -2, -3]

    # Retrieves all the fasta files and their path
    for fasta_file in os.listdir(input_fasta_isos):
        if fasta_file.endswith(".fasta"):
            input_path = os.path.join(input_fasta_isos, fasta_file)
            output_path = os.path.join(output_translated_fasta, fasta_file)

            with open(output_path, "w") as outfile:
                for frame in frames:
                    # Translates in each frame with SeqKit [3]
                    cmd = f"seqkit translate --frame {frame} --line-width 0 {input_path}"
                    # Processs the code line and waits until each run finishes (subprocess)
                    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

                    # Tag headers with frame information
                    frame_tag = f"_frame{frame}"
                    for line in result.stdout.splitlines():
                        if line.startswith(">"):
                            base_header = line[1:]
                            new_header = f">{base_header}{frame_tag}"
                            outfile.write(new_header + "\n")
                        else:
                            outfile.write(line + "\n")

output_longest_translated_fasta = "translated_isoforms_longest"

def retrieve_sensible_fasta(sequence):
    ''' Returns the sequence before
    it hits the stop codon.'''
    if '*' in sequence:
        return sequence.split('*')[0]
    return sequence

# Ensures not to repeat writing the files for each run of the code (debugging)
if not os.path.exists(output_longest_translated_fasta):
    os.makedirs(output_longest_translated_fasta)

    for tr_fasta in os.listdir(output_translated_fasta):
        if tr_fasta.endswith(".fasta"):
            # Creates the paths for fasta files
            input_path = os.path.join(output_translated_fasta, tr_fasta)
            output_path = os.path.join(output_longest_translated_fasta, tr_fasta)

            # Iterate the frames into a list [4] [5]
            frames = list(SeqIO.parse(input_path, "fasta"))
            best_frame = None
            longest_seq = -1

            for frame in frames:
                shortened_seq = retrieve_sensible_fasta(str(frame.seq))
                # Iterates over all frames and keeps the longest seq
                if len(shortened_seq) > longest_seq:
                    longest_seq = len(shortened_seq)
                    best_frame = frame
                    best_frame.seq = best_frame.seq[:len(shortened_seq)]

            # Writes the longest frame
            if best_frame:
                with open(output_path, "w") as out_f:
                    writer = FastaWriter(out_f, wrap=None)
                    writer.write_file([best_frame])

# Assigns isoform ID, interactor ID and interaction information
iso_id = df_mmc3.iloc[:, 2]
inter_id = df_mmc3.iloc[:, 6]
inter_info = df_mmc3.iloc[:, 7]

inter_id[inter_id == "LOC100288797"] = "TMEM239"

# Maps all three information for pairs
map_iso_inter_excel = [[iso, inter, info] for iso, inter, info in
                     zip(iso_id, inter_id, inter_info)]

# Creates the dataframe for the pair map
df_pairs = pd.DataFrame(map_iso_inter_excel,
                        columns=['Iso_ID', 'Inter_ID', 'Interaction_Found'])
print("Interaction info in dataframe:")
print(df_pairs)


# Saves only the interaction data
df_pairs.to_csv("df_strctmn.csv", index=False)

# Gets rid of the missing values
df_pairs = df_pairs.dropna(subset=['Interaction_Found'])

# Adds a Base_Iso column to make it easy to select
df_pairs['Base_Iso'] = df_pairs['Iso_ID'].str.split('_').str[0]

print(
    f"\nDataframe after NAN values and missing inter proteins "
    f"are removed and Base_Iso column is added:\n {df_pairs}\n")

# Groups by Base_Iso name and Inter_ID to get
# cases where an isoform interacted the same with all
# interactor proteins [6]
conflicting_interactions = df_pairs.groupby(['Base_Iso', 'Inter_ID'])[
                               'Interaction_Found'].nunique() > 1

# Filters based on conflicting interactions
df_pairs_filtered = df_pairs[df_pairs.set_index(['Base_Iso', 'Inter_ID']).index.isin(
    conflicting_interactions[conflicting_interactions].index)]

print(f"Filtered dataframe:\n {df_pairs_filtered}\n")

# Some statistics:
print(f"\nUnique Iso Number: {df_pairs_filtered['Iso_ID'].nunique()}")
print(f"Unique Gene Number: {df_pairs_filtered['Base_Iso'].nunique()}")
print(f"Unique Inter Protein Number: {df_pairs_filtered['Inter_ID'].nunique()}\n")

# Merging and writing the pairs step
inter_dir = "fasta_interactors"
output_pairs_dir = "fasta_pairs"
iso_filtered_dir = "fasta_isoforms_filtered"
inter_filtered_dir = "fasta_interactors_filtered"

# Ensures not to repeat writing the files for each run of the code (debugging)
for path in [output_pairs_dir, iso_filtered_dir, inter_filtered_dir]:
    os.makedirs(path, exist_ok=True)

for _, row in df_pairs.iterrows():
    iso_id = row["Iso_ID"]
    inter_id = row["Inter_ID"]

    # File paths for the output files
    iso_fasta = os.path.join(output_longest_translated_fasta, f"{iso_id}.fasta")
    inter_fasta = os.path.join(inter_dir, f"{inter_id}.fasta")
    output_fasta = os.path.join(output_pairs_dir, f"{iso_id}__{inter_id}.fasta")
    iso_out = os.path.join(iso_filtered_dir, f"{iso_id}.fasta")
    inter_out = os.path.join(inter_filtered_dir, f"{inter_id}.fasta")

    if os.path.exists(iso_fasta) and os.path.exists(inter_fasta):
        with open(output_fasta, "w") as outfile:
            # Writes isoform FASTA
            with open(iso_fasta, "r") as infile:
                iso_data = infile.read().strip()
                outfile.write(iso_data + "\n")
            # Writes interactor FASTA
            with open(inter_fasta, "r") as infile:
                inter_data = infile.read().strip()
                outfile.write(inter_data + "\n")

        # Also writes individual isoform FASTA files
        if not os.path.exists(iso_out):
            with open(iso_out, "w") as f:
                f.write(iso_data + "\n")

        # And, writes individual interactor FASTA files
        if not os.path.exists(inter_out):
            with open(inter_out, "w") as f:
                f.write(inter_data + "\n")

    else:
        print(f"{iso_id} and {inter_id} are missing")

all_fasta = list()
for dir in [iso_filtered_dir, inter_filtered_dir]:
    for file in os.listdir(dir):
        if file.endswith(".fasta"):
            all_fasta.append(os.path.join(dir, file))

with open("structman_input.fasta", "w") as fw:
    for path in all_fasta:
        with open(path, "r") as fr:
            fw.write(fr.read())

structman_output_dir = "structman_output"
def swap_colon(pdb):
    '''Swaps chains in Structure Recommendation
    column in .protein_protein_interactions'''
    if ":" in pdb:
        l, r = pdb.rsplit(":", 1)
        if l and r:
            l, r = l[:-1] + r[0], l[-1] + r[1:]
        return l + ":" + r
    return pdb

# Needs the output file
if os.path.exists(structman_output_dir):
    temp = pd.read_csv(
        f'{structman_output_dir}/structman_input.protein_protein_interactions'
                                '.tsv', sep='\t')

    # Selects only the relevant columns
    df_strctmn = temp[['Input Protein ID A', 'Input Protein ID B', 'Structure Recommendation',
                            'Interaction Score']].copy()

    # Changes the column names
    df_strctmn.columns= ['Iso_ID', 'Inter_ID', 'Structure_Recommendation',
                               'Interaction_Score']

    # Mask that looks for cloned isoforms in the interactor column [7]
    swap_mask = df_strctmn['Inter_ID'].str.contains('frame')

    # Swaps the values when it finds an isoform in the interactor column
    df_strctmn.loc[swap_mask, ['Iso_ID', 'Inter_ID']] = (
        df_strctmn.loc[swap_mask, ['Inter_ID', 'Iso_ID']].values)

    # Change the Structure Recommendation structure (7ANK:B:C would be changed to
    # 7ANK:C:B if the column were swapped)
    df_strctmn.loc[swap_mask, 'Structure_Recommendation'] = \
        df_strctmn.loc[swap_mask, 'Structure_Recommendation'].apply(swap_colon)

    # Gets rid of the _frame* suffices [8]
    df_strctmn['Iso_ID'] = df_strctmn['Iso_ID'].astype(str).str.replace(
        r'_frame\d+$', '', regex=True)
    df_strctmn['Inter_ID'] = df_strctmn['Inter_ID'].astype(str).str.replace(
        r'_frame\d+$', '', regex=True)

    # Filters by the selected df_pairs [9]
    filter_pairs = set(
        df_pairs_filtered.apply(lambda row: frozenset([row['Iso_ID'], row['Inter_ID']]),
                                axis=1))
    ppi_filtered = df_strctmn[df_strctmn.apply(lambda row: frozenset(
        [row['Iso_ID'], row['Inter_ID']]) in filter_pairs,
                                                 axis=1)]

    # Merges the Y2H information
    sorted_ppi = ppi_filtered.merge(
        df_pairs_filtered[['Iso_ID', 'Inter_ID', 'Interaction_Found']],
        on=['Iso_ID', 'Inter_ID'],
        how='left'
    )

    # Reformats the order of the columns
    ppi_condensed = sorted_ppi[["Iso_ID", "Inter_ID", "Interaction_Found",
                                "Interaction_Score",
                                "Structure_Recommendation"]]

    ppi_condensed = ppi_condensed.sort_values(by='Iso_ID')

    ppi_condensed = ppi_condensed.drop_duplicates(subset=['Iso_ID', 'Inter_ID',
                                                    'Interaction_Score',
                                                    'Structure_Recommendation'])

    iso_id = ppi_condensed['Iso_ID'].tolist()

    # Gets the prefixes
    prefix_iso = [x[:-2] for x in iso_id]

    # Adds the protein name without the isoform-number suffix
    structman_df = ppi_condensed.copy()
    structman_df["Base_Iso"] = (
        structman_df["Iso_ID"].str.rsplit("_", n=1).str[0]
    )

    # Keeps interactors represented by at least two isoforms of the same protein
    non_strict_comparisons = structman_df.groupby(
        ["Base_Iso", "Inter_ID"]
    )

    non_strict_iso_count = (
        non_strict_comparisons["Iso_ID"].transform("nunique")
    )

    non_strict_df = structman_df[
        non_strict_iso_count >= 2
        ].copy()

    # Saves the broader result for possible later analyses
    non_strict_df.to_csv(
        f"{structman_output_dir}/structman_result_non_strict.tsv",
        sep="\t",
        index=False,
    )

    # The strict comparison also requires the same template and chain assignment
    comparison_columns = [
        "Base_Iso",
        "Inter_ID",
        "Structure_Recommendation",
    ]

    strict_comparisons = non_strict_df.groupby(comparison_columns)

    strict_iso_count = (
        strict_comparisons["Iso_ID"].transform("nunique")
    )

    inter_count = (
        strict_comparisons["Interaction_Found"].transform("nunique")
    )

    # Keeps templates shared by different isoforms with opposite Y2H results
    final_df = non_strict_df[
        (strict_iso_count >= 2) &
        (inter_count >= 2)
        ].copy()

    print(f"Final result:\n{final_df}")

    # Saves the strict result used for downstream analyses and modelling
    final_df.to_csv(
        f"{structman_output_dir}/structman_result.tsv",
        sep="\t",
        index=False,
    )
# ============= GLOBAL SETTINGS FOR DOWNSTREAM ANALYSIS ==============

# Set for M1 Macbook Air
threads = 6


# ============= MMSEQS2 TEMPLATE SEARCH ===============
tpl_dir = "template_search"
# Output file that's going to be used for downstream analysis
mmseqs_output = os.path.join(tpl_dir, "mmseqs_pairs.tsv")

# Files that are produced by the MMseqs2 search
pdb_fasta = os.path.join(tpl_dir, "pdb_seqres.fasta")
iso_hits_file = os.path.join(tpl_dir, "iso_vs_pdb.tsv")
inter_hits_file = os.path.join(tpl_dir, "int_vs_pdb.tsv")

# Column order used when the MMseqs2 results are written
mmseqs_cols = [
    "query", "target", "evalue", "bits", "pident", "alnlen",
    "qcov", "tcov", "qstart", "qend", "tstart", "tend",
]

# Performs the template search if the template_search folder doesn't exist
if not os.path.exists(tpl_dir):
    os.makedirs(tpl_dir)

    # Combines isoform FASTA files for template search [5]
    all_iso_out = os.path.join(tpl_dir, "all_isoforms.fasta")
    with open(all_iso_out, "w") as output:
        for filename in sorted(os.listdir(iso_filtered_dir)):
            iso_name = os.path.splitext(filename)[0]
            rec = SeqIO.read(
                    os.path.join(iso_filtered_dir, filename), "fasta"
            )
            output.write(f">{iso_name}\n{rec.seq}\n")

    # Combines interactor FASTA files for template search [5]
    all_inter_out = os.path.join(tpl_dir, "all_interactors.fasta")
    with open(all_inter_out, "w") as output:
        for filename in sorted(os.listdir(inter_filtered_dir)):
            inter_name = os.path.splitext(filename)[0]
            rec = SeqIO.read(
                os.path.join(inter_filtered_dir, filename), "fasta"
            )
            output.write(f">{inter_name}\n{rec.seq}\n")

    # Downloads the PDB chain sequences [15]
    url = "https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt"
    urllib.request.urlretrieve(url, pdb_fasta)

    # Creates the three MMseqs2 databases
    iso_db = os.path.join(tpl_dir, "isoformsDB")
    inter_db = os.path.join(tpl_dir, "interactorsDB")
    pdb_db = os.path.join(tpl_dir, "pdbChainsDB")

    # Creates MMseqs2 databases [17]. For example:
    # mmseqs createdb template_search/all_isoforms.fasta template_search/isoformsDB
    subprocess.run(
        ["mmseqs", "createdb", all_iso_out, iso_db], check=True
    )
    subprocess.run(
        ["mmseqs", "createdb", all_inter_out, inter_db], check=True
    )
    subprocess.run(
        ["mmseqs", "createdb", pdb_fasta, pdb_db], check=True
    )

    # Searches the isoforms and interactors against the PDB chain database [16] [17]
    tmp_dir = os.path.join(tpl_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)
    iso_results = os.path.join(tpl_dir, "iso_vs_pdb")
    inter_results = os.path.join(tpl_dir, "int_vs_pdb")

    # Example: mmseqs search template_search/isoformsDB template_search/pdbChainsDB
    # template_search/iso_vs_pdb template_search/tmp --threads 6 -s 7.5 -e 1e-3
    # --max-seqs 200
    subprocess.run([
        "mmseqs", "search", iso_db, pdb_db,
        iso_results, tmp_dir,
        "--threads", str(threads), "-s", "7.5", "-e", "1e-3", "--max-seqs", "200",
    ], check=True)

    # Example: mmseqs search template_search/interactorsDB template_search/pdbChainsDB
    # template_search/int_vs_pdb template_search/tmp --threads 6 -s 7.5 -e 1e-3
    # --max-seqs 200
    subprocess.run([
        "mmseqs", "search", inter_db, pdb_db,
        inter_results, tmp_dir,
        "--threads", str(threads), "-s", "7.5", "-e", "1e-3", "--max-seqs", "200",
    ], check=True)

    # Saves the search results as tables with the given column names (order matters)
    format_output = ",".join(mmseqs_cols)

    # Converts the MMseqs2 results into tables with the selected alignment fields [17]
    # Example: mmseqs convertalis template_search/isoformsDB template_search/pdbChainsDB
    # template_search/iso_vs_pdb template_search/iso_vs_pdb.tsv --format-output "query, target, etc."
    subprocess.run([
        "mmseqs", "convertalis", iso_db, pdb_db,
        iso_results, iso_hits_file,
        "--format-output", format_output,
    ], check=True)

    # Example: mmseqs convertalis template_search/interactorsDB template_search/pdbChainsDB
    # template_search/int_vs_pdb template_search/int_vs_pdb.tsv --format-output "query, target, etc."
    subprocess.run([
        "mmseqs", "convertalis", inter_db, pdb_db,
        inter_results, inter_hits_file,
        "--format-output", format_output,
    ], check=True)

# Reads the two result tables
iso_hits = pd.read_csv(
    iso_hits_file, sep="\t", names=mmseqs_cols
)
inter_hits = pd.read_csv(
    inter_hits_file, sep="\t", names=mmseqs_cols
)

# Separates identifiers such as 2r3v_A into a PDB ID and chain ID
# So that it'd be: Shared_PDB = 8RVE and Iso/Inter_Chain = AE
iso_tpl = iso_hits["target"].str.split("_", n=1)
inter_tpl = inter_hits["target"].str.split("_", n=1)
iso_hits["Shared_PDB"] = iso_tpl.str[0].str.upper()
iso_hits["Iso_Chain"] = iso_tpl.str[1]
inter_hits["Shared_PDB"] = inter_tpl.str[0].str.upper()
inter_hits["Inter_Chain"] = inter_tpl.str[1]
iso_hits["Iso_ID"] = iso_hits["query"]
inter_hits["Inter_ID"] = inter_hits["query"]

# Thresholds for relatively meaningful results
iso_hits = iso_hits[
    (iso_hits["pident"] >= 25) &
    (iso_hits["qcov"] >= 0.60) &
    (iso_hits["tcov"] >= 0.60)
].copy()
inter_hits = inter_hits[
    (inter_hits["pident"] >= 25) &
    (inter_hits["qcov"] >= 0.60) &
    (inter_hits["tcov"] >= 0.60)
].copy()

# Keeps only the best 20 hits for each query [19]
iso_hits = iso_hits.sort_values(
    ["Iso_ID", "bits", "pident", "qcov"],
    ascending=[True, False, False, False],
).groupby("Iso_ID", group_keys=False).head(20)

inter_hits = inter_hits.sort_values(
    ["Inter_ID", "bits", "pident", "qcov"],
    ascending=[True, False, False, False],
).groupby("Inter_ID", group_keys=False).head(20)

# Finds PDB entries shared by each experimental isoform-interactor pair [18]
selected_pairs = df_pairs_filtered[
    ["Iso_ID", "Inter_ID", "Interaction_Found", "Base_Iso"]
].drop_duplicates().copy()
# Merges the isoforms by Iso_ID
pair_tpls = selected_pairs.merge(
    iso_hits, on="Iso_ID", how="inner"
)
# Merges the interactors by Inter_ID and Shared_PDB
pair_tpls = pair_tpls.merge(
    inter_hits, on=["Inter_ID", "Shared_PDB"],
    how="inner", suffixes=("_iso", "_inter"),
)

# The isoform and interactor must use two different chains
pair_tpls = pair_tpls[
    pair_tpls["Iso_Chain"] != pair_tpls["Inter_Chain"]
].copy()

# Uses both MMseqs2 bit scores for ranking after the identity
# and coverage thresholds have been applied [17]
pair_tpls["Combined_Bit_Score"] = (
        pair_tpls["bits_iso"] +
        pair_tpls["bits_inter"]
)

# Used as a tie-breaker so that both sides of the template have good coverage
pair_tpls["Minimum_qcov"] = pair_tpls[
    ["qcov_iso", "qcov_inter"]
].min(axis=1)

# Treats A-B and B-A as the same physical chain pair during ranking
# to also preserve homomeric information if there is any
pair_tpls["Unordered_Chain_Pair"] = pair_tpls.apply(
    lambda row: ":".join(sorted([
        str(row["Iso_Chain"]), str(row["Inter_Chain"])
    ])),
    axis=1,
)

# 2XRF chain A needs a symmetry-generated copy to form its dimer
# Chains B and C are the biological dimer present in the PDB file [36] [37]
# This is applied before ranking so that A-B is not selected from tied hits
pair_tpls = pair_tpls[
    (pair_tpls["Shared_PDB"] != "2XRF") |
    (pair_tpls["Unordered_Chain_Pair"] == "B:C")
].copy()

# Checks whether both template chains have the same sequence for homomeric information
pdb_seqs = {}
for entry in SeqIO.parse(pdb_fasta, "fasta"):
    target_parts = entry.id.split("_", 1)
    if len(target_parts) == 2:
        pdb_seqs[(target_parts[0].upper(), target_parts[1])] = str(entry.seq)

homomeric_tpl = []
for _, row in pair_tpls.iterrows():
    iso_tpl_seq = pdb_seqs.get(
        (row["Shared_PDB"], row["Iso_Chain"])
    )
    inter_tpl_seq = pdb_seqs.get(
        (row["Shared_PDB"], row["Inter_Chain"])
    )
    homomeric_tpl.append(
        iso_tpl_seq is not None and
        # Both template sequences need to be equal for homomeric status
        iso_tpl_seq == inter_tpl_seq
    )
pair_tpls["Homomeric_Tpl"] = homomeric_tpl

# Finds template-chain assignments shared by different isoforms
# with opposite experimental interaction results
comparison_columns = [
    "Base_Iso", "Inter_ID", "Shared_PDB",
    "Iso_Chain", "Inter_Chain",
]

comparisons = pair_tpls.groupby(comparison_columns)
iso_count = comparisons["Iso_ID"].transform("nunique")
inter_count = comparisons["Interaction_Found"].transform("nunique")

# This filtering is performed before choosing the best template
# Otherwise, a common second-best template could be removed
comparable_tpls = pair_tpls[
    (iso_count >= 2) & (inter_count >= 2)
    ].copy()

# Scores every common template according to its weakest isoform result
# this prevents one very strong alignment from hiding a weaker alignment
# for the other isoform
comparisons = comparable_tpls.groupby(comparison_columns)

comparable_tpls["Comparison_Bit_Score"] = (
    comparisons["Combined_Bit_Score"].transform("min")
)

comparable_tpls["Comparison_qcov"] = (
    comparisons["Minimum_qcov"].transform("min")
)

# Ranks templates that can actually be used for the isoform comparison
# Bit score is used first, while query coverage is used as a tie-breaker
comparable_tpls = comparable_tpls.sort_values(
    [
        "Base_Iso", "Inter_ID",
        "Comparison_Bit_Score", "Comparison_qcov",
    ],
    ascending=[True, True, False, False],
)

# Selects the best shared PDB and physical chain pair
# A-B and B-A have the same Unordered_Chain_Pair
best_physical_pairs = comparable_tpls.drop_duplicates(
    ["Base_Iso", "Inter_ID"]
)[[
    "Base_Iso", "Inter_ID", "Shared_PDB",
    "Unordered_Chain_Pair", "Homomeric_Tpl",
]]

# Retrieves all isoform rows belonging to the selected physical chain pair
best_pair_tpls = comparable_tpls.merge(
    best_physical_pairs,
    on=[
        "Base_Iso", "Inter_ID", "Shared_PDB",
        "Unordered_Chain_Pair", "Homomeric_Tpl",
    ],
    how="inner",
)

# Both orientations are retained when the template-chain sequences
# are identical
homomeric_tpls = best_pair_tpls[
    best_pair_tpls["Homomeric_Tpl"]
].copy()

# For other templates, only the best directed orientation is retained
other_tpls = best_pair_tpls[
    ~best_pair_tpls["Homomeric_Tpl"]
].copy()

best_other_orientation = (
    other_tpls
    .sort_values(
        [
            "Base_Iso", "Inter_ID",
            "Comparison_Bit_Score", "Comparison_qcov",
        ],
        ascending=[True, True, False, False],
    )
    .drop_duplicates(["Base_Iso", "Inter_ID"])
    [comparison_columns]
)

other_tpls = other_tpls.merge(
    best_other_orientation,
    on=comparison_columns,
    how="inner",
)

# Combines both types of result
mmseqs_df = pd.concat(
    [homomeric_tpls, other_tpls],
    ignore_index=True,
).drop_duplicates([
    "Iso_ID", "Inter_ID", "Shared_PDB",
    "Iso_Chain", "Inter_Chain",
])


# Add similar column name to StructMAn : Structure_Recommendation
mmseqs_df["Structure_Recommendation"] = (
    mmseqs_df["Shared_PDB"] + ":" +
    mmseqs_df["Iso_Chain"] + ":" +
    mmseqs_df["Inter_Chain"]
)
# General table
output_columns = [
    "Iso_ID", "Inter_ID", "Interaction_Found", "Base_Iso",
    "Shared_PDB", "Iso_Chain", "Inter_Chain",
    "Iso_pident", "Iso_qcov", "Iso_tcov", "Iso_bits",
    "Inter_pident", "Inter_qcov", "Inter_tcov", "Inter_bits",
    "Combined_Bit_Score", "Comparison_Bit_Score",
    "Iso_target", "Inter_target", "Iso_evalue", "Inter_evalue",
    "Homomeric_Tpl", "Structure_Recommendation",
]
# Rename the column names for better representability
mmseqs_df = mmseqs_df.rename(columns={
    "pident_iso": "Iso_pident", "qcov_iso": "Iso_qcov",
    "tcov_iso": "Iso_tcov", "bits_iso": "Iso_bits",
    "pident_inter": "Inter_pident", "qcov_inter": "Inter_qcov",
    "tcov_inter": "Inter_tcov", "bits_inter": "Inter_bits",
    "target_iso": "Iso_target", "target_inter": "Inter_target",
    "evalue_iso": "Iso_evalue", "evalue_inter": "Inter_evalue",
})[output_columns]
mmseqs_df = mmseqs_df.sort_values([
    "Base_Iso", "Inter_ID", "Shared_PDB",
    "Iso_Chain", "Inter_Chain", "Iso_ID",
]).reset_index(drop=True)

mmseqs_df.to_csv(mmseqs_output, sep="\t", index=False)

print("MMseqs2 pairs:")
print(mmseqs_df)


# ============= CANDIDATE TEMPLATES ===============
# Prepares StructMAn candidates
structman_candidates = final_df.copy()

# Separates the template and two chains
structman_parts = structman_candidates[
    "Structure_Recommendation"
].str.rsplit(":", n=2, expand=True)

structman_candidates["Base_Tpl"] = (
    structman_parts[0].str.replace("_AU", "", regex=False)
)
structman_candidates["Iso_Chain"] = structman_parts[1]
structman_candidates["Inter_Chain"] = structman_parts[2]

# Uses the same template format for both sources
structman_candidates["Structure_Recommendation"] = (
    structman_candidates["Base_Tpl"] + ":" +
    structman_candidates["Iso_Chain"] + ":" +
    structman_candidates["Inter_Chain"]
)

structman_candidates["Source"] = "StructMAn"


# Prepares MMseqs2 candidates
mmseqs_candidates = mmseqs_df.copy()
mmseqs_candidates["Base_Tpl"] = mmseqs_candidates["Shared_PDB"]
mmseqs_candidates["Source"] = "MMseqs2"


# Columns needed for modelling
candidate_cols = [
    "Iso_ID", "Inter_ID", "Interaction_Found", "Source",
    "Iso_Chain", "Inter_Chain", "Structure_Recommendation",
    "Base_Iso", "Base_Tpl",
]

# Combines candidates from both sources
candidate_df = pd.concat(
    [
        structman_candidates[candidate_cols],
        mmseqs_candidates[candidate_cols],
    ],
    ignore_index=True,
)

# Columns that define one exact modelling candidate
same_candidate = [
    "Iso_ID", "Inter_ID", "Interaction_Found",
    "Iso_Chain", "Inter_Chain", "Structure_Recommendation",
    "Base_Iso", "Base_Tpl",
]

# Marks a candidate if the same template and chain assignment came from both sources
source_count = candidate_df.groupby(
    same_candidate
)["Source"].transform("nunique")

candidate_df.loc[source_count == 2, "Source"] = "Both"

# Removes repeated rows for the same modelling candidate [20]
candidate_df = candidate_df.drop_duplicates(same_candidate)

candidate_df = candidate_df.sort_values([
    "Base_Iso", "Inter_ID", "Base_Tpl",
    "Iso_Chain", "Inter_Chain", "Iso_ID",
]).reset_index(drop=True)

print("Candidate tpls:")
print(candidate_df)

# Removes selected chain pairs that do not represent biological interfaces [21] [22]
nefl_mask = (
    (candidate_df["Base_Iso"] == "NEFL") &
    (candidate_df["Inter_ID"] == "DES") &
    (candidate_df["Base_Tpl"] == "8RVE")
)

mef2a_mask = (
    (candidate_df["Base_Iso"] == "MEF2A") &
    (candidate_df["Inter_ID"] == "MEOX1") &
    (candidate_df["Base_Tpl"] == "6WC2") &
    (candidate_df["Iso_Chain"] == "D") &
    (candidate_df["Inter_Chain"] == "O")
)

# Stores the removed candidates for the supplementary tables
excluded_candidate_df = candidate_df[
    nefl_mask | mef2a_mask
].copy()

# Explains why the selected NEFL chains were not modelled
excluded_candidate_df.loc[
    excluded_candidate_df["Base_Iso"] == "NEFL",
    "Exclusion_Reason",
] = "selected chains do not form a physical interface"

# Explains why the selected MEF2A chains were not modelled
excluded_candidate_df.loc[
    excluded_candidate_df["Base_Iso"] == "MEF2A",
    "Exclusion_Reason",
] = "selected chains belong to different biological assemblies"

candidate_df = candidate_df[
    ~(nefl_mask | mef2a_mask)
].copy()

print("Candidate tpls after filtering:")
print(candidate_df)

candidate_df.to_csv(
    "candidate_templates.tsv",
    sep="\t",
    index=False,
)


# ============= TEMPLATE STRUCTURES ===============

pdb_dir = "pdb_files"
af_models_dir = "af_models"

os.makedirs(pdb_dir, exist_ok=True)

# Downloads the templates that are not already present [14]
for tpl_id in candidate_df["Base_Tpl"].drop_duplicates():
    if len(tpl_id) == 4:
        tpl_path = os.path.join(
            pdb_dir, f"{tpl_id}.pdb"
        )

        if not os.path.exists(tpl_path):
            urllib.request.urlretrieve(
                f"https://files.rcsb.org/download/{tpl_id}.pdb",
                tpl_path,
            )


# ============= MODELLER ===============


# Output directory
modeller_output_dir = "modeller_complex_runs"
os.makedirs(modeller_output_dir, exist_ok=True)

# The same isoform, template and chains only need to be modelled once
model_jobs = candidate_df.drop_duplicates([
    "Iso_ID", "Base_Tpl",
    "Iso_Chain", "Inter_Chain",
]).copy()

# Modified amino acids are mapped to one letter because of
# the usage of coordinate-derived sequences
aa_map = {
    "MSE": "M", "SEC": "U", "PYL": "O", "SEP": "S", "TPO": "T",
    "PTR": "Y", "CSO": "C", "CME": "C", "HYP": "P",
}

# Iterates through models and assigns isoform-interactor
# information to variables
for _, model_row in model_jobs.iterrows():
    iso_name = model_row["Iso_ID"]
    tpl_id = model_row["Base_Tpl"]
    iso_chain = model_row["Iso_Chain"]
    inter_chain = model_row["Inter_Chain"]

    # Sets and creates the folder names for each model
    # Example:
    # Isoform: UPP2_2 Template: 2XRF
    # Isoform Chain: A Interactor Chain: B
    # Run: UPP2_2__2XRF_A__partner_B
    run_name = (
        f"{iso_name}__{tpl_id}_{iso_chain}"
        f"__partner_{inter_chain}"
    )
    # "abspath" used because of changing working directory later on
    run_dir = os.path.abspath(
        os.path.join(modeller_output_dir, run_name)
    )

    # Skips the job if MODELLER models are already present
    # to save time for debugging
    if os.path.exists(run_dir) and any(
            name.startswith(f"{iso_name}.B9999")
            for name in os.listdir(run_dir)
    ):
        print(f"{run_name} already modelled")
        continue

    os.makedirs(run_dir, exist_ok=True)

    # Reads the full isoform sequence
    iso_path = os.path.join(
        output_longest_translated_fasta,
        f"{iso_name}.fasta",
    )
    # Reads the sequence from FASTA and converts it to a string
    iso_seq = str(
        SeqIO.read(iso_path, "fasta").seq
    )

    # Selects the experimental or AlphaFold template file
    if len(tpl_id) == 4:
        # Experimental PDB files are all 4 letters
        tpl_path = os.path.join(
            pdb_dir, f"{tpl_id}.pdb"
        )
    else:
        # For AlphaFold models
        tpl_path = os.path.join(
            af_models_dir, f"{tpl_id}.pdb"
        )

    # Reads the template structure [13]
    # Strucuter -> Model -> Chain -> Residue -> Atom
    parser = PDBParser(QUIET=True)
    tpl_structure = parser.get_structure(
        tpl_id, tpl_path
    )
    # Selects the first structural model
    tpl_model = tpl_structure[0]

    # Places the isoform chain first and partner chain second
    prepared_model = Model(0) # Bio.PDB.Model.Model
    # Will store coordinate derived sequences
    tpl_sequences = []

    # First isoform-side chain then partner chain
    for chain_id in [iso_chain, inter_chain]:
        chain = tpl_model[chain_id].copy()
        sequence = ""

        # Uses only amino acids that have coordinates in the template
        for res in list(chain):
            # list(chain) creates a fixed snapshot of the residues
            # so that it protects against removal later on
            if is_aa(res, standard=False):
                # Accepts standard and recognized modified AAs
                # Waters, ions, ligands and other non-amino-acid
                # residues are rejected
                sequence += seq1(
                    # Converts 3 letter residue name to one-letter
                    res.resname,
                    custom_map=aa_map,
                )
            else:
                # Removes non-AA residues
                chain.detach_child(res.id)

        # Adds cleaned chain to the new model and stores its sequence
        # After two iteratins:
        # tpl_sequences[0] would be the isoform-side template sequence
        # tpl_sequences[1] would be the partner template sequence
        prepared_model.add(chain)
        tpl_sequences.append(sequence)

    # Labels
    tpl_label = (
        f"{tpl_id}_{iso_chain}_{inter_chain}"
    )
    prepared_label = f"TPL_{tpl_label}"

    # To read from beginning of the first chain to the end of the
    # second chain:
    iso_first = (

        prepared_model[iso_chain].get_list()[0].id[1]
    )
    # get_list()[0] is the first amino-acid residue in the isoform-side chain
    # Biopython residue ID structure: (hetero flag, residue number, insertion code)
    # id[1] represents residue number
    inter_last = (
        prepared_model[inter_chain].get_list()[-1].id[1]
    )

    prepared_path = os.path.join(
        run_dir, f"{prepared_label}.pdb"
    )

    # Saves the prepared template
    writer = PDBIO()
    writer.set_structure(prepared_model)
    writer.save(prepared_path)

    # Aligns the full isoform to the residues present in the template [10]
    aligner = PairwiseAligner()
    aligner.mode = "global"
    # BLOSUM62 is used
    aligner.substitution_matrix = substitution_matrices.load(
        "BLOSUM62"
    )
    # Opening a new gap costs -10
    aligner.open_gap_score = -10
    # Extending it costs -0.5
    aligner.extend_gap_score = -0.5

    # Selects the best (first) alignment
    alignment = aligner.align(
        tpl_sequences[0],
        iso_seq,
    )[0]

    # The two rows contain the gapped template and isoform sequences [10]
    # The partner chain stays the same as the template partner
    # alignment[0]  # gapped template-chain sequence
    # alignment[1]  # gapped full isoform sequence
    aligned_tpls = [
        alignment[0],
        tpl_sequences[1],
    ]
    aligned_targets = [
        alignment[1],
        tpl_sequences[1],
    ]

    # Writes the two-chain MODELLER alignment (.ali) [11]
    alignment_name = (
        f"ALN_COMPLEX__{iso_name}__{tpl_id}_"
        f"{iso_chain}__partner_{inter_chain}.ali"
    )
    alignment_path = os.path.join(run_dir, alignment_name)

    with open(alignment_path, "w") as file:
        file.write(f">P1;{tpl_label}\n")
        file.write(
            f"structureX:{prepared_label}:"
            f"{iso_first}:{iso_chain}:"
            f"{inter_last}:{inter_chain}::::\n"
        )
        file.write("/".join(aligned_tpls) + "*\n\n")
        file.write(f">P1;{iso_name}\n")
        file.write(f"sequence:{iso_name}::::::::\n")
        file.write("/".join(aligned_targets) + "*\n")

    # MODELLER writes its files into the current directory
    previous_dir = os.getcwd()
    os.chdir(run_dir)

    # Returns to the previous directory even if MODELLER raises an error
    try:
        # Minimal output from the MODELLER
        log.minimal()
        # Initializes MODELLER's libraries and settigs
        environment = Environ()
        # Location for prepared template PDB files
        environment.io.atom_files_directory = [run_dir]

        # Uses six MODELLER workers to build models in parallel [31]
        parallel_job = job()
        for _ in range(6):
            parallel_job.append(local_slave())

        model = AutoModel(
            environment,
            # PIR alignment
            alnfile=alignment_path,
            # Code after >P1; for the known template
            knowns=tpl_label,
            # Code after >P1; for the target
            sequence=iso_name,
            # Gets DOPE and Normalized DOPE scores
            assess_methods=(
                assess.DOPE,
                assess.normalized_dope,
            ),
        )

        model.use_parallel_job(parallel_job)

        # Uses the highest optimization settings with 100 models for each pair [27][28][29][30]
        model.starting_model = 1
        model.ending_model = 100
        model.library_schedule = autosched.slow
        model.max_var_iterations = 500
        model.md_level = refine.very_slow
        model.repeat_optimization = 3

        # Stops models whose objective function remains above this limit [24]
        model.max_molpdf = 1e6

        model.make()

        print(f"Finished: {run_name}")

    finally:
        os.chdir(previous_dir)


# =========== MODEL QC AND SELECTION ============

# Will store the best model from each modelling job
best_model_rows = []

# Creates a result tsv file for each modelling job
for _, model_row in model_jobs.iterrows():
    iso_name = model_row["Iso_ID"]
    tpl_id = model_row["Base_Tpl"]
    iso_chain = model_row["Iso_Chain"]
    inter_chain = model_row["Inter_Chain"]

    run_name = (
        f"{iso_name}__{tpl_id}_{iso_chain}"
        f"__partner_{inter_chain}"
    )
    run_dir = os.path.abspath(
        os.path.join(modeller_output_dir, run_name)
    )

    # Location of the score table for this modelling job
    results_path = os.path.join(
        run_dir, "model_results.tsv"
    )

    # Uses the scores calculated during the previous MODELLER run
    if os.path.exists(results_path):
        results = pd.read_csv(results_path, sep="\t")

    # Calculates the scores when the result table is not present
    else:
        # Selects the models generated by MODELLER
        model_files = sorted([
            name for name in os.listdir(run_dir)
            if name.startswith(f"{iso_name}.B9999")
        ])

        # Will store the record for each candidate model
        model_rows = []

        # Reads the MODELLER scores stored in each model PDB file
        for model_name in model_files:
            model_path = os.path.join(run_dir, model_name)

            with open(model_path, "r") as file:
                for line in file:
                    # Looks for MODELLER entries in the model files
                    if line.startswith(
                        "REMARK   6 MODELLER OBJECTIVE FUNCTION:"
                    ):
                        molpdf = float(line.split(":")[-1])

                    elif line.startswith("REMARK   6 DOPE score:"):
                        dope = float(line.split(":")[-1])

                    elif line.startswith(
                        "REMARK   6 Normalized DOPE score:"
                    ):
                        norm_dope = float(line.split(":")[-1])

            model_rows.append({
                "name": model_name,
                "molpdf": molpdf,
                "DOPE score": dope,
                "Normalized DOPE score": norm_dope,
            })

        results = pd.DataFrame(model_rows)

        # Sorts the models using the same selection order
        # Normalized DOPE score gets the highest priority
        # followed by DOPE score, molpdf and name [25][26]
        results = results.sort_values([
            "Normalized DOPE score",
            "DOPE score",
            "molpdf",
            "name",
        ]).reset_index(drop=True)

        # Writes all model scores for this modelling job
        results.to_csv(
            results_path,
            sep="\t",
            index=False,
        )

    # Selects the first model in the ranked table
    best_model = results.iloc[0]

    # Stores the selected model for the combined result table
    best_model_rows.append({
        "Model_Job": run_name,
        "Iso_ID": iso_name,
        "Base_Tpl": tpl_id,
        "Iso_Chain": iso_chain,
        "Inter_Chain": inter_chain,
        # MODELLER assigns output chains alphabetically [23]
        "Model_Iso_Chain": "A",
        "Model_Inter_Chain": "B",
        "Model_Name": best_model["name"],
        "molpdf": best_model["molpdf"],
        "DOPE_Score": best_model["DOPE score"],
        "Normalized_DOPE_Score": (
            best_model["Normalized DOPE score"]
        ),
    })

    # Writes the selected model inside its modelling directory as .txt file
    best_path = os.path.join(
        run_dir, "best_model.txt"
    )

    with open(best_path, "w") as file:
        file.write(
            "Best model by normalized DOPE: "
            f"{best_model['name']}\n"
        )
        file.write(
            "Normalized DOPE score: "
            f"{best_model['Normalized DOPE score']}\n"
        )
        file.write(
            f"DOPE score: {best_model['DOPE score']}\n"
        )
        file.write(
            f"molpdf: {best_model['molpdf']}\n"
        )

    print(
        f"Selected model for {run_name}: "
        f"{best_model['name']}"
    )

# Creates one row for each modelling job
best_model_job_df = pd.DataFrame(best_model_rows)

# Adds the candidate and interactor information to the models [18]
best_model_df = candidate_df.merge(
    best_model_job_df,
    on=[
        "Iso_ID",
        "Base_Tpl",
        "Iso_Chain",
        "Inter_Chain",
    ],
    how="left", # preserve every row from candidate_df
)

# Separates the original template chains from the model output chains
best_model_df = best_model_df.rename(columns={
    "Iso_Chain": "Tpl_Iso_Chain",
    "Inter_Chain": "Tpl_Inter_Chain",
})

best_model_df = best_model_df[[
    "Iso_ID",
    "Inter_ID",
    "Interaction_Found",
    "Source",
    "Model_Job",
    "Base_Tpl",
    "Tpl_Iso_Chain",
    "Model_Iso_Chain",
    "Tpl_Inter_Chain",
    "Model_Inter_Chain",
    "Model_Name",
    "molpdf",
    "DOPE_Score",
    "Normalized_DOPE_Score",
]]

# Writes the selected model for each candidate
best_model_df.to_csv(
    os.path.join(
        modeller_output_dir,
        "best_modeller_models.tsv",
    ),
    sep="\t",
    index=False,
)

print("Best MODELLER models:")
print(best_model_df)


# ================ INTERFACE ANALYSIS ===================

# Output folder for interface results
interface_output_dir = "interface_analysis"
os.makedirs(interface_output_dir, exist_ok=True)

# Maximum heavy atom distance used to define a contact
contact_cutoff = 6.0

# Reads the selected MODELLER structures [13]
parser = PDBParser(QUIET=True)

# Aligns corresponding residues between two model chains [10]
interface_aligner = PairwiseAligner()
interface_aligner.mode = "global"
interface_aligner.substitution_matrix = substitution_matrices.load(
    "BLOSUM62"
)
interface_aligner.open_gap_score = -10
interface_aligner.extend_gap_score = -0.5


# The following functions are used for each model comparison
# Keeping them as functions prevents repeating the same alignment
# and contact calculations for model 1 and model 2


def residue_id(res):
    '''Creates a readable ID for one residue.'''
    # Bio.PDB stores an optional insertion code with the residue number
    insertion = res.id[2].strip()
    # Example output: A:25:GLY
    return (
        f"{res.get_parent().id}:{res.id[1]}"
        f"{insertion}:{res.resname}"
    )


def get_chain_residues(model, chain_id):
    '''Gets the amino acids that have CA coordinates from one chain.'''
    residues = []

    # Selects the requested chain from the MODELLER structure
    for res in model[chain_id]:
        # Non-standard amino acids are allowed but residues without
        # a CA atom cannot be used for structural alignment
        if is_aa(res, standard=False) and "CA" in res:
            residues.append(res)

    return residues


def align_residues(residues_1, residues_2, iso_2):
    '''Aligns two model chains and matches their residue numbers.'''
    # Converts three-letter amino acid names into sequences
    sequence_1 = ""
    for res in residues_1:
        sequence_1 += seq1(res.resname)

    sequence_2 = ""
    for res in residues_2:
        sequence_2 += seq1(res.resname)

    # Selects the highest-scoring pairwise alignment [10]
    alignment = interface_aligner.align(
        sequence_1, sequence_2
    )[0]

    # The alignment rows include gaps shown with a dash
    aligned_1 = alignment[0]
    aligned_2 = alignment[1]

    # Keeps track of positions in the original ungapped sequences
    index_1 = 0
    index_2 = 0
    aligned_pairs = []

    for aa_1, aa_2 in zip(aligned_1, aligned_2):
        # Stores positions only when both sequences have a residue
        if aa_1 != "-" and aa_2 != "-":
            aligned_pairs.append((index_1, index_2))

        if aa_1 != "-":
            index_1 += 1
        if aa_2 != "-":
            index_2 += 1

    # Model 1 is used as the common numbering system
    map_1 = {}
    for res in residues_1:
        map_1[residue_id(res)] = residue_id(res)

    # Unaligned model 2 residues first receive a unique label
    # This prevents them from being mistaken for model 1 residues
    # that happen to have the same residue number
    map_2 = {}
    for res in residues_2:
        map_2[residue_id(res)] = (
            f"{iso_2}_{residue_id(res)}"
        )

    # Aligned model 2 residues receive the equivalent model 1 ID
    for index_1, index_2 in aligned_pairs:
        res_1 = residues_1[index_1]
        res_2 = residues_2[index_2]
        map_2[residue_id(res_2)] = residue_id(res_1)

    # Counts identical amino acids among the aligned positions
    identical = 0
    for index_1, index_2 in aligned_pairs:
        if sequence_1[index_1] == sequence_2[index_2]:
            identical += 1

    # These values describe the reliability of the residue mapping
    identity = identical / len(aligned_pairs)
    coverage_1 = len(aligned_pairs) / len(residues_1)
    coverage_2 = len(aligned_pairs) / len(residues_2)

    return (
        aligned_pairs, map_1, map_2,
        identity, coverage_1, coverage_2,
    )


def get_contacts(iso_residues, inter_residues):
    '''Finds residue pairs within the interface cutoff [13].'''
    iso_atoms = []
    # Collects all non-hydrogen atoms from the isoform chain
    for res in iso_residues:
        for atom in res:
            if atom.element != "H":
                iso_atoms.append(atom)

    inter_atoms = []
    # Collects all non-hydrogen atoms from the interactor chain
    for res in inter_residues:
        for atom in res:
            if atom.element != "H":
                inter_atoms.append(atom)

    iso_chain = iso_residues[0].get_parent().id
    inter_chain = inter_residues[0].get_parent().id

    contacts = {}

    # NeighborSearch finds atom pairs within 6 Angstroms
    atom_search = NeighborSearch(iso_atoms + inter_atoms)
    close_atoms = atom_search.search_all(
        contact_cutoff, level="A"
    )

    for atom_1, atom_2 in close_atoms:
        chain_1 = atom_1.get_parent().get_parent().id
        chain_2 = atom_2.get_parent().get_parent().id

        # Places the isoform atom first if the pair connects
        # the isoform and interactor chains
        if chain_1 == iso_chain and chain_2 == inter_chain:
            iso_atom = atom_1
            inter_atom = atom_2
        elif chain_2 == iso_chain and chain_1 == inter_chain:
            iso_atom = atom_2
            inter_atom = atom_1
        else:
            # Ignores atom pairs that are within the same chain
            continue

        iso_res = residue_id(iso_atom.get_parent())
        inter_res = residue_id(inter_atom.get_parent())
        pair = (iso_res, inter_res)

        # Bio.PDB calculates the distance between the two atoms
        distance = iso_atom - inter_atom

        # Several atoms may connect the same two residues
        # Only their shortest heavy atom distance is kept
        if pair not in contacts:
            contacts[pair] = distance
        elif distance < contacts[pair]:
            contacts[pair] = distance

    return contacts


def map_contact_pairs(contacts, iso_map, inter_map):
    '''Converts contact pairs into the model 1 numbering system.'''
    mapped_contacts = {}

    for pair, distance in contacts.items():
        iso_res = pair[0]
        inter_res = pair[1]

        mapped_iso_res = iso_map[iso_res]
        mapped_inter_res = inter_map[inter_res]
        mapped_pair = (mapped_iso_res, mapped_inter_res)

        # Keeps the shortest distance if multiple contacts map
        # to the same residue pair
        if mapped_pair not in mapped_contacts:
            mapped_contacts[mapped_pair] = distance
        elif distance < mapped_contacts[mapped_pair]:
            mapped_contacts[mapped_pair] = distance

    return mapped_contacts


def jaccard(set_1, set_2):
    '''Calculates the Jaccard similarity between two sets [33].'''
    shared = set_1 & set_2
    combined = set_1 | set_2
    return len(shared) / len(combined)


# Adds the protein name shared by the isoforms
best_model_df["Base_Iso"] = (
    best_model_df["Iso_ID"].str.rsplit("_", n=1).str[0]
)

# Columns that define a comparable group of models
interface_group_cols = [
    "Base_Iso",
    "Inter_ID",
    "Base_Tpl",
    "Tpl_Iso_Chain",
    "Tpl_Inter_Chain",
]

# Stores the identity of each model comparison
comparison_rows = []

# Stores the alignment and structural similarity measurements
alignment_rows = []

# Stores the interface residue measurements
residue_rows = []

# Stores the mapped contact summary measurements
contact_summary_rows = []

# Stores every individual mapped contact
contact_rows = []

for _, interface_group in best_model_df.groupby(
        interface_group_cols, sort=False
):
    # Starts the row numbers again within each group
    interface_group = interface_group.reset_index(drop=True)

    # Compares every unique pair within the same template group
    for i in range(len(interface_group)):
        for j in range(i + 1, len(interface_group)):
            model_row_1 = interface_group.iloc[i]
            model_row_2 = interface_group.iloc[j]

            # Keeps only isoforms with opposite interaction results
            if (
                model_row_1["Interaction_Found"] ==
                model_row_2["Interaction_Found"]
            ):
                continue

            # Model 1 is negative and model 2 is positive
            if model_row_1["Interaction_Found"] == "positive":
                model_row_1, model_row_2 = (
                    model_row_2, model_row_1
                )

            # Reads the names used throughout this comparison
            iso_1 = model_row_1["Iso_ID"]
            iso_2 = model_row_2["Iso_ID"]
            inter_name = model_row_1["Inter_ID"]
            tpl_id = model_row_1["Base_Tpl"]
            tpl_iso_chain = model_row_1[
                "Tpl_Iso_Chain"
            ]
            tpl_inter_chain = model_row_1[
                "Tpl_Inter_Chain"
            ]

            # Makes a short label for figures and tables
            iso_comparison = f"{iso_1}__vs__{iso_2}"

            # Makes a unique label that also includes the template
            comparison_id = (
                f"{iso_1}__vs__{iso_2}__{inter_name}__"
                f"{tpl_id}__{tpl_iso_chain}__"
                f"{tpl_inter_chain}"
            )

            # Locates the two selected MODELLER structures
            model_path_1 = os.path.join(
                modeller_output_dir,
                model_row_1["Model_Job"],
                model_row_1["Model_Name"],
            )
            model_path_2 = os.path.join(
                modeller_output_dir,
                model_row_2["Model_Job"],
                model_row_2["Model_Name"],
            )

            # Reads the first structural model from each PDB
            model_1 = parser.get_structure(
                iso_1, model_path_1
            )[0]
            model_2 = parser.get_structure(
                iso_2, model_path_2
            )[0]

            # Gets the modelled isoform and interactor residues
            iso_residues_1 = get_chain_residues(
                model_1, model_row_1["Model_Iso_Chain"]
            )
            iso_residues_2 = get_chain_residues(
                model_2, model_row_2["Model_Iso_Chain"]
            )
            inter_residues_1 = get_chain_residues(
                model_1, model_row_1["Model_Inter_Chain"]
            )
            inter_residues_2 = get_chain_residues(
                model_2, model_row_2["Model_Inter_Chain"]
            )

            # Aligns the common interactor chains
            (
                inter_pairs, inter_map_1, inter_map_2,
                _, _, _,
            ) = align_residues(
                inter_residues_1, inter_residues_2, iso_2
            )

            # Aligns the two isoform chains for contact mapping
            (
                _, iso_map_1, iso_map_2,
                iso_identity, iso_coverage_1, iso_coverage_2,
            ) = align_residues(
                iso_residues_1, iso_residues_2, iso_2
            )

            # Reads the complete processed interactor sequence [5]
            inter_path = os.path.join(
                inter_filtered_dir,
                f"{inter_name}.fasta",
            )
            canonical_inter_seq = str(
                SeqIO.read(inter_path, "fasta").seq
            ).rstrip("*")

            # Converts the model 1 interactor chain into a sequence
            model_inter_seq = ""
            for res in inter_residues_1:
                model_inter_seq += seq1(res.resname)

            # Aligns the coordinate-derived partner chain to the
            # complete processed interactor sequence [10]
            canonical_alignment = interface_aligner.align(
                model_inter_seq,
                canonical_inter_seq,
            )[0]

            # Will connect model 1 interactor residues to the
            # corresponding positions in the processed interactor FASTA
            canonical_inter_map = {}
            model_index = 0
            canonical_index = 0

            for model_aa, canonical_aa in zip(
                    canonical_alignment[0],
                    canonical_alignment[1],
            ):
                # Stores a mapping only when both sequences have a residue
                if model_aa != "-" and canonical_aa != "-":
                    model_res = inter_residues_1[model_index]
                    canonical_inter_map[residue_id(model_res)] = (
                        f"{inter_name}:{canonical_index + 1}:"
                        f"{seq3(canonical_aa).upper()}"
                    )

                # Moves through the two original ungapped sequences
                if model_aa != "-":
                    model_index += 1
                if canonical_aa != "-":
                    canonical_index += 1

            # Superposes the interactor CA atoms and calculates RMSD [32]
            super_imposer = Superimposer()
            super_imposer.set_atoms(
                [
                    inter_residues_1[x]["CA"]
                    for x, _ in inter_pairs
                ],
                [
                    inter_residues_2[y]["CA"]
                    for _, y in inter_pairs
                ],
            )

            # Finds the interface contacts in each model
            contacts_1 = get_contacts(
                iso_residues_1, inter_residues_1
            )
            contacts_2 = get_contacts(
                iso_residues_2, inter_residues_2
            )

            # Maps both sides of every contact to a shared numbering frame
            mapped_contacts_1 = map_contact_pairs(
                contacts_1, iso_map_1, inter_map_1
            )
            mapped_contacts_2 = map_contact_pairs(
                contacts_2, iso_map_2, inter_map_2
            )

            # Reverses the maps so each common residue can be traced
            # back to its original residue in both selected models
            iso_reverse_1 = {
                mapped: original
                for original, mapped in iso_map_1.items()
            }
            iso_reverse_2 = {
                mapped: original
                for original, mapped in iso_map_2.items()
            }
            inter_reverse_1 = {
                mapped: original
                for original, mapped in inter_map_1.items()
            }
            inter_reverse_2 = {
                mapped: original
                for original, mapped in inter_map_2.items()
            }

            # Converts the contact dictionaries into sets
            contact_set_1 = set(mapped_contacts_1)
            contact_set_2 = set(mapped_contacts_2)

            # Interactor residues contacted by each isoform
            inter_set_1 = {
                inter_map_1[inter_res]
                for _, inter_res in contacts_1
            }
            inter_set_2 = {
                inter_map_2[inter_res]
                for _, inter_res in contacts_2
            }

            # Inter residues contacted in both models
            shared_inter = inter_set_1 & inter_set_2

            # Interactor contacts lost in the positive isoform
            lost_inter = inter_set_1 - inter_set_2

            # Interactor contacts gained in the positive isoform
            gained_inter = inter_set_2 - inter_set_1

            # Mapped residue pairs found in both models
            shared_contacts = contact_set_1 & contact_set_2

            # Mapped contacts lost in the positive isoform
            lost_contacts = contact_set_1 - contact_set_2

            # Mapped contacts gained in the positive isoform
            gained_contacts = contact_set_2 - contact_set_1

            # Stores the model and template information
            comparison_rows.append({
                "Comparison_ID": comparison_id,
                "Iso_Comparison": iso_comparison,
                "Base_Iso": model_row_1["Base_Iso"],
                "Iso_1": iso_1,
                "Iso_2": iso_2,
                "Inter_ID": inter_name,
                "Interaction_Found_1": model_row_1[
                    "Interaction_Found"
                ],
                "Interaction_Found_2": model_row_2[
                    "Interaction_Found"
                ],
                "Source": model_row_1["Source"],
                "Base_Tpl": tpl_id,
                "Tpl_Iso_Chain": tpl_iso_chain,
                "Tpl_Inter_Chain": tpl_inter_chain,
                "Model_Name_1": model_row_1["Model_Name"],
                "Model_Name_2": model_row_2["Model_Name"],
            })

            # Stores the alignment quality measurements
            alignment_rows.append({
                "Comparison_ID": comparison_id,
                "Inter_CA_Pairs": len(inter_pairs),
                "Inter_CA_RMSD": round(
                    super_imposer.rms, 4
                ),
                "Iso_Sequence_Identity": round(
                    iso_identity, 4
                ),
                "Iso_Alignment_Coverage_1": round(
                    iso_coverage_1, 4
                ),
                "Iso_Alignment_Coverage_2": round(
                    iso_coverage_2, 4
                ),
            })

            # Stores the interface residue changes
            residue_rows.append({
                "Comparison_ID": comparison_id,
                "Iso_1_Interface_Residue_Count": len({
                    iso_res for iso_res, _ in contacts_1
                }),
                "Iso_2_Interface_Residue_Count": len({
                    iso_res for iso_res, _ in contacts_2
                }),
                "Shared_Inter_Interface_Residue_Count": len(
                    shared_inter
                ),
                "Lost_Inter_Interface_Residue_Count": len(
                    lost_inter
                ),
                "Gained_Inter_Interface_Residue_Count": len(
                    gained_inter
                ),
                "Inter_Interface_Jaccard": round(
                    jaccard(inter_set_1, inter_set_2), 4
                ),
            })

            # Stores the mapped contact changes
            contact_summary_rows.append({
                "Comparison_ID": comparison_id,
                "Shared_Mapped_Contact_Count": len(
                    shared_contacts
                ),
                "Lost_Mapped_Contact_Count": len(
                    lost_contacts
                ),
                "Gained_Mapped_Contact_Count": len(
                    gained_contacts
                ),
                "Mapped_Contact_Jaccard": round(
                    jaccard(contact_set_1, contact_set_2), 4
                ),
            })

            # Stores the mapped contacts for later plots and tables
            for change_type, contact_set in [
                ("shared", shared_contacts),
                ("lost_in_positive", lost_contacts),
                ("gained_in_positive", gained_contacts),
            ]:
                for iso_res, inter_res in sorted(contact_set):
                    pair = (iso_res, inter_res)
                    contact_rows.append({
                        "Comparison_ID": comparison_id,
                        "Change_Type": change_type,
                        "Iso_Residue": iso_res,
                        "Inter_Residue": inter_res,
                        # Original residues in the two selected models
                        "Iso_1_Residue": iso_reverse_1.get(iso_res),
                        "Iso_2_Residue": iso_reverse_2.get(iso_res),
                        "Inter_1_Residue": inter_reverse_1.get(inter_res),
                        "Inter_2_Residue": inter_reverse_2.get(inter_res),
                        # Position in the complete processed interactor FASTA
                        "Canonical_Inter_Residue": (
                            canonical_inter_map.get(inter_res)
                        ),
                        "Distance_Iso_1": (
                            round(mapped_contacts_1[pair], 3)
                            if pair in mapped_contacts_1 else None
                        ),
                        "Distance_Iso_2": (
                            round(mapped_contacts_2[pair], 3)
                            if pair in mapped_contacts_2 else None
                        ),
                    })


# Makes one DataFrame for each analysis stage
comparison_df = pd.DataFrame(comparison_rows)
alignment_df = pd.DataFrame(alignment_rows)
residue_df = pd.DataFrame(residue_rows)
contact_summary_df = pd.DataFrame(contact_summary_rows)
contact_df = pd.DataFrame(contact_rows)

# Combines the comparison-level results using their unique ID [18]
interface_df = comparison_df.merge(
    alignment_df, on="Comparison_ID"
)
interface_df = interface_df.merge(
    residue_df, on="Comparison_ID"
)
interface_df = interface_df.merge(
    contact_summary_df, on="Comparison_ID"
)

# Writes the mapped residue contacts
contact_df.to_csv(
    os.path.join(
        interface_output_dir,
        "interface_contact_changes.tsv",
    ),
    sep="\t",
    index=False,
)

print("Interface comparison results:")
print(interface_df)
print(f"Mapped contact rows: {len(contact_df)}")

# =============== FOLDX ==================

# Folder containing the FoldX files for every selected model
foldx_output_dir = "foldx_results"
os.makedirs(foldx_output_dir, exist_ok=True)

# FoldX is not installed through the conda environment
# The separately downloaded program must be available as "foldx"
foldx_path = "foldx"

# Some MODELLER structures are used for more than one interactor row
# Each structure only needs to be analysed by FoldX once
foldx_models_df = best_model_df.drop_duplicates(
    ["Model_Job", "Model_Name"]
).copy()

foldx_rows = []

for _, foldx_model in foldx_models_df.iterrows():
    # Reads the model information
    model_job = foldx_model["Model_Job"]
    iso_name = foldx_model["Iso_ID"]
    model_name = foldx_model["Model_Name"]

    # Makes one FoldX folder for each model
    foldx_run_dir = os.path.join(
        foldx_output_dir, model_job
    )
    os.makedirs(foldx_run_dir, exist_ok=True)

    # Locates the selected MODELLER structure
    model_path = os.path.join(
        modeller_output_dir,
        model_job,
        model_name,
    )

    # Uses the same simple filenames for every FoldX run
    foldx_input = os.path.join(
        foldx_run_dir, "complex.pdb"
    )
    repaired_model = os.path.join(
        foldx_run_dir, "complex_Repair.pdb"
    )
    interaction_file = os.path.join(
        foldx_run_dir,
        "Interaction_complex_Repair_AC.fxout",
    )

    # Separates the isoform and interactor sides of the complex
    chain_group = (
        f"{foldx_model['Model_Iso_Chain']},"
        f"{foldx_model['Model_Inter_Chain']}"
    )

    # Skips FoldX when its final interaction output is already present
    if not os.path.exists(interaction_file):
        # Copies the selected MODELLER structure only once
        if not os.path.exists(foldx_input):
            shutil.copyfile(model_path, foldx_input)

        # Repairs the structure before the energy analysis [34]
        if not os.path.exists(repaired_model):
            with open(
                os.path.join(foldx_run_dir, "repair_stdout.txt"),
                "w",
            ) as output:
                with open(
                    os.path.join(foldx_run_dir, "repair_stderr.txt"),
                    "w",
                ) as error:
                    subprocess.run(
                        [
                            foldx_path,
                            "--command=RepairPDB",
                            "--pdb=complex.pdb",
                        ],
                        cwd=foldx_run_dir,
                        stdout=output,
                        stderr=error,
                        check=True,
                    )

        # Calculates the interaction energy between the chains [35]
        with open(
            os.path.join(foldx_run_dir, "analyse_stdout.txt"),
            "w",
        ) as output:
            with open(
                os.path.join(foldx_run_dir, "analyse_stderr.txt"),
                "w",
            ) as error:
                subprocess.run(
                    [
                        foldx_path,
                        "--command=AnalyseComplex",
                        "--pdb=complex_Repair.pdb",
                        f"--analyseComplexChains={chain_group}",
                    ],
                    cwd=foldx_run_dir,
                    stdout=output,
                    stderr=error,
                    check=True,
                )

    # Reads the interaction energy from the FoldX table [35]
    with open(interaction_file, "r") as result:
        for line in result:
            if line.startswith("Pdb\t"):
                foldx_header = line.strip().split("\t")
                foldx_values = next(result).strip().split("\t")
                break

    energy_column = foldx_header.index("Interaction Energy")
    interaction_energy = float(foldx_values[energy_column])

    # Stores one result for each unique model
    foldx_rows.append({
        "Model_Job": model_job,
        "Iso_ID": iso_name,
        "Base_Tpl": foldx_model["Base_Tpl"],
        "Tpl_Iso_Chain": foldx_model[
            "Tpl_Iso_Chain"
        ],
        "Tpl_Inter_Chain": foldx_model[
            "Tpl_Inter_Chain"
        ],
        "Model_Iso_Chain": foldx_model[
            "Model_Iso_Chain"
        ],
        "Model_Inter_Chain": foldx_model[
            "Model_Inter_Chain"
        ],
        "Model_Name": model_name,
        "FoldX_Interaction_Energy": interaction_energy,
    })


# Makes the model-level FoldX result table
foldx_df = pd.DataFrame(foldx_rows)

# Writes one FoldX energy for each unique model
foldx_df.to_csv(
    os.path.join(
        foldx_output_dir,
        "foldx_model_results.tsv",
    ),
    sep="\t",
    index=False,
)

print("FoldX model results:")
print(foldx_df)

# Connects each model name with its FoldX energy
foldx_energy = foldx_df.set_index(
    "Model_Name"
)["FoldX_Interaction_Energy"]

# Keeps only the information needed for the FoldX comparison
foldx_summary_df = interface_df[[
    "Comparison_ID",
    "Iso_Comparison",
    "Base_Iso",
    "Inter_ID",
    "Source",
    "Base_Tpl",
    "Tpl_Iso_Chain",
    "Tpl_Inter_Chain",
    "Iso_1",
    "Interaction_Found_1",
    "Model_Name_1",
    "Iso_2",
    "Interaction_Found_2",
    "Model_Name_2",
]].copy()

# Uses Y2H names in the final table
foldx_summary_df = foldx_summary_df.rename(columns={
    "Interaction_Found_1": "Y2H_Result_1",
    "Interaction_Found_2": "Y2H_Result_2",
})

# Adds the FoldX energy for the negative isoform model
foldx_summary_df["FoldX_Interaction_Energy_1"] = (
    foldx_summary_df["Model_Name_1"].map(foldx_energy)
)

# Adds the FoldX energy for the positive isoform model
foldx_summary_df["FoldX_Interaction_Energy_2"] = (
    foldx_summary_df["Model_Name_2"].map(foldx_energy)
)

# Model 1 is Y2H negative and model 2 is Y2H positive
# A negative difference means more favorable binding for model 2
foldx_summary_df[
    "FoldX_Energy_Difference_Positive_Minus_Negative"
] = (
    foldx_summary_df["FoldX_Interaction_Energy_2"] -
    foldx_summary_df["FoldX_Interaction_Energy_1"]
)

# Describes which isoform has the more favorable FoldX energy
foldx_summary_df["FoldX_Energy_Direction"] = (
    "same_energy"
)
foldx_summary_df.loc[
    foldx_summary_df[
        "FoldX_Energy_Difference_Positive_Minus_Negative"
    ] > 0,
    "FoldX_Energy_Direction",
] = (
    "negative_iso_more_favorable"
)
foldx_summary_df.loc[
    foldx_summary_df[
        "FoldX_Energy_Difference_Positive_Minus_Negative"
    ] < 0,
    "FoldX_Energy_Direction",
] = "positive_iso_more_favorable"

# Marks whether the FoldX direction agrees with the Y2H result
foldx_summary_df["FoldX_Direction_Matches_Y2H"] = "no"
foldx_summary_df.loc[
    foldx_summary_df[
        "FoldX_Energy_Difference_Positive_Minus_Negative"
    ] < 0,
    "FoldX_Direction_Matches_Y2H",
] = "yes"

# ================== SUMMARIES ======================

# Output folder for the complete result table
summary_output_dir = "summary_results"
os.makedirs(summary_output_dir, exist_ok=True)

# Adds the FoldX results to the interface results [18]
all_results_df = interface_df.merge(
    foldx_summary_df[[
        "Comparison_ID",
        "FoldX_Interaction_Energy_1",
        "FoldX_Interaction_Energy_2",
        "FoldX_Energy_Difference_Positive_Minus_Negative",
        "FoldX_Energy_Direction",
        "FoldX_Direction_Matches_Y2H",
    ]],
    on="Comparison_ID",
)

# Writes every comparison into one complete table
all_results_df.to_csv(
    os.path.join(
        summary_output_dir,
        "all_results.tsv",
    ),
    sep="\t",
    index=False,
)

# Prints candidate and modelling counts for the thesis text
print("Template selection summary:")

for source_name in ["StructMAn", "MMseqs2"]:
    source_df = candidate_df[
        candidate_df["Source"] == source_name
    ]
    source_jobs = source_df.drop_duplicates([
        "Iso_ID",
        "Base_Tpl",
        "Iso_Chain",
        "Inter_Chain",
    ])
    pdb_count = (
        source_df["Base_Tpl"].str.len() == 4
    ).sum()
    predicted_count = len(source_df) - pdb_count
    base_count = source_df["Base_Iso"].nunique()
    iso_count = source_df["Iso_ID"].nunique()
    inter_count = source_df["Inter_ID"].nunique()
    tpl_count = source_df["Base_Tpl"].nunique()

    print(
        f"{source_name}: {len(source_df)} candidate rows, "
        f"{len(source_jobs)} modelling jobs"
    )
    print(
        f"{base_count} base isoforms, {iso_count} isoforms, "
        f"{inter_count} interactors, {tpl_count} templates"
    )
    print(
        f"{pdb_count} PDB rows and "
        f"{predicted_count} predicted rows"
    )

all_jobs = candidate_df.drop_duplicates([
    "Iso_ID",
    "Base_Tpl",
    "Iso_Chain",
    "Inter_Chain",
])
all_pdb_count = (
    candidate_df["Base_Tpl"].str.len() == 4
).sum()
all_base_count = candidate_df["Base_Iso"].nunique()
all_iso_count = candidate_df["Iso_ID"].nunique()
all_inter_count = candidate_df["Inter_ID"].nunique()
all_tpl_count = candidate_df["Base_Tpl"].nunique()

print(
    f"Combined: {len(candidate_df)} candidate rows, "
    f"{len(all_jobs)} modelling jobs"
)
print(
    f"{all_base_count} base isoforms, {all_iso_count} isoforms, "
    f"{all_inter_count} interactors, "
    f"{all_tpl_count} templates"
)
print(
    f"{all_pdb_count} PDB rows and "
    f"{len(candidate_df) - all_pdb_count} predicted rows"
)

# Keeps only comparisons below the final RMSD limit
reliable_results_df = all_results_df[
    all_results_df["Inter_CA_RMSD"] <= 2
].copy()

# Gives every selected model pair one name
reliable_results_df["Model_Pair"] = (
    reliable_results_df["Model_Name_1"] + "__" +
    reliable_results_df["Model_Name_2"]
)

# Groups templates and chain orientations from the same comparison [19]
foldx_result_df = reliable_results_df.groupby([
    "Base_Iso",
    "Inter_ID",
    "Iso_Comparison",
    "Iso_1",
    "Interaction_Found_1",
    "Iso_2",
    "Interaction_Found_2",
], as_index=False).agg(
    Model_Pair_Set=(
        "Model_Pair",
        lambda values: ", ".join(sorted(set(values))),
    ),
    Median_FoldX_Difference=(
        "FoldX_Energy_Difference_Positive_Minus_Negative",
        "median",
    ),
)

# Removes repeated labels that used the same structural models [20]
foldx_result_df = foldx_result_df.drop_duplicates([
    "Base_Iso",
    "Iso_Comparison",
    "Iso_1",
    "Interaction_Found_1",
    "Iso_2",
    "Interaction_Found_2",
    "Model_Pair_Set",
])

# Counts FoldX results that agree and disagree with Y2H
match_count = (
    foldx_result_df["Median_FoldX_Difference"] < 0
).sum()
opposite_count = (
    foldx_result_df["Median_FoldX_Difference"] > 0
).sum()
comparison_count = match_count + opposite_count

# Calculates the average and median energy differences
mean_foldx_difference = foldx_result_df[
    "Median_FoldX_Difference"
].mean()
median_foldx_difference = foldx_result_df[
    "Median_FoldX_Difference"
].median()

print("Complete result summary:")
print(f"All structural comparisons: {len(all_results_df)}")
print(
    "Comparisons with inter RMSD at or below 2 A: "
    f"{len(reliable_results_df)}"
)
print(
    "FoldX direction matching Y2H: "
    f"{match_count}/{comparison_count} "
    f"({100 * match_count / comparison_count:.2f}%)"
)
print(f"FoldX direction opposite to Y2H: {opposite_count}")
print(
    "Mean FoldX difference, positive minus negative (kcal/mol): "
    f"{mean_foldx_difference:.4f}"
)
print(
    "Median FoldX difference, positive minus negative (kcal/mol): "
    f"{median_foldx_difference:.4f}"
)


# ======== CHIMERAX FOLDER CREATION ==========

# Main folder for the ChimeraX images and sessions
chimerax_dir = os.path.join("figures", "chimerax")

# Creates the main folder if it is not already present
if not os.path.exists(chimerax_dir):
    os.makedirs(chimerax_dir)

# Stores the model information needed for ChimeraX
chimerax_rows = []

# Creates one folder for every comparison below the 2 A limit
for _, result_row in reliable_results_df.iterrows():
    source_name = result_row["Source"].lower()

    # Shortens the two Y2H results for the folder name
    interaction_1 = "pos" if result_row[
        "Interaction_Found_1"
    ] == "positive" else "neg"
    interaction_2 = "pos" if result_row[
        "Interaction_Found_2"
    ] == "positive" else "neg"

    # Combines the source, proteins, template and compared isoforms
    chimerax_name = (
        f"{source_name}__"
        f"{result_row['Base_Iso']}-{result_row['Inter_ID']}__"
        f"tmpl-{result_row['Base_Tpl']}-"
        f"{result_row['Tpl_Iso_Chain']}-"
        f"{result_row['Tpl_Inter_Chain']}__"
        f"{result_row['Iso_1']}-{interaction_1}_vs_"
        f"{result_row['Iso_2']}-{interaction_2}"
    )

    chimerax_result_dir = os.path.join(
        chimerax_dir,
        chimerax_name,
    )

    # Skips folders that were already created
    if not os.path.exists(chimerax_result_dir):
        os.makedirs(chimerax_result_dir)

    # Model 1 is Y2H negative and model 2 is Y2H positive
    negative_job = (
        f"{result_row['Iso_1']}__{result_row['Base_Tpl']}_"
        f"{result_row['Tpl_Iso_Chain']}__partner_"
        f"{result_row['Tpl_Inter_Chain']}"
    )
    positive_job = (
        f"{result_row['Iso_2']}__{result_row['Base_Tpl']}_"
        f"{result_row['Tpl_Iso_Chain']}__partner_"
        f"{result_row['Tpl_Inter_Chain']}"
    )

    # Keeps the exact model names and paths for later figure preparation
    negative_path = os.path.abspath(os.path.join(
        modeller_output_dir,
        negative_job,
        result_row["Model_Name_1"],
    ))
    positive_path = os.path.abspath(os.path.join(
        modeller_output_dir,
        positive_job,
        result_row["Model_Name_2"],
    ))

    chimerax_rows.append({
        "Case_ID": chimerax_name,
        "Comparison_ID": result_row["Comparison_ID"],
        "Source": result_row["Source"],
        "Base_Iso": result_row["Base_Iso"],
        "Inter_ID": result_row["Inter_ID"],
        "Base_Tpl": result_row["Base_Tpl"],
        "Tpl_Iso_Chain": result_row["Tpl_Iso_Chain"],
        "Tpl_Inter_Chain": result_row["Tpl_Inter_Chain"],
        "Negative_Iso": result_row["Iso_1"],
        "Negative_Model": result_row["Model_Name_1"],
        "Negative_PDB_Path": negative_path,
        "Positive_Iso": result_row["Iso_2"],
        "Positive_Model": result_row["Model_Name_2"],
        "Positive_PDB_Path": positive_path,
        "Inter_CA_RMSD": round(
            result_row["Inter_CA_RMSD"], 4
        ),
        "Inter_Interface_Jaccard": round(
            result_row["Inter_Interface_Jaccard"], 4
        ),
        "Mapped_Contact_Jaccard": round(
            result_row["Mapped_Contact_Jaccard"], 4
        ),
        "FoldX_Energy_Difference_Positive_Minus_Negative": (
            round(result_row[
                "FoldX_Energy_Difference_Positive_Minus_Negative"
            ], 4)
        ),
    })

# Writes the model and comparison information for the ChimeraX figures
chimerax_df = pd.DataFrame(chimerax_rows)
chimerax_df.to_csv(
    os.path.join(chimerax_dir, "chimerax_manifest.tsv"),
    sep="\t",
    index=False,
)

print(f"ChimeraX folders created: {len(chimerax_df)}")


# =============== SUPPLEMENTARY TABLES =================

# Output folder for the supplementary tables
supp_table_dir = "supplementary_tables"
os.makedirs(supp_table_dir, exist_ok=True)


# Keeps the candidate template information used before modelling
# The structure recommendation already contains the template and chains
# StructMAn and MMseqs2 scores are not mixed because they are different scores
candidate_table_df = candidate_df[[
    "Base_Iso",
    "Iso_ID",
    "Inter_ID",
    "Interaction_Found",
    "Source",
    "Structure_Recommendation",
]].copy()

# Uses full column names for the supplementary table
candidate_table_df = candidate_table_df.rename(columns={
    "Base_Iso": "Base_Isoform",
    "Iso_ID": "Isoform_ID",
    "Inter_ID": "Interactor_ID",
    "Interaction_Found": "Y2H_Result",
})

# Orders related isoforms and interactors together
candidate_table_df = candidate_table_df.sort_values([
    "Base_Isoform",
    "Isoform_ID",
    "Interactor_ID",
    "Source",
])

candidate_table_df.to_csv(
    os.path.join(
        supp_table_dir,
        "table_s1_candidate_templates.tsv",
    ),
    sep="\t",
    index=False,
)


# Starts the MODELLER table from the selected best models
model_table_df = best_model_df.copy()

# Combines the template name and original template chains
model_table_df["Template"] = (
    model_table_df["Base_Tpl"] + " (" +
    model_table_df["Tpl_Iso_Chain"] + ":" +
    model_table_df["Tpl_Inter_Chain"] + ")"
)

# Shows the chain names used in the produced model
model_table_df["Model_Chains"] = (
    model_table_df["Model_Iso_Chain"] + ":" +
    model_table_df["Model_Inter_Chain"]
)

# Keeps only the model identity and its selection scores
model_table_df = model_table_df[[
    "Iso_ID",
    "Inter_ID",
    "Interaction_Found",
    "Source",
    "Template",
    "Model_Chains",
    "Model_Name",
    "molpdf",
    "DOPE_Score",
    "Normalized_DOPE_Score",
]].copy()

# Uses readable names in the exported table
model_table_df = model_table_df.rename(columns={
    "Iso_ID": "Isoform_ID",
    "Inter_ID": "Interactor_ID",
    "Interaction_Found": "Y2H_Result",
    "Model_Name": "Selected_Model",
    "molpdf": "MOLPDF",
})

# Orders the models by isoform, interactor and template
model_table_df = model_table_df.sort_values([
    "Isoform_ID",
    "Interactor_ID",
    "Template",
])

model_table_df.to_csv(
    os.path.join(
        supp_table_dir,
        "table_s2_modeller_models.tsv",
    ),
    sep="\t",
    index=False,
)


# Starts the structural tables from the final comparisons
# reliable_results_df already contains only comparisons at or below 2 A
supp_result_df = reliable_results_df.copy()

# Model 1 is Y2H negative and model 2 is Y2H positive
supp_result_df["Isoform_Comparison"] = (
    supp_result_df["Iso_1"] + " (Y2H-) vs " +
    supp_result_df["Iso_2"] + " (Y2H+)"
)

# Combines the template name and directed chain pair
supp_result_df["Template"] = (
    supp_result_df["Base_Tpl"] + " (" +
    supp_result_df["Tpl_Iso_Chain"] + ":" +
    supp_result_df["Tpl_Inter_Chain"] + ")"
)

# Gives the shared columns their table names
supp_result_df = supp_result_df.rename(columns={
    "Inter_ID": "Interactor_ID",
    "Inter_CA_Pairs": "Interactor_CA_Pairs",
    "Inter_CA_RMSD": "Interactor_CA_RMSD_Angstrom",
})

# Orders the same comparisons equally in all three tables
supp_result_df = supp_result_df.sort_values([
    "Base_Iso",
    "Interactor_ID",
    "Iso_1",
    "Iso_2",
    "Base_Tpl",
    "Tpl_Iso_Chain",
    "Tpl_Inter_Chain",
])

# Columns shared by the interface, contact and FoldX tables
common_cols = [
    "Isoform_Comparison",
    "Interactor_ID",
    "Source",
    "Template",
    "Interactor_CA_Pairs",
    "Interactor_CA_RMSD_Angstrom",
]


# Keeps the interactor interface residue changes
interface_table_df = supp_result_df[
    common_cols + [
        "Shared_Inter_Interface_Residue_Count",
        "Lost_Inter_Interface_Residue_Count",
        "Gained_Inter_Interface_Residue_Count",
        "Inter_Interface_Jaccard",
    ]
].copy()

# Describes loss and gain relative to the Y2H-positive isoform
interface_table_df = interface_table_df.rename(columns={
    "Shared_Inter_Interface_Residue_Count": "Shared_Residues",
    "Lost_Inter_Interface_Residue_Count": "Lost_in_Y2H_Positive",
    "Gained_Inter_Interface_Residue_Count": "Gained_in_Y2H_Positive",
    "Inter_Interface_Jaccard": "Interface_Jaccard",
})

interface_table_df.to_csv(
    os.path.join(
        supp_table_dir,
        "table_s3_interface_residues.tsv",
    ),
    sep="\t",
    index=False,
)


# Keeps the mapped residue-pair contact changes
contact_table_df = supp_result_df[
    common_cols + [
        "Shared_Mapped_Contact_Count",
        "Lost_Mapped_Contact_Count",
        "Gained_Mapped_Contact_Count",
        "Mapped_Contact_Jaccard",
    ]
].copy()

# Uses the same loss and gain direction as the interface table
contact_table_df = contact_table_df.rename(columns={
    "Shared_Mapped_Contact_Count": "Shared_Contacts",
    "Lost_Mapped_Contact_Count": "Lost_in_Y2H_Positive",
    "Gained_Mapped_Contact_Count": "Gained_in_Y2H_Positive",
})

contact_table_df.to_csv(
    os.path.join(
        supp_table_dir,
        "table_s4_mapped_contacts.tsv",
    ),
    sep="\t",
    index=False,
)


# Keeps the FoldX energies for the same reliable comparisons
foldx_table_df = supp_result_df[
    common_cols + [
        "FoldX_Interaction_Energy_1",
        "FoldX_Interaction_Energy_2",
        "FoldX_Energy_Difference_Positive_Minus_Negative",
        "FoldX_Direction_Matches_Y2H",
    ]
].copy()

# Model 1 is Y2H negative and model 2 is Y2H positive
# The FoldX column names also include their energy units
foldx_table_df = foldx_table_df.rename(columns={
    "FoldX_Interaction_Energy_1": (
        "FoldX_Energy_Y2H_Negative_kcal_mol"
    ),
    "FoldX_Interaction_Energy_2": (
        "FoldX_Energy_Y2H_Positive_kcal_mol"
    ),
    "FoldX_Energy_Difference_Positive_Minus_Negative": (
        "FoldX_Difference_Y2H_Positive_Minus_Negative_kcal_mol"
    ),
})

foldx_table_df.to_csv(
    os.path.join(
        supp_table_dir,
        "table_s5_foldx.tsv",
    ),
    sep="\t",
    index=False,
)

# =============== SEQUENCE CHECK =====================

# Keeps every isoform and interactor sequence pair only once
sequence_pair_df = candidate_df[[
    "Base_Iso",
    "Iso_ID",
    "Inter_ID",
]].drop_duplicates()

# Will store the processed sequence information
sequence_rows = []

for _, sequence_row in sequence_pair_df.iterrows():
    iso_name = sequence_row["Iso_ID"]
    inter_name = sequence_row["Inter_ID"]

    # Paths of the processed sequences used in the analysis
    iso_path = os.path.join(
        output_longest_translated_fasta,
        f"{iso_name}.fasta",
    )
    inter_path = os.path.join(
        inter_filtered_dir,
        f"{inter_name}.fasta",
    )

    # Reads the complete processed FASTA sequences [5]
    iso_seq = str(
        SeqIO.read(iso_path, "fasta").seq
    ).rstrip("*")
    inter_seq = str(
        SeqIO.read(inter_path, "fasta").seq
    ).rstrip("*")

    # Same gene pairs have the same base isoform and interactor names
    same_gene = sequence_row["Base_Iso"] == inter_name

    # Complete sequence equality also requires the same length
    exact_match = iso_seq == inter_seq

    sequence_rows.append({
        "Iso_ID": iso_name,
        "Inter_ID": inter_name,
        "Iso_Length_AA": len(iso_seq),
        "Inter_Length_AA": len(inter_seq),
        "Same_Gene_Pair": "yes" if same_gene else "no",
        "Exact_Sequence_Match": "yes" if exact_match else "no",
    })

# Creates the processed sequence audit table
sequence_check_df = pd.DataFrame(sequence_rows)
sequence_check_df = sequence_check_df.sort_values([
    "Iso_ID",
    "Inter_ID",
])

sequence_check_df.to_csv(
    os.path.join(
        summary_output_dir,
        "sequence_pair_audit.tsv",
    ),
    sep="\t",
    index=False,
)

# Selects same-gene pairs with completely identical sequences
exact_sequence_df = sequence_check_df[
    (sequence_check_df["Same_Gene_Pair"] == "yes") &
    (sequence_check_df["Exact_Sequence_Match"] == "yes")
].copy()

print(f"Sequence pairs checked: {len(sequence_check_df)}")
print(f"Exact same-gene sequence matches: {len(exact_sequence_df)}")
print(exact_sequence_df)


# Finds the selected models belonging to the exact sequence matches
exact_model_df = best_model_df.merge(
    exact_sequence_df[["Iso_ID", "Inter_ID"]],
    on=["Iso_ID", "Inter_ID"],
)

# Will store the sequences found in the actual selected model chains
model_sequence_rows = []

for _, model_row in exact_model_df.iterrows():
    model_path = os.path.join(
        modeller_output_dir,
        model_row["Model_Job"],
        model_row["Model_Name"],
    )

    # Reads the selected MODELLER structure [13]
    model_structure = parser.get_structure(
        model_row["Model_Name"],
        model_path,
    )[0]

    # Reads the isoform chain first and interactor chain second
    model_sequences = []
    for chain_id in [
        model_row["Model_Iso_Chain"],
        model_row["Model_Inter_Chain"],
    ]:
        chain_seq = ""

        # Converts every amino acid in the model chain to one letter
        for res in model_structure[chain_id]:
            if is_aa(res, standard=False):
                chain_seq += seq1(
                    res.resname,
                    custom_map=aa_map,
                )

        model_sequences.append(chain_seq)

    # Checks whether the final model contains identical chains
    model_chain_match = model_sequences[0] == model_sequences[1]

    model_sequence_rows.append({
        "Iso_ID": model_row["Iso_ID"],
        "Inter_ID": model_row["Inter_ID"],
        "Model_Job": model_row["Model_Job"],
        "Selected_Model": model_row["Model_Name"],
        "Model_Iso_Chain": model_row["Model_Iso_Chain"],
        "Model_Inter_Chain": model_row["Model_Inter_Chain"],
        "Model_Iso_Length_AA": len(model_sequences[0]),
        "Model_Inter_Length_AA": len(model_sequences[1]),
        "Exact_Model_Chain_Match": (
            "yes" if model_chain_match else "no"
        ),
    })

# Creates the selected model chain audit table
model_sequence_check_df = pd.DataFrame(model_sequence_rows)
model_sequence_check_df = model_sequence_check_df.sort_values([
    "Iso_ID",
    "Inter_ID",
    "Model_Job",
])

model_sequence_check_df.to_csv(
    os.path.join(
        summary_output_dir,
        "selected_model_chain_audit.tsv",
    ),
    sep="\t",
    index=False,
)

print("Selected model chains for the exact sequence matches:")
print(model_sequence_check_df)





'''
######### REFERENCES ######### 
1 - https://www.uniprot.org/help/api_queries
2 - https://requests.readthedocs.io/en/latest/user/quickstart/
3 - https://bioinf.shenwei.me/seqkit/usage/
4 - https://stackoverflow.com/questions/66923711/how-to-parse-fasta-string-using-seq-io-from-biopython
5 - https://biopython.org/wiki/SeqIO
6 - https://stackoverflow.com/questions/72651199/finding-non-unique-rows-in-pandas-dataframe
7 - https://stackoverflow.com/questions/69663149/how-to-change-values-of-masked-column-in-pandas
8 - https://stackoverflow.com/questions/70079432/remove-suffix-if-string-matches-regular-expression-in-pandas
9 - https://stackoverflow.com/questions/55425324/pandas-drop-duplicates-based-on-2-columns-sometimes-reversed
10 - https://biopython.org/docs/latest/Tutorial/chapter_pairwise.html
11 - https://salilab.org/modeller/manual/node502.html
13 - https://biopython.org/docs/latest/Tutorial/chapter_pdb.html
14 - https://www.rcsb.org/docs/programmatic-access/file-download-services
15 - https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt
16 - https://doi.org/10.1038/nbt.3988
17 - https://www.mmseqs.com/latest/userguide.pdf
18 - https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.merge.html
19 - https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.groupby.html
20 - https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.drop_duplicates.html
21 - https://www.rcsb.org/structure/8RVE
22 - https://www.rcsb.org/structure/6WC2
23 - https://salilab.org/modeller/manual/node30.html
24 - https://salilab.org/modeller/manual/node52.html
25 - https://salilab.org/modeller/manual/node261.html
26 - https://salilab.org/modeller/manual/node206.html
27 - https://salilab.org/modeller/manual/node45.html
28 - https://salilab.org/modeller/manual/node46.html
29 - https://salilab.org/modeller/manual/node50.html
30 - https://salilab.org/modeller/manual/node51.html
31 - https://salilab.org/modeller/manual/node77.html
32 - https://biopython.org/docs/latest/api/Bio.PDB.Superimposer.html
33 - https://scikit-learn.org/stable/modules/generated/sklearn.metrics.jaccard_score.html
34 - https://foldxsuite.crg.eu/command/RepairPDB
35 - https://foldxsuite.crg.eu/command/AnalyseComplex
36 - https://www.wwpdb.org/documentation/file-format-content/format33/remarks2.html#REMARK%20350
37 - https://www.rcsb.org/structure/2XRF
38 - https://ars.els-cdn.com/content/image/1-s2.0-S0092867416300435-mmc2.xlsx
39 - https://ars.els-cdn.com/content/image/1-s2.0-S0092867416300435-mmc3.xlsx
'''
