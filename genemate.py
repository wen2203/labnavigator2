from Bio import SeqIO
import os
import requests
from Bio import Entrez

def fetch_gene_data(gene_name, organism):
    try:
        handle = Entrez.esearch(db="nucleotide", term=f"{gene_name}[Gene] AND {organism}[Organism]", retmax=1)
        record = Entrez.read(handle)
        ids = record.get("IdList", [])
        if not ids:
            return None
        fetch_handle = Entrez.efetch(db="nucleotide", id=ids[0], rettype="fasta", retmode="text")
        fasta_data = fetch_handle.read()
        return fasta_data if fasta_data.strip() else None
    except Exception as e:
        return None
