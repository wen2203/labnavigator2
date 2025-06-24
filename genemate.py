from Bio import Entrez # https://biopython.org/docs/1.76/api/Bio.Entrez.html

# haalt gen uit ncbi database
def fetch_gene_data(gene_name, organism):
    try:
        # zoekt het gen in database
        handle = Entrez.esearch(db="nucleotide", term=f"{gene_name}[Gene] AND {organism}[Organism]", retmax=1)
        record = Entrez.read(handle)
        ids = record.get("IdList", [])
        if not ids:
            return None
        
        # haalt fasta data op voor het eerste resultaat uit database
        # dit kan nog later uitgebreid worden zodat gebruikers uit genen kunnen zoeken
        fetch_handle = Entrez.efetch(db="nucleotide", id=ids[0], rettype="fasta", retmode="text")

        
        fasta_data = fetch_handle.read()
        return fasta_data if fasta_data.strip() else None
    except Exception as e:
        return None
