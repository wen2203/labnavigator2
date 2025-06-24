import streamlit as st  # import Streamlit om een webapp te maken https://docs.streamlit.io/develop/api-reference
import os  # Import os voor bestands- en mapbewerkingen
from labmate import connect_db, calculate_tm  # zelfgemaakte functies: verbinden met DB, smelttemp bereken
from genemate import fetch_gene_data

from Bio import SeqIO
from Bio import Entrez
import io
import csv

# Verbind met de database en krijg cursor om sql queries uit te voeren
# conn is de verbinding en c is de pen n queries uit te voeren
conn, c = connect_db()

st.image("Untitled_Artwork.png", width=500) # voeg een foto toe als logo en


# titel bovenaan de webapp dit is met Chatgpt gedaan voor de kleur
st.markdown("<h1 style='color:hotpink;'>✧˖✧ Lab Navigator ✧˖✧</h1>", unsafe_allow_html=True)

# kleine alinea over Lab Navigator
st.write("Welkom bij Lab Navigator! " \
"Dit is een hulpmiddel om eenvoudig het lab te navigeren. " \
"Maak een keuze welke tool je wilt gebruiken, zoals een planner om experimenten in te plannen," \
" Gene Fetcher om gensequenties te downloaden, een FASTQ naar FASTA converter en de smelttemperatuur-rekenmachine.")


# een dropdownmenu/selectbox waarin de gebruiker een optie kiest
menu = st.sidebar.radio("Maak een keuze uit de toolbox:", [
    "Nieuw experiment",
    "Bekijk experimenten",
    "Rond experiment af",
    "Verwijder experiment",
    "Exporteer CSV",
    "Smelttemperatuur berekenen",
    "Convert FASTQ → FASTA",
    "Gen downloaden (NCBI database)"
])

# voor elke optie wordt iets anders getoond en uitgevoerd:

# wanneer er nieuw exp word geselecteerd
if menu == "Nieuw experiment":
    st.markdown("<h3 style='color:deeppink;'>Nieuw experiment toevoegen</h3>", unsafe_allow_html=True)
    # Hier vraagt de app om gegevens van het experiment
    name = st.text_input("Experimentnaam")
    date = st.date_input("Datum:")
    time = st.time_input("Tijd:")
    
    duration = st.number_input("Duur in minuten", min_value=1)  # getal voor duur
    user = st.text_input("Gebruiker")  # naam van de persoon

    
    # als gebruiker op knop toevoegen klikt:
    if st.button("Toevoegen"):
        # Zet date en time om naar string formats
        date_str = date.strftime("%Y-%m-%d")
        time_str = time.strftime("%H:%M")

        # voeg de ingevoerde data toe in de database
        c.execute('''
            INSERT INTO experiments (name, date, start_time, duration, user)
            VALUES (?, ?, ?, ?, ?)
        ''', (name, date_str, time_str, duration, user))
        # sla de wijziging op in DB
        conn.commit()
        # laat een succesbericht zien
        st.success("Experiment toegevoegd!")


elif menu == "Bekijk experimenten":
    st.markdown("<h3 style='color:deeppink;'>Experimentenlijst</h3>", unsafe_allow_html=True)
    c.execute("SELECT * FROM experiments") #selecteer  alles van experimenten tabel
    # fetch alle rows
    rows = c.fetchall()
    if rows:
        # loop gaat door alle rows
        for row in rows:
            # gaat door alle rows en schrijft info per row op
            st.write(f"ID: {row[0]} | Naam: {row[1]} | Datum: {row[2]} | Start: {row[3]} | "
                     f"Duur: {row[4]} min | Gebruiker: {row[5]} | Status: {row[6]}")
    # geen rows is geen experimenten gevonden
    else:
        st.info("Geen experimenten gevonden.")


elif menu == "Rond experiment af":
    st.markdown("<h3 style='color:deeppink;'>Markeer experiment als afgerond</h3>", unsafe_allow_html=True)
    c.execute("SELECT id, name, date, start_time, duration, user, status FROM experiments WHERE status != 'done'")
    rows = c.fetchall()

    if rows:
        opties = [f"{row[0]} - {row[1]} ({row[2]}) Status: {row[6]}" for row in rows]
        keuze = st.selectbox("Selecteer experiment om af te ronden", opties)
        exp_id = int(keuze.split(" - ")[0])

        if st.button("Afronden"):
            c.execute("UPDATE experiments SET status='done' WHERE id=?", (exp_id,))
            conn.commit()
            st.success(f"Experiment {exp_id} gemarkeerd als afgerond.")
    else:
        st.info("Geen openstaande experimenten gevonden.")


elif menu == "Verwijder experiment":
    st.markdown("<h3 style='color:deeppink;'>Verwijder experiment</h3>", unsafe_allow_html=True)
    c.execute("SELECT id, name, date, start_time, duration, user, status FROM experiments")
    rows = c.fetchall()
    if rows:
        opties = [f"{row[0]} - {row[1]} ({row[2]}) Status: {row[6]}" for row in rows]
        keuze = st.selectbox("Selecteer experiment om te verwijderen", opties)
        exp_id = int(keuze.split(" - ")[0])
    
        if st.button("Verwijderen"):
            c.execute("DELETE FROM experiments WHERE id=?", (exp_id,))
            conn.commit()
            st.success(f"Experiment {exp_id} verwijderd.")
    else:
        st.info("Geen openstaande experimenten gevonden.")


elif menu == "Exporteer CSV":
    
    st.markdown("<h3 style='color:deeppink;'>Exporteer experimenten naar CSV</h3>", unsafe_allow_html=True)

    if st.button("Exporteer CSV"):
        # Data ophalen uit de database
        c.execute("SELECT * FROM experiments")
        rows = c.fetchall()

        # Bestand op schijf schrijven
        output_dir = "output"
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, "experiments_export.csv")
        with open(output_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(['ID','Naam','Datum','Starttijd','Duur','Gebruiker','Materialen','Locatie','Status'])
            writer.writerows(rows)

        st.success("✅ CSV-bestand gemaakt!")

        # Lees het bestand als bytes voor download
        with open(output_path, 'rb') as f:
            file_bytes = f.read()

        # Downloadknop met de echte file-inhoud
        st.download_button(
            label="⬇️ Download CSV-bestand",
            data=file_bytes,
            file_name="experiments.csv",
            mime="text/csv"
        )


elif menu == "Smelttemperatuur berekenen":
    st.markdown("<h3 style='color:deeppink;'>Smelttemperatuur berekenen</h3>", unsafe_allow_html=True)

    seq = st.text_input("Voer DNA-sequentie in (A,T,G,C):")
    if st.button("Bereken"):
        resultaat = calculate_tm(seq)
        st.write(resultaat)


elif menu == "Convert FASTQ to FASTA":
    #st.markdown("<h3 style='color:deeppink;'>Convert FASTQ to FASTA</h3>", unsafe_allow_html=True)
    uploaded_file = st.file_uploader("Upload FASTQ-bestand", type=["fastq"])
    
    if uploaded_file is not None:
        out_name = st.text_input("Naam outputbestand (zonder .fasta):", value="output")
    
        if st.button("Converteren"):
            if not out_name.endswith(".fasta"):
                out_name += ".fasta"
    
            # use BytesIO for uploaded file
            fastq_io = io.BytesIO(uploaded_file.getvalue())
            fasta_io = io.StringIO()
    
            SeqIO.convert(fastq_io, "fastq", fasta_io, "fasta")
            fasta_str = fasta_io.getvalue()
    
            st.success("✅ Conversie gelukt!")
            st.download_button(
                label="⬇️ Download FASTA bestand",
                data=fasta_str,
                file_name=out_name,
                mime="text/plain"
            )


elif menu == "Gen downloaden (NCBI database)":
    Entrez.email = "wendypan22@hotmail.com"  

    st.markdown("<h3 style='color:deeppink;'>Gen downloaden van NCBI</h3>", unsafe_allow_html=True)
    gene = st.text_input("Gennaam:")
    organism = st.selectbox("Organisme:", [
        "Homo sapiens",
        "Mus musculus",
        "Rattus norvegicus",
        "Danio rerio",
        "Gallus gallus domesticus",
        "Oryctolagus cuniculus",
        "Sus scrofa domesticus",
        "Macaca mulatta",
        "Cavia porcellus",
        "Drosophila melanogaster",
        "Gorilla gorilla"
    ])

    
    if st.button("Downloaden"): # als user op de download knop drukt
        data = fetch_gene_data(gene, organism) # roept eigen functie aan 
        if data: # als er data is 
            st.success(f"✅ Gen {gene} voor {organism} gevonden!") # geef aan dat het gen is gevonden aan gebruiker
            # dit zorgt ervoor dat het gen gedownload word
            st.download_button(
                label="⬇️ Download FASTA bestand",
                data=data,
                file_name=f"{gene}.fasta",
                mime="text/plain"
            )
        # als er geen gen is gevonden
        else:
            st.error(f"❌ Gen '{gene}' niet gevonden of download mislukt.") # geeft aan dat het niet is gelukt

