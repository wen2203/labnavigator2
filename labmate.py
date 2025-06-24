# om met sql te werken
import sqlite3
# voor regex 
import re

def connect_db():
    # connect met databese
    conn = sqlite3.connect("labplanner.db")
    # maakt een cursor aan om sql queries uit te voeren
    c = conn.cursor()
    c.execute('''
        CREATE TABLE IF NOT EXISTS experiments (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT NOT NULL,
            date TEXT NOT NULL,
            start_time TEXT NOT NULL,
            duration INTEGER NOT NULL,
            user TEXT NOT NULL,
            status TEXT NOT NULL DEFAULT 'pending'
        )
    ''')
    conn.commit()
    return conn, c

def calculate_tm(sequence):
    # maakt sequentie uppercase
    sequence = sequence.upper()
    # als het niet matcht met A, T, G OF C dan fout
    if not re.fullmatch(r"[ATGC]+", sequence):
        return("Fout: de sequentie mag alleen A, T, G en C bevatten.")

    # tel de nucleotiden
    a = sequence.count("A")
    t = sequence.count("T")
    g = sequence.count("G")
    c = sequence.count("C")

    # smelttemperatuur formule
    tm = 2 * (a + t) + 4 * (g + c)

    # returnt temperatuur met string
    return(f"Geschatte smelttemperatuur: {tm}°C\n")
