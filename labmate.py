import sqlite3
import csv
import re

DB_NAME = "labplanner.db"

def connect_db():
    conn = sqlite3.connect(DB_NAME)
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
    sequence = sequence.upper()
    if not re.fullmatch(r"[ATGC]+", sequence):
        return("Fout: de sequentie mag alleen A, T, G en C bevatten.")
    a = sequence.count("A")
    t = sequence.count("T")
    g = sequence.count("G")
    c = sequence.count("C")
    tm = 2 * (a + t) + 4 * (g + c)
    return(f"Geschatte smelttemperatuur: {tm}°C\n")
