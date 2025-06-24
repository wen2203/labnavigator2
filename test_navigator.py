# importeer functies die je wilt testen
from labmate import connect_db, calculate_tm
from genemate import fetch_gene_data
from unittest.mock import patch, MagicMock # https://docs.python.org/3/library/unittest.mock.html


def test_connect_db_creates_table(): # deze unittest is door Chatgpt geschreven
    # Roep de connect_db() functie aan, deze maakt (indien nodig) de database en tabel aan
    conn, c = connect_db()
    # Controleer of de tabel 'experiments' bestaat in de SQLite-database
    c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='experiments'")
    # Haal het resultaat op van de query
    result = c.fetchone()
    # Sluit de database
    conn.close()
    # Test slaagt als result niet None is (de tabel bestaat dus)
    assert result is not None


def test_calculate_tm_valid_sequence():
    # checkt sequentie of er een correcte output is
    assert calculate_tm("ATGC") == (f"Geschatte smelttemperatuur: {12}°C\n")
    assert calculate_tm("atgc") == (f"Geschatte smelttemperatuur: {12}°C\n")
    # check of de nucleotiden de juiste temp geeft en of de formule klopt
    assert calculate_tm("a") == (f"Geschatte smelttemperatuur: {2}°C\n")
    assert calculate_tm("t") == (f"Geschatte smelttemperatuur: {2}°C\n")
    assert calculate_tm("g") == (f"Geschatte smelttemperatuur: {4}°C\n")
    assert calculate_tm("c") == (f"Geschatte smelttemperatuur: {4}°C\n")


def test_calculate_tm_invalid_sequence():
    # checkt of de verkeerde input een fout geeft
    assert calculate_tm("AXTG") == "Fout: de sequentie mag alleen A, T, G en C bevatten."
    assert calculate_tm("hello") == "Fout: de sequentie mag alleen A, T, G en C bevatten."


# deze code is door Chatgpt geschreven
@patch("genemate.Entrez.esearch")
@patch("genemate.Entrez.read")
@patch("genemate.Entrez.efetch")
def test_fetch_gene_data_found(mock_efetch, mock_read, mock_esearch):
    # Simuleer een gevonden gen
    mock_read.return_value = {"IdList": ["123456"]}
    mock_fetch = MagicMock()
    mock_fetch.read.return_value = ">MockGen\nATGC"
    mock_efetch.return_value = mock_fetch

    result = fetch_gene_data("BRCA1", "Homo sapiens")
    assert result is not None
    assert "ATGC" in result


# deze code is door Chatgpt geschreven
@patch("genemate.Entrez.esearch")
@patch("genemate.Entrez.read")
def test_fetch_gene_data_not_found(mock_read, mock_esearch):
    # Simuleer geen resultaten
    mock_read.return_value = {"IdList": []}
    result = fetch_gene_data("FakeGene", "Homo sapiens")
    assert result is None
