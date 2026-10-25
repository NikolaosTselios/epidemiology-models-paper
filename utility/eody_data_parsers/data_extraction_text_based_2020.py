import pdfplumber
import re
import glob
import os
import pandas as pd

# Ο φάκελος που περιέχει τα PDF του 2020
folder_path = r"C:\Users\doublefool\Pandemic_Modelling\SIR_Model_Project_Matlab\eody\eody_pdf_2020"
pdf_files = glob.glob(os.path.join(folder_path, "*.pdf"))

# Regex για τη μορφή του 2020
patterns = {
    "Data_UpTo": r"μέχρι\s*τις\s*([\d]+\s*\S+\s*[\d]{4}\s*\(ώρα\s*\d+:\d+\))",
    "Confirmed_Cases": r"είναι\s*([\d\.]+)[,\.]?",
    "Total_Cases": r"ανέρχεται\s*σε\s*([\d\.]+)",
    "Daily_Change": r"ημερήσια\s*μεταβολή\s*\+([\d\.]+)%",
    "New_Deaths": r"COVID-19\s*είναι\s*([\d\.]+)",
    "Total_Deaths": r"συνολικά\s*([\d\.]+)\s*θάνατοι",
    "Median_Age_Deaths": r"είναι\s*([\d\.]+)\s*έτη"
}

# Λίστα για να συγκεντρωθούν τα αποτελέσματα
results = []

for pdf_file in pdf_files:
    print(f"Επεξεργασία αρχείου: {pdf_file}")
    try:
        with pdfplumber.open(pdf_file) as pdf:
            page_text = pdf.pages[0].extract_text()
            text = page_text if page_text else ""

            # Καθαρισμός κειμένου 
            text = text.replace('\n', '').replace('\r', '').replace('\x0c', '').strip()

            # Εξαγωγή των δεδομένων 
            record = {"File": os.path.basename(pdf_file)}
            for key, pattern in patterns.items():
                match = re.search(pattern, text, flags=re.IGNORECASE)
                record[key] = match.group(1).strip().rstrip(',.') if match else ""
            results.append(record)

    except Exception as e:
        print(f"Σφάλμα κατά την επεξεργασία του {pdf_file}: {e}")

# Αποθήκευση των αποτελεσμάτων σε CSV
df = pd.DataFrame(results)
output_csv = os.path.join(folder_path, "thefinal_extracted_data_specific2020.csv")
df.to_csv(output_csv, index=False)
print(f"Τα δεδομένα αποθηκεύτηκαν στο: {output_csv}")
