import pdfplumber
import re
import glob
import os
import pandas as pd

# Ο φάκελος που περιέχει τα PDF του 2022
folder_path = r"C:\Users\doublefool\Pandemic_Modelling\SIR_Model_Project_Matlab\eody\eody_pdf_2022"
pdf_files = glob.glob(os.path.join(folder_path, "*.pdf"))

# Regex για τη μορφή του 2022
patterns = {
    "Data_UpTo": r"μέχριτις([\d]{1,2}[^\d\(]+[\d]{4}\(ώρα\d+:\d+\))",
    "Confirmed_Cases": r"είναι([\d\.]+)",
    "Total_Cases":  r"ανέρχεται\s*σε\s*([\d\.]+)",
    "Daily_Change": r"ημερήσια\s*μεταβολή\s*\+([\d\.]+)%?",
    "New_Deaths": r"COVID[-]?[ ]?19\s*είναι\s*([\d\.]+)",
    "Total_Deaths": r"([\d\.]+)θάνατοι",
    "Median_Age_Deaths": r"([\d\.]+)ετών"
}

# Λίστα για να συγκεντρωθούν τα αποτελέσματα
results = []

for pdf_file in pdf_files:
    print("Επεξεργασία αρχείου:", pdf_file)
    try:
        with pdfplumber.open(pdf_file) as pdf:
            text = ""
            for i in range(min(3, len(pdf.pages))):
                page_text = pdf.pages[i].extract_text()
                if page_text:
                    text += page_text + "\n"

            # Καθαρισμός κειμένου 
            text = text.replace('\n', '').replace('\r', '').replace('\x0c', '').strip()

            # Εξαγωγή των δεδομένων
            record = {"File": os.path.basename(pdf_file)}
            for key, pattern in patterns.items():
                match = re.search(pattern, text, flags=re.IGNORECASE)
                record[key] = match.group(1).strip().rstrip(',.') if match else ""

            results.append(record)
    except Exception as e:
        print(f"Σφάλμα στο {pdf_file}: {e}")

# Αποθήκευση των αποτελεσμάτων σε CSV
df = pd.DataFrame(results)
output_csv = os.path.join(folder_path, "final_extracted_data_specific2022.csv")
df.to_csv(output_csv, index=False)
print("Τα δεδομένα αποθηκεύτηκαν στο", output_csv)
