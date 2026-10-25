function periods = generate_fitting_periods(start_date_str, end_date_str, step_days)
% generate_fitting_periods - Δημιουργεί διαδοχικά χρονικά παράθυρα fitting με βάση αρχική ημερομηνία και βήμα.

% Μετατροπή εισόδων σε datetime
start_date = datetime(start_date_str);
end_date   = datetime(end_date_str);

% Αρχικοποίηση μεταβλητών
current_start = start_date;
periods = [];

% Δημιουργία περιόδων χωρίς overlapping
while current_start <= end_date
    current_end = current_start + days(step_days - 1);
    if current_end > end_date
        current_end = end_date;
    end

    % Προσθήκη τρέχουσας περιόδου στο αποτέλεσμα
    periods = [periods; string(current_start), string(current_end)];

    % Ενημέρωση για το επόμενο χρονικό παράθυρο
    current_start = current_end + days(1);
end

end