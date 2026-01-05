#!/usr/bin/awk -f

BEGIN {
    FS = ","; OFS = "\t";
}

NR == 1 {
    # Process header to dynamically find variable and R-variable pairs
    for (i = 1; i <= NF; i++) {
        if ($i ~ /^R-/) {
            var = substr($i, 3); # Extract variable name (e.g., Milk from R-Milk)
            r_index[var] = i;    # Store index of R-variable
        } else {
            main_index[$i] = i;  # Store index of main variable
        }
    }
    # Prepare the output header
    header = "#IID";
    for (var in r_index) {
        header = header OFS var;
    }
    print header;
}

NR > 1 {
    output = $1; # Start output with the ID
    for (var in r_index) {
        # Calculate value / (R * (1 - R)) if R is valid
        r_val = $(r_index[var]);
        if (r_val > 0 && r_val < 1) {
            calc = $(main_index[var]) / (r_val * (1 - r_val));
        } else {
            calc = "N/A"; # Handle invalid R values (e.g., 0 or 1)
        }
        output = output OFS calc;
    }
    print output;
}