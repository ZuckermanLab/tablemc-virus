($0!~/max_total_iterations/) {print $0; next}
($0~/max_total_iterations/) {loc=index($0,$NF); $0=substr($0,1,loc-1) ($NF+nadd); print $0}
