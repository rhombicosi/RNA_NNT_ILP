import RNA

# Example RNA sequence
sequence = "GCGCUUCGCCG"

# Fold it
(fc, mfe) = RNA.fold(sequence)

print("Sequence:", sequence)
print("Structure:", fc)
print("Minimum Free Energy:", mfe)