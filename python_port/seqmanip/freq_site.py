import argparse
import sys
from Bio import SeqIO

def run_freq_site(args):
    """Replicates FreqOfSite.pl"""
    print(f"Calculating frequency at site {args.site} relative to reference '{args.ref}' in {args.infile}...")
    
    try:
        sequences = list(SeqIO.parse(args.infile, "fasta"))
        
        # 1. Find Reference
        ref_seq = None
        for seq in sequences:
            if args.ref in seq.id:
                ref_seq = seq
                print(f"Found reference: {seq.id}")
                break
        
        if not ref_seq:
            print(f"Error: Reference '{args.ref}' not found in alignment.")
            return

        # 2. Map Ungapped Site to Aligned Column
        # FreqOfSite.pl logic:
        # Ungapped site 1 = First non-gap character.
        
        req_site = args.site
        aligned_col = -1
        
        ungapped_count = 0
        ref_str = str(ref_seq.seq)
        
        for i, char in enumerate(ref_str):
            if char != '-':
                ungapped_count += 1
            
            if ungapped_count == req_site:
                aligned_col = i
                break
                
        if aligned_col == -1:
            print(f"Error: Site {req_site} exceeds reference length (ungapped length: {ungapped_count}).")
            return
            
        print(f"Ungapped Site {req_site} maps to Alignment Column {aligned_col + 1} (1-based index: {aligned_col})")
        
        # 3. Count Nucleotides at Column
        counts = {'A': 0, 'C': 0, 'T': 0, 'G': 0}
        
        for seq in sequences:
            s_str = str(seq.seq)
            if aligned_col < len(s_str):
                base = s_str[aligned_col].upper()
                if base in counts:
                    counts[base] += 1
                # Optional: Count gaps or Ns if needed? Original script prints As, Cs, Ts, Gs using a matrix hash.
                
        # 4. Output
        print(f"Site {req_site} Reference {ref_str[aligned_col]}")
        print(f"Site with gaps {aligned_col + 1}")
        print(f"As: {counts['A']}")
        print(f"Cs: {counts['C']}")
        print(f"Ts: {counts['T']}")
        print(f"Gs: {counts['G']}")
            
    except Exception as e:
        print(f"Error: {e}")

def main():
    parser = argparse.ArgumentParser(description="Frequency of Site Tool")
    parser.add_argument('-in', '--input', dest='infile', required=True, help="Input Alignment FASTA")
    parser.add_argument('-ref', dest='ref', required=True, help="Reference Sequence ID (partial match)")
    parser.add_argument('-site', type=int, required=True, help="Site of interest (ungapped, 1-based)")
    
    args = parser.parse_args()
    
    run_freq_site(args)

if __name__ == "__main__":
    main()
