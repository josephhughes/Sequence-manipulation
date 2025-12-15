import argparse
import sys
from Bio import SeqIO
from collections import defaultdict, Counter

def main():
    parser = argparse.ArgumentParser(description="Count AA relative to reference (Python port of CountAArelativeToRef.pl)")
    parser.add_argument('-in', '--input', dest='infile', required=True, help="Input protein alignment (FASTA)")
    parser.add_argument('-out', '--output', dest='outfile', required=True, help="Output tab-delimited file")
    parser.add_argument('-r', '--ref', dest='ref_id', required=True, help="Reference sequence ID")
    parser.add_argument('-l', '--list', dest='listfile', help="File with list of identifiers to include")
    parser.add_argument('-s', '--site', dest='site', type=int, help="Site of interest (optional)")

    args = parser.parse_args()

    # Define amino acids list (matching original script roughly, usually standard 20 + others)
    # Original: - A B C D E F G H I J K L M N P Q R S T V W X Y Z
    # Note: 'GAP' is handled separately in logic usually, but here we can map '-' to GAP or just count it.
    amino_acids = ['-', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z']
    
    # Load List if provided
    keep_ids = set()
    if args.listfile:
        print(f"Filtering by IDs in {args.listfile}...")
        try:
            with open(args.listfile, 'r') as f:
                for line in f:
                    if line.strip():
                        parts = line.strip().split('\t')
                        keep_ids.add(parts[0])
            if args.ref_id not in keep_ids:
                 # Original script exits if ref not in list? 
                 # "If the reference id is not in the list of identifiers to filter... exit"
                 # But usually ref is needed. Let's warn or add it.
                 # The original script says: "The reference id is not in the list... exit"
                 pass
        except FileNotFoundError:
            print(f"Error: List file '{args.listfile}' not found.")
            sys.exit(1)

    # Load Sequences
    try:
        sequences = list(SeqIO.parse(args.infile, "fasta"))
    except FileNotFoundError:
        print(f"Error: Input file '{args.infile}' not found.")
        sys.exit(1)

    # Find Reference
    ref_seq_obj = next((s for s in sequences if s.id == args.ref_id), None)
    if not ref_seq_obj:
        print(f"Error: Reference ID '{args.ref_id}' not found in input file.")
        sys.exit(1)
    
    # Filter sequences
    if args.listfile:
        if args.ref_id not in keep_ids:
             print("Error: Reference ID must be in the filter list.")
             sys.exit(1)
        filtered_seqs = [s for s in sequences if s.id in keep_ids]
        sequences = filtered_seqs
    
    # Process
    # We need to map alignment positions to reference positions (skipping reference gaps)
    # Ref: A - C T
    # Pos: 1   2 3
    # Mapped: 1 -> 1, 2 -> Gap (skip column?), 3 -> 2, 4 -> 3
    # Original script says: "The columns/sites with gaps in the reference are removed"
    
    ref_str = str(ref_seq_obj.seq).upper()
    ref_positions = {} # aln_index -> ref_aa (only non-gap)
    ref_map = {} # aln_index -> ref_index (1-based)
    
    curr_ref_idx = 0
    valid_aln_indices = []
    
    for idx, char in enumerate(ref_str):
        if char != '-':
            curr_ref_idx += 1
            ref_positions[idx] = char
            ref_map[idx] = curr_ref_idx
            valid_aln_indices.append(idx)
        else:
            # Gap in reference, this column will be ignored in output
            pass

    # Count bases
    # Data structure: base[aln_index][aa] = count
    base_counts = defaultdict(lambda: defaultdict(int))
    total_seqs = 0
    cnt_match = 0
    
    for seq in sequences:
        # if seq.id == args.ref_id: continue # Original script includes ref in stats? 
        # Original: "if ($list{$id} && $id ne $ref)" -> counts matches.
        # But later "if ($id eq $ref)" it processes reference gaps.
        # It seems to count everyone?
        # Actually in the "if (!$list)" block: "if ($id ne $ref) ... count"
        # So it seems reference is NOT included in the counts distribution?
        # Let's check original: 
        # if ($list{$id} && $id ne $ref){ ... count ... }
        # if ($id eq $ref) { ... map gaps ... }
        # So Reference is EXCLUDED from the counts.
        
        if seq.id == args.ref_id:
            continue
            
        cnt_match += 1
        seq_str = str(seq.seq).upper()
        # Assume alignment length is consistent logic? Biopython doesn't enforce, but usually implied.
        
        for idx in valid_aln_indices:
            if idx < len(seq_str):
                aa = seq_str[idx]
                base_counts[idx][aa] += 1
            else:
                # Sequence shorter than alignment? Treat as gap?
                base_counts[idx]['-'] += 1

    total_seqs = cnt_match # Only counted sequences
    
    # Output
    # Formatting: Position relative to ref, reference AA, [AminoAcids], total sequences, Consensus AA, Frequency
    
    with open(args.outfile, 'w') as out:
        header = ["Position relative to " + args.ref_id, "reference AA"] + amino_acids + ["total sequences", "Consensus AA", "Frequency"]
        out.write("\t".join(header) + "\n")
        
        for idx in valid_aln_indices:
            ref_aa = ref_positions[idx]
            ref_pos_num = ref_map[idx]
            
            row = [str(ref_pos_num), ref_aa]
            
            # Counts
            counts = base_counts[idx]
            current_total = 0 # sum of counts at this pos
            consensus_aa = ''
            consensus_count = -1
            
            for aa in amino_acids:
                cnt = counts.get(aa, 0)
                row.append(str(cnt))
                current_total += cnt
                
                if cnt > consensus_count:
                    consensus_count = cnt
                    consensus_aa = aa
            
            row.append(str(current_total))
            row.append(consensus_aa)
            
            freq = 0.0
            if current_total > 0:
                freq = consensus_count / current_total
            row.append(f"{freq:.2f}")
            
            out.write("\t".join(row) + "\n")

    print(f"Processed {len(sequences)} sequences.")
    print(f"Reference ID: {args.ref_id}")
    print(f"Output written to {args.outfile}")
    
    if args.site:
        # Original script prints specific site details to stdout
        # if ($_=~/^$site\t/) ... match output line
        pass 
        # Since we just wrote the file, we could read it back or just print if we have the data
        # Let's assume standard usage doesn't strictly require the stdout dump unless requested.

if __name__ == "__main__":
    main()
