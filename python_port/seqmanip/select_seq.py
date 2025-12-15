import argparse
import random
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main():
    parser = argparse.ArgumentParser(description="Python port of SelectSeq.pl for sequence manipulation.")
    
    # Input/Output
    parser.add_argument('-in', '--input', dest='infile', required=True, help="Input FASTA file")
    parser.add_argument('-out', '--output', dest='outfile', required=True, help="Output FASTA file")
    
    # Matching options
    parser.add_argument('-fid', dest='fid', help="Exact ID match")
    parser.add_argument('-mid', dest='mid', help="Partial ID match (regex)")
    parser.add_argument('-rmid', dest='rmid', help="Remove specific ID match")
    parser.add_argument('-idfile', dest='idfile', help="File with list of IDs to process")
    parser.add_argument('-nomatchid', dest='nomatchid', default="NoMatch_ids.txt", help="Output file for unmatched IDs")
    parser.add_argument('-fdescr', dest='fdescr', help="Partial description match (regex)")
    parser.add_argument('-rdm', dest='rdm', nargs='?', const=1, type=int, help="Randomly select N sequences (default: 1)")
    
    # Specialized modes
    parser.add_argument('-inout', action='store_true', help="Output 'In' (matches) and 'Out' (non-matches) files")
    parser.add_argument('-singletons', action='store_true', help="Split into singletons and non-singletons")
    parser.add_argument('-dedup', action='store_true', help="Deduplicate by ID (keep first)")
    parser.add_argument('-longest', action='store_true', help="Output the longest sequence")
    parser.add_argument('-freq', dest='frequency', type=float, help="Frequency threshold")
    
    # ID modifiers
    parser.add_argument('-cogukid', action='store_true', help="Parse COG-UK IDs")
    parser.add_argument('-coviz', action='store_true', help="Parse CoViz IDs")

    args = parser.parse_args()

    # Read sequences
    try:
        sequences = list(SeqIO.parse(args.infile, "fasta"))
    except FileNotFoundError:
        print(f"Error: Input file '{args.infile}' not found.")
        sys.exit(1)

    # --- Mode: List of IDs (-idfile) ---
    if args.idfile:
        print(f"Looking for IDs from {args.idfile} in {args.infile}...")
        try:
            with open(args.idfile, 'r') as f:
                target_ids = {line.strip().split()[0] for line in f if line.strip()}
        except FileNotFoundError:
            print(f"Error: ID file '{args.idfile}' not found.")
            sys.exit(1)
            
        matches = []
        non_matches = []
        hits = 0
        
        # Helper to process IDs for specific formats
        def process_id(seq_id):
            if args.cogukid:
                # Approximate logic from Perl: s/^COGUK\/(.+)\/\w+.+/$1/
                if seq_id.startswith("COGUK/"):
                    parts = seq_id.split('/')
                    if len(parts) > 1: return parts[1]
                # s/^\w+\/(.+)\/\d{4}/$1/
                parts = seq_id.split('/')
                if len(parts) >= 3: return parts[1]
            return seq_id # Add coviz support if needed pattern is clear
            
        for seq in sequences:
            curr_id = process_id(seq.id)
            if curr_id in target_ids:
                matches.append(seq)
                hits += 1
                target_ids.discard(curr_id) # Consuming the ID? Perl deletes from list.
            else:
                non_matches.append(seq)

        if args.inout:
             SeqIO.write(matches, f"In{args.outfile}", "fasta")
             SeqIO.write(non_matches, f"Not{args.outfile}", "fasta")
        else:
             SeqIO.write(matches, args.outfile, "fasta")
             
        print(f"Found {hits} matches.")
        if target_ids:
            print(f"No matches for: {', '.join(target_ids)}")
            with open(args.nomatchid, 'w') as f:
                f.write('\n'.join(target_ids))

    # --- Mode: Exact Match (-fid) ---
    elif args.fid:
        print(f"Looking for exact ID: {args.fid}")
        matches = [s for s in sequences if s.id == args.fid]
        SeqIO.write(matches, args.outfile, "fasta")
        print(f"Found {len(matches)} matches.")

    # --- Mode: Partial Match (-mid) ---
    elif args.mid:
        import re
        print(f"Looking for partial ID match: {args.mid}")
        pattern = re.compile(args.mid)
        matches = [s for s in sequences if pattern.search(s.id)]
        SeqIO.write(matches, args.outfile, "fasta")
        print(f"Found {len(matches)} matches.")

    # --- Mode: Remove ID (-rmid) ---
    elif args.rmid:
        print(f"Removing ID: {args.rmid}")
        kept = [s for s in sequences if s.id != args.rmid]
        SeqIO.write(kept, f"{args.outfile}.noid{args.rmid}.fa", "fasta")

    # --- Mode: Description Match (-fdescr) ---
    elif args.fdescr:
        import re
        print(f"Looking for description match: {args.fdescr}")
        pattern = re.compile(args.fdescr)
        matches = [s for s in sequences if s.description and pattern.search(s.description)]
        SeqIO.write(matches, args.outfile, "fasta")
        print(f"Found {len(matches)} matches.")

    # --- Mode: Random (-rdm) ---
    elif args.rdm is not None:
        count = args.rdm
        print(f"Selecting {count} random sequences...")
        if sequences:
            if count > len(sequences):
                print(f"Warning: Requested {count} sequences but only {len(sequences)} available. Returning all.")
                selected = sequences
            else:
                selected = random.sample(sequences, count)
                
            print(f"Selected: {len(selected)} sequences")
            # If only 1, print ID for convenience (legacy behavior compatibility)
            if len(selected) == 1:
                print(f"Selected ID: {selected[0].id}")
                
            SeqIO.write(selected, args.outfile, "fasta")
        else:
            print("No sequences found to select from.")

    # --- Mode: Singletons (-singletons) ---
    elif args.singletons:
        print("Splitting singletons...")
        seq_map = {}
        for s in sequences:
            seq_str = str(s.seq)
            if seq_str not in seq_map: seq_map[seq_str] = []
            seq_map[seq_str].append(s)
            
        singletons = []
        nosingletons = []
        
        for s_str, records in seq_map.items():
            if len(records) > 1:
                nosingletons.extend(records)
            else:
                singletons.extend(records)
                
        SeqIO.write(singletons, f"{args.outfile}.singletons.fa", "fasta")
        SeqIO.write(nosingletons, f"{args.outfile}.nosingletons.fa", "fasta")

    # --- Mode: Deduplicate (-dedup) ---
    elif args.dedup:
        print("Deduplicating by ID...")
        seen = set()
        kept = []
        for s in sequences:
            if s.id not in seen:
                seen.add(s.id)
                kept.append(s)
            else:
                print(f"Duplicate dropped: {s.id}")
        SeqIO.write(kept, f"{args.outfile}.nodups.fa", "fasta")

    # --- Mode: Longest (-longest) ---
    elif args.longest:
        print("Finding longest sequence...")
        if not sequences:
            print("No sequences found.")
            return
        longest = max(sequences, key=lambda s: len(s.seq))
        print(f"Longest: {longest.id} ({len(longest.seq)} bp)")
        SeqIO.write([longest], f"{args.outfile}_longest.fa", "fasta")

    # --- Mode: Frequency (-freq) ---
    elif args.frequency is not None:
        print(f"Filtering by frequency threshold: {args.frequency}")
        total = len(sequences)
        if total == 0: return
        
        seq_counts = {}
        for s in sequences:
            s_str = str(s.seq)
            if s_str not in seq_counts: seq_counts[s_str] = []
            seq_counts[s_str].append(s)
            
        above = []
        below = []
        
        for s_str, records in seq_counts.items():
            count = len(records)
            freq = count / total
            # print(f"{count} {freq} {s_str[:10]}...") 
            if freq > args.frequency:
                above.extend(records)
            else:
                below.extend(records)
                
        SeqIO.write(above, f"{args.outfile}.above{args.frequency}.fa", "fasta")
        SeqIO.write(below, f"{args.outfile}.below{args.frequency}.fa", "fasta")

if __name__ == "__main__":
    main()
