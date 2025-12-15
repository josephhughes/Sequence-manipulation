import argparse
import sys
from Bio import SeqIO
from Bio.Seq import Seq

def mismatch_count(s1, s2):
    """Calculates number of mismatched characters between two strings of equal length."""
    if len(s1) != len(s2):
        raise ValueError(f"Strings must be same length for mismatch count ({len(s1)} vs {len(s2)})")
    return sum(1 for a, b in zip(s1, s2) if a != b)

def remove_chars(s1, s2, chars_to_remove):
    """Removes characters from both strings if EITHER contains the char at that position."""
    # Perl logic: "if($char2=~/N/i || $char1=~/N/i) { ... } else { append }"
    # So if index i has N in s1 OR s2, drop from both.
    
    r1 = []
    r2 = []
    chars_set = set(chars_to_remove.lower())
    
    for a, b in zip(s1, s2):
        if a.lower() in chars_set or b.lower() in chars_set:
            continue
        r1.append(a)
        r2.append(b)
        
    return "".join(r1), "".join(r2)

def run_hamming(args):
    """Replicates hammingdistanceToRef.pl"""
    print(f"Calculating Hamming distances aligned to reference...")
    
    try:
        sequences = SeqIO.to_dict(SeqIO.parse(args.infile, "fasta"))
        
        # Determine Reference
        ref_seq = None
        if args.rfile:
            ref_recs = list(SeqIO.parse(args.rfile, "fasta"))
            if ref_recs:
                ref_seq = ref_recs[0]
                # Merge into main dict if not present (logic from Perl script)
                sequences[ref_seq.id] = ref_seq
        elif args.refname:
            if args.refname in sequences:
                ref_seq = sequences[args.refname]
            else:
                # Try finding by partial match?
                # Perl script uses hash lookup directly usually.
                print(f"Error: Reference '{args.refname}' not found in input file.")
                return
        else:
            print("Error: Reference must be specified via -ref or -rfile")
            return
            
        print(f"Reference: {ref_seq.id}")
        
        total_bases = 0
        total_seqs = 0
        total_mismatches = 0
        seq_with_mismatch = 0
        
        ref_str = str(ref_seq.seq)
        
        for seq_id, seq_rec in sequences.items():
            seq_str = str(seq_rec.seq)
            
            # Ensure same length or truncate/warn?
            # Perl script assumes same length.
            if len(seq_str) != len(ref_str):
                print(f"Warning: {seq_id} length {len(seq_str)} != Ref {len(ref_str)}. Truncating to min.")
                l = min(len(seq_str), len(ref_str))
                s1 = ref_str[:l]
                s2 = seq_str[:l]
            else:
                s1 = ref_str
                s2 = seq_str
                
            mm = mismatch_count(s1, s2)
            
            total_bases += len(s2)
            total_mismatches += mm
            if mm > 0:
                seq_with_mismatch += 1
            total_seqs += 1
            
        print(f"Total number of sequences = {total_seqs}")
        print(f"Total number of bases = {total_bases}")
        print(f"Total number of mismatches = {total_mismatches}")
        print(f"Total sequences with mismatches = {seq_with_mismatch}")
        
    except Exception as e:
        print(f"Error: {e}")

def run_compare(args):
    """Replicates compareRefs.pl"""
    print(f"Comparing sequences in {args.infile}...")
    
    try:
        sequences = list(SeqIO.parse(args.infile, "fasta"))
        if len(sequences) < 2:
            print("Error: Need at least 2 sequences to compare.")
            return
            
        s1 = sequences[0]
        s2 = sequences[1]
        
        str1 = str(s1.seq)
        str2 = str(s2.seq)
        
        # Assuming equal length
        l = min(len(str1), len(str2))
        
        diff = 0
        for i in range(l):
            if str1[i] != str2[i]:
                diff += 1
                
        print(f"Comparing {s1.id} vs {s2.id}")
        print(f"Total number of differences {diff}")
        
    except Exception as e:
        print(f"Error: {e}")

def run_sliding(args):
    """Replicates SlidingWindowNuc.pl"""
    print(f"Running Sliding Window Analysis on {args.infile} (Ref: {args.refname})...")
    
    try:
        if not args.outfile:
            print("Error: Output file required.")
            return
            
        sequences = SeqIO.to_dict(SeqIO.parse(args.infile, "fasta"))
        
        if args.refname not in sequences:
            print(f"Error: Reference {args.refname} not found.")
            return
            
        ref_rec = sequences[args.refname]
        
        # Region of interest
        # Perl start/end seems to apply to the subseq extraction BEFORE translation?
        # "read in the codon alignment... subseq($start, $end)"
        # Note: Perl is 1-based.
        
        start = args.start if args.start else 1
        end = args.end if args.end else len(ref_rec)
        
        # Adjust to 0-based
        s_idx = start - 1
        e_idx = end
        
        print(f"Region: {start}-{end}")
        
        # Store Nuc and AA sequences
        data = {}
        
        for sid, srec in sequences.items():
            # Extract region
            nuc_seq = srec.seq[s_idx:e_idx]
            # Translate
            # Assuming standard table, handling gaps? 
            # Perl script does: $orf->translate
            # If alignment is codon aligned, gaps should be in triplets usually?
            # We'll use gap='-'
            prot_seq = nuc_seq.translate(gap='-')
            
            data[sid] = {
                'nuc': str(nuc_seq),
                'aa': str(prot_seq)
            }
            
        # Analysis loop
        # For each seq vs ref
        with open(args.outfile, 'w') as out:
            out.write("Type\tQuery\tRef\tWindowIdx\tPosition\tMismatches\n")
            
            ref_data = data[args.refname]
            
            for sid, sdata in data.items():
                if sid == args.refname: continue
                
                # AA Sliding Window
                # Window size is in AA
                aa_seq = sdata['aa']
                ref_aa = ref_data['aa']
                w_size = args.window
                
                # Length check
                length_aa = min(len(aa_seq), len(ref_aa))
                
                for i in range(length_aa - w_size + 1):
                    w1 = aa_seq[i : i+w_size]
                    w2 = ref_aa[i : i+w_size]
                    
                    if args.remove:
                        w1, w2 = remove_chars(w1, w2, ['X', '-', '*']) # Remove X, gap, stop? Perl: X, -, N
                        
                    mm = mismatch_count(w1, w2)
                    pos = i * 3 + start # Mapping back to Nuc position?
                    out.write(f"AA\t{sid}\t{args.refname}\t{i}\t{pos}\t{mm}\n")
                    
                # Nuc Sliding Window
                # Window size for Nuc is 3 * AA_Window? 
                # Perl script: my $nucwindow=3*$window;
                nuc_seq = sdata['nuc']
                ref_nuc = ref_data['nuc']
                nw_size = 3 * w_size
                
                length_nuc = min(len(nuc_seq), len(ref_nuc))
                
                for j in range(length_nuc - nw_size + 1):
                    w1 = nuc_seq[j : j+nw_size]
                    w2 = ref_nuc[j : j+nw_size]
                    
                    if args.remove:
                        w1, w2 = remove_chars(w1, w2, ['N', '-', 'X'])
                        
                    mm = mismatch_count(w1, w2)
                    # Perl prints: $j and $j (start of window relative to extraction?)
                    # Let's print absolute pos?
                    pos = j + start
                    out.write(f"Nuc\t{sid}\t{args.refname}\t{j}\t{pos}\t{mm}\n")
                    
        print(f"Analysis written to {args.outfile}")

    except Exception as e:
        print(f"Error: {e}")

def main():
    parser = argparse.ArgumentParser(description="Distance and Comparison Tools")
    subparsers = parser.add_subparsers(dest='command', help='Sub-commands')
    
    # Hamming
    parser_ham = subparsers.add_parser('hamming', help='Calculate Hamming distance to reference')
    parser_ham.add_argument('-in', '--input', dest='infile', required=True, help="Input Alignment")
    parser_ham.add_argument('-ref', dest='refname', help="Reference Sequence ID")
    parser_ham.add_argument('-rfile', dest='rfile', help="Reference Sequence File")

    # Compare
    parser_comp = subparsers.add_parser('compare', help='Compare two sequences')
    parser_comp.add_argument('-in', '--input', dest='infile', required=True, help="Input File (2 sequences)")
    parser_comp.add_argument('-out', '--output', dest='outfile', help="Output File")

    # Sliding Window
    parser_slide = subparsers.add_parser('sliding', help='Sliding Window Analysis')
    parser_slide.add_argument('-in', '--input', dest='infile', required=True, help="Input Alignment")
    parser_slide.add_argument('-ref', dest='refname', required=True, help="Reference Sequence ID")
    parser_slide.add_argument('-w', '--window', dest='window', type=int, required=True, help="Window Size (AA)")
    parser_slide.add_argument('-start', type=int, help="Start position (1-based)")
    parser_slide.add_argument('-end', type=int, help="End position (1-based)")
    parser_slide.add_argument('-out', '--output', dest='outfile', required=True, help="Output TSV")
    parser_slide.add_argument('-remove', action='store_true', help="Remove Gaps/Ns/Xs before counting")

    args = parser.parse_args()
    
    if args.command == 'hamming':
        run_hamming(args)
    elif args.command == 'compare':
        run_compare(args)
    elif args.command == 'sliding':
        run_sliding(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
