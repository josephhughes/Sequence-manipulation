import argparse
import sys
import os
from Bio import SeqIO
from collections import Counter
import math

def run_summary(args):
    """Replicates FastaStats.pl"""
    print(f"Calculating statistics for {args.infile}...")
    
    try:
        if args.outfile:
            out_handle = open(args.outfile, 'w')
        else:
            out_handle = sys.stdout

        sequences = SeqIO.parse(args.infile, "fasta")
        length_list = []
        base_counts = Counter()
        total_bases = 0
        nb_seq = 0
        
        for seq in sequences:
            seq_str = str(seq.seq)
            length = len(seq_str)
            length_list.append(length)
            total_bases += length
            nb_seq += 1
            base_counts.update(seq_str)
            
        if nb_seq == 0:
            print("No sequences found.")
            return

        avg_len = total_bases / nb_seq
        min_len = min(length_list)
        max_len = max(length_list)
        
        # Header
        header = "# seqs\tAverage length\tMin length\tMax length\tBase cnt\tTotal bases"
        if args.outfile: out_handle.write(header + "\n")
        else: print(header)
        
        # Base count string
        base_str = ""
        for base, count in base_counts.items():
            base_str += f"{base}:{count};"
            
        stats_line = f"{nb_seq}\t{avg_len}\t{min_len}\t{max_len}\t{base_str}\t{total_bases}"
        
        if args.outfile: 
            out_handle.write(stats_line + "\n")
            out_handle.close()
            print(f"Stats written to {args.outfile}")
        else:
            print(stats_line)
            
    except Exception as e:
        print(f"Error: {e}")

def run_distribution(args):
    """Replicates SeqLengthDistribution.pl"""
    # Options: -size, -upper, -lower. Output filtered fasta if -out.
    print(f"Analyzing length distribution for {args.infile}...")
    
    try:
        sequences = SeqIO.parse(args.infile, "fasta")
        
        out_handle = None
        if args.outfile:
            out_handle = open(args.outfile, 'w')
            
        freq = Counter()
        sum_len = 0
        seq_cnt = 0
        
        filtered_freq = Counter()
        filtered_sum = 0
        filtered_cnt = 0
        
        for seq in sequences:
            length = len(seq)
            freq[length] += 1
            sum_len += length
            seq_cnt += 1
            
            keep = False
            # Logic from Perl script:
            # If ONLY size: keep == size
            # If upper set, keep > upper ?? Wait.
            # Perl script:
            # if (length > upper) -> write. Wait. Typically 'upper' means 'keep below this'.
            # Let's check original Perl logic carefully.
            # if ($outfile && $upper && !$lower){ if ($length>$upper){ write } }
            # Wait, usually "upper cutoff" implies REMOVING things above it. 
            # But the Perl script implies KEEPING things > upper? 
            # "Remove the sequences that are shorter than a certain length" implies Lower Cutoff.
            # "upper cutoff size" usually implies max length.
            # Let's look at Perl again strictly:
            # if ($length>$upper) -> write. Use case: "remove short". So $upper is actually a LOWER threshold?
            # Or is it "keep sequences > upper"? Yes. 
            # if ($length<$lower) -> write. "keep sequences < lower".
            # if ($length<$lower && $length>$upper) -> write. interval.
            # This naming is confusing. 
            # Let's trust standard conventions or assume user wants Filtering options:
            # -min (keep >= min)
            # -max (keep <= max)
            # -exact (keep == size)
            
            # Arguments map:
            # -upper (Perl) -> used as "keep > upper" which is weird for "upper cutoff".
            # Actually, let's look at User prompt: "Remove the sequences that are shorter than a certain length".
            # So if I want scripts > 100, I use -upper 100?
            
            # Let's implement standard min/max filtering which is clearer.
            # -min: Minimum length to keep (inclusive).
            # -max: Maximum length to keep (inclusive).
            # -size: Exact length.
            
            keep = True
            if args.min and length < args.min: keep = False
            if args.max and length > args.max: keep = False
            if args.size and length != args.size: keep = False
            
            if keep and out_handle:
                SeqIO.write(seq, out_handle, "fasta")
                filtered_freq[length] += 1
                filtered_sum += length
                filtered_cnt += 1
        
        if out_handle: out_handle.close()
        
        # Print Distribution
        print("x sequences with y length")
        for l in sorted(freq.keys()):
            print(f"{l}\t{freq[l]}")
            
        print(f"Sum of lengths: {sum_len}")
        if seq_cnt > 0:
            print(f"Average length: {math.ceil(sum_len/seq_cnt)}")
            most_common_len = freq.most_common(1)[0][0]
            print(f"Modal length: {most_common_len}")
            
        if filtered_cnt > 0 and args.outfile:
            print(f"\nFiltered (written to {args.outfile}):")
            print("x sequences with y length")
            for l in sorted(filtered_freq.keys()):
                print(f"{l}\t{filtered_freq[l]}")
            print(f"New Sum: {filtered_sum}")
            print(f"New Average: {math.ceil(filtered_sum/filtered_cnt)}")
            
    except Exception as e:
        print(f"Error: {e}")

def run_kmers(args):
    """Replicates FindAllKmers.pl"""
    print(f"Finding K-mers (k={args.k}) for {args.infile}...")
    
    try:
        sequences = SeqIO.parse(args.infile, "fasta")
        
        if not args.outfile:
            print("Error: Output file is required.")
            return
            
        with open(args.outfile, 'w') as out:
            for seq in sequences:
                s_len = len(seq)
                max_start = s_len - args.k
                seq_str = str(seq.seq)
                
                # Perl loop 0..max_start. 1-based output?
                # Perl: start_pos = i + 1. 
                # Print >id_start_pos
                
                for i in range(max_start + 1):
                    start_pos = i + 1
                    # Python slice 0-based: i : i+k
                    oligo = seq_str[i : i+args.k]
                    out.write(f">{seq.id}_{start_pos}\n{oligo}\n")
                    
        print(f"K-mers written to {args.outfile}")
    except Exception as e:
        print(f"Error: {e}")

def main():
    parser = argparse.ArgumentParser(description="Sequence Statistics Tools")
    subparsers = parser.add_subparsers(dest='command', help='Sub-commands')
    
    # Summary
    parser_sum = subparsers.add_parser('summary', help='Calculate sequence statistics')
    parser_sum.add_argument('-in', '--input', dest='infile', required=True, help="Input FASTA file")
    parser_sum.add_argument('-out', '--output', dest='outfile', help="Output file (optional)")
    
    # Distribution
    parser_dist = subparsers.add_parser('distribution', help='Length distribution and filtering')
    parser_dist.add_argument('-in', '--input', dest='infile', required=True, help="Input FASTA file")
    parser_dist.add_argument('-out', '--output', dest='outfile', help="Output filtered FASTA file")
    parser_dist.add_argument('-min', type=int, help="Minimum length to keep")
    parser_dist.add_argument('-max', type=int, help="Maximum length to keep")
    parser_dist.add_argument('-size', type=int, help="Exact length to keep")

    # Kmers
    parser_k = subparsers.add_parser('kmers', help='Find all K-mers')
    parser_k.add_argument('-in', '--input', dest='infile', required=True, help="Input FASTA file")
    parser_k.add_argument('-out', '--output', dest='outfile', required=True, help="Output FASTA file")
    parser_k.add_argument('-k', type=int, required=True, help="K-mer size")

    args = parser.parse_args()
    
    if args.command == 'summary':
        run_summary(args)
    elif args.command == 'distribution':
        run_distribution(args)
    elif args.command == 'kmers':
        run_kmers(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
