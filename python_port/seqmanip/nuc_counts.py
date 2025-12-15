import argparse
import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict

def codon_table(code_type):
    # Standard genetic code (NCBI Table 1)
    # The Perl script defines its own hash. Biopython has standard tables.
    # The Perl script allows 'dna' or 'rna'.
    # We can use Biopython's CodonTable.
    from Bio.Data import CodonTable
    if code_type.lower() == 'rna':
        return CodonTable.unambiguous_rna_by_name["Standard"]
    else:
        return CodonTable.unambiguous_dna_by_name["Standard"]

def run_ambig(args):
    """Replicates CountAmbig.pl"""
    print(f"Counting ambiguities in {args.infile}...")
    
    # Definitions from Perl script
    # All chars including N
    all_ambig_chars = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N']
    # exclude N
    strict_ambig_chars = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V']
    
    target_chars = strict_ambig_chars if args.ambigonly else all_ambig_chars
    
    try:
        if args.outfile:
             out_handle = open(args.outfile, 'w')
        else:
             out_handle = sys.stdout # Default to stdout if not implicitly redirectable in original? 
             # Original script prints to stdout. Python port can do same.
             
        sequences = SeqIO.parse(args.infile, "fasta")
        
        for seq in sequences:
            seq_str = str(seq.seq).upper()
            total_count = 0
            details = defaultdict(int)
            
            for char in target_chars:
                # Regex count or simple count?
                # Perl uses: $str =~ /\Q$c/gi
                # String count is faster
                cnt = seq_str.count(char)
                total_count += cnt
                details[char] = cnt
                
            if args.detail:
                # Summary of each char across the sequence?
                # Perl: print "$char\t$details{$char}\n" (This seems to print for EVERY sequence? Wait)
                # Perl loop: foreach sequence -> calculate -> if (detail) print stats ELSE print id sum
                # So if detail is on, it prints stats for EACH sequence? That's verbose. 
                # Let's inspect Perl again.
                # Yes: inside while(next_seq), if (detail) foreach key print.
                # It doesn't print the ID in detail mode? "print "$char\t$details{$char}\n";"
                # That seems odd, multiple seqs would mix. But let's replicate logic.
                if args.outfile: out_handle.write(f"Sequence: {seq.id}\n")
                else: print(f"Sequence: {seq.id}")
                
                for char in sorted(details.keys()):
                    line = f"{char}\t{details[char]}"
                    if args.outfile: out_handle.write(line + "\n")
                    else: print(line)
            else:
                line = f"{seq.id}\t{total_count}"
                if args.outfile: out_handle.write(line + "\n")
                else: print(line)
                
        if args.outfile: out_handle.close()
            
    except FileNotFoundError:
        print(f"Error: Input file '{args.infile}' not found.")
        sys.exit(1)

def run_syn_nonsyn(args):
    """Replicates CountSynNonSyn.pl"""
    print(f"Counting Syn/NonSyn mutations relative to {args.ref_file}...")
    
    # Load Reference
    try:
        ref_seq_obj = next(SeqIO.parse(args.ref_file, "fasta"))
    except StopIteration:
        print("Error: Empty reference file.")
        sys.exit(1)
    except FileNotFoundError:
        print(f"Error: Reference file '{args.ref_file}' not found.")
        sys.exit(1)
        
    ref_seq_str = str(ref_seq_obj.seq).upper()
    
    # Validate length
    if len(ref_seq_str) % 3 != 0:
        print("Warning: Reference sequence length is not a multiple of 3.")

    # Split ref into codons
    ref_codons = [ref_seq_str[i:i+3] for i in range(0, len(ref_seq_str), 3)]
    
    # Get Genetic Code
    # Biopython table
    table = codon_table(args.code)
    # Helper to translate
    def translate_codon(codon):
        # Handle gaps or partials
        if len(codon) < 3 or '-' in codon or 'N' in codon:
            return 'X'
        try:
            return table.forward_table.get(codon, table.stop_codons.get(codon, '*')) # Check stop codons manually if not in forward
            # Actually forward_table only has amino acids. stop_codons list has stops.
            if codon in table.stop_codons: return '*'
            if codon in table.forward_table: return table.forward_table[codon]
            return 'X' 
        except:
            return 'X'

    # Stats containers
    syn_dist = defaultdict(int)     # Count of seqs with X syn mutations
    nonsyn_dist = defaultdict(int)  # Count of seqs with X nonsyn mutations
    nuc_dist = defaultdict(int)     # Count of seqs with X nuc mutations
    
    try:
        sequences = SeqIO.parse(args.infile, "fasta")
        out_handle = open(args.outfile, 'w') if args.outfile else sys.stdout 
        # Actually original script requires -out argument.
        
        for seq in sequences:
            seq_str = str(seq.seq).upper()
            seq_codons = [seq_str[i:i+3] for i in range(0, len(seq_str), 3)]
            
            nonsyn_count = 0
            syn_count = 0
            nuc_count = 0
            
            # Iterate up to min length (ignore extra length?)
            limit = min(len(ref_codons), len(seq_codons))
            
            for i in range(limit):
                ref_c = ref_codons[i]
                seq_c = seq_codons[i]
                
                # Check nucleotide diffs
                codon_diffs = 0
                for j in range(len(ref_c)):
                    if j < len(seq_c) and ref_c[j] != seq_c[j]:
                        codon_diffs += 1
                        nuc_count += 1
                
                if codon_diffs > 0:
                    ref_aa = translate_codon(ref_c)
                    seq_aa = translate_codon(seq_c)
                    
                    if ref_aa != seq_aa:
                        nonsyn_count += codon_diffs # Original script adds 'codoncount' (number of nuc mutations in this codon)
                    else:
                        syn_count += codon_diffs
            
            # Store stats
            syn_dist[syn_count] += 1
            nonsyn_dist[nonsyn_count] += 1
            nuc_dist[nuc_count] += 1
            
            if nonsyn_count > 4:
                # Original output: print "$seq_id has $nonsyncount nonsynonymous mutations\n";
                print(f"{seq.id} has {nonsyn_count} nonsynonymous mutations")

        # Output Summary
        out_handle.write("Number of sequences with x synonymous mutations\n")
        for k in sorted(syn_dist.keys()):
            out_handle.write(f"{k} {syn_dist[k]}\n")
            
        out_handle.write("Number of sequences with x nonsynonymous mutations\n")
        for k in sorted(nonsyn_dist.keys()):
            out_handle.write(f"{k} {nonsyn_dist[k]}\n")
            
        out_handle.write("Number of sequences with x nucleotide mutations\n")
        for k in sorted(nuc_dist.keys()):
            out_handle.write(f"{k} {nuc_dist[k]}\n")
            
        if args.outfile: out_handle.close()
        
    except FileNotFoundError:
        print(f"Error: Input file '{args.infile}' not found.")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Nucleotide Sequence counts (Ambig & Syn/NonSyn)")
    subparsers = parser.add_subparsers(dest='command', help='Sub-commands')
    
    # Ambig Subcommand
    parser_ambig = subparsers.add_parser('ambig', help='Count ambiguous characters')
    parser_ambig.add_argument('-in', '--input', dest='infile', required=True, help="Input FASTA file")
    parser_ambig.add_argument('-out', '--output', dest='outfile', help="Output file (optional, defaults to stdout)")
    parser_ambig.add_argument('-detail', action='store_true', help="Provide summary of each ambiguous IUPAC char")
    parser_ambig.add_argument('-ambigonly', action='store_true', help="Count only ambiguous sites (exclude N)")
    
    # Syn-NonSyn Subcommand
    parser_sns = subparsers.add_parser('syn-nonsyn', help='Count Synonymous/Non-synonymous mutations')
    parser_sns.add_argument('-in', '--input', dest='infile', required=True, help="Aligned FASTA file")
    parser_sns.add_argument('-out', '--output', dest='outfile', required=True, help="Output data file")
    parser_sns.add_argument('-ref', '--reference', dest='ref_file', required=True, help="Reference sequence file")
    parser_sns.add_argument('-code', dest='code', choices=['dna', 'rna'], default='dna', help="Genetic code (dna/rna)")

    args = parser.parse_args()
    
    if args.command == 'ambig':
        run_ambig(args)
    elif args.command == 'syn-nonsyn':
        run_syn_nonsyn(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
