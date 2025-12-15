import argparse
import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable

def run_standard(args):
    """Replicates translate_seq.pl and translate_seq_noIUPAC.pl"""
    print(f"Translating {args.infile} to {args.outfile}...")
    
    try:
        if args.outfile:
            out_handle = open(args.outfile, 'w')
        else:
            out_handle = sys.stdout

        sequences = SeqIO.parse(args.infile, "fasta")
        count = 0
        skipped = 0
        
        # Regex for strict mode (allow A,C,G,T,U,N case insensitive)
        strict_pattern = re.compile(r'^[ACGTUN]+$', re.IGNORECASE)
        
        for seq in sequences:
            seq_str = str(seq.seq)
            
            if args.strict:
                if not strict_pattern.match(seq_str):
                    print(f"Skipping {seq.id}: contains non-standard characters")
                    skipped += 1
                    continue
            
            # Translate
            try:
                # Biopython translate. 
                # Original script translate_seq.pl uses: $seq->translate(-unknown => $unknown);
                # BioPerl translate handles gaps by making them X usually? Or preserves them?
                # translate_seq.pl seems to be for standard sequences.
                # translate_seq_noIUPAC.pl implies strictness.
                
                # Check if gap='-' needed?
                # If standard tool, maybe just use default behavior (gaps -> X or Error).
                # But BioPerl often flexible.
                # Let's use 'X' for unknown/gap in standard mode unless it's a gap?
                protein_seq = seq.seq.translate(table="Standard", to_stop=False, cds=False, gap=None) 
                
                # Wait, gap=None means gaps raise error?
                # If we want to mimic Perl: unknown => 'X'.
                # Can simple .translate method handle unknown?
                # Biopython: "Ambiguous codons... 'X'".
                
                # Construct new record
                new_rec = seq
                new_rec.seq = protein_seq
                SeqIO.write(new_rec, out_handle, "fasta")
                count += 1
                
            except Exception as e:
                print(f"Error translating {seq.id}: {e}")
                skipped += 1
                
        if args.outfile: out_handle.close()
        print(f"Translated {count} sequences. Skipped {skipped}.")

    except Exception as e:
        print(f"Error: {e}")

def run_gapped(args):
    """Replicates TranslateWithGaps.pl"""
    print(f"Translating gapped alignment {args.infile} (Code: {args.gencode})...")
    
    # Custom mapping for simple universal gapped
    # '---' -> '-'
    # Biopython translate(gap='-') handles this if gaps are multiple of 3.
    
    try:
        sequences = SeqIO.parse(args.infile, "fasta")
        
        if args.outfile:
            out_handle = open(args.outfile, 'w')
        else:
            out_handle = sys.stdout
            
        count = 0
        
        for seq in sequences:
            # Check length % 3
            if len(seq) % 3 != 0:
                print(f"Warning: {seq.id} length {len(seq)} is not divisible by 3. Result may be truncated/aligned poorly.")
            
            # Use Biopython's gap handling
            # If gap='-', it preserves '---' as '-'
            # It puts 'X' for 'N--' or other frameshifts?
            try:
                p_seq = seq.seq.translate(table="Standard", gap='-')
                seq.seq = p_seq
                SeqIO.write(seq, out_handle, "fasta")
                count += 1
            except Exception as e:
                print(f"Error translating {seq.id}: {e}")
        
        if args.outfile: out_handle.close()
        print(f"Translated {count} sequences.")
        
    except Exception as e:
        print(f"Error: {e}")

def main():
    parser = argparse.ArgumentParser(description="Translation Tools")
    subparsers = parser.add_subparsers(dest='command', help='Sub-commands')
    
    # Standard
    parser_std = subparsers.add_parser('standard', help='Standard translation')
    parser_std.add_argument('-in', '--input', dest='infile', required=True, help="Input FASTA file")
    parser_std.add_argument('-out', '--output', dest='outfile', help="Output file")
    parser_std.add_argument('-u', '--unknown', dest='unknown', default='X', help="Character for unknown AA [default: X]")
    parser_std.add_argument('-strict', action='store_true', help="Skip sequences with non-standard characters (replicates translate_seq_noIUPAC)")

    # Gapped
    parser_gap = subparsers.add_parser('gapped', help='Translate alignment with gaps')
    parser_gap.add_argument('-in', '--input', dest='infile', required=True, help="Input FASTA file")
    parser_gap.add_argument('-out', '--output', dest='outfile', help="Output file")
    parser_gap.add_argument('-gc', dest='gencode', default='universal', help="Genetic Code (universal/flu/degen) [default: universal]")

    args = parser.parse_args()
    
    if args.command == 'standard':
        run_standard(args)
    elif args.command == 'gapped':
        run_gapped(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
