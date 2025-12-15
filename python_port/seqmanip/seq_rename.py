import argparse
import sys
import os
from Bio import SeqIO

def run_table_rename(args):
    """Replicates RenameFastaWithTable.pl"""
    print(f"Renaming IDs in {args.infile} using table {args.table}...")
    
    # Parse columns (1-based to 0-based)
    col_idx = args.col - 1
    newcol_idx = args.newcol - 1
    
    # Load Table
    id_map = {}
    try:
        with open(args.table, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) > max(col_idx, newcol_idx):
                    old_id = parts[col_idx]
                    new_id = parts[newcol_idx]
                    # Common cleanup from original script
                    new_id = new_id.replace(' ', '_')
                    id_map[old_id] = new_id
    except FileNotFoundError:
        print(f"Error: Table file '{args.table}' not found.")
        sys.exit(1)
        
    # Process Sequences
    try:
        if args.outfile:
             out_handle = open(args.outfile, 'w')
        else:
             out_handle = sys.stdout

        sequences = SeqIO.parse(args.infile, "fasta")
        count = 0
        renamed = 0
        
        for seq in sequences:
            original_id = seq.id
            # Original script had specific cleaning: s/cds://g; s/\/\d+-\d+//g;
            # We will try exact match first.
            
            if original_id in id_map:
                seq.id = id_map[original_id]
                seq.description = seq.description.replace(original_id, seq.id) # Update desc if it contains ID
                renamed += 1
            else:
                 # Check if cleanup helps matching?
                 # For now, stick to exact match to be safe.
                 pass
                 
            SeqIO.write(seq, out_handle, "fasta")
            count += 1
            
        if args.outfile: out_handle.close()
        print(f"Processed {count} sequences. Renamed {renamed}.")
        
    except FileNotFoundError:
        print(f"Error: Input file '{args.infile}' not found.")
        sys.exit(1)

def run_prefix_rename(args):
    """Replicates SeqRenamer.pl (prefix mode)"""
    print(f"Renaming sequences in {args.infile} with prefix '{args.prefix}'...")
    
    try:
        sequences = SeqIO.parse(args.infile, args.fmt)
        
        if args.outfile:
             out_handle = open(args.outfile, 'w')
        else:
             out_handle = sys.stdout
             
        # Map file
        map_filename = f"{args.prefix}_nametable.txt"
        map_handle = open(map_filename, 'w')
        
        count = 0
        for seq in sequences:
            old_id = seq.id
            new_id = f"{args.prefix}{count}" # Original: $pre$cnt (0-based) or 1-based?
            # Perl: $pre$cnt then $cnt++. So 0-based.
            
            seq.id = new_id
            seq.description = f"{new_id} {old_id}" # Original appends old ID to desc
            
            SeqIO.write(seq, out_handle, args.fmt)
            map_handle.write(f"{old_id}\t{new_id}\n")
            count += 1
            
        if args.outfile: out_handle.close()
        map_handle.close()
        print(f"Renamed {count} sequences. Mapping saved to {map_filename}")

    except FileNotFoundError:
        print(f"Error: Input file '{args.infile}' not found.")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Sequence renaming tools")
    subparsers = parser.add_subparsers(dest='command', help='Sub-commands')
    
    # Table Rename
    parser_table = subparsers.add_parser('table', help='Rename IDs using a tab-delimited table')
    parser_table.add_argument('-in', '--input', dest='infile', required=True, help="Input FASTA file")
    parser_table.add_argument('-out', '--output', dest='outfile', help="Output file")
    parser_table.add_argument('-table', dest='table', required=True, help="Tab-delimited table")
    parser_table.add_argument('-col', dest='col', type=int, default=1, help="Column with old name (1-based) [default: 1]")
    parser_table.add_argument('-newcol', dest='newcol', type=int, default=2, help="Column with new name (1-based) [default: 2]")
    
    # Prefix Rename
    parser_pre = subparsers.add_parser('prefix', help='Rename with prefix and counter')
    parser_pre.add_argument('-in', '--input', dest='infile', required=True, help="Input file")
    parser_pre.add_argument('-out', '--output', dest='outfile', help="Output file")
    parser_pre.add_argument('-pre', '--prefix', dest='prefix', required=True, help="Prefix string")
    parser_pre.add_argument('-fmt', dest='fmt', default='fasta', help="Format (fasta/fastq) [default: fasta]")

    args = parser.parse_args()
    
    if args.command == 'table':
        run_table_rename(args)
    elif args.command == 'prefix':
        run_prefix_rename(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
