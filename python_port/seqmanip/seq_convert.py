import argparse
import sys
import os
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator

def run_phylip(args):
    """Fasta <-> Phylip conversion"""
    fmt_in = 'phylip' if args.reverse else 'fasta'
    fmt_out = 'fasta' if args.reverse else 'phylip'
    
    print(f"Converting {args.infile} ({fmt_in}) to {args.outfile} ({fmt_out})...")
    
    try:
        # Biopython's 'phylip' is strict (10 char names). 'phylip-relaxed' allows longer.
        # Original Perl script padded to 10 chars, implying strict.
        # However, for modern usage, 'phylip-relaxed' is often preferred if names are long.
        # Let's try Strict first to match original intent, but warn or fallback if issues?
        # Actually, using 'phylip-relaxed' is safer for general usage unless strict is required.
        # But if the user wants to mimic the Perl script exactly...
        # The Perl script explicitly calculated padding for 10 chars. 
        # So we should probably use 'phylip' (strict).
        
        # NOTE: Biopython might truncate long names silently or throw error in strict mode.
        count = SeqIO.convert(args.infile, fmt_in, args.outfile, fmt_out)
        print(f"Converted {count} sequences.")
        
    except Exception as e:
        print(f"Error converting: {e}")
        # If strict failed, maybe suggest relaxed?
        if fmt_out == 'phylip':
             print("Tip: If names are too long for strict Phylip (10 chars), Biopython enforces truncation/error.")

def run_gb2fa(args):
    """GenBank -> Fasta"""
    print(f"Converting {args.infile} (GenBank) to {args.outfile} (FASTA)...")
    
    try:
        if args.gi:
            # Custom ID formatting matching original script: gi|primary_id|gb|accession
            records = []
            for rec in SeqIO.parse(args.infile, "genbank"):
                # Clean/extract IDs
                # Biopython rec.id is usually Accession.Version
                # rec.name is ACCESSION
                # rec.annotations['gi'] might exist if old style GB file
                
                gi = rec.annotations.get('gi', 'NA')
                acc = rec.annotations.get('accessions', ['NA'])[0]
                
                # Original: my $gi = $inseq->primary_id(); my $acc = $inseq->accession_number();
                # Biopython mappings vary by parser version/file type.
                # Let's try to construct a compatible ID.
                new_id = f"gi|{gi}|gb|{acc}"
                rec.id = new_id
                rec.description = rec.description # Keep description
                records.append(rec)
            
            count = SeqIO.write(records, args.outfile, "fasta")
        else:
            # Standard conversion
            count = SeqIO.convert(args.infile, "genbank", args.outfile, "fasta")
            
        print(f"Converted {count} sequences.")
        
    except Exception as e:
        print(f"Error converting: {e}")

def run_fastq_qual(args):
    """Fastq <-> Fasta+Qual"""
    
    base_name = os.path.splitext(args.infile)[0]
    
    if args.mode == 'to_fastq':
        # FASTA + QUAL -> FASTQ
        # Input: args.infile (fasta). Assumes .qual file exists with same base name OR specified?
        # Original script: takes fasta, looks for .qual.
        
        qual_file = base_name + ".qual"
        if not os.path.exists(qual_file):
             # Try replacing extension
             if args.infile.endswith(".fasta") or args.infile.endswith(".fna") or args.infile.endswith(".fa"):
                 base = args.infile.rsplit('.', 1)[0]
                 qual_file = base + ".qual"
        
        print(f"Merging {args.infile} and {qual_file} into {args.outfile}...")
        
        try:
            with open(args.infile) as fa_handle, open(qual_file) as qual_handle, open(args.outfile, "w") as out_handle:
                records = PairedFastaQualIterator(fa_handle, qual_handle)
                count = SeqIO.write(records, out_handle, "fastq")
            print(f"Converted {count} sequences.")
        except FileNotFoundError as e:
            print(f"Error: {e}")
            sys.exit(1)
        except Exception as e:
            print(f"Error: {e}")
            sys.exit(1)

    elif args.mode == 'from_fastq':
        # FASTQ -> FASTA + QUAL
        print(f"Splitting {args.infile} into FASTA + QUAL...")
        
        out_fasta = args.outfile if args.outfile else base_name + ".fasta"
        out_qual = base_name + ".qual"
        if args.outfile:
             # If user provided output name, use it for fasta, append .qual for qualitative?
             # Or assume user provided base? 
             # Let's stick to: Output is the FASTA output. Qual is derived.
             base_out = os.path.splitext(args.outfile)[0]
             out_qual = base_out + ".qual"
             
        try:
            # Biopython convert doesn't split automatically to Qual.
            # Convert to Fasta is easy.
            # Convert to Qual requires 'qual' format writer.
            
            SeqIO.convert(args.infile, "fastq", out_fasta, "fasta")
            SeqIO.convert(args.infile, "fastq", out_qual, "qual")
            
            print(f"Created {out_fasta} and {out_qual}.")
            
        except Exception as e:
            print(f"Error: {e}")

def run_nexus(args):
    """Fasta -> Nexus"""
    # Original used Bio::AlignIO, implying it expects an alignment.
    # Nexus format generally requires alignment (same length seqs).
    print(f"Converting {args.infile} (FASTA) to {args.outfile} (NEXUS)...")
    try:
        # Check alignment or just convert?
        # AlignIO.convert is robust.
        count = AlignIO.convert(args.infile, "fasta", args.outfile, "nexus", molecule_type=args.type)
        print(f"Converted {count} alignments.")
    except Exception as e:
        print(f"Error converting: {e}")
        print("Note: Nexus conversion usually requires all sequences to be the same length (an alignment).")

def run_split(args):
    """Replicates SplitFasta.pl"""
    print(f"Splitting {args.infile} based on metadata {args.meta}...")
    
    try:
        # Load Metadata
        group_map = {}
        col_idx = -1
        
        with open(args.meta, 'r') as f:
            header = f.readline().strip().split('\t')
            try:
                col_idx = header.index(args.colname)
            except ValueError:
                print(f"Error: Column '{args.colname}' not found in metadata header.")
                return

            for line in f:
                parts = line.strip().split('\t')
                if len(parts) > col_idx:
                    seq_id = parts[0] # Assumes ID is first col
                    group = parts[col_idx]
                    group_map[seq_id] = group
                    
        # Process Fasta
        sequences = SeqIO.parse(args.infile, "fasta")
        groups = {}
        
        for seq in sequences:
            if seq.id in group_map:
                g = group_map[seq.id]
                if g not in groups: groups[g] = []
                groups[g].append(seq)
        
        # Write files
        outval = args.outdir if args.outdir else "."
        if not os.path.exists(outval): os.makedirs(outval)
        
        for g, seqs in groups.items():
            # Clean group name for filename
            fname = "".join(x for x in g if x.isalnum() or x in "._-")
            outpath = os.path.join(outval, f"{fname}.fasta")
            count = SeqIO.write(seqs, outpath, "fasta")
            print(f"Wrote {count} sequences to {outpath}")
            
    except Exception as e:
        print(f"Error: {e}")

def run_oneline(args):
    """Replicates fa2oneline.pl"""
    print(f"Linearizing {args.infile} to {args.outfile}...")
    try:
        # 'fasta-2line' format available in newer Biopython (1.71+)
        # If not, we can fall back or user manual write.
        # Check by trying?
        # Actually simplest is to just use it, catch error.
        try:
            count = SeqIO.convert(args.infile, "fasta", args.outfile, "fasta-2line")
            print(f"Converted {count} sequences.")
        except Exception:
            # Fallback for older biopython or just manual
            records = SeqIO.parse(args.infile, "fasta")
            with open(args.outfile, "w") as f:
                c = 0
                for rec in records:
                    f.write(f">{rec.description}\n{str(rec.seq)}\n")
                    c += 1
            print(f"Converted {c} sequences (manual).")
            
    except Exception as e:
        print(f"Error: {e}")

def main():
    parser = argparse.ArgumentParser(description="Sequence Format Conversion Tools")
    subparsers = parser.add_subparsers(dest='command', help='Sub-commands')
    
    # Phylip
    parser_phy = subparsers.add_parser('phylip', help='Convert between FASTA and PHYLIP')
    parser_phy.add_argument('-in', '--input', dest='infile', required=True, help="Input file")
    parser_phy.add_argument('-out', '--output', dest='outfile', required=True, help="Output file")
    parser_phy.add_argument('-reverse', action='store_true', help="Convert PHYLIP to FASTA (default is FASTA to PHYLIP)")
    
    # GenBank
    parser_gb = subparsers.add_parser('gb2fa', help='Convert GenBank to FASTA')
    parser_gb.add_argument('-in', '--input', dest='infile', required=True, help="Input GenBank file")
    parser_gb.add_argument('-out', '--output', dest='outfile', required=True, help="Output FASTA file")
    parser_gb.add_argument('-gi', action='store_true', help="Use custom GI|Accession ID format")
    
    # Fastq/Qual
    parser_fq = subparsers.add_parser('fastq-qual', help='Convert between FASTQ and FASTA+QUAL')
    parser_fq.add_argument('-in', '--input', dest='infile', required=True, help="Input file")
    parser_fq.add_argument('-out', '--output', dest='outfile', help="Output file (FASTA or FASTQ depending on mode)")
    parser_fq.add_argument('-mode', dest='mode', choices=['to_fastq', 'from_fastq'], required=True, help="Conversion direction")

    # Nexus
    parser_nex = subparsers.add_parser('nexus', help='Convert FASTA to NEXUS')
    parser_nex.add_argument('-in', '--input', dest='infile', required=True, help="Input FASTA file")
    parser_nex.add_argument('-out', '--output', dest='outfile', required=True, help="Output NEXUS file")
    parser_nex.add_argument('-type', dest='type', default='DNA', choices=['DNA', 'RNA', 'protein'], help="Molecule type (DNA/RNA/protein) [default: DNA]")

    # Split
    parser_split = subparsers.add_parser('split', help='Split FASTA by metadata group')
    parser_split.add_argument('-in', '--input', dest='infile', required=True, help="Input FASTA file")
    parser_split.add_argument('-meta', dest='meta', required=True, help="Metadata table (tab-delim)")
    parser_split.add_argument('-col', dest='colname', required=True, help="Column name to group by")
    parser_split.add_argument('-outdir', dest='outdir', help="Output directory [default: current]")

    # Oneline
    parser_one = subparsers.add_parser('oneline', help='Convert FASTA to single-line sequence')
    parser_one.add_argument('-in', '--input', dest='infile', required=True, help="Input FASTA file")
    parser_one.add_argument('-out', '--output', dest='outfile', required=True, help="Output FASTA file")

    args = parser.parse_args()
    
    if args.command == 'phylip':
        run_phylip(args)
    elif args.command == 'gb2fa':
        run_gb2fa(args)
    elif args.command == 'fastq-qual':
        run_fastq_qual(args)
    elif args.command == 'nexus':
        run_nexus(args)
    elif args.command == 'split':
        run_split(args)
    elif args.command == 'oneline':
        run_oneline(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
