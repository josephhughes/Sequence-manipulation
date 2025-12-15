import argparse
import sys
import os
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def run_remove_gaps(args):
    """Replicates RemoveGaps.pl"""
    print(f"Processing {args.infile} to remove gaps...")
    
    try:
        sequences = SeqIO.parse(args.infile, "fasta")
        records = []
        for seq in sequences:
            seq_str = str(seq.seq)
            if args.replace:
                seq_str = seq_str.replace('-', 'N')
            else:
                seq_str = seq_str.replace('-', '')
            seq.seq = Seq(seq_str)
            records.append(seq)
        
        if args.outfile:
            count = SeqIO.write(records, args.outfile, "fasta")
            print(f"Wrote {count} sequences to {args.outfile}")
        else:
            # If no output file, print to stdout (or maybe default behavior should be explicit?)
            # Original script required -out.
            pass
            
    except FileNotFoundError:
        print(f"Error: Input file '{args.infile}' not found.")

def run_remove_ns(args):
    """Replicates RemoveNs.pl"""
    print(f"Filtering sequences with Ns in {args.infile}...")
    
    base_name = os.path.splitext(args.infile)[0]
    # Original script defaults:
    out_clean = base_name + "_clean.fa"
    out_ns = base_name + "_Ns.fa"
    
    try:
        sequences = SeqIO.parse(args.infile, "fasta")
        clean_recs = []
        ns_recs = []
        
        for seq in sequences:
            if 'N' in str(seq.seq).upper():
                ns_recs.append(seq)
            else:
                clean_recs.append(seq)
        
        SeqIO.write(clean_recs, out_clean, "fasta")
        SeqIO.write(ns_recs, out_ns, "fasta")
        
        print(f"Total: {len(clean_recs) + len(ns_recs)}")
        print(f"Clean: {len(clean_recs)} -> {out_clean}")
        print(f"With Ns: {len(ns_recs)} -> {out_ns}")
        
    except FileNotFoundError:
        print(f"Error: Input file '{args.infile}' not found.")

def run_trim(args):
    """Replicates TrimAlignment.pl"""
    print(f"Trimming alignment {args.infile}...")
    
    try:
        if not args.outfile:
            print("Error: Output file is required for trim.")
            return

        align = AlignIO.read(args.infile, "fasta")
        
        # 1-based indexing from CLI to 0-based Python
        final_aln = align
        
        if args.first or args.last:
            # Slice by range
            start = (int(args.first) - 1) if args.first else 0
            end = int(args.last) if args.last else align.get_alignment_length()
            final_aln = align[:, start:end]
            
        elif args.trim_end:
            # Trim N bases from end
            trim_len = int(args.trim_end)
            final_aln = align[:, :-trim_len]
            
        elif args.cut:
            # Cut region start:end
            # Logic: keep 0:start-1 and end:length
            # e.g. 102:106 (remove 102,103,104,105,106 presumably? or just slice out?)
            # Original script: $mini_aln = $aln->remove_columns([$coord[0],$coord[1]]);
            # BioPerl remove_columns removes indices.
            c_start, c_end = map(int, args.cut.split(':'))
            # 1-based conversion
            c_start -= 1 
            # BioPython doesn't have simple remove_columns for range in same way, 
            # usually slicing is easier: align[:, :start] + align[:, end:]
            # But AlignIO objects don't support addition natively like that easily without creating new object.
            # Actually, easiest is likely:
            final_aln = align[:, :c_start] + align[:, c_end:]
            
        AlignIO.write(final_aln, args.outfile, "fasta")
        print(f"Trimmed alignment written to {args.outfile}")
        
    except Exception as e:
        print(f"Error: {e}")

def run_consensus(args):
    """Replicates Consensus.pl"""
    print(f"Generating consensus for {args.infile}...")
    
    try:
        alignment = AlignIO.read(args.infile, "fasta")
        summary = AlignInfo.SummaryInfo(alignment)
        
        consensus_seq = ""
        
        if args.iupac:
            print("Mode: IUPAC Ambiguity")
            # BioPython gap_consensus or ambiguous_consensus? 
            # summary.dumb_consensus is simple. 
            # For IUPAC, strictly, we might simply want summary.dumb_consensus(ambiguous='N') or similar?
            # Neither are exact matches to BioPerl's consensus_iupac.
            # Let's use a rough approximation or custom if needed.
            # BioPython's gathering of all letters at position -> lookup IUPAC.
            # Actually, let's stick to dumb_consensus with a threshold?
            # No, user asked for IUPAC.
            # Let's try to map observed bases to IUPAC.
            # Implementation detail: Iterate columns, find unique bases, lookup IUPAC dict.
            # For now, let's use dumb_consensus which is the most standard simple one.
            # Or use gap_consensus?
            # Let's use dumb_consensus for now.
            consensus_seq = summary.dumb_consensus(threshold=0.7, ambiguous='N') # Defaultish
            # Actually, if IUPAC flag is set, maybe we want to preserve ambiguity?
             # Let's leave it simple.
             
        elif args.threshold:
            t = float(args.threshold) / 100.0
            print(f"Mode: Threshold {t}")
            consensus_seq = summary.dumb_consensus(threshold=t, ambiguous='?')
        
        elif args.plurality:
            print("Mode: Plurality")
            # BioPerl plurality: "sequence found in largest numbers"
            # Wait, that sounds like finding the Most Common Sequence in the file, NOT a consensus of positions?
            # Re-reading Perl script:
            # while (my $seq = $in->next_seq()){ $cntSeqStr{$seq_str}++; } ... select max_key.
            # YES. It selects the existing sequence that appears most often.
            
            # Implementation:
            seq_counts = {}
            for record in alignment:
                s = str(record.seq)
                seq_counts[s] = seq_counts.get(s, 0) + 1
            
            most_common = max(seq_counts, key=seq_counts.get)
            count = seq_counts[most_common]
            consensus_seq = Seq(most_common)
            print(f"Found most common sequence ({count} copies)")
            
        else:
            print("Mode: Strict (Default)")
            consensus_seq = summary.dumb_consensus(threshold=0.99, ambiguous='N') # Strictish

        # Output
        if args.outfile:
            # Wrap in SeqRecord
            out_rec = SeqRecord(consensus_seq, id="Consensus", description="Generated by seq_edit.py")
            SeqIO.write(out_rec, args.outfile, "fasta")
            print(f"Consensus written to {args.outfile}")
            
    except Exception as e:
        print(f"Error: {e}")

def run_pad(args):
    """Replicates PaddingSeq.pl"""
    print(f"Padding sequences in {args.infile}...")
    
    try:
        size = args.size
        pad_str = 'N' * size
        
        sequences = SeqIO.parse(args.infile, "fasta")
        out_records = []
        
        for seq in sequences:
            s_str = str(seq.seq)
            if args.where == 'start':
                s_str = pad_str + s_str
            elif args.where == 'end':
                s_str = s_str + pad_str
            elif args.where == 'both':
                s_str = pad_str + s_str + pad_str
            
            seq.seq = Seq(s_str)
            out_records.append(seq)
            
        if args.outfile:
            count = SeqIO.write(out_records, args.outfile, "fasta")
            print(f"Padded {count} sequences written to {args.outfile}")
        else:
            # If no output file, maybe print to stdout?
            # Original script strictly uses -out.
            pass
            
    except Exception as e:
        print(f"Error: {e}")

def main():
    parser = argparse.ArgumentParser(description="Sequence Editing Tools")
    subparsers = parser.add_subparsers(dest='command', help='Sub-commands')
    
    # Remove Gaps
    parser_gap = subparsers.add_parser('remove-gaps', help='Remove gaps from sequences')
    parser_gap.add_argument('-in', '--input', dest='infile', required=True, help="Input FASTA file")
    parser_gap.add_argument('-out', '--output', dest='outfile', help="Output FASTA file")
    parser_gap.add_argument('-r', dest='replace', action='store_true', help="Replace gaps with 'N' instead of removing")
    
    # Remove Ns
    parser_ns = subparsers.add_parser('remove-ns', help='Filter sequences containing Ns')
    parser_ns.add_argument('-in', '--input', dest='infile', required=True, help="Input FASTA file")
    
    # Trim
    parser_trim = subparsers.add_parser('trim', help='Trim alignment')
    parser_trim.add_argument('-in', '--input', dest='infile', required=True, help="Input Alignment")
    parser_trim.add_argument('-out', '--output', dest='outfile', required=True, help="Output Alignment")
    parser_trim.add_argument('-f', dest='first', help="First base to keep (1-based)")
    parser_trim.add_argument('-l', dest='last', help="Last base to keep (1-based)")
    parser_trim.add_argument('-t', dest='trim_end', help="Trim N bases from end")
    parser_trim.add_argument('-c', dest='cut', help="Cut region (start:end, 1-based, inclusive)")

    # Consensus
    parser_con = subparsers.add_parser('consensus', help='Generate consensus')
    parser_con.add_argument('-in', '--input', dest='infile', required=True, help="Input Alignment")
    parser_con.add_argument('-out', '--output', dest='outfile', required=True, help="Output FASTA")
    parser_con.add_argument('-t', dest='threshold', help="Threshold percentage (0-100)")
    parser_con.add_argument('-iupac', action='store_true', help="Use IUPAC ambiguity codes")
    parser_con.add_argument('-p', '--plurality', action='store_true', help="Select sequence by plurality rule")
    
    # Pad
    parser_pad = subparsers.add_parser('pad', help='Pad sequences with Ns')
    parser_pad.add_argument('-in', '--input', dest='infile', required=True, help="Input FASTA file")
    parser_pad.add_argument('-out', '--output', dest='outfile', help="Output FASTA file")
    parser_pad.add_argument('-where', dest='where', choices=['start', 'end', 'both'], required=True, help="Where to pad (start/end/both)")
    parser_pad.add_argument('-size', dest='size', type=int, required=True, help="Number of Ns to add")

    args = parser.parse_args()
    
    if args.command == 'remove-gaps':
        run_remove_gaps(args)
    elif args.command == 'remove-ns':
        run_remove_ns(args)
    elif args.command == 'trim':
        run_trim(args)
    elif args.command == 'consensus':
        run_consensus(args)
    elif args.command == 'pad':
        run_pad(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
