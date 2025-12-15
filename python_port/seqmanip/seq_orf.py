import argparse
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def find_longest_orf(seq_obj, table=1, min_len=0, strict=True):
    """
    Finds the longest ORF in a sequence.
    strict=True: ORF must start with ATG.
    strict=False: ORF can start at beginning of frame.
    Returns: (start, end, strand, length, seq_str, protein_str)
    """
    best_orf = None
    max_len = -1
    
    seq_len = len(seq_obj)
    
    # 3 Forward frames, 3 Reverse frames
    directions = [(1, seq_obj.seq), (-1, seq_obj.seq.reverse_complement())]
    
    for strand, nuc_seq in directions:
        for frame in range(3):
            # Translate frame
            # We can use Biopython's translate, splitting by stop codons.
            # But we need coordinates. 
            # Logic:
            # 1. Translate entire frame
            # 2. Find longest stretch between stops (or start-stop if strict)
            
            trans = nuc_seq[frame:].translate(table=table)
            prot_str = str(trans)
            
            # This splits by stop (*), doesn't handle ATG start strictly easily without iteration.
            # Let's iterate manually or split.
            
            # Simple approach: Split by *. Iterate segments. 
            # If strict, find first M in segment.
            
            parts = prot_str.split('*')
            
            # Track position in protein string
            current_aa_idx = 0
            
            # Identify full protein sequence length relative to frame start
            
            for i, part in enumerate(parts):
                # part is AA sequence (excluding stop)
                # If strict, we need an M.
                start_aa_rel = 0 # Relative to part start
                
                if strict:
                    m_index = part.find('M')
                    if m_index == -1:
                        # No Methionine, no valid ORF in this segment if strict
                        current_aa_idx += len(part) + 1 # +1 for stop
                        continue
                    start_aa_rel = m_index
                
                # ORF content
                orf_prot = part[start_aa_rel:]
                orf_len_aa = len(orf_prot)
                orf_len_nt = orf_len_aa * 3
                
                if orf_len_nt > max_len:
                    # Calculate coordinates
                    # Start index in nuc_seq (0-based)
                    # frame offset + (current_aa_idx + start_aa_rel) * 3
                    
                    start_nt_local = frame + (current_aa_idx + start_aa_rel) * 3
                    end_nt_local = start_nt_local + orf_len_nt
                    
                    # Coordinates in original sequence?
                    # If strand 1: start_nt_local is absolute (0-based)
                    # If strand -1: seq was rev_comp. 
                    # start_nt_local is from 5' of rev_comp.
                    # Original index: len - end_nt_local?
                    
                    if strand == 1:
                        final_start = start_nt_local
                        final_end = end_nt_local # Python slice exclusive
                    else:
                        # RevComp logic fix
                        # nuc_seq is rev_comp. indices start_nt_local, end_nt_local are on rev_comp.
                        # Original sequence index:
                        # 5' RC ------ 3' RC
                        # 3' Fwd ------ 5' Fwd
                        # Index i in RC corresponds to index (len - 1 - i) in Forward?
                        # No.
                        # RC = Fwd[::-1].complement()
                        # RC[0] is Compl(Fwd[len-1])
                        # RC[i] is Compl(Fwd[len-1-i])
                        # The segment RC[start:end] corresponds to Fwd[len-end : len-start] (reverse complemented)
                        # So original coords are len-end_nt_local to len-start_nt_local.
                        
                        final_start = seq_len - end_nt_local
                        final_end = seq_len - start_nt_local
                        
                    max_len = orf_len_nt
                    
                    # Get DNA sequence logic fix
                    # We output the DNA sequence of the ORF. 
                    # If on reverse strand, it should be the SEQUENCE of the ORF (e.g. TTA...), 
                    # which is what `nuc_seq[start:end]` contains because nuc_seq IS the reverse complement.
                    # So using nuc_seq slice is correct for content.
                    
                    # BUT wait, the test failed?
                    # Test Find: Expected ATGTTTTAA. Got ATGTTT.
                    # My split('*') logic: "current_aa_idx += len(part) + 1".
                    # But the 'part' excludes the stop codon.
                    # So 'len(part) + 1' accounts for length + stop.
                    # BUT orf_prot is part[start_aa_rel:].
                    # And orf_len_nt = orf_len_aa * 3.
                    # Does orf_prot include '*' ? NO. split remove it.
                    # So my calculated length DOES NOT include the stop codon.
                    # But normally an ORF includes the stop codon?
                    # The Perl script regex: `m/^.{$frame}(taa|tga|tag)/i`. Finding stop.
                    # Then: `push @ends,pos($strand{$direction})-2;` (position of last base BEFORE stop?)
                    # `subseq($s,$e)`.
                    # Let's check if Perl script includes stop.
                    # Usually "ORF" implies Start -> Stop (inclusive).
                    # My `part` split REMOVES stop.
                    # If I want to include stop, I must add 3 nt to length.
                    
                    final_seq_str = str(nuc_seq[start_nt_local : end_nt_local + 3]) # +3 for stop
                    # Be careful if stop was at very end of frame?
                    
                    # Re-calculating coords to include stop
                    orf_len_nt += 3
                    
                    if strand == 1:
                         final_end += 3
                    else:
                         final_start -= 3 # expanding left side in Fwd coords
                    
                    # Boundaries check
                    final_start = max(0, final_start)
                    final_end = min(seq_len, final_end)

                    best_orf = {
                        'start': final_start,
                        'end': final_end,
                        'strand': strand,
                        'len': orf_len_nt,
                        'dna': final_seq_str,
                        'prot': orf_prot + '*'
                    }
                
                current_aa_idx += len(part) + 1 # +1 for stop
                
    return best_orf

def run_find(args):
    """Replicates longestORF.pl"""
    print(f"Finding longest ORFs in {args.infile}...")
    
    try:
        sequences = SeqIO.parse(args.infile, "fasta")
        
        if args.outfile:
            out_handle = open(args.outfile, 'w')
        else:
            out_handle = sys.stdout
            
        count = 0
        strict_mode = not args.notstrict
        
        for seq in sequences:
            result = find_longest_orf(seq, strict=strict_mode)
            
            if result:
                # Output format:
                # >ID_Start_End_Strand
                strand_sym = '+' if result['strand'] == 1 else '-'
                # 1-based coords for display in ID
                disp_start = result['start'] + 1
                disp_end = result['end']
                
                new_id = f"{seq.id}_{disp_start}_{disp_end}_{strand_sym}"
                new_desc = f"Length={result['len']} Strand={strand_sym} Frame=?"
                
                out_rec = SeqRecord(Seq(result['dna']), id=new_id, description=new_desc)
                SeqIO.write(out_rec, out_handle, "fasta")
                count += 1
            else:
                print(f"No ORF found for {seq.id}")
        
        if args.outfile: out_handle.close()
        print(f"Processed {count} sequences.")
        
    except Exception as e:
        print(f"Error: {e}")

def run_extract(args):
    """Replicates ORFfromORFpredictor.pl - Extracts regions based on coords."""
    print(f"Extracting regions from {args.infile} using {args.coords}...")
    
    try:
        # Load Coords
        # Format expected: ID Start End [Strand]
        # Tab or space delimited
        
        coords_map = {}
        with open(args.coords, 'r') as f:
            for line in f:
                if not line.strip(): continue
                parts = line.split()
                if len(parts) >= 3:
                    sid = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    strand = '+'
                    if len(parts) > 3: strand = parts[3]
                    
                    coords_map[sid] = (start, end, strand)
                    
        sequences = SeqIO.to_dict(SeqIO.parse(args.infile, "fasta"))
        
        if args.outfile:
            out_handle = open(args.outfile, 'w')
        else:
            out_handle = sys.stdout
            
        count = 0
        
        for sid, (start, end, strand) in coords_map.items():
            if sid in sequences:
                seq_obj = sequences[sid]
                
                # Check bounds
                # Assuming 1-based inclusive input (standard for GFF/Perl scripts)
                # Convert to 0-based exclusive
                s_idx = start - 1
                e_idx = end
                
                sub_seq = seq_obj.seq[s_idx:e_idx]
                
                if strand == '-' or strand == '-1':
                    try:
                        sub_seq = sub_seq.reverse_complement()
                    except Exception as e:
                       print(f"Error extracting RC: {e}")
                    
                new_id = f"{sid}_{start}_{end}_{strand}"
                out_rec = SeqRecord(sub_seq, id=new_id, description="")
                SeqIO.write(out_rec, out_handle, "fasta")
                count += 1
            else:
                print(f"Warning: ID {sid} not found in sequence file.")
                
        if args.outfile: out_handle.close()
        print(f"Extracted {count} regions.")

    except Exception as e:
        print(f"Error: {e}")

def main():
    parser = argparse.ArgumentParser(description="ORF Tools")
    subparsers = parser.add_subparsers(dest='command', help='Sub-commands')
    
    # Find
    parser_find = subparsers.add_parser('find', help='Find longest ORF in each sequence')
    parser_find.add_argument('-in', '--input', dest='infile', required=True, help="Input FASTA")
    parser_find.add_argument('-out', '--output', dest='outfile', help="Output FASTA")
    parser_find.add_argument('-notstrict', action='store_true', help="Allow ORFs to start without ATG (starts at beginning of reading frame)")
    
    # Extract
    parser_ext = subparsers.add_parser('extract', help='Extract CDS/Regions from coordinates')
    parser_ext.add_argument('-in', '--input', dest='infile', required=True, help="Input FASTA")
    parser_ext.add_argument('-coords', dest='coords', required=True, help="Coordinates file (ID Start End [Strand])")
    parser_ext.add_argument('-out', '--output', dest='outfile', help="Output FASTA")

    args = parser.parse_args()
    
    if args.command == 'find':
        run_find(args)
    elif args.command == 'extract':
        run_extract(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
