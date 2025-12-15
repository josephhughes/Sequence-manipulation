import os
import subprocess
import sys

def create_dummy_fasta(filename):
    with open(filename, 'w') as f:
        f.write(">test1\nATGC\n")
        f.write(">test2\nCGTA\n")
        f.write(">dup\nAAAA\n")
        f.write(">dup\nAAAA\n")
        f.write(">long\nATGCATGCATGC\n")

def run_test(args, expected_file, check_content=None):
    cmd = [sys.executable, "python_port/seqmanip/select_seq.py"] + args
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    if os.path.exists(expected_file):
        print(f"PASS: {expected_file} created.")
        if check_content:
             with open(expected_file, 'r') as f:
                 content = f.read()
                 if check_content in content:
                     print(f"PASS: Content matched '{check_content}'")
                 else:
                     print(f"FAIL: Content mismatch. Found:\n{content}")
    else:
        print(f"FAIL: {expected_file} not created.")

if __name__ == "__main__":
    create_dummy_fasta("test_input.fasta")
    
    # Test 1: Exact Match
    run_test(["-in", "test_input.fasta", "-out", "out_exact.fasta", "-fid", "test1"], 
             "out_exact.fasta", ">test1")

    # Test 2: Longest
    run_test(["-in", "test_input.fasta", "-out", "out_long", "-longest"],
             "out_long_longest.fa", ">long")
             
    # Test 3: Deduplicate
    run_test(["-in", "test_input.fasta", "-out", "out_dedup", "-dedup"],
             "out_dedup.nodups.fa")

    # Clean up
    # os.remove("test_input.fasta")
