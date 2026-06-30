#!/usr/bin/env python3
import argparse
import struct
import json
import os
import sys
import time
import numpy as np
from collections import defaultdict
import re # Import re at the top level
from concurrent.futures import ProcessPoolExecutor
# Maybe import cpu_count to default workers
# from os import cpu_count

# ─────────────────────────────────────────────────────────────────────────────
# Constants (moved globals here or made local where appropriate)
RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}

# Global for worker process state (file handle)
# Each worker process will have its own instance of this global
worker_dat_file = None

# ─────────────────────────────────────────────────────────────────────────────
# Helper Functions (remain mostly unchanged, ensure defined at top level)

def reverse_complement(sequence):
    complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(complement_map)[::-1]

def parse_idx_file(idx_path):
    node_index = {}
    try:
        with open(idx_path, 'rb') as f:
            # Check file size for basic validation
            file_size = os.fstat(f.fileno()).st_size
            if file_size < 4:
                 print(f"❌ Error: Index file {idx_path} is too small.", file=sys.stderr)
                 sys.exit(1)

            num_nodes_bytes = f.read(4)
            if len(num_nodes_bytes) < 4:
                print(f"❌ Error: Could not read number of nodes from {idx_path}.", file=sys.stderr)
                sys.exit(1)
            num_nodes = struct.unpack('<I', num_nodes_bytes)[0]

            # Check if file size matches expected size based on num_nodes
            expected_min_size = 4 + num_nodes * 22
            if file_size < expected_min_size:
                print(f"❌ Error: Index file {idx_path} appears truncated. Expected at least {expected_min_size} bytes for {num_nodes} nodes, found {file_size}.", file=sys.stderr)
                # Continue cautiously or exit, depending on desired robustness
                # sys.exit(1) # Or just warn

            print(f"🔹 Reading index for {num_nodes} nodes...")
            for i in range(num_nodes):
                record_bytes = f.read(22)
                if len(record_bytes) < 22:
                     print(f"❌ Error: Index file ended prematurely while reading record {i+1}/{num_nodes}.", file=sys.stderr)
                     break # Stop processing further records
                node_id, offset, block_size, n_records, _ = struct.unpack('<I Q I I H', record_bytes)
                node_index[node_id] = (offset, n_records)

        print(f"✔ Parsed {len(node_index)} nodes from index file.")
        if len(node_index) != num_nodes:
             print(f"⚠️ Warning: Expected {num_nodes} nodes based on header, but parsed {len(node_index)}. File might be corrupt or truncated.", file=sys.stderr)
        return node_index
    except FileNotFoundError:
        print(f"❌ Error: Index file not found at {idx_path}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error parsing index file {idx_path}: {e}", file=sys.stderr)
        sys.exit(1)


def load_node_sequences_from_gfa(gfa_path, target_node_ids):
    node_sequences = {}
    target_node_set = set(target_node_ids) # Use a set for faster lookups
    try:
        with open(gfa_path, 'r') as f:
            print(f"🔹 Reading GFA file: {gfa_path}")
            line_counter = 0
            parsed_count = 0
            for line in f:
                line_counter += 1
                if line_counter % 10_000_000 == 0:
                    print(f"  Checked {line_counter:,} lines in GFA file...")

                if not line.startswith('S\t'):
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 3:
                    # print(f"⚠️ Warning: Skipping malformed S line {line_counter}: {line.strip()}", file=sys.stderr)
                    continue # Skip malformed lines

                try:
                    # Handle potential non-integer node IDs gracefully if needed
                    # For now, assume they should be integers
                    nid = int(parts[1])
                except ValueError:
                    # print(f"⚠️ Warning: Skipping S line {line_counter} with non-integer ID: {parts[1]}", file=sys.stderr)
                    continue # Skip nodes with non-integer IDs

                if nid in target_node_set:
                    # Store sequence, remove potential trailing newline just in case
                    node_sequences[nid] = parts[2]
                    parsed_count += 1
                    # Optional: remove from set to stop searching if all found
                    # target_node_set.remove(nid)
                    # if not target_node_set:
                    #    print(f"✔ Found all {len(target_node_ids)} target sequences.")
                    #    break

            print(f"✔ Checked {line_counter:,} lines in GFA.")
            print(f"✔ Loaded {len(node_sequences)} target sequences from GFA.")
            if len(node_sequences) != len(target_node_ids):
                 print(f"⚠️ Warning: Found sequences for {len(node_sequences)} out of {len(target_node_ids)} target nodes.", file=sys.stderr)

    except FileNotFoundError:
        print(f"❌ Error: GFA file not found at {gfa_path}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error reading GFA file {gfa_path}: {e}", file=sys.stderr)
        sys.exit(1)

    return node_sequences

def decode_cigar(cigar_string):
    # Using pre-compiled regex might be slightly faster if called extremely often,
    # but re.findall caches implicitly, so improvement might be negligible.
    # CIGAR_REGEX = re.compile(r'(\d+)([MXDI])')
    # return CIGAR_REGEX.findall(cigar_string)
    if not cigar_string or cigar_string == '*': # Handle empty or null CIGAR
        return []
    try:
        return re.findall(r'(\d+)([MIDX=])', cigar_string) # Include M, =, X
    except Exception as e:
        print(f"⚠️ Warning: Could not parse CIGAR string '{cigar_string}': {e}", file=sys.stderr)
        return []


def detect_variants_from_cigar(offset, cigar_string):
    variants = []
    pos = offset # Position on the reference node sequence
    cigar_ops = decode_cigar(cigar_string)

    for length_str, op in cigar_ops:
        try:
            length = int(length_str)
        except ValueError:
            print(f"⚠️ Warning: Invalid length in CIGAR operation '{length_str}{op}' in string '{cigar_string}'", file=sys.stderr)
            continue # Skip this operation

        if op == 'M' or op == '=': # Match or sequence match
            pos += length
        elif op == 'X': # Sequence mismatch
            # Report one variant event per base for simplicity here
            for i in range(length):
                variants.append((pos + i, 'X'))
            pos += length
        elif op == 'I': # Insertion to the reference
            # Variant occurs *after* the reference position `pos - 1`
            variants.append((pos - 1, 'I'))
            # Insertions do not consume reference bases, so pos does not change here.
        elif op == 'D': # Deletion from the reference
            # Deletion starts *at* the reference position `pos`
            variants.append((pos, 'D'))
            pos += length # Deletions consume reference bases
        # Other CIGAR ops like S, H, P, N are ignored for variant detection here

    return variants

# ─────────────────────────────────────────────────────────────────────────────
# Worker Process Initialization and Target Function

def init_worker(dat_file_path):
    """Initializer for each worker process: opens the .dat file."""
    global worker_dat_file
    # print(f"[Worker {os.getpid()}] Initializing and opening {dat_file_path}") # Debug print
    try:
        worker_dat_file = open(dat_file_path, 'rb')
    except FileNotFoundError:
         print(f"❌ Error [Worker {os.getpid()}]: DAT file not found at {dat_file_path}", file=sys.stderr)
         # Exit the worker process cleanly if the file isn't found
         sys.exit(1)
    except Exception as e:
        print(f"❌ Error [Worker {os.getpid()}] opening DAT file {dat_file_path}: {e}", file=sys.stderr)
        sys.exit(1)


def process_node_parallel(task_args):
    """
    Function executed by each worker process.
    Processes a single node to find variants and generate pileups.
    """
    node_id, offset, n_records, sequence = task_args
    global worker_dat_file

    if worker_dat_file is None:
         # This shouldn't happen if initializer worked, but good safeguard
         print(f"❌ Error [Worker {os.getpid()}]: Worker DAT file handle not initialized for node {node_id}.", file=sys.stderr)
         # Return an empty result or raise an error, depending on desired handling
         return node_id, {} # Return empty dict for this node

    # Make sure sequence is not None or empty if needed later
    if not sequence:
         # print(f"⚠️ Warning [Worker {os.getpid()}]: Empty sequence provided for node {node_id}. Skipping pileup generation.", file=sys.stderr)
         return node_id, {} # Cannot process without sequence

    node_len = len(sequence)
    segments = []

    try:
        # Seek to the correct position for this node's records
        # The first 10 bytes (node_id, offset, block_size) are metadata in the block
        worker_dat_file.seek(offset + 10)

        for record_idx in range(n_records):
            try:
                data = worker_dat_file.read(RECORD_SIZE)
                if len(data) < RECORD_SIZE:
                    # print(f"⚠️ Warning [Worker {os.getpid()}]: Short read ({len(data)} bytes) for node {node_id}, record {record_idx+1}/{n_records}. Stopping reads for this node.", file=sys.stderr)
                    break # Stop reading if data is incomplete

                # Unpack record data
                off, raw_seq, raw_bq, raw_cigar, mapq, strand_byte = RECORD_STRUCT.unpack(data)

                # Basic filtering
                if mapq < 10:
                    continue

                # Decode fields, handling potential errors
                try:
                    # Use 'replace' or 'ignore' for bytes that aren't valid ASCII
                    seq = raw_seq.rstrip(b'\x00').decode('ascii', errors='ignore')
                    cigar = raw_cigar.rstrip(b'\x00').decode('ascii', errors='ignore')
                    strand_char = strand_byte.decode('ascii')
                except UnicodeDecodeError as ude:
                    # print(f"⚠️ Warning [Worker {os.getpid()}]: Unicode decode error in record {record_idx+1} for node {node_id}: {ude}. Skipping record.", file=sys.stderr)
                    continue # Skip record if decoding fails

                read_len = len(seq)

                # Adjust offset and reverse complement sequence for reverse strand reads
                # Handle potential edge case where node_len might be 0 or off/read_len calculation is invalid
                if strand_char == '-' and node_len > 0:
                    # Original calculation: off = node_len - off - read_len
                    # Let's verify this logic. If a read maps to the reverse strand,
                    # the 'offset' usually refers to the *start* of the alignment on the
                    # forward strand representation of the node.
                    # If the alignment starts at `off` (0-based) on forward strand and has length `read_len`,
                    # on the reverse complement node sequence, the alignment effectively starts at
                    # `node_len - (off + read_len)`. Let's stick to the original calculation for now,
                    # assuming it matches the data's convention. Need example data to be 100% sure.
                     adj_off = node_len - off - read_len # Calculate adjusted offset first
                     if adj_off < 0:
                        # print(f"⚠️ Warning [Worker {os.getpid()}]: Negative adjusted offset ({adj_off}) for node {node_id}, record {record_idx+1}. Original off={off}, node_len={node_len}, read_len={read_len}. Using original offset.", file=sys.stderr)
                        # Decide how to handle: skip, use original `off`, use 0? Using original `off` might be wrong.
                        # Using 0 might be safer if the coordinate system is unclear. Let's use 0 for now.
                        # adj_off = 0
                        # Or perhaps better to skip this read if coordinates seem inconsistent
                        # print(f"   Skipping read due to inconsistent coordinates.")
                        continue # Skip this read

                     seq = reverse_complement(seq) # Reverse complement the sequence
                     off = adj_off # Use the adjusted offset for further calculations


                segments.append((off, seq, cigar))

            except struct.error as se:
                 print(f"❌ Error [Worker {os.getpid()}]: Failed to unpack record {record_idx+1} for node {node_id}: {se}. Stopping reads for this node.", file=sys.stderr)
                 break # Stop processing this node if unpacking fails
            except Exception as e_inner:
                 print(f"❌ Error [Worker {os.getpid()}]: Unexpected error processing record {record_idx+1} for node {node_id}: {e_inner}", file=sys.stderr)
                 # Continue to next record or break depending on severity
                 continue


    except IOError as ioe:
         print(f"❌ Error [Worker {os.getpid()}]: I/O error reading data for node {node_id} at offset {offset}: {ioe}", file=sys.stderr)
         return node_id, {} # Return empty result for this node
    except Exception as e_outer:
        # Catch any other unexpected error during the reading loop
        print(f"❌ Error [Worker {os.getpid()}]: Unexpected error seeking/reading node {node_id}: {e_outer}", file=sys.stderr)
        return node_id, {}


    # --- Variant Detection and Pileup Generation ---
    reads_by_variant = defaultdict(list)
    for off, seq, cigar in segments:
        # Ensure offset is non-negative before passing to variant detection
        if off < 0:
            # print(f"⚠️ Warning [Worker {os.getpid()}]: Skipping segment with negative offset ({off}) for node {node_id}.", file=sys.stderr)
            continue

        # Detect variants based on CIGAR string
        variants = detect_variants_from_cigar(off, cigar)
        for vpos, vtype in variants:
            # Ensure variant position is within node boundaries
            # For Insertions ('I'), vpos is the base *before* the insertion. Check vpos >= -1?
            # For Deletions ('D') and Mismatches ('X'), vpos is the first affected base. Check vpos >= 0.
            is_valid_pos = False
            if vtype == 'I' and vpos >= -1 and vpos < node_len: # Insertion after vpos
                 is_valid_pos = True
            elif vtype in ('D', 'X') and vpos >= 0 and vpos < node_len: # Deletion/Mismatch starting at vpos
                 is_valid_pos = True

            if is_valid_pos:
                reads_by_variant[(vpos, vtype)].append((off, seq))
            #else:
            #    print(f"Debug: Variant ({vpos}, {vtype}) outside node bounds [0, {node_len}) for node {node_id}")


    # --- Generate Pileups ---
    pileups = {}
    window = 60
    half = window // 2
    for (vpos, vtype), reads in reads_by_variant.items():
        if not reads: # Skip if no reads support this variant somehow
             continue

        mat = np.full((len(reads), window), BASE_TO_INDEX['N'], dtype=np.uint8) # Fill with 'N' index
        for i, (read_offset, read_seq) in enumerate(reads):
            read_len = len(read_seq)
            # Calculate the start position in the read sequence corresponding to the
            # start of the pileup window (vpos - half).
            # pileup_window_start_on_node = vpos - half
            # pileup_window_start_in_read = pileup_window_start_on_node - read_offset
            start_in_read = vpos - read_offset - half

            for j in range(window): # j is the index within the pileup window (0 to window-1)
                # Calculate the corresponding index in the read sequence
                read_idx = start_in_read + j
                if 0 <= read_idx < read_len:
                    base = read_seq[read_idx].upper()
                    mat[i, j] = BASE_TO_INDEX.get(base, BASE_TO_INDEX['N'])
                # else: base is outside the read, defaults to 'N' (index 4)

        # Store pileup matrix as list of lists for JSON serialization
        pileups[f"{vpos}_{vtype}"] = mat.tolist()

    # print(f"[Worker {os.getpid()}] Finished processing node {node_id}. Found {len(pileups)} variants.") # Debug print
    return node_id, pileups

# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Parallel variant-centered pileup generation.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter # Show defaults
    )
    parser.add_argument("dat", help=".dat file path (read alignment data)")
    parser.add_argument("idx", help=".idx file path (index for .dat file)")
    parser.add_argument("output", help="JSON output file path for pileups")
    parser.add_argument("--gfa", help="GFA graph file path (needed if node sequence cache is not used/built)")
    parser.add_argument("--load-cache", help="Load node sequences from this JSON cache file")
    parser.add_argument("--save-cache", help="Save node sequences to this JSON cache file (used if --gfa is provided)")
    parser.add_argument("-w", "--workers", type=int, default=8,
                        help="Number of worker processes to use")
    parser.add_argument("-c", "--chunksize", type=int, default=200,
                        help="Number of nodes processed by a worker before returning results (approx)")
    args = parser.parse_args()

    # --- Input Validation ---
    if not os.path.isfile(args.dat):
         print(f"❌ Error: DAT file not found: {args.dat}", file=sys.stderr)
         sys.exit(1)
    if not os.path.isfile(args.idx):
         print(f"❌ Error: Index file not found: {args.idx}", file=sys.stderr)
         sys.exit(1)

    if not args.load_cache and not args.gfa:
        print("❌ Error: You must provide either a GFA file (`--gfa`) to build node sequences or a pre-built cache file (`--load-cache`).", file=sys.stderr)
        sys.exit(1)
    if args.gfa and not args.save_cache and not args.load_cache:
        # If using GFA and no cache exists, force saving it for future runs.
        # Determine a default cache name or require --save-cache.
        # Let's require it for clarity.
        print("❌ Error: When providing a GFA file (`--gfa`) without loading a cache (`--load-cache`), you must also specify where to save the new cache using `--save-cache`.", file=sys.stderr)
        sys.exit(1)
    if args.load_cache and not os.path.isfile(args.load_cache):
        print(f"❌ Error: Specified cache file to load does not exist: {args.load_cache}", file=sys.stderr)
        sys.exit(1)
    if args.gfa and not os.path.isfile(args.gfa):
         print(f"❌ Error: GFA file not found: {args.gfa}", file=sys.stderr)
         sys.exit(1)

    # --- Load Index ---
    print("🔹 Parsing index file...")
    start_time = time.time()
    node_index = parse_idx_file(args.idx)
    if not node_index:
         print("❌ Error: Failed to parse node index. Exiting.", file=sys.stderr)
         sys.exit(1)
    print(f"✔ Index parsing took {time.time() - start_time:.2f} seconds.")


    # --- Load or Build Node Sequences ---
    node_sequences = {}
    cache_path = args.load_cache or args.save_cache # Determine the relevant cache path

    if args.load_cache and os.path.isfile(args.load_cache):
        print(f"🔹 Loading node sequences from cache: {args.load_cache}...")
        start_time = time.time()
        try:
            with open(args.load_cache, 'r') as cf:
                # Load keys as integers directly
                loaded_data = json.load(cf)
                node_sequences = {int(k): v for k, v in loaded_data.items()}
            print(f"✔ Loaded {len(node_sequences)} sequences from cache in {time.time() - start_time:.2f} seconds.")
        except json.JSONDecodeError as jde:
             print(f"❌ Error decoding JSON from cache file {args.load_cache}: {jde}", file=sys.stderr)
             sys.exit(1)
        except Exception as e:
             print(f"❌ Error loading cache file {args.load_cache}: {e}", file=sys.stderr)
             sys.exit(1)
    elif args.gfa:
        print(f"🔹 Building node sequence cache from GFA: {args.gfa}...")
        start_time = time.time()
        # We only need sequences for nodes present in the index file
        node_sequences = load_node_sequences_from_gfa(args.gfa, node_index.keys())
        print(f"✔ Sequence loading from GFA took {time.time() - start_time:.2f} seconds.")
        if args.save_cache:
            print(f"🔹 Saving node sequence cache to: {args.save_cache}...")
            start_time = time.time()
            try:
                with open(args.save_cache, 'w') as cf:
                    # Ensure keys are saved as strings in JSON
                    json.dump({str(k): v for k, v in node_sequences.items()}, cf)
                print(f"✔ Saved cache in {time.time() - start_time:.2f} seconds.")
            except IOError as ioe:
                 print(f"❌ Error saving cache file {args.save_cache}: {ioe}", file=sys.stderr)
                 # Continue execution even if saving cache fails? Or exit? Let's continue.
            except Exception as e:
                 print(f"❌ Error during cache saving to {args.save_cache}: {e}", file=sys.stderr)

    else:
        # This case should have been caught by earlier validation, but double-check
        print("❌ Internal Error: No sequence source available. Exiting.", file=sys.stderr)
        sys.exit(1)

    if not node_sequences:
         print("⚠️ Warning: No node sequences were loaded or built. Output will likely be empty.", file=sys.stderr)
         # Decide whether to exit or continue
         # sys.exit(1)


    # --- Prepare Tasks for Parallel Processing ---
    print("🔹 Preparing tasks for parallel processing...")
    tasks = []
    nodes_missing_sequence = 0
    for node_id, (offset, n_records) in node_index.items():
        sequence = node_sequences.get(node_id)
        if sequence is not None:
            # Add task only if sequence exists
            tasks.append((node_id, offset, n_records, sequence))
        else:
            # print(f"Node {node_id} found in index but not in sequence data. Skipping.")
            nodes_missing_sequence += 1

    if nodes_missing_sequence > 0:
        print(f"⚠️ Warning: Skipped {nodes_missing_sequence} nodes present in index but missing from sequence data.")

    if not tasks:
        print("❌ Error: No tasks to process (no nodes with available sequence data found). Exiting.", file=sys.stderr)
        sys.exit(1)

    total_tasks = len(tasks)
    num_workers = min(args.workers, total_tasks) # Don't need more workers than tasks
    if num_workers <= 0:
        num_workers = 1 # Ensure at least one worker

    print(f"🔹 Processing {total_tasks} nodes using {num_workers} workers (chunksize: {args.chunksize})...")

    # --- Execute in Parallel ---
    results = {}
    start_proc_time = time.time()

    # Use ProcessPoolExecutor with context manager for proper cleanup
    try:
        with ProcessPoolExecutor(max_workers=num_workers,
                                 initializer=init_worker,
                                 initargs=(args.dat,)) as executor:
            # map processes tasks in order and returns results as they complete
            # Pass the list of task tuples
            future_results = executor.map(process_node_parallel, tasks, chunksize=args.chunksize)

            processed_count = 0
            # Iterate through results as they become available
            for node_id, pileup_dict in future_results:
                 if pileup_dict is not None: # Store result if valid
                     results[node_id] = pileup_dict
                 # else: Error likely occurred in worker, message printed there.

                 processed_count += 1
                 # Print progress update periodically
                 if processed_count % 10000 == 0 or processed_count == total_tasks:
                     elapsed = time.time() - start_proc_time
                     rate = processed_count / elapsed if elapsed > 0 else 0
                     print(f"✔ {processed_count}/{total_tasks} nodes processed ({rate:.1f} nodes/s) — Elapsed: {elapsed:.2f}s")

    except Exception as pool_exc:
         print(f"\n❌ An error occurred during parallel processing: {pool_exc}", file=sys.stderr)
         # Depending on the error, 'results' might be partially filled.
         # Decide if you want to save partial results or exit.
         # For simplicity, we'll proceed to save whatever results we got.
         print("⚠️ Attempting to write any partial results obtained...", file=sys.stderr)


    total_elapsed_time = time.time() - start_proc_time
    print(f"✔ Parallel processing finished in {total_elapsed_time:.2f} seconds.")

    # --- Write Output ---
    print(f"🔹 Writing {len(results)} results to JSON output: {args.output}")
    start_write_time = time.time()
    try:
        with open(args.output, 'w') as out_f:
            # Save node IDs as strings for compatibility, consistent with cache saving
            json.dump({str(k): v for k, v in results.items()}, out_f, indent=2)
        write_elapsed_time = time.time() - start_write_time
        print(f"✔ Output written in {write_elapsed_time:.2f} seconds.")
        print(f"✅ Done. Output saved to {args.output}")
    except IOError as ioe:
         print(f"❌ Error writing output JSON to {args.output}: {ioe}", file=sys.stderr)
         sys.exit(1)
    except Exception as e:
         print(f"❌ Unexpected error writing output JSON: {e}", file=sys.stderr)
         sys.exit(1)


if __name__ == '__main__':
    # Good practice: protect the main execution block
    # This is necessary for multiprocessing on some platforms (like Windows)
    main()