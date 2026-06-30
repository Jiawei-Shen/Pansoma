#!/usr/bin/env python3
import argparse
import struct
import json
import os
import sys
import time
import math  # For math.ceil
# import numpy as np # Not strictly needed if outputting lists of lists for JSON
from collections import defaultdict
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
import psutil  # For memory usage reporting

# --- Constants ---
RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4,
                 'a': 0, 'c': 1, 'g': 2, 't': 3, 'n': 4}
N_INDEX = BASE_TO_INDEX['N']
DEFAULT_WINDOW_SIZE = 100

worker_dat_file = None  # Global for worker


# --- Helper Functions (assumed correct from previous versions) ---
def reverse_complement(sequence):
    if not isinstance(sequence, str): return ""
    try:
        complement_map = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
        return sequence.translate(complement_map)[::-1]
    except Exception:
        return ""


def parse_cigar(cigar_string):
    if not cigar_string or cigar_string == '*': return []
    try:
        return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHPX=])', cigar_string)]
    except Exception as e:
        print(f"⚠️ Warning: Could not parse CIGAR string '{cigar_string}': {e}", file=sys.stderr)
        return []


def detect_variants_core(ref_seq, read_seq, start_offset, cigar_ops):
    ref_len = len(ref_seq);
    read_len = len(read_seq)
    ref_ptr = start_offset;
    read_ptr = 0
    last_valid_ref_ptr = start_offset - 1
    try:
        for length, op in cigar_ops:
            if op in ('M', '=', 'X'):
                for i in range(length):
                    current_ref_pos = ref_ptr + i;
                    current_read_pos = read_ptr + i
                    if current_ref_pos >= ref_len or current_read_pos >= read_len: return
                    ref_base = ref_seq[current_ref_pos].upper();
                    read_base = read_seq[current_read_pos].upper()
                    last_valid_ref_ptr = current_ref_pos
                    if op != '=' and ref_base != 'N' and read_base != 'N' and ref_base != read_base:
                        yield (current_ref_pos, 'X', ref_base, read_base)
                ref_ptr += length;
                read_ptr += length
            elif op == 'I':
                if read_ptr + length > read_len: return
                inserted_bases = read_seq[read_ptr: read_ptr + length].upper()
                yield (last_valid_ref_ptr, 'I', '-', inserted_bases)
                read_ptr += length
            elif op == 'D':
                start_del_pos = ref_ptr
                if ref_ptr + length > ref_len: return
                deleted_bases = ref_seq[ref_ptr: ref_ptr + length].upper()
                yield (start_del_pos, 'D', deleted_bases, '-')
                ref_ptr += length;
                last_valid_ref_ptr = ref_ptr - 1
            elif op == 'S':
                read_ptr += length
            elif op == 'N':
                ref_ptr += length; last_valid_ref_ptr = ref_ptr - 1
            elif op in ('H', 'P'):
                pass
    except IndexError:
        print(
            f"⚠️ Warning: IndexError during variant detection. RefLen={ref_len}, ReadLen={read_len}, Offset={start_offset}, CIGAR='{''.join(map(str, cigar_ops))}', RefPtr={ref_ptr}, ReadPtr={read_ptr}",
            file=sys.stderr)
    except Exception as e:
        print(f"❌ Error: Unexpected error in detect_variants_core: {e}", file=sys.stderr)


def create_pileup_row(read_seq, read_offset, cigar_ops, window_start, window_end):
    window_len = window_end - window_start
    if window_len <= 0: return None
    pileup_row = [N_INDEX] * window_len;
    read_len = len(read_seq)
    ref_ptr = read_offset;
    read_ptr = 0;
    row_filled = False
    try:
        for length, op in cigar_ops:
            if ref_ptr >= window_end and not row_filled: return None
            if op in ('M', '=', 'X'):
                for i in range(length):
                    current_ref_pos = ref_ptr + i;
                    current_read_pos = read_ptr + i
                    if current_read_pos >= read_len: return pileup_row if row_filled else None
                    if window_start <= current_ref_pos < window_end:
                        pileup_idx = current_ref_pos - window_start
                        pileup_row[pileup_idx] = BASE_TO_INDEX.get(read_seq[current_read_pos].upper(), N_INDEX)
                        row_filled = True
                ref_ptr += length;
                read_ptr += length
            elif op == 'I':
                if read_ptr + length > read_len: return pileup_row if row_filled else None
                read_ptr += length
            elif op == 'D':
                for i in range(length):
                    current_ref_pos = ref_ptr + i
                    if window_start <= current_ref_pos < window_end:
                        pileup_row[current_ref_pos - window_start] = N_INDEX;
                        row_filled = True
                ref_ptr += length
            elif op == 'S':
                read_ptr += length
            elif op == 'N':
                for i in range(length):
                    current_ref_pos = ref_ptr + i
                    if window_start <= current_ref_pos < window_end:
                        pileup_row[current_ref_pos - window_start] = N_INDEX;
                        row_filled = True
                ref_ptr += length
            elif op in ('H', 'P'):
                pass
        return pileup_row if row_filled else None
    except IndexError:
        print(
            f"⚠️ Warning: IndexError during pileup creation. ReadLen={read_len}, Offset={read_offset}, CIGAR='{''.join(map(str, cigar_ops))}', RefPtr={ref_ptr}, ReadPtr={read_ptr}, Win=[{window_start}:{window_end}]",
            file=sys.stderr); return pileup_row if row_filled else None
    except Exception as e:
        print(f"❌ Error: Unexpected error in create_pileup_row: {e}", file=sys.stderr); return None


def parse_idx_file(idx_path):
    node_index = {}
    try:
        with open(idx_path, 'rb') as f:
            file_size = os.fstat(f.fileno()).st_size
            if file_size < 4: print(f"❌ Error: Index file {idx_path} is too small.", file=sys.stderr); sys.exit(1)
            num_nodes_bytes = f.read(4)
            if len(num_nodes_bytes) < 4: print(f"❌ Error: Could not read number of nodes from {idx_path}.",
                                               file=sys.stderr); sys.exit(1)
            num_nodes_header = struct.unpack('<I', num_nodes_bytes)[0]
            expected_min_size = 4 + num_nodes_header * 22
            if file_size < expected_min_size: print(
                f"⚠️ Warning: Index file {idx_path} may be truncated. Expected {expected_min_size}, found {file_size}.",
                file=sys.stderr)
            print(f"🔹 Reading index for up to {num_nodes_header} nodes (from header)...")
            nodes_parsed_count = 0
            for i in range(num_nodes_header):
                record_bytes = f.read(22)
                if len(record_bytes) < 22: print(
                    f"❌ Error: Index file ended prematurely reading record {i + 1}/{num_nodes_header}.",
                    file=sys.stderr); break
                node_id, offset, _, n_records, _ = struct.unpack('<I Q I I H', record_bytes)
                if n_records > 0: node_index[node_id] = (offset, n_records); nodes_parsed_count += 1
            print(f"✔ Parsed {nodes_parsed_count} nodes with records from index file.")
            if nodes_parsed_count == 0 and num_nodes_header > 0:
                print(f"⚠️ Warning: Header {num_nodes_header} nodes, but 0 with records found.", file=sys.stderr)
            elif len(node_index) != num_nodes_header and nodes_parsed_count > 0 and file_size >= expected_min_size:
                print(f"⚠️ Warning: Header {num_nodes_header} nodes, parsed {nodes_parsed_count} with records.",
                      file=sys.stderr)
            return node_index
    except FileNotFoundError:
        print(f"❌ Error: Index file not found at {idx_path}", file=sys.stderr); sys.exit(1)
    except Exception as e:
        print(f"❌ Error parsing index file {idx_path}: {e}", file=sys.stderr); sys.exit(1)


def load_node_sequences_from_gfa(gfa_path, target_node_ids_set):
    node_sequences = {};
    nodes_found_count = 0
    current_target_ids_set = set(target_node_ids_set)
    try:
        with open(gfa_path, 'r') as f:
            print(f"🔹 Reading GFA file: {gfa_path} for up to {len(current_target_ids_set)} target sequences...")
            line_counter = 0
            for line in f:
                line_counter += 1
                if not current_target_ids_set: break
                if line_counter % 1_000_000 == 0: print(
                    f"  GFA Checked {line_counter:,} lines... found {nodes_found_count}/{len(target_node_ids_set)} target sequences.",
                    end='\r')
                if not line.startswith('S\t'): continue
                parts = line.strip().split('\t')
                if len(parts) < 3: continue
                try:
                    nid = int(parts[1])
                except ValueError:
                    continue
                if nid in current_target_ids_set:
                    seq = parts[2].upper()
                    if seq != '*' and re.match(r'^[ACGTN]+$', seq):
                        node_sequences[nid] = seq;
                        nodes_found_count += 1
                        current_target_ids_set.remove(nid)
            print(f"\n✔ Finished GFA. Checked {line_counter:,} lines.")
            print(f"✔ Loaded sequences for {len(node_sequences)} target nodes from GFA.")
            if current_target_ids_set: print(
                f"⚠️ Warning (GFA): Did not find sequences for {len(current_target_ids_set)} target node(s), e.g., {list(current_target_ids_set)[:5]}.",
                file=sys.stderr)
    except FileNotFoundError:
        print(f"❌ Error: GFA file not found at {gfa_path}", file=sys.stderr); sys.exit(1)
    except Exception as e:
        print(f"❌ Error reading GFA file {gfa_path}: {e}", file=sys.stderr); sys.exit(1)
    return node_sequences


def init_worker(dat_file_path):
    global worker_dat_file
    try:
        worker_dat_file = open(dat_file_path, 'rb')
    except FileNotFoundError:
        print(f"❌ Error [Worker {os.getpid()}]: DAT file not found: {dat_file_path}", file=sys.stderr); sys.exit(1)
    except Exception as e:
        print(f"❌ Error [Worker {os.getpid()}] opening DAT: {e}", file=sys.stderr); sys.exit(1)


def process_node_parallel(task_args):
    node_id, file_read_offset, n_records, node_sequence, window_size = task_args
    global worker_dat_file
    if worker_dat_file is None: return node_id, {"error": "worker_dat_file not initialized"}
    if not node_sequence: return node_id, {"error": "node_sequence is empty"}

    node_len = len(node_sequence);
    half_window = window_size // 2
    reads_by_variant = defaultdict(list);
    mapq_filter = 10
    final_pileups = {}
    try:
        worker_dat_file.seek(file_read_offset)
        for _ in range(n_records):
            try:
                data = worker_dat_file.read(RECORD_SIZE)
                if len(data) < RECORD_SIZE: break
                read_offset_on_node, raw_seq, _, raw_cigar, mapq, strand_byte = RECORD_STRUCT.unpack(data)
                if mapq < mapq_filter: continue
                seq = raw_seq.rstrip(b'\x00').decode('ascii', errors='ignore')
                cigar_str = raw_cigar.rstrip(b'\x00').decode('ascii', errors='ignore')
                strand_char = strand_byte.decode('ascii', errors='ignore')
                if not seq or not cigar_str or cigar_str == '*': continue
                oriented_seq = seq;
                current_offset = read_offset_on_node
                if strand_char == '-':
                    oriented_seq = reverse_complement(seq)
                    if not oriented_seq: continue
                cigar_ops = parse_cigar(cigar_str)
                if not cigar_ops: continue
                variants_found_in_read = list(
                    detect_variants_core(node_sequence, oriented_seq, current_offset, cigar_ops))
                for vpos, vtype, ref_base, alt_base in variants_found_in_read:
                    if (vtype == 'I' and -1 <= vpos < node_len) or \
                            (vtype in ('D', 'X') and 0 <= vpos < node_len):
                        variant_key = f"{vpos}_{vtype}_{ref_base}_{alt_base}"
                        reads_by_variant[variant_key].append((current_offset, oriented_seq, cigar_ops))
            except struct.error:
                break
            except UnicodeDecodeError:
                continue
            except Exception:
                continue
    except IOError as e:
        return node_id, {"error": f"IOError: {e}"}
    except Exception as e:
        return node_id, {"error": f"Outer error: {e}"}

    if not reads_by_variant: return node_id, {}
    for variant_key, supporting_reads in reads_by_variant.items():
        if not supporting_reads: continue
        try:
            vpos = int(variant_key.split('_')[0])
        except (ValueError, IndexError):
            continue
        window_start, window_end = 0, 0
        if node_len <= window_size:
            window_start, window_end = 0, node_len
        else:
            center_pos = vpos;
            window_start = max(0, center_pos - half_window)
            window_end = min(node_len, window_start + window_size)
            window_start = max(0, window_end - window_size)
        if window_start >= window_end: continue
        pileup_rows = []
        for ro, rs, rco in supporting_reads:
            row = create_pileup_row(rs, ro, rco, window_start, window_end)
            if row: pileup_rows.append(row)
        if pileup_rows: final_pileups[variant_key] = pileup_rows
    return node_id, final_pileups


def generate_task_args(node_info_list, node_sequences_dict, window_size_arg):
    tasks_yielded_count = 0
    for node_id, file_offset, n_records in node_info_list:
        sequence = node_sequences_dict.get(node_id)
        if sequence:
            tasks_yielded_count += 1
            yield (node_id, file_offset, n_records, sequence, window_size_arg)
    print(f"\nℹ️ Task generator finished. Yielded {tasks_yielded_count} tasks for processing.")


def write_results_to_jsonl(filepath, results_batch_dict):
    if not results_batch_dict: return 0
    count_written_this_batch = 0
    try:
        with open(filepath, 'a') as outfile:
            for node_id, pileup_data in results_batch_dict.items():
                record_to_write = {str(node_id): pileup_data}
                json.dump(record_to_write, outfile)
                outfile.write('\n')
                count_written_this_batch += 1
    except Exception as e:
        print(f"\n❌ Error writing results batch to {filepath}: {e}", file=sys.stderr)
    return count_written_this_batch


# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Pileup generation.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("dat", help=".dat file")
    parser.add_argument("idx", help=".idx file")
    parser.add_argument("output", help="JSON Lines output file (.jsonl recommended)")
    parser.add_argument("--gfa", help="GFA file (if not using cache)")
    parser.add_argument("--load-cache", help="JSON cache file for node sequences")
    parser.add_argument("--save-cache", help="JSON cache file to save node sequences")
    parser.add_argument("--window", type=int, default=DEFAULT_WINDOW_SIZE, help="Pileup window size")
    parser.add_argument("-w", "--workers", type=int, default=os.cpu_count(), help="Number of worker processes")
    parser.add_argument("--progress-interval", type=int, default=1000,
                        help="Write results to file and print progress every N completed tasks.")
    parser.add_argument("--max-active-futures", type=int, default=0,
                        help="Max futures in flight (0 for num_workers * 2, min 1)")

    args = parser.parse_args()

    if not os.path.isfile(args.dat): print(f"❌ DAT not found: {args.dat}", file=sys.stderr); sys.exit(1)
    if not os.path.isfile(args.idx): print(f"❌ IDX not found: {args.idx}", file=sys.stderr); sys.exit(1)
    if not args.load_cache and not args.gfa: print("❌ Must provide --gfa or --load-cache.", file=sys.stderr); sys.exit(
        1)
    if args.load_cache and not os.path.isfile(args.load_cache) and not args.gfa: print(
        f"❌ Cache {args.load_cache} not found and no GFA fallback.", file=sys.stderr); sys.exit(1)
    if args.gfa and not os.path.isfile(args.gfa) and not (args.load_cache and os.path.isfile(args.load_cache)): print(
        f"❌ GFA {args.gfa} not found and no valid cache.", file=sys.stderr); sys.exit(1)
    if args.window <= 0: print("❌ Window must be > 0.", file=sys.stderr); sys.exit(1)

    print(f"Script started at {time.strftime('%Y-%m-%d %H:%M:%S')}")
    overall_start_time = time.time()
    main_process = psutil.Process(os.getpid())

    try:
        with open(args.output, 'w') as f:
            f.write("")
        print(f"🔹 Output file {args.output} initialized/cleared.")
    except Exception as e:
        print(f"❌ Error initializing output file {args.output}: {e}", file=sys.stderr); sys.exit(1)

    node_index_dict = parse_idx_file(args.idx)
    if not node_index_dict: print("❌ No nodes in index. Exiting.", file=sys.stderr); sys.exit(1)
    all_indexed_nodes_info = [(nid, data[0], data[1]) for nid, data in node_index_dict.items()]
    if not all_indexed_nodes_info: print("❌ No processable nodes in index. Exiting.", file=sys.stderr); sys.exit(1)
    target_ids_for_sequences = {info[0] for info in all_indexed_nodes_info}

    all_node_sequences = {}
    loaded_from_cache = False
    if args.load_cache and os.path.isfile(args.load_cache):
        print(f"🔹 Loading node sequences from cache: {args.load_cache}...")
        cache_load_start = time.time()
        try:
            with open(args.load_cache, 'r') as cf:
                loaded_data = json.load(cf)
            for k_str, seq_val in loaded_data.items():
                try:
                    k_int = int(k_str)
                except ValueError:
                    continue
                if k_int in target_ids_for_sequences: all_node_sequences[k_int] = seq_val
            print(
                f"✔ Loaded {len(all_node_sequences)} relevant sequences from cache in {time.time() - cache_load_start:.2f}s.")
            loaded_from_cache = True
        except Exception as e:
            print(f"❌ Error loading/parsing cache {args.load_cache}: {e}. Will try GFA if provided.", file=sys.stderr)
            all_node_sequences = {}
    if not loaded_from_cache or len(all_node_sequences) < len(target_ids_for_sequences):
        if args.gfa:
            print(
                f"🔹 {'Attempting to load remaining' if loaded_from_cache and len(all_node_sequences) > 0 else 'Loading ALL'} sequences from GFA: {args.gfa}...")
            still_needed_ids = target_ids_for_sequences - set(
                all_node_sequences.keys()) if loaded_from_cache else target_ids_for_sequences
            if still_needed_ids:
                gfa_sequences = load_node_sequences_from_gfa(args.gfa, still_needed_ids)
                all_node_sequences.update(gfa_sequences)
            elif not still_needed_ids and loaded_from_cache:
                print("ℹ️ All required sequences were already found in cache.")
            if args.save_cache and (not loaded_from_cache or (still_needed_ids and gfa_sequences)):
                print(f"🔹 Saving {len(all_node_sequences)} sequences to cache: {args.save_cache}...")
                try:
                    with open(args.save_cache, 'w') as cf:
                        json.dump({str(k): v for k, v in all_node_sequences.items()}, cf, indent=2)
                    print(f"✔ Saved cache to {args.save_cache}")
                except Exception as e:
                    print(f"❌ Error saving cache: {e}", file=sys.stderr)
        elif not all_node_sequences:
            print("❌ No sequence source for required nodes. Exiting.", file=sys.stderr); sys.exit(1)
    if not all_node_sequences: print("❌ Failed to load any required sequences. Exiting.", file=sys.stderr); sys.exit(1)

    process_node_info_list = [info for info in all_indexed_nodes_info if info[0] in all_node_sequences]
    total_tasks_to_process = len(process_node_info_list)
    if not process_node_info_list: print(
        "❌ No nodes to process after matching index with sequences. Exiting."); sys.exit(1)
    if total_tasks_to_process < len(all_indexed_nodes_info): print(
        f"ℹ️ Will process {total_tasks_to_process} nodes with sequences (out of {len(all_indexed_nodes_info)} indexed).")

    task_arg_generator = generate_task_args(process_node_info_list, all_node_sequences, args.window)

    total_nodes_with_pileups_written = 0
    results_this_interval = {}
    nodes_processed_since_last_write = 0
    future_to_node_id = {}
    num_workers = min(args.workers, os.cpu_count() or 1, total_tasks_to_process);
    num_workers = max(1, num_workers)
    max_active_futures = args.max_active_futures if args.max_active_futures > 0 else num_workers * 2;
    max_active_futures = max(1, max_active_futures)
    progress_write_interval = args.progress_interval

    print(f"🔹 Concurrently submitting and processing {total_tasks_to_process} tasks using {num_workers} workers.")
    print(f"   Max active futures to manage: {max_active_futures}")
    print(f"   Will write to output file and print memory approx every {progress_write_interval} completed tasks.")

    nodes_completed_count = 0
    tasks_submitted_count = 0

    with ProcessPoolExecutor(max_workers=num_workers, initializer=init_worker, initargs=(args.dat,)) as executor:
        for _ in range(max_active_futures):  # Submit initial batch
            try:
                task_args = next(task_arg_generator)
                future = executor.submit(process_node_parallel, task_args)
                future_to_node_id[future] = task_args[0]
                tasks_submitted_count += 1
            except StopIteration:
                break
        if tasks_submitted_count > 0: print(f"  Initial submission: {tasks_submitted_count} tasks in flight.")

        while future_to_node_id:
            active_futures_list = list(future_to_node_id.keys())
            if not active_futures_list: break
            for future in as_completed(active_futures_list):
                nodes_completed_count += 1
                nodes_processed_since_last_write += 1
                original_node_id = future_to_node_id.pop(future)
                try:
                    _ret_node_id, pileup_dict = future.result()
                    if _ret_node_id != original_node_id: print(
                        f"⚠️ Main: Node ID mismatch! Original: {original_node_id}, Returned: {_ret_node_id}",
                        file=sys.stderr)
                    if pileup_dict and not (isinstance(pileup_dict, dict) and "error" in pileup_dict) and pileup_dict:
                        results_this_interval[original_node_id] = pileup_dict
                    elif isinstance(pileup_dict, dict) and "error" in pileup_dict:
                        print(f"  ℹ️ Worker for node {original_node_id} reported error: {pileup_dict['error']}",
                              file=sys.stderr)
                except Exception as exc:
                    print(f"\n❌ Task for node {original_node_id} generated exception during result retrieval: {exc}",
                          file=sys.stderr)

                if len(future_to_node_id) < max_active_futures:
                    try:
                        if tasks_submitted_count < total_tasks_to_process:
                            new_task_args = next(task_arg_generator)
                            new_future = executor.submit(process_node_parallel, new_task_args)
                            future_to_node_id[new_future] = new_task_args[0]
                            tasks_submitted_count += 1
                    except StopIteration:
                        pass

                should_write_now = (nodes_processed_since_last_write >= progress_write_interval) or \
                                   (tasks_submitted_count == total_tasks_to_process and not future_to_node_id) or \
                                   (nodes_completed_count == total_tasks_to_process)

                if results_this_interval and should_write_now:
                    num_written_this_batch = write_results_to_jsonl(args.output, results_this_interval)
                    total_nodes_with_pileups_written += num_written_this_batch
                    results_this_interval.clear()
                    nodes_processed_since_last_write = 0

                if nodes_completed_count % progress_write_interval == 0 or \
                        nodes_completed_count == total_tasks_to_process or \
                        not future_to_node_id:
                    elapsed = time.time() - overall_start_time
                    rate = nodes_completed_count / elapsed if elapsed > 0 else 0
                    mem_rss_mb = main_process.memory_info().rss / (1024 * 1024)
                    print(
                        f"  Progress: {nodes_completed_count}/{total_tasks_to_process} done. {len(future_to_node_id)} active. Rate: {rate:.1f}/s. Mem: {mem_rss_mb:.1f}MB. Written: {total_nodes_with_pileups_written}")
                break
            else:
                if not future_to_node_id: break

    if results_this_interval:  # Final write
        num_written_this_batch = write_results_to_jsonl(args.output, results_this_interval)
        total_nodes_with_pileups_written += num_written_this_batch
        results_this_interval.clear()

    print()
    total_script_elapsed = time.time() - overall_start_time
    mem_rss_mb = main_process.memory_info().rss / (1024 * 1024)
    print(
        f"\n🏁 Finished processing in {total_script_elapsed:.2f} seconds. Final main process memory: {mem_rss_mb:.1f}MB.")
    print(f"🔹 Found and wrote pileups for {total_nodes_with_pileups_written} nodes overall to {args.output}")
    print(f"✅ Done. Total time: {time.time() - overall_start_time:.2f} seconds.")


if __name__ == '__main__':
    main()