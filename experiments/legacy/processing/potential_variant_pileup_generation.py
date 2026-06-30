#!/usr/bin/env python3
import json
import numpy as np
import argparse
from collections import Counter

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def get_max_width(segments):
    """Determine the width of the pileup matrix needed."""
    max_end = 0
    for segment in segments:
        if segment['read_quality'] <= 10:
            continue
        seq = segment['sequence']
        offset = segment['offset']
        if segment['strand'] == "-":
            offset = offset  # keep raw offset; corrected in pileup
        end = offset + len(seq)
        max_end = max(max_end, end)
    return max_end


def pileup_char_matrix(data):
    """
    Create an adaptive character-based pileup matrix.
    Returns a matrix (numpy array) and the final width used.
    """
    valid_segments = [s for s in data if s['read_quality'] > 10]
    if not valid_segments:
        return np.array([[]], dtype='<U1')  # empty

    max_width = get_max_width(valid_segments)
    pileup = []

    for segment in valid_segments:
        offset = segment['offset']
        seq = segment['sequence']

        if segment['strand'] == "-":
            seq = reverse_complement(seq)
            offset = max_width - offset - len(seq)

        for i, base in enumerate(seq):
            x = offset + i
            if 0 <= x < max_width and base.upper() in "ACGT":
                placed = False
                for row in pileup:
                    if x >= len(row):
                        row.extend(['.'] * (x - len(row) + 1))
                    if row[x] == '.':
                        row[x] = base.upper()
                        placed = True
                        break
                if not placed:
                    new_row = ['.'] * max_width
                    new_row[x] = base.upper()
                    pileup.append(new_row)

    # Normalize row lengths
    for row in pileup:
        if len(row) < max_width:
            row.extend(['.'] * (max_width - len(row)))

    return np.array(pileup, dtype='<U1')


def calculate_allele_frequencies(pileup_matrix):
    if pileup_matrix.size == 0:
        return {}

    height, width = pileup_matrix.shape
    af_dict = {}

    for col in range(width):
        bases = [pileup_matrix[row][col] for row in range(height) if pileup_matrix[row][col] in "ACGT"]
        if bases:
            count = Counter(bases)
            total = sum(count.values())
            af = {base: round(freq / total, 4) for base, freq in count.items()}
            af_dict[str(col)] = af

    return af_dict


def process_json(input_path, output_path):
    with open(input_path) as f:
        raw_data = json.load(f)

    result = {}

    for node_id, segments in raw_data.items():
        matrix = pileup_char_matrix(segments)
        af = calculate_allele_frequencies(matrix)
        result[node_id] = {
            "pileup": ["".join(row) for row in matrix.tolist()],
            "allele_frequencies": af
        }

    with open(output_path, 'w') as f:
        json.dump(result, f, indent=2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Adaptive pileup with allele frequency per column.")
    parser.add_argument("input_json", help="Input JSON file")
    parser.add_argument("output_json", help="Output JSON file")
    args = parser.parse_args()

    process_json(args.input_json, args.output_json)
