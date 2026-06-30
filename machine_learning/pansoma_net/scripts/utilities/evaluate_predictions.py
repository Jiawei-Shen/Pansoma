import os
import json
import argparse

def generate_full_report(json_path, data_folder_path, comparison_file_path=None):
    """
    Calculates full classification metrics, lists all items, performs an optional
    comparison for false positives, and saves all results to output.json.

    Args:
        json_path (str): File path for the JSON with prediction lists.
        data_folder_path (str): File path for the root data folder.
        comparison_file_path (str, optional): Path to a .txt file for comparison.
    """
    # --- Initialize the output dictionary ---
    output_data = {}

    # --- 1. Load Model Predictions ---
    try:
        with open(json_path, 'r') as f:
            predictions = json.load(f)
        predicted_true = set(predictions.get('true', []))
        predicted_false = set(predictions.get('false', []))
    except FileNotFoundError:
        print(f"Error: The file '{json_path}' was not found.")
        return
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from '{json_path}'.")
        return

    # --- 2. Get Ground Truth Files ---
    actual_true_dir = os.path.join(data_folder_path, 'val', 'true')
    actual_false_dir = os.path.join(data_folder_path, 'val', 'false')
    if not os.path.isdir(actual_true_dir) or not os.path.isdir(actual_false_dir):
        print(f"Error: Ensure both '{actual_true_dir}' and '{actual_false_dir}' exist.")
        return
    actual_true = {os.path.basename(f) for f in os.listdir(actual_true_dir)}
    actual_false = {os.path.basename(f) for f in os.listdir(actual_false_dir)}

    # --- 3. Identify Classification Categories ---
    true_positives = sorted(list(predicted_true.intersection(actual_true)))
    true_negatives = sorted(list(predicted_false.intersection(actual_false)))
    false_positives = sorted(list(predicted_true.intersection(actual_false)))
    false_negatives = sorted(list(predicted_false.intersection(actual_true)))

    # --- 4. Calculate Metrics ---
    tp_count, tn_count = len(true_positives), len(true_negatives)
    fp_count, fn_count = len(false_positives), len(false_negatives)

    precision = tp_count / (tp_count + fp_count) if (tp_count + fp_count) > 0 else 0
    recall = tp_count / (tp_count + fn_count) if (tp_count + fn_count) > 0 else 0
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

    output_data['performance_metrics'] = {
        'precision': round(precision, 4),
        'recall': round(recall, 4),
        'f1_score': round(f1_score, 4)
    }
    output_data['confusion_matrix_counts'] = {
        'true_positives': tp_count,
        'true_negatives': tn_count,
        'false_positives': fp_count,
        'false_negatives': fn_count
    }
    output_data['detailed_lists'] = {
        'true_positives': true_positives,
        'true_negatives': true_negatives,
        'false_positives': false_positives,
        'false_negatives': false_negatives
    }

    # --- 5. (OPTIONAL) Compare False Positives ---
    if comparison_file_path:
        if not os.path.isfile(comparison_file_path):
            print(f"\nWarning: Comparison file not found at '{comparison_file_path}'")
            output_data['false_positive_comparison'] = {'error': 'Comparison file not found'}
        else:
            with open(comparison_file_path, 'r') as f:
                comparison_files = {line.strip() for line in f if line.strip()}

            overlapped_items = sorted(list(set(false_positives).intersection(comparison_files)))
            not_overlapped_items = sorted(list(set(false_positives).difference(comparison_files)))
            overlap_rate = len(overlapped_items) / fp_count if fp_count > 0 else 0

            output_data['false_positive_comparison'] = {
                'comparison_file': os.path.basename(comparison_file_path),
                'total_false_positives': fp_count,
                'total_items_in_comparison_file': len(comparison_files),
                'overlap_rate': round(overlap_rate, 4),
                'overlapped_items_count': len(overlapped_items),
                'not_overlapped_items_count': len(not_overlapped_items),
                'overlapped_items': overlapped_items,
                'not_overlapped_items': not_overlapped_items
            }

    # --- 6. Save Results to JSON File ---
    output_filename = 'output.json'
    with open(output_filename, 'w') as f:
        json.dump(output_data, f, indent=4)

    print(f"\n✅ Successfully generated report. All results saved to '{output_filename}'.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Run a full model evaluation and save the report to output.json."
    )
    parser.add_argument("json_file", help="Path to the JSON file with prediction results.")
    parser.add_argument("data_folder", help="Path to the root data folder.")
    parser.add_argument(
        "--compare_file",
        help="Optional: Path to a .txt file with filenames to compare False Positives against."
    )

    args = parser.parse_args()
    generate_full_report(args.json_file, args.data_folder, args.compare_file)