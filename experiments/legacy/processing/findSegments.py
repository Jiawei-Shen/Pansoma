import subprocess
import json
import sys
import argparse


def run_vg_find_position(position, graph_path):
    """Find the node corresponding to a given position."""
    # cmd = f'/home/jiawei/anaconda3/envs/pange/bin/vg find -x ./Pangenome_Graph/hprc-v1.1-mc-grch38.d9.gbz -p {position} | vg view -j -'
    cmd = f'vg find -x {graph_path} -p {position} | vg view -j -'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print("Error running vg find for position:", result.stderr)
        return None
    return json.loads(result.stdout)


def run_vg_find_neighbors(node_id, graph_path, context=5):
    """Find neighboring nodes within the given context range."""
    # cmd = f'/home/jiawei/anaconda3/envs/pange/bin/vg find -x {graph_path} -n {node_id} -c {context} | /home/jiawei/anaconda3/envs/pange/bin/vg view -j -'
    cmd = f'vg find -x {graph_path} -n {node_id} -c {context} | vg view -j -'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print("Error running vg find for neighbors:", result.stderr)
        return None
    return json.loads(result.stdout)


def main():
    parser = argparse.ArgumentParser(description="Find segments in a graph based on a given position.")
    parser.add_argument('position', type=str, help='Position to search in the graph')
    parser.add_argument('graph_path', default="./Pangenome_Graph/hprc-v1.1-mc-grch38.d9.gbz", type=str, help='Path to the graph file')
    args = parser.parse_args()

    position_data = run_vg_find_position(args.position, args.graph_path)
    if not position_data or "node" not in position_data:
        print("No node found for given position.")
        sys.exit(1)

    node_id = position_data["node"][0]["id"]
    print(f"Found node ID: {node_id}")

    neighbors_data = run_vg_find_neighbors(node_id, args.graph_path)
    if not neighbors_data:
        print("No neighboring nodes found.")
        sys.exit(1)

    print(json.dumps(neighbors_data, indent=4))


if __name__ == "__main__":
    main()
