import re
from os import fdopen
from Bio import Phylo

def to_newick(tree):
    """
    Convert a tree dictionary (as returned by parse) back to a Newick-format string.
    """

    def recurse(node):
        # If there are child nodes, enclose them in parentheses
        if node['children']:
            children_str = ",".join(recurse(child) for child in node['children'])
            name_part = node['name'] if node['name'] else ""
            length_part = f":{node['length']}" if node['length'] is not None else ""
            return f"({children_str}){name_part}{length_part}"
        else:
            # This is a leaf node
            name_part = node['name'] if node['name'] else ""
            length_part = f":{node['length']}" if node['length'] is not None else ""
            return f"{name_part}{length_part}"

    # Build the string for the entire tree and add a semicolon
    return recurse(tree) + ";"

def parse(newick):
    tokens = re.finditer(r"([^:;,()\s]*)(?:\s*:\s*([\d.]+)\s*)?([,);])|(\S)", newick+";")

    def recurse(nextid = 0, parentid = -1): # one node
        thisid = nextid#;
        children = []

        name, length, delim, ch = next(tokens).groups(0)
        if ch == "(":
            while ch in "(,":
                node, ch, nextid = recurse(nextid+1, thisid)
                children.append(node)
            name, length, delim, ch = next(tokens).groups(0)
        return {"id": thisid, "name": name, "length": float(length) if length else None,
                "parentid": parentid, "children": children, 'vector':None}, delim, nextid
    return recurse()[0]

# sample_vectors = None
def do_and(result, thr=0.5):
    res_vec = []
    n = len(list(result.values())[0])
    for i in range(n):
        res = 1
        for vec in result.values():
            if vec[i] == 0:
                res = 0
                break
        res_vec.append(res)
    return res_vec

def do_and_percent(result, thr=0.5):
    res_vec = []
    n = len(list(result.values())[0])
    for i in range(n):
        ones, zeros = 0, 0
        for vec in result.values():
            if vec[i] == 0:
                zeros += 1
            else:
                ones += 1
        if ones/(ones+zeros) > thr:
            res_vec.append(1)
        else:
            res_vec.append(0)
    return res_vec



def read_matrix(path):
    # Read the presence/absence file
    with open(path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    # First line is header: variant_id, sample1, sample2, ...
    header = lines[0].split('\t')
    sample_names = header[1:]  # Exclude 'variant_id'
    num_samples = len(sample_names)

    # Initialize a list (or dict) to hold each sample's presence/absence vector
    # Key: sample_name, Value: list of 0/1 for each variant
    samp_vectors = {sample: [] for sample in sample_names}

    # Read each subsequent line (one per variant)
    for line in lines[1:]:
        fields = line.split('\t')
        variant_id = fields[0]  # not needed for distance calculations, but present
        calls = fields[1:]  # the 0/1 calls for each sample

        # Store the presence/absence in each sample's vector
        for i, sample in enumerate(sample_names):
            samp_vectors[sample].append(int(calls[i]))
    return samp_vectors

def distance(v1, v2):
    dist = 0
    for i in range(len(v1)):
        if v1[i] != v2[i]:
            dist += 1
    return dist


def improve_tree(donkey_monkey, presence_matrix):
    if not donkey_monkey["children"]:
        #change to a global variable later so it will not run 6o times
        samp_vectors = read_matrix(presence_matrix)
        return samp_vectors[donkey_monkey["name"]]
    result = {}
    for child in donkey_monkey["children"]:
        result[child['id']] = improve_tree(child, presence_matrix)
    my_vector = do_and_percent(result)
    for child in donkey_monkey["children"]:
        child['length'] = distance(my_vector, result[child['id']])
        child['vector'] = my_vector
    return my_vector

def read_newick_file(file_path):
    with open(file_path, 'r') as file:
        newick_string = file.read().strip()
    return newick_string

def main(newick_tree_path, presence_matrix_path, improved_tree_path):
    tree = read_newick_file(newick_tree_path)
    tree = parse(tree)
    improve_tree(tree, presence_matrix_path)
    new_tree = to_newick(tree)
    print(new_tree)
    # Phylo.write(new_tree, improved_tree_path, "newick")

if __name__ == "__main__":
    main()