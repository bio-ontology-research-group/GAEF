from pythoncyc import PGDB
from pythoncyc.PToolsFrame import PFrame
import pythoncyc.config
import re
import pandas as pd
import sys

# Set the hostname and port for the Pathway Tools API server
pythoncyc.config.set_host_name('localhost')
pythoncyc.config.set_host_port(5008)
go_version = sys.argv[1]

# Initialize the PGDB object for the 'meta' organism
meta = PGDB('meta')

# Define the pathway ID
pathway_id = 'PLPSAL-PWY'

# Retrieve the pathway PFrame object
pathway = PFrame(pathway_id, meta, getFrameData=True)

# Function to get reactions in a pathway, including handling subpathways
def get_pathway_reactions(pathway):
    reactions = set()
    if hasattr(pathway, 'reaction_list') and pathway.reaction_list is not None:
        reactions.update(pathway.reaction_list)
    if hasattr(pathway, 'predecessors') and pathway.predecessors is not None:
        for pred in pathway.predecessors:
            if isinstance(pred, PFrame) and pred.classname == 'Pathway':
                subpathway = PFrame(pred.frameid, meta, getFrameData=True)
                reactions.update(get_pathway_reactions(subpathway))
            else:
                reactions.add(pred)
    return reactions

# Expand subpathways and get all reactions
def expand_subpathways(pathway_id):
    reactions = set()
    pathway = PFrame(pathway_id, meta, getFrameData=True)
    if hasattr(pathway, 'reaction_list'):
        for reaction_id in pathway.reaction_list:
            reaction = PFrame(reaction_id, meta, getFrameData=True)
            if 'PWY' in reaction_id:
                reactions.update(expand_subpathways(reaction_id))
            else:
                reactions.add(reaction_id)
    return reactions

# Get all reactions in the pathway, including those in subpathways
reactions_in_pathway = expand_subpathways(pathway_id)

# Function to get predecessors of a reaction
def get_reaction_predecessors(reaction_id, pathway):
    try:
        predecessors = meta.sendPgdbFnCallList('get-predecessors', reaction_id, pathway.frameid)
        return [PFrame(pred, meta, getFrameData=True) for pred in predecessors]
    except Exception as e:
        print(f"Error getting predecessors for {reaction_id}: {e}")
        return []

# Convert the predecessors data into a dictionary for easy lookup
predecessors_dict = {}
for reaction_id in reactions_in_pathway:
    predecessors = get_reaction_predecessors(reaction_id, pathway)
    predecessors_dict[reaction_id] = [pred.frameid for pred in predecessors]

# Function to find all possible paths with cycle detection
def find_all_paths(predecessors_dict):
    def find_paths(reaction, visited):
        if reaction in visited:
            return []  # Avoid cycles
        if reaction not in predecessors_dict or not predecessors_dict[reaction]:
            return [[reaction]]
        paths = []
        visited.add(reaction)
        for pred in predecessors_dict[reaction]:
            for path in find_paths(pred, visited):
                paths.append(path + [reaction])
        visited.remove(reaction)
        return paths

    # Find all starting reactions (reactions with no predecessors)
    all_predecessors = {pred for preds in predecessors_dict.values() for pred in preds}
    starting_reactions = set(predecessors_dict.keys()) - all_predecessors

    # If no starting reactions are found, we have a cycle or disjoint graph
    if not starting_reactions:
        starting_reactions = set(predecessors_dict.keys())

    # Find all paths from each starting reaction
    all_paths = []
    for start_reaction in starting_reactions:
        all_paths.extend(find_paths(start_reaction, set()))
    return all_paths

# Get all possible paths
all_paths = find_all_paths(predecessors_dict)

# #Print all possible paths with EC numbers
# print("\nAll Paths with EC Numbers:")
# for path in all_paths:
#     ec_path = []
#     for reaction_id in path:
#         reaction = PFrame(reaction_id, meta, getFrameData=True)
#         if hasattr(reaction, 'ec_number') and reaction.ec_number:
#             ec_numbers = reaction.ec_number
#         else:
#             ec_numbers = ['No EC number']
#         ec_path.append(f"{reaction_id} (EC: {', '.join(ec_numbers)})")
#     print(" -> ".join(ec_path))



# Extract EC numbers and format them
def format_ec_numbers(ec_numbers):
    formatted_ecs = []
    for ec in ec_numbers:
        # Updated regex pattern to match EC numbers with digits or letters in each part
        matches = re.findall(r'\|EC-([A-Za-z0-9]+(?:\.[A-Za-z0-9]+){0,3})\|', ec)
        formatted_ecs.extend(matches)
    return formatted_ecs

# Function to generate all combinations of paths for multiple EC numbers
def generate_ec_combinations(paths_with_ecs):
    if not paths_with_ecs:
        return [[]]
    first, rest = paths_with_ecs[0], paths_with_ecs[1:]
    combinations = []
    for ec in first:
        for combination in generate_ec_combinations(rest):
            combinations.append([ec] + combination)
    return combinations

def process_pathway(pathway_id):
    try:
        reactions_in_pathway = expand_subpathways(pathway_id)
        predecessors_dict = {}
        for reaction_id in reactions_in_pathway:
            predecessors = get_reaction_predecessors(reaction_id, PFrame(pathway_id, meta, getFrameData=True))
            predecessors_dict[reaction_id] = [pred.frameid for pred in predecessors]
        all_paths = find_all_paths(predecessors_dict)
        paths_with_ecs = []
        for path in all_paths:
            path_ecs = []
            for reaction_id in path:
                reaction = PFrame(reaction_id, meta, getFrameData=True)
                if hasattr(reaction, 'ec_number') and reaction.ec_number:
                    ec_numbers = reaction.ec_number
                    formatted_ecs = format_ec_numbers(ec_numbers)
                    path_ecs.append(formatted_ecs)
            paths_with_ecs.append(path_ecs)
        all_ec_paths = []
        for path_ecs in paths_with_ecs:
            all_ec_paths.extend(generate_ec_combinations(path_ecs))
        return all_ec_paths
    except Exception as e:
        print(f"Error processing pathway {pathway_id}: {e}")
        return None  # Return None to indicate an error occurred

# Read the metacyc_GO.tsv file
df = pd.read_csv(f'metacyc_GO_v{go_version}.tsv', sep='\t')

# Prepare to collect results and skipped pathways
results = []
skipped_pathways = []

# Process each pathway
for index, row in df.iterrows():
    pathway_id = row['MetaCyc_Pathway']
    ec_combinations = process_pathway(pathway_id)
    if ec_combinations is not None:
        results.append({
            'GO_Term': row['GO_Term'],
            'MetaCyc_Pathway': pathway_id,
            'EC_Combinations': ec_combinations
        })
    else:
        skipped_pathways.append(pathway_id)

# Convert results to DataFrame and save to new file
results_df = pd.DataFrame(results)
results_df.to_csv(f'metacyc_GO_v{go_version}_with_EC.tsv', sep='\t', index=False)

# Save skipped pathways to a file
with open('skipped.txt', 'w') as f:
    for pathway in skipped_pathways:
        f.write(pathway + '\n')
