import json
import pandas as pd

# Load the TSV file
tsv_file = 'non_wgd_event_distributions.tsv'
df = pd.read_csv(tsv_file, sep='\t')

# Load the JSON file
json_file = '../configs/all_events_noWGD.json'
with open(json_file, 'r') as f:
    json_data = json.load(f)

# Create a mapping of event types to probabilities from the TSV file
prob_mapping = {}
for _, row in df.iterrows():
    string = row['loc']
    string = "".join([c.upper() if i == 0 else c for i, c in enumerate(string)])
    event_type = f"{string}{'Duplication' if row['type'] == 'gain' else 'Deletion'}"
    prob_mapping[event_type] = row['prob']

# Update the probabilities in the JSON file
for event in json_data['Signatures']['CNVs']['Events']:
    event_type = event['Type']
    if event_type in prob_mapping:
        event['Prob'] = prob_mapping[event_type]

# Save the updated JSON file
updated_json_file = 'updated_all_events_noWGD.json'
with open(updated_json_file, 'w') as f:
    json.dump(json_data, f, indent=4)

print(f"Updated JSON saved to {updated_json_file}")