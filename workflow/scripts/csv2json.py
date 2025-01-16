import csv
import json

# Read the CSV file and extract values
def csv_to_json(csv_file, json_file):
    with open(csv_file, mode='r', encoding='utf-8') as file:
        reader = csv.DictReader(file)
        json_data = {}
        
        for row in reader:
            attribute = row['Attribute']
            value = row['Value']
            json_data[attribute] = value
            

    # Write to JSON file
    with open(json_file, 'w', encoding='utf-8') as outfile:
        json.dump(json_data, outfile, ensure_ascii=False, indent=4)

# Usage
csv_file = "/home/shamlym/Downloads/Test case 7_ Metadata fields - 30-year mean_ diff-norheatwave.csv"  # Replace with your uploaded CSV file path
json_file = "/home/shamlym/workspace/klima-kverna/Klimakverna-Pilot1/config/testcase_7/30_year_mean_diff_norheatwave.json"
csv_to_json(csv_file, json_file)

print(f"JSON has been written to {json_file}")
