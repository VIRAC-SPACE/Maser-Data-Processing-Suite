import json

filename = "results/cepa.json"

with open(filename) as result_data:
    results = json.load(result_data)

for experiment in results:
    results[experiment]["flag"] = False

with open(filename, "w") as result_data:
    result_data.write(json.dumps(results, indent=2))