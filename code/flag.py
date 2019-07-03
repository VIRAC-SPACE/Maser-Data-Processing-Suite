import json

with open("results/cepa.json") as result_data:
    results = json.load(result_data)

for experiment in results:
    results[experiment]["flag"] = False

with open("results/cepa.json", "w") as result_data:
    result_data.write(json.dumps(results, indent=2))
