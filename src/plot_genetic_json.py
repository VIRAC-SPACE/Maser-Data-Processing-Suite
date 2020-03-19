import sys
import os
import json


def get_all_json_file():
    path = "."
    return [file for file in os.listdir(path) if ".json" in file]


def main():
    for file in get_all_json_file():
        with open(file) as data:
            genetic_results = json.load(data)
    sys.exit(0)


if __name__ == "__main__":
    main()