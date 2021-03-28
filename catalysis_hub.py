#!/bin/env python

"""
A script for querying the Catalysis-hub (https://catalysis-hub.org) database.
The data comprises heterogeneous catalytic reactions with atomic structures.
"""

import os
import json
import requests


# Define API root path and selected key values
# for the reactions
ROOT = 'http://api.catalysis-hub.org/graphql'
KEY_VALUES = [
    "chemicalComposition",
    "surfaceComposition",
    "facet",
    "sites",
    "coverages",
    "reactants",
    "products",
    "Equation",
    "reactionEnergy",
    "activationEnergy",
    "dftCode",
    "dftFunctional",
    "username",
    "pubId",
    "id",
]


def query_reactions(endcursor):
    """
    Run a batch query of 50 reactions. Searching the database in batches
    prevents the backend from overloading. The endcursor defines the starting
    position for each batch query. The query string is a structured string
    that contains all the key values that are included in the search.

    Parameters:
      endcursor (string):  Cursor position for the batch
    Returns:
      data (dict):         The reaction data for the current batch
    """

    # Define the query string
    # (for more info, see: http://docs.catalysis-hub.org/en/latest/)
    query_string = "{"
    query_string += 'reactions('
    query_string += 'first: 50,'
    query_string += f'after: "{endcursor}"'
    query_string += ', activationEnergy: 10, op: "<"'
    query_string += ") {"
    query_string += "  totalCount"
    query_string += "  pageInfo {"
    query_string += "    endCursor"
    query_string += "  }"
    query_string += "  edges {"
    query_string += "    node {"

    # Append the key values with correct structure
    for key_value in KEY_VALUES:
        query_string += str("\n" + " "*6 + key_value)

    query_string += "\n"
    query_string += "      systems {"
    query_string += "        id"
    query_string += "        Trajdata"
    query_string += "        energy"
    query_string += '        InputFile(format: "xyz")'
    query_string += "        keyValuePairs"
    query_string += "      }"
    query_string += "    }"
    query_string += "  }"
    query_string += "}}"

    # Perform the batch query
    data = requests.post(ROOT, {"query": query_string})
    try:
        data = data.json()["data"]
    except Exception as e:
        print(e)
        print("Error: Something went wrong.")
    return data


def parse_reaction(reaction):
    """
    Parse the reaction data so that all missing values are
    assigned a string value "None". The data for each reaction
    is saved in a dictionary.

    Parameters:
      reaction (dict):       Raw data describing a single reaction
    Returns:
      reaction_dict (dict):  Parsed reaction data
    """

    # "reaction_dict" has two keys: "key_value_pairs" and "structures"
    reaction_dict = {}

    # Record the reaction key values into a dictionary
    key_value_pairs = {}
    for key_value in KEY_VALUES:
        try:
            key_value_pairs[key_value] = reaction[key_value]
        except ValueError:
            key_value_pairs[key_value] = "None"
    # Some null values need to be parsed explicitly
    for key_value in KEY_VALUES:
        if not key_value_pairs[key_value]:
            key_value_pairs[key_value] = "None"
    reaction_dict["key_value_pairs"] = key_value_pairs

    # Record the structural data into a list of dictionaries
    structures = []
    for structure in reaction["systems"]:
        struct = {}
        struct["energy"] = structure["energy"]
        struct["InputFile"] = structure["InputFile"]
        struct["keyValuePairs"] = structure["keyValuePairs"]
        structures.append(struct)
    reaction_dict["structures"] = structures

    return reaction_dict


def main():
    """
    Reaction dictionaries are saved into a larger dictionary using
    the unique ids provided in the reaction data.
    """

    reaction_list = {}
    N_fetched = 0
    endcursor = ""
    N = 0
    totalcount = 100000

    # The database is queried in batches.
    # The variable totalcount is updated from the arbitrary
    # initial value to the actual number of found reactions.
    while N * 50 + 1 < totalcount:
        # Fetch a reaction data batch
        data = query_reactions(endcursor)

        # Record the reaction data into a dictionary
        for reaction in data["reactions"]["edges"]:
            reaction = reaction["node"]
            reaction_list[reaction["id"]] = parse_reaction(reaction)
        endcursor = data["reactions"]["pageInfo"]["endCursor"]
        totalcount = data["reactions"]["totalCount"]
        N_fetched += totalcount
        count = 50 * (N + 1)
        if count >= totalcount:
            count = totalcount
        print(f"Fetched reactions {50 * N + 1}-{count}/{totalcount}")
        N += 1
    print("Done!")

    cwd = os.getcwd()
    with open(f"{cwd}/data/reactions_cathub.json", "w") as outfile:
        json.dump(reaction_list, outfile)


if __name__ == '__main__':
    main()
