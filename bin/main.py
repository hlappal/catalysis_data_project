#!/usr/bin/env python3

import os
import sys
from catalysis_hub import catalysis_hub
from catapp import catapp
from preprocess import preprocess, construct_struct_descriptors


def main(arg):
    """
    """

    # Get the root directory
    os.chdir("../")
    root_dir = os.getcwd()

    # Fetch and save the data
    if arg in ["-r", "--re-fetch"]:
        # Re-fetch the data even if it exists
        fetch_data("Catalysis-hub", root_dir, re_fetch=True)
        fetch_data("CatApp", root_dir, re_fetch=True)
    else:
        # Check if Catalysis-hub data exists
        if not os.path.isfile(f"{root_dir}/data/reactions_cathub.json"):
            # Fetch the data if it does not exist
            fetch_data("Catalysis-hub", root_dir, re_fetch=False)
        else:
            print("Catalysis-hub data already exists at data/reactions_cathub.json")
        # Check if CatApp data exists
        if not os.path.isfile(f"{root_dir}/data/reactions_catapp.json"):
            # Fetch the data if it does not exist
            fetch_data("CatApp", root_dir, re_fetch=False)
        else:
            print("CatApp data already exists at data/reactions_catapp.json")

    # Load and preprocess the data
    print("Processing data...")
    df = preprocess(root_dir)
    print(df.head())

    # Construct structural descriptors
    construct_struct_descriptors(root_dir)


def fetch_data(database, root_dir, re_fetch=False):
    """ """
    # Define the fetch message
    if re_fetch:
        msg = "Re-fetching "
    else:
        msg = "Fetching "
    msg += "data from " + database + " database..."

    print(msg)
    if database == "Catalysis-hub":
        catalysis_hub(root_dir)
    elif database == "CatApp":
        catapp(root_dir)
    else:
        print(f"Database error: No such database as {database}")


def print_help():
    print("Usage: python main.py [options]")
    print("\nOptions:")
    print(f"{'  -r, --re-fetch':20} Re-fetch the data from both databases.")
    print(f"{'  -h, --help':20} Print this help.")


if __name__ == "__main__":
    if len(sys.argv) > 2:
        print("Usage: python main.py [options]")
    else:
        if len(sys.argv) == 2:
            OPT = sys.argv[1]
        else:
            OPT = None
    # Parse options
    if OPT in ["-h", "--help"]:
        print_help()
    elif OPT in ["-r", "--re-fetch"]:
        main(OPT)
    elif OPT is None:
        main(OPT)
    else:
        print(f"Invalid option '{OPT}'\n")
        print_help()
