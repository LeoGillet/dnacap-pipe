import os

import pandas as pd

from src import bed
from src import common


PYLORI_COLUMNS = [
    "16S_1 Mutations",
    "16S_2 Mutations",
    "23S_1 Mutations",
    "23S_2 Mutations",
    "gyrA Mutations",
    "rpoB Mutations",
    "cagA",
    "MLST best match",
    "MLST match 2",
    "MLST match 3",
    "MLST match 4",
    "MLST match 5",
]


def get_or_create(pipeline_type):
    if os.path.isfile("output/data.pickle"):
        return pd.read_pickle("output/data.pickle")
    if pipeline_type == "pylori":
        regions = []
        for name in bed._load_names(common.get_genome_regions(pipeline_type) + ".names.list"):
            regions.extend([name+" mean depth", name+" coverage"])
        columns = regions.extend(PYLORI_COLUMNS)
    return pd.DataFrame(
        index=[basename for (_, basename, _) in common.get_bam_files()],
        columns=columns,
    )


def save(dataframe: pd.DataFrame):
    dataframe.sort_index(inplace=True)
    dataframe.to_pickle("output/data.pickle")


def export_to_csv():
    if not os.path.isfile("output/data.pickle"):
        raise FileNotFoundError(
            "Pickled dataframe was not found under output/data.pickle"
        )
    df: pd.DataFrame = pd.read_pickle("output/data.pickle")
    df.to_csv("data.csv", sep=",")

def get_size_of_pickled_dataframe():
    if not os.path.isfile("output/data.pickle"):
        return "not found"
    return f"{os.path.getsize('output/data.pickle')/1000:.2f} KB"