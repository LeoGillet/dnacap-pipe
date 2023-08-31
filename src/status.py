"""
Contains all classes and functions used to track all steps
that were already executed during the run
Uses a JSON file <status.json> in the output directory to save information
"""

import json
import os


class Status:
    """
    Used by pipeline helper to track status of workflow,
    which tasks have been done, and which have yet to be completed
    """

    def __init__(self, json_path="output/status.json"):
        self.json_path = json_path
        self.statuses = {
            "json_created": 0,
            "prerequisites": 0,
            "fastqc_reports": 0,
            "cleaning": 0,
            "mapping": 0,
            "flagstats": 0,
            "coverage": 0,
            "depth_txt": 0,
            "depth_plots": 0,
            "genome_extracted": 0,
            "genes_extracted": 0,
            "stat_seq_n": "?",
            "stat_ref_gen_n": "?",
            "pipeline_type": "unset",
        }
        self.get_or_create()
        self.set("json_created", 1)

    def to_dict(self) -> dict:
        """
        Returns statuses dict
        """
        return self.statuses

    def dump(self):
        """
        Dumps object to JSON file by converting it to dict
        """
        with open(self.json_path, "w", encoding="UTF-8") as json_file:
            json.dump(self.to_dict(), json_file)

    def set(self, checkpoint, value):
        """
        Sets checkpoint status to provided value
        """
        self.statuses[checkpoint] = value
        self.dump()

    def get(self, checkpoint):
        """
        Gets value of checkpoint
        """
        try:
            return self.statuses[checkpoint]
        except KeyError as exc:
            raise KeyError(
                "Unknown checkpoint", checkpoint, ", impossible to get"
            ) from exc

    def from_dict(self, data: dict):
        """
        Reads dict and updates checkpoint statuses
        """
        for checkpoint, value in data.items():
            if checkpoint in self.to_dict():
                self.set(checkpoint, value)

    def get_or_create(self):
        """
        Called at instanciation, find status JSON in directory and loads it or creates it
        """
        if not os.path.isfile(self.json_path):
            if not os.path.isdir("output"):
                os.mkdir("output")
            choice = "unset"
            while choice not in ("pylori", "genitalium"):
                choice = input(
                    "Please set type of pipeline [pylori|genitalium] : "
                ).lower()
            self.statuses["pipeline_type"] = choice
            with open(self.json_path, "w", encoding="UTF-8") as json_file:
                json.dump(self.to_dict(), json_file)
            return
        with open(self.json_path, "r", encoding="UTF-8") as json_file:
            data: dict = json.load(json_file)
            self.from_dict(data)
