import json
import os


class Status:
    """
    Used by pipeline helper to track status of workflow,
    which tasks have been done, and which have yet to be completed
    """

    def __init__(self, json_path="output/status.json"):
        self.json_path = json_path
        self.json_created = 0
        self.prerequisites = 0
        self.fastqc_reports = 0
        self.cleaning = 0
        self.mapping = 0
        self.flagstats = 0
        self.coverage = 0
        self.depth_txt = 0
        self.depth_plots = 0
        self.stat_seq_n = '?'
        self.stat_ref_gen_n = '?'
        self.get_or_create()
        self.set('json_created', 1)

    def to_dict(self) -> dict:
        """
        Converts object to dict
        """
        return {
            "json_created": self.json_created,
            "prerequisites": self.prerequisites,
            "fastqc_reports": self.fastqc_reports,
            "cleaning": self.cleaning,
            "mapping": self.mapping,
            "flagstats": self.flagstats,
            "coverage": self.coverage,
            "depth_txt": self.depth_txt,
            "depth_plots": self.depth_plots,
            "stat_seq_n": self.stat_seq_n,
            "stat_ref_gen_n": self.stat_ref_gen_n,
        }

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
        match checkpoint:
            case "json_created":
                self.json_created = value
            case "prerequisites":
                self.prerequisites = value
            case "fastqc_reports":
                self.fastqc_reports = value
            case "cleaning":
                self.cleaning = value
            case "mapping":
                self.mapping = value
            case "flagstats":
                self.flagstats = value
            case "coverage":
                self.coverage = value
            case "depth_txt":
                self.depth_txt = value
            case "depth_plots":
                self.depth_plots = value
            case "stat_seq_n":
                self.stat_seq_n = value
            case "stat_ref_gen_n":
                self.stat_ref_gen_n = value
        self.dump()

    def from_dict(self, data: dict):
        """
        Reads dict and updates checkpoint statuses
        """
        for checkpoint, value in data.items():
            if checkpoint in self.to_dict() and isinstance(value, int):
                self.set(checkpoint, value)

    def get_or_create(self):
        """
        Called at instanciation, find status JSON in directory and loads it or creates it
        """
        if not os.path.isfile(self.json_path):
            with open(self.json_path, "w", encoding="UTF-8") as json_file:
                print(self.to_dict())
                json.dump(self.to_dict(), json_file)
            return
        with open(self.json_path, "r", encoding="UTF-8") as json_file:
            data: dict = json.load(json_file)
            self.from_dict(data)
