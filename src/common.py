# -*- coding: utf-8 -*-
import json

from src.logger import *


def get_tool_path(tool: str):
    with open('config.json', 'r') as config_file:
        config = json.load(config_file)
        if tool not in config['paths']:
            log(f"Tool {tool} not found in config file. Stopping", level='ERROR')
            raise KeyError(f"Tool {tool} not found in config file.")
        return config['paths'][tool]
