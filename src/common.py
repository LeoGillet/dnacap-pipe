# -*- coding: utf-8 -*-
"""
This module stores any miscellaneous function that is used throughout multiple modules
"""
import json

from src.logger import log


def get_tool_path(tool: str):
    """
    Reads config.json file to fetch tool paths
    :param tool: name of the tool - key in 'paths' object of config.json
    :return: string containing path of the tool
    """
    with open('config.json', 'r', encoding='UTF-8') as config_file:
        config = json.load(config_file)
        if tool not in config['paths']:
            log(f"Tool {tool} not found in config file. Stopping", level='ERROR')
            raise KeyError(f"Tool {tool} not found in config file.")
        return config['paths'][tool]
