# -*- coding: utf-8 -*-
from src import common


def _sickle_command(pe_file1, pe_file2):
    return f"""{common.get_tool_path('sickle')} pe -f input/{pe_file1} -r input/{pe_file2} -o output/intermediate_files/sickled_{pe_file1} -p output/intermediate_files/sickled_{pe_file2} -s output/intermediate_files/sickled_singles.fastq"""
