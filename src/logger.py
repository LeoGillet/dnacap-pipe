# -*- coding: utf-8 -*-
"""
This module contains all functions relating to the logging system implemented
and used in all other modules.
"""
from datetime import datetime

from colorama import init, Fore, Style

init()


def str_err(msg):
    """Formats string to appear like an error"""
    return Fore.RED + msg + Style.RESET_ALL


def str_success(msg):
    """Formats string to appear like a success"""
    return Fore.GREEN + msg + Style.RESET_ALL


def str_info(msg):
    """Formats string to appear like an information"""
    return Fore.CYAN + msg + Style.RESET_ALL


def str_warn(msg):
    """Formats string to appear like a warning"""
    return Fore.YELLOW + msg + Style.RESET_ALL


# ---------------------------------------------

def _fmt_log_line(info: str, lvl: str) -> str:
    """Formats string to add conventional pipe-separated log info"""
    now = datetime.now().strftime("%d/%m/%Y|%H:%M:%S.%f")
    return f"{now}|{lvl.upper()}|{info}\n"


def log(msg, level='INFO'):
    """Writes formatted line to log file"""
    with open("log/latest.log", 'a', encoding='UTF-8') as logfile:
        logfile.write(_fmt_log_line(msg, level))
