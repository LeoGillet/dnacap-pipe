# -*- coding: utf-8 -*-
from datetime import datetime

from colorama import init, Fore, Fore, Style

init()


def str_err(msg):
    return Fore.RED + msg + Style.RESET_ALL


def str_success(msg):
    return Fore.GREEN + msg + Style.RESET_ALL


def str_info(msg):
    return Fore.CYAN + msg + Style.RESET_ALL


def str_warn(msg):
    return Fore.YELLOW + msg + Style.RESET_ALL


# ---------------------------------------------

def _fmt_log_line(s: str, lvl: str) -> str:
    now = datetime.now().strftime("%d/%m/%Y|%H:%M:%S.%f")
    return f"{now}|{lvl.upper()}|{s}\n"


def log(msg, level='INFO'):
    with open("log/latest.log", 'a') as logfile:
        logfile.write(_fmt_log_line(msg, level))
