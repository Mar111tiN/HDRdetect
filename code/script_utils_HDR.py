import os
from subprocess import PIPE, run
from io import StringIO
from datetime import datetime as dt
import pandas as pd


ansii_colors = {
    "magenta": "[1;35;2m",
    "green": "[1;9;2m",
    "red": "[1;31;1m",
    "cyan": "[1;36;1m",
    "gray": "[1;30;1m",
    "orange": "[1;33;1m",
    "black": "[0m",
}

colors = {
    "process": ansii_colors["green"],
    "time": ansii_colors["magenta"],
    "normal": ansii_colors["orange"],
    "warning": ansii_colors["red"],
    "success": ansii_colors["cyan"],
}


def show_output(text, color="normal", multi=False, time=True, **kwargs):
    """
    get colored output to the terminal
    """
    time = (
        f"\033{colors['time']}{dt.now().strftime('%H:%M:%S')}\033[0m : " if time else ""
    )
    proc = f"\033{colors['process']}Process {os.getpid()}\033[0m : " if multi else ""
    text = f"\033{colors[color]}{text}\033[0m"
    print(time + proc + text, **kwargs)


def show_command(command, list=False, time=False, multi=True, **kwargs):
    """
    prints the command line if debugging is active
    """
    time = (
        f"\033{colors['time']}{dt.now().strftime('%H:%M:%S')}\033[0m : " if time else ""
    )
    proc = f"\033[92mProcess {os.getpid()}\033[0m : " if multi else ""
    if list:
        command = f"\033[1m$ {' '.join(command)}\033[0m"
    else:
        command = f"\033[1m$ {command}\033[0m"
    print(time + proc + command, **kwargs)
    return


def run_cmd(cmd, **kwargs):
    show_command(cmd, **kwargs)
    exit = run(cmd, shell=True, check=True)
    return exit == 0


def cmd2df(cmd, show=False, **kwargs):
    '''
    wrapper for running shell commands directly into a dataframe
    optional output with show argument that passes kwargs to show_command
    '''

    if show:
        show_command(cmd, **kwargs)
    
    # certain parsings (eg. bam_csvs) fail because of erroneous quoting in Qual fields
    try:
        cmd_df = pd.read_csv(
            StringIO(run(cmd, stdout=PIPE, check=True, shell=True).stdout.decode("utf-8")),
            sep="\t"
        )
    except:
        show_output(f"Parsing error in cmd:<{cmd}> falling back to non-quoting", color="warning")
        cmd_df = pd.read_csv(
            StringIO(run(cmd, stdout=PIPE, check=True, shell=True).stdout.decode("utf-8")),
            sep="\t",
            quoting=3
        )
    return cmd_df
