#!/usr/bin/env python3
"""

@author: valengo
"""
import argparse
from typing import Tuple

import pandas
from enum import Enum

import sys

from cnvfinder.wrapper.argdesc import Strings


def create_parser(description: str, command: str = Strings.command.value, usage: str = None) -> argparse.ArgumentParser:
    if usage is None:
        usage = Strings.usage.value.format(command)
    return argparse.ArgumentParser(description=description, usage=usage,
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)


def parse_sub_command(parser: argparse.ArgumentParser) -> Tuple:
    return parser.parse_known_args(sys.argv[2:])


def get_arg_name_from_enum(arg: Enum):
    if type(arg.value) == dict:
        return '--' + list(arg.value.keys())[-1]
    return '--' + arg.name.replace('_', '-')


def get_arg_help_from_enum(arg: Enum):
    if type(arg.value) == dict:
        return list(arg.value.items())[-1][-1]
    return arg.value


def getattr_by(arg: Enum, args: Tuple):
    if type(args) == tuple:  # parse_known_args returns tuple
        args = args[0]

    name = arg.name
    if type(arg.value) == dict:
        name = list(arg.value.keys())[-1]
    try:
        return getattr(args, name.replace('-', '_'))
    except AttributeError:
        sys.exit(Strings.parsing_error.value.format(name, Strings.issue_url.value))


def load_table(filename: str):
    return pandas.read_table(filename)


def parse_legend(legends: list):
    parsed = []
    for legend in legends:
        parsed.append(dict(label=legend[0], color='#{}'.format(legend[-1])))
    return parsed
