#!/usr/bin/env python

import errno
import git
import pathlib
import os

from pathlib import Path
from typing import Any, Dict
from yaml import YAMLError, safe_load


def check_file(
    filename: str
) -> str:
    """Checks if file exists

    Args:
        filename (str): Complete path to file
    Raises:
        FileNotFoundError: If file does not exist
    Returns:
        str: Original file name
    """

    if not Path(filename).is_file():
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), filename
        )

    return filename


def check_path(
    filename: pathlib.PosixPath
) -> str:
    """Checks the path of a given file and recursively creates new directories
       if necessary
    Args:
        filename (pathlib.PosixPath): Complete path to file
    Returns:
        str: Original file name
    """

    if not Path(filename).parents[0].is_dir():
        filename.parents[0].mkdir(parents=True)

    return str(filename)


def get_path() -> str:
    """Get project path

    Returns:
        str: Full path to the project
    """

    git_repo = git.Repo(search_parent_directories=True)

    return str(Path(git_repo.working_tree_dir))


def read_yaml(
    config: str
) -> Dict[str, Any]:
    """Read .yaml file into dict

    Args:
        config (str): .yaml configuration file

    Returns:
        Dict[str, Any]: .yaml file content
    """

    with open(check_file(config), "r") as stream:
        try:
            data = safe_load(stream)
        except YAMLError as exc:
            print(exc)

    return data
