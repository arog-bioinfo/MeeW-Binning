"""Thin source-cache-safe CLI boundary for the shared CheckM 1 runtime."""

from __future__ import print_function

import argparse
import importlib.util
import sys
from pathlib import Path


def main(arguments=None):
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--helper-source", required=True)
    parsed, remaining = parser.parse_known_args(arguments)
    source = Path(parsed.helper_source)
    spec = importlib.util.spec_from_file_location(
        "_meew_binning_checkm1_runtime", str(source)
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load CheckM 1 runtime helper from %s" % source)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    return helper.main(remaining)


if __name__ == "__main__":
    sys.exit(main())
