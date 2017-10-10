#!/usr/bin/env python3
#-*- coding: utf-8 -*-

"""
"""

import argparse
import sys
try:
    from numpy import random
except ImportError:
    import random


def get_args():  # pragma: no cover
    parser = argparse.ArgumentParser(
        description="Given a seed, " +
        " return a pseudrando integer between 1 and 9999, separated by " +
        "newlines, to stdout.\n"
        "usage\n: " +
        "`seedRand.py 27 10` would return 10 random numbers seeded with 27")
    parser.add_argument("seed", action="store",
                        help="seed", type=int)
    parser.add_argument("n", action="store",
                        help="number of random numbers to return, " +
                        "must be > 0", type=int)
    args = parser.parse_args()
    return args


def main(args=None):
    """Generate a certain amount of pseudorandom numbers

    Args:
        args: argparse Namespace
    Returns:
        (str): 0 if all is well, 1 if error

    """
    try:
        random.seed(args.seed)
        for i in range(0, args.n):
            sys.stdout.write(str(random.randint(1, 9999)) + "\n")
        return 0
    except Exception as e:
        sys.stderr.write(e)
        return 1

if __name__ == "__main__":
    args = get_args()
    main(args)
