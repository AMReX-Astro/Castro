#!/usr/bin/env python

import argparse
import sys


def doit():

    parser = argparse.ArgumentParser()
    parser.add_argument("network", type=str, nargs=1)
    parser.add_argument("auto_source_dir", type=str, nargs=1)

    args = parser.parse_args()

    net = args.network[0]
    auto_dir = args.auto_source_dir[0]

    try:
        with open(f"{auto_dir}/NETWORK_USED") as f:
            stored_net = f.readlines()[0].strip()
            if stored_net != net:
                sys.exit("network inconsistent.  do 'make clean' and retry")

    except OSError:
        pass


if __name__ == "__main__":
    doit()
