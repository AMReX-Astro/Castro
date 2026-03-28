import argparse
import re

PR_URL_BASE = r"https://github.com/AMReX-Astro/Castro/pull/"

pr = re.compile(r"(\#)(\d+)")


def doit(clfile):

    with open(clfile) as cl:
        for line in cl:
            new_line = re.sub(pr, rf"[\g<0>]({PR_URL_BASE}\g<2>)", line)
            print(new_line, end="")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("changelog", type=str, nargs=1,
                        help="ChangeLog file")

    args = parser.parse_args()

    doit(args.changelog[0])
