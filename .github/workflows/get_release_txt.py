#!/usr/bin/env python

"""
Get the text for the release from CHANGES.md
"""

import sys
from pathlib import Path

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("No version provided!")
    else:
        version = sys.argv[1]
        lines = Path("CHANGES.md").read_text().splitlines()

        release_lines = []
        in_release = False

        for line in lines:
            if line.startswith("## "):
                heading = line.removeprefix("## ").strip()

                if in_release:
                    break

                in_release = heading == version
                continue

            if in_release:
                release_lines.append(line)

        if not in_release:
            sys.exit(f"Could not find release notes for {version} in CHANGES.md")

        while release_lines and not release_lines[0].strip():
            release_lines.pop(0)

        while release_lines and not release_lines[-1].strip():
            release_lines.pop()

        print("\n".join(release_lines))
