"""
Given a backtrace file, with lots of ??'s, this will run addr2line on all
addresses and replace the question marks with the human-readable file names
and line numbers.
"""

import re
import sys
import subprocess

def convert_backtrace(filename):

    txt = ""

    with open(filename) as file:
        txt = file.read()

    new_text = ""

    a = re.compile(r"(\.\/.*\.ex)\((\+0x[0-9a-f]+)\).*$\n(^[\t ]*)\?\?\n[\t ]*\?\?:0", re.M)

    b = re.compile(r".*$", re.M)

    # if there are no matches, do nothing
    if a.search(txt) is None:
        return

    last_pos = 0

    for m in a.finditer(txt):
        # executable
        executable = m.group(1)
        # address to need to look up
        address = m.group(2)
        # indentation to use
        indentation = m.group(3)

        firstline = b.match(m.group(0)).group(0)

        bashCommand = f"addr2line -Cfie {executable} {address}"

        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()

        new_text += txt[last_pos:m.start()] + firstline + '\n' + output.decode("utf-8")
        last_pos = m.end()

    new_text += txt[last_pos:]
    print(new_text)

    with open(filename, 'w') as file:
        file.write(new_text)


if __name__ == "__main__":

    if len(sys.argv) < 2:
        sys.exit("Please provide a Backtrace file to read!")

    convert_backtrace(sys.argv[1])
