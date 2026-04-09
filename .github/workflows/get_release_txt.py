#!/usr/bin/env python

"""
Get the text for the release from CHANGES.md
"""

import re
import sys

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print('No version provided!')
    else:
        gen_version_re = re.compile(r'#\s(\d\d\.\d\d)')
        this_version_re = re.compile(fr'#\s{sys.argv[1]}')

        with open('CHANGES.md') as file:
            txt = file.read()
            m = re.search(this_version_re, txt)
            if m:
                # find next date
                m_next = re.search(gen_version_re, txt[m.end():])
                if m_next:
                    txt = txt[m.end():m.end()+m_next.start()].strip()
                else:
                    txt = txt[m.end():].strip()
            else:
                txt = ""

            # we now need to substitute characters in the string so that
            # the action can deal with line breaks
            txt = txt.replace('%', '%25')
            txt = txt.replace('\n', '%0A')
            txt = txt.replace('\r', '%0D')
            txt = txt.replace('%0A   *', '%0A*')

            print(f'"RELEASE_TXT=${{{txt}}}"')