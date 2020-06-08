#!/usr/bin/env python3

"""
Get the text for the release from CHANGES.md
"""

import re
import sys

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit('No version provided!')

    gen_version_re = re.compile(r'#\s(\d\d\.\d\d)')
    this_version_re = re.compile(f'#\s{sys.argv[1]}')

    with open('CHANGES.md', 'r') as file:
        txt = file.read()
        m = re.search(this_version_re, txt)
        if m:
            # find next date
            m_next = re.search(gen_version_re, txt[m.end():])
            if m_next:
                txt = '   ' + txt[m.end():m.end()+m_next.start()].strip()
            else:
                txt = '   ' + txt[m.end():].strip()
        else:
            txt = ""
                
        print(f'::set-env name=RELEASE_TXT::{txt}')