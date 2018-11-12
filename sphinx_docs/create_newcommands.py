import re

mathjax_config = {'TeX': {'Macros': {}}}

with open('mathsymbols.tex', 'r') as f:
    for line in f:
        macros = re.findall(r'\\newcommand{\\(.*?)}(\[(\d)\])?{(.+)}', line)
        for macro in macros:
            if len(macro[1]) == 0:
                mathjax_config['TeX']['Macros'][macro[0]] = "{"+macro[3]+"}"
            else:
                mathjax_config['TeX']['Macros'][macro[0]] = ["{"+macro[3]+"}", int(macro[2])]




