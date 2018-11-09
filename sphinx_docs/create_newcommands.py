import re

mathjax_config = """
MathJax.Hub.Config({{
    TeX: {{
        Macros: {{
{}
        }}
    }}
}});
"""

commands = {}

with open('mathsymbols.tex', 'r') as f:
    for line in f:
        macros = re.findall(r'\\newcommand{\\(.*?)}(\[(\d)\])?{(.+)}', line)
        for macro in macros:
            if len(macro[1]) == 0:
                commands[macro[0]] = "{"+macro[3]+"}"
            else:
                commands[macro[0]] = ["{"+macro[3]+"}", int(macro[2])]

cmd_str = ""
for k, v in commands.items():
    cmd_str += "            {}: {},\n".format(k, repr(v))

print(mathjax_config.format(cmd_str))


