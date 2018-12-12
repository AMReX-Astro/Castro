#!/usr/bin/python

import re
import sys


def make_class_header(class_name, description):

    # remove // from description
    description = re.sub(r"//", "", description).strip()
    description = re.sub(r"\n", "\n///", description)
    class_name = re.sub(r"{", "", class_name).strip()
    class_name = class_name.split(':')[0].strip()

    boilerplate = f"""
///
/// @class {class_name}
/// @brief {description}
///"""

    return boilerplate


def make_method_header(description="", parameters=[]):
    # remove // from description
    description = re.sub(r"//", "", description).strip()
    description = re.sub(r"\n", "\n///", description)

    boilerplate = ""

    if description != "":

        boilerplate += f"""
///
/// {description}
///
"""
    elif parameters != []:
        boilerplate += """
///
"""

    if parameters != []:
        for param in parameters:
            boilerplate += f"""/// @param {(param.split('=')[0].strip()).split(' ')[-1]}
"""

        boilerplate += r"""///
"""

    return boilerplate


def make_method_doxycomment(description=""):
    # remove // from description
    description = re.sub(r"//", "", description).strip()
    description = re.sub(r"\n", "\n///", description)

    if description == "":
        return ""

    else:

        return f"""
///
/// @note
/// {description}
///
"""


def make_variable_docstring(description):
    description = re.sub(r"//", "", description).strip()
    description = re.sub(r"\n", "\n///", description)

    if description == "":
        return ""
    else:
        return f"""
///
/// {description}
///
"""


def process_header_file(filename):

    output_data = ""

    # find comments in lines above
    re_comments = re.compile(
        r"[ \t]*\/\/\s*\n[ \t]*(\/\/[ \t]*[\S ^\n]*?)\n[ \t]*\/\/")

    with open(filename) as input_file:
        data = input_file.read()

        # find classes
        re_class_name = re.compile(r"\n[^\n\s]*class ([^\n]+)")

        last_index = 0

        for m in re.finditer(re_class_name, data):

            comments = None
            for comments in re.finditer(re_comments,
                                        data[last_index:m.start()]):
                pass

            if comments and (m.start() - comments.end() - last_index) < 2:
                output_data += data[last_index:last_index + comments.start()]
                class_header = make_class_header(m.group(1), comments.group(1))
            else:
                output_data += data[last_index:m.start()]
                class_header = make_class_header(m.group(1), "")

            output_data += class_header

            last_index = m.start()

    output_data += data[last_index:]

    data = output_data
    output_data = ""
    last_index = 0

    re_prototype = re.compile(
        r"(?:^[\w&:*\t ]+\n)*^[ \t]*[~\w:*& <>]+\(([*\w\: \,&\n\t_=\<>\-.]*)\)", flags=re.MULTILINE)

    # markup methods
    for m in re.finditer(re_prototype, data):
        # print("match = ", m.group(1))

        parameters = m.group(1).split(",")
        parameters = [param.strip() for param in parameters]
        parameters = [param for param in parameters if param != ""]

        comments = None
        for comments in re.finditer(re_comments,
                                    data[last_index:m.start()]):
            pass

        if comments and (m.start() - comments.end() - last_index) < 2:
            # print(comments.span())
            output_data += data[last_index:last_index + comments.start()]
            method_header = make_method_header(comments.group(1), parameters)
            last_index = m.start()

        else:
            output_data += data[last_index:m.start()]
            method_header = make_method_header("", parameters)
            last_index = m.start()

        output_data += method_header

    output_data += data[last_index:]

    data = output_data
    output_data = ""
    last_index = 0

    re_comments = re.compile(
        r"^[ \t]*(\/\/[ \t]*[\S \n]*?)\n^(?![ \t]*\/\/)", flags=re.MULTILINE)

    re_variable = re.compile(
        r"^[ \t]*[~\w:*& <>\[\]]+;", flags=re.MULTILINE)

    # markup variables
    for m in re.finditer(re_variable, data):

        # print("match =", m.group(0))

        if " return " in m.group(0):
            continue

        comments = None
        for comments in re.finditer(re_comments,
                                    data[last_index:m.start() + 1]):
            pass

        # print(data[last_index:m.start()-1])

        # if comments:
        #     print(comments.group(0), m.start(), comments.end()+last_index)

        if comments and (m.start() - comments.end() - last_index) < 1:
            output_data += data[last_index:last_index + comments.start()]
            variable_header = make_variable_docstring(comments.group(1))
            output_data += variable_header
            last_index = m.start()

        else:
            output_data += data[last_index:m.start()]
            last_index = m.start()

    output_data += data[last_index:]

    output_filename = filename + ".doxygen"

    # print(output_data)

    with open(output_filename, 'w+') as output_file:
        output_file.write(output_data)


def process_cpp_file(filename):

    output_data = ""

    # find comments in lines above
    re_comments = re.compile(
        r"[ \t]*(\/\/[ \t]*[\S ^\n]*?)\n^[ \t]*[^\/]", flags=re.MULTILINE)

    re_prototype = re.compile(
        r"^\w*\n^\w[~\w:*& ]+\([\w\: \,&\n\t_=<>.]*\)\n?[\s\S]*?\n?{", flags=re.MULTILINE)

    with open(filename) as input_file:
        data = input_file.read()

        last_index = 0

        for m in re.finditer(re_prototype, data):

            comments = None

            for comments in re.finditer(re_comments,
                                        data[last_index:m.start() + 1]):
                pass

            if comments and (m.start() - comments.end() - last_index) < 3:
                output_data += data[last_index:last_index + comments.start()]
                method_header = make_method_doxycomment(comments.group(1))
                last_index = m.start()

            else:
                output_data += data[last_index:m.start()]
                method_header = make_method_doxycomment("")
                last_index = m.start()

            output_data += method_header

        output_data += data[last_index:]

    output_filename = filename + ".doxygen"

    with open(output_filename, 'w+') as output_file:
        output_file.write(output_data)


if __name__ == "__main__":
    filename = sys.argv[1]
    # print(filename[])

    if filename[-2:] == ".H":
        process_header_file(filename)
    elif filename[-4:] == ".cpp":
        process_cpp_file(filename)
