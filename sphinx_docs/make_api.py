import os
import re

rootdir = "../Source"

outfile_path = "source/filelist.rst"

with open(outfile_path, 'w') as outfile:

    output_data = """File list
=========

.. toctree::
   :maxdepth: 1

   """
    for subdir in sorted(os.listdir(rootdir)):
        if not os.path.isdir(os.path.join(rootdir, subdir)):
            continue

        output_data += f"""{subdir}_files
   """
        subdir_file_name = f"source/{subdir}_files.rst"

        with open(subdir_file_name, 'w') as subdir_file:

            subdir_output_data = f"{subdir.capitalize()}\n"
            subdir_output_data += "=" * len(subdir)
            subdir_output_data += """

.. toctree::
   :maxdepth: 2

   """

            for f in sorted(os.listdir(os.path.join(rootdir, subdir))):
                if (f[-4:] != ".cpp" and f[-4:] != ".F90" and f[-4:] != ".f90" and f[-2:] != ".H"):
                    continue

                rst_name = re.sub("_", "__", f)
                rst_name = re.sub("\.", "_8", rst_name)

                subdir_output_data += f"""file/{rst_name}
   """

            subdir_file.write(subdir_output_data)

    outfile.write(output_data)
