import os
import sys

# ONLINE means file is only on disk
# NEARLINE mean file is only on tape
# ONLINE_and_NEARLINE means the file is available both ways

def run_cat_command(file_path):
    directory = os.path.dirname(file_path)
    file_name = os.path.basename(file_path)

    command = f"cat {directory}/\".(get)({file_name})(locality)\""

    os.system(command)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_file>")
        sys.exit(1)

    file_path = sys.argv[1]
    run_cat_command(file_path)
