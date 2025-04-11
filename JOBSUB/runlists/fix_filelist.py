def filter_root_files(input_file, output_file):
    """
    Filters lines ending with '.root' from the input file and writes them to the output file.

    :param input_file: Path to the input text file.
    :param output_file: Path to the output text file.
    """
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                if line.strip().endswith('.root'):
                    outfile.write(line)
        print(f"Filtered lines written to {output_file}")
    except FileNotFoundError:
        print(f"Error: File {input_file} not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Example usage
if __name__ == "__main__":
    input_file = "29218_filepaths.txt"  # Replace with your input file path
    output_file = "29218_filepaths_fixed.txt"  # Replace with your desired output file path
    filter_root_files(input_file, output_file)
