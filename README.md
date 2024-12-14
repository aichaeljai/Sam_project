# Sam_project

# SAM File Analysis Script

This script analyzes SAM files, providing insights into read mapping, quality scores, and reference alignment. It performs various tasks such as counting mapped reads, analyzing alignment coverage, and filtering low-quality reads.

## Features

- Count the number of mapped and unmapped reads.
- Distinguish between fully and partially mapped reads based on CIGAR strings.
- Analyze mapping quality and visualize results.
- Determine read distribution across chromosomes.
- Filter and export reads based on quality thresholds.
- Provide graphical representation of coverage and mapping quality.

## Requirements

The script requires the following Python libraries:

- `csv` (built-in): For reading and writing structured data.
- `matplotlib`: For generating plots.
- `argparse` (built-in): For parsing command-line arguments.
- `os` (built-in): For file system interactions.
- `collections.Counter` (built-in): For counting elements efficiently.

## Installation

To install the required libraries, run:

for linux/unix systems or environnements like conda:

```bash
pip install matplotlib
```
for windows systems:

```bash
py -m pip install matplotlib
```

## Usage

Run the script from the command line with the following options:

```bash
python script_name.py -i <input.sam> -o <output.sam> [-fs <flag_size>]
```

### Arguments

- `-i, --input`: **Required**. Path to the input SAM file.
- `-o, --output`: **Required**. Path to save the filtered output SAM file.
- `-fs, --flag_size`: Optional. Size of the flag in bits (default: 12).

### Example

```bash
python script_name.py -i sample.sam -o filtered_output.sam -fs 12
```

## What the Script Does

1. **SAM File Reading**:
   - Reads and processes the SAM file, extracting flags, mapping quality scores, CIGAR strings, and reference sequences.

2. **Mapped Reads Analysis**:
   - Counts the number of mapped reads and calculates the percentage of fully and partially mapped reads.
   - Identifies reads mapped to the complementary strand.

3. **Paired-End Analysis**:
   - Counts the number of paired reads and determines which are first or second in the pair.

4. **Chromosome Distribution**:
   - Analyzes read distribution across chromosomes and generates a coverage map.

5. **Mapping Quality Analysis**:
   - Visualizes mapping quality and filters low-quality reads (threshold: 30).

6. **Filtered Output**:
   - Outputs filtered reads to the specified file.

## Output

The script generates the following outputs:

- A filtered SAM file containing high-quality reads.
- Graphical plots for:
  - Mapping quality distribution.
  - Coverage across reference sequences.

## Limitations

- Only works with properly formatted SAM files.
- Requires sufficient memory for large files.

## Notes

Ensure the input SAM file path is correct and accessible. The script validates the input file and outputs an error message if the file is not found.

## License

This project is licensed under the MIT License.

---

For further questions or contributions, feel free to open an issue or a pull request.

