import csv  # comma-separated values to manipulate structured files (read & write)
import matplotlib.pyplot as plt  # Library for plotting figures
from collections import Counter  # Import Counter function from collections to count in a dictionary
# import pysam
import argparse  # Library that automatically creates a help message and manages arguments
import os

def script_call():
    # Create a parser for arguments
    parser = argparse.ArgumentParser(
        description="This script analyzes SAM files and provides information about the reads."
    )
    
    # Add arguments 
    # 1st Arg: used to specify the path of the SAM file to analyze
    parser.add_argument(
        '-i', '--input', 
        type=str,
        required=True, 
        help="Specify the path of the SAM file."
    )
    # 2nd Arg: used to specify the output file containing filtered reads
    parser.add_argument(
        '-o', '--output', 
        type=str, 
        required=True,
        help="Specify the output SAM file."
    )
    # 3rd Arg: used to specify the flag size
    parser.add_argument(
        '-fs', '--flag_size', 
        type=int,
        default=12, 
        required=False,
        help="Specify the size of FLAGS in bits."
    )
    
    # Parse all arguments
    args = parser.parse_args()
    
    sam_file_path = args.input
    if os.path.exists(sam_file_path):
        print("SAM file found")
    else:
        print("The path is incorrect: File not found")
        print("\n")
    
    output_sam_file = args.output
    flag_size = args.flag_size

    return sam_file_path, output_sam_file, flag_size


# Reading the SAM file
# Define a function named sam_reading that takes the SAM file path as an argument/parameter
def sam_reading(sam_file_path):
    # Open the SAM file in read mode "r"
    with open(sam_file_path, "r") as sam_file:
        # Initialize the CSV reader with tab delimiter as SAM files are tab-separated
        # sam_reader is an array (matrix)
        sam_reader = csv.reader(sam_file, delimiter='\t')
        # Replace with r split
        # Iterate through each line of the SAM file in a loop
        # Create (declare) empty lists and dictionaries to populate as the lines are processed
        flags = []  # list
        quals = []
        maps_score = []
        cigars = []
        references = []
        start_positions = []
        seq_lengths = []
        # For each line in the SAM file, iterate in a loop
        for row in sam_reader:
            # Ignore header lines that start with '@'
            if row[0].startswith("@"):
                continue
            
            # Access or extract information from each column
            qname = row[0]  # Read name
            flag = int(row[1])  # Flag as an integer
            flags.append(flag)  # Add each line's flag to the list
            rname = row[2]  # Reference sequence name
            references.append(rname)
            start_pos = int(row[3])  # Start position of alignment
            start_positions.append(start_pos)
            mapq = int(row[4])  # Mapping quality
            maps_score.append(mapq)  # Add each line's quality to the list
            cigar = row[5]  # CIGAR string (describes how a read aligns to the reference sequence)
            cigars.append(cigar)
            rnext = row[6]  # Reference of the next read in case of pairs
            pnext = int(row[7])  # Position of the next read in case of pairs
            tlen = int(row[8])  # Fragment length for pairs
            seq = row[9]  # DNA sequence
            seq_length = len(seq)
            seq_lengths.append(seq_length)
            qual = row[10]  # Base quality for each sequence
            quals.append(qual)  # Add each line's quality to the list

    return flags, quals, maps_score, cigars, references, start_positions, seq_lengths


#Question 1: How many reads are mapped? 
# Count the number of mapped reads by analyzing the SAM flag field, which encodes various attributes of a read.
# In sequencing, the SAM format is used to store alignment information, including whether a read is successfully aligned (mapped) or not.
# Define a function to convert SAM flag into binary representation, as flags are stored in a bit-packed integer
# Each bit in the flag represents a specific characteristic of the read (e.g., mapping status, strand orientation, etc.).
def flags_to_binary(flag_size, flags):
    # Loop through each flag in the input list
    for i in range(len(flags)):
       # Convert each flag into binary
        flags[i] = bin(flags[i])
        # Remove the first two characters (0 and b) & pad with zeros on the left to reach the desired size, which is 12
        flags[i] = flags[i][2:].zfill(flag_size)
        # Return the list of flags converted to binary
    return flags


# Take the binary flags and test the value of bit 2, which indicates if the read is mapped or not
# Define a function named number_of_mapped_reads that takes the flag size and binary flags as parameters
def number_of_mapped_reads(flag_size, binary_flags):
    # Declare an empty list to store mapped reads
    mapped_reads = []
    flag_size = flag_size - 1 # (0 --> 11)
    nbr = 0
    for i in range(len(binary_flags)):
        # Store the binary flag of the current line (i) in the variable flag
        flag = binary_flags[i]
        # Bit 2 indicates if the read is mapped (0 = mapped), so subtract 2 from the flag size (12-1=11)
        if (flag[flag_size-2] == "0"):
            # Add 1 to count the number of reads
            nbr = nbr + 1
            mapped_reads.append(i)
    # Calculate the % of mapped reads compared to the total number of reads
    nbr_prcnt = (nbr / len(binary_flags)) * 100
    print("First question:")
    print(f"There are {len(binary_flags)} reads.")
    print(f"There are {nbr} mapped reads, representing {nbr_prcnt:.2f}% of the reads.")
    print(f"There are {len(binary_flags)-nbr} unmapped reads, representing {100-nbr_prcnt:.2f}% of the reads.\n")
    # The percentage of mapped reads is a useful metric in sequencing and genome alignment studies
    # Mapped reads are crucial in sequencing studies as they indicate which parts of the genome are represented in the sequencing data.
    return mapped_reads, nbr


# Question 2: How are the reads and read pairs mapped? Count the number of reads for each flag
# For mapped reads, how many are partially/fully/incorrectly mapped?
# We analyze how well the reads are mapped to the reference genome, 
# and how the mapping status is encoded through CIGAR strings.
def fully_or_partially_mapped_reads(nbr, cigars, mapped_reads):
    mapped_reads = set(mapped_reads)
    fully_mapped = 0
    partially_mapped = 0
    incorrectly_mapped = 0

# Iterate over the CIGAR strings of all reads
    for i in range(len(cigars)):
        if i in mapped_reads:
            # A CIGAR string indicates the alignment of a read to the reference.
            # If a soft clipping ('S') or hard clipping ('H') is present, the read is partially mapped.
            if "S" in cigars[i] or "H" in cigars[i]:
                partially_mapped = partially_mapped + 1
            # If insertions ('I'), deletions ('D'), skips ('N'), mismatches ('X'), or padding ('P') are present,
            # the read is considered incorrectly mapped.    
            elif "I" in cigars[i] or "D" in cigars[i] or "N" in cigars[i] or "X" in cigars[i] or "P" in cigars[i]:
                incorrectly_mapped = incorrectly_mapped + 1
            else:
                # Otherwise, the read is considered fully mapped (perfect match with the reference).
                fully_mapped = fully_mapped + 1
    # Calculate the percentage of each mapping type
    fully_mapped_prcnt = (fully_mapped / nbr) * 100
    partially_mapped_prcnt = (partially_mapped / nbr) * 100
    incorrectly_mapped_prcnt = (incorrectly_mapped / nbr) * 100

    # Print results for mapping status 
    print("Second question:")
    print(f"There are {fully_mapped} fully mapped reads, representing {fully_mapped_prcnt:.2f}% of the mapped reads.")
    print(f"There are {partially_mapped} partially mapped reads, representing {partially_mapped_prcnt:.2f}% of the mapped reads.")
    print(f"There are {incorrectly_mapped} incorrectly mapped reads, representing {incorrectly_mapped_prcnt:.2f}% of the mapped reads.")
    print("\n")


# Calculate the number of paired reads and how many are mapped; 
# Test the value of bit 6 (indicating the first read of the pair) and bit 0 (indicating if it's part of a pair)
# Paired-end sequencing generates two reads per fragment, and it’s essential to track both parts of the pair for alignment quality.
def paired_reads_analysis(flag_size, binary_flags):
    # Adjust flag size to analyze specific bits related to read pairs
    flag_size = flag_size - 1
    # Initialize counters for the first and second read in each pair
    nbr_1 = 0 # First of the pair
    nbr_2 = 0 # Second of the pair
    # Iterate over all binary flags to identify paired reads
    for i in range(len(binary_flags)):
        flag = binary_flags[i]
        # Bit 6 indicates if it's the first read of the pair (1 = first read),
        # and bit 0 indicates whether the read is part of a pair (1 = part of a pair).
        if ((flag[flag_size-6] == "1") and flag[flag_size-0] == "1"):
            nbr_1 = nbr_1 + 1 # Increment count for the first read of the pair
        elif ((flag[flag_size-7] == "1") and flag[flag_size-0] == "1"):
            nbr_2 = nbr_2 + 1 # Increment count for the second read of the pair
    # Calculate the total number of paired reads
    nbr_de_paires = nbr_1 + nbr_2
    nbr_1_prcnt = (nbr_1 / nbr_de_paires) * 100 # Percentage of first reads
    nbr_2_prcnt = (nbr_2 / nbr_de_paires) * 100 # Percentage of second reads        
 
    # Print results for paired reads 
    print(f"There are {nbr_de_paires} paired reads.")
    print(f"There are {nbr_1} reads that are the first of the pair, representing {nbr_1_prcnt:.2f}%.")
    print(f"There are {nbr_2} reads that are the second of the pair, representing {nbr_2_prcnt:.2f}%.")
    print("\n")


# Question 3: Where are the reads mapped? Is the alignment uniform along the reference sequence? 
# Count the number of reads per chromosome. 
# In sequencing, it's important to know which chromosomes or regions of the genome the reads are aligning to, 
# and whether the depth is uniform or biased in certain areas.

def chromosome_distribution(references):
    # Initialize a dictionary to store the number of reads mapped to each chromosome (reference).
    chrom_counts = {}
    for i in range(len(references)):
        chrom = references[i]
        # Count the number of reads mapped to each chromosome.
        if chrom in chrom_counts:
            chrom_counts[chrom] += 1 # Increment the count for an already existing chromosome.
        else:
            chrom_counts[chrom] = 1 # Initialize the count for a new chromosome.
    print("3rd question:")
    print("Distribution by chromosome:")
    # Iterate over the chromosomes and display the number of reads mapped to each.
    for chrom in chrom_counts.keys():
        print(f"I have {chrom_counts[chrom]} reads on {chrom}")
    print("\n")

# In the next function, we calculate the read depth at specific positions along the chromosomes
# to check if the reads are uniformly distributed or if there are regions with higher or lower depth.
def read_positions(start_positions, seq_lengths, references):
    # Initialize a dictionary to store the depth (depth) for each chromosome and its positions.
    depth_by_ref = {}
    for i in range(len(references)):
        chrom = references[i]
        # If the chromosome is not already in the dictionary, initialize it with an empty dictionary.
        if chrom not in depth_by_ref:
            depth_by_ref[chrom] = {}
        # Iterate through the positions covered by the current read
        for pos in range(start_positions[i], start_positions[i] + seq_lengths[i]):
            # If the position exists in the second nested dictionary, increment its value by 1, otherwise initialize it to 1.
            if pos in depth_by_ref[chrom]:
                depth_by_ref[chrom][pos] += 1
            else:
                depth_by_ref[chrom][pos] = 1

    # Loop with 2 variables: reference and depth (how reads cover the reference).
    for ref, depth in depth_by_ref.items():
        # Sort positions on the chromosome to visualize depth in order.
        positions = sorted(depth.keys())
        # Extract the depth values for the positions.
        depth_values = list(depth.values())

        # Plot the depth along the reference sequence.
        plt.figure(figsize=(10, 5))  # Create a figure with a specified size (in inches).
        plt.plot(positions, depth_values, label=f"Depth on {ref}")  # Plot the depth for the current chromosome.
        plt.xlabel("Position on the reference sequence")  # Label the X-axis with the genomic positions.
        plt.ylabel("Number of reads (depth)")  # Label the Y-axis with the number of reads covering each position.
        plt.title("Read depth along the reference sequence")  # Title of the plot.
        plt.legend()  # Display the legend with the reference name.
        plt.show()  # Show the plot.

# Question 4: With which quality are the reads mapped? 
# Mapping quality is a critical metric in sequencing that indicates the confidence of a read's alignment to the reference genome.
# Higher mapping quality values reflect greater confidence, while lower values may indicate ambiguous or unreliable alignments.
def mapping_quality(maps_score):
    # The list 'maps_score' contains the mapping quality values for each read. 
    # Mapping quality scores help to evaluate the reliability of the alignment for each read. 
    # These scores are particularly useful in downstream analyses, such as variant calling, where accurate mapping is essential.
    
    # Use Python's Counter to count the number of reads for each mapping quality score.
    mapq_counts = Counter(maps_score) # Returns a dictionary with mapping quality scores as keys and counts as values.

    print("4th question:")
    # Visualize the distribution of mapping quality scores as a bar plot.
    plt.figure(figsize=(10, 5))  # Initialize a figure with the specified size.
    plt.bar(list(mapq_counts.keys()), list(mapq_counts.values()), label="Number of reads by quality value") 
    plt.xlabel("Quality")  # Label the X-axis to indicate quality scores.
    plt.ylabel("Number of reads by quality")  # Label the Y-axis to indicate the count of reads for each quality score.
    plt.title("Number of reads according to quality values")  # Title of the plot.
    plt.legend()  # Display the legend explaining the plot.
    plt.show()  # Show the bar plot.

    # Define a threshold of mapping quality at 30.
    # Reads with mapping quality below 30 are considered less reliable based on common bioinformatics standards.
    seuil_mapq = 30
    read_quality_lower_30 = []

    # Identify reads with mapping quality below the threshold (30).
    # Such reads may represent ambiguous alignments or reads mapped to multiple regions.
    for i in range(len(maps_score)):
        if maps_score[i] < seuil_mapq:
            read_quality_lower_30.append(i)

    return read_quality_lower_30

def filtred_reads(mapped_reads, read_quality_lower_30, sam_file_path, filtred_sam_file):
    # Filter out reads with low mapping quality (<30) and preserve the remaining high-confidence reads.
    # Filtering helps ensure that downstream analyses are based on reliable data, improving the robustness of results.

    # Combine the indices of reads that are mapped (fully or partially) and those with low quality.
    filtered_indices = set(mapped_reads).union(read_quality_lower_30)

    # Open the input SAM file and create a new file for the filtered reads.
    with open(sam_file_path, "r") as sam_file, open(filtred_sam_file, "w", newline='') as output_file:
        # Use a CSV reader to process the SAM file line-by-line (SAM files are tab-delimited).
        sam_reader = csv.reader(sam_file, delimiter='\t')
        output_writer = csv.writer(output_file, delimiter='\t')

        for index, row in enumerate(sam_reader):
            # Preserve header lines (starting with '@') in the output SAM file, as they contain metadata.
            if row[0].startswith("@"):
                output_writer.writerow(row)
                continue

            # Write the reads to the output file if their index matches the filtered criteria.
            if index in filtered_indices:
                output_writer.writerow(row)

# Biologically informed explanation:
# - High mapping quality reads are critical for accurate variant detection and downstream genomic analyses.
# - Reads with low mapping quality (e.g., <30) may be misaligned or map to multiple locations in the genome, 
#   making them less reliable for identifying true genomic variants.

# Function calls for processing the data and filtering based on mapping quality.
sam_file_path, sam_output_file, flag_size = script_call()
flags, quals, maps_score, cigars, references, start_positions, seq_lengths = sam_reading(sam_file_path)
# Convert the flag values into binary representation.
binary_flags = flags_to_binary(flag_size, flags)
# Determine the mapped reads and their count.
mapped_reads, nbr = number_of_mapped_reads(flag_size, binary_flags)
# Analyze fully, partially, and incorrectly mapped reads.
fully_or_partially_mapped_reads(nbr, cigars, mapped_reads)
# Analyze paired reads (first or second of the pair).
paired_reads_analysis(flag_size, binary_flags)
# Examine the distribution of reads across chromosomes.
chromosome_distribution(references)
# Analyze read depth along the reference genome.
read_positions(start_positions, seq_lengths, references)
# Evaluate the mapping quality and identify low-quality reads.
read_quality_greater_30 = mapping_quality(maps_score)
# Filter reads based on quality and save the filtered reads into a new SAM file.
filtred_reads(mapped_reads, read_quality_greater_30, sam_file_path, sam_output_file)
