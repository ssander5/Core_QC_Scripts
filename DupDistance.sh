#######################
# Written by Sheri Sanders
# Jun 21, 2024
# with the help of ChatGPT :)
########################

import re
from collections import defaultdict
import math

def parse_fastq(file_path):
    """Parse the FASTQ file and return a list of reads with their headers."""
    reads = []
    with open(file_path, 'r') as file:
        while True:
            header = file.readline().strip()
            if not header:
                break
            sequence = file.readline().strip()
            file.readline()  # Plus line
            quality = file.readline().strip()
            reads.append((header, sequence, quality))
    return reads

def find_duplicates(reads):
    """Find and return duplicate reads along with their headers."""
    read_dict = defaultdict(list)
    for header, sequence, quality in reads:
        read_dict[sequence].append((header, quality))
    
    duplicates = {seq: headers for seq, headers in read_dict.items() if len(headers) > 1}
    return duplicates

def extract_coordinates(header):
    """Extract and return the coordinates from the header."""
    match = re.search(r':(\d+):(\d+):(\d+) ', header)
    if match:
        zone = int(match.group(1))
        x_coord = int(match.group(2))
        y_coord = int(match.group(3))
        return (zone, x_coord, y_coord)
    return None

def calculate_distance(coord1, coord2):
    """Calculate the Euclidean distance between two coordinates."""
    _, x1, y1 = coord1
    _, x2, y2 = coord2
    return math.sqrt((x2 - x1)**2 + (y2 - y1)**2)

def report_proximity_duplicates(duplicates, proximity_threshold=100):
    """Report duplicate reads that are in physical proximity based on coordinates."""
    proximity_duplicates = []
    for sequence, headers in duplicates.items():
        coords = [extract_coordinates(header) for header, _ in headers]
        for i in range(len(coords)):
            for j in range(i + 1, len(coords)):
                if coords[i] and coords[j] and coords[i][0] == coords[j][0]:  # Same zone
                    distance = calculate_distance(coords[i], coords[j])
                    if distance <= proximity_threshold:
                        proximity_duplicates.append((sequence, headers, coords, distance))
                        break  # Once a close pair is found, no need to check further
    return proximity_duplicates

def main(file_path, proximity_threshold=100):
    reads = parse_fastq(file_path)
    duplicates = find_duplicates(reads)
    proximity_duplicates = report_proximity_duplicates(duplicates, proximity_threshold)
    
    # Print the results
    for sequence, headers, coords, distance in proximity_duplicates:
        print(f"Sequence: {sequence}")
        print(f"Distance: {distance}")
        for (header, quality), coord in zip(headers, coords):
            print(f"Header: {header}, Coordinates: {coord}, Quality: {quality}")
        print()

# Replace 'your_fastq_file.fastq' with the actual path to your FASTQ file
main('forward.fq', proximity_threshold=10000)
