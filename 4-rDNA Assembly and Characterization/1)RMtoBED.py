import sys
import argparse
import os

def convert_rm_to_bed(rm_out, bed_out):
    """
    Parses a RepeatMasker .out file and writes a 6-column BED file.
    BED Format: chrom, chromStart(0-based), chromEnd, name, score, strand
    """
    converted_count = 0
    
    with open(rm_out, 'r') as infile, open(bed_out, 'w') as outfile:
        for line in infile:
            line = line.strip()
            
            # Skip empty lines
            if not line:
                continue
                
            # Skip the 3 standard header lines in a RepeatMasker .out file
            if "SW" in line and "perc" in line:
                continue
            if "score" in line and "div." in line:
                continue
            if not line[0].isdigit(): # Actual data lines always start with the SW score (a number)
                continue
            
            parts = line.split()
            
            # Ensure the line has enough columns to parse
            if len(parts) < 11:
                continue
            
            # Extract relevant columns
            # RM format: 0:SW_score, 1:perc_div, 2:perc_del, 3:perc_ins, 4:query_sequence, 
            # 5:query_start, 6:query_end, 7:(left), 8:strand, 9:repeat_name, 10:repeat_class
            
            chrom = parts[4]
            # Convert 1-based RepeatMasker start to 0-based BED start
            start = int(parts[5]) - 1 
            end = parts[6]
            
            # RepeatMasker uses 'C' or '+' for strand. BED uses '-' or '+'
            strand_char = parts[8]
            if strand_char.upper() == 'C':
                strand = '-'
            elif strand_char == '+':
                strand = '+'
            else:
                strand = '.' # Fallback if undefined
            
            repeat_name = parts[9]
            sw_score = parts[0] # Using Smith-Waterman score as the BED score column
            
            # Write out standard 6-column BED format
            outfile.write(f"{chrom}\t{start}\t{end}\t{repeat_name}\t{sw_score}\t{strand}\n")
            converted_count += 1

    return converted_count

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert RepeatMasker .out to .bed for IGV")
    parser.add_argument("-i", "--input", required=True, help="Input RepeatMasker .out file")
    parser.add_argument("-o", "--output", required=True, help="Output .bed file")
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' not found.")
        sys.exit(1)
        
    print(f"Converting '{args.input}' to BED format...")
    count = convert_rm_to_bed(args.input, args.output)
    print(f"Success! {count} repeat features written to '{args.output}'.")
