import pysam
import pandas as pd
import numpy as np
from collections import defaultdict
from sklearn.preprocessing import LabelEncoder
from sklearn.ensemble import RandomForestClassifier
from sklearn.cluster import AgglomerativeClustering
from river import tree, preprocessing
import gzip
import argparse
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Load VCF file and extract SNP information
def load_vcf(vcf_file):
    snps = {}
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            snps[(chrom, pos)] = (ref, alt)
    logging.info(f"Loaded {len(snps)} SNPs from VCF file.")
    return snps

# Load barcode list
def load_barcodes(barcode_file):
    barcodes = []
    with gzip.open(barcode_file, 'rt') as f:
        for line in f:
            barcodes.append(line.strip())
    logging.info(f"Loaded {len(barcodes)} barcodes.")
    return barcodes

# Extract reads overlapping SNPs from BAM file in chunks
def extract_reads(bam_file, snps, chunk_size=1000000):
    reads = defaultdict(list)
    bam = pysam.AlignmentFile(bam_file, "rb")
    total_reads = 0
    
    for chrom, pos in snps.keys():
        for read in bam.fetch(chrom, pos-1, pos):
            if not read.is_unmapped and not read.is_duplicate:
                reads[read.query_name].append((chrom, pos, read.query_sequence, read.get_tag('CB') if read.has_tag('CB') else ''))
                total_reads += 1
                if total_reads % chunk_size == 0:
                    yield reads
                    reads = defaultdict(list)
    
    bam.close()
    if reads:
        yield reads

# Create feature matrix for machine learning
def create_feature_matrix(reads, snps, barcodes):
    data = []
    labels = []
    
    for read_name, read_data in reads.items():
        read_snps = {pos: seq for chrom, pos, seq, barcode in read_data}
        if any(barcode in barcodes for chrom, pos, seq, barcode in read_data):
            features = [read_snps.get(pos, 'N') for chrom, pos in snps.keys()]
            data.append(features)
            labels.append(read_data[0][3])  # Use barcode as label
    
    logging.info(f"Created feature matrix with {len(data)} samples.")
    return pd.DataFrame(data), labels

# Initialize model
model = tree.HoeffdingTreeClassifier()
scaler = preprocessing.StandardScaler()

# Process BAM chunk
def process_bam_chunk(reads_chunk, snps, barcodes, n_clusters):
    try:
        # Create feature matrix
        data, labels = create_feature_matrix(reads_chunk, snps, barcodes)
        
        # Encode features and labels for sklearn
        encoded_data, encoded_labels, label_encoder = encode_data(data, labels)
        
        # Initialize RandomForestClassifier
        rf_classifier = RandomForestClassifier(n_estimators=100, random_state=42)
        
        # Train the classifier
        rf_classifier.fit(encoded_data, encoded_labels)
        
        # Cluster the barcodes
        clusters = AgglomerativeClustering(n_clusters=n_clusters).fit_predict(encoded_data)
        
        # Annotate data
        annotated_data = annotate_data(data, clusters)
        
        # Update the model and predict
        results = []
        for x, y in zip(data.values, labels):
            x = scaler.learn_one(x)
            model.learn_one(x, y)
            prediction = model.predict_one(x)
            
            # Store the prediction in the results
            results.append({
                'features': x,
                'label': y,
                'prediction': prediction
            })
        
        return results
    except Exception as e:
        logging.error(f"Error processing BAM chunk: {e}")
        return []

# Encode features and labels for sklearn
def encode_data(data, labels):
    label_encoder = LabelEncoder()
    encoded_labels = label_encoder.fit_transform(labels)
    
    encoded_data = pd.get_dummies(data)
    
    return encoded_data, encoded_labels, label_encoder

# Annotate data with classification results
def annotate_data(data, clusters):
    data['cluster'] = clusters
    data['classification'] = 'unassigned'
    singlets = [idx for idx, count in enumerate(np.bincount(clusters)) if count == 1]
    doublets = [idx for idx, count in enumerate(np.bincount(clusters)) if count == 2]
    
    data.loc[data['cluster'].isin(singlets), 'classification'] = 'singlet'
    data.loc[data['cluster'].isin(doublets), 'classification'] = 'doublet'
    return data

# Process incoming messages
def process_messages(vcf_file, barcode_file, bam_file, n_clusters, output_file):
    snps = load_vcf(vcf_file)
    barcodes = load_barcodes(barcode_file)

    # Initialize empty DataFrame for final output
    final_output = pd.DataFrame()

    # Process the BAM file in chunks
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = []
        for reads_chunk in extract_reads(bam_file, snps):
            futures.append(executor.submit(process_bam_chunk, reads_chunk, snps, barcodes, n_clusters))
        
        for future in as_completed(futures):
            results = future.result()
            if results:
                final_output = final_output.append(results, ignore_index=True)
                logging.info(f"Appended results to final output. Current size: {len(final_output)}")

    # Save final results to CSV
    final_output.to_csv(output_file, index=False)
    logging.info(f"Final results saved to {output_file}")

# Argument parser
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Real-Time Demultiplexing of Pooled scRNA-seq or scATAC-seq Data")
    parser.add_argument("--vcf_file", required=True, help="Path to the VCF file")
    parser.add_argument("--barcode_file", required=True, help="Path to the gzipped barcode file")
    parser.add_argument("--bam_file", required=True, help="Path to the BAM file")
    parser.add_argument("--n_clusters", type=int, required=True, help="Number of expected donors/clusters")
    parser.add_argument("--output_file", required=True, help="Path to the output CSV file")

    args = parser.parse_args()

    process_messages(args.vcf_file, args.barcode_file, args.bam_file, args.n_clusters, args.output_file)
