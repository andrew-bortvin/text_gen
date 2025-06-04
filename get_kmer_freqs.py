import re
import sys
import random
import pandas as pd
from scipy.stats import zscore
import numpy as np

k = 4
author = "willy"
target_length = 100
scaling_factor = 4


if len(sys.argv) > 1:
	k = int(sys.argv[1])
if len(sys.argv) > 2:
	author = sys.argv[2]
if len(sys.argv) > 3:
	target_length = int(sys.argv[3])
if len(sys.argv) > 4:
	scaling_factor = int(sys.argv[4])

###########################################
#Functions for constructing k-mer freq table#
###########################################


def fill_kmer_dict_construct_kf(in_string, k, in_dict, in_list, n_authors, current_author):
	# Args:
	# in_string: str, the text to be k-mer counted
	# k: int, length of k-mer
	# in_dict: dict, dictionary of k-mer counts across all authors
	# in_list: list, list of total k-mer count across all authors
	# n_authors: int, number of authors in the dataset
	# current_author: int, position in the list of authors 
	### Note that this is ZERO INDEXED
	for i in range(len(in_string) - k + 1):
		# Get k-mer
		current_kmer = in_string[i:i+k]
		# Default is a list, with one entry per author
		in_dict.setdefault(current_kmer, [0] * n_authors)
		# iterate the current author's k-mer count by 1
		in_dict[current_kmer][current_author] += 1
		# iterate the current author's total k-mer number by 1 
		in_list[current_author] += 1

	return (in_dict, in_list)

def parse_text_construct_kf(in_dict, k, in_file, in_list, n_authors, current_author):
	#Args:
	# in_dict: dict, the dictionary to be Modified
	# in_list: list, list of total k-mer count across all authors
	# k: int, length of k-mer
	# in_file: an opened file stream; the current book to be read
	# n_authors: int, number of authors in the dataset
	# current_author: the index of the current author to be read. Determines which
	### element of in_dict is modified
	### Note that this is ZERO INDEXED

	# Skip opening text: 
	l = in_file.readline()
	while not l.startswith("*** START OF THE PROJECT GUTENBERG EBOOK"):
		l = in_file.readline()

	text = ""
	for l in in_file:
		if l.startswith("*** END OF THE PROJECT GUTENBERG EBOOK"):
			break
		# Remove empty lines
		if l != "\n":
			# Convert newlines and tabs to spaces
			l = re.sub('\n', ' ', l)
			l = re.sub('\t', ' ', l)
			# Remove non alphanumerical characters
			l = re.sub(r'[^a-zA-Z\ ]', '', l)
			# Remove consecutive spaces
			l = re.sub(' +', ' ', l)
			# Remove case
			l = l.upper()

			text += l

	in_dict, in_list = fill_kmer_dict_construct_kf(text, k, in_dict, in_list, n_authors, current_author)

	return in_dict, in_list

###########################################
#Functions for constructing NEXT k-mer table#
###########################################

def fill_next_kmer_dict(in_string, k, in_dict):
	for i in range(len(in_string) - k - k + 1):
		current_kmer = in_string[i:i+k]
		next_kmer = in_string[i + k:i + k + k]
		in_dict.setdefault(current_kmer, dict())
		in_dict[current_kmer].setdefault(next_kmer, 0)
		in_dict[current_kmer][next_kmer] += 1

	return(in_dict)

def parse_text(in_dict, in_file):
	# Skip opening text: 
	l = in_file.readline()
	while not l.startswith("*** START OF THE PROJECT GUTENBERG EBOOK"):
		l = in_file.readline()

	text = ""
	for l in in_file:
		if l.startswith("*** END OF THE PROJECT GUTENBERG EBOOK"):
			break
		# Remove empty lines
		if l != "\n":
			# Convert newlines and tabs to spaces
			l = re.sub('\n', ' ', l)
			l = re.sub('\t', ' ', l)
			# Remove non alphanumerical characters
			l = re.sub(r'[^a-zA-Z\'\ ]', '', l)
			# Remove consecutive spaces
			l = re.sub(' +', ' ', l)
			# Remove case
			l = l.upper()

			text += l

	in_dict = fill_next_kmer_dict(text, k, in_dict)

	return kmer_count_dict

##############################
#Functions for generating text#
##############################

def add_sampled_kmer(in_dict, in_kmer, weights_df, current_author):
	# helper function for generate_text()
	# Args:
	# in_dict: dict, dictionary of next k-mer frequencies
	# in_kmer: the preceeding k-mer
	# weights_df: pd df, the data frame of z scores for all k-mers
	# current_author: the position of the current author,
	#### used to index the weights

	# Waht are the potential next k-mers
	potential_next_kmers = list(in_dict[in_kmer].keys())

	# For each of the potential NEXT k-mers, get the 
	# corresponding weight for this author
	weights = np.zeros(len(potential_next_kmers))
	for i in range(len(potential_next_kmers)):
		weights[i] =  weights_df[i][current_author]

	# adjust weights by frequency in text
	weights = weights**scaling_factor
	weights = weights * np.array(list(in_dict[in_kmer].values()))

	out_kmer = random.choices(potential_next_kmers,
			weights = weights)[0]

	return out_kmer

def generate_text(in_dict, target_length, weights_df, current_author):
	# Args:
	# in_dict: dict, dictionary of NEXT k-mers
	# target_length: int, length of desired sequence
	# weights_df: pd. df, df of z scores
	# current_author: int, number of current author  
	sampled_kmer = random.choices(list(in_dict.keys()))[0]
	out_string = sampled_kmer

	while (out_string[-1] != " ") or (len(out_string) < target_length):
		sampled_kmer = add_sampled_kmer(in_dict, sampled_kmer, weights_df, current_author)
		out_string += sampled_kmer
	
	return out_string

##############################
#Constructing k-mer freq table#
##############################

print("Identifying k-mers that are enriched for each author.")

# Generate a list of all authors and their books
authors_book_dict = {
	'willy' : ["/Users/andrew/texts/willy.txt"],
	'dickens' : ["/Users/andrew/texts/acc.txt",
	"/Users/andrew/texts/bleak.txt",
	"/Users/andrew/texts/dc.txt",
	"/Users/andrew/texts/great.txt",
	"/Users/andrew/texts/pick.txt",
	"/Users/andrew/texts/times.txt"],
	'beo' : ["/Users/andrew/texts/beo.txt"],
	'moby' : ["/Users/andrew/texts/moby.txt"],
	'joyce' : ["/Users/andrew/texts/ulysees.txt"],
	'grimm' : ["/Users/andrew/texts/grimm.txt"],
	'august' : ["/Users/andrew/texts/august.txt"],
	'bronte': ["/Users/andrew/texts/wuthering.txt"],
	'cant': ["/Users/andrew/texts/cant.txt"],
	'austen' : ["/Users/andrew/texts/northanger.txt",
	"/Users/andrew/texts/pride.txt",
	"/Users/andrew/texts/sense.txt",
	"/Users/andrew/texts/emma.txt"]
}

n_authors = len(authors_book_dict.keys())

# Initialize empty dictionary
kmer_count_dict = {}

# Initialize a list to store the total number of k-mers per author 
total_kmer_count_list = [0] * n_authors

author_names = list(authors_book_dict.keys())
# Point to a specific author
for i in range(len(author_names)):
	print(f"Reading in {author_names[i]}")
	# Pull out all books by this author
	for file_name in authors_book_dict[author_names[i]]:
		# Read in file and parse
		f = open(file_name)
		# Update k-mer counts and total number of k-mers for the author
		kmer_count_dict,  total_kmer_count_list = parse_text_construct_kf(kmer_count_dict, k, f, total_kmer_count_list, n_authors, i)
		f.close()

# Convert the values to an array, for easier operations:
kmer_count_dict_pd = pd.DataFrame.from_dict(kmer_count_dict, orient='index')

# Normalize k-mer count by the total size of each text 
kmer_count_dict_pd = kmer_count_dict_pd.div(total_kmer_count_list, axis=1) * 1000000

# Get the z scores 
kmer_count_z_scores = kmer_count_dict_pd.apply(zscore, axis=1)

# Shift the z scores so that they are non-negative
# Get min value
all_z_scores = np.concatenate(kmer_count_z_scores.values)
minZ = np.min(all_z_scores)
kmer_count_z_scores = kmer_count_z_scores - minZ


#####################################
#Reading in target text and generating#
#####################################



for author in author_names:
	print(f"Generating text for: {author}")

	current_author =  author_names.index(author)

	kmer_count_dict = {}
	if author == "willy":
		f = open("/Users/andrew/texts/willy.txt")
		kmer_count_dict = parse_text(kmer_count_dict, f)
		f.close()
	elif author == "dickens":
		f = open("/Users/andrew/texts/acc.txt")
		kmer_count_dict = parse_text(kmer_count_dict, f)
		f.close()
		f = open("/Users/andrew/texts/bleak.txt")
		kmer_count_dict = parse_text(kmer_count_dict, f)
		f.close()
		f = open("/Users/andrew/texts/dc.txt")
		kmer_count_dict = parse_text(kmer_count_dict, f)
		f.close()
		f = open("/Users/andrew/texts/great.txt")
		kmer_count_dict = parse_text(kmer_count_dict, f)
		f.close()
		f = open("/Users/andrew/texts/pick.txt")
		kmer_count_dict = parse_text(kmer_count_dict, f)
		f.close()
		f = open("/Users/andrew/texts/times.txt")
		kmer_count_dict = parse_text(kmer_count_dict, f)
		f.close()
	elif author == "beo":
		f = open("/Users/andrew/texts/beo.txt")
		kmer_count_dict = parse_text(kmer_count_dict, f)
		f.close()
	elif author == "moby":
		f = open("/Users/andrew/texts/moby.txt")
		kmer_count_dict = parse_text(kmer_count_dict, f)
		f.close()
	elif author == "joyce":
		f = open("/Users/andrew/texts/ulysees.txt")
		kmer_count_dict = parse_text(kmer_count_dict, f)
		f.close()
	elif author == "grimm":
		f = open("/Users/andrew/texts/grimm.txt")
		kmer_count_dict = parse_text(kmer_count_dict, f)
		f.close()
	elif author == "august":
		f = open("/Users/andrew/texts/august.txt")
		kmer_count_dict = parse_text(kmer_count_dict, f)
		f.close()
	elif author == "bronte":
		f = open("/Users/andrew/texts/wuthering.txt")
		kmer_count_dict = parse_text(kmer_count_dict, f)
		f.close()
	elif author == "cant":
		f = open("/Users/andrew/texts/cant.txt")
		kmer_count_dict = parse_text(kmer_count_dict, f)
		f.close()
	elif author == "austen":
		f = open("/Users/andrew/texts/northanger.txt")
		kmer_count_dict = parse_text(kmer_count_dict, f)
		f.close()
		f = open("/Users/andrew/texts/pride.txt")
		kmer_count_dict = parse_text(kmer_count_dict, f)
		f.close()
		f = open("/Users/andrew/texts/sense.txt")
		kmer_count_dict = parse_text(kmer_count_dict, f)
		f.close()
		f = open("/Users/andrew/texts/emma.txt")
		kmer_count_dict = parse_text(kmer_count_dict, f)
		f.close()
	print(generate_text(kmer_count_dict, target_length, kmer_count_z_scores, current_author))