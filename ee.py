import re
import sys
import random
import pandas as pd

k = 4
stanzas = 4

if len(sys.argv) > 1:
	k = int(sys.argv[1])

text = pd.read_csv("poems.csv")

def count_words_per_line(in_poem,in_list):
	# Returns a list of how many words there are per line
	# Args:
	# in_poem: str, the full text of the poem to analyze
	in_poem = in_poem.split("\n")
	# Dealing with multiple spaces:
	# Compress multiple spaces to one
	in_poem = [re.sub(' +', ' ', x) for x in in_poem]
	# Remove left and right spaces
	in_poem = [x.lstrip() for x in in_poem]
	in_poem = [x.rstrip() for x in in_poem]
	# Count spaces
	counts = [x.count(' ') + 1 for  x in in_poem if x != ""]

	return in_list + counts

def count_lines_per_stanza(in_poem, in_list):
	# Returns a list of how many lines there are in a stanza
	# Split poem into lines
	# Args:
	# in_poem: str, the poem to parse
	in_poem = in_poem.split("\n")
	# Get the positions of the stanza breaks
	line_breaks = [i for i, val in enumerate(in_poem) if val == ""]
	line_breaks = [-1] + line_breaks # negative 1 is the index of the "first" line break before line 0
	# Get wait time (in units of lines)
	
	return in_list + [line_breaks[x] - (line_breaks[x-1] + 1) for x in range(1, len(line_breaks))]
	
def get_word_start_kmers(in_poem, k, in_list):
	# Get the k-mers that start all of our words.
	# Convert newlines and tabs to spaces
	in_poem = re.sub('\n', '  ', in_poem)
	in_poem = re.sub('\t', ' ', in_poem)
	# Remove non alphanumerical characters
	in_poem = re.sub(r'[^a-zA-Z\ ]', '', in_poem)
	# Remove consecutive spaces
	in_poem = re.sub(' +', ' ', in_poem)

	# Separate into words
	words = in_poem.split(" ")

	# Get starting k-mers
	in_list += [x[0:k] for x in words if len(x) > k]

	return in_list

def parse_poem(in_poem, k, in_dict):
	# Function takes a poem and parses it,
	# removing new lines, punctuation, spaces
	# Args:
	# in_poem: string, the entire text of the poem to analyze
	# in_dict: dict, the dictionary with next k-mers

	# Convert newlines and tabs to spaces
	in_poem = re.sub('\n', '  ', in_poem)
	in_poem = re.sub('\t', ' ', in_poem)
	# Remove non alphanumerical characters
	in_poem = re.sub(r'[^a-zA-Z\ ]', '', in_poem)
	# Remove consecutive spaces
	in_poem = re.sub(' +', ' ', in_poem)
	# Remove case
	in_poem = in_poem.lower()

	in_dict = fill_next_kmer_dict(in_poem, k, in_dict)

	return in_dict

def fill_next_kmer_dict(in_string, k, in_dict):
	for i in range(len(in_string) - k - k + 1):
		current_kmer = in_string[i:i+k]
		next_kmer = in_string[i + k:i + k + k]
		in_dict.setdefault(current_kmer, dict())
		in_dict[current_kmer].setdefault(next_kmer, 0)
		in_dict[current_kmer][next_kmer] += 1

	return(in_dict)

def add_sampled_kmer(in_dict, in_kmer):
	# helper function for generate_text()
	try:
  		out_kmer = random.choices(list(in_dict[in_kmer].keys()),
			weights = list(in_dict[in_kmer].values()))[0]
	except:
  		out_kmer = random.choices(start_kmer_list)[0]
	
	return out_kmer

def generate_text(start_kmer_list, in_dict):
	sampled_kmer = random.choices(start_kmer_list)[0]
	out_string = sampled_kmer
	n_stanzas = 0
	
	n_target_stanzas = 4

	# write poem
	while (out_string[-1] != "&") or (n_stanzas < n_target_stanzas):
		n_target_words = random.choices(line_lengths)[0]
		n_lines = 0
		n_target_lines = random.choices(stanza_lengths)[0]
		# Write stanza
		while (out_string[-1] != "\n") or n_lines < n_target_lines:
			n_words = 0
			n_target_words = random.choices(line_lengths)[0]
			while (out_string[-1] != " ") or n_words < n_target_words:
				sampled_kmer = add_sampled_kmer(in_dict, sampled_kmer)
				out_string += sampled_kmer
				n_words = out_string.split("\n")[-1].count(' ')
				#print(out_string)
			out_string += "\n"
			n_lines = out_string.split("&")[-1].count("\n")
		out_string += "&"
		n_stanzas = out_string.count("&")

	
	return out_string

next_kmer_dict = {}
starting_kmers = []
stanza_lengths = []
line_lengths = []
for i in range(text.shape[0]):
	title = text.loc[:,"first_line"].iloc[i]
	print(f"Reading {title}")

	poem = text.loc[:,"text"].iloc[i]
	
	next_kmer_dict = parse_poem(poem, k, next_kmer_dict)
	starting_kmers = get_word_start_kmers(poem, k, starting_kmers)
	stanza_lengths = count_lines_per_stanza(poem, stanza_lengths)
	line_lengths = count_words_per_line(poem, line_lengths)


print("\n\n")
unformatted_poem = generate_text(starting_kmers, next_kmer_dict)
out_poem = re.sub('&', '\n\n', unformatted_poem)
print(out_poem)
#print(repr(poem))
#print(count_lines_per_stanza(poem))
#print(count_words_per_line(poem))
#print(parse_poem(poem))
#print(get_word_start_kmers(poem, k))