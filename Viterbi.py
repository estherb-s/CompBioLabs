import numpy as np

neginf = float("-inf")

# Hidden markov model for detecting CG rich regions
# The alphabet of the nucleotide sites
nucleotide_alphabet = "ACGT" 
# The alphabet of the states (start, CG High, CG Low)
state_alphabet = "SHL" 

# Emissions probability matrix, the rows = the three states (start, CG-rich, CG-poor/AT-rich)
# Columns = probability of emitting [A, C, G, T] depending on the state
emiss_probs = np.array([
	[0.00, 0.00, 0.00, 0.00],
	[0.2459, 0.2079, 0.2478, 0.2984],
	[0.2698, 0.3237, 0.2080, 0.1985]
])

# Transition probabilities matrix, the rows = state at site i - 1,
# The columns = probability of each state at site i
trans_probs = np.array([
	[0.00, 0.50, 0.50],
	[0.00, 0.9997, 0.0003],
	[0.00, 0.0002, 0.9998]
])

emiss_log_probs = np.log(emiss_probs)
trans_log_probs = np.log(trans_probs)

sequence = open('lambda.txt').read().replace('\n', '')
# sequence = "GCGGCCGAGGT"

# Length of the Viterbi matrix
length = len(sequence) + 1 

# Log-probabilities
v_matrix = np.zeros((length, 3)) 
# Pointers
p_matrix = np.zeros((length, 3), dtype = "int8") 

v_matrix.fill(neginf)
p_matrix.fill(-1)
v_matrix[0, 0] = 0.0

for i in range(1, length):
	nucleotide = nucleotide_alphabet.index(sequence[i - 1])
	for j in range(3):
		for k in range(3):
			vij = emiss_log_probs[j, nucleotide] + v_matrix[i - 1, k] + trans_log_probs[k, j]			
			if vij > v_matrix[i, j]:
				v_matrix[i, j] = vij
				p_matrix[i, j] = k

next_map_k = np.argmax(v_matrix[length - 1])
map_states = ""
for i in reversed(range(1, length)):
	map_states = state_alphabet[next_map_k] + map_states
	next_map_k = p_matrix[i, next_map_k]

# print the map states:
# First seq = true simulated path
# Second seq = estimated path using Viterbi algorithm (using above code)
print(map_states)




