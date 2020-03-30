import numpy as np

# Hidden markov model for detecting CG rich regions
# The alphabet of the nucleotide sites
nucleotide_alphabet = "ACGT" 
# The alphabet of the states (start, CG High, CG Low)
state_alphabet = "SHL" 

# Emissions probability matrix, the rows = the three states (start, CG-rich, CG-poor/AT-rich)
# Columns = probability of emitting [A, C, G, T] depending on the state
emission_probs = np.array([
	[0.00, 0.00, 0.00, 0.00],
	[0.2459, 0.2079, 0.2478, 0.2984],
	[0.2698, 0.3237, 0.2080, 0.1985]
])

# Transition probabilities matrix, the rows = state at site i - 1,
# The columns = probability of each state at site i
transition_probs = np.array([
	[0.00, 0.50, 0.50],
	[0.00, 0.9997, 0.0003],
	[0.00, 0.0002, 0.9998]
])

# Start at the start state
current_state = 0 
# Start with an empty sequence
sim_sites = "" 
# Record the sequence of states
sim_states = "" 

# Length of simulated sequence
length = 100

for i in range(length):
	# Choose a new state at site i, weighted by the transitional probability matrix
	current_state = np.random.choice(3, p = transition_probs[current_state])
	# Choose a new nucleotide at site i, weighted by the emission probability matrix
	nucleotide = np.random.choice(4, p = emission_probs[current_state])

	# Append the letter codes of the chosen state and nucleotide
	sim_states += state_alphabet[current_state]
	sim_sites += nucleotide_alphabet[nucleotide]

print(sim_states)
print(sim_sites)
