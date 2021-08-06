def degen(input):
	dictionary = {'A': 'GCN', 'R': 'NNN', 'N': 'AAY', 'D': 'GAY', 'B': 'RAY', 
	'C': 'UGY', 'Q': 'CAR', 'E': 'GAR', 'Z': 'SAR', 'G': 'GGN', 'H': 'CAY', 'I': 'AUH',
	'L': 'NNN', 'K': 'AAR', 'F': 'UUY', 'P': 'CCN', 'S': 'NNN', 'T': 'ACN',
	'Y': 'UAY', 'V': 'GUN', 'M': 'AUG', 'W': 'UGG'}
	s = ''
	for char in input:
		s += dictionary[char]
	return s







