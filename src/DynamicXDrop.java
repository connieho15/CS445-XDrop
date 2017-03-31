
public class DynamicXDrop {
	public static void main(String[] args) {
		dynamicXDrop();
	}

	public static void dynamicXDrop() {
		// best overall score
		double T_prime = 0;
		
		// best alignment score
		double T = 0;
		
		// counter
		double k = 0;
		
		// lower bound for x-coordinate
		double L = 0;
		
		// upper bound for x-coordinate 
		double U = 0; 
		double X = 3;
		
		// match score
		double mat = 3;
		
		// mismatch score
		double mis = -2;
		
		// insertion/deletion score
		double ind = mis - (mat/2);
		double half = 0.5;
		
		// MRSA strain MSTB8 mecA gene sequence
		String MSTB8 = "tggggtggttacaacgttacaagatatgaagtggtaaatggtaatatcgacttaa"
				+ "aacaagcaatagaatcatcagataacattttctttgctagagtagcactcgaatta"
				+ "ggcagtaagaaatttgaaaaaggcatgaaaaaactaggtgttggtgaagatatacc"
				+ "aagtgattatccattttataatgctcaaatttcaaacaaaaatttagataatgaaa"
				+ "tattattagctgattcaggttacggacaaggtgaaatactgattaacccagtacag"
				+ "atcctttcaatctatagcgcattagaaaataatggcaatattaacgcacctcactt"
				+ "attaaaagacacgaaaaacaaagtttggaagaaaaatattatttccaaagaaaata"
				+ "tcaatctattaactgatggtatgcaacaagtcgtaaataaaacacataaagaagat"
				+ "atttatagatcttatgca";
		// MRSA strain MSTC7 mecA gene sequence
		String MSTC7 = "tggggtggttacaacgttacaagatatgaagtggtaaatggtaatatcgacttaa"
				+ "aacaagcaatagaatcatcagataacattttctttgctagagtagcactcgaatta"
				+ "ggcagtaagaaatttgaaaaaggcatgaaaaaactaggtgttggtgaagatatacc"
				+ "aagtgattatccattttataatgctcaaatttcaaacaaaaatttagataatgaaa"
				+ "tattattagctgattcaggttacggacaaggtgaaatactgattaacccagtacag"
				+ "atcctttcaatctatagcgcattagaaaataatggcaatattaacgcacctcactt"
				+ "attaaaagacacgaaaaacaaagtttggaagaaaaatattatttccaaagaaaata"
				+ "tcaatctattaactgatggtatgcaacaagtcgtaaataaaacacataaagaag";
		int N = MSTB8.length();
		int M = MSTC7.length();
		
		// Scoring matrix s 
		double[][] s = new double[N+1][M+1];
		s[0][0] = 0;

		long startTime = System.nanoTime();
		
		while (L <= U+1) {
			k++;
			for (double i = Math.ceil(L); i <= Math.floor(U)+1; i = i + half) {
				double j = k-i; 
				
				if ((i+half < N) && (j+half < M)) {
					// Do this if values are not halves 
					if (i == (int) i) {
						double s1 = Double.NEGATIVE_INFINITY;
						double s2 = Double.NEGATIVE_INFINITY;
						double s3 = Double.NEGATIVE_INFINITY;
						double s4 = Double.NEGATIVE_INFINITY;					
						
						if (L <= i-half && U >= i-half && MSTB8.charAt((int) i) == (MSTC7.charAt((int) j))) {
							s1 = s[(int) (i-half)][(int) (j-half)] + (mat/2);
						} if (L <= i-half && U >= i-half && MSTB8.charAt((int) i) != (MSTC7.charAt((int) j))) {
							s2 = s[(int) (i-half)][(int) (j-half)] + (mis/2);
						} if (i <= U) {
							s3 = s[(int) i][(int) (j-1)] + ind;
						} if (L <= i-1) {
							s4 = s[(int) (i-1)][(int) j] + ind;
						}
		
						s[(int) i][(int) j] = Math.max(s1, Math.max(s2, Math.max(s3, s4)));
					
					// Do this if values are half values
					} else {
						double score = 0;
						if (MSTB8.charAt((int) (i+half)) == MSTC7.charAt((int) (j+half))) {
							score = mat/2;
						} else if (MSTB8.charAt((int) (i+half)) != MSTC7.charAt((int) (j+half))) {
							score = mis/2;
						}
						s[(int) i][(int) j] = s[(int) (i-half)][(int) (j-half)] + score;
					}
					
					T_prime = Math.max(T_prime, s[(int) i][(int) j]);
					
					// X-drop: Don't consider extensions of alignment 
					// more than X below best score seen so far
					if (s[(int) i][(int) j] < T-X) {
						s[(int) i][(int) j] = Double.NEGATIVE_INFINITY;
					}
				}
			}

			for (double i = Math.ceil(L); i < N; i++) {
				if ((k - i > 0) && s[(int) i][(int) k - (int) i] > Double.NEGATIVE_INFINITY) {
					L = i;
					break;
				}
			}

			int maxSoFar = 0;
			for (double i = Math.ceil(L); i < N; i++) {
				if ((k - i > 0) && s[(int) i][(int) k - (int) i] > Double.NEGATIVE_INFINITY) {
					if (i > maxSoFar) {
						maxSoFar = (int) i;
					}
				}
			} U = maxSoFar;

			L = Math.max(L, k+1-N);
			U = Math.min(U, M-1);
			T = T_prime; // Update best score seen so far
		}

		System.out.println("Alignment Score: " + T_prime);
		
		long endTime = System.nanoTime();
		long duration = (endTime - startTime) / 1000000;  //divide by 1000000 to get milliseconds.
		System.out.println("Runtime: " + duration + "ms");
	}
}
