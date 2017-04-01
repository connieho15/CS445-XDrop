
public class DynamicXDrop {
	public static void main(String[] args) {
		dynamicXDrop();
	}

	public static void dynamicXDrop() {
		double T_prime = 0;			// best overall score
		double T = 0;				// best alignment score
		double k = 0;				// counter 
		double L = 0;				// lower bound for x-coordinate
		double U = 0; 				// upper bound for x-coordinate 
		double X = 10;				// X for X drop approach
		
		double mat = 3;				// match score
		double mis = -2; 			// mismatch score
		double ind = mis - (mat/2);	// insertion/deletion score
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
		
		// scoring matrix s 
		double[][] s = new double[N+1][M+1];
		s[0][0] = 0;

		// start runtime timer 
		long startTime = System.nanoTime();
		
		while (L <= U+1) {
			k++;
			for (double i = Math.ceil(L); i <= Math.floor(U)+1; i = i + half) {
				double j = k-i; 
				
				if ((i+half < N) && (j+half < M)) {
					// do this if values are not halves 
					if (i == (int) i) {
						double s1 = Double.NEGATIVE_INFINITY;
						double s2 = Double.NEGATIVE_INFINITY;
						double s3 = Double.NEGATIVE_INFINITY;
						double s4 = Double.NEGATIVE_INFINITY;					
						
						// match case:
						if (L <= i-half && U >= i-half && MSTB8.charAt((int) i) == (MSTC7.charAt((int) j))) {
							s1 = s[(int) (i-half)][(int) (j-half)] + (mat/2);
						// mismatch case:
						} if (L <= i-half && U >= i-half && MSTB8.charAt((int) i) != (MSTC7.charAt((int) j))) {
							s2 = s[(int) (i-half)][(int) (j-half)] + (mis/2);
						// deletion case:
						} if (i <= U) {
							s3 = s[(int) i][(int) (j-1)] + ind;
						// insertion case:
						} if (L <= i-1) {
							s4 = s[(int) (i-1)][(int) j] + ind;
						}
		
						s[(int) i][(int) j] = Math.max(s1, Math.max(s2, Math.max(s3, s4)));
					
					// do this if values are half values
					} else {
						double score = 0;
						// match case: 
						if (MSTB8.charAt((int) (i+half)) == MSTC7.charAt((int) (j+half))) {
							score = mat/2;
						// mismatch case: 
						} else if (MSTB8.charAt((int) (i+half)) != MSTC7.charAt((int) (j+half))) {
							score = mis/2;
						}
						s[(int) i][(int) j] = s[(int) (i-half)][(int) (j-half)] + score;
					}
					
					// set best overall score
					T_prime = Math.max(T_prime, s[(int) i][(int) j]);
					
					// X-drop: Don't consider extensions of alignment 
					// more than X below best score seen so far
					if (s[(int) i][(int) j] < T-X) {
						s[(int) i][(int) j] = Double.NEGATIVE_INFINITY;
					}
				}
			}

			// trim upper and lower bounds:
			for (double p = Math.ceil(L); p < N; p++) {
				if ((k-p > 0 && k-p < N && k-p < M) && s[(int) p][(int) k - (int) p] > Double.NEGATIVE_INFINITY) {
					L = p;
					break;
				}
			}
			
			double maxSoFar = Double.NEGATIVE_INFINITY;
			for (double q = Math.ceil(L); q < N; q++) {
				if ((k-q > 0 && k-q < N && k-q < M) && s[(int) q][(int) k - (int) q] > Double.NEGATIVE_INFINITY) {
					if (q > maxSoFar) {
						maxSoFar = (int) q;
					}
				}
			} 
			
			U = maxSoFar;
			L = Math.max(L, k+1-N);
			U = Math.min(U, M-1);
			T = T_prime;	// update best score seen so far
		}

		// return the best overall score
		System.out.println("Alignment Score: " + T_prime);
		
		// stop runtime timer
		long endTime = System.nanoTime();
		long duration = (endTime - startTime) / 1000000;  //divide by 1000000 to get milliseconds.
		System.out.println("Runtime: " + duration + "ms");
	}
}
