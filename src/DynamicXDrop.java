
public class DynamicXDrop {
	static // MRSA strain MSTB8 mecA gene sequence
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
			static String MSTC7 = "tggggtggttacaacgttacaagatatgaagtggtaaatggtaatatcgacttaa"
					+ "aacaagcaatagaatcatcagataacattttctttgctagagtagcactcgaatta"
					+ "ggcagtaagaaatttgaaaaaggcatgaaaaaactaggtgttggtgaagatatacc"
					+ "aagtgattatccattttataatgctcaaatttcaaacaaaaatttagataatgaaa"
					+ "tattattagctgattcaggttacggacaaggtgaaatactgattaacccagtacag"
					+ "atcctttcaatctatagcgcattagaaaataatggcaatattaacgcacctcactt"
					+ "attaaaagacacgaaaaacaaagtttggaagaaaaatattatttccaaagaaaata"
					+ "tcaatctattaactgatggtatgcaacaagtcgtaaataaaacacataaagaag";
			
	public static void main(String[] args) {
		dynamicXDrop(MSTB8, MSTC7);
	}

	public static void dynamicXDrop(String a, String b) {
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
		
		double best_i = 0;
		double best_j = 0;
		
		int N = a.length();
		int M = b.length();
		
		// scoring matrix s 
		CustomTwoDArray s = new CustomTwoDArray(N+1, M+1);
		s.set(0, 0, 0);
		
		// pointer matrix
		CustomTwoDArray ptr = new CustomTwoDArray(N+1, M+1);
		
		// start runtime timer 
		long startTime = System.nanoTime();
		
		// last state visited (match, mismatch, insertion or deletion)
		double ptrp = 0.0;
		double best_p = 0;
		
		while (L <= U+1) {
			k++;
			for (double i = Math.ceil(L); i <= Math.floor(U)+1; i = i + half) {
				double j = k-i; 
				
				if ((i+half < N) && (j+half < M)) {
					System.out.println("i: " + i + "j: " + j);
					// do this if values are not halves 
					if (i == (int) i) {
						
						double s1 = Double.NEGATIVE_INFINITY;
						double s2 = Double.NEGATIVE_INFINITY;
						double s3 = Double.NEGATIVE_INFINITY;
						double s4 = Double.NEGATIVE_INFINITY;					
						
						// match case:
						if (L <= i-half && U >= i-half && a.charAt((int) i) == (b.charAt((int) j))) {
							s1 = s.get(i-half,  j-half) + (mat/2);
							System.out.println("	MAT " + s1);
						// mismatch case:
						} if (L <= i-half && U >= i-half && a.charAt((int) i) != (b.charAt((int) j))) {
							s2 = s.get(i-half,  j-half) + (mis/2);
							System.out.println("	MIS " + s2);
						// deletion case:
						} if (i <= U) {
							s3 = s.get(i, j-1) + ind;
							System.out.println("	DEL " + s3);
						// insertion case:
						} if (L <= i-1) {
							s4 = s.get(i-1, j) + ind;
							System.out.println("	INS " + s4);
						}
						if (s1 == Math.max(s1, Math.max(s2, Math.max(s3, s4)))) {ptrp = 1;} else {
							if (s2 == Math.max(s1, Math.max(s2, Math.max(s3, s4)))) {ptrp = 2;} else {
								if (s3 == Math.max(s1, Math.max(s2, Math.max(s3, s4)))) {ptrp = 3;} else {
									if (s4 == Math.max(s1, Math.max(s2, Math.max(s3, s4)))) {
										ptrp = 4;}}}}
						ptr.set(i, j, ptrp);
						s.set(i, j, Math.max(s1, Math.max(s2, Math.max(s3, s4))));
						
					
					// do this if values are half values
					} else {
						double score = 0;
						// match case: 
						if (a.charAt((int) (i+half)) == b.charAt((int) (j+half))) {
							score = mat/2;
							System.out.println("	MAT " + score);
							ptrp = 1;
						// mismatch case: 
						} else if (a.charAt((int) (i+half)) != b.charAt((int) (j+half))) {
							score = mis/2;
							System.out.println("	MIS " + score);
							ptrp = 2;
						}
						s.set(i, j, s.get(i-half, j-half) + score);
						ptr.set(i, j, ptrp);
					}
					
					// set best overall score
					if (s.get(i,  j) > T_prime){
						best_i = i;
						best_j = j;
						best_p = ptrp;
						System.out.println(ptr.get(best_i, best_j));
						System.out.println(s.get(i, j));
					}
					
					T_prime = Math.max(T_prime, s.get(i, j));
					
					// X-drop: Don't consider extensions of alignment 
					// more than X below best score seen so far
					if (s.get(i, j) < T-X) {
						s.set(i, j, Double.NEGATIVE_INFINITY);
					}
				}
			}

			// trim upper and lower bounds:
			for (double p = Math.ceil(L); p < N; p++) {
				if ((k-p > 0 && k-p < N && k-p < M) && s.get(p, k - p) > Double.NEGATIVE_INFINITY) {
					L = p;
					break;
				}
			}
			

			for (double q = Math.ceil(L); q < N; q++) {
				if ((k-q > 0 && k-q < N && k-q < M) && s.get(q, k - q) > Double.NEGATIVE_INFINITY) {
						U = (int) q;
				}
			} 
			
			L = Math.max(L, k+1-N);
			U = Math.min(U, M-1);
			T = T_prime;	// update best score seen so far
		}

		// return the best overall score
		System.out.println("Alignment Score: " + T_prime);
		
		System.out.println("besti: " + best_i + "bestj: " + best_j);
		
		// Traceback
		double pi = best_i;
		double pj = best_j;
		String x = "";
		String y = "";
		String ptrs = "";
		while (pi > 0 && pj > 0 && ptr.get(pi, pj) > 0){
			int p = (int) ptr.get(pi, pj);
			ptrs = p + ptrs;
			// end on a match or mismatch
			if (p==1 || p==2) {
				x = a.charAt((int)pi) + x;
				y = b.charAt((int)pj) + y;
				pi=pi-1;
				pj=pj-1;
			} else 
				// insertion
				if (p==4){
					System.out.println("INSERTION AT" + pi + pj);
					x = a.charAt((int)pi) + x;
					y = "-" + y;
					pi--;
					// deletion
				} else if (p==3){
					x = "-" + x;
					y = b.charAt((int)pj) + y;
					pj--;
				}

		}
	
			System.out.println(x);
			System.out.println(y);
			System.out.println(ptrs);
		
		
		// stop runtime timer
		long endTime = System.nanoTime();
		long duration = (endTime - startTime) / 1000000;  
		System.out.println("Runtime: " + duration + "ms");
	}
	
	
}


