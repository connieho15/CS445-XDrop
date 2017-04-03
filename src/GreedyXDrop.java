import java.util.ArrayList;

public class GreedyXDrop {
	
	// Scoring
	protected static final double MISMATCH = -2;
	protected static final double MATCH = 3;
	protected static final double X = 3;
	
	public static void main(String[] args) {
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

		GreedyXDropAlg(MSTC7, MSTB8);
	}

	public static void GreedyXDropAlg(String x, String y){
		
		// Seq Lengths
		int M = x.length();
		int N = y.length();

		double[][] R = new double[N+M+1][M+1];
		double[] T = new double[N+1];			

		double D_prime = 0;			// Optimal distance
		

		double i = 0;
		double j = 0;

		// start runtime timer 
		long startTime = System.nanoTime();

		// Skips similar chars
		while ((i < Math.min(M, N)) && x.charAt((int) i) == y.charAt((int) i)){
			i++;
		}

		
		R[0][0] = i;
		double S_prime = computeS(i,i,0);
		double T_prime = T[0] = S_prime;

		double D = 0;			//Distance
		double L = 0;			//Lower bound
		double U = 0;			//Upper bound
		
		double c1 = Double.NEGATIVE_INFINITY;
		double c2 = Double.NEGATIVE_INFINITY;
		double c3 = Double.NEGATIVE_INFINITY;
		
		while (L <= U+2){
			if (D >M) {
			break;
		}
			D++;
			D_prime = D - Math.floor((X + (MATCH/2))/(MATCH-MISMATCH)) - 1;
			
			for (double k = L; k < U+1; k++){
				if (L<k){
					c1 = R[(int) (D-1)][(int) (k-1)] + 1;
				}
				if (L<=k && k<=U){
					c2 = R[(int) (D-1)][(int) k] + 1;	
				}
				if (k<U){
					c3 = R[(int) (D-1)][(int) (k+1)];
				}

				i = Math.max(c1, Math.max(c2, c3));

				j = i - k;
				
				if (i> Double.NEGATIVE_INFINITY && computeS(i,j,D)>=T[(int) D_prime]-X){
					while (i<M && j<N && x.charAt((int) i) == y.charAt((int) j)) {
						i++;
						j++;
					}
					R[(int) D][(int) k] = i;
					T_prime = Math.max(T_prime,computeS(i,j,D));
				} else {
					R[(int) D][(int) k] = Double.NEGATIVE_INFINITY;
				}
			}

			T[(int) (D)] = T_prime;

			for (double k = 0; k < M; k++) {
				if( R[(int) D][(int) k] > Double.NEGATIVE_INFINITY){
					L = Math.min(L, R[(int) D][(int) k]);
				}
			}

			for (double k = 0; k < M; k++) {
				if(R[(int) D][(int) k] > Double.NEGATIVE_INFINITY){
					U = Math.max(U, R[(int) D][(int) k]);
				}
			}

			for (double k = 0; k < M; k++) {
				System.out.println(R[(int) D][(int) k]);
				if(R[(int) D][(int) k] == N+k){
					L = Math.max(L, R[(int) D][(int) k] + 2);
				}
			}

			for (double k = 0; k < M; k++) {
				if(R[(int) D][(int) k] == M){
					U = Math.min(U, R[(int) D][(int) k] - 2);
				}
			}


		}
		// return the best overall score
		System.out.println("Alignment Score: " + T_prime);

		// stop runtime timer
		long endTime = System.nanoTime();
		long duration = (endTime - startTime) / 1000000;  //divide by 1000000 to get milliseconds.
		System.out.println("Runtime: " + duration + "ms");

	}
	
	private static double computeS(double i, double j, double d){
		return (i+j)*(MATCH/2)-d*(MATCH-MISMATCH);
	}

}
