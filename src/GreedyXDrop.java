import java.util.ArrayList;

public class GreedyXDrop {
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

		GreedyXDropAlg(MSTB8, MSTC7);
	}

	public static void GreedyXDropAlg(String x, String y){

		// Scoring
		double MISMATCH = -2;
		double MATCH = 3;

		// Seq Lengths
		int M = x.length();
		int N = y.length();

		double[][] R = new double[M+N][N];
		double[] T = new double[N];

		double D = 0;
		double L = 0;
		double U = 0;
		double X = 0;

		double D_prime = 0;
		double[][] S_prime = new double[M+N][N];

		double i = 0;
		double j = 0;

		// start runtime timer 
		long startTime = System.nanoTime();

		while ((i < Math.min(M, N)) && x.charAt((int) i) == y.charAt((int) i)){
			i++;
		}

		double T_prime = T[0] = S_prime[(int) (i+i)][0];

		D++;
		D_prime = D - (X + (MATCH/2)/(MATCH-MISMATCH)) - 1;

		double c1 = Double.NEGATIVE_INFINITY;
		double c2 = Double.NEGATIVE_INFINITY;
		double c3 = Double.NEGATIVE_INFINITY;

		while (L <= U+2){
			for (double k = L-1; k < U+1; k++){
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
				if (i> Double.NEGATIVE_INFINITY && S_prime[(int) (i+j)][(int) D]>=T[(int) D_prime]-X){
					while (i<M && j<N && x.charAt((int) i) == y.charAt((int) j)) {
						i++;
						j++;
					}
					R[(int) D][(int) k] = i;
					T_prime = Math.max(T_prime,S_prime[(int) (i+j)][(int) D]);
				} else {
					R[(int) D][(int) k] = Double.NEGATIVE_INFINITY;
				}
			}

			T[(int) (D)] = T_prime;

			for (double k = 0; k < N; k++) {
				if(R[(int) D][(int) k] > Double.NEGATIVE_INFINITY){
					L = Math.min(L, R[(int) D][(int) k]);
				}
			}

			for (double k = 0; k < N; k++) {
				if(R[(int) D][(int) k] > Double.NEGATIVE_INFINITY){
					U = Math.max(L, R[(int) D][(int) k]);
				}
			}

			for (double k = 0; k < N; k++) {
				if(R[(int) D][(int) k] == N+k){
					U = Math.max(L, R[(int) D][(int) k] + 2);
				}
			}

			for (double k = 0; k < N; k++) {
				if(R[(int) D][(int) k] == M){
					L = Math.min(L, R[(int) D][(int) k] - 2);
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

}
