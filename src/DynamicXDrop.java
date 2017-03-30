
public class DynamicXDrop {
	public static void main(String[] args) {
		dynamicXDrop();
	}

	public static void dynamicXDrop() {
		double T_prime = 0;
		double T = 0;
		double k = 0;
		double L = 0;
		double U = 0; 
		double X = 10;
		
		double mat = +3;
		double mis = -2;
		double ind = mis - (mat/2);
		double half = 0.5;
		
		String a = "ABC";
		String b = "ABC";
		int N = a.length();
		int M = b.length();
		
		// Scoring matrix s 
		double[][] s = new double[N][M];
		s[0][0] = 0;
		
		while (L <= U+1) {
			k++;
			for (double i = Math.ceil(L); i <= Math.floor(U)+1; i = i + half) {
				double j = k-i; 
				
				if ((i+half < N) && (j+half < M)) {
					if (i == (int) i) {
						double s1 = Double.NEGATIVE_INFINITY;
						double s2 = Double.NEGATIVE_INFINITY;
						double s3 = Double.NEGATIVE_INFINITY;
						double s4 = Double.NEGATIVE_INFINITY;					
						
						if (L <= i-half && U > i-half && a.charAt((int) i) == (b.charAt((int) j))) {
							s1 = s[(int) (i-half)][(int) (j-half)] + (mat/2);
						} if (L <= i-half && U > i-half && a.charAt((int) i) != (b.charAt((int) j))) {
							s2 = s[(int) (i-half)][(int) (j-half)] + (mis/2);
						} if (i <= U) {
							s3 = s[(int) i][(int) (j-half)] + ind;
						} if (L <= i-1) {
							s4 = s[(int) (i-half)][(int) j] + ind;
						}
		
						s[(int) i][(int) j] = Math.max(s1, Math.max(s2, Math.max(s3, s4)));
						
					} else {
						double score = 0;
						if (a.charAt((int) (i+half)) == b.charAt((int) (j+half))) {
							score = mat/2;
						} else if (a.charAt((int) (i+half)) != b.charAt((int) (j+half))) {
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
			
			L = Math.max(L, k+1-N);
			U = Math.min(U, M-1);
			T = T_prime; // Update best score seen so far
		}
		System.out.println(T_prime);
	}
}
