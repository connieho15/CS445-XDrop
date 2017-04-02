
public class CustomTwoDArray {
	private double[][] array;
	
	// custom two d array that works like a regular array except
	// double in length and width to accommodate 
	// S[i-1/2][j-1/2] "half nodes"
	public CustomTwoDArray(int N, int M){
		this.array = new double[2*N][2*M];
	}
	
	// gets value from array
	// s[i][j]
	public double get(double i, double j){
		return array[(int)(2*i)][(int)(2*j)];
	}
	
	// sets a new value in array
	// S[i][j] = n
	public void set(double i, double j, double n){
		array[(int)(2*i)][(int)(2*j)] = n;
	}
}
