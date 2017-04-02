public class Test {

	public static void main(String[] args) {
		
		// Testing DynamicXDrop
		System.out.println("Testing...");
		
		// Two completely different sequences should return a score of 0
		// Because this algorithm is meant for closely related sequences, no alignment will be produced
		DynamicXDrop.dynamicXDrop("AAAAAAAA", "GGGGGGGG");
		
		// Two completely same sequences should return a high score and all-match alignment
		DynamicXDrop.dynamicXDrop("CAAAAAAAAACAAAA", "CAAAAAAAAACAAAA");
		
		// TODO: Tests for Greedy method 
		
	} 

}
