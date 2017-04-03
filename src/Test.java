public class Test {

	public static void main(String[] args) {
		
		// Two completely different sequences should return a score of 0
		// Because this algorithm is meant for closely related sequences, no alignment will be produced	
		System.out.println("Aligning 2 identical sequences:");
		DynamicXDrop.dynamicXDrop("GGGGGGAAAGGGGGGG", "GGGGGGGAAAGGGGGG");
		
		System.out.println("\nAligning 2 identical sequences with 1 deletion:");
		DynamicXDrop.dynamicXDrop("GGGGGGGGGGGG", "GGGGGGGGAGGG");
		
		System.out.println("\nAligning 2 identical sequences with 1 insertion:");
		DynamicXDrop.dynamicXDrop("GGGGGGAGGG", "GGGGGGGGGG");
		
		System.out.println("\nAligning 2 different length sequences:");
		DynamicXDrop.dynamicXDrop("GGGGGCCCGGGGGGGG", "GGGGGCCCGGGGGGGGGGAAAAAA");	
		
		System.out.println("\nAligning 2 completely different sequences:");
		DynamicXDrop.dynamicXDrop("AAAAAATTTAAAAAAAAA", "GGGGGCCCGGGGGGGGGG");	
		// TODO: Tests for Greedy method 
		System.out.println("GREEDY ALG:");
		System.out.println("Aligning 2 identical sequences:");
		GreedyXDrop.GreedyXDropAlg("GGGGGGAAAGGGGGGG", "GGGGGGGAAAGGGGGG");
		
		System.out.println("\nAligning 2 identical sequences with 1 deletion:");
		GreedyXDrop.GreedyXDropAlg("GGGGGGGGGGGG", "GGGGGGGGAGGG");
		
		System.out.println("\nAligning 2 identical sequences with 1 insertion:");
		GreedyXDrop.GreedyXDropAlg("GGGGGGAGGG", "GGGGGGGGGG");
		
		System.out.println("\nAligning 2 different length sequences:");
		GreedyXDrop.GreedyXDropAlg("GGGGGCCCGGGGGGGG", "GGGGGCCCGGGGGGGGGGAAAAAA");	
		
		System.out.println("\nAligning 2 completely different sequences:");
		GreedyXDrop.GreedyXDropAlg("AAAAAATTTAAAAAAAAA", "GGGGGCCCGGGGGGGGGG");	
	 
	} 

}
