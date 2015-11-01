import java.util.HashSet;


public class Predictor {
	final int interested_label = 1; // looking at first ">"
	final String input = "GCF_000091665.1_ASM9166v1_genomic.fna";
	final String output = "GCF_000091665.1_ASM9166v1_genomic.gbff";
	
	HashSet<Integer> A_indices = new HashSet<Integer>();
	HashSet<Integer> C_indices = new HashSet<Integer>();
	HashSet<Integer> G_indices = new HashSet<Integer>();
	HashSet<Integer> T_indices = new HashSet<Integer>(); // will be empty
	
	public static void main(String[] args) {
		
	}

}
