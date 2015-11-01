import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;


public class Predictor {
	final int interested_label = 1; // looking at first ">"
	final String input = "GCF_000091665.1_ASM9166v1_genomic.fna";
	final String output = "GCF_000091665.1_ASM9166v1_genomic.gbff";
	
	static HashSet<Integer> A = new HashSet<Integer>();
	static HashSet<Integer> C = new HashSet<Integer>();
	static HashSet<Integer> G = new HashSet<Integer>();
	static HashSet<Integer> T = new HashSet<Integer>();
	
	// mod, orfts, short_orfs, all share same reference
	static ArrayList<ORF> mod0 = new ArrayList<ORF>();
	static ArrayList<ORF> mod1 = new ArrayList<ORF>();
	static ArrayList<ORF> mod2 = new ArrayList<ORF>();
	
	static ArrayList<ORF> orfs = new ArrayList<ORF>();
	static ArrayList<ORF> short_orfs = new ArrayList<ORF>();
	static ArrayList<ORF> long_orfs = new ArrayList<ORF>();
	
	static ArrayList<ORF> cds = new ArrayList<ORF>();

	// String is three letters (k), and score contains count of each letter k+1 
	static HashMap<String, Score> P = new HashMap<String, Score>(); 
	static HashMap<String, Score> Q = new HashMap<String, Score>();
	
	public static void main(String[] args) {

		// Part 1: Reading files and parsing
		// parse and find all ORFS in three reading frames
		ORFParser.findORF(A, T, C, G, mod0, mod1, mod2, orfs);
		report(false, false, true);
		
		// process ORFs
		ORFParser.shorterORFs(orfs, short_orfs, 50);
		ORFParser.longerORFs(orfs, long_orfs, 1400);
		System.out.println("shorter: " + short_orfs.size() + ", longer: " + long_orfs.size());
		
		// find all CDS & matching
		CDSParser.parse(cds);
		int totalCDS = CDSParser.scoreORF(orfs, cds);
		System.out.println("score: " + totalCDS + "/" + cds.size());
		
		calculateProb(short_orfs, Q);
		calculateProb(long_orfs, P);
		
		verifyProb(long_orfs, P);
		verifyProb(short_orfs, Q);
		
		String[] temp = new String[]{"A", "C", "G", "T"};
		System.out.println("\tA\t\tC\t\tG\t\tT");
		for (int i = 0; i < 4; i++) {
			System.out.print(temp[i] + "\t");
			for (int j = 0; j < 4; j++) {
				String key = "A" + temp[i] + temp[j];
				System.out.print(P.get(key).A);
				System.out.print("|");
				System.out.print(Q.get(key).A);
				System.out.print("\t");
			}
			System.out.println();
		}
		
	}
	
	private static void verifyProb(ArrayList<ORF> long_orfs, HashMap<String, Score> p2) {
		// sum of sum = length of the orfs
		for (int i = 0; i < orfs.size(); i++) {
			int length = orfs.get(i).count - 1;
			int count = 0;
			for (String key : p2.keySet()) {
				count += p2.get(key).sum;
			}
			assert(length == count);
		}
	}

	public static void calculateProb(ArrayList<ORF> orfs, HashMap<String, Score> result){
		for (int i = 0; i < orfs.size(); i++) {
			int start = orfs.get(i).start;
			int finish = orfs.get(i).finish;
			
			for (int j = start; j + 3 < finish - 2; j++) {
				String triplet = getNT(j) + getNT(j+1) + getNT(j+2);
				String next = getNT(j+3);
				if (result.containsKey(triplet)) {
					Score each = result.get(triplet);
					each.sum += 1;
					if (next.equals("A")) {
						each.A += 1;
					} else if (next.equals("C")) {
						each.C += 1;
					} else if (next.equals("G")) {
						each.G += 1;
					} else {
						each.T += 1;
					}	
				} else {
					Score each = new Score(0, 0, 0, 0, 1);
					if (next.equals("A")) {
						each.A += 1;
					} else if (next.equals("C")) {
						each.C += 1;
					} else if (next.equals("G")) {
						each.G += 1;
					} else {
						each.T += 1;
					}
					result.put(triplet, each);					
				}
			}
		}
	}
	
	public static String getNT(int j) {
		if (A.contains(j)) {
			return "A"; 
		} else if (C.contains(j)) {
			return "C";
		} else if (G.contains(j)) {
			return "G";
		} else {
			return "T";
		}
	}
	
	public static void report(boolean eachNuc, boolean ORF, boolean ORFsize) {
		if (eachNuc) {
			System.out.println("A: " + A.toString());
			System.out.println("C: " + C.toString());
			System.out.println("G: " + G.toString());
			System.out.println("T: " + T.toString());
		}
		if (ORF) {
			System.out.println("mod1: " + mod1.toString());
			System.out.println("mod2: " + mod2.toString());
			System.out.println("mod0: " + mod0.toString());
		}
		if (ORFsize) {
			System.out.println("mod1: " + mod1.size() + ", mod2: " + mod2.size() + ", mod0: " + mod0.size());
		}
	}

}
