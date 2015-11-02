import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;


public class Predictor {
	
	final int interested_label = 1; // looking at first ">"
	
	public static void main(String[] args) {
		
		// MASSIVE PARAMETERS
		// mod, orfts, short_orfs, all share same reference
		HashSet<Integer> A = new HashSet<Integer>();
		HashSet<Integer> C = new HashSet<Integer>();
		HashSet<Integer> G = new HashSet<Integer>();
		HashSet<Integer> T = new HashSet<Integer>();
		
		ArrayList<ORF> mod0 = new ArrayList<ORF>();
		ArrayList<ORF> mod1 = new ArrayList<ORF>();
		ArrayList<ORF> mod2 = new ArrayList<ORF>();
		
		ArrayList<ORF> orfs = new ArrayList<ORF>();
		ArrayList<ORF> short_orfs = new ArrayList<ORF>();
		ArrayList<ORF> long_orfs = new ArrayList<ORF>();
		ArrayList<ORF> cds = new ArrayList<ORF>();

		// String is three letters (k), and score contains count of each letter k+1 
		HashMap<String, Score> P = new HashMap<String, Score>(); 
		HashMap<String, Score> Q = new HashMap<String, Score>();
		
		
		// PART 1: Reading files and parsing
		ORFParser.findORF(A, C, G, T, mod1, mod2, mod0, orfs); // find all ORFS in three reading frames

		System.out.println("mod1: " + mod1.size() + ", mod2: " + mod2.size() + ", mod0: " + mod0.size());
		
		ORFParser.shorterORFs(orfs, short_orfs, 50); // process ORFs
		ORFParser.longerORFs(orfs, long_orfs, 1400);
		System.out.println("shorter: " + short_orfs.size() + ", longer: " + long_orfs.size());
		
		CDSParser.parseCDSFile(cds); // find all CDS & matching
		int totalCDS = CDSParser.simpleScoreORF(orfs, cds);
		System.out.println("score: " + totalCDS);
		
		
		// PART 1: trying all threshold
//		int max_threshold = 1400; // skip by 100
//		
//		// Length threshold
//		for (int i = 100; i < 1500; i+= 100) {
//			int[] counts = new int[4];
//			lengthThresholdORF(orfs, i);
//		}
//		

//		 PART 2: Calculating transition & emission probabilities

		setupProb(long_orfs, P, A, C, G, T);
		setupProb(short_orfs, Q, A, C, G, T);
		
		
		//answerQuestion1e(P, Q); // For question 1e
//		
		Collections.sort(short_orfs, start_comparator);
		Collections.sort(long_orfs, start_comparator);
		CDSParser.markovScoreORF(short_orfs, P, Q, A, C, G, T);	
		CDSParser.markovScoreORF(long_orfs, P, Q, A, C, G, T);	

		answerQuestion1f(short_orfs, long_orfs);
//			
//		for (double i = -10; i < 10; i+= 1.0) {
//			markovThresholdORF(orfs, i);
//		}
	}
	
	private static void markovThresholdORF(ArrayList<ORF> orfs, double threshold) {
		int true_positive = 0;
		int false_positive = 0;
		int true_negative = 0;
		int false_negative = 0;
		
		for (int i = 0; i < orfs.size(); i++) {
			if (orfs.get(i).markov_score > threshold && orfs.get(i).isCDS) {
				true_positive += 1;
			} else if (orfs.get(i).markov_score > threshold && !orfs.get(i).isCDS) {
				false_positive += 1;
			} else if (orfs.get(i).markov_score <= threshold && orfs.get(i).isCDS) {
				false_negative += 1;
			} else if (orfs.get(i).markov_score <= threshold && !orfs.get(i).isCDS) {
				true_negative += 1;
			}
		}
		System.out.println("false negative; " + false_negative);
		System.out.println("++: " + 1.0 * true_positive / (true_positive + false_negative) +
						  " +-: " + 1.0 * false_positive/ (false_positive + true_negative));
	}
	
	private static void lengthThresholdORF(ArrayList<ORF> orfs, int threshold) {
		int true_positive = 0;
		int false_positive = 0;
		int true_negative = 0;
		int false_negative = 0;
		
		for (int i = 0; i < orfs.size(); i++) {
			if (orfs.get(i).length > threshold && orfs.get(i).isCDS) {
				true_positive += 1;
			} else if (orfs.get(i).length > threshold && !orfs.get(i).isCDS) {
				false_positive += 1;
			} else if (orfs.get(i).length <= threshold && orfs.get(i).isCDS) {
				false_negative += 1;
			} else if (orfs.get(i).length <= threshold && !orfs.get(i).isCDS) {
				true_negative += 1;
			}
		}
		System.out.println("++: " + 1.0 * true_positive / (true_positive + false_negative) +
						  " +-: " + 1.0 * false_positive/ (false_positive + true_negative));
	}
	
	private static void answerQuestion1f(ArrayList<ORF> short_orfs, ArrayList<ORF> long_orfs) {
		Collections.sort(short_orfs, start_comparator);
		Collections.sort(long_orfs, start_comparator);
		//its start coordinate, length, log-base-2 Markov model score
		// and a flag indicating whether this ORF was/was not among the "simple forward strand CDSs" in GenBan
		
		System.out.println("Start\tLength\tIS_CDS\tScore");
		for (int i = 0; i < 5; i++) {
			System.out.println(short_orfs.get(i).toString());
		}
		for (int i = 0; i < 5; i++) {
			System.out.println(long_orfs.get(i).toString());
		}
	}
	
	static Comparator<ORF> start_comparator = new Comparator<ORF>() {
		@Override
		public int compare(ORF o1, ORF o2) {
			if (o1.start > o2.start) {
				return 1;
			} else if (o1.finish == o2.finish) {
				return 0;
			} else {
				return -1;
			}
		}
	};
	/* HELPER METHODS */
	/* methods for creating/using Markov model */
	// Getting corresponding nucleotide at index j (
	private static String getNT(int j, HashSet<Integer> A, HashSet<Integer> C, HashSet<Integer> G, HashSet<Integer> T) {
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
	
	// Calculate the P/Q by counting
	public static void setupProb(ArrayList<ORF> orfs, HashMap<String, Score> result,
								HashSet<Integer> A, HashSet<Integer> C, HashSet<Integer> G, HashSet<Integer> T){
		for (int i = 0; i < orfs.size(); i++) {
			int start = orfs.get(i).start;
			int finish = orfs.get(i).finish;
			
			for (int j = start; j + 3 < finish - 2; j++) {
				String triplet = getNT(j, A, C, G, T) + getNT(j+1, A, C, G, T) + getNT(j+2, A, C, G, T);
				String next = getNT(j+3, A, C, G, T);
				if (!result.containsKey(triplet)) {
					result.put(triplet, new Score(0,0,0,0));
				}
				if (next.equalsIgnoreCase("A")) {
					result.get(triplet).A += 1;
				} else if (next.equalsIgnoreCase("C")) {
					result.get(triplet).C += 1;
				} else if (next.equalsIgnoreCase("G")) {
					result.get(triplet).G += 1;
				} else {
					result.get(triplet).T += 1;
				}	
			}
		}
	}
	
	
	/* PRIVATE METHODS */
	/* Differs from public methods; not essential part of analysis */
	// printing out P & Q for Question 1e
	private static void answerQuestion1e(HashMap<String, Score> P, HashMap<String, Score> Q) {
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
}
