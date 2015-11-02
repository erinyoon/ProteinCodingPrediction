import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;


public class Predictor {
	
	final static int interested_label = 1; // looking at first ">"
	final static String input = "GCF_000091665.1_ASM9166v1_genomic.fna";
	final static String output = "GCF_000091665.1_ASM9166v1_genomic.gbff";
	
	public static void main(String[] args) {
		
		DNASequence seq = new DNASequence(input, interested_label);
		ORFCollection orfs = new ORFCollection(seq, 50, 1400);
		ArrayList<ORF> cds = new ArrayList<ORF>();

		// PART 1: Reading files and parsing
		CDSParser.parseCDSFile(cds, seq, output, interested_label); // find all CDS & matching
		int totalCDS = CDSParser.markCdsOrfs(orfs, cds);
		System.out.println("score: " + totalCDS);
		
		// Part 2: Calculating Markov scores
		MarkovModel markovScores = new MarkovModel(orfs, seq);
		//answerQuestion1e(markovScores); // For question 1e, printing out some scores with A
		answerQuestion1f(orfs);			// For question 1f, printing out long & short 5 ORFs
		
		// PART 3.1: Length threshold
		//for (int i = 100; i <= 1500; i+= 100) {
		//	lengthThresholdORF(orfs, i);
		//}
		// PART 3.2: Markov threshold
		for (double i = 10; i <= 100; i+= 10.0) {
			markovThresholdORF(orfs, i);
		}
	}
	
	private static void markovThresholdORF(ORFCollection orfs, double threshold) {
		int true_positive = 0;
		int false_positive = 0;
		int true_negative = 0;
		int false_negative = 0;
		
		for (int i = 0; i < orfs.orfs.size(); i++) {
			if (orfs.orfs.get(i).markov_score > threshold && orfs.orfs.get(i).isCDS) {
				true_positive += 1;
			} else if (orfs.orfs.get(i).markov_score > threshold && !orfs.orfs.get(i).isCDS) {
				false_positive += 1;
			} else if (orfs.orfs.get(i).markov_score <= threshold && orfs.orfs.get(i).isCDS) {
				false_negative += 1;
			} else if (orfs.orfs.get(i).markov_score <= threshold && !orfs.orfs.get(i).isCDS) {
				true_negative += 1;
			}
		}
		System.out.println("MARKOV: "  + threshold +
						   ", tpr: " + 1.0 * true_positive / (true_positive + false_negative) +
						   ", fpr: " + 1.0 * false_positive/ (false_positive + true_negative));
	}
	
	private static void lengthThresholdORF(ORFCollection orfs, int threshold) {
		int true_positive = 0;
		int false_positive = 0;
		int true_negative = 0;
		int false_negative = 0;
		
		for (int i = 0; i < orfs.orfs.size(); i++) {
			if (orfs.orfs.get(i).length > threshold && orfs.orfs.get(i).isCDS) {
				true_positive += 1;
			} else if (orfs.orfs.get(i).length > threshold && !orfs.orfs.get(i).isCDS) {
				false_positive += 1;
			} else if (orfs.orfs.get(i).length <= threshold && orfs.orfs.get(i).isCDS) {
				false_negative += 1;
			} else if (orfs.orfs.get(i).length <= threshold && !orfs.orfs.get(i).isCDS) {
				true_negative += 1;
			}
		}
		System.out.println("LENGTH: "  + threshold +
				   		   ", tpr: " + 1.0 * true_positive / (true_positive + false_negative) +
				           ", fpr: " + 1.0 * false_positive/ (false_positive + true_negative));
	}
	
	private static void answerQuestion1f(ORFCollection orfs) {
		orfs.sortByStart();
		System.out.println("Start\tLength\tIS_CDS\tScore");
		for (int i = 0; i < 5; i++) {
			System.out.println(orfs.short_orfs.get(i).toString());
		}
		for (int i = 0; i < 5; i++) {
			System.out.println(orfs.long_orfs.get(i).toString());
		}
	}
	

	/* HELPER METHODS */
	/* methods for creating/using Markov model */
	// Calculate the P/Q by counting
	
	
	
	/* PRIVATE METHODS */
	/* Differs from public methods; not essential part of analysis */
	// printing out P & Q for Question 1e
	private static void answerQuestion1e(MarkovModel mm) {
		String[] temp = new String[]{"A", "C", "G", "T"};
		System.out.println("\tA\t\tC\t\tG\t\tT");
		for (int i = 0; i < 4; i++) {
			System.out.print(temp[i] + "\t");
			for (int j = 0; j < 4; j++) {
				String key = "A" + temp[i] + temp[j];
				System.out.print(mm.P.get(key).A);
				System.out.print("|");
				System.out.print(mm.Q.get(key).A);
				System.out.print("\t");
			}
			System.out.println();
		}
	}
}
