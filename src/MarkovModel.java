import java.util.ArrayList;
import java.util.HashMap;

public class MarkovModel {

	HashMap<String, Score> P = new HashMap<String, Score>(); 
	HashMap<String, Score> Q = new HashMap<String, Score>();
	
	public MarkovModel(ORFCollection orfs, DNASequence seq) {	
		setupProbHelp(orfs.long_orfs, seq, P);
		setupProbHelp(orfs.short_orfs, seq, Q);
		markovScoreORF(orfs, seq);
	}	
	
	private void markovScoreORF(ORFCollection orfs, DNASequence seq) {

		int lengthP = 0;
		int lengthQ = 0;

		for (int i = 0; i < orfs.orfs.size(); i++) {
			if (orfs.orfs.get(i).length < 50) {
				lengthQ += orfs.orfs.get(i).length - 5;
			} else if (orfs.orfs.get(i).length > 1400) {
				lengthP += orfs.orfs.get(i).length - 5;
			}
		}

		System.out.println("P " + lengthP + " Q " + lengthQ);

		for (int i = 0; i < orfs.orfs.size(); i++) {
			// for each sequence, compute below
			// i.e) P(x) = P(x1 x2 x3) * P(x4 | x1 x2 x3) * P(x5 | x2 x3 x4) ... P(xn | xn-3 xn-2 xn-1)
			ORF each = orfs.orfs.get(i);

			if (each.length <= 3) {
				each.markov_score = 0;
				continue;
			}
			
			String x1x2x3 = seq.getNT(each.start) + seq.getNT(each.start + 1) + seq.getNT(each.start + 2); // P(x1, x2, x3)			
			each.markov_score += Math.log(1.0 * P.get(x1x2x3).getTotal() / lengthP)
			- Math.log(1.0 * Q.get(x1x2x3).getTotal() / lengthQ);

			for (int j = each.start; j < each.finish - 5 ; j++) {
				String conditional = seq.getNT(j) + seq.getNT(j + 1) + seq.getNT(j + 2);
				String emitted = seq.getNT(j + 3);
				each.markov_score += (Math.log(P.get(conditional).getProb(emitted)) - Math.log(Q.get(conditional).getProb(emitted)));
			}
			each.markov_score /= Math.log(2);
		}
	}

	private void setupProbHelp(ArrayList<ORF> orfs, DNASequence seq, HashMap<String, Score> result){
		for (int i = 0; i < orfs.size(); i++) {
			int start = orfs.get(i).start;
			int finish = orfs.get(i).finish;
	
			for (int j = start; j + 3 < finish - 2; j++) {
				String triplet = seq.getNT(j) + seq.getNT(j+1) + seq.getNT(j+2);
				String next = seq.getNT(j+3);
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
	
	public class Score {
		int A; // each seen
		int C;
		int G;
		int T;
		
		public Score(int a, int c, int g, int t) {
			A = a; C = c; G = g; T = t;
		}
		
		public int getTotal() {
			return A + C + G + T;
		}
		
		public double getProb(String temp) {
			if (temp.equalsIgnoreCase("a")) {
				return 1.0 * A / getTotal();
			} else if (temp.equalsIgnoreCase("c")) {
				return 1.0 * C / getTotal();
			} else if (temp.equalsIgnoreCase("g")) {
				return 1.0 * G / getTotal();
			} else {
				return 1.0 * T / getTotal();
			}
		}
		
		public String toString() {
			return "A:" + A + " C:" + C + " G:" + G + " T:" + T;
		}
	}

}
