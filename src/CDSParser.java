import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;


public class CDSParser {
	final static int interested_label = 1; // looking at first ">"
	final static String output = "GCF_000091665.1_ASM9166v1_genomic.gbff";

	public static void parseCDSFile(ArrayList<ORF> result) {
		BufferedReader br = null;
		int newline = 0;
		try {
			String currentLine;
			br = new BufferedReader(new FileReader("data/" + File.separator + output));
			
			while ((currentLine = br.readLine()) != null) {
				if (currentLine.startsWith("LOCUS")) {
					if (newline >= interested_label) {
						break;
					}
					newline += 1;
					continue;
				}
				if (newline == interested_label) {
					currentLine = currentLine.trim();
					if (currentLine.startsWith("CDS") && !currentLine.contains("complement")) {
						String[] arr = currentLine.split("\\s+");
						String[] pos = arr[1].split("\\.\\.");
						ORF eachORF = new ORF(Integer.parseInt(pos[0]), Integer.parseInt(pos[1]));
						result.add(eachORF);
					}
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (br != null)br.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
	}

	public static int simpleScoreORF(ArrayList<ORF> mod, ArrayList<ORF> cds) {
		int score = 0;
		int index = 0;
		for (int i = 0; i < mod.size(); i++) {
			if (mod.get(i).finish > cds.get(index).finish) {
				while(mod.get(i).finish > cds.get(index).finish) {
					index += 1;
					if (index >= cds.size()) return score;
				}
			}
			if (mod.get(i).finish == cds.get(index).finish) {
				mod.get(i).isCDS = true;
				score += 1;
				index += 1;
				if (index >= cds.size()) return score;
			}
		}
		return score;
	}
	
	public static void markovScoreORF(ArrayList<ORF> orfs, HashMap<String, Score> P, HashMap<String, Score> Q, 
									  HashSet<Integer> A, HashSet<Integer> C, HashSet<Integer> G, HashSet<Integer> T) {

		for (int i = 0; i < orfs.size(); i++) {
			// for each sequence, compute below
			// i.e) P(x) = P(x1 x2 x3) * P(x4 | x1 x2 x3) * P(x5 | x2 x3 x4) ... P(xn | xn-3 xn-2 xn-1)
			ORF each = orfs.get(i);
			
			if (each.length <= 3) {
				each.markov_score = 0;
				continue;
			}
			String x1x2x3 = getNT(each.start, A, C, G, T) + getNT(each.start + 1, A, C, G, T) + getNT(each.start + 2, A, C, G, T); // P(x1, x2, x3)			
			each.markov_score += Math.log(1.0 * P.get(x1x2x3).getTotal() / (each.length - 5)) 
								- Math.log(1.0 * Q.get(x1x2x3).getTotal() / (each.length - 5));

			for (int j = each.start; j < each.finish - 5 ; j++) {
				String conditional = getNT(j, A, C, G, T) + getNT(j + 1, A, C, G, T) + getNT(j + 2, A, C, G, T);
				String emitted = getNT(j + 3, A, C, G, T);
				each.markov_score += (Math.log(P.get(conditional).getProb(emitted)) - Math.log(Q.get(conditional).getProb(emitted)));
			}
			each.markov_score /= Math.log(2);
		}
	}
	
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
}
