import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;


/**
 * ORFParser takes in file of Methanocaldococcus jannaschii
 *  
 *  NOTE:
 *  Ignore the first line (the one starting with ">") and all whitespace
 *  Treat every other character as a nucleotide, case-insensitive
 *  For simplicity, just treat anything other than aAcCgG as a T.
 *  Again in the interest of simplicity, we are only going to look at the long chromosome
 *
 */
public class ORFParser {
	final static String input = "GCF_000091665.1_ASM9166v1_genomic.fna";
	//final static String input = "dumb";
	
	public static int parse(HashSet<Integer> A, HashSet<Integer> C, HashSet<Integer> G, HashSet<Integer> T) {
		BufferedReader br = null;
		int index = 0;
		int newline = 0;
		try {
			String currentLine;
			br = new BufferedReader(new FileReader("data/" + input));
			
			while ((currentLine = br.readLine()) != null) {
				if (currentLine.startsWith(">")) {
					if (newline >= 1) {
						break;
					}
					newline += 1;
					continue;
				}
				String[] arr = currentLine.split("");
				// ignore first one
				for (int i = 1; i < arr.length; i++) {
					String value = arr[i].toLowerCase();
					if (value.equals("a")) {
						A.add(i + index);
					} else if (value.equals("c")) {
						C.add(i + index);
					} else if (value.equals("g")) {
						G.add(i + index);
					} else { // all others are T
						T.add(i + index);
					}
				}
				index += arr.length - 1;
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
		return index;
	}
	
	public static void findORF(HashSet<Integer> A, HashSet<Integer> C, HashSet<Integer> G, HashSet<Integer> T,
								ArrayList<ORF> mod1, ArrayList<ORF> mod2, ArrayList<ORF> mod0, ArrayList<ORF> orfs) {
		int size = parse(A, C, G, T);
		for (int i = 1; i <= 3; i++) {
			// mod 1, 2, 0
			// stop codon (TAA, TAG, or TGA).
			boolean started = false;
			int temp = i;
			for (int j = i; j + 2 <= size; j+= 3) {
				if ((T.contains(j) && A.contains(j+1) && (A.contains(j+2) || G.contains(j+2))) ||
					(T.contains(j) && G.contains(j+1) && A.contains(j+2))) {
					if (started) { // seeing more than once
						ORF eachORF = new ORF(temp, j + 2);
						if (i == 1) {
							mod1.add(eachORF);
						} else if (i == 2) {
							mod2.add(eachORF);
						} else {
							mod0.add(eachORF);
						}
						orfs.add(eachORF);
						temp = j + 3;
					} else { // seeing first time
						ORF eachORF = new ORF(i, j + 2);
						if (i == 1) {
							mod1.add(eachORF);
						} else if (i == 2) {
							mod2.add(eachORF);
						} else {
						    mod0.add(eachORF);
						}
						orfs.add(eachORF);
						temp = j + 3;
						started = true;
					}
				}
			}
		}
		Collections.sort(orfs, orf_comparator);
	}
	
	static Comparator<ORF> orf_comparator = new Comparator<ORF>() {
		@Override
		public int compare(ORF o1, ORF o2) {
			if (o1.finish > o2.finish) {
				return 1;
			} else if (o1.finish == o2.finish) {
				return 0;
			} else {
				return -1;
			}
		}
	};

	public static void shorterORFs(ArrayList<ORF> arr, ArrayList<ORF> res, int max) {
		for (int i = 0; i < arr.size(); i++) {
			if (arr.get(i).length < max) {
				res.add(arr.get(i));
			}
		}
	}

	public static void longerORFs(ArrayList<ORF> arr, ArrayList<ORF> res, int min) {
		for (int i = 0; i < arr.size(); i++) {
			if (arr.get(i).length > min) {
				res.add(arr.get(i));
			}
		}
	}
}
