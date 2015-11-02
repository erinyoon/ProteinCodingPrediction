import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;

/*  
 *  NOTE:
 *  Ignore the first line (the one starting with ">") and all whitespace
 *  Treat every other character as a nucleotide, case-insensitive
 *  For simplicity, just treat anything other than aAcCgG as a T.
 *  Again in the interest of simplicity, we are only going to look at the long chromosome
 *  */
public class DNASequence {

	HashSet<Integer> A = new HashSet<Integer>();
	HashSet<Integer> C = new HashSet<Integer>();
	HashSet<Integer> G = new HashSet<Integer>();
	HashSet<Integer> T = new HashSet<Integer>();
	int totalLength;
	
	public DNASequence (String input, int interested) {
		BufferedReader br = null;
		int index = 0;
		int newline = 0;
		try {
			String currentLine;
			br = new BufferedReader(new FileReader("data/" + input));
			
			while ((currentLine = br.readLine()) != null) {
				if (currentLine.startsWith(">")) {
					if (newline == interested) {
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
		totalLength = index;
	}
	
	public boolean isA(int j) {
		return A.contains(j);
	}
	
	public boolean isC(int j) {
		return C.contains(j);
	}
	
	public boolean isG(int j) {
		return G.contains(j);

	}
	
	public boolean isT(int j) {
		return T.contains(j);
	}
	
	// Getting corresponding nucleotide at index j (
	public String getNT(int j) {
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
