import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;


public class CDSParser {
	
	public static void parseCDSFile(ArrayList<ORF> result, DNASequence seq, String output, int interested) {
		BufferedReader br = null;
		int newline = 0;
		try {
			String currentLine;
			br = new BufferedReader(new FileReader("data/" + File.separator + output));
			
			while ((currentLine = br.readLine()) != null) {
				if (currentLine.startsWith("LOCUS")) {
					if (newline == interested) {
						break;
					}
					newline += 1;
					continue;
				}
				if (newline == interested) {
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

	public static int markCdsOrfs(ORFCollection mod, ArrayList<ORF> cds) {
		mod.sortByFinish();
		int score = 0;
		int index = 0;
		for (int i = 0; i < mod.orfs.size(); i++) {
			if (mod.orfs.get(i).finish > cds.get(index).finish) {
				while(mod.orfs.get(i).finish > cds.get(index).finish) {
					index += 1;
					if (index >= cds.size()) return score;
				}
			}
			if (mod.orfs.get(i).finish == cds.get(index).finish) {
				mod.orfs.get(i).isCDS = true;
				score += 1;
				index += 1;
				if (index >= cds.size()) return score;
			}
		}
		return score;
	}
}
