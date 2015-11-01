import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;


public class CDSParser {
	final static int interested_label = 1; // looking at first ">"
	final String input = "GCF_000091665.1_ASM9166v1_genomic.fna";
	final static String output = "GCF_000091665.1_ASM9166v1_genomic.gbff";
	
	public static void parse(ArrayList<ORF> result) {
		BufferedReader br = null;
		int index = 0;
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

	public static int scoreORF(ArrayList<ORF> mod, ArrayList<ORF> cds) {
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
}
