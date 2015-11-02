import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;


public class ORFCollection {
	
	ArrayList<ORF> orfs = new ArrayList<ORF>();
	ArrayList<ORF> short_orfs = new ArrayList<ORF>();
	ArrayList<ORF> long_orfs = new ArrayList<ORF>();
	int[] mods = new int[3];
	
	public ORFCollection (DNASequence seq, int min, int max) { 
		for (int i = 1; i <= 3; i++) {
			// stop codon (TAA, TAG, or TGA).
			boolean started = false;
			int temp = i;
			for (int j = i; j + 2 <= seq.totalLength; j+= 3) {
				if ((seq.isT(j) && seq.isA(j+1) && (seq.isA(j+2) || seq.isG(j+2))) ||
					(seq.isT(j) && seq.isG(j+1) && seq.isA(j+2))) {
					if (started) { // seeing more than once
						ORF eachORF = new ORF(temp, j + 2);
						if (i == 1) {
							mods[1] += 1;
						} else if (i == 2) {
							mods[2] += 1;
						} else {
							mods[0] += 1;
						}
						orfs.add(eachORF);
						temp = j + 3;
					} else { // seeing first time
						ORF eachORF = new ORF(i, j + 2);
						if (i == 1) {
							mods[1] += 1;
						} else if (i == 2) {
							mods[2] += 1;
						} else {
							mods[0] += 1;
						}
						orfs.add(eachORF);
						temp = j + 3;
						started = true;
					}
				}
			}
		}
		System.out.println("mod1: " + mods[1] + ", mod2: " + mods[2] + ", mod0: " + mods[0]);
		
		for (int i = 0; i < orfs.size(); i++) {
			if (orfs.get(i).length < min) {
				short_orfs.add(orfs.get(i));
			}
		}
		
		for (int i = 0; i < orfs.size(); i++) {
			if (orfs.get(i).length > max) {
				long_orfs.add(orfs.get(i));
			}
		}
		System.out.println("shorter: " + short_orfs.size() + ", longer: " + long_orfs.size());
	}
	
	public void sortByFinish() {
		Collections.sort(short_orfs, finish_comparator);
		Collections.sort(long_orfs, finish_comparator);
		Collections.sort(orfs, finish_comparator);
	}
	
	public void sortByStart() {
		Collections.sort(short_orfs, start_comparator);
		Collections.sort(long_orfs, start_comparator);
		Collections.sort(orfs, start_comparator);	}
	
	/** INNER CLASSES **/
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
	
	static Comparator<ORF> finish_comparator = new Comparator<ORF>() {
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

}
