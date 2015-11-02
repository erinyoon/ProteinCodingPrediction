
public class ORF {
	int start;
	int finish;
	int length;
	boolean isCDS = false;
	boolean markov_isCDS = false;
	double markov_score = 0;
	
	public ORF(int start, int finish) {
		this.start = start;
		this.finish = finish;
		this.length = (finish - start + 1);
	}
	
	public String toString() {
		//its start coordinate, length, log-base-2 Markov model score
		// and a flag indicating whether this ORF was/was not among the "simple forward strand CDSs" in GenBan
		return start + "\t" + length + "\t" + isCDS + "\t" + markov_score;
	}
}
